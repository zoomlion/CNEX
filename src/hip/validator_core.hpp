// MIT License
//
// Copyright (c) 2024 JiangminZheng
//
// Core validation logic - shared between Python module and standalone binary.
// No pybind11 dependency.

#pragma once
#include <tsl/robin_map.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Optional counters for the segmented genome-mode validator. Enabled by
// setting the VALIDATE_PROFILE environment variable to a non-zero value.
struct ValidateProfile {
    std::atomic<long long> windows{0};
    std::atomic<long long> gathers{0};
    std::atomic<long long> emits{0};
    std::atomic<long long> pushed{0};
    std::atomic<long long> tally_calls{0};
    std::atomic<long long> chain_calls{0};
    std::atomic<long long> prep_ns{0};
    std::atomic<long long> win_ns{0};
    std::atomic<long long> gather_ns{0};
    std::atomic<long long> tally_ns{0};
};
inline ValidateProfile g_val_profile;

inline bool validate_profile_enabled() {
    static const bool on = [] {
        const char* v = std::getenv("VALIDATE_PROFILE");
        return v != nullptr && v[0] != '\0' && v[0] != '0';
    }();
    return on;
}

// ─── 2-bit DNA encoding ───

inline uint32_t dna_encoder(const std::string& dna) {
    static const std::array<uint32_t, 256> lookup = [] {
        std::array<uint32_t, 256> table = {};
        table['A'] = 0; table['C'] = 1; table['G'] = 2; table['T'] = 3;
        return table;
    }();

    if (dna.length() > 16) {
        throw std::overflow_error("Input DNA string is too long.");
    }

    uint32_t result = 0;
    for (char c : dna) {
        if (lookup[(unsigned char)c] == 0 && c != 'A') {
            throw std::invalid_argument("Invalid DNA character: " + std::string(1, c));
        }
        result = (result << 2) | lookup[(unsigned char)c];
    }
    return result;
}


inline bool encode_mer_at(const std::string& seq, size_t pos, size_t k, uint32_t& result) {
    static constexpr std::array<int8_t, 256> lookup = [] {
        std::array<int8_t, 256> table = {};
        for (auto& v : table) v = -1;
        table['A'] = 0; table['C'] = 1; table['G'] = 2; table['T'] = 3;
        return table;
    }();

    result = 0;
    for (size_t i = 0; i < k; ++i) {
        int8_t bits = lookup[static_cast<unsigned char>(seq[pos + i])];
        if (bits < 0) return false;
        result = (result << 2) | static_cast<uint32_t>(bits);
    }
    return true;
}


inline void reverse_complement(const std::string& dna, std::string& out) {
    static const std::array<char, 256> complement = [] {
        std::array<char, 256> table = {};
        for (auto& v : table) v = 'N';
        table['A'] = 'T'; table['T'] = 'A';
        table['C'] = 'G'; table['G'] = 'C';
        table['a'] = 't'; table['t'] = 'a';
        table['c'] = 'g'; table['g'] = 'c';
        table['N'] = 'N'; table['n'] = 'n';
        return table;
    }();

    out.clear();
    out.reserve(dna.length());
    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        out.push_back(complement[(unsigned char)*it]);
    }
}




struct VoteTally {
    int best_id = -1;
    int best_cnt = 0;
    int second_cnt = 0;
    int n = 0;
};

// O(m) tally with epoch stamps: element ids are bounded, no sort/clear per window.
template <typename IdAt>
inline VoteTally tally_votes_impl(int n, IdAt&& id_at) {
    VoteTally r;
    r.n = n;
    if (n == 0) return r;

    thread_local std::vector<unsigned> stamps;
    thread_local std::vector<int> counts;
    thread_local std::vector<int> seen;
    thread_local unsigned epoch = 0;

    ++epoch;
    if (epoch == 0) {
        std::fill(stamps.begin(), stamps.end(), 0u);
        epoch = 1;
    }
    seen.clear();

    for (int i = 0; i < n; ++i) {
        int id = id_at(i);
        if (id >= static_cast<int>(stamps.size())) {
            stamps.resize(id + 1, 0u);
            counts.resize(id + 1, 0);
        }
        if (stamps[id] != epoch) {
            stamps[id] = epoch;
            counts[id] = 1;
            seen.push_back(id);
        } else {
            ++counts[id];
        }
    }

    for (int id : seen) {
        int c = counts[id];
        if (c > r.best_cnt) {
            r.second_cnt = r.best_cnt;
            r.best_cnt = c;
            r.best_id = id;
        } else if (c > r.second_cnt) {
            r.second_cnt = c;
        }
    }
    return r;
}

inline VoteTally tally_votes(const std::vector<int>& nums) {
    return tally_votes_impl(static_cast<int>(nums.size()),
                            [&](int i) { return nums[i]; });
}

// Gate matrix: T1 = legacy defaults (from the flags); T2, T3, ... are extra
// OR-combined rows. A window is accepted when ANY row passes. Rows share one
// tally and one chain evaluation per window.
struct GateTier {
    double vote_ratio = 3.0;
    int min_c = 7;
    int min_span = 25;
};

// In-validate rescue: windows/reads that fail every strict tier but reach the
// slightly looser T3 row are verified with SIMD Smith-Waterman against the
// element's panel(s); a normalized score >= tau emits them as validated
// evidence (same downstream assembly path as the strict output).
struct VerifyConfig {
    bool enabled = false;
    GateTier t3{1.5, 5, 25};
    double tau = 0.9;
};

struct ValidatorConfig {
    int window_size = 150;
    int step_size = 25;
    double vote_frac = 0.1;
    int max_diff = 2;
    std::vector<GateTier> tiers = {GateTier{}};   // T1 OR T2 ...
    VerifyConfig verify;
};

// chain_fn(best_id) -> (passing, span); invoked at most once per window.
inline bool gate_ratio_precheck(const VoteTally& t, const ValidatorConfig& cfg) {
    if (t.best_id < 0 || cfg.tiers.empty()) return false;
    float threshold = static_cast<float>(t.n) * static_cast<float>(cfg.vote_frac);
    if (t.best_cnt <= threshold) return false;
    for (const auto& tier : cfg.tiers) {
        if (t.best_cnt >= static_cast<int>(tier.vote_ratio * t.second_cnt)) return true;
    }
    return false;
}

template <typename ChainFn>
inline bool gate_pass(const VoteTally& t, const ValidatorConfig& cfg, ChainFn&& chain_fn,
                      int& passing_out, int& span_out) {
    passing_out = 0;
    span_out = 0;
    if (!gate_ratio_precheck(t, cfg)) return false;

    // Integer truncation matches the historical T1 semantics; keep it even
    // for fractional ratios (e.g. 1.5) to avoid changing accepted windows.
    auto ratio_ok = [&](const GateTier& tier) {
        return t.best_cnt >= static_cast<int>(tier.vote_ratio * t.second_cnt);
    };

    auto cs = chain_fn(t.best_id);
    passing_out = cs.first;
    span_out = cs.second;

    for (const auto& tier : cfg.tiers) {
        if (ratio_ok(tier) && passing_out >= tier.min_c && span_out >= tier.min_span) {
            return true;
        }
    }
    return false;
}

// T3 rescue gate: evaluated only when the strict rows failed. Uses the strict
// concentration threshold (vote_frac) and the relaxed T3 chain/dominance row.
template <typename ChainFn>
inline std::pair<int, int> rescue_gate(const VoteTally& t, const ValidatorConfig& cfg,
                                       ChainFn&& chain_fn) {
    if (!cfg.verify.enabled || t.best_id < 0) return {-1, 0};
    const double threshold = static_cast<double>(cfg.vote_frac) * t.n;
    if (!(static_cast<double>(t.best_cnt) > threshold)) return {-1, 0};
    if (t.best_cnt < static_cast<int>(cfg.verify.t3.vote_ratio * t.second_cnt)) return {-1, 0};
    auto cs = chain_fn(t.best_id);
    if (cs.first < cfg.verify.t3.min_c || cs.second < cfg.verify.t3.min_span) return {-1, 0};
    return {t.best_id, t.best_cnt};
}


// ─── MerQueryManager ───

class MerQueryManager {
public:
    // Packed assignment ranges into assign_arena:
    //   high 32 bits = offset, low 32 bits = number of (element, locus) assignments.
    // All rows of one k-mer are contiguous (mertable dump order), so loading keeps
    // a flat arena instead of one heap vector per k-mer.
    tsl::robin_map<uint32_t, uint64_t> compressed_mer_query;
    std::vector<std::pair<int, int>> assign_arena;
    int mer_size = 0;

    // Optional direct-index accelerator for the genome-mode prep pass
    // (k <= 13): maps each 2-bit k-mer code to a slot holding rank+1, with
    // the packed range stored in direct_rng[rank]. Costs 4^k * 4 bytes.
    std::vector<uint32_t> direct_index;
    std::vector<uint64_t> direct_rng;

    bool has_direct_index() const { return !direct_index.empty(); }

    void build_direct_index() {
        if (mer_size <= 0 || mer_size > 13) return;
        direct_index.assign(size_t(1) << (2 * mer_size), 0);
        direct_rng.resize(compressed_mer_query.size());
        uint32_t rank = 0;
        for (const auto& kv : compressed_mer_query) {
            direct_index[kv.first] = rank + 1;
            direct_rng[rank] = kv.second;
            ++rank;
        }
    }

    inline uint64_t direct_query(uint32_t code) const {
        uint32_t r = direct_index[code];
        return r ? direct_rng[r - 1] : 0;
    }

    static inline uint64_t pack_range(uint64_t off, uint64_t cnt) {
        return (off << 32) | (cnt & 0xffffffffULL);
    }
    static inline uint32_t range_off(uint64_t v) { return static_cast<uint32_t>(v >> 32); }
    static inline uint32_t range_cnt(uint64_t v) { return static_cast<uint32_t>(v & 0xffffffffULL); }

    inline bool has_id(uint64_t range, int id) const {
        uint32_t off = range_off(range);
        uint32_t cnt = range_cnt(range);
        for (uint32_t i = 0; i < cnt; ++i) {
            if (assign_arena[off + i].first == id) return true;
        }
        return false;
    }

    void add_mer(const std::string& mer, int id, int loci) {
        uint32_t encoded_mer = dna_encoder(mer);
        auto it = compressed_mer_query.find(encoded_mer);
        if (it == compressed_mer_query.end()) {
            compressed_mer_query[encoded_mer] = pack_range(assign_arena.size(), 1);
        } else {
            compressed_mer_query[encoded_mer] =
                pack_range(assign_arena.size(), range_cnt(it->second) + 1);
        }
        assign_arena.push_back({id, loci});
        if (mer_size == 0) mer_size = static_cast<int>(mer.length());
    }

    // Load a mers table. Rows are grouped per k-mer (dump order). When
    // max_assign > 0, groups with more than max_assign (element, locus)
    // assignments are dropped at load time (repeat/paralog filter); 0 = keep all.
    void load_from_file(const std::string& path, int max_assign = 0) {
        std::ifstream file(path);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open mers file: " + path);
        }
        file.seekg(0, std::ios::end);
        size_t file_size = static_cast<size_t>(file.tellg());
        file.seekg(0, std::ios::beg);
        compressed_mer_query.reserve(file_size / 25);
        assign_arena.reserve(file_size / 20);

        std::string line;
        uint32_t cur_key = 0;
        bool have_cur = false;
        std::vector<std::pair<int, int>> group;
        group.reserve(64);
        long long kept_mers = 0, kept_rows = 0, dropped_mers = 0, dropped_rows = 0;
        long long shared_mers = 0;

        auto flush_group = [&]() {
            if (!have_cur) return;
            const int sz = static_cast<int>(group.size());
            if (sz > 1) ++shared_mers;
            if (max_assign == 0 || sz <= max_assign) {
                const size_t off = assign_arena.size();
                assign_arena.insert(assign_arena.end(), group.begin(), group.end());
                compressed_mer_query[cur_key] = pack_range(off, group.size());
                ++kept_mers;
                kept_rows += sz;
            } else {
                ++dropped_mers;
                dropped_rows += sz;
            }
            group.clear();
        };

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            const char* p = line.data();
            const char* e = p + line.size();
            const char* t1 = static_cast<const char*>(memchr(p, '\t', e - p));
            if (t1 == nullptr) continue;
            const char* t2 = static_cast<const char*>(memchr(t1 + 1, '\t', e - t1 - 1));
            if (t2 == nullptr) continue;
            const char* t3 = static_cast<const char*>(memchr(t2 + 1, '\t', e - t2 - 1));
            const char* loci_end = (t3 != nullptr) ? t3 : e;
            if (mer_size == 0) mer_size = static_cast<int>(t1 - p);

            uint32_t encoded_mer = 0;
            bool valid = true;
            for (const char* q = p; q < t1; ++q) {
                switch (*q) {
                    case 'A': encoded_mer <<= 2; break;
                    case 'C': encoded_mer = (encoded_mer << 2) | 1u; break;
                    case 'G': encoded_mer = (encoded_mer << 2) | 2u; break;
                    case 'T': encoded_mer = (encoded_mer << 2) | 3u; break;
                    default: valid = false; break;
                }
                if (!valid) break;
            }
            if (!valid) continue;
            int id = 0;
            for (const char* q = t1 + 1; q < t2; ++q) id = id * 10 + (*q - '0');
            int loci = 0;
            for (const char* q = t2 + 1; q < loci_end; ++q) loci = loci * 10 + (*q - '0');
            if (!have_cur || encoded_mer != cur_key) {
                flush_group();
                cur_key = encoded_mer;
                have_cur = true;
            }
            group.emplace_back(id, loci);
        }
        flush_group();

        if (max_assign > 0) {
            const long long total_mers = kept_mers + dropped_mers;
            std::cerr << "  mer load filter: max-assign=" << max_assign
                      << "  kept " << kept_mers << " mers/" << kept_rows << " rows"
                      << "  dropped " << dropped_mers << " mers/" << dropped_rows << " rows";
            if (total_mers > 0) {
                std::cerr << " (" << (100.0 * static_cast<double>(dropped_mers) /
                                      static_cast<double>(total_mers)) << "% mers)";
            }
            std::cerr << "  shared(>1 assign): " << shared_mers << "\n";
        }
    }

    int get_mer_size() const { return mer_size; }
    size_t size() const { return compressed_mer_query.size(); }
};


// ─── validate_read ───

// Reads-mode validation with optional T3 rescue verification. `verify` is
// invoked only when both strands fail the strict gates and T3 passes for a
// strand; it receives the read (forward orientation) and the strand to test.
template <typename Verify>
inline std::pair<int, int> validate_read_verified(
    const std::string& seq,
    const MerQueryManager& compressed_mer_query,
    int mer_size,
    int target_ele_id,
    const ValidatorConfig& cfg,
    Verify&& verify,
    double* rescue_score = nullptr)
{
    thread_local std::vector<int> match_ids;
    thread_local std::vector<int> match_starts;
    thread_local std::vector<int> match_loci;

    auto validate_chain = [&](const std::string& chain, int& r_ele, int& r_votes)
        -> std::pair<int, int> {
        match_ids.clear();
        match_starts.clear();
        match_loci.clear();

        int chain_len = static_cast<int>(chain.size());
        int mer_size_local = mer_size;

        if (chain_len < mer_size_local) {
            return {-1, 0};
        }

        // Rolling hash: O(1) per position
        static auto base_lookup = [] {
            std::array<int8_t, 256> t{};
            t.fill(-1);
            t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
            return t;
        }();

        match_ids.reserve(chain_len);
        match_starts.reserve(chain_len);
        match_loci.reserve(chain_len);

        uint32_t rolling = 0;
        uint32_t mer_mask = static_cast<uint32_t>((1ULL << (2 * mer_size_local)) - 1);
        int valid_run = 0;

        for (int i = 0; i < chain_len; ++i) {
            int8_t bits = base_lookup[static_cast<unsigned char>(chain[i])];
            if (bits < 0) {
                valid_run = 0;
                rolling = 0;
                continue;
            }
            rolling = ((rolling << 2) | static_cast<uint32_t>(bits)) & mer_mask;
            valid_run++;
            if (valid_run >= mer_size_local) {
                int start = i - mer_size_local + 1;
                auto it = compressed_mer_query.compressed_mer_query.find(rolling);
                if (it != compressed_mer_query.compressed_mer_query.end()) {
                    uint32_t off = MerQueryManager::range_off(it->second);
                    uint32_t cnt = MerQueryManager::range_cnt(it->second);
                    for (uint32_t a = 0; a < cnt; ++a) {
                        int id = compressed_mer_query.assign_arena[off + a].first;
                        int loci = compressed_mer_query.assign_arena[off + a].second;
                        if (target_ele_id == -1 || id == target_ele_id) {
                            match_ids.push_back(id);
                            match_starts.push_back(start);
                            match_loci.push_back(loci);
                        }
                    }
                }
            }
        }

        if (match_ids.empty()) return {-1, 0};

        // ─── Vote + gate matrix (T1 OR T2 ...) ───
        VoteTally tally = tally_votes(match_ids);

        auto chain_of = [&](int confi_id) -> std::pair<int, int> {
            int first = -1, last = -1, prev_x = -1, prev_l = -1, passing = 0, m = 0;
            for (size_t i = 0; i < match_ids.size(); ++i) {
                if (match_ids[i] != confi_id) continue;
                int x = match_starts[i];
                int l = match_loci[i];
                if (first < 0) first = x;
                else if (std::abs((x - prev_x) - (l - prev_l)) <= cfg.max_diff) ++passing;
                prev_x = x; prev_l = l; last = x; ++m;
            }
            if (m == 0) return {0, 0};
            return {passing, last - first + mer_size};
        };

        int passing = 0, span = 0;
        if (!gate_pass(tally, cfg, chain_of, passing, span)) {
            auto rg = rescue_gate(tally, cfg, chain_of);
            if (rg.first >= 0) {
                r_ele = rg.first;
                r_votes = rg.second;
            }
            return {-1, 0};
        }

        return {tally.best_id, static_cast<int>(tally.n)};
    };

    int f_ele = -1, f_votes = 0, r_ele = -1, r_votes = 0;
    auto [forward_id, forward_count] = validate_chain(seq, f_ele, f_votes);

    thread_local std::string _rc_buf;
    reverse_complement(seq, _rc_buf);
    auto [reverse_id, reverse_count] = validate_chain(_rc_buf, r_ele, r_votes);

    if (forward_id == -1 && reverse_id == -1) {
        if (cfg.verify.enabled) {
            int ele = -1, strand = 0;
            if (f_ele >= 0 && (r_ele < 0 || f_votes >= r_votes)) {
                ele = f_ele;
                strand = 1;
            } else if (r_ele >= 0) {
                ele = r_ele;
                strand = -1;
            }
            if (ele >= 0) {
                double sc = verify(ele, seq, strand);
                if (sc > 0.0) {
                    if (rescue_score) *rescue_score = sc;
                    return {ele, strand};
                }
            }
        }
        return {-1, 0};
    } else if (forward_id == -1) {
        return {reverse_id, -1};
    } else if (reverse_id == -1) {
        return {forward_id, 1};
    } else {
        return (forward_count >= reverse_count) ?
               std::make_pair(forward_id, 1) :
               std::make_pair(reverse_id, -1);
    }
}

// Strict-only wrapper (rescue verification disabled).
inline std::pair<int, int> validate_read(
    const std::string& seq,
    const MerQueryManager& compressed_mer_query,
    int mer_size,
    int target_ele_id = -1,
    const ValidatorConfig& cfg = ValidatorConfig{})
{
    return validate_read_verified(seq, compressed_mer_query, mer_size, target_ele_id, cfg,
                                  [](int, const std::string&, int) { return 0.0; });
}


// ─── Segmented sliding-window validation (genome mode) ───
//
// Validates every window [s, s+window_size) with s = 0, step_size, 2*step_size, ...
// exactly like validate_read(seq.substr(s, window_size)), but each k-mer position is
// looked up only once (its forward and reverse-complement ranges are cached), instead
// of re-hashing the 125/150 overlapping positions of neighbouring windows.
// `emit(start, strand, element_id)` is called for every accepted window.

inline uint32_t rc_code_2bit(uint32_t code, int k) {
    uint32_t res = 0;
    for (int i = 0; i < k; ++i) {
        uint32_t b = (code >> (2 * i)) & 0x3u;
        res = (res << 2) | (b ^ 0x3u);
    }
    return res;
}

template <typename Emit, typename Verify>
inline void validate_sequence_segmented(
    const std::string& seq,
    const MerQueryManager& mqm,
    int start_offset,
    const ValidatorConfig& cfg,
    Emit&& emit,
    Verify&& verify)
{
    const int k = mqm.mer_size;
    const int window_size = cfg.window_size;
    const int step_size = std::max(1, cfg.step_size);
    const int L = static_cast<int>(seq.size());
    if (k <= 0 || L < window_size) return;
    const int nspan = window_size - k + 1;   // k-mer start positions per window
    if (nspan <= 0) return;

    static const std::array<int8_t, 256> base_lookup = [] {
        std::array<int8_t, 256> t{};
        t.fill(-1);
        t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
        return t;
    }();
    const uint32_t mer_mask = static_cast<uint32_t>((1ULL << (2 * k)) - 1);

    const int SEG = 1 << 18;   // 256k bases per segment

    thread_local std::vector<uint64_t> fwd_ranges;
    thread_local std::vector<uint64_t> rc_ranges;
    thread_local std::vector<uint64_t> upper_bits;
    thread_local std::vector<uint32_t> codes_fwd;
    thread_local std::vector<uint32_t> codes_rc;
    // Per-segment match cache: each position's assignment range is expanded
    // once, and all overlapping windows read contiguous slices of it.
    thread_local std::vector<uint32_t> off_fwd, off_rc;    // position -> offset
    thread_local std::vector<uint64_t> pack_fwd, pack_rc;  // (id << 32) | locus
    int cur_o = 0;   // current window origin within the segment
    const bool vp = validate_profile_enabled();
    thread_local long long tl_prep_ns = 0;
    thread_local long long tl_win_ns = 0;
    thread_local long long tl_gather_ns = 0;
    thread_local long long tl_tally_ns = 0;

    // Tally updated directly from the range arrays, without materializing the
    // match list. The list is only built for windows that pass the ratio
    // pre-check and therefore need chain evaluation.
    struct TallyAcc {
        std::vector<uint64_t> sc;   // (stamp << 32) | count
        std::vector<int> seen;
        unsigned epoch = 0;
        int total = 0;
        void begin() {
            ++epoch;
            if (epoch == 0) {
                std::fill(sc.begin(), sc.end(), 0ull);
                epoch = 1;
            }
            seen.clear();
            total = 0;
        }
        inline void add(int id) {
            ++total;
            if (id >= static_cast<int>(sc.size())) sc.resize(id + 1, 0ull);
            uint64_t v = sc[id];
            if (static_cast<unsigned>(v >> 32) != epoch) {
                sc[id] = (static_cast<uint64_t>(epoch) << 32) | 1ull;
                seen.push_back(id);
            } else {
                sc[id] = v + 1;
            }
        }
        VoteTally finish() const {
            VoteTally r;
            r.n = total;
            for (int id : seen) {
                int c = static_cast<int>(sc[id] & 0xffffffffu);
                if (c > r.best_cnt) {
                    r.second_cnt = r.best_cnt;
                    r.best_cnt = c;
                    r.best_id = id;
                } else if (c > r.second_cnt) {
                    r.second_cnt = c;
                }
            }
            return r;
        }
    };
    thread_local TallyAcc tally_acc;

    auto make_chain = [&](const std::vector<uint32_t>* off,
                          const std::vector<uint64_t>* pack, bool rc_strand) {
        return [&, off, pack, rc_strand](int confi_id) -> std::pair<int, int> {
            if (vp) ++g_val_profile.chain_calls;
            int first = -1, last = -1, prev_x = -1, prev_l = -1, passing = 0, m = 0;
            for (int x = 0; x < nspan; ++x) {
                // The RC strand enumerates k-mer starts backwards (its first
                // match is the window's last position), so map x to the cache
                // position in that order.
                int p = rc_strand ? (cur_o + (nspan - 1 - x)) : (cur_o + x);
                uint32_t b = (*off)[p];
                uint32_t e = (*off)[p + 1];
                for (uint32_t j = b; j < e; ++j) {
                    uint64_t v = (*pack)[j];
                    int id = static_cast<int>(v >> 32);
                    if (id != confi_id) continue;
                    int l = static_cast<int>(static_cast<uint32_t>(v));
                    if (first < 0) first = x;
                    else if (std::abs((x - prev_x) - (l - prev_l)) <= cfg.max_diff) ++passing;
                    prev_x = x; prev_l = l; last = x; ++m;
                }
            }
            if (m == 0) return {0, 0};
            return {passing, last - first + k};
        };
    };
    auto chain_fwd = make_chain(&off_fwd, &pack_fwd, false);
    auto chain_rc = make_chain(&off_rc, &pack_rc, true);

    const int seg_stride = SEG - (window_size - 1);
    for (int seg_start = 0; seg_start < L; seg_start += seg_stride) {
        const int seg_end = std::min(L, seg_start + SEG);
        const int seg_len = seg_end - seg_start;

        fwd_ranges.assign(seg_len, 0);
        rc_ranges.assign(seg_len, 0);
        upper_bits.assign((seg_len + 63) >> 6, 0ULL);

        auto _t_prep = std::chrono::steady_clock::now();
        if (mqm.has_direct_index()) {
            // Two passes: compute k-mer codes sequentially (no dependent
            // memory stalls), then run the direct-index lookups with software
            // prefetching to hide the random-access latency.
            codes_fwd.assign(seg_len, UINT32_MAX);
            codes_rc.assign(seg_len, UINT32_MAX);
            uint32_t rolling = 0;
            uint32_t rc_rolling = 0;
            int valid_run = 0;
            for (int i = seg_start; i < seg_end; ++i) {
                unsigned char raw = static_cast<unsigned char>(seq[i]);
                char c = static_cast<char>(raw);
                if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
                if (raw >= 'A' && raw <= 'Z')
                    upper_bits[(i - seg_start) >> 6] |= (1ULL << ((i - seg_start) & 63));
                int8_t bits = base_lookup[static_cast<unsigned char>(c)];
                if (bits < 0) {
                    valid_run = 0;
                    rolling = 0;
                    rc_rolling = 0;
                    continue;
                }
                rolling = ((rolling << 2) | static_cast<uint32_t>(bits)) & mer_mask;
                rc_rolling = (rc_rolling >> 2) |
                             (static_cast<uint32_t>(3 - bits) << (2 * k - 2));
                ++valid_run;
                if (valid_run >= k) {
                    int st = i - k + 1 - seg_start;
                    codes_fwd[st] = rolling;
                    codes_rc[st] = rc_rolling;
                }
            }
            constexpr int kPrefetchDist = 32;
            for (int st = 0; st < seg_len; ++st) {
                int pf = st + kPrefetchDist;
                if (pf < seg_len && codes_fwd[pf] != UINT32_MAX) {
                    __builtin_prefetch(&mqm.direct_index[codes_fwd[pf]]);
                    __builtin_prefetch(&mqm.direct_index[codes_rc[pf]]);
                }
                uint32_t cf = codes_fwd[st];
                if (cf == UINT32_MAX) continue;
                uint32_t r = mqm.direct_index[cf];
                if (r) fwd_ranges[st] = mqm.direct_rng[r - 1];
                uint32_t rr = mqm.direct_index[codes_rc[st]];
                if (rr) rc_ranges[st] = mqm.direct_rng[rr - 1];
            }
        } else {
            uint32_t rolling = 0;
            int valid_run = 0;
            for (int i = seg_start; i < seg_end; ++i) {
                unsigned char raw = static_cast<unsigned char>(seq[i]);
                char c = static_cast<char>(raw);
                if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
                if (raw >= 'A' && raw <= 'Z')
                    upper_bits[(i - seg_start) >> 6] |= (1ULL << ((i - seg_start) & 63));
                int8_t bits = base_lookup[static_cast<unsigned char>(c)];
                if (bits < 0) {
                    valid_run = 0;
                    rolling = 0;
                    continue;
                }
                rolling = ((rolling << 2) | static_cast<uint32_t>(bits)) & mer_mask;
                ++valid_run;
                if (valid_run >= k) {
                    int st = i - k + 1 - seg_start;   // k-mer start index within the segment
                    auto it = mqm.compressed_mer_query.find(rolling);
                    if (it != mqm.compressed_mer_query.end()) fwd_ranges[st] = it->second;
                    uint32_t rc = rc_code_2bit(rolling, k);
                    auto it2 = mqm.compressed_mer_query.find(rc);
                    if (it2 != mqm.compressed_mer_query.end()) rc_ranges[st] = it2->second;
                }
            }
        }
        if (vp) {
            tl_prep_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - _t_prep).count();
        }

        // Expand every position's assignment range once for the whole segment;
        // all overlapping windows then share contiguous slices of this cache.
        off_fwd.assign(seg_len + 1, 0);
        off_rc.assign(seg_len + 1, 0);
        pack_fwd.clear();
        pack_rc.clear();
        {
            auto _t_cache = vp ? std::chrono::steady_clock::now()
                               : std::chrono::steady_clock::time_point{};
            for (int p = 0; p < seg_len; ++p) {
                constexpr int kCachePrefetch = 32;
                if (p + kCachePrefetch < seg_len) {
                    uint64_t rp = fwd_ranges[p + kCachePrefetch];
                    if (rp)
                        __builtin_prefetch(&mqm.assign_arena[MerQueryManager::range_off(rp)]);
                }
                uint64_t rf = fwd_ranges[p];
                if (rf) {
                    uint32_t off = MerQueryManager::range_off(rf);
                    uint32_t cnt = MerQueryManager::range_cnt(rf);
                    for (uint32_t a = 0; a < cnt; ++a) {
                        const auto& pr = mqm.assign_arena[off + a];
                        pack_fwd.push_back(
                            (static_cast<uint64_t>(static_cast<uint32_t>(pr.first)) << 32) |
                            static_cast<uint32_t>(pr.second));
                    }
                }
                off_fwd[p + 1] = static_cast<uint32_t>(pack_fwd.size());
            }
            for (int p = 0; p < seg_len; ++p) {
                constexpr int kCachePrefetch = 32;
                if (p + kCachePrefetch < seg_len) {
                    uint64_t rp = rc_ranges[p + kCachePrefetch];
                    if (rp)
                        __builtin_prefetch(&mqm.assign_arena[MerQueryManager::range_off(rp)]);
                }
                uint64_t rr = rc_ranges[p];
                if (rr) {
                    uint32_t off = MerQueryManager::range_off(rr);
                    uint32_t cnt = MerQueryManager::range_cnt(rr);
                    for (uint32_t a = 0; a < cnt; ++a) {
                        const auto& pr = mqm.assign_arena[off + a];
                        pack_rc.push_back(
                            (static_cast<uint64_t>(static_cast<uint32_t>(pr.first)) << 32) |
                            static_cast<uint32_t>(pr.second));
                    }
                }
                off_rc[p + 1] = static_cast<uint32_t>(pack_rc.size());
            }
            if (vp) {
                tl_gather_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - _t_cache).count();
                g_val_profile.pushed += (long long)(pack_fwd.size() + pack_rc.size());
            }
        }
        auto _t_win = std::chrono::steady_clock::now();

        int s0 = seg_start;
        if (step_size > 1) {
            int rem = (s0 + start_offset) % step_size;
            if (rem != 0) s0 += step_size - rem;
        }
        const int last_start = seg_end - window_size;
        for (int s = s0; s <= last_start; s += step_size) {
            if (vp) ++g_val_profile.windows;
            // soft-mask skip: window must contain at least one uppercase base
            int o = s - seg_start;
            int e = o + window_size;
            bool has_upper = false;
            for (int w = o >> 6; w <= (e - 1) >> 6; ++w) {
                uint64_t word = upper_bits[w];
                if (w == (o >> 6)) word &= (~0ULL << (o & 63));
                if (w == ((e - 1) >> 6)) {
                    int hi = (e - 1) & 63;
                    word &= (hi == 63) ? ~0ULL : ((1ULL << (hi + 1)) - 1);
                }
                if (word) { has_upper = true; break; }
            }
            if (!has_upper) continue;

            // forward strand: tally a contiguous slice of the segment cache;
            // chain evaluation runs only when the ratio pre-check can pass.
            cur_o = o;
            if (vp) ++g_val_profile.tally_calls;
            auto _tt = vp ? std::chrono::steady_clock::now()
                          : std::chrono::steady_clock::time_point{};
            tally_acc.begin();
            for (uint32_t j = off_fwd[o]; j < off_fwd[o + nspan]; ++j) {
                tally_acc.add(static_cast<int>(pack_fwd[j] >> 32));
            }
            if (vp) {
                tl_tally_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - _tt).count();
            }
            int fwd_id = -1, fwd_count = 0;
            bool fwd_ok = false;
            VoteTally fwd_t = tally_acc.finish();
            if (gate_ratio_precheck(fwd_t, cfg)) {
                int passing = 0, span = 0;
                if (gate_pass(fwd_t, cfg, chain_fwd, passing, span)) {
                    fwd_ok = true;
                    fwd_id = fwd_t.best_id;
                    fwd_count = fwd_t.n;
                }
            }

            // reverse-complement strand: same slice, but matches are enumerated
            // backwards (window's last k-mer position first) to preserve the
            // legacy match order for tally tie-breaking.
            if (vp) ++g_val_profile.tally_calls;
            _tt = vp ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
            tally_acc.begin();
            for (int p = o + nspan - 1; p >= o; --p) {
                for (uint32_t j = off_rc[p]; j < off_rc[p + 1]; ++j) {
                    tally_acc.add(static_cast<int>(pack_rc[j] >> 32));
                }
            }
            if (vp) {
                tl_tally_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - _tt).count();
            }
            int rev_id = -1, rev_count = 0;
            bool rev_ok = false;
            VoteTally rev_t = tally_acc.finish();
            if (gate_ratio_precheck(rev_t, cfg)) {
                int passing = 0, span = 0;
                if (gate_pass(rev_t, cfg, chain_rc, passing, span)) {
                    rev_ok = true;
                    rev_id = rev_t.best_id;
                    rev_count = rev_t.n;
                }
            }
            if (vp) {
                tl_tally_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - _tt).count();
            }

            if (fwd_ok || rev_ok) {
                if (vp) ++g_val_profile.emits;
                if (!fwd_ok) emit(s, -1, rev_id, 0.0);
                else if (!rev_ok) emit(s, 1, fwd_id, 0.0);
                else if (fwd_count >= rev_count) emit(s, 1, fwd_id, 0.0);
                else emit(s, -1, rev_id, 0.0);
            } else if (cfg.verify.enabled) {
                // T3 near-miss: SIMD SW verification decides (both strands).
                auto rf = rescue_gate(fwd_t, cfg, chain_fwd);
                auto rr = rescue_gate(rev_t, cfg, chain_rc);
                int ele = -1, strand = 0;
                if (rf.first >= 0 && (rr.first < 0 || rf.second >= rr.second)) {
                    ele = rf.first;
                    strand = 1;
                } else if (rr.first >= 0) {
                    ele = rr.first;
                    strand = -1;
                }
                if (ele >= 0) {
                    std::string w(seq, s, window_size);
                    double sc = verify(ele, w, strand);
                    if (sc > 0.0) {
                        if (vp) ++g_val_profile.emits;
                        emit(s, strand, ele, sc);
                    }
                }
            }
        }
        if (vp) {
            tl_win_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - _t_win).count();
        }
    }
    if (vp) {
        g_val_profile.prep_ns += tl_prep_ns;
        g_val_profile.win_ns += tl_win_ns;
        g_val_profile.gather_ns += tl_gather_ns;
        g_val_profile.tally_ns += tl_tally_ns;
        tl_prep_ns = 0;
        tl_win_ns = 0;
        tl_gather_ns = 0;
        tl_tally_ns = 0;
    }
}
