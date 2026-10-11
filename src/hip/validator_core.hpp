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
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

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
inline VoteTally tally_votes(const std::vector<int>& nums) {
    VoteTally r;
    r.n = static_cast<int>(nums.size());
    if (nums.empty()) return r;

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

    for (int id : nums) {
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

// Gate matrix: T1 = legacy defaults (from the flags); T2, T3, ... are extra
// OR-combined rows. A window is accepted when ANY row passes. Rows share one
// tally and one chain evaluation per window.
struct GateTier {
    double vote_ratio = 3.0;
    int min_c = 7;
    int min_span = 25;
};

struct ValidatorConfig {
    int window_size = 150;
    int step_size = 25;
    double vote_frac = 0.1;
    int max_diff = 2;
    std::vector<GateTier> tiers = {GateTier{}};   // T1 OR T2 ...
};

// chain_fn(best_id) -> (passing, span); invoked at most once per window.
template <typename ChainFn>
inline bool gate_pass(const VoteTally& t, const ValidatorConfig& cfg, ChainFn&& chain_fn,
                      int& passing_out, int& span_out) {
    passing_out = 0;
    span_out = 0;
    if (t.best_id < 0 || cfg.tiers.empty()) return false;

    float threshold = static_cast<float>(t.n) * static_cast<float>(cfg.vote_frac);
    if (t.best_cnt <= threshold) return false;

    // Integer truncation matches the historical T1 semantics; keep it even
    // for fractional ratios (e.g. 1.5) to avoid changing accepted windows.
    auto ratio_ok = [&](const GateTier& tier) {
        return t.best_cnt >= static_cast<int>(tier.vote_ratio * t.second_cnt);
    };

    bool any_ratio = false;
    for (const auto& tier : cfg.tiers) {
        if (ratio_ok(tier)) {
            any_ratio = true;
            break;
        }
    }
    if (!any_ratio) return false;

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

    void load_from_file(const std::string& path) {
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
        size_t cur_off = 0;
        uint32_t cur_cnt = 0;
        bool have_cur = false;
        auto flush_range = [&]() {
            if (have_cur && cur_cnt > 0) {
                compressed_mer_query[cur_key] = pack_range(cur_off, cur_cnt);
            }
        };

        while (std::getline(file, line)) {
            if (line.empty()) continue;
            size_t t1 = line.find('\t');
            if (t1 == std::string::npos) continue;
            size_t t2 = line.find('\t', t1 + 1);
            if (t2 == std::string::npos) continue;
            size_t t3 = line.find('\t', t2 + 1);
            std::string loc_s = (t3 == std::string::npos)
                                ? line.substr(t2 + 1)
                                : line.substr(t2 + 1, t3 - t2 - 1);
            int id = std::stoi(line.substr(t1 + 1, t2 - t1 - 1));
            int loci = std::stoi(loc_s);
            if (mer_size == 0) mer_size = static_cast<int>(t1);
            uint32_t encoded_mer;
            try {
                encoded_mer = dna_encoder(line.substr(0, t1));
            } catch (const std::invalid_argument&) {
                continue;
            }
            if (!have_cur || encoded_mer != cur_key) {
                flush_range();
                cur_key = encoded_mer;
                cur_off = assign_arena.size();
                cur_cnt = 0;
                have_cur = true;
            }
            assign_arena.push_back({id, loci});
            ++cur_cnt;
        }
        flush_range();
    }

    int get_mer_size() const { return mer_size; }
    size_t size() const { return compressed_mer_query.size(); }
};


// ─── validate_read ───

inline std::pair<int, int> validate_read(
    const std::string& seq,
    const MerQueryManager& compressed_mer_query,
    int mer_size,
    int target_ele_id = -1,
    const ValidatorConfig& cfg = ValidatorConfig{})
{
    thread_local std::vector<int> match_ids;
    thread_local std::vector<int> match_starts;
    thread_local std::vector<int> match_loci;

    auto validate_chain = [&](const std::string& chain) -> std::pair<int, int> {
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
        if (!gate_pass(tally, cfg, chain_of, passing, span)) return {-1, 0};

        return {tally.best_id, static_cast<int>(tally.n)};
    };

    auto [forward_id, forward_count] = validate_chain(seq);

    thread_local std::string _rc_buf;
    reverse_complement(seq, _rc_buf);
    auto [reverse_id, reverse_count] = validate_chain(_rc_buf);

    if (forward_id == -1 && reverse_id == -1) {
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

template <typename Emit>
inline void validate_sequence_segmented(
    const std::string& seq,
    const MerQueryManager& mqm,
    int start_offset,
    const ValidatorConfig& cfg,
    Emit&& emit)
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
    thread_local std::vector<int> match_ids;
    thread_local std::vector<int> match_starts;
    thread_local std::vector<int> match_loci;

    auto chain_stats = [&](int confi_id) -> std::pair<int, int> {
        int first = -1, last = -1, prev_x = -1, prev_l = -1, passing = 0, m = 0;
        for (size_t idx = 0; idx < match_ids.size(); ++idx) {
            if (match_ids[idx] != confi_id) continue;
            int x = match_starts[idx];
            int l = match_loci[idx];
            if (first < 0) first = x;
            else if (std::abs((x - prev_x) - (l - prev_l)) <= cfg.max_diff) ++passing;
            prev_x = x; prev_l = l; last = x; ++m;
        }
        if (m == 0) return {0, 0};
        return {passing, last - first + k};
    };

    const int seg_stride = SEG - (window_size - 1);
    for (int seg_start = 0; seg_start < L; seg_start += seg_stride) {
        const int seg_end = std::min(L, seg_start + SEG);
        const int seg_len = seg_end - seg_start;

        fwd_ranges.assign(seg_len, 0);
        rc_ranges.assign(seg_len, 0);
        upper_bits.assign((seg_len + 63) >> 6, 0ULL);

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

        int s0 = seg_start;
        if (step_size > 1) {
            int rem = (s0 + start_offset) % step_size;
            if (rem != 0) s0 += step_size - rem;
        }
        const int last_start = seg_end - window_size;
        for (int s = s0; s <= last_start; s += step_size) {
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

            // forward chain: k-mer starts increasing
            match_ids.clear(); match_starts.clear(); match_loci.clear();
            for (int i = s; i < s + nspan; ++i) {
                uint64_t r = fwd_ranges[i - seg_start];
                if (!r) continue;
                uint32_t off = MerQueryManager::range_off(r);
                uint32_t cnt = MerQueryManager::range_cnt(r);
                int x = i - s;
                for (uint32_t a = 0; a < cnt; ++a) {
                    match_ids.push_back(mqm.assign_arena[off + a].first);
                    match_starts.push_back(x);
                    match_loci.push_back(mqm.assign_arena[off + a].second);
                }
            }
            int fwd_id = -1, fwd_count = 0;
            bool fwd_ok = false;
            if (!match_ids.empty()) {
                VoteTally t = tally_votes(match_ids);
                int passing = 0, span = 0;
                if (gate_pass(t, cfg, chain_stats, passing, span)) {
                    fwd_ok = true;
                    fwd_id = t.best_id;
                    fwd_count = t.n;
                }
            }

            // reverse-complement chain: k-mer starts decreasing => RC positions increasing
            match_ids.clear(); match_starts.clear(); match_loci.clear();
            for (int i = s + nspan - 1; i >= s; --i) {
                uint64_t r = rc_ranges[i - seg_start];
                if (!r) continue;
                uint32_t off = MerQueryManager::range_off(r);
                uint32_t cnt = MerQueryManager::range_cnt(r);
                int x = (nspan - 1) - (i - s);
                for (uint32_t a = 0; a < cnt; ++a) {
                    match_ids.push_back(mqm.assign_arena[off + a].first);
                    match_starts.push_back(x);
                    match_loci.push_back(mqm.assign_arena[off + a].second);
                }
            }
            int rev_id = -1, rev_count = 0;
            bool rev_ok = false;
            if (!match_ids.empty()) {
                VoteTally t = tally_votes(match_ids);
                int passing = 0, span = 0;
                if (gate_pass(t, cfg, chain_stats, passing, span)) {
                    rev_ok = true;
                    rev_id = t.best_id;
                    rev_count = t.n;
                }
            }

            if (!fwd_ok && !rev_ok) continue;
            if (!fwd_ok) emit(s, -1, rev_id);
            else if (!rev_ok) emit(s, 1, fwd_id);
            else if (fwd_count >= rev_count) emit(s, 1, fwd_id);
            else emit(s, -1, rev_id);
        }
    }
}
