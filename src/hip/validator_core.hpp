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
#include <cerrno>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// Optional counters for the segmented genome-mode validator. Enabled by
// setting the VALIDATE_PROFILE environment variable to a non-zero value.
struct ValidateProfile {
    std::atomic<long long> windows{0};
    std::atomic<long long> emits{0};
    std::atomic<long long> tally_calls{0};
    std::atomic<long long> chain_calls{0};
    std::atomic<long long> prep_ns{0};
    std::atomic<long long> win_ns{0};
    std::atomic<long long> tally_ns{0};
    // Reads mode breakdown.
    std::atomic<long long> rd_reads{0};
    std::atomic<long long> rd_scan_ns{0};
    std::atomic<long long> rd_tally_ns{0};
    std::atomic<long long> rd_gate_ns{0};
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

// Reverse complement of a 2-bit-encoded k-mer (A=00,C=01,G=10,T=11).
inline uint32_t rc_code_2bit(uint32_t code, int k) {
    uint32_t res = 0;
    for (int i = 0; i < k; ++i) {
        uint32_t b = (code >> (2 * i)) & 0x3u;
        res = (res << 2) | (b ^ 0x3u);
    }
    return res;
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
    std::vector<uint32_t> present_codes;
    size_t kept_mers_ = 0;
    int mer_size = 0;

    // Optional direct-index accelerator (k <= 13): maps each 2-bit k-mer code
    // to the packed assignment range (count in the top 3 bits, arena offset in
    // the low 29). Storing the range inline lets a scan resolve code -> range in
    // ONE random access instead of the previous code -> rank -> range chain
    // (which needed a second direct_rng array). Costs 4^k * 4 bytes.
    std::vector<uint32_t> direct_index;

    static constexpr int kDirCntShift = 29;                 // 3 count bits, 29 offset bits
    static constexpr uint32_t kDirOffMask = (1u << kDirCntShift) - 1;
    static inline uint32_t pack_dir(uint32_t off, uint32_t cnt) { return (cnt << kDirCntShift) | off; }
    static inline uint32_t dir_off(uint32_t v) { return v & kDirOffMask; }
    static inline uint32_t dir_cnt(uint32_t v) { return v >> kDirCntShift; }

    bool has_direct_index() const { return !direct_index.empty(); }

    void build_direct_index() {
        if (mer_size <= 0 || mer_size > 13) return;
        if (!direct_index.empty()) return;
        direct_index.assign(size_t(1) << (2 * mer_size), 0);
        for (const auto& kv : compressed_mer_query) {
            uint32_t off = range_off(kv.second);
            uint32_t cnt = range_cnt(kv.second);
            if (off > kDirOffMask || cnt > 7) { direct_index.clear(); return; }
            direct_index[kv.first] = pack_dir(off, cnt);
        }
    }

    void ensure_hashmap() {
        if (!compressed_mer_query.empty() || kept_mers_ == 0) return;
        compressed_mer_query.reserve(kept_mers_);
        if (!direct_index.empty() && !present_codes.empty()) {
            for (uint32_t code : present_codes) {
                uint32_t v = direct_index[code];
                if (v) compressed_mer_query[code] = pack_range(dir_off(v), dir_cnt(v));
            }
        }
    }

    void drop_raw_after_canon() {
        std::vector<uint32_t>().swap(direct_index);
        std::vector<uint32_t>().swap(present_codes);
        std::vector<std::pair<int, int>>().swap(assign_arena);
        compressed_mer_query.clear();
    }

    inline bool has_id_fast(uint32_t code, int id) const {
        if (!direct_index.empty()) {
            uint32_t v = direct_index[code];
            if (!v) return false;
            uint32_t off = dir_off(v);
            uint32_t cnt = dir_cnt(v);
            for (uint32_t i = 0; i < cnt; ++i) {
                if (assign_arena[off + i].first == id) return true;
            }
            return false;
        }
        auto it = compressed_mer_query.find(code);
        return it != compressed_mer_query.end() && has_id(it->second, id);
    }

    // ─── Canonical (strand-merged) index ───
    // Keyed by min(code, rc(code)); a canonical key aggregates the assignments
    // of both a k-mer and its reverse complement. Each row carries a form bit
    // (0 = this assignment came from the canonical-form k-mer, 1 = from its
    // reverse complement), so a single forward scan can resolve both strands.
    // Row layout: form in bit 31, element id in bits 16..30, locus in bits 0..15.
    std::vector<uint32_t> dir_canon;
    std::vector<uint32_t> arena_canon;   // packed form|id|loci
    bool has_canon_index() const { return !dir_canon.empty(); }
    static inline uint32_t pack_crow(int id, int loci, uint32_t form) {
        return (form << 31) | (static_cast<uint32_t>(id) << 16) | (static_cast<uint32_t>(loci) & 0xffffu);
    }
    static inline int crow_id(uint32_t v) { return static_cast<int>((v >> 16) & 0x7fffu); }
    static inline int crow_loci(uint32_t v) { return static_cast<int>(v & 0xffffu); }
    static inline uint32_t crow_form(uint32_t v) { return v >> 31; }

    // Build the canonical index from the loaded (raw) table. Returns false (and
    // leaves dir_canon empty) if the id/loci ranges do not fit the row packing;
    // callers then fall back to the raw direct index.
    bool build_canonical_index() {
        if (mer_size <= 0 || mer_size > 13) return false;
        const int k = mer_size;

        // Zero-hash array path when direct_index + present_codes are already populated.
        if (!direct_index.empty() && !present_codes.empty()) {
            dir_canon.assign(size_t(1) << (2 * k), 0);
            std::vector<uint32_t> canon_keys;
            canon_keys.reserve(present_codes.size());
            size_t total = 0;
            for (uint32_t code : present_codes) {
                uint32_t v = direct_index[code];
                uint32_t n = dir_cnt(v);
                if (n == 0) continue;
                uint32_t c = std::min(code, rc_code_2bit(code, k));
                if (dir_canon[c] == 0) canon_keys.push_back(c);
                dir_canon[c] += n;
                total += n;
            }
            if (total > kDirOffMask) { dir_canon.clear(); return false; }
            uint32_t off = 0;
            for (uint32_t c : canon_keys) {
                uint32_t cnt = dir_canon[c];
                if (cnt > 7) { dir_canon.clear(); return false; }
                dir_canon[c] = pack_dir(off, cnt);
                off += cnt;
            }
            arena_canon.resize(total);
            for (uint32_t code : present_codes) {
                uint32_t v = direct_index[code];
                uint32_t o = dir_off(v), n = dir_cnt(v);
                if (n == 0) continue;
                uint32_t c = std::min(code, rc_code_2bit(code, k));
                uint32_t form = (code == c) ? 0u : 1u;
                uint32_t cur_off = dir_off(dir_canon[c]);
                for (uint32_t a = 0; a < n; ++a) {
                    int id = assign_arena[o + a].first;
                    int loci = assign_arena[o + a].second;
                    if (id < 0 || id > 0x7fff || loci < 0 || loci > 0xffff) {
                        dir_canon.clear();
                        arena_canon.clear();
                        return false;
                    }
                    arena_canon[cur_off + a] = pack_crow(id, loci, form);
                }
                dir_canon[c] += n;
            }
            for (uint32_t c : canon_keys) {
                dir_canon[c] -= dir_cnt(dir_canon[c]);
            }
            return true;
        }

        tsl::robin_map<uint32_t, uint32_t> cnt;   // canonical -> row count
        cnt.reserve(compressed_mer_query.size());
        size_t total = 0;
        for (const auto& kv : compressed_mer_query) {
            uint32_t code = kv.first;
            uint32_t c = std::min(code, rc_code_2bit(code, k));
            uint32_t n = range_cnt(kv.second);
            if (range_off(kv.second) + n > 0xFFFFFFFFu) return false;
            cnt[c] += n;
            total += n;
        }
        if (total > kDirOffMask) return false;
        dir_canon.assign(size_t(1) << (2 * k), 0);
        tsl::robin_map<uint32_t, uint32_t> cursor;
        cursor.reserve(cnt.size());
        uint32_t off = 0;
        bool ok = true;
        for (const auto& kv : cnt) {
            if (kv.second > 7) { ok = false; break; }
            dir_canon[kv.first] = pack_dir(off, kv.second);
            cursor[kv.first] = off;
            off += kv.second;
        }
        cnt.clear();
        if (!ok) { dir_canon.clear(); return false; }
        arena_canon.resize(total);
        for (const auto& kv : compressed_mer_query) {
            uint32_t code = kv.first;
            uint32_t c = std::min(code, rc_code_2bit(code, k));
            uint32_t form = (code == c) ? 0u : 1u;
            uint32_t o = range_off(kv.second), n = range_cnt(kv.second);
            for (uint32_t a = 0; a < n; ++a) {
                int id = assign_arena[o + a].first;
                int loci = assign_arena[o + a].second;
                if (id > 0x7fff || loci > 0xffff) { dir_canon.clear(); arena_canon.clear(); return false; }
                arena_canon[cursor[c]++] = pack_crow(id, loci, form);
            }
        }
        return true;
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
    // When defer_hashmap is true and the table fits direct_index (k<=13, cnt<=7),
    // skip populating the 12M-entry robin_map and build direct_index in place.
    void load_from_file(const std::string& path, int max_assign = 0, bool defer_hashmap = false) {
        struct BinHeader {
            char magic[8];
            uint32_t version;
            int32_t mer_size;
            int32_t max_assign;
            uint32_t fits_direct;
            int64_t src_size;
            int64_t src_mtime_sec;
            int64_t src_mtime_nsec;
            uint64_t kept_mers;
            uint64_t arena_size;
            uint64_t dropped_mers;
            uint64_t dropped_rows;
            uint64_t shared_mers;
        };

        auto read_exact = [](int fd, void* dst, size_t nbytes) -> bool {
            char* p = static_cast<char*>(dst);
            size_t got = 0;
            while (got < nbytes) {
                ssize_t r = ::read(fd, p + got, nbytes - got);
                if (r < 0) {
                    if (errno == EINTR) continue;
                    return false;
                }
                if (r == 0) return false;
                got += static_cast<size_t>(r);
            }
            return true;
        };

        auto write_exact = [](int fd, const void* src, size_t nbytes) -> bool {
            const char* p = static_cast<const char*>(src);
            size_t sent = 0;
            while (sent < nbytes) {
                ssize_t w = ::write(fd, p + sent, nbytes - sent);
                if (w < 0) {
                    if (errno == EINTR) continue;
                    return false;
                }
                sent += static_cast<size_t>(w);
            }
            return true;
        };

        auto print_filter_stats = [&](long long km, long long kr, long long dm,
                                      long long dr, long long sm) {
            if (max_assign > 0) {
                const long long total_mers = km + dm;
                std::cerr << "  mer load filter: max-assign=" << max_assign
                          << "  kept " << km << " mers/" << kr << " rows"
                          << "  dropped " << dm << " mers/" << dr << " rows";
                if (total_mers > 0) {
                    std::cerr << " (" << (100.0 * static_cast<double>(dm) /
                                          static_cast<double>(total_mers)) << "% mers)";
                }
                std::cerr << "  shared(>1 assign): " << sm << "\n";
            }
        };

        struct stat src_st{};
        if (::stat(path.c_str(), &src_st) != 0) {
            throw std::runtime_error("Cannot open mers file: " + path);
        }

        const char* no_cache_env = std::getenv("CNEX_NO_BIN_CACHE");
        const bool use_cache = (no_cache_env == nullptr || no_cache_env[0] == '\0' || no_cache_env[0] == '0');
        const std::string cache_path = path + ".ma" + std::to_string(max_assign) + ".cnexidx";

        if (use_cache) {
            struct stat cst{};
            if (::stat(cache_path.c_str(), &cst) == 0 &&
                cst.st_size >= static_cast<off_t>(sizeof(BinHeader))) {
                int cfd = ::open(cache_path.c_str(), O_RDONLY);
                if (cfd >= 0) {
                    BinHeader hdr{};
                    bool ok = read_exact(cfd, &hdr, sizeof(hdr)) &&
                              std::memcmp(hdr.magic, "CNEXIDX1", 8) == 0 &&
                              hdr.version == 1 &&
                              hdr.max_assign == max_assign &&
                              hdr.src_size == static_cast<int64_t>(src_st.st_size) &&
                              hdr.src_mtime_sec == static_cast<int64_t>(src_st.st_mtim.tv_sec) &&
                              hdr.src_mtime_nsec == static_cast<int64_t>(src_st.st_mtim.tv_nsec);
                    uint64_t expected_sz = sizeof(BinHeader) +
                        hdr.kept_mers * sizeof(uint32_t) +
                        hdr.kept_mers * sizeof(uint64_t) +
                        hdr.arena_size * sizeof(std::pair<int, int>);
                    if (ok && static_cast<uint64_t>(cst.st_size) == expected_sz) {
                        present_codes.resize(hdr.kept_mers);
                        std::vector<uint64_t> packed_ranges(hdr.kept_mers);
                        assign_arena.resize(hdr.arena_size);
                        bool body_ok =
                            read_exact(cfd, present_codes.data(), hdr.kept_mers * sizeof(uint32_t)) &&
                            read_exact(cfd, packed_ranges.data(), hdr.kept_mers * sizeof(uint64_t)) &&
                            read_exact(cfd, assign_arena.data(), hdr.arena_size * sizeof(std::pair<int, int>));
                        ::close(cfd);
                        if (body_ok) {
                            mer_size = hdr.mer_size;
                            kept_mers_ = hdr.kept_mers;
                            if (defer_hashmap && hdr.fits_direct && mer_size >= 7 && mer_size <= 13) {
                                direct_index.assign(size_t(1) << (2 * mer_size), 0);
                                for (size_t i = 0; i < hdr.kept_mers; ++i) {
                                    uint64_t r = packed_ranges[i];
                                    direct_index[present_codes[i]] = pack_dir(range_off(r), range_cnt(r));
                                }
                            } else {
                                compressed_mer_query.reserve(hdr.kept_mers);
                                for (size_t i = 0; i < hdr.kept_mers; ++i) {
                                    compressed_mer_query[present_codes[i]] = packed_ranges[i];
                                }
                            }
                            print_filter_stats(static_cast<long long>(hdr.kept_mers),
                                               static_cast<long long>(hdr.arena_size),
                                               static_cast<long long>(hdr.dropped_mers),
                                               static_cast<long long>(hdr.dropped_rows),
                                               static_cast<long long>(hdr.shared_mers));
                            return;
                        }
                        present_codes.clear();
                        assign_arena.clear();
                    } else {
                        ::close(cfd);
                    }
                }
            }
        }

        const size_t file_size = static_cast<size_t>(src_st.st_size);
        assign_arena.reserve(file_size / 20);
        present_codes.reserve(file_size / 25);

        uint32_t cur_key = 0;
        bool have_cur = false;
        std::vector<std::pair<int, int>> group;
        group.reserve(64);
        long long kept_mers = 0, kept_rows = 0, dropped_mers = 0, dropped_rows = 0;
        long long shared_mers = 0;
        bool dir_inited = false;
        bool can_direct = false;

        auto switch_to_hashmap = [&]() {
            if (!can_direct) return;
            compressed_mer_query.reserve(present_codes.size() + 1024);
            for (uint32_t code : present_codes) {
                uint32_t v = direct_index[code];
                if (v) compressed_mer_query[code] = pack_range(dir_off(v), dir_cnt(v));
            }
            std::vector<uint32_t>().swap(direct_index);
            can_direct = false;
        };

        auto flush_group = [&]() {
            if (!have_cur) return;
            const int sz = static_cast<int>(group.size());
            if (sz > 1) ++shared_mers;
            if (max_assign == 0 || sz <= max_assign) {
                const size_t off = assign_arena.size();
                assign_arena.insert(assign_arena.end(), group.begin(), group.end());
                if (can_direct && (sz > 7 || off > kDirOffMask)) {
                    switch_to_hashmap();
                }
                if (can_direct) {
                    if (direct_index[cur_key] == 0) present_codes.push_back(cur_key);
                    direct_index[cur_key] = pack_dir(static_cast<uint32_t>(off), static_cast<uint32_t>(sz));
                } else {
                    if (compressed_mer_query.find(cur_key) == compressed_mer_query.end()) {
                        present_codes.push_back(cur_key);
                    }
                    compressed_mer_query[cur_key] = pack_range(off, group.size());
                }
                ++kept_mers;
                kept_rows += sz;
            } else {
                ++dropped_mers;
                dropped_rows += sz;
            }
            group.clear();
        };

        auto process_line = [&](const char* p, size_t len) {
            if (len == 0) return;
            const char* e = p + len;
            const char* t1 = static_cast<const char*>(memchr(p, '\t', len));
            if (t1 == nullptr) return;
            const char* t2 = static_cast<const char*>(memchr(t1 + 1, '\t', static_cast<size_t>(e - t1 - 1)));
            if (t2 == nullptr) return;
            const char* t3 = static_cast<const char*>(memchr(t2 + 1, '\t', static_cast<size_t>(e - t2 - 1)));
            const char* loci_end = (t3 != nullptr) ? t3 : e;
            if (mer_size == 0) {
                mer_size = static_cast<int>(t1 - p);
                if (!dir_inited) {
                    dir_inited = true;
                    if (defer_hashmap && mer_size >= 7 && mer_size <= 13) {
                        direct_index.assign(size_t(1) << (2 * mer_size), 0);
                        can_direct = true;
                    } else {
                        compressed_mer_query.reserve(file_size / 25);
                    }
                }
            }

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
            if (!valid) return;
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
        };

        bool mapped_ok = false;
        if (file_size > 0) {
            int fd = ::open(path.c_str(), O_RDONLY);
            if (fd >= 0) {
                void* map_ptr = ::mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
                if (map_ptr != MAP_FAILED) {
                    ::madvise(map_ptr, file_size, MADV_SEQUENTIAL);
                    const char* cur = static_cast<const char*>(map_ptr);
                    const char* end = cur + file_size;
                    while (cur < end) {
                        const char* nl = static_cast<const char*>(
                            memchr(cur, '\n', static_cast<size_t>(end - cur)));
                        size_t line_len = nl ? static_cast<size_t>(nl - cur)
                                             : static_cast<size_t>(end - cur);
                        process_line(cur, line_len);
                        cur = nl ? (nl + 1) : end;
                    }
                    ::munmap(map_ptr, file_size);
                    mapped_ok = true;
                }
                ::close(fd);
            }
        }
        if (!mapped_ok && file_size > 0) {
            std::ifstream file(path);
            if (!file.is_open()) {
                throw std::runtime_error("Cannot open mers file: " + path);
            }
            std::string line;
            while (std::getline(file, line)) {
                process_line(line.data(), line.size());
            }
        }
        flush_group();
        kept_mers_ = present_codes.size();

        if (use_cache && kept_mers_ > 0 && mer_size > 0) {
            bool fits_dir = can_direct;
            if (!fits_dir && mer_size >= 7 && mer_size <= 13) {
                fits_dir = true;
                for (uint32_t code : present_codes) {
                    uint64_t r = compressed_mer_query[code];
                    if (range_off(r) > kDirOffMask || range_cnt(r) > 7) {
                        fits_dir = false;
                        break;
                    }
                }
            }
            std::vector<uint64_t> packed_ranges;
            packed_ranges.reserve(kept_mers_);
            if (can_direct) {
                for (uint32_t code : present_codes) {
                    uint32_t v = direct_index[code];
                    packed_ranges.push_back(pack_range(dir_off(v), dir_cnt(v)));
                }
            } else {
                for (uint32_t code : present_codes) {
                    packed_ranges.push_back(compressed_mer_query[code]);
                }
            }
            const std::string tmp_path = cache_path + ".tmp." + std::to_string(::getpid());
            int wfd = ::open(tmp_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (wfd >= 0) {
                BinHeader hdr{};
                std::memcpy(hdr.magic, "CNEXIDX1", 8);
                hdr.version = 1;
                hdr.mer_size = mer_size;
                hdr.max_assign = max_assign;
                hdr.fits_direct = fits_dir ? 1u : 0u;
                hdr.src_size = static_cast<int64_t>(src_st.st_size);
                hdr.src_mtime_sec = static_cast<int64_t>(src_st.st_mtim.tv_sec);
                hdr.src_mtime_nsec = static_cast<int64_t>(src_st.st_mtim.tv_nsec);
                hdr.kept_mers = static_cast<uint64_t>(kept_mers_);
                hdr.arena_size = static_cast<uint64_t>(assign_arena.size());
                hdr.dropped_mers = static_cast<uint64_t>(dropped_mers);
                hdr.dropped_rows = static_cast<uint64_t>(dropped_rows);
                hdr.shared_mers = static_cast<uint64_t>(shared_mers);
                bool wok = write_exact(wfd, &hdr, sizeof(hdr)) &&
                           write_exact(wfd, present_codes.data(), kept_mers_ * sizeof(uint32_t)) &&
                           write_exact(wfd, packed_ranges.data(), kept_mers_ * sizeof(uint64_t)) &&
                           write_exact(wfd, assign_arena.data(), assign_arena.size() * sizeof(std::pair<int, int>));
                if (::close(wfd) != 0) wok = false;
                if (!wok || ::rename(tmp_path.c_str(), cache_path.c_str()) != 0) {
                    ::unlink(tmp_path.c_str());
                }
            }
        }

        print_filter_stats(kept_mers, kept_rows, dropped_mers, dropped_rows, shared_mers);
    }

    int get_mer_size() const { return mer_size; }
    size_t size() const { return kept_mers_ > 0 ? kept_mers_ : compressed_mer_query.size(); }
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
    if (static_cast<int>(seq.size()) < mer_size) return {-1, 0};

    static auto base_lookup = [] {
        std::array<int8_t, 256> t{};
        t.fill(-1);
        t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
        return t;
    }();

    const int k = mer_size;
    const int chain_len = static_cast<int>(seq.size());
    const uint32_t mer_mask = static_cast<uint32_t>((1ULL << (2 * k)) - 1);
    const int n_codes = chain_len - k + 1;

    thread_local std::vector<std::array<int, 3>> f_hits, r_hits;   // (pos, id, loci)
    f_hits.clear();
    r_hits.clear();
    thread_local std::vector<uint32_t> scan_codes, scan_rcodes;
    scan_codes.assign(n_codes, UINT32_MAX);
    scan_rcodes.assign(n_codes, UINT32_MAX);

    const bool _vp = validate_profile_enabled();
    std::chrono::steady_clock::time_point _t0, _t1;
    if (_vp) _t0 = std::chrono::steady_clock::now();

    // Pass 1: rolling 2-bit codes of every forward k-mer start position, plus
    // each k-mer's reverse-complement code (also via a rolling update), which
    // the canonical index needs to form min(code, rc).
    {
        uint32_t rolling = 0, rc_rolling = 0;
        int valid_run = 0;
        for (int i = 0; i < chain_len; ++i) {
            int8_t bits = base_lookup[static_cast<unsigned char>(seq[i])];
            if (bits < 0) { valid_run = 0; rolling = 0; rc_rolling = 0; continue; }
            rolling = ((rolling << 2) | static_cast<uint32_t>(bits)) & mer_mask;
            rc_rolling = (rc_rolling >> 2) |
                         (static_cast<uint32_t>(3 - bits) << (2 * k - 2));
            ++valid_run;
            if (valid_run >= k) {
                scan_codes[i - k + 1] = rolling;
                scan_rcodes[i - k + 1] = rc_rolling;
            }
        }
    }

    // Legacy per-strand scan (raw direct index, or hash) of an already-encoded
    // code array, appending (position, element, locus) to `out`.
    auto scan_raw = [&](const std::vector<uint32_t>& codes,
                        std::vector<std::array<int, 3>>& out) {
        const int nc = static_cast<int>(codes.size());
        if (compressed_mer_query.has_direct_index()) {
            constexpr int kPf = 32;
            thread_local std::vector<uint32_t> scan_pk;
            scan_pk.assign(nc, 0);
            const int ramp = std::min(kPf, nc);
            for (int s = 0; s < ramp; ++s) {
                uint32_t c = codes[s];
                if (c != UINT32_MAX)
                    __builtin_prefetch(&compressed_mer_query.direct_index[c]);
            }
            for (int s = 0; s < nc; ++s) {
                int pf = s + kPf;
                if (pf < nc) {
                    uint32_t c = codes[pf];
                    if (c != UINT32_MAX)
                        __builtin_prefetch(&compressed_mer_query.direct_index[c]);
                }
                uint32_t c = codes[s];
                if (c == UINT32_MAX) continue;
                scan_pk[s] = compressed_mer_query.direct_index[c];
            }
            for (int s = 0; s < ramp; ++s) {
                uint32_t v = scan_pk[s];
                if (v) __builtin_prefetch(
                    &compressed_mer_query.assign_arena[MerQueryManager::dir_off(v)]);
            }
            for (int s = 0; s < nc; ++s) {
                int pf = s + kPf;
                if (pf < nc) {
                    uint32_t v = scan_pk[pf];
                    if (v) __builtin_prefetch(
                        &compressed_mer_query.assign_arena[MerQueryManager::dir_off(v)]);
                }
                uint32_t v = scan_pk[s];
                if (!v) continue;
                uint32_t off = MerQueryManager::dir_off(v);
                uint32_t cnt = MerQueryManager::dir_cnt(v);
                for (uint32_t a = 0; a < cnt; ++a) {
                    int id = compressed_mer_query.assign_arena[off + a].first;
                    if (target_ele_id != -1 && id != target_ele_id) continue;
                    out.push_back({s, id, compressed_mer_query.assign_arena[off + a].second});
                }
            }
        } else {
            for (int s = 0; s < nc; ++s) {
                uint32_t code = codes[s];
                if (code == UINT32_MAX) continue;
                auto it = compressed_mer_query.compressed_mer_query.find(code);
                if (it == compressed_mer_query.compressed_mer_query.end()) continue;
                uint32_t off = MerQueryManager::range_off(it->second);
                uint32_t cnt = MerQueryManager::range_cnt(it->second);
                for (uint32_t a = 0; a < cnt; ++a) {
                    int id = compressed_mer_query.assign_arena[off + a].first;
                    if (target_ele_id != -1 && id != target_ele_id) continue;
                    out.push_back({s, id, compressed_mer_query.assign_arena[off + a].second});
                }
            }
        }
    };

    if (compressed_mer_query.has_canon_index()) {
        // Canonical (strand-merged): a single forward pass resolves both
        // strands. Code -> canonical key (min(code, rc)) -> range in one random
        // access; each row's form bit says whether the hit is on the forward
        // strand (append at s) or the reverse strand (append at n_codes-1-s).
        constexpr int kPf = 32;
        thread_local std::vector<uint32_t> scan_pk;
        scan_pk.assign(n_codes, 0);
        const int ramp = std::min(kPf, n_codes);
        for (int s = 0; s < ramp; ++s) {
            if (scan_codes[s] == UINT32_MAX) continue;
            __builtin_prefetch(&compressed_mer_query.dir_canon[
                std::min(scan_codes[s], scan_rcodes[s])]);
        }
        for (int s = 0; s < n_codes; ++s) {
            int pf = s + kPf;
            if (pf < n_codes && scan_codes[pf] != UINT32_MAX)
                __builtin_prefetch(&compressed_mer_query.dir_canon[
                    std::min(scan_codes[pf], scan_rcodes[pf])]);
            if (scan_codes[s] == UINT32_MAX) continue;
            scan_pk[s] = compressed_mer_query.dir_canon[
                std::min(scan_codes[s], scan_rcodes[s])];
        }
        for (int s = 0; s < ramp; ++s) {
            uint32_t v = scan_pk[s];
            if (v) __builtin_prefetch(
                &compressed_mer_query.arena_canon[MerQueryManager::dir_off(v)]);
        }
        for (int s = 0; s < n_codes; ++s) {
            int pf = s + kPf;
            if (pf < n_codes) {
                uint32_t v = scan_pk[pf];
                if (v) __builtin_prefetch(
                    &compressed_mer_query.arena_canon[MerQueryManager::dir_off(v)]);
            }
            uint32_t v = scan_pk[s];
            if (!v) continue;
            uint32_t fm = (scan_codes[s] <= scan_rcodes[s]) ? 0u : 1u;
            uint32_t off = MerQueryManager::dir_off(v);
            uint32_t cnt = MerQueryManager::dir_cnt(v);
            for (uint32_t a = 0; a < cnt; ++a) {
                uint32_t row = compressed_mer_query.arena_canon[off + a];
                int id = MerQueryManager::crow_id(row);
                if (target_ele_id != -1 && id != target_ele_id) continue;
                int loci = MerQueryManager::crow_loci(row);
                if (MerQueryManager::crow_form(row) == fm)
                    f_hits.push_back({s, id, loci});
                else
                    r_hits.push_back({n_codes - 1 - s, id, loci});
            }
        }
        // Reverse hits were appended with position decreasing; the legacy
        // reverse scan visits them ascending. Reverse the list (the chain only
        // depends on across-position order, which is now ascending).
        std::reverse(r_hits.begin(), r_hits.end());
    } else {
        scan_raw(scan_codes, f_hits);
        // Reverse strand: encode the reverse-complement codes of the forward
        // k-mers, then scan in ascending rc-read position order.
        thread_local std::vector<uint32_t> rc_codes;
        rc_codes.assign(n_codes, UINT32_MAX);
        {
            uint32_t rolling = 0, rc_rolling = 0;
            int valid_run = 0;
            for (int i = 0; i < chain_len; ++i) {
                int8_t bits = base_lookup[static_cast<unsigned char>(seq[i])];
                if (bits < 0) { valid_run = 0; rolling = 0; rc_rolling = 0; continue; }
                rolling = ((rolling << 2) | static_cast<uint32_t>(bits)) & mer_mask;
                rc_rolling = (rc_rolling >> 2) |
                             (static_cast<uint32_t>(3 - bits) << (2 * k - 2));
                ++valid_run;
                if (valid_run >= k) rc_codes[chain_len - k - (i - k + 1)] = rc_rolling;
            }
        }
        scan_raw(rc_codes, r_hits);
    }

    if (_vp) {
        g_val_profile.rd_scan_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t0).count();
        _t1 = std::chrono::steady_clock::now();
    }

    // Vote + gate matrix (T1 OR T2 ...) for one strand frame.
    thread_local std::vector<int> frame_ids;
    auto eval_frame = [&](const std::vector<std::array<int, 3>>& hits,
                          int& r_ele, int& r_votes) -> std::pair<int, int> {
        if (hits.empty()) return {-1, 0};
        frame_ids.clear();
        frame_ids.reserve(hits.size());
        for (const auto& h : hits) frame_ids.push_back(h[1]);
        VoteTally tally = tally_votes(frame_ids);
        std::pair<int, int> cached_chain{0, 0};
        int cached_chain_id = -1;
        auto chain_of = [&](int confi_id) -> std::pair<int, int> {
            if (confi_id == cached_chain_id) return cached_chain;
            int first = -1, last = -1, prev_x = -1, prev_l = -1, passing = 0, m = 0;
            for (const auto& h : hits) {
                if (h[1] != confi_id) continue;
                int x = h[0], l = h[2];
                if (first < 0) first = x;
                else if (std::abs((x - prev_x) - (l - prev_l)) <= cfg.max_diff) ++passing;
                prev_x = x; prev_l = l; last = x; ++m;
            }
            cached_chain_id = confi_id;
            cached_chain = (m == 0)
                ? std::pair<int,int>{0, 0}
                : std::pair<int,int>{passing, last - first + mer_size};
            if (m == 0) return cached_chain;
            return cached_chain;
        };
        int passing = 0, span = 0;
        const bool gate_ok = gate_pass(tally, cfg, chain_of, passing, span);
        if (!gate_ok) {
            auto rg = rescue_gate(tally, cfg, chain_of);
            if (rg.first >= 0) { r_ele = rg.first; r_votes = rg.second; }
            return {-1, 0};
        }
        return {tally.best_id, static_cast<int>(tally.n)};
    };

    int f_ele = -1, f_votes = 0, r_ele = -1, r_votes = 0;
    auto [forward_id, forward_count] = eval_frame(f_hits, f_ele, f_votes);
    auto [reverse_id, reverse_count] = eval_frame(r_hits, r_ele, r_votes);

    if (_vp) {
        g_val_profile.rd_gate_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - _t1).count();
        g_val_profile.rd_reads += 1;
    }

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

    thread_local std::vector<uint32_t> fwd_dir;
    thread_local std::vector<uint32_t> rc_dir;
    thread_local std::vector<uint64_t> fwd_ranges;
    thread_local std::vector<uint64_t> rc_ranges;
    thread_local std::vector<uint64_t> upper_bits;
    thread_local std::vector<uint32_t> codes_fwd;
    thread_local std::vector<uint32_t> codes_rc;
    int cur_o = 0;   // current window origin within the segment
    const bool vp = validate_profile_enabled();
    const bool use_dir = mqm.has_direct_index();
    thread_local long long tl_prep_ns = 0;
    thread_local long long tl_win_ns = 0;
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

    // Incremental tally for sliding windows. Counts are updated by adding the
    // new right edge and removing the old left edge. The exact legacy tally is
    // recomputed only when the top count is tied (the only case where the
    // first-seen order can change the winning element id).
    // Uses a contiguous power-of-2 ring buffer instead of std::deque.
    struct SlidingTally {
        std::vector<int> cnt;
        std::vector<int> pos;
        std::vector<int> active;
        std::vector<int> pos_ring;
        uint32_t ring_head = 0;
        uint32_t ring_tail = 0;
        uint32_t ring_mask = 0;
        int total = 0;

        void ensure_ring(int nspan_val) {
            uint32_t need = static_cast<uint32_t>(nspan_val) + 4u;
            if (ring_mask + 1 >= need) return;
            uint32_t cap = 16;
            while (cap < need) cap <<= 1;
            pos_ring.resize(cap);
            ring_mask = cap - 1;
        }
        void clear() {
            for (int id : active) { cnt[id] = 0; pos[id] = -1; }
            active.clear();
            ring_head = 0;
            ring_tail = 0;
            total = 0;
        }
        inline void push_hit(int p) {
            pos_ring[ring_tail++ & ring_mask] = p;
        }
        inline void pop_hit_if(int p) {
            if (ring_head != ring_tail && pos_ring[ring_head & ring_mask] == p) {
                ++ring_head;
            }
        }
        inline void add(int id) {
            if (id >= static_cast<int>(cnt.size())) {
                cnt.resize(id + 1, 0);
                pos.resize(id + 1, -1);
            }
            if (cnt[id]++ == 0) {
                pos[id] = static_cast<int>(active.size());
                active.push_back(id);
            }
            ++total;
        }
        inline void remove(int id) {
            --total;
            if (--cnt[id] == 0) {
                int p = pos[id];
                int last = active.back();
                active[p] = last;
                pos[last] = p;
                active.pop_back();
                pos[id] = -1;
            }
        }
        VoteTally finish(bool& tied) const {
            VoteTally r;
            r.n = total;
            int best = 0, second = 0, best_id = -1, ties = 0;
            for (int id : active) {
                int c = cnt[id];
                if (c > best) {
                    second = best;
                    best = c;
                    best_id = id;
                    ties = 1;
                } else if (c == best) {
                    ++ties;
                    second = best;
                } else if (c > second) {
                    second = c;
                }
            }
            r.best_id = best_id;
            r.best_cnt = best;
            r.second_cnt = second;
            tied = ties > 1;
            return r;
        }
    };
    thread_local SlidingTally slide_f, slide_r;
    slide_f.ensure_ring(nspan);
    slide_r.ensure_ring(nspan);

    auto get_range = [&](int p, bool rc, uint32_t& off, uint32_t& cnt) -> bool {
        if (use_dir) {
            uint32_t v = (rc ? rc_dir : fwd_dir)[p];
            if (!v) return false;
            off = MerQueryManager::dir_off(v);
            cnt = MerQueryManager::dir_cnt(v);
            return true;
        }
        uint64_t r = (rc ? rc_ranges : fwd_ranges)[p];
        if (!r) return false;
        off = MerQueryManager::range_off(r);
        cnt = MerQueryManager::range_cnt(r);
        return true;
    };

    auto make_chain = [&](const SlidingTally* slide, bool rc_strand) {
        return [&, slide, rc_strand](int confi_id) -> std::pair<int, int> {
            if (vp) ++g_val_profile.chain_calls;
            int first = -1, last = -1, prev_x = -1, prev_l = -1, passing = 0, m = 0;
            const uint32_t npos = slide->ring_tail - slide->ring_head;
            for (uint32_t qi = 0; qi < npos; ++qi) {
                uint32_t slot = rc_strand
                    ? (slide->ring_tail - 1 - qi)
                    : (slide->ring_head + qi);
                int p = slide->pos_ring[slot & slide->ring_mask];
                int x = rc_strand ? (nspan - 1 - (p - cur_o)) : (p - cur_o);
                uint32_t off = 0, cnt = 0;
                if (!get_range(p, rc_strand, off, cnt)) continue;
                for (uint32_t a = 0; a < cnt; ++a) {
                    const auto& pr = mqm.assign_arena[off + a];
                    int id = pr.first;
                    if (id != confi_id) continue;
                    int l = pr.second;
                    if (first < 0) first = x;
                    else if (std::abs((x - prev_x) - (l - prev_l)) <= cfg.max_diff) ++passing;
                    prev_x = x; prev_l = l; last = x; ++m;
                }
            }
            if (m == 0) return {0, 0};
            return {passing, last - first + k};
        };
    };
    auto chain_fwd = make_chain(&slide_f, false);
    auto chain_rc = make_chain(&slide_r, true);

    int min_req_c = INT_MAX;
    for (const auto& tier : cfg.tiers) min_req_c = std::min(min_req_c, tier.min_c);
    if (cfg.verify.enabled) min_req_c = std::min(min_req_c, cfg.verify.t3.min_c);
    const int min_any_hits = (min_req_c == INT_MAX) ? 0 : (min_req_c + 1);

    const int seg_stride = SEG - (window_size - 1);
    for (int seg_start = 0; seg_start < L; seg_start += seg_stride) {
        const int seg_end = std::min(L, seg_start + SEG);
        const int seg_len = seg_end - seg_start;

        if (use_dir) {
            fwd_dir.assign(seg_len, 0);
            rc_dir.assign(seg_len, 0);
        } else {
            fwd_ranges.assign(seg_len, 0);
            rc_ranges.assign(seg_len, 0);
        }
        upper_bits.assign((seg_len + 63) >> 6, 0ULL);

        auto _t_prep = std::chrono::steady_clock::now();
        if (use_dir) {
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
                fwd_dir[st] = mqm.direct_index[cf];
                rc_dir[st] = mqm.direct_index[codes_rc[st]];
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

        auto _t_win = std::chrono::steady_clock::now();

        int s0 = seg_start;
        if (step_size > 1) {
            int rem = (s0 + start_offset) % step_size;
            if (rem != 0) s0 += step_size - rem;
        }
        const int last_start = seg_end - window_size;
        double min_gate_ratio = 1e9;
        for (const auto& tier : cfg.tiers)
            min_gate_ratio = std::min(min_gate_ratio, tier.vote_ratio);
        if (cfg.verify.enabled)
            min_gate_ratio = std::min(min_gate_ratio, cfg.verify.t3.vote_ratio);
        const bool tie_needs_exact = min_gate_ratio <= 1.0;

        auto add_pos = [&](int p, bool rc, SlidingTally& st) {
            uint32_t off = 0, cnt = 0;
            if (!get_range(p, rc, off, cnt)) return;
            for (uint32_t a = 0; a < cnt; ++a)
                st.add(mqm.assign_arena[off + a].first);
            st.push_hit(p);
        };
        auto remove_pos = [&](int p, bool rc, SlidingTally& st) {
            uint32_t off = 0, cnt = 0;
            if (!get_range(p, rc, off, cnt)) return;
            for (uint32_t a = 0; a < cnt; ++a)
                st.remove(mqm.assign_arena[off + a].first);
            st.pop_hit_if(p);
        };
        auto rebuild_slide = [&](int o, bool rc, SlidingTally& st) {
            st.clear();
            for (int p = o; p < o + nspan; ++p) add_pos(p, rc, st);
        };
        auto full_tally = [&](int o, bool rc) {
            tally_acc.begin();
            if (!rc) {
                for (int p = o; p < o + nspan; ++p) {
                    uint32_t off = 0, cnt = 0;
                    if (!get_range(p, false, off, cnt)) continue;
                    for (uint32_t a = 0; a < cnt; ++a)
                        tally_acc.add(mqm.assign_arena[off + a].first);
                }
            } else {
                for (int p = o + nspan - 1; p >= o; --p) {
                    uint32_t off = 0, cnt = 0;
                    if (!get_range(p, true, off, cnt)) continue;
                    for (uint32_t a = 0; a < cnt; ++a)
                        tally_acc.add(mqm.assign_arena[off + a].first);
                }
            }
            return tally_acc.finish();
        };

        bool slide_init = false;
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
            if (!has_upper) { slide_init = false; continue; }

            cur_o = o;
            if (vp) ++g_val_profile.tally_calls;
            auto _tt = vp ? std::chrono::steady_clock::now()
                          : std::chrono::steady_clock::time_point{};
            if (!slide_init || step_size >= window_size) {
                rebuild_slide(o, false, slide_f);
                rebuild_slide(o, true, slide_r);
                slide_init = true;
            } else {
                const int prev_o = o - step_size;
                for (int p = prev_o; p < prev_o + step_size; ++p) {
                    remove_pos(p, false, slide_f);
                    remove_pos(p, true, slide_r);
                }
                for (int p = prev_o + nspan; p < prev_o + nspan + step_size; ++p) {
                    add_pos(p, false, slide_f);
                    add_pos(p, true, slide_r);
                }
            }
            int fwd_id = -1, fwd_count = 0;
            bool fwd_ok = false;
            VoteTally fwd_t;
            if (slide_f.total >= min_any_hits) {
                bool tie_f = false;
                fwd_t = slide_f.finish(tie_f);
                if (tie_f && tie_needs_exact) fwd_t = full_tally(o, false);
                if (gate_ratio_precheck(fwd_t, cfg)) {
                    int passing = 0, span = 0;
                    if (gate_pass(fwd_t, cfg, chain_fwd, passing, span)) {
                        fwd_ok = true;
                        fwd_id = fwd_t.best_id;
                        fwd_count = fwd_t.n;
                    }
                }
            }

            int rev_id = -1, rev_count = 0;
            bool rev_ok = false;
            VoteTally rev_t;
            if (slide_r.total >= min_any_hits) {
                bool tie_r = false;
                rev_t = slide_r.finish(tie_r);
                if (tie_r && tie_needs_exact) rev_t = full_tally(o, true);
                if (gate_ratio_precheck(rev_t, cfg)) {
                    int passing = 0, span = 0;
                    if (gate_pass(rev_t, cfg, chain_rc, passing, span)) {
                        rev_ok = true;
                        rev_id = rev_t.best_id;
                        rev_count = rev_t.n;
                    }
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
        g_val_profile.tally_ns += tl_tally_ns;
        tl_prep_ns = 0;
        tl_win_ns = 0;
        tl_tally_ns = 0;
    }
}
