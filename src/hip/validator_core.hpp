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




inline int findFrequentWithMap(const std::vector<int>& nums, double vote_frac=0.1, double vote_ratio=3.0) {
    if (nums.empty()) return -1;

    int n = static_cast<int>(nums.size());

    thread_local std::vector<int> sorted;
    sorted.assign(nums.begin(), nums.end());
    std::sort(sorted.begin(), sorted.end());

    int best_id = sorted[0];
    int best_cnt = 1;
    int second_cnt = 0;
    int cur_id = sorted[0];
    int cur_cnt = 1;

    for (int i = 1; i < n; ++i) {
        if (sorted[i] == cur_id) {
            ++cur_cnt;
        } else {
            if (cur_cnt > best_cnt) {
                second_cnt = best_cnt;
                best_cnt = cur_cnt;
                best_id = cur_id;
            } else if (cur_cnt > second_cnt) {
                second_cnt = cur_cnt;
            }
            cur_id = sorted[i];
            cur_cnt = 1;
        }
    }
    if (cur_cnt > best_cnt) {
        second_cnt = best_cnt;
        best_cnt = cur_cnt;
        best_id = cur_id;
    } else if (cur_cnt > second_cnt) {
        second_cnt = cur_cnt;
    }

    float threshold = static_cast<float>(n) * static_cast<float>(vote_frac);
    if (best_cnt > threshold && best_cnt >= static_cast<int>(vote_ratio * second_cnt)) {
        return best_id;
    }
    return -1;
}


// ─── MerQueryManager ───

class MerQueryManager {
public:
    tsl::robin_map<uint32_t, std::pair<int, int>> compressed_mer_query;
    int mer_size = 0;

    void add_mer(const std::string& mer, int id, int loci) {
        uint32_t encoded_mer = dna_encoder(mer);
        compressed_mer_query[encoded_mer] = {id, loci};
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

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            size_t t1 = line.find('\t');
            if (t1 == std::string::npos) continue;
            size_t t2 = line.find('\t', t1 + 1);
            if (t2 == std::string::npos) continue;
            size_t t3 = line.find('\t', t2 + 1);
            if (t3 == std::string::npos) continue;
            std::string mer = line.substr(0, t1);
            int id = std::stoi(line.substr(t1 + 1, t2 - t1 - 1));
            int loci = std::stoi(line.substr(t2 + 1, t3 - t2 - 1));
            if (mer_size == 0) mer_size = static_cast<int>(mer.length());
            try {
                uint32_t encoded_mer = dna_encoder(mer);
                compressed_mer_query[encoded_mer] = {id, loci};
            } catch (const std::invalid_argument&) {
                continue;
            }
        }
    }

    int get_mer_size() const { return mer_size; }
    size_t size() const { return compressed_mer_query.size(); }
};


// ─── validate_read ───

inline std::pair<int, int> validate_read(
    const std::string& seq,
    const MerQueryManager& compressed_mer_query,
    int mer_size,
    int min_c = 7,
    int max_diff = 2,
    int min_span = 25,
    int target_ele_id = -1,
    double vote_frac = 0.1,
    double vote_ratio = 3.0)
{
    thread_local std::vector<int> match_ids;
    thread_local std::vector<int> match_starts;
    thread_local std::vector<int> match_loci;
    thread_local std::vector<int> seq_positions;
    thread_local std::vector<int> loci_positions;

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
                    auto [id, loci] = it->second;
                    if (target_ele_id == -1 || id == target_ele_id) {
                        match_ids.push_back(id);
                        match_starts.push_back(start);
                        match_loci.push_back(static_cast<int>(loci));
                    }
                }
            }
        }

        if (match_ids.empty()) return {-1, 0};

        // ─── Vote + filtered check ───
        int confi_id = findFrequentWithMap(match_ids, vote_frac, vote_ratio);
        if (confi_id == -1) return {-1, 0};

        seq_positions.clear();
        loci_positions.clear();
        size_t n = match_ids.size();
        for (size_t i = 0; i < n; ++i) {
            if (match_ids[i] == confi_id) {
                seq_positions.push_back(match_starts[i]);
                loci_positions.push_back(match_loci[i]);
            }
        }

        size_t m = seq_positions.size();
        int passing = 0;
        for (size_t i = 1; i < m; ++i) {
            int gap_read = seq_positions[i] - seq_positions[i - 1];
            int gap_loci = loci_positions[i] - loci_positions[i - 1];
            if (std::abs(gap_read - gap_loci) <= max_diff)
                passing++;
        }

        int span = seq_positions.back() - seq_positions.front() + mer_size;

        if (passing >= min_c && span >= min_span)
            return {confi_id, static_cast<int>(n)};

        return {-1, 0};
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
