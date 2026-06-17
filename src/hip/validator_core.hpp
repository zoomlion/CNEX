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


inline std::string reverse_complement(const std::string& dna) {
    static const std::array<char, 256> complement = [] {
        std::array<char, 256> table = {};
        table['A'] = 'T'; table['T'] = 'A';
        table['C'] = 'G'; table['G'] = 'C';
        return table;
    }();

    std::string result;
    result.reserve(dna.length());
    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        result.push_back(complement[(unsigned char)*it]);
    }
    return result;
}


inline int findFrequentWithMap(const std::vector<int>& nums) {
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

    float threshold = static_cast<float>(n) / 10.0f;
    if (best_cnt > threshold && best_cnt >= 3 * second_cnt) {
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
    int min_c)
{
    auto is_sorted = [](const std::vector<int>& lst) {
        return std::is_sorted(lst.begin(), lst.end());
    };

    auto is_pattern_gapped = [](
        const std::vector<int>& lst1,
        const std::vector<int>& lst2,
        const int min_c,
        const int max_g = 1) {
        if (static_cast<int>(lst1.size()) < min_c) {
            return false;
        }
        int count = 0;
        for (size_t i = 0; i < lst1.size() - 1; ++i) {
            if (std::abs(std::abs(lst1[i + 1] - lst1[i]) - std::abs(lst2[i + 1] - lst2[i])) <= max_g) {
                if (++count >= min_c) return true;
            }
        }
        return false;
    };

    auto validate_chain = [&](const std::string& chain) -> std::pair<int, int> {
        tsl::robin_map<int, std::vector<std::pair<int, int>>> ordinals;
        std::vector<int> candidates;
        int chain_len = static_cast<int>(chain.size());
        int mer_size_local = mer_size;

        if (chain_len < mer_size_local) {
            return {-1, 0};
        }

        // Rolling hash: O(1) per position instead of O(k) per position
        static auto base_lookup = [] {
            std::array<int8_t, 256> t{};
            t.fill(-1);
            t['A'] = 0; t['C'] = 1; t['G'] = 2; t['T'] = 3;
            return t;
        }();

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
                    candidates.push_back(id);
                    ordinals[id].emplace_back(start, static_cast<int>(loci));
                }
            }
        }

        int confi_id = findFrequentWithMap(candidates);
        if (confi_id == -1) {
            return {-1, 0};
        }

        const auto& positions = ordinals[confi_id];
        std::vector<int> seq_positions;
        std::vector<int> loci_positions;
        for (const auto& pos : positions) {
            seq_positions.push_back(pos.first);
            loci_positions.push_back(pos.second);
        }

        if (is_pattern_gapped(seq_positions, loci_positions, min_c) && is_sorted(loci_positions)) {
            return {confi_id, static_cast<int>(candidates.size())};
        }

        return {-1, 0};
    };

    auto [forward_id, forward_count] = validate_chain(seq);

    std::string rev_comp = reverse_complement(seq);
    auto [reverse_id, reverse_count] = validate_chain(rev_comp);

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
