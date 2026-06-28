#pragma once
#include <tsl/robin_map.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

inline float mer_entropy(const std::string& seq) {
    for (char c : seq) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return 0.0f;
    }
    tsl::robin_map<char, int> freq;
    for (char c : seq) freq[c]++;
    float entropy = 0.0f;
    float len = static_cast<float>(seq.size());
    for (const auto& p : freq) {
        float prob = static_cast<float>(p.second) / len;
        entropy -= prob * std::log2(prob);
    }
    return entropy;
}

class TableMer {
public:
    tsl::robin_map<uint32_t, std::vector<std::pair<int, int>>> merTable;
    tsl::robin_map<uint32_t, int> merCounter;
    static uint8_t base_to_int(char c) {
        switch (c) {
            case 'A': return 0b00;
            case 'C': return 0b01;
            case 'G': return 0b10;
            case 'T': return 0b11;
            default: return 0xFF;
        }
    }
    size_t k;
    float min_entropy;

    TableMer(size_t input_k = 13, float input_min_entropy = 1.4f)
        : k(input_k), min_entropy(input_min_entropy) {}

    void set_min_entropy(float v) { min_entropy = v; }

    uint32_t mer2int(const std::string& mer) const {
        uint32_t result = 0;
        for (char base : mer) {
            uint8_t bits = base_to_int(base);
            if (bits == 0xFF) return 0;
            result = (result << 2) | bits;
        }
        return result;
    }

    std::string int2mer(uint32_t value) const {
        static const char INT_TO_BASE[4] = {'A', 'C', 'G', 'T'};
        std::string mer;
        for (size_t i = 0; i < k; ++i) {
            mer.push_back(INT_TO_BASE[value & 0b11]);
            value >>= 2;
        }
        std::reverse(mer.begin(), mer.end());
        return mer;
    }

    void add(const std::string& dna, int bunch_id, int loci, int local_count) {
        if (mer_entropy(dna) < min_entropy) return;
        uint32_t key = mer2int(dna);
        auto it = merTable.find(key);
        if (it == merTable.end()) {
            merCounter[key] = local_count;
            merTable[key] = {{bunch_id, loci}};
        } else {
            if (local_count > merCounter[key]) {
                merCounter[key] = local_count;
                merTable[key] = {{bunch_id, loci}};
            } else if (local_count == merCounter[key]) {
                merTable[key].push_back({bunch_id, loci});
            }
        }
    }

    std::tuple<int, int, int> get(uint32_t key) const {
        auto it = merTable.find(key);
        if (it == merTable.end() || it->second.empty()) return {-1, -1, 0};
        std::mt19937 generator(42);
        std::uniform_int_distribution<int> dist(0, static_cast<int>(it->second.size()) - 1);
        auto p = it->second[dist(generator)];
        auto cnt = merCounter.find(key);
        return {p.first, p.second, cnt != merCounter.end() ? cnt->second : 0};
    }

    void dump(const std::string& filename) const {
        std::ofstream out(filename);
        for (const auto& [key, _] : merTable) {
            auto [bunch_id, loci, count] = get(key);
            out << int2mer(key) << "\t" << bunch_id << "\t" << loci << "\t" << count << "\n";
        }
    }

    size_t size() const { return merTable.size(); }
};
