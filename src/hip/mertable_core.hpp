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
    struct MerAssign {
        int bunch_id;
        int loci;
        int count;
    };
    tsl::robin_map<uint32_t, std::vector<MerAssign>> merTable;
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
    int max_assign = 5;

    TableMer(size_t input_k = 13, float input_min_entropy = 1.4f)
        : k(input_k), min_entropy(input_min_entropy) {}

    void set_min_entropy(float v) { min_entropy = v; }
    void set_max_assign(int v) { max_assign = v; }

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
        auto& assigns = merTable[key];
        for (auto& a : assigns) {
            if (a.bunch_id == bunch_id && a.loci == loci) {
                if (local_count > a.count) a.count = local_count;
                return;
            }
        }
        assigns.push_back({bunch_id, loci, local_count});
    }

    void dump(const std::string& filename) {
        std::ofstream out(filename);
        for (const auto& [key, assigns] : merTable) {
            // A k-mer with too many loci/elements is a repeat signal: drop it entirely
            if (static_cast<int>(assigns.size()) > max_assign) continue;
            std::string mer = int2mer(key);
            for (const auto& a : assigns) {
                out << mer << "\t" << a.bunch_id << "\t" << a.loci << "\t" << a.count << "\n";
            }
        }
    }

    size_t size() const { return merTable.size(); }
};
