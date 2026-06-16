// MIT License
//
// Copyright (c) 2024 JiangminZheng
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tsl/robin_map.h>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>
#include <utility>
#include <fstream>

namespace py = pybind11;

uint32_t dna_encoder(const std::string& dna) {
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
        if (lookup[c] == 0 && c != 'A') {
            throw std::invalid_argument("Invalid DNA character: " + std::string(1, c));
        }
        result = (result << 2) | lookup[c];
    }
    return result;
}

bool encode_mer_at(const std::string& seq, size_t pos, size_t k, uint32_t& result) {
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

std::string reverse_complement(const std::string& dna) {
    static const std::array<char, 256> complement = [] {
        std::array<char, 256> table = {};
        table['A'] = 'T'; table['T'] = 'A';
        table['C'] = 'G'; table['G'] = 'C';
        return table;
    }();

    std::string result;
    result.reserve(dna.length());
    for (auto it = dna.rbegin(); it != dna.rend(); ++it) {
        result.push_back(complement[*it]);
    }
    return result;
}

int findFrequentWithMap(const std::vector<int>& nums) {
    tsl::robin_map<int, int> frequency;
    // Count frequency of each number
    for (int num : nums) {
        frequency[num]++;
    }
    
    // Convert map to a vector of pairs and sort by frequency (descending)
    std::vector<std::pair<int, int>> freq_vec(frequency.begin(), frequency.end());
    std::sort(freq_vec.begin(), freq_vec.end(), [](const auto& a, const auto& b) {
        return a.second > b.second; // Sort by frequency in descending order
    });
    
    float threshold = (float)nums.size() / 10;
    
    if (freq_vec.size() <= 1) {
        if (!freq_vec.empty() && freq_vec[0].second > threshold) {
            return freq_vec[0].first;
        }
        return -1;
    }
    
    int max_freq = freq_vec[0].second;
    int second_freq = freq_vec[1].second;
    
    if (max_freq > threshold && max_freq >= 3 * second_freq) {
        return freq_vec[0].first;
    }
    
    return -1;
}

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
        // Estimate line count from file size for pre-allocation
        file.seekg(0, std::ios::end);
        size_t file_size = static_cast<size_t>(file.tellg());
        file.seekg(0, std::ios::beg);
        compressed_mer_query.reserve(file_size / 25);

        std::string line;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            // Format: mer\tid\tloci\tcount
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

std::pair<int, int> validate_read(
    const std::string& seq, 
    const MerQueryManager& compressed_mer_query, 
    int mer_size, 
    int min_c) {
    auto is_sorted = [](const std::vector<int>& lst) {
        return std::is_sorted(lst.begin(), lst.end());
    };

    auto is_pattern_gapped = [](
        const std::vector<int>& lst1, 
        const std::vector<int>& lst2, 
        const int min_c, 
        const int max_g = 1) {
        std::vector<int> patterned_gaps;
        if (static_cast<int>(lst1.size()) < min_c) {
            return false;
        }
        for (size_t i = 0; i < lst1.size() - 1; ++i) {
            if (std::abs(std::abs(lst1[i + 1] - lst1[i]) - std::abs(lst2[i + 1] - lst2[i])) <= max_g) {
                patterned_gaps.push_back(i);
            }
        }
        return static_cast<int>(patterned_gaps.size()) >= min_c;
    };

    auto validate_chain = [&](const std::string& chain) -> std::pair<int, int> {
        tsl::robin_map<int, std::vector<std::pair<int, int>>> ordinals;
        std::vector<int> candidates;
        uint32_t encoded_mer;

        for (size_t i = 0; i <= chain.size() - static_cast<size_t>(mer_size); ++i) {
            if (!encode_mer_at(chain, i, static_cast<size_t>(mer_size), encoded_mer)) {
                continue;
            }
            auto it = compressed_mer_query.compressed_mer_query.find(encoded_mer);
            if (it != compressed_mer_query.compressed_mer_query.end()) {
                auto [id, loci] = it->second;
                candidates.push_back(id);
                ordinals[id].emplace_back(static_cast<int>(i), static_cast<int>(loci));
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

    // Validate forward chain
    auto [forward_id, forward_count] = validate_chain(seq);
    
    // Validate reverse complement chain
    std::string rev_comp = reverse_complement(seq);
    auto [reverse_id, reverse_count] = validate_chain(rev_comp);

    // Compare results and return the best match
    if (forward_id == -1 && reverse_id == -1) {
        return {-1, 0}; // No valid match found
    } else if (forward_id == -1) {
        return {reverse_id, -1}; // Reverse complement match
    } else if (reverse_id == -1) {
        return {forward_id, 1}; // Forward match
    } else {
        // Both chains match, return the one with more matching k-mers
        return (forward_count >= reverse_count) ? 
               std::make_pair(forward_id, 1) : 
               std::make_pair(reverse_id, -1);
    }
}

PYBIND11_MODULE(validator, m) {
    m.doc() = "A Python module for validating DNA reads using compressed mer query";

    py::class_<MerQueryManager>(m, "MerQueryManager")
        .def(py::init<>())
        .def("add_mer", &MerQueryManager::add_mer, "Add a mer to the query map",
             py::arg("mer"), py::arg("id"), py::arg("loci"))
        .def("load_from_file", &MerQueryManager::load_from_file, 
             "Load mers from TSV file directly (C++ fast path)",
             py::arg("path"))
        .def("get_mer_size", &MerQueryManager::get_mer_size, "Get k-mer size")
        .def("size", &MerQueryManager::size, "Number of mers in the table");

    m.def("validate_read", &validate_read, "Validate a read sequence using compressed mer query",
          py::arg("seq"), py::arg("compressed_mer_query"), py::arg("mer_size"), py::arg("min_c"));
    m.def("dna_encoder", &dna_encoder, "Encode a DNA string to an integer",
          py::arg("dna"));
    m.def("encode_mer_at", &encode_mer_at, "Encode a k-mer at a position in a string (zero-alloc)",
          py::arg("seq"), py::arg("pos"), py::arg("k"), py::arg("result"));
}