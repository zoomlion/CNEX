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
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <queue>
#include <set>
#include <random>

namespace py = pybind11;

float mer_entropy(std::string seq) {
    // if non-ACGT characters are found, return 0.0
    for (char nucleotide : seq) {
        if (nucleotide!= 'A' && nucleotide!= 'C' && nucleotide!= 'G' && nucleotide!= 'T') {
            return 0.0;
        }
    }

    // Calculate the entropy of a DNA sequence
    tsl::robin_map<char, int> frequency;
    int length = seq.length();

    // Count the frequency of each nucleotide
    for (char nucleotide : seq) {
        frequency[nucleotide]++;
    }

    // Calculate entropy
    float entropy = 0.0;
    for (const auto& pair : frequency) {
        float probability = static_cast<float>(pair.second) / length;
        entropy -= probability * log2(probability);
    }

    return entropy;
}

// set parameter k for k-mer when building the merTable
class TableMer{
public:
    tsl::robin_map<u_int32_t, std::vector<std::pair<int, int>> > merTable; 
    tsl::robin_map<u_int32_t, int> merCounter; 
    const tsl::robin_map<char, uint8_t> BASE_TO_INT = {
        {'A', 0b00}, {'C', 0b01}, {'G', 0b10}, {'T', 0b11}
    };

    size_t k;
    float min_entropy; // entropy threshold for filtering out low-entropy kmers

    TableMer(size_t input_k = 13, float input_min_entropy = 1.4) 
        : k(input_k), min_entropy(input_min_entropy) {
        // Validate kmer_size
        if (k > 16) {
            throw std::invalid_argument("k must be <= 16 due to 32-bit integer limitations");
        }
        if (k < 7) {
            throw std::invalid_argument("k must be >= 7");
        }
        // Validate entropy threshold
        if (min_entropy <= 0) {
            throw std::invalid_argument("min_entropy must be > 0");
        }
    }

    // Setter for k
    void set_k(size_t k) {
        if (k > 16) {
            throw std::invalid_argument("k must be <= 16 due to 32-bit integer limitations");
        }
        if (k < 7) {
            throw std::invalid_argument("k must be >= 7");
        }
        k = k;
    }

    // Setter for min_entropy
    void set_min_entropy(float min_entropy) {
        if (min_entropy <= 0) {
            throw std::invalid_argument("min_entropy must be > 0");
        }
        min_entropy = min_entropy;
    }

    uint32_t mer2int(const std::string &mer){
        uint32_t result = 0;
        for (char base : mer)
        {
            auto it = BASE_TO_INT.find(base);
            if (it == BASE_TO_INT.end())
            {
                return 0; // Handle invalid characters
            }
            result = (result << 2) | it->second;
        }
        return result;
    }

    std::string int2mer(uint32_t value, size_t k){
        std::string mer;
        for (size_t i = 0; i < k; ++i)
        {
            uint8_t base_value = value & 0b11;
            mer.push_back(INT_TO_BASE[base_value]);
            value >>= 2;
        }
        reverse(mer.begin(), mer.end());
        return mer;
    }

    const std::vector<char> INT_TO_BASE = {'A', 'C', 'G', 'T'};

    // add func: add a new key-value(pair) pair to the merTable
    void add(std::string dna, int bunch_id, int loci, int local_count){
        // check entropy over min_entropy
        if (float(mer_entropy(dna)) < min_entropy) {
            return;
        }
        // convert string to uint32_t
        u_int32_t key = mer2int(dna);
        // not in the merTable, add a new key-value pair
        if (merTable.find(key) == merTable.end()){
            merCounter[key] = local_count;
            std::vector<std::pair<int, int>> temp;
            temp.push_back(std::make_pair(bunch_id, loci));
            merTable[key] = temp;
        }
        else
        {
            // already in the merTable, add a new value to the existing key
            if (local_count > merCounter[key]){
                merCounter[key] = local_count;
                std::vector<std::pair<int, int>> temp;
                temp.push_back(std::make_pair(bunch_id, loci));
                merTable[key] = temp;
            }
            else if (local_count == merCounter[key]){
                std::vector<std::pair<int, int>> temp = merTable[key];
                temp.push_back(std::make_pair(bunch_id, loci));
                merTable[key] = temp;
            }
        }
    }

    std::tuple<int, int, int> get(u_int32_t key) {
        // get random and stable one (bunch_id, loci) pair from the merTable
        // set seed to 42
        short seed = 42; 
        std::mt19937 generator(seed); 
        std::uniform_int_distribution<int> distribution(0, merTable[key].size() - 1);
        int randomIndex = distribution(generator);
        std::pair<int, int> infos = merTable[key][randomIndex];
        int count = merCounter[key];

        return std::make_tuple(infos.first, infos.second, count);
    }

    void dump(std::string filename) {
        // dump random and stable one (bunch_id, loci) pair from the merTable to a file
        std::ofstream outfile(filename);
        for (const auto& pair : merTable) {
            std::string mer = int2mer(pair.first, k);
            std::tuple<int, int, int> infos = get(pair.first);
            // sep with tab
            outfile << mer << "\t" << std::get<0>(infos) << "\t" << std::get<1>(infos) << "\t" << std::get<2>(infos) << std::endl;
        }
    }
}; 


PYBIND11_MODULE(tablemer, m) {
    m.doc() = "A cpp based mertable data structure"; 

    py::class_<TableMer>(m, "TableMer")
        .def(py::init<size_t>())
        .def("add", &TableMer::add, "add a new key-value(pair) pair to the merTable", 
             py::arg("dna"), py::arg("bunch_id"), py::arg("loci"), py::arg("local_count"))
        .def("set_k", &TableMer::set_k, "Set the k-mer size", py::arg("k"))
        .def("set_min_entropy", &TableMer::set_min_entropy, "Set the minimum entropy threshold", py::arg("min_entropy"))
        .def("get", &TableMer::get, "get a random pair from the merTable", 
             py::arg("dna"))
        .def("dump", &TableMer::dump, "dump the merTable to a file", 
             py::arg("filename")); 

    m.def("mer_entropy", &mer_entropy, "Calculate the entropy of a DNA sequence", 
          py::arg("seq"));
}