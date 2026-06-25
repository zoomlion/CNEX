#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <cstdint>
#include <vector>
#include <string>

#include <cmath>

namespace py = pybind11;

#ifdef __GNUC__
#define POPCOUNT64(x) __builtin_popcountll(x)
#else
inline uint64_t popcount64_fallback(uint64_t x) {
    x = x - ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x + (x >> 4)) & 0x0F0F0F0F0F0F0F0FULL;
    return (x * 0x0101010101010101ULL) >> 56;
}
#define POPCOUNT64(x) popcount64_fallback(x)
#endif

py::bytes build_bitmap(const std::vector<uint32_t>& codes, int max_code) {
    size_t n_words = (static_cast<size_t>(max_code) + 63) / 64;
    std::string data(n_words * 8, '\0');
    auto* words = reinterpret_cast<uint64_t*>(data.data());
    for (uint32_t code : codes) {
        size_t word_idx = static_cast<size_t>(code) >> 6;
        size_t bit_pos  = static_cast<size_t>(code) & 63;
        words[word_idx] |= (1ULL << bit_pos);
    }
    return py::bytes(data);
}

uint64_t popcount(py::bytes bitmap) {
    std::string data = bitmap;
    auto* words = reinterpret_cast<const uint64_t*>(data.data());
    size_t n = data.size() / 8;
    uint64_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        count += POPCOUNT64(words[i]);
    }
    return count;
}

double jaccard_distance(py::bytes a, py::bytes b) {
    std::string da = a, db = b;
    auto* wa = reinterpret_cast<const uint64_t*>(da.data());
    auto* wb = reinterpret_cast<const uint64_t*>(db.data());
    size_t n = da.size() / 8;
    uint64_t inter = 0, union_cnt = 0;
    for (size_t i = 0; i < n; ++i) {
        inter += POPCOUNT64(wa[i] & wb[i]);
        union_cnt += POPCOUNT64(wa[i] | wb[i]);
    }
    if (union_cnt == 0) return 1.0;
    return 1.0 - static_cast<double>(inter) / static_cast<double>(union_cnt);
}

double bray_curtis_distance(py::bytes a, py::bytes b) {
    std::string da = a, db = b;
    auto* wa = reinterpret_cast<const uint64_t*>(da.data());
    auto* wb = reinterpret_cast<const uint64_t*>(db.data());
    size_t n = da.size() / 8;
    uint64_t inter = 0, cnt_a = 0, cnt_b = 0;
    for (size_t i = 0; i < n; ++i) {
        inter += POPCOUNT64(wa[i] & wb[i]);
        cnt_a += POPCOUNT64(wa[i]);
        cnt_b += POPCOUNT64(wb[i]);
    }
    uint64_t sum = cnt_a + cnt_b;
    if (sum == 0) return 1.0;
    return 1.0 - 2.0 * static_cast<double>(inter) / static_cast<double>(sum);
}

double mash_distance(py::bytes a, py::bytes b, int k) {
    std::string da = a, db = b;
    auto* wa = reinterpret_cast<const uint64_t*>(da.data());
    auto* wb = reinterpret_cast<const uint64_t*>(db.data());
    size_t n = da.size() / 8;
    uint64_t inter = 0, cnt_a = 0, cnt_b = 0;
    for (size_t i = 0; i < n; ++i) {
        inter += POPCOUNT64(wa[i] & wb[i]);
        cnt_a += POPCOUNT64(wa[i]);
        cnt_b += POPCOUNT64(wb[i]);
    }
    if (inter == 0) return 1.0;
    double J = static_cast<double>(inter) / static_cast<double>(cnt_a + cnt_b - inter);
    double ratio = 2.0 * J / (1.0 + J);
    return -std::log(ratio) / static_cast<double>(k);
}

double containment_distance(py::bytes a, py::bytes b, int k) {
    std::string da = a, db = b;
    auto* wa = reinterpret_cast<const uint64_t*>(da.data());
    auto* wb = reinterpret_cast<const uint64_t*>(db.data());
    size_t n = da.size() / 8;
    uint64_t inter = 0, cnt_a = 0, cnt_b = 0;
    for (size_t i = 0; i < n; ++i) {
        inter += POPCOUNT64(wa[i] & wb[i]);
        cnt_a += POPCOUNT64(wa[i]);
        cnt_b += POPCOUNT64(wb[i]);
    }
    uint64_t min_cnt = cnt_a < cnt_b ? cnt_a : cnt_b;
    if (min_cnt == 0 || inter == 0) return 1.0;
    double C = static_cast<double>(inter) / static_cast<double>(min_cnt);
    return -std::log(C) / static_cast<double>(k);
}

py::bytes build_bitmap_from_contigs(
    const std::vector<std::vector<uint32_t>>& contig_codes,
    const std::vector<int>& indices,
    int max_code)
{
    size_t n_words = (static_cast<size_t>(max_code) + 63) / 64;
    std::string data(n_words * 8, '\0');
    auto* words = reinterpret_cast<uint64_t*>(data.data());
    for (int idx : indices) {
        for (uint32_t code : contig_codes[static_cast<size_t>(idx)]) {
            size_t word_idx = static_cast<size_t>(code) >> 6;
            size_t bit_pos  = static_cast<size_t>(code) & 63;
            words[word_idx] |= (1ULL << bit_pos);
        }
    }
    return py::bytes(data);
}

PYBIND11_MODULE(popcount, m) {
    m.doc() = "Fast bitmap popcount & Jaccard for align-free phylogenetics";
    m.def("build_bitmap", &build_bitmap,
          py::arg("codes"), py::arg("max_code"));
    m.def("build_bitmap_from_contigs", &build_bitmap_from_contigs,
          py::arg("contig_codes"), py::arg("indices"), py::arg("max_code"));
    m.def("popcount", &popcount,
          py::arg("bitmap"));
    m.def("jaccard_distance", &jaccard_distance,
          py::arg("a"), py::arg("b"));
    m.def("bray_curtis_distance", &bray_curtis_distance,
          py::arg("a"), py::arg("b"));
    m.def("mash_distance", &mash_distance,
          py::arg("a"), py::arg("b"), py::arg("k"));
    m.def("containment_distance", &containment_distance,
          py::arg("a"), py::arg("b"), py::arg("k"));
}
