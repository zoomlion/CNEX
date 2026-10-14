// panel_core.hpp - panel loading for the in-validate T3 rescue.
//
// Loads the top-k longest non-ref sequences per element (MSA block) from the
// same leave-one-out block file used to build the mers table. The panel is
// consumed by the SIMD Smith-Waterman verification of near-miss windows/reads.

#pragma once

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace panel_core {

inline char upper(char c) {
    return (c >= 'a' && c <= 'z') ? static_cast<char>(c - 'a' + 'A') : c;
}

// All elements are loaded (the caller cannot know in advance which windows
// will reach T3); memory is roughly k * avg_panel_len * n_elements.
inline void load_panels(const std::string& msa, int panel_k,
                        std::unordered_map<int, std::vector<std::string>>& panel) {
    std::ifstream f(msa);
    if (!f) throw std::runtime_error("Cannot open panel MSA: " + msa);
    std::string line;
    std::vector<std::string> cur;
    bool in_block = false;
    int bunch_id = 0;
    auto flush = [&]() {
        if (in_block && !cur.empty()) {
            auto& v = panel[bunch_id];
            v = std::move(cur);
            std::sort(v.begin(), v.end(),
                      [](const std::string& a, const std::string& b) {
                          return a.size() > b.size();
                      });
            if (static_cast<int>(v.size()) > panel_k) v.resize(panel_k);
        }
        ++bunch_id;
        cur.clear();
        in_block = false;
    };
    std::string hdr;
    std::string seq;
    bool wanted = false;
    auto flush_seq = [&]() {
        if (wanted && !hdr.empty() && !seq.empty() && hdr.rfind("ref", 0) != 0) {
            std::string pure;
            pure.reserve(seq.size());
            for (char c : seq) {
                if (c != '-' && c != ' ') pure.push_back(upper(c));
            }
            if (!pure.empty()) cur.push_back(std::move(pure));
        }
        hdr.clear();
        seq.clear();
        wanted = false;
    };
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("###", 0) == 0) {
            flush_seq();
            flush();
            continue;
        }
        if (!line.empty() && line[0] == '>') {
            flush_seq();
            in_block = true;
            hdr = line.substr(1);
            wanted = true;
        } else if (in_block) {
            seq += line;
        }
    }
    flush_seq();
    flush();
}

}  // namespace panel_core
