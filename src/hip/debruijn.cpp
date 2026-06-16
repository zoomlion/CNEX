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
#include <string_view>
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <cstdint>
#include <climits>
#include <algorithm>

namespace py = pybind11;

// Graph is std::unordered_map for pybind11 compatibility
// Internal computations use tsl::robin_map for speed
using Graph = std::unordered_map<std::string, std::vector<std::string>>;

Graph de_bruijn_graph(
    const std::vector<std::string>& reads, uint64_t k, int min_count = 1
) {
    // Count (k+1)-mer occurrences to filter rare k-mers (likely sequencing errors)
    tsl::robin_map<std::string_view, int> kmer_counts;
    kmer_counts.reserve(reads.size() * 10);
    for (const auto& read : reads) {
        if (read.length() < (k + 1)) continue;
        auto sv = std::string_view(read);
        for (size_t i = 0; i <= read.length() - (k + 1); ++i) {
            kmer_counts[sv.substr(i, k + 1)]++;
        }
    }

    // Build edges only from (k+1)-mers meeting minimum count threshold
    std::set<std::pair<std::string_view, std::string_view>> valid_edges;
    for (const auto& [kmer, count] : kmer_counts) {
        if (count < min_count) continue;
        auto prefix = kmer.substr(0, k);
        auto suffix = kmer.substr(1, k);
        valid_edges.insert({prefix, suffix});
    }

    Graph graph;
    graph.reserve(valid_edges.size());
    for (const auto& edge : valid_edges) {
        graph[std::string(edge.first)].push_back(std::string(edge.second));
    }
    return graph;
}


tsl::robin_map<std::string, int> bfs_distances(
    const Graph& graph, const std::string& start
) {
    tsl::robin_map<std::string, int> distances;
    std::queue<std::string> q;
    distances[start] = 0;
    q.push(start);

    while (!q.empty()) {
        auto current = q.front(); q.pop();
        int d = distances[current];
        auto it = graph.find(current);
        if (it == graph.end()) continue;
        for (const auto& neighbor : it->second) {
            if (distances.find(neighbor) == distances.end()) {
                distances[neighbor] = d + 1;
                q.push(neighbor);
            }
        }
    }
    return distances;
}


std::pair<std::string, std::string> find_farthest_start_end_nodes(
    const Graph& graph
) {
    tsl::robin_map<std::string, int> in_degree, out_degree;

    for (const auto& [node, neighbors] : graph) {
        out_degree[node] = static_cast<int>(neighbors.size());
        for (const auto& neighbor : neighbors) {
            in_degree[neighbor]++;
        }
    }

    // Collect top-3 start and end candidates by degree difference
    std::vector<std::pair<std::string, int>> start_scores, end_scores;
    for (const auto& [node, _] : graph) {
        int diff = out_degree[node] - in_degree[node];
        if (diff > 0) start_scores.push_back({node, diff});
        else if (diff < 0) end_scores.push_back({node, -diff});
    }
    std::sort(start_scores.begin(), start_scores.end(),
        [](auto& a, auto& b) { return a.second > b.second; });
    std::sort(end_scores.begin(), end_scores.end(),
        [](auto& a, auto& b) { return a.second > b.second; });

    std::vector<std::string> start_cands, end_cands;
    for (int i = 0; i < std::min(3, (int)start_scores.size()); i++)
        start_cands.push_back(start_scores[i].first);
    for (int i = 0; i < std::min(3, (int)end_scores.size()); i++)
        end_cands.push_back(end_scores[i].first);

    if (start_cands.empty() || end_cands.empty()) {
        for (const auto& [node, _] : graph) {
            return {node, node};
        }
    }

    // BFS from each top start candidate to find farthest end
    std::string best_start, best_end;
    int max_dist = -1;

    for (const auto& start : start_cands) {
        auto dists = bfs_distances(graph, start);
        for (const auto& end : end_cands) {
            auto it = dists.find(end);
            if (it != dists.end() && it->second > max_dist) {
                max_dist = it->second;
                best_start = start;
                best_end = end;
            }
        }
    }

    if (best_start.empty() || best_end.empty()) {
        return {start_cands[0], end_cands[0]};
    }
    return {best_start, best_end};
}


// Greedy trace from start following highest start-distance neighbors
std::vector<std::string> greedy_trace(
    const Graph& graph,
    const std::string& start,
    const std::string& end,
    const tsl::robin_map<std::string, int>& start_dists,
    uint64_t max_path_length
) {
    std::vector<std::string> path;
    std::string current = start;

    while (true) {
        path.push_back(current);
        if (path.size() > max_path_length) break;
        if (current == end) break;

        auto it = graph.find(current);
        if (it == graph.end() || it->second.empty()) break;

        std::string best_next;
        int best_dist = -1;
        for (const auto& neighbor : it->second) {
            auto dit = start_dists.find(neighbor);
            if (dit != start_dists.end() && static_cast<int>(dit->second) > best_dist) {
                best_dist = static_cast<int>(dit->second);
                best_next = neighbor;
            }
        }
        if (best_next.empty()) break;
        current = best_next;
    }
    return path;
}

// Extend path forward from last node (follow most-covered branch)
void extend_forward(
    const Graph& graph,
    std::vector<std::string>& path,
    uint64_t max_path_length
) {
    tsl::robin_map<std::string, int> visited;
    for (const auto& n : path) visited[n] = 1;

    std::string current = path.back();
    while (path.size() < max_path_length) {
        auto it = graph.find(current);
        if (it == graph.end() || it->second.empty()) break;

        std::string best_next;
        int max_out = 0;
        for (const auto& neighbor : it->second) {
            if (visited.find(neighbor) != visited.end()) continue;
            auto nit = graph.find(neighbor);
            int out_deg = (nit != graph.end()) ? static_cast<int>(nit->second.size()) : 0;
            if (out_deg > max_out) {
                max_out = out_deg;
                best_next = neighbor;
            }
        }
        if (best_next.empty()) break;
        visited[best_next] = 1;
        path.push_back(best_next);
        current = best_next;
    }
}

// Build reverse graph for backward extension
tsl::robin_map<std::string, std::vector<std::string>> build_reverse_graph(const Graph& graph) {
    tsl::robin_map<std::string, std::vector<std::string>> rev;
    for (const auto& [node, neighbors] : graph) {
        for (const auto& neighbor : neighbors) {
            rev[neighbor].push_back(node);
        }
    }
    return rev;
}

// Extend path backward from first node (prepend nodes)
void extend_backward(
    const Graph& graph,
    std::vector<std::string>& path,
    uint64_t max_path_length
) {
    auto rev = build_reverse_graph(graph);
    tsl::robin_map<std::string, int> visited;
    for (const auto& n : path) visited[n] = 1;

    std::string current = path.front();
    while (path.size() < max_path_length) {
        auto it = rev.find(current);
        if (it == rev.end() || it->second.empty()) break;

        std::string best_prev;
        int max_in = 0;
        for (const auto& prev : it->second) {
            if (visited.find(prev) != visited.end()) continue;
            auto pit = rev.find(prev);
            int in_deg = (pit != rev.end()) ? static_cast<int>(pit->second.size()) : 0;
            if (in_deg > max_in) {
                max_in = in_deg;
                best_prev = prev;
            }
        }
        if (best_prev.empty()) break;
        visited[best_prev] = 1;
        path.insert(path.begin(), best_prev);
        current = best_prev;
    }
}


std::vector<std::string> trace_longest_path(
    const Graph& graph,
    const std::string& start,
    const std::string& end,
    uint64_t max_path_length = 2500
) {
    auto start_dists = bfs_distances(graph, start);

    auto path = greedy_trace(graph, start, end, start_dists, max_path_length);

    // Extend forward/backward to capture more sequence
    if (path.size() < max_path_length) {
        extend_forward(graph, path, max_path_length);
    }
    if (path.size() < max_path_length) {
        extend_backward(graph, path, max_path_length);
    }

    return path;
}


std::string assemble_sequence(const Graph& graph) {
    if (graph.empty()) {
        return "";
    }

    auto [start_node, end_node] = find_farthest_start_end_nodes(graph);
    auto path = trace_longest_path(graph, start_node, end_node);

    if (path.empty()) {
        return "";
    }

    std::string result = path[0];
    for (size_t i = 1; i < path.size(); ++i) {
        result += path[i].back();
    }

    return result;
}


void output_gml(const Graph& graph, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::ios_base::failure("Failed to open file: " + filename);
    }

    file << "graph [\n";

    int node_id = 0;
    tsl::robin_map<std::string, int> node_map;

    for (const auto& [node, neighbors] : graph) {
        if (node_map.find(node) == node_map.end()) {
            node_map[node] = node_id++;
            file << "  node [\n    id " << node_map[node] << "\n    label \"" << node << "\"\n  ]\n";
        }
        for (const auto& neighbor : neighbors) {
            if (node_map.find(neighbor) == node_map.end()) {
                node_map[neighbor] = node_id++;
                file << "  node [\n    id " << node_map[neighbor] << "\n    label \"" << neighbor << "\"\n  ]\n";
            }
        }
    }

    for (const auto& [node, neighbors] : graph) {
        for (const auto& neighbor : neighbors) {
            file << "  edge [\n    source " << node_map[node] << "\n    target " << node_map[neighbor] << "\n  ]\n";
        }
    }

    file << "]\n";
}


PYBIND11_MODULE(debruijn, m) {
    m.doc() = "De Bruijn graph and sequence assembler";

    m.def("de_bruijn_graph", &de_bruijn_graph, "Build a De Bruijn graph from reads",
          py::arg("reads"), py::arg("k"), py::arg("min_count") = 1);
    m.def("assemble_sequence", &assemble_sequence, "Assemble sequence from De Bruijn graph",
          py::arg("graph"));
    m.def("output_gml", &output_gml, "Output De Bruijn graph as GML file",
          py::arg("graph"), py::arg("filename"));
}
