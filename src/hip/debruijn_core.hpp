#pragma once
#include <tsl/robin_map.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

// Graph type: node -> list of successor nodes
using Graph = std::unordered_map<std::string, std::vector<std::string>>;

inline Graph de_bruijn_graph(
    const std::vector<std::string>& reads, uint64_t k, int min_count = 1
) {
    tsl::robin_map<std::string_view, int> kmer_counts;
    kmer_counts.reserve(reads.size() * 10);
    for (const auto& read : reads) {
        if (read.length() < (k + 1)) continue;
        auto sv = std::string_view(read);
        for (size_t i = 0; i <= read.length() - (k + 1); ++i) {
            kmer_counts[sv.substr(i, k + 1)]++;
        }
    }

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


inline tsl::robin_map<std::string, int> bfs_distances(
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


inline std::pair<std::string, std::string> find_farthest_start_end_nodes(
    const Graph& graph
) {
    tsl::robin_map<std::string, int> in_degree, out_degree;

    for (const auto& [node, neighbors] : graph) {
        out_degree[node] = static_cast<int>(neighbors.size());
        for (const auto& neighbor : neighbors) {
            in_degree[neighbor]++;
        }
    }

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


inline std::vector<std::string> greedy_trace(
    const Graph& graph,
    const std::string& start,
    const std::string& end,
    const tsl::robin_map<std::string, int>& start_dists,
    uint64_t max_path_length,
    const std::unordered_map<std::string, int>* node_scores = nullptr
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
        int best_score = -1;
        for (const auto& neighbor : it->second) {
            auto dit = start_dists.find(neighbor);
            if (dit == start_dists.end()) continue;
            int dist = static_cast<int>(dit->second);
            int score = 0;
            if (node_scores) {
                auto sit = node_scores->find(neighbor);
                if (sit != node_scores->end()) score = sit->second;
            }
            if (dist > best_dist || (dist == best_dist && score > best_score)) {
                best_dist = dist;
                best_score = score;
                best_next = neighbor;
            }
        }
        if (best_next.empty()) break;
        current = best_next;
    }
    return path;
}


inline void extend_forward(
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


inline tsl::robin_map<std::string, std::vector<std::string>> build_reverse_graph(const Graph& graph) {
    tsl::robin_map<std::string, std::vector<std::string>> rev;
    for (const auto& [node, neighbors] : graph) {
        for (const auto& neighbor : neighbors) {
            rev[neighbor].push_back(node);
        }
    }
    return rev;
}


inline void extend_backward(
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


inline std::vector<std::string> trace_longest_path(
    const Graph& graph,
    const std::string& start,
    const std::string& end,
    uint64_t max_path_length = 2500,
    const std::unordered_map<std::string, int>* node_scores = nullptr
) {
    auto start_dists = bfs_distances(graph, start);
    auto path = greedy_trace(graph, start, end, start_dists, max_path_length, node_scores);

    if (path.size() < max_path_length) {
        extend_forward(graph, path, max_path_length);
    }
    if (path.size() < max_path_length) {
        extend_backward(graph, path, max_path_length);
    }

    return path;
}


inline std::string assemble_sequence(const Graph& graph,
                               const std::unordered_map<std::string, int>& node_scores) {
    if (graph.empty()) {
        return "";
    }

    auto [start_node, end_node] = find_farthest_start_end_nodes(graph);
    const std::unordered_map<std::string, int>* ns_ptr = node_scores.empty() ? nullptr : &node_scores;
    auto path = trace_longest_path(graph, start_node, end_node, 2500, ns_ptr);

    if (path.empty()) {
        return "";
    }

    std::string result = path[0];
    for (size_t i = 1; i < path.size(); ++i) {
        result += path[i].back();
    }

    return result;
}
