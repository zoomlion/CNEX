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
using KmerCounts = tsl::robin_map<std::string, int>;

using EdgeCounts = tsl::robin_map<std::string, int>;  // "src\ttgt" → read count

inline std::tuple<Graph, KmerCounts, EdgeCounts> de_bruijn_graph(
    const std::vector<std::string>& reads, uint64_t k, int min_count = 1
) {
    tsl::robin_map<std::string_view, int> sv_raw;
    KmerCounts k_counts;
    EdgeCounts e_counts;
    sv_raw.reserve(reads.size() * 10);
    k_counts.reserve(reads.size() * 2);
    e_counts.reserve(reads.size());

    for (const auto& read : reads) {
        if (read.length() < (k + 1)) continue;
        auto sv = std::string_view(read);
        for (size_t i = 0; i <= read.length() - k; ++i)
            k_counts[std::string(sv.substr(i, k))]++;
        for (size_t i = 0; i <= read.length() - (k + 1); ++i)
            sv_raw[sv.substr(i, k + 1)]++;
    }

    std::set<std::pair<std::string_view, std::string_view>> valid_edges;
    for (const auto& [kmer, count] : sv_raw) {
        if (count < min_count) continue;
        auto prefix = kmer.substr(0, k);
        auto suffix = kmer.substr(1, k);
        valid_edges.insert({prefix, suffix});
    }

    // Build graph and count edge traversals (per read, dedup)
    Graph graph;
    graph.reserve(valid_edges.size());
    for (const auto& edge : valid_edges) {
        graph[std::string(edge.first)].push_back(std::string(edge.second));
    }

    // Count edge traversals: each read counts each edge at most once
    for (const auto& read : reads) {
        if (read.length() < (k + 1)) continue;
        auto sv = std::string_view(read);
        std::set<std::string> seen_edges;
        for (size_t i = 0; i <= read.length() - (k + 1); ++i) {
            auto kmer = sv.substr(i, k + 1);
            auto prefix = std::string(kmer.substr(0, k));
            auto suffix = std::string(kmer.substr(1, k));
            auto git = graph.find(prefix);
            if (git == graph.end()) continue;
            bool has_edge = false;
            for (const auto& tgt : git->second) {
                if (tgt == suffix) { has_edge = true; break; }
            }
            if (has_edge) {
                std::string ekey = prefix + "\t" + suffix;
                if (seen_edges.find(ekey) == seen_edges.end()) {
                    seen_edges.insert(ekey);
                    e_counts[ekey]++;
                }
            }
        }
    }

    return {std::move(graph), std::move(k_counts), std::move(e_counts)};
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
                               const std::unordered_map<std::string, int>& node_scores,
                               std::vector<std::string>* out_path = nullptr) {
    if (graph.empty()) {
        return "";
    }

    auto [start_node, end_node] = find_farthest_start_end_nodes(graph);
    const std::unordered_map<std::string, int>* ns_ptr = node_scores.empty() ? nullptr : &node_scores;
    auto path = trace_longest_path(graph, start_node, end_node, 2500, ns_ptr);

    if (path.empty()) {
        return "";
    }

    if (out_path) *out_path = path;

    std::string result = path[0];
    for (size_t i = 1; i < path.size(); ++i) {
        result += path[i].back();
    }

    return result;
}


// ─── Export path GFA for a single element (only contig path + branches) ───

inline void export_path_gfa(
    const Graph& graph,
    const std::unordered_map<std::string, int>& node_scores,
    const std::vector<std::string>& contig_path,
    int ele_id, int k,
    std::ofstream& out,
    const KmerCounts* kmer_counts = nullptr)
{
    if (contig_path.empty()) return;

    // Collect nodes in contig path
    tsl::robin_map<std::string, int> in_path;
    for (size_t i = 0; i < contig_path.size(); ++i)
        in_path[contig_path[i]] = 1;

    // Collect all nodes: path + immediate branches
    tsl::robin_map<std::string, int> used_set = in_path;
    for (const auto& node : contig_path) {
        auto git = graph.find(node);
        if (git == graph.end()) continue;
        for (const auto& next : git->second) {
            if (in_path.find(next) == in_path.end())
                used_set[next] = 1;
        }
    }

    // Assign IDs: {ele_id}_{n}
    tsl::robin_map<std::string, int> node_id;
    int n = 0;
    for (const auto& [node, _] : used_set)
        node_id[node] = n++;

    // S lines
    for (const auto& [node, _] : used_set) {
        int kc = 0;
        if (kmer_counts) {
            auto kit = kmer_counts->find(node);
            if (kit != kmer_counts->end()) kc = kit->second;
        }
        auto sit = node_scores.find(node);
        int ms = (sit != node_scores.end()) ? sit->second : 0;
        out << "S\t" << ele_id << "_" << node_id[node]
            << "\t" << node
            << "\tKC:i:" << kc
            << "\tMS:i:" << ms
            << "\tEL:Z:" << ele_id << "\n";
    }

    // L lines (only where both nodes are used)
    for (const auto& [node, _] : used_set) {
        auto git = graph.find(node);
        if (git == graph.end()) continue;
        int src = node_id[node];
        for (const auto& next : git->second) {
            auto nit = node_id.find(next);
            if (nit == node_id.end()) continue;
            out << "L\t" << ele_id << "_" << src << "\t+\t"
                << ele_id << "_" << nit->second << "\t+\t"
                << (k - 1) << "M\n";
        }
    }

    // P line: contig path
    out << "P\tcontig_" << ele_id;
    for (const auto& node : contig_path) {
        auto it = node_id.find(node);
        if (it != node_id.end())
            out << "\t" << ele_id << "_" << it->second;
    }
    out << "\t" << std::string(contig_path.size(), '+') << "\n";
}


// ─── SNP scan (De Bruijn bubble detection) ───

inline bool scan_snps(  // returns true if any variant was found
    const Graph& graph,
    const EdgeCounts& e_counts,
    const std::vector<std::string>& contig_path,
    const std::string& contig_seq,
    int ele_id, int k,
    std::ofstream& out,
    double min_af = 0.05,
    int max_branch_depth = 5)
{
    bool has_var = false;
    if (contig_path.empty()) return false;

    // Build position lookup: contig_path node → its index in path
    tsl::robin_map<std::string, int> path_pos;
    for (size_t i = 0; i < contig_path.size(); ++i)
        path_pos[contig_path[i]] = static_cast<int>(i);

    for (size_t i = 0; i + 1 < contig_path.size(); ++i) {
        const std::string& node = contig_path[i];
        const std::string& next_main = contig_path[i + 1];

        auto git = graph.find(node);
        if (git == graph.end() || git->second.size() < 2) continue;

        for (const auto& branch : git->second) {
            if (branch == next_main) continue;  // main path, not a branch

            // BFS from branch node: look for reconnection to main path
            std::vector<std::string> frontier = {branch};
            tsl::robin_map<std::string, int> visited;
            visited[branch] = 1;
            std::string reconnect_node;
            int reconnect_pos = -1;

            for (int depth = 1; depth <= max_branch_depth && !frontier.empty(); ++depth) {
                std::vector<std::string> next_frontier;
                for (const auto& cur : frontier) {
                    auto cgit = graph.find(cur);
                    if (cgit == graph.end()) continue;
                    for (const auto& cn : cgit->second) {
                        auto pp = path_pos.find(cn);
                        if (pp != path_pos.end() && pp->second > static_cast<int>(i)) {
                            // Found reconnection
                            reconnect_node = cn;
                            reconnect_pos = pp->second;
                            break;
                        }
                        if (visited.find(cn) == visited.end()) {
                            visited[cn] = depth + 1;
                            next_frontier.push_back(cn);
                        }
                    }
                    if (reconnect_pos >= 0) break;
                }
                frontier = std::move(next_frontier);
                if (reconnect_pos >= 0) break;
            }

            if (reconnect_pos < 0) continue;

            int branch_len = reconnect_pos - static_cast<int>(i);
            if (branch_len > 10) continue;

            // Edge count: how many reads traverse this specific edge
            std::string main_key = node + "\t" + next_main;
            std::string alt_key  = node + "\t" + branch;

            auto rit = e_counts.find(main_key);
            auto ait = e_counts.find(alt_key);
            int ref_cov = (rit != e_counts.end()) ? rit->second : 0;
            int alt_cov = (ait != e_counts.end()) ? ait->second : 0;
            if (ref_cov < 2 || alt_cov < 2) continue;

            double af = static_cast<double>(alt_cov) / (ref_cov + alt_cov);
            if (af < min_af) continue;

            int contig_pos = static_cast<int>(i + k - 1);
            const char* type = (branch_len == 1) ? "SNP" : "INDEL";
            char ref_base = next_main.back();
            char alt_base = branch.back();

            out << ele_id << "\t" << contig_pos
                << "\t" << ref_base << "\t" << alt_base
                << "\t" << ref_cov << "\t" << alt_cov
                << "\t" << af << "\t" << branch_len
                << "\t" << type << "\n";
            has_var = true;
        }
    }
    return has_var;
}
