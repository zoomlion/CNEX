#pragma once
#include <tsl/robin_map.h>
#include <tsl/robin_set.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
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
    const std::vector<std::string>& reads, uint64_t k, int min_count = 1,
    bool need_kcounts = false, bool need_ecounts = false
) {
    tsl::robin_map<std::string_view, int> sv_raw;
    KmerCounts k_counts;
    EdgeCounts e_counts;
    sv_raw.reserve(reads.size() * 10);
    if (need_kcounts) k_counts.reserve(reads.size() * 2);
    if (need_ecounts) e_counts.reserve(reads.size());

    for (const auto& read : reads) {
        if (read.length() < (k + 1)) continue;
        auto sv = std::string_view(read);
        if (need_kcounts) {
            for (size_t i = 0; i <= read.length() - k; ++i)
                k_counts[std::string(sv.substr(i, k))]++;
        }
        for (size_t i = 0; i <= read.length() - (k + 1); ++i)
            sv_raw[sv.substr(i, k + 1)]++;
    }

    std::vector<std::pair<std::string_view, std::string_view>> valid_edges;
    valid_edges.reserve(sv_raw.size());
    for (const auto& [kmer, count] : sv_raw) {
        if (count < min_count) continue;
        valid_edges.emplace_back(kmer.substr(0, k), kmer.substr(1, k));
    }
    // Same ordering as std::set over the (prefix, suffix) pairs.
    std::sort(valid_edges.begin(), valid_edges.end());
    valid_edges.erase(std::unique(valid_edges.begin(), valid_edges.end()),
                      valid_edges.end());

    // Build graph and count edge traversals (per read, dedup)
    Graph graph;
    graph.reserve(valid_edges.size());
    for (const auto& edge : valid_edges) {
        graph[std::string(edge.first)].push_back(std::string(edge.second));
    }

    // Count edge traversals: each read counts each edge at most once.
    // Only needed for SNP calling; skipped when not requested.
    if (need_ecounts) {
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
    tsl::robin_set<std::string> visited;
    visited.insert(current);

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
            if (visited.count(neighbor)) continue;
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
        visited.insert(current);
    }
    return path;
}


// Extend past the tail of `path`. Neighbors are ranked by element score first,
// out-degree only as a tie-break.
inline void extend_forward(
    const Graph& graph,
    std::vector<std::string>& path,
    uint64_t max_path_length,
    const std::unordered_map<std::string, int>* node_scores
) {
    tsl::robin_map<std::string, int> visited;
    for (const auto& n : path) visited[n] = 1;

    std::string current = path.back();
    while (path.size() < max_path_length) {
        auto it = graph.find(current);
        if (it == graph.end() || it->second.empty()) break;

        std::string best_next;
        int best_score = -1;
        int best_out = 0;
        for (const auto& neighbor : it->second) {
            if (visited.find(neighbor) != visited.end()) continue;
            int score = 0;
            if (node_scores) {
                auto sit = node_scores->find(neighbor);
                if (sit != node_scores->end()) score = sit->second;
            }
            auto nit = graph.find(neighbor);
            int out_deg = (nit != graph.end()) ? static_cast<int>(nit->second.size()) : 0;
            if (score > best_score || (score == best_score && out_deg > best_out)) {
                best_score = score;
                best_out = out_deg;
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


// Mirror of extend_forward for the head of `path`.
inline void extend_backward(
    const Graph& graph,
    std::vector<std::string>& path,
    uint64_t max_path_length,
    const std::unordered_map<std::string, int>* node_scores
) {
    auto rev = build_reverse_graph(graph);
    tsl::robin_map<std::string, int> visited;
    for (const auto& n : path) visited[n] = 1;

    std::vector<std::string> bwd;
    std::string current = path.front();
    while (path.size() + bwd.size() < max_path_length) {
        auto it = rev.find(current);
        if (it == rev.end() || it->second.empty()) break;

        std::string best_prev;
        int best_score = -1;
        int best_in = 0;
        for (const auto& prev : it->second) {
            if (visited.find(prev) != visited.end()) continue;
            int score = 0;
            if (node_scores) {
                auto sit = node_scores->find(prev);
                if (sit != node_scores->end()) score = sit->second;
            }
            auto pit = rev.find(prev);
            int in_deg = (pit != rev.end()) ? static_cast<int>(pit->second.size()) : 0;
            if (score > best_score || (score == best_score && in_deg > best_in)) {
                best_score = score;
                best_in = in_deg;
                best_prev = prev;
            }
        }
        if (best_prev.empty()) break;
        visited[best_prev] = 1;
        bwd.push_back(best_prev);
        current = best_prev;
    }
    if (!bwd.empty()) {
        std::reverse(bwd.begin(), bwd.end());
        path.insert(path.begin(),
                    std::make_move_iterator(bwd.begin()),
                    std::make_move_iterator(bwd.end()));
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

    // Extension is evidence-ranked (score first, degree only as tie-break) but
    // deliberately does not stop at the heuristic start/end nodes: those are
    // topology guesses and legitimate element sequence often continues past
    // them.
    if (path.size() < max_path_length) {
        extend_forward(graph, path, max_path_length, node_scores);
    }
    if (path.size() < max_path_length) {
        extend_backward(graph, path, max_path_length, node_scores);
    }

    return path;
}


inline std::string path_to_string(const std::vector<std::string>& path) {
    if (path.empty()) return "";
    std::string s = path[0];
    for (size_t i = 1; i < path.size(); ++i) s += path[i].back();
    return s;
}


// Evidence-anchored walk: seeded at the highest-scoring node, extends both
// directions preferring element-scored neighbors, bridging short unscored gaps.
// Used as a second pass when the topology-driven path carries no element evidence.
inline std::vector<std::string> anchored_trace(
    const Graph& graph,
    const std::unordered_map<std::string, int>& node_scores,
    uint64_t max_path_length = 2500,
    int gap_limit = 3)
{
    if (graph.empty() || node_scores.empty()) return {};

    std::unordered_map<std::string, std::vector<std::string>> rev;
    rev.reserve(graph.size());
    for (const auto& [node, nexts] : graph)
        for (const auto& nx : nexts) rev[nx].push_back(node);

    auto get_score = [&](const std::string& n) {
        auto it = node_scores.find(n);
        return it == node_scores.end() ? 0 : it->second;
    };
    auto degree = [&](const std::string& n) {
        auto it = graph.find(n);
        size_t out = (it == graph.end()) ? 0 : it->second.size();
        auto rit = rev.find(n);
        size_t in = (rit == rev.end()) ? 0 : rit->second.size();
        return out + in;
    };

    std::string seed;
    int best_sc = 0;
    size_t best_deg = 0;
    for (const auto& [node, _] : graph) {
        int sc = get_score(node);
        if (sc <= 0) continue;
        size_t deg = degree(node);
        if (sc > best_sc || (sc == best_sc && deg > best_deg)) {
            best_sc = sc; best_deg = deg; seed = node;
        }
    }
    if (seed.empty()) return {};

    auto walk = [&](const Graph& g, const std::string& start) {
        std::vector<std::string> path;
        tsl::robin_set<std::string> visited;
        std::string cur = start;
        visited.insert(cur);
        int gap = 0;
        while (path.size() + 1 < max_path_length) {
            path.push_back(cur);
            auto it = g.find(cur);
            if (it == g.end() || it->second.empty()) break;
            std::string chosen;
            int chosen_sc = -1;
            size_t chosen_deg = 0;
            for (const auto& nx : it->second) {
                if (visited.count(nx)) continue;
                int sc = get_score(nx);
                if (sc <= 0) continue;
                size_t deg = degree(nx);
                if (sc > chosen_sc || (sc == chosen_sc && deg > chosen_deg)) {
                    chosen = nx; chosen_sc = sc; chosen_deg = deg;
                }
            }
            if (chosen.empty()) {
                if (gap >= gap_limit) break;
                int allowed = gap_limit - gap;
                int best_reach = -1;
                for (const auto& nx : it->second) {
                    if (visited.count(nx)) continue;
                    tsl::robin_set<std::string> seen = visited;
                    seen.insert(nx);
                    std::vector<std::string> frontier{nx};
                    int reach = -1;
                    for (int d = 1; d <= allowed && reach < 0; ++d) {
                        std::vector<std::string> next_frontier;
                        for (const auto& f : frontier) {
                            auto fit = g.find(f);
                            if (fit == g.end()) continue;
                            for (const auto& f2 : fit->second) {
                                if (seen.count(f2)) continue;
                                int sc2 = get_score(f2);
                                if (sc2 > 0) { reach = sc2; break; }
                                seen.insert(f2);
                                next_frontier.push_back(f2);
                            }
                            if (reach >= 0) break;
                        }
                        frontier = std::move(next_frontier);
                    }
                    if (reach > best_reach) { best_reach = reach; chosen = nx; }
                }
                if (chosen.empty()) break;
                ++gap;
            } else {
                gap = 0;
            }
            visited.insert(chosen);
            cur = chosen;
        }
        return path;
    };

    auto fwd = walk(graph, seed);
    auto bwd = walk(rev, seed);

    std::vector<std::string> path;
    for (auto rit = bwd.rbegin(); rit != bwd.rend(); ++rit)
        if (*rit != seed) path.push_back(*rit);
    path.insert(path.end(), fwd.begin(), fwd.end());
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

    // Collect all nodes: path + branches (BFS to merge point, max depth 10)
    // Also record bubble paths for P-line export
    constexpr int MAX_BRANCH_DEPTH = 10;
    tsl::robin_map<std::string, int> used_set = in_path;
    std::vector<std::vector<std::string>> bubble_paths;
    for (const auto& node : contig_path) {
        auto git = graph.find(node);
        if (git == graph.end()) continue;
        for (const auto& next : git->second) {
            if (in_path.find(next) != in_path.end()) continue;
            // BFS with parent tracking to reconstruct bubble path
            struct QItem { std::string node; int depth; };
            std::queue<QItem> q;
            tsl::robin_map<std::string, std::string> parent;
            q.push({next, 0});
            parent[next] = "";
            while (!q.empty()) {
                auto [cur, depth] = q.front(); q.pop();
                used_set[cur] = 1;
                if (depth >= MAX_BRANCH_DEPTH) continue;
                auto cit = graph.find(cur);
                if (cit == graph.end()) continue;
                for (const auto& cn : cit->second) {
                    if (in_path.find(cn) != in_path.end()) {
                        used_set[cn] = 1;
                        // Reconstruct bubble path from branch start to merge
                        std::vector<std::string> bp;
                        std::string p = cur;
                        while (!p.empty()) {
                            bp.push_back(p);
                            auto it = parent.find(p);
                            if (it != parent.end()) p = it->second;
                            else break;
                        }
                        std::reverse(bp.begin(), bp.end());
                        if (bp.size() >= 1) bubble_paths.push_back(bp);
                        continue;
                    }
                    if (parent.find(cn) != parent.end()) continue;
                    parent[cn] = cur;
                    q.push({cn, depth + 1});
                }
            }
        }
    }

    // Assign IDs: {ele_id}_{n}
    tsl::robin_map<std::string, int> node_id;
    int n = 0;
    for (const auto& [node, _] : used_set)
        node_id[node] = n++;

    // S lines: GFA 1.0 format with KC, MS, EL tags
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

    // L lines: GFA 1.0 format
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

    // P line: main path
    out << "P\tcontig_" << ele_id;
    for (const auto& node : contig_path) {
        auto it = node_id.find(node);
        if (it != node_id.end())
            out << "\t" << ele_id << "_" << it->second;
    }
    out << "\t" << std::string(contig_path.size(), '+') << "\n";

    // P lines: bubble paths
    for (size_t bi = 0; bi < bubble_paths.size(); ++bi) {
        out << "P\tcontig_" << ele_id << "_alt" << std::to_string(bi);
        for (const auto& node : bubble_paths[bi]) {
            auto it = node_id.find(node);
            if (it != node_id.end())
                out << "\t" << ele_id << "_" << it->second;
        }
        out << "\t" << std::string(bubble_paths[bi].size(), '+') << "\n";
    }
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


// ─── Indexed integer-ID De Bruijn graph (zero string allocations during traversal) ───

struct IndexedGraph {
    int V = 0;
    std::vector<std::string_view> nodes;   // node id -> k-mer view into caller's reads
    std::vector<int> adj_off;              // size V + 1
    std::vector<int> adj_dst;              // size E (sorted ascending per node)
    std::vector<int> rev_off;              // size V + 1
    std::vector<int> rev_dst;              // size E
    std::vector<int> in_deg;               // size V
    std::vector<int> out_deg;              // size V
    std::vector<int> kmer_counts;          // populated when need_kcounts
    std::vector<int> edge_counts;          // populated when need_ecounts (indexed by adj_dst slot)

    bool empty() const { return V == 0; }

    int find_edge(int u, int v) const {
        int b = adj_off[u], e = adj_off[u + 1];
        for (int ei = b; ei < e; ++ei) {
            if (adj_dst[ei] == v) return ei;
        }
        return -1;
    }
};

inline IndexedGraph build_indexed_graph(
    const std::vector<std::string>& reads, uint64_t k, int min_count = 1,
    bool need_kcounts = false, bool need_ecounts = false
) {
    IndexedGraph g;
    tsl::robin_map<std::string_view, int> sv_raw;
    sv_raw.reserve(reads.size() * 10);

    tsl::robin_map<std::string_view, int> raw_kcounts;
    if (need_kcounts) raw_kcounts.reserve(reads.size() * 2);

    for (const auto& read : reads) {
        if (read.length() < (k + 1)) continue;
        std::string_view sv(read);
        if (need_kcounts) {
            for (size_t i = 0; i <= read.length() - k; ++i) {
                raw_kcounts[sv.substr(i, k)]++;
            }
        }
        for (size_t i = 0; i <= read.length() - (k + 1); ++i) {
            sv_raw[sv.substr(i, k + 1)]++;
        }
    }

    std::vector<std::pair<std::string_view, std::string_view>> valid_edges;
    valid_edges.reserve(sv_raw.size());
    for (const auto& kv : sv_raw) {
        if (kv.second < min_count) continue;
        valid_edges.emplace_back(kv.first.substr(0, k), kv.first.substr(1, k));
    }
    if (valid_edges.empty()) return g;

    std::sort(valid_edges.begin(), valid_edges.end());
    valid_edges.erase(std::unique(valid_edges.begin(), valid_edges.end()),
                      valid_edges.end());

    std::vector<std::string_view> all_nodes;
    all_nodes.reserve(valid_edges.size() * 2);
    for (const auto& e : valid_edges) {
        all_nodes.push_back(e.first);
        all_nodes.push_back(e.second);
    }
    std::sort(all_nodes.begin(), all_nodes.end());
    all_nodes.erase(std::unique(all_nodes.begin(), all_nodes.end()), all_nodes.end());

    g.V = static_cast<int>(all_nodes.size());
    tsl::robin_map<std::string_view, int> id_of;
    id_of.reserve(static_cast<size_t>(g.V) * 2);

    static const bool use_lex_order = (std::getenv("CNEX_LEX_GRAPH") != nullptr);
    if (use_lex_order) {
        g.nodes = std::move(all_nodes);
        for (int i = 0; i < g.V; ++i) {
            id_of.emplace(g.nodes[i], i);
        }
    } else {
        // Assign active source IDs in std::unordered_map<string_view> bucket order
        // so tie-breaks in find_farthest_start_end_nodes / anchored_trace / extend_backward
        // match baseline's Graph without any std::string or per-node vector allocation.
        std::unordered_map<std::string_view, int> ukeys;
        ukeys.reserve(valid_edges.size());
        for (const auto& e : valid_edges) {
            ukeys.try_emplace(e.first, 0);
        }
        g.nodes.reserve(g.V);
        for (const auto& kv : ukeys) {
            int id = static_cast<int>(g.nodes.size());
            g.nodes.push_back(kv.first);
            id_of.emplace(kv.first, id);
        }
        for (const auto& s : all_nodes) {
            if (id_of.find(s) == id_of.end()) {
                int id = static_cast<int>(g.nodes.size());
                g.nodes.push_back(s);
                id_of.emplace(s, id);
            }
        }
    }

    const size_t E = valid_edges.size();
    std::vector<std::pair<int, int>> int_edges;
    int_edges.reserve(E);
    g.out_deg.assign(g.V, 0);
    g.in_deg.assign(g.V, 0);
    for (const auto& e : valid_edges) {
        int u = id_of.find(e.first)->second;
        int v = id_of.find(e.second)->second;
        int_edges.emplace_back(u, v);
        ++g.out_deg[u];
        ++g.in_deg[v];
    }

    g.adj_off.resize(g.V + 1);
    g.rev_off.resize(g.V + 1);
    g.adj_off[0] = 0;
    g.rev_off[0] = 0;
    for (int i = 0; i < g.V; ++i) {
        g.adj_off[i + 1] = g.adj_off[i] + g.out_deg[i];
        g.rev_off[i + 1] = g.rev_off[i] + g.in_deg[i];
    }
    g.adj_dst.resize(E);
    g.rev_dst.resize(E);
    std::vector<int> adj_cur = g.adj_off;
    for (const auto& uv : int_edges) {
        g.adj_dst[adj_cur[uv.first]++] = uv.second;
    }
    std::vector<int> rev_cur = g.rev_off;
    for (int u = 0; u < g.V; ++u) {
        for (int ei = g.adj_off[u]; ei < g.adj_off[u + 1]; ++ei) {
            g.rev_dst[rev_cur[g.adj_dst[ei]]++] = u;
        }
    }

    if (need_kcounts) {
        g.kmer_counts.assign(g.V, 0);
        for (int i = 0; i < g.V; ++i) {
            auto kit = raw_kcounts.find(g.nodes[i]);
            if (kit != raw_kcounts.end()) g.kmer_counts[i] = kit->second;
        }
    }

    if (need_ecounts) {
        g.edge_counts.assign(E, 0);
        std::vector<uint32_t> edge_epoch(E, 0);
        uint32_t epoch = 0;
        for (const auto& read : reads) {
            if (read.length() < (k + 1)) continue;
            if (++epoch == 0) {
                std::fill(edge_epoch.begin(), edge_epoch.end(), 0);
                epoch = 1;
            }
            std::string_view sv(read);
            for (size_t i = 0; i <= read.length() - (k + 1); ++i) {
                auto uit = id_of.find(sv.substr(i, k));
                if (uit == id_of.end()) continue;
                auto vit = id_of.find(sv.substr(i + 1, k));
                if (vit == id_of.end()) continue;
                int ei = g.find_edge(uit->second, vit->second);
                if (ei >= 0 && edge_epoch[ei] != epoch) {
                    edge_epoch[ei] = epoch;
                    ++g.edge_counts[ei];
                }
            }
        }
    }

    return g;
}

inline void bfs_distances_idx(const IndexedGraph& g, int start, std::vector<int>& dist,
                              std::vector<int>& q) {
    dist.assign(g.V, -1);
    q.clear();
    dist[start] = 0;
    q.push_back(start);
    for (size_t qi = 0; qi < q.size(); ++qi) {
        int u = q[qi];
        int d = dist[u] + 1;
        for (int ei = g.adj_off[u]; ei < g.adj_off[u + 1]; ++ei) {
            int v = g.adj_dst[ei];
            if (dist[v] < 0) {
                dist[v] = d;
                q.push_back(v);
            }
        }
    }
}

inline std::pair<int, int> find_farthest_start_end_nodes_idx(const IndexedGraph& g) {
    std::vector<std::pair<int, int>> start_scores, end_scores;
    int first_active = -1;
    for (int u = 0; u < g.V; ++u) {
        if (g.out_deg[u] == 0) continue;   // matches Graph iteration (keys have out_deg > 0)
        if (first_active < 0) first_active = u;
        int diff = g.out_deg[u] - g.in_deg[u];
        if (diff > 0) start_scores.emplace_back(u, diff);
        else if (diff < 0) end_scores.emplace_back(u, -diff);
    }
    auto by_score_desc = [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
        return a.second > b.second;
    };
    std::sort(start_scores.begin(), start_scores.end(), by_score_desc);
    std::sort(end_scores.begin(), end_scores.end(), by_score_desc);

    int n_start = std::min(3, static_cast<int>(start_scores.size()));
    int n_end = std::min(3, static_cast<int>(end_scores.size()));
    if (n_start == 0 || n_end == 0) {
        int fallback = (first_active >= 0) ? first_active : 0;
        return {fallback, fallback};
    }

    int best_start = -1, best_end = -1, max_dist = -1;
    std::vector<int> dist, q;
    q.reserve(g.V);
    for (int si = 0; si < n_start; ++si) {
        int st = start_scores[si].first;
        bfs_distances_idx(g, st, dist, q);
        for (int ei = 0; ei < n_end; ++ei) {
            int en = end_scores[ei].first;
            if (dist[en] > max_dist) {
                max_dist = dist[en];
                best_start = st;
                best_end = en;
            }
        }
    }
    if (best_start < 0 || best_end < 0) {
        return {start_scores[0].first, end_scores[0].first};
    }
    return {best_start, best_end};
}

inline std::vector<int> trace_longest_path_idx(
    const IndexedGraph& g, int start, int end, uint64_t max_path_length,
    const std::vector<int>* node_scores
) {
    std::vector<int> start_dists, q;
    q.reserve(g.V);
    bfs_distances_idx(g, start, start_dists, q);

    std::vector<int> path;
    std::vector<uint8_t> visited(g.V, 0);
    int current = start;
    visited[current] = 1;

    while (true) {
        path.push_back(current);
        if (path.size() > max_path_length) break;
        if (current == end) break;
        int b = g.adj_off[current], e = g.adj_off[current + 1];
        if (b == e) break;

        int best_next = -1, best_dist = -1, best_score = -1;
        for (int ei = b; ei < e; ++ei) {
            int nb = g.adj_dst[ei];
            if (visited[nb]) continue;
            int dist = start_dists[nb];
            if (dist < 0) continue;
            int sc = node_scores ? (*node_scores)[nb] : 0;
            if (dist > best_dist || (dist == best_dist && sc > best_score)) {
                best_dist = dist;
                best_score = sc;
                best_next = nb;
            }
        }
        if (best_next < 0) break;
        current = best_next;
        visited[current] = 1;
    }

    // extend_forward
    current = path.back();
    while (path.size() < max_path_length) {
        int b = g.adj_off[current], e = g.adj_off[current + 1];
        if (b == e) break;
        int best_next = -1, best_score = -1, best_out = 0;
        for (int ei = b; ei < e; ++ei) {
            int nb = g.adj_dst[ei];
            if (visited[nb]) continue;
            int sc = node_scores ? (*node_scores)[nb] : 0;
            int od = g.out_deg[nb];
            if (sc > best_score || (sc == best_score && od > best_out)) {
                best_score = sc;
                best_out = od;
                best_next = nb;
            }
        }
        if (best_next < 0) break;
        visited[best_next] = 1;
        path.push_back(best_next);
        current = best_next;
    }

    // extend_backward
    if (path.size() < max_path_length) {
        std::vector<int> bwd;
        current = path.front();
        while (path.size() + bwd.size() < max_path_length) {
            int b = g.rev_off[current], e = g.rev_off[current + 1];
            if (b == e) break;
            int best_prev = -1, best_score = -1, best_in = 0;
            for (int ei = b; ei < e; ++ei) {
                int pv = g.rev_dst[ei];
                if (visited[pv]) continue;
                int sc = node_scores ? (*node_scores)[pv] : 0;
                int idg = g.in_deg[pv];
                if (sc > best_score || (sc == best_score && idg > best_in)) {
                    best_score = sc;
                    best_in = idg;
                    best_prev = pv;
                }
            }
            if (best_prev < 0) break;
            visited[best_prev] = 1;
            bwd.push_back(best_prev);
            current = best_prev;
        }
        if (!bwd.empty()) {
            std::reverse(bwd.begin(), bwd.end());
            bwd.insert(bwd.end(), path.begin(), path.end());
            path = std::move(bwd);
        }
    }

    return path;
}

inline std::vector<int> anchored_trace_idx(
    const IndexedGraph& g,
    const std::vector<int>& node_scores,
    uint64_t max_path_length = 2500,
    int gap_limit = 3
) {
    if (g.empty()) return {};
    int seed = -1, best_sc = 0, best_deg = 0;
    for (int u = 0; u < g.V; ++u) {
        if (g.out_deg[u] == 0) continue;   // matches Graph iteration
        int sc = node_scores[u];
        if (sc <= 0) continue;
        int deg = g.out_deg[u] + g.in_deg[u];
        if (sc > best_sc || (sc == best_sc && deg > best_deg)) {
            best_sc = sc;
            best_deg = deg;
            seed = u;
        }
    }
    if (seed < 0) return {};

    auto walk = [&](const std::vector<int>& off, const std::vector<int>& dst) {
        std::vector<int> path;
        std::vector<uint8_t> visited(g.V, 0);
        std::vector<uint8_t> seen(g.V, 0);
        std::vector<int> frontier, next_frontier;
        int cur = seed;
        visited[cur] = 1;
        int gap = 0;
        while (path.size() + 1 < max_path_length) {
            path.push_back(cur);
            int b = off[cur], e = off[cur + 1];
            if (b == e) break;
            int chosen = -1, chosen_sc = -1, chosen_deg = 0;
            for (int ei = b; ei < e; ++ei) {
                int nx = dst[ei];
                if (visited[nx]) continue;
                int sc = node_scores[nx];
                if (sc <= 0) continue;
                int deg = g.out_deg[nx] + g.in_deg[nx];
                if (sc > chosen_sc || (sc == chosen_sc && deg > chosen_deg)) {
                    chosen = nx; chosen_sc = sc; chosen_deg = deg;
                }
            }
            if (chosen < 0) {
                if (gap >= gap_limit) break;
                int allowed = gap_limit - gap;
                int best_reach = -1;
                for (int ei = b; ei < e; ++ei) {
                    int nx = dst[ei];
                    if (visited[nx]) continue;
                    seen = visited;
                    seen[nx] = 1;
                    frontier.clear();
                    frontier.push_back(nx);
                    int reach = -1;
                    for (int d = 1; d <= allowed && reach < 0; ++d) {
                        next_frontier.clear();
                        for (int f : frontier) {
                            for (int fe = off[f]; fe < off[f + 1]; ++fe) {
                                int f2 = dst[fe];
                                if (seen[f2]) continue;
                                int sc2 = node_scores[f2];
                                if (sc2 > 0) { reach = sc2; break; }
                                seen[f2] = 1;
                                next_frontier.push_back(f2);
                            }
                            if (reach >= 0) break;
                        }
                        frontier.swap(next_frontier);
                    }
                    if (reach > best_reach) { best_reach = reach; chosen = nx; }
                }
                if (chosen < 0) break;
                ++gap;
            } else {
                gap = 0;
            }
            visited[chosen] = 1;
            cur = chosen;
        }
        return path;
    };

    auto fwd = walk(g.adj_off, g.adj_dst);
    auto bwd = walk(g.rev_off, g.rev_dst);
    std::vector<int> path;
    path.reserve(bwd.size() + fwd.size());
    for (auto rit = bwd.rbegin(); rit != bwd.rend(); ++rit) {
        if (*rit != seed) path.push_back(*rit);
    }
    path.insert(path.end(), fwd.begin(), fwd.end());
    return path;
}

inline std::string path_to_string_idx(const IndexedGraph& g, const std::vector<int>& path) {
    if (path.empty()) return "";
    std::string s;
    s.reserve(g.nodes[path[0]].size() + path.size() - 1);
    s.append(g.nodes[path[0]].data(), g.nodes[path[0]].size());
    for (size_t i = 1; i < path.size(); ++i) {
        s.push_back(g.nodes[path[i]].back());
    }
    return s;
}

inline std::string assemble_sequence_idx(const IndexedGraph& g,
                                         const std::vector<int>& node_scores,
                                         bool has_any_score,
                                         std::vector<int>* out_path = nullptr) {
    if (g.empty()) return "";
    auto [start_id, end_id] = find_farthest_start_end_nodes_idx(g);
    const std::vector<int>* ns_ptr = has_any_score ? &node_scores : nullptr;
    auto path = trace_longest_path_idx(g, start_id, end_id, 2500, ns_ptr);
    if (path.empty()) return "";
    if (out_path) *out_path = path;
    return path_to_string_idx(g, path);
}

inline void export_path_gfa_idx(
    const IndexedGraph& g,
    const std::vector<int>& node_scores,
    const std::vector<int>& contig_path,
    int ele_id, int k,
    std::ofstream& out)
{
    if (contig_path.empty()) return;
    std::vector<uint8_t> in_path(g.V, 0);
    for (int u : contig_path) in_path[u] = 1;

    constexpr int MAX_BRANCH_DEPTH = 10;
    std::vector<uint8_t> used = in_path;
    std::vector<std::vector<int>> bubble_paths;
    std::vector<int> parent(g.V, -1);
    std::vector<uint8_t> seen(g.V, 0);

    for (int u : contig_path) {
        for (int ei = g.adj_off[u]; ei < g.adj_off[u + 1]; ++ei) {
            int next = g.adj_dst[ei];
            if (in_path[next]) continue;
            std::fill(seen.begin(), seen.end(), 0);
            std::queue<std::pair<int, int>> q;
            q.push({next, 0});
            parent[next] = -1;
            seen[next] = 1;
            while (!q.empty()) {
                auto [cur, depth] = q.front(); q.pop();
                used[cur] = 1;
                if (depth >= MAX_BRANCH_DEPTH) continue;
                for (int ce = g.adj_off[cur]; ce < g.adj_off[cur + 1]; ++ce) {
                    int cn = g.adj_dst[ce];
                    if (in_path[cn]) {
                        used[cn] = 1;
                        std::vector<int> bp;
                        int p = cur;
                        while (p >= 0) {
                            bp.push_back(p);
                            p = parent[p];
                        }
                        std::reverse(bp.begin(), bp.end());
                        if (!bp.empty()) bubble_paths.push_back(std::move(bp));
                        continue;
                    }
                    if (seen[cn]) continue;
                    seen[cn] = 1;
                    parent[cn] = cur;
                    q.push({cn, depth + 1});
                }
            }
        }
    }

    std::vector<int> node_id(g.V, -1);
    int n = 0;
    for (int u = 0; u < g.V; ++u) {
        if (used[u]) node_id[u] = n++;
    }

    for (int u = 0; u < g.V; ++u) {
        if (!used[u]) continue;
        int kc = g.kmer_counts.empty() ? 0 : g.kmer_counts[u];
        int ms = node_scores.empty() ? 0 : node_scores[u];
        out << "S\t" << ele_id << "_" << node_id[u]
            << "\t" << g.nodes[u]
            << "\tKC:i:" << kc
            << "\tMS:i:" << ms
            << "\tEL:Z:" << ele_id << "\n";
    }
    for (int u = 0; u < g.V; ++u) {
        if (!used[u]) continue;
        int src = node_id[u];
        for (int ei = g.adj_off[u]; ei < g.adj_off[u + 1]; ++ei) {
            int v = g.adj_dst[ei];
            if (node_id[v] < 0) continue;
            out << "L\t" << ele_id << "_" << src << "\t+\t"
                << ele_id << "_" << node_id[v] << "\t+\t"
                << (k - 1) << "M\n";
        }
    }
    out << "P\tcontig_" << ele_id;
    for (int u : contig_path) {
        if (node_id[u] >= 0) out << "\t" << ele_id << "_" << node_id[u];
    }
    out << "\t" << std::string(contig_path.size(), '+') << "\n";

    for (size_t bi = 0; bi < bubble_paths.size(); ++bi) {
        out << "P\tcontig_" << ele_id << "_alt" << bi;
        for (int u : bubble_paths[bi]) {
            if (node_id[u] >= 0) out << "\t" << ele_id << "_" << node_id[u];
        }
        out << "\t" << std::string(bubble_paths[bi].size(), '+') << "\n";
    }
}

inline bool scan_snps_idx(
    const IndexedGraph& g,
    const std::vector<int>& contig_path,
    int ele_id, int k,
    std::ofstream& out,
    double min_af = 0.05,
    int max_branch_depth = 5)
{
    bool has_var = false;
    if (contig_path.empty()) return false;
    std::vector<int> path_pos(g.V, -1);
    for (size_t i = 0; i < contig_path.size(); ++i)
        path_pos[contig_path[i]] = static_cast<int>(i);

    std::vector<uint8_t> visited(g.V, 0);
    std::vector<int> frontier, next_frontier, touched;

    for (size_t i = 0; i + 1 < contig_path.size(); ++i) {
        int u = contig_path[i];
        int next_main = contig_path[i + 1];
        if (g.out_deg[u] < 2) continue;
        int main_ei = g.find_edge(u, next_main);

        for (int ei = g.adj_off[u]; ei < g.adj_off[u + 1]; ++ei) {
            int branch = g.adj_dst[ei];
            if (branch == next_main) continue;

            frontier.clear();
            touched.clear();
            frontier.push_back(branch);
            visited[branch] = 1;
            touched.push_back(branch);
            int reconnect_pos = -1;

            for (int depth = 1; depth <= max_branch_depth && !frontier.empty(); ++depth) {
                next_frontier.clear();
                for (int cur : frontier) {
                    for (int ce = g.adj_off[cur]; ce < g.adj_off[cur + 1]; ++ce) {
                        int cn = g.adj_dst[ce];
                        if (path_pos[cn] > static_cast<int>(i)) {
                            reconnect_pos = path_pos[cn];
                            break;
                        }
                        if (!visited[cn]) {
                            visited[cn] = 1;
                            touched.push_back(cn);
                            next_frontier.push_back(cn);
                        }
                    }
                    if (reconnect_pos >= 0) break;
                }
                frontier.swap(next_frontier);
                if (reconnect_pos >= 0) break;
            }
            for (int t : touched) visited[t] = 0;
            if (reconnect_pos < 0) continue;

            int branch_len = reconnect_pos - static_cast<int>(i);
            if (branch_len > 10) continue;

            int ref_cov = (main_ei >= 0 && !g.edge_counts.empty()) ? g.edge_counts[main_ei] : 0;
            int alt_cov = (!g.edge_counts.empty()) ? g.edge_counts[ei] : 0;
            if (ref_cov < 2 || alt_cov < 2) continue;

            double af = static_cast<double>(alt_cov) / (ref_cov + alt_cov);
            if (af < min_af) continue;

            int contig_pos = static_cast<int>(i + k - 1);
            const char* type = (branch_len == 1) ? "SNP" : "INDEL";
            char ref_base = g.nodes[next_main].back();
            char alt_base = g.nodes[branch].back();

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
