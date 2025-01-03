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
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <queue>
#include <set>

namespace py = pybind11;
using Graph = std::unordered_map<std::string, std::vector<std::string>>;

Graph de_bruijn_graph(const std::vector<std::string>& reads, u_int64_t k) {
    Graph graph;
    std::set<std::pair<std::string, std::string>> valid_edges;
    
    // 首先收集所有有效的边
    for (const auto& read : reads) {
        if (read.length() < (k+1)) continue;
        
        for (size_t i = 0; i <= read.length() - (k+1); ++i) {
            std::string kmer = read.substr(i, (k+1));
            std::string prefix = kmer.substr(0, k);
            std::string suffix = kmer.substr(1, k);
            valid_edges.insert({prefix, suffix});
        }
    }
    
    // 只添加有效的边到图中
    for (const auto& edge : valid_edges) {
        graph[edge.first].push_back(edge.second);
    }
    
    return graph;
}


std::unordered_map<std::string, int> dijkstra(const Graph& graph, const std::string& start) {
    std::unordered_map<std::string, int> distances;

    // Initialize distances for all nodes in the graph
    for (const auto& [node, neighbors] : graph) {
        distances[node] = std::numeric_limits<int>::max();
        for (const auto& neighbor : neighbors) {
            if (distances.find(neighbor) == distances.end()) {
                distances[neighbor] = std::numeric_limits<int>::max();
            }
        }
    }
    distances[start] = 0;

    // Priority queue for selecting the next node with the smallest tentative distance
    std::priority_queue<std::pair<int, std::string>, std::vector<std::pair<int, std::string>>, std::greater<>> pq;
    pq.push({0, start});

    while (!pq.empty()) {
        auto [current_distance, current_node] = pq.top();
        pq.pop();

        if (current_distance > distances[current_node]) continue;

        // Process all neighbors of the current node
        if (graph.find(current_node) != graph.end()) {
            for (const auto& neighbor : graph.at(current_node)) {
                int new_distance = current_distance + 1; // Assuming unweighted edges

                if (new_distance < distances[neighbor]) {
                    distances[neighbor] = new_distance;
                    pq.push({new_distance, neighbor});
                }
            }
        }
    }

    return distances;
}

std::pair<std::string, std::string> find_farthest_start_end_nodes(const Graph& graph) {
    std::unordered_map<std::string, int> in_degree, out_degree;
    // Calculate in-degree and out-degree for each node
    for (const auto& node : graph) {
        out_degree[node.first] = node.second.size();
        for (const auto& neighbor : node.second) {
            in_degree[neighbor]++;
            if (out_degree.find(neighbor) == out_degree.end()) {
                out_degree[neighbor] = 0;
            }
        }
    }

    // Identify candidate start and end nodes based on degree differences
    std::set<std::string> start_candidates, end_candidates;
    for (const auto& node : graph) {
        if (out_degree[node.first] - in_degree[node.first] > 0) {
            start_candidates.insert(node.first);
        }
        if (in_degree[node.first] - out_degree[node.first] > 0) {
            end_candidates.insert(node.first);
        }
        // check neighbours by the way
        for (const auto& neighbor : node.second) {
            if (out_degree[neighbor] - in_degree[neighbor] > 0) {
                start_candidates.insert(neighbor);
            }
            if (in_degree[neighbor] - out_degree[neighbor] > 0) {
                end_candidates.insert(neighbor);
            }
        }
    }

    // If no candidates are found, assume any node can be a start or end
    if (start_candidates.empty() || end_candidates.empty()) {
        for (const auto& node : graph) {
            start_candidates.insert(node.first);
            end_candidates.insert(node.first);
            std::cout << "Warning: No start or end candidates found, assuming any node can be a start or end.\n";
            break;
        }
    }

    std::pair<std::string, std::string> farthest_pair;
    int max_distance = -1;

    // Compute the distance between all pairs of start and end candidates
    for (const auto& start : start_candidates) {
        auto distances = dijkstra(graph, start);

        for (const auto& end : end_candidates) {
            if (distances[end] != std::numeric_limits<int>::max() && distances[end] > max_distance) {
                max_distance = distances[end];
                farthest_pair = {start, end};
            }
            // std::cout << "Distance from " << start << " to " << end << ": " << distances[end] << '\n';
        }
    }

    // // print selected start and end nodes
    // std::cout << "Selected start node: " << farthest_pair.first << '\n';
    // std::cout << "Selected end node: " << farthest_pair.second << '\n';
    return farthest_pair;
}

std::vector<std::string> find_longest_path(
    const Graph& graph,
    const std::string& start,
    const std::string& end,
    u_int64_t max_path_length = 2500  // 添加最大路径长度限制
) {
    Graph paths;
    std::unordered_map<std::string, int> visit_count;
    std::vector<std::string> longest_path;
    
    std::function<void(const std::string&, std::vector<std::string>&)> dfs = 
    [&](const std::string& current, std::vector<std::string>& current_path) {
        // 添加访问次数限制
        if (visit_count[current] >= 2) return;
        if (current_path.size() > max_path_length) return;
        
        visit_count[current]++;
        current_path.push_back(current);
        
        if (current_path.size() > longest_path.size()) {
            longest_path = current_path;
        }
        
        if (graph.count(current)) {
            for (const auto& next : graph.at(current)) {
                dfs(next, current_path);
            }
        }
        
        current_path.pop_back();
        visit_count[current]--;
    };
    
    std::vector<std::string> current_path;
    dfs(start, current_path);
    
    return longest_path;
}

std::string assemble_sequence(const Graph& graph) {
    if (graph.empty()) {
        return "";
    }
    
    auto [start_node, end_node] = find_farthest_start_end_nodes(graph);
    auto path = find_longest_path(graph, start_node, end_node);
    
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
    std::unordered_map<std::string, int> node_map;

    // allocate ID
    for (const auto& node : graph) {
        if (node_map.find(node.first) == node_map.end()) {
            node_map[node.first] = node_id++;
            file << "  node [\n    id " << node_map[node.first] << "\n    label \"" << node.first << "\"\n  ]\n";
        }
        for (const auto& neighbor : node.second) {
            if (node_map.find(neighbor) == node_map.end()) {
                node_map[neighbor] = node_id++;
                file << "  node [\n    id " << node_map[neighbor] << "\n    label \"" << neighbor << "\"\n  ]\n";
            }
        }
    }

    for (const auto& node : graph) {
        for (const auto& neighbor : node.second) {
            file << "  edge [\n    source " << node_map[node.first] << "\n    target " << node_map[neighbor] << "\n  ]\n";
            // std::cout << "Edge: " << node.first << " -> " << neighbor << '\n';
            // std::cout << "Mapped IDs: " << node_map[node.first] << " -> " << node_map[neighbor] << '\n';
        }
    }

    file << "]\n";
}

PYBIND11_MODULE(debruijn, m) {
    m.doc() = "De Bruijn graph and sequence assembler";

    m.def("de_bruijn_graph", &de_bruijn_graph, "Build a De Bruijn graph from reads", py::arg("reads"), py::arg("k"));
    m.def("assemble_sequence", &assemble_sequence, "Assemble sequence from De Bruijn graph", py::arg("graph"));
    m.def("output_gml", &output_gml, "Output De Bruijn graph as GML file", py::arg("graph"), py::arg("filename"));
}