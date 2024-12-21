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

std::unordered_map<std::string, std::vector<std::string>> de_bruijn_graph(const std::vector<std::string>& reads, int k) {
    std::unordered_map<std::string, std::vector<std::string>> graph;
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

std::pair<std::string, std::string> find_start_end_nodes(const std::unordered_map<std::string, std::vector<std::string>>& graph) {
    std::unordered_map<std::string, int> in_degree, out_degree;
    
    for (const auto& node : graph) {
        out_degree[node.first] = node.second.size();
        for (const auto& neighbor : node.second) {
            in_degree[neighbor]++;
            out_degree[neighbor]; // 确保节点存在于out_degree中
        }
    }
    
    std::string start_node = "", end_node = "";
    int max_degree_diff = 0;
    
    // 寻找入度出度差最大的节点作为起点和终点
    for (const auto& node : graph) {
        int in = in_degree[node.first];
        int out = out_degree[node.first];
        int diff = out - in;
        
        if (diff > max_degree_diff) {
            max_degree_diff = diff;
            start_node = node.first;
        }
        if (-diff > max_degree_diff) {
            max_degree_diff = -diff;
            end_node = node.first;
        }
    }
    
    // 如果没有找到明显的起点或终点，使用度数最大的节点
    if (start_node.empty() || end_node.empty()) {
        int max_degree = 0;
        std::string max_degree_node;
        
        for (const auto& node : graph) {
            int total_degree = in_degree[node.first] + out_degree[node.first];
            if (total_degree > max_degree) {
                max_degree = total_degree;
                max_degree_node = node.first;
            }
        }
        
        if (start_node.empty()) start_node = max_degree_node;
        if (end_node.empty()) end_node = max_degree_node;
    }
    
    return {start_node, end_node};
}

std::vector<std::string> find_longest_path(
    const std::unordered_map<std::string, std::vector<std::string>>& graph,
    const std::string& start,
    const std::string& end,
    int max_path_length = 1500  // 添加最大路径长度限制
) {
    std::unordered_map<std::string, std::vector<std::string>> paths;
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

std::string assemble_sequence(const std::unordered_map<std::string, std::vector<std::string>>& graph) {
    if (graph.empty()) {
        return "";
    }
    
    auto [start_node, end_node] = find_start_end_nodes(graph);
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

void output_gml(const std::unordered_map<std::string, std::vector<std::string>>& graph, const std::string& filename) {
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