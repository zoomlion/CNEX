#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>

namespace py = pybind11;

// 修改 De Bruijn 图构建函数
std::unordered_map<std::string, std::vector<std::string>> de_bruijn_graph(const std::vector<std::string>& reads, int k) {
    std::unordered_map<std::string, std::vector<std::string>> graph;
    
    for (const auto& read : reads) {
        if (read.length() < k) continue;
        
        // 处理每个k-mer
        for (size_t i = 0; i <= read.length() - k; ++i) {
            std::string kmer = read.substr(i, k);
            std::string prefix = kmer.substr(0, k-1);
            std::string suffix = kmer.substr(1, k-1);
            graph[prefix].push_back(suffix);
        }
    }
    
    return graph;
}

// 查找起始节点
std::string find_start_node(const std::unordered_map<std::string, std::vector<std::string>>& graph) {
    std::unordered_map<std::string, int> in_degree, out_degree;
    
    // 计算入度和出度
    for (const auto& node : graph) {
        out_degree[node.first] = node.second.size();
        for (const auto& neighbor : node.second) {
            in_degree[neighbor]++;
        }
    }
    
    // 首先尝试找到入度为0的节点
    for (const auto& node : graph) {
        if (in_degree[node.first] == 0) {
            return node.first;
        }
    }
    
    // 如果没有入度为0的节点，返回第一个节点
    return graph.begin()->first;
}

// 修改后的序列组装函数
std::string assemble_sequence(const std::unordered_map<std::string, std::vector<std::string>>& original_graph) {
    if (original_graph.empty()) {
        return "";
    }
    
    // 创建图的可修改副本
    auto graph = original_graph;
    
    // 找到起始节点
    std::string start = find_start_node(graph);
    
    // 存储最终路径
    std::string result = start;
    std::string current = start;
    
    // 贪婪策略：总是选择可用的下一个节点
    while (true) {
        auto it = graph.find(current);
        if (it == graph.end() || it->second.empty()) {
            break;
        }
        
        // 获取下一个节点
        std::string next = it->second.back();
        it->second.pop_back();
        
        // 只添加新节点的最后一个字符
        result += next.back();
        current = next;
    }
    
    return result;
}

// output_gml 函数保持不变
void output_gml(const std::unordered_map<std::string, std::vector<std::string>>& graph, const std::string& filename) {
    std::ofstream file(filename);
    file << "graph [\n";
    
    int node_id = 0;
    std::unordered_map<std::string, int> node_map;
    for (const auto& node : graph) {
        if (node_map.find(node.first) == node_map.end()) {
            node_map[node.first] = node_id++;
            file << "  node [\n    id " << node_map[node.first] << "\n    label \"" << node.first << "\"\n  ]\n";
        }
    }
    
    for (const auto& node : graph) {
        for (const auto& neighbor : node.second) {
            file << "  edge [\n    source " << node_map[node.first] << "\n    target " << node_map[neighbor] << "\n  ]\n";
        }
    }
    
    file << "]\n";
}

PYBIND11_MODULE(debruijn, m) {
    m.def("de_bruijn_graph", &de_bruijn_graph, "Build a De Bruijn graph from reads", py::arg("reads"), py::arg("k"));
    m.def("assemble_sequence", &assemble_sequence, "Assemble sequence from De Bruijn graph", py::arg("graph"));
    m.def("output_gml", &output_gml, "Output De Bruijn graph as GML file", py::arg("graph"), py::arg("filename"));
}