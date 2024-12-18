#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: JiangminZheng
Date: 2024-12-16
'''
import debruijn

# 示例读取数据
reads = []
with open('./Assemble.reads') as f:
    for line in f:
        line = line.strip()
        eles = line.split('\t')
        if len(eles) > 2:
            continue
        id, read = eles
        if int(id) == 17843:
            reads.append(read)

reads = reads[::-1]

print(len(reads))
k = 21

graph = debruijn.de_bruijn_graph(reads, k)

assembled_sequence = debruijn.assemble_sequence(graph)
print(f"Assembled sequence: {assembled_sequence}")

debruijn.output_gml(graph, "debruijn_graph.gml")
print("GML file has been saved as 'debruijn_graph.gml'.")
