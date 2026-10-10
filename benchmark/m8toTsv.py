#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: JiangminZheng
Date: 2025-7-18
Description: From m8 out in a folder to tsv matrix
'''


import os, sys, re
import glob, tqdm
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict


def picker_prune(sp_name, raw_query, multi_Q_out):
    length_info = {}
    with open(raw_query) as query_f:
        for line in query_f.readlines():
            if line.startswith('>') and not line.startswith('>ref.'):
                CNE_id, s, e = re.search(r">(\S+)@\S+:(\d+)-(\d+)", line).groups()
                s, e = int(s), int(e)
                if e - s > length_info.get(CNE_id, 0):
                    length_info[CNE_id] = e - s
    local_len, local_loci = {}, {}

    with open(multi_Q_out, 'r') as blastn_f:
        for line in blastn_f.readlines():
            query, chr, s, e, strand, supp = line.split()
            s, e = int(s), int(e)
            query_len  = length_info[query.split('@')[0]]
            target_len = e - s + 1
            cov = target_len / query_len
            if cov > 1.5 or cov < 0.05:
                continue
            local_len[query.split('@')[0]] = min(target_len, query_len)
            local_loci[query.split('@')[0]] = f"{sp_name}.{chr}{strand}:{s}-{e}"

    return local_len, local_loci, length_info


def main():
    seed_fas = "/data01/zhengjiangmin/project/CNEfinding/resultsHub/06.nonLinearLoss/03.Rosales/data/Rosales.seed"
    total_cov = {}
    total_loci = {}
    lenght_info = {}
    for file in tqdm.tqdm(glob.glob('./test/*.out')):
        sp_name = os.path.basename(file).replace('.out', '')
        local_len, local_loci, length_info = picker_prune(sp_name, seed_fas, file)
        total_cov[sp_name] = local_len
        total_loci[sp_name] = local_loci

    total_cov_df = pd.DataFrame.from_dict(total_cov)
    total_cov_df = total_cov_df.fillna(0)
    total_cov_df[total_cov_df.mean(axis=0).sort_values().index].to_csv(
        'mmseqs@len.tsv', sep='\t'
    )
    print(total_cov_df)

    total_loci_df = pd.DataFrame.from_dict(total_loci)
    total_loci_df = total_loci_df.fillna('NA')
    total_loci_df.to_csv('mmseqs@loci.tsv', sep='\t')

    length_series = pd.Series(length_info).reindex(total_cov_df.index)
    df_normed = total_cov_df.div(length_series, axis=0)
    df_sampled = df_normed.sample(frac=0.05, random_state=42)

    fig = sns.clustermap(
        df_sampled, 
        figsize=(80, 20), 
    )

    plt.tight_layout()
    plt.savefig('test.png', format='png', dpi=300)

if __name__ == '__main__':
    main()
