#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: JiangminZheng
Date: 2024-12-16
"""

import os, re, sys
import glob
import argparse
import tqdm
import subprocess
import shutil
import string
import random
from collections import defaultdict
from hip import debruijn
from tempfile import NamedTemporaryFile


def reverse_complement(seq):
    """
    Return the reverse complement of a DNA sequence.
    with maketrans method
    """
    reverse_map = str.maketrans("ATCG", "TAGC")
    return seq.translate(reverse_map)[::-1]


def debruijn_assembler(reads, k):
    """
    Assemble a longest DNA sequence from a set of reads using a De Bruijn Graph.
    """
    graph = debruijn.de_bruijn_graph(reads, k)
    return debruijn.assemble_sequence(graph)


def reads_generator(input_f: str):
    """
    get a generator for reads from input file
    """
    prev_ele_id = None
    temp_reads = {}
    with open(input_f) as f:
        for line in f.readlines():
            reads_id, strand, ele_id, seq = line.strip().split("\t")
            seq = seq.upper()
            strand = int(strand)
            if prev_ele_id is None:
                prev_ele_id = ele_id
            if ele_id != prev_ele_id:
                yield prev_ele_id, temp_reads
                temp_reads = {}
                temp_reads[reads_id] = (strand, seq)
                prev_ele_id = ele_id
            else:
                temp_reads[reads_id] = (strand, seq)
        yield prev_ele_id, temp_reads


def bunch_generator(fasta_file):
    """
    Build a generator to yield bunched FASTA sequences.
    """
    buffer = ''
    with open(fasta_file, 'r') as file:
        for line in file:
            buffer += line
            while '###\n' in buffer:
                part, buffer = buffer.split('###\n', 1)
                if part.strip():
                    yield part.strip()
        # Yield any remaining part after the loop
        if buffer.strip():
            yield buffer.strip()


def bunch2conservedfas(bunch: str) -> dict:
    """
    Convert a bunch string to longest fasta, with sequence IDs as keys and sequences as values.
    Sequences with only '-' are ignored.
    """
    def get_most_frequent_base(base_list):
        return max(set(base_list), key=base_list.count)

    fas_dict = {}
    for line in bunch.split('\n'):
        if not line:
            continue
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            fas_dict[seq_id] = ''
        else:
            if set(line.strip()) - set('-'):  # not all '-'
                fas_dict[seq_id] += line.strip()
            else:
                fas_dict.pop(seq_id, None)

    seq_len = len(list(fas_dict.values())[0])
    #  fas_dict.values()
    most_cons = ''.join(
        [get_most_frequent_base([seq[i] for seq in fas_dict.values()])
          for i in range(seq_len)]
    )

    # select the longest sequence as ref
    conserv = 0
    ref_seq_id = None
    for seq_id, seq in fas_dict.items():
        local_conserv = 0
        if seq_id.startswith('ref'):
            continue
        for i in range(seq_len):
            if seq[i] == most_cons[i]:
                local_conserv += 1
        if local_conserv > conserv:
            conserv = local_conserv
            ref_seq_id = seq_id

    if ref_seq_id is None:
        raise ValueError("No reference sequence found.")

    return f">{ref_seq_id}\n{fas_dict[ref_seq_id].replace('-', '')}\n"


def merge_loci(locus):
    """
    merge loci by their range: return a dict of affiliation
    e.g. 16:93246851-93247000
         16:93246901-93247050
         16:93246951-93247100
         16:93247001-93247150
         16:93247051-93247200 --> {'16:93246851-93247200': [16:93246851-93247000, .., 16:93247051-93247200]}
    """
    intervals = []
    loci_dict = defaultdict(list)
    for loci in locus:
        chrom, s, e = re.search(r'^(\S+):(\d+)-(\d+)$', loci).groups()
        intervals.append((chrom, int(s), int(e)))
    intervals.sort(key=lambda x: (x[0], x[1]))
    if len(intervals) == 1:
        info = f"{intervals[0][0]}:{intervals[0][1]}-{intervals[0][2]}"
        loci_dict[info] = [info]
    else:
        curr_chrom, curr_start, curr_end = intervals[0]
        current_group = [f"{curr_chrom}:{curr_start}-{curr_end}"]
        
        for i in range(1, len(intervals)):
            chrom, start, end = intervals[i]
            # If different chromosome or no overlap
            if chrom != curr_chrom or start > curr_end:
                # Save current group
                merged = f"{curr_chrom}:{curr_start}-{curr_end}"
                loci_dict[merged] = current_group
                # Start new group
                curr_chrom, curr_start, curr_end = chrom, start, end
                current_group = [f"{chrom}:{start}-{end}"]
            else:
                # Extend current range
                curr_end = max(curr_end, end)
                current_group.append(f"{chrom}:{start}-{end}")
        
        # Add the last group
        merged = f"{curr_chrom}:{curr_start}-{curr_end}"
        loci_dict[merged] = current_group
    
    return loci_dict


def axt2rank(
        axt_str: str, 
        index2id: dict, 
        int_assembled_seqs: dict, 
        mode: str, 
        min_percent=0.3
) -> list:
    """
    parse axt alignment string and return a list of (id, seq) pairs 
    in order of alignment scores, tagged with score rank
    """
    ranked_seqs = []
    index2score = {}
    aln_info = defaultdict(list)
    max_score = 0
    for line in axt_str.split('\n'):
        if line.startswith('#') or not line:
            continue
        if not line[0].isdigit():
            continue
        fields = line.strip().split()
        # all 1-based coordinates
        _, ref_id, ref_s, ref_e, query_index, query_s, query_e, _, score = fields
        if float(score) > max_score:
            max_score = float(score)
        ref_s, ref_e, query_s, query_e, score = map(int, [ref_s, ref_e, query_s, query_e, score])
        aln_info[int(query_index)].append((query_s, query_e, score))

    if len(aln_info) == 0:
        return []

    for query_index in aln_info:
        for s, e, score in sorted(aln_info[query_index], key=lambda x: (x[2], -x[0]), reverse=True):
            if score < min_percent * max_score:
                continue
            if query_index not in index2score:
                index2score[query_index] = (s, e, score)
            else:
                prev_s, prev_e, prev_score = index2score[query_index]
                if s > prev_e:
                    index2score[query_index] = (prev_s, e, score)
                elif e < prev_s:
                    index2score[query_index] = (s, prev_e, score)
    
    if mode=='genome':
        # wabble query_id and query_seq by query_s and query_e
        for rank, (index, (s_bias, e_bias, score)) in enumerate(
            sorted(index2score.items(), key=lambda x: x[1][2], reverse=True)
        ):
            s_bias, e_bias, score = map(int, [s_bias, e_bias, score])
            query_id = index2id[index]
            query_seq = int_assembled_seqs[index]
            chrom_info, s, e = re.search(
                r'^(\S+[+-]):(\d+)-(\d+)$', query_id
            ).groups()
            s, e = map(int, [s, e])
            query_id = f"{chrom_info}:{s+s_bias-1}-{s+e_bias-1}#{rank}"
            query_seq = query_seq[s_bias-1:e_bias]
            ranked_seqs.append((query_id, query_seq))
    elif mode=='fq':
        for rank, (index, (s_bias, e_bias, score)) in enumerate(
            sorted(index2score.items(), key=lambda x: x[1][2], reverse=True)
        ):
            s_bias, e_bias, score = map(int, [s_bias, e_bias, score])
            query_id = index2id[index]
            query_seq = int_assembled_seqs[index]
            query_id = f"{query_id}#{rank}"
            query_seq = query_seq[s_bias-1:e_bias]
            ranked_seqs.append((query_id, query_seq))

    return ranked_seqs


def assemble_reads(ele_id: str, raw_seq: str, temp_reads: dict, mode: str, k=35):
    """
    assemble reads using de bruijn graph from hip
    """
    if mode == 'fq':
        reads = []
        for reads_id, (strand, seq) in temp_reads.items():
            local_seq = seq if strand == 1 else reverse_complement(seq)
            reads.append(local_seq)
        assembled = debruijn_assembler(reads, k)
        return {ele_id: assembled}
    elif mode == 'genome':
        locus = []
        reads = {}
        for reads_id, (strand, seq) in temp_reads.items():
            local_seq = seq if strand == 1 else reverse_complement(seq)
            strand = '+' if strand == 1 else '-'
            chrom, s, e = re.search(r'^(\S+):(\d+)-(\d+)$', reads_id).groups()
            locus.append(f"{chrom}{strand}:{s}-{e}")
            reads[f"{chrom}{strand}:{s}-{e}"] = local_seq
        if len(locus) > 100:  # too many loci, skip
            return {}
        loci_dict = merge_loci(locus)
        all_assembled = []
        for merged_locus, reads_list in loci_dict.items():
            assembled = debruijn_assembler(
                [reads[loci] for loci in reads_list], k
            )
            assembled_id = f"{ele_id}@{merged_locus}"
            all_assembled.append((assembled_id, assembled))
        return {id: seq for id, seq in all_assembled}


def main():
    parser = argparse.ArgumentParser(description="De Bruijn Graph Assembly")
    parser.add_argument("inputs_dir", type=str, help="input dirs holding files")
    parser.add_argument("-r", "--raw_msa_file", type=str, help="raw msa file", required=True)
    parser.add_argument("-k", "--kmer", type=int, default=35, help="k-mer length")
    parser.add_argument("-l", "--lastz_exec", type=str, default="lastz", help="path to lastz")
    parser.add_argument("-o", "--output", type=str, default="assembled.fasta", help="output file")
    parser.add_argument("-t", "--threads", type=int, default=8, help="number of threads")
    args = parser.parse_args()

    if sys.version_info < (3, 7):
        raise ValueError("Python version must be >= 3.7 for ordered dict support.")

    # Check if input directory exists
    if not os.path.exists(args.inputs_dir):
        raise FileNotFoundError(f"Input directory not found: {args.inputs_dir}")

    # Check if Lastz executable exists and is executable
    if not os.path.exists(args.lastz_exec):
        raise FileNotFoundError(f"Lastz executable not found: {args.lastz_exec}")
    if not os.access(args.lastz_exec, os.X_OK):
        raise PermissionError(f"Lastz executable not executable: {args.lastz_exec}")

    # cat and sort all reads to tempfile
    input_f = f"{''.join([random.choice(string.ascii_letters) for _ in range(8)])}.reads"
    os.system(f"cat {args.inputs_dir}/*.reads | sort -k3n > {input_f}")
    if not os.path.exists(input_f):
        raise FileNotFoundError(f"File {input_f} not found.")
    
    # cache_dir in /dev/shm with 8 chars
    cache_dir = f"/dev/shm/{''.join([random.choice(string.ascii_letters) for _ in range(8)])}"
    os.makedirs(cache_dir, exist_ok=True)
    
    refs = {}
    for bunch_id, bunch in enumerate(bunch_generator(args.raw_msa_file)):
        ref_seq = bunch2conservedfas(bunch)
        refs[f"{bunch_id}"] = ref_seq

    # get reads type by first line
    # if starts with '@' is read mode, elif : in first line is fasta mode
    type = None
    with open(input_f) as f:
        first_line = f.readline().strip()
        if first_line.startswith('@'):
            type = 'fq'
        elif ':' in first_line:
            type = 'genome'
        else:
            raise ValueError(f"Unknown file type: {first_line}")

    results = []
    for ele_id, temp_reads in reads_generator(input_f):
        if not ele_id:
            continue
        print(f"\rAssembling {int(ele_id)+1}", end="", flush=True)
        # step 1: assemble all possible reads
        assembled_seqs = assemble_reads(ele_id, args.raw_msa_file, temp_reads, mode=type, k=args.kmer)
        if not assembled_seqs:
            continue
        int_assembled_seqs = {index: seq for index, seq in enumerate(assembled_seqs.values())}
        index2id = {index: id for index, id in enumerate(assembled_seqs.keys())}

        # step 2: align to reference by lastz
        #  a. filter noise sequences
        #  b. get paralogs
        filtered_seqs = []
        # pack seqeucnes into ref and query files
        with NamedTemporaryFile(mode='w', delete=True, dir=cache_dir) as ref_file, \
             NamedTemporaryFile(mode='w', delete=True, dir=cache_dir) as query_file:
            ref_file.write(refs[ele_id])
            ref_file.flush()
            query_file.write(''.join([f">{index}\n{seq}\n" for index, seq in int_assembled_seqs.items()]))
            query_file.flush()
            lastz_cmd = [
                args.lastz_exec, 
                ref_file.name, 
                query_file.name, 
                '--hspthresh=500', 
                '--gappedthresh=2000', 
                '--strand=forward', 
                '--ambiguous=iupac', 
                '--chain', 
                '--axt'
            ]
            align_out = subprocess.run(lastz_cmd, check=True, stdout=subprocess.PIPE).stdout.decode()
            filtered_seqs = axt2rank(align_out, index2id, int_assembled_seqs, mode=type)

        results.extend(filtered_seqs)
    print()
    
    # write to fasta file
    with open(args.output, "w") as f:
        for (id, seq) in results:
            f.write(f">{id}\n{seq}\n")

    # clean up input_f
    os.remove(input_f)
    shutil.rmtree(cache_dir)
            
if __name__ == "__main__":
    main()
