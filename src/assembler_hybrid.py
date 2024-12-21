#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: JiangminZheng
Date: Use reads to assemble cne
"""

import os, re, sys
import mappy
from collections import defaultdict
import gzip
import argparse
import tqdm
import multiprocessing
from hip import debruijn
from hip.validator import validate_read, MerQueryManager


def fq_generator(file: str, divider: int = 1, mark: int = 0):
    """
    Generate reads from a FASTQ file (supports both gzipped and plain text files).
    Only yields entries where the line number is divisible by the divider parameter.

    Args:
        file (str): Path to the FASTQ file (can be gzipped with .gz or plain text).
        divider (int): Only yield entries where line number is divisible by this value.
                      Default is 1 (yield all entries).

    Yields:
        tuple: A tuple containing:
            - seq_id (str): The sequence ID (from the first line, starting with @).
            - sequence (str): The nucleotide sequence (from the second line).
            - quality (str): The quality scores (from the fourth line).
    """
    # Open file (supports both gzipped and plain text)
    open_func = gzip.open if file.endswith(".gz") else open

    with open_func(file, "rt") as f:
        entry_count = 0  # Counter for FASTQ entries

        while True:
            try:
                # Read the 4 lines of each FASTQ entry
                seq_id = f.readline().strip()  # 1st line: sequence identifier
                if not seq_id:  # EOF
                    break

                sequence = f.readline().strip()  # 2nd line: nucleotide sequence
                _ = f.readline().strip()  # 3rd line: separator (usually "+")
                quality = f.readline().strip()  # 4th line: quality scores

                entry_count += 1

                # Only yield entries where entry_count is divisible by divider
                if entry_count % divider == mark:
                    yield (seq_id, sequence, quality)

            except Exception as e:
                print(f"Error while reading file: {e}")
                break


def out_write(file, infos):
    with open(file, "a") as f:
        f.write("\n".join(["\t".join([str(ele) for ele in info]) for info in infos]))


def assembler(reads_files, mers, assemble_out, depth=20_000_000, threads=4, thread_id=0):
    # read in mers
    mer_query = MerQueryManager()
    mer_size = 0
    with open(mers) as f:
        for line in tqdm.tqdm(f, unit="mer"):
            mer, id, loci, count = line.strip("\n").split("\t")
            mer_size = len(mer)
            try:
                mer_query.add_mer(mer, int(id), int(loci))
            except ValueError:
                continue
    if mer_size == 0:
        raise ValueError("Confident mer size suggested to be over 7")
    # read in reads
    confi_reads = []
    cache_size = 5_000
    open(assemble_out, "w").close()
    for file in reads_files:
        fqs = fq_generator(file, divider=threads, mark=thread_id)
        for seq_id, seq, qua in tqdm.tqdm(fqs, total=depth, unit="reads"):
            depth -= 1
            # forward
            confi_id = validate_read(seq, mer_query, mer_size)
            if confi_id > -1:
                confi_reads.append((confi_id, seq))
            if len(confi_reads) > cache_size:
                # write
                out_write(assemble_out, confi_reads)
                confi_reads = []
            if depth < 0:
                break

    if confi_reads:
        out_write(assemble_out, confi_reads)


def main():
    parser = argparse.ArgumentParser(
        description="Assemble element based on confi mers."
    )
    parser.add_argument("reads_files", nargs="+", help="input reads file in gz or not")
    parser.add_argument("--mers", type=str, required=True, help="mers table")
    parser.add_argument(
        "--depth", type=int, default=200_000_000, help="reads number to assemble"
    )
    parser.add_argument("--threads", type=int, default=4, help="threads number to use")
    parser.add_argument("--output_dir", type=str, default="out", help="output directory")
    args = parser.parse_args()

    average_depth = args.depth // (len(args.reads_files) * args.threads)

    pool = multiprocessing.Pool(processes=args.threads)
    os.makedirs(args.output_dir, exist_ok=True)
    for i in range(args.threads):
        pool.apply_async(
            assembler,
            args=(
                args.reads_files,
                args.mers,
                os.path.join(args.output_dir, f"Assemble.{i}.reads"),
                average_depth,
                args.threads,
                i,
            ),
        )
    pool.close()
    pool.join()

    # assembler(
    #     args.reads_files, 
    #     args.mers, 
    #     args.assemble_out, 
    #     args.depth, 
    #     args.threads
    #     )


if __name__ == "__main__":
    main()
