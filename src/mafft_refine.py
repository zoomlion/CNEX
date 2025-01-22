#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: JiangminZheng
Date: 2025-1-9
Description: Refine MSA by realign with mafft
'''
import os, re, sys
import multiprocessing
import subprocess
import itertools
import argparse
from tempfile import NamedTemporaryFile
import shutil


def bunch2fas(bunch: str) -> dict:
    """
    Convert a bunch string to a FASTA dictionary, with sequence IDs as keys and sequences as values.
    Sequences with only '-' are ignored.
    """
    fas_dict = {}
    for line in bunch.split('\n'):
        if not line:
            continue
        if line.startswith('>'):
            seq_id = line.strip('>').strip()
            fas_dict[seq_id] = ''
        else:
            fas_dict[seq_id] += line.strip().upper()

    # Remove sequences with no valid characters (e.g., only '-')
    return {
        seq_id: seq
        for seq_id, seq in fas_dict.items()
        if len(set(seq) - {'-'}) > 0 and not seq_id.startswith('ref')  # Exclude sequences with only '-' and 'ref' prefix
    }


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


def align_bunch(args):
    """
    Align a bunch of sequences using MAFFT.
    
    Args:
        fas_dict (dict): A dictionary with sequence IDs as keys and sequences as values.
        mafft (str): Path to the MAFFT executable.

    Returns:
        str: Aligned sequences in FASTA format.
    """
    fas_dict, mafft = args
    input_string = '\n'.join([f">{seq_id}\n{seq}" for seq_id, seq in fas_dict.items()])

    # save to temporary file at '/dev/shm' to avoid OOM
    with NamedTemporaryFile(mode='w', dir='/dev/shm', delete=True) as temp_file:
        temp_file.write(input_string)
        temp_file.flush()
        mafft_cmd = [mafft, '--preservecase', '--quiet', '--retree', '1', '--thread', '1', temp_file.name]
        try:
            result = subprocess.run(
                mafft_cmd,
                text=True,           # Handle input/output as text
                capture_output=True, # Capture stdout and stderr
                check=True           # Raise an error if the command fails
            )
            # Return the aligned sequences from stdout and set chop false
            # to avoid n.a. chop in fasta output
            out_dict = bunch2fas(result.stdout)
            return '\n'.join([f">{seq_id}\n{seq}" for seq_id, seq in out_dict.items()]) + '\n'
        except subprocess.CalledProcessError as e:
            # Handle errors (e.g., if MAFFT fails)
            raise RuntimeError(f"MAFFT alignment failed: {e.stderr}") from e


def process_in_chunks(generator, output_file, mafft, chunk_size, threads=4):
    """
    Process a generator in chunks, and write the output to a file to save ram usage.
    """
    open(output_file, 'w').close()
    with multiprocessing.Pool(processes=threads) as pool:
        while True:
            chunk = list(itertools.islice(generator, chunk_size))
            if not chunk:
                break
            # Align the chunk using mafft with mp.Pool
            results = pool.map(align_bunch, [(bunch2fas(item), mafft) for item in chunk])
            # Write the results to the output file
            with open(output_file, 'a') as file:
                for result in results:
                    file.write('###\n')
                    file.write(result)


def main():
    parser = argparse.ArgumentParser(description='Refine MSA by realign with mafft')
    parser.add_argument('input', help='input MSA file')
    parser.add_argument('output', help='output MSA file')
    parser.add_argument('--mafft', type=str, required=True, help='mafft executable')
    parser.add_argument('--min_c', type=int, default=3, help='minimum number of records in a bunch')
    parser.add_argument('-t', '--threads', type=int, default=4, help='number of threads to use')
    args = parser.parse_args()

    # check diff in input and output
    if os.path.exists(args.output):
        if os.path.samefile(args.input, args.output):
            raise ValueError('input and output are the same file')

    # check if mafft exists and executable
    if not os.path.isfile(args.mafft) or not os.access(args.mafft, os.X_OK):
        raise ValueError('mafft not found or not executable')
    
    process_in_chunks(bunch_generator(args.input), args.output, args.mafft, 1000, args.threads)


if __name__ == '__main__':
    main()
