#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: JiangminZheng
Description: Assemble CNE using sequencing reads with progress tracking
"""

import os
import re
import sys
import mappy
import gzip
import argparse
import tqdm
import multiprocessing as mp
import shutil
from collections import defaultdict
from queue import Empty
from hip import debruijn
from hip.validator import validate_read, MerQueryManager


def chunk_generator(file: str, chunk_size: int = 30000):
    """
    Generate chunks of reads from a FASTQ file or a GENOME fasta file.
    
    Args:
        file (str): Path to file (gzipped or plain text)
        chunk_size (int): Number of reads per chunk
        
    Yields:
        list: List of tuples containing (seq_id, sequence, quality) for each chunk
    """
    open_func = gzip.open if file.endswith(".gz") else open
    chunk = []
    
    with open_func(file, "rt") as f:
        first_line = f.readline()
        if first_line.startswith('>'):
            type = 'genome'
            f.seek(0)  # Rewind to include the first line in processing
        elif first_line.startswith('@'):
            type = 'fq'
            f.seek(0)
        else:
            raise ValueError("Unsupported file type")
        window_size = 150
        step_size = 50 
        if type == 'fq':
            lines = f
            while True:
                seq_id = next(lines, None)
                if seq_id is None:
                    if chunk:
                        yield chunk
                    break
                seq_id = seq_id.strip()
                
                sequence = next(lines, None)
                if sequence is None:
                    break
                sequence = sequence.strip()
                
                plus_line = next(lines, None)
                if plus_line is None:
                    break
                plus_line = plus_line.strip()
                
                quality = next(lines, None)
                if quality is None:
                    break
                quality = quality.strip()
                
                chunk.append((seq_id, sequence, '')) # strip quality
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
        elif type == 'genome':
            seq_id = ''
            seq = ''
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if seq_id != '' and seq != '':
                        # Generate pseudo-reads for the previous sequence
                        seq_len = len(seq)
                        for start in range(0, seq_len, step_size):
                            end = start + window_size
                            if end > seq_len:
                                end = seq_len
                            if end - start < window_size:
                                break  # Skip incomplete reads
                            read_id = f"{seq_id}:{start+1}-{end}"
                            # if all lowercase, skip
                            if not re.search('[A-Z]', seq[start:end]):
                                continue
                            chunk.append((read_id, seq[start:end].upper(), ''))
                            if len(chunk) >= chunk_size:
                                yield chunk
                                chunk = []
                        seq = ''
                    # Start new sequence
                    seq_id = line[1:].split()[0]  # Use the first word after '>'
                else:
                    seq += line
            # Generate pseudo-reads for the last sequence
            if seq_id != '' and seq != '':
                
                seq_len = len(seq)
                for start in range(0, seq_len, step_size):
                    end = start + window_size
                    if end > seq_len:
                        end = seq_len
                    if end - start < window_size:
                        break  # Skip incomplete reads
                    read_id = f"{seq_id}:{start+1}-{end}"
                    # if all lowercase, skip
                    if not re.search('[A-Z]', seq[start:end]):
                        continue
                    chunk.append((read_id, seq[start:end].upper(), ''))
                    if len(chunk) >= chunk_size:
                        yield chunk
                        chunk = []
            if chunk:
                yield chunk


def reader_process(files, queue, chunk_size, total_depth, progress_queue):
    """
    Reader process that reads files and puts data chunks into queue.
    
    Args:
        files (list): List of input FASTQ files
        queue (Queue): Queue for data chunks
        chunk_size (int): Size of each chunk
        total_depth (int): Maximum number of reads to process
        progress_queue (Queue): Queue for progress updates
    """
    reads_count = 0
    for file in files:
        for chunk in chunk_generator(file, chunk_size):
            if reads_count >= total_depth:
                break
            queue.put(chunk)
            reads_count += len(chunk)
            progress_queue.put(len(chunk))
    
    # Send end signals
    for _ in range(mp.cpu_count()):
        queue.put(None)
    progress_queue.put(None)


def worker_process(queue, mer_file, output_file, thread_id):
    """
    Worker process that processes data chunks from queue.
    
    Args:
        queue (Queue): Queue containing data chunks
        mer_file (str): Path to mer file
        output_file (str): Path to output file
        thread_id (int): Worker thread ID
    """
    # Initialize mer query manager
    mer_query = MerQueryManager()
    with open(mer_file) as f:
        for line in f:
            mer, id, loci, count = line.strip("\n").split("\t")
            mer_size = len(mer)
            try:
                mer_query.add_mer(mer, int(id), int(loci))
            except ValueError:
                continue
    
    confi_reads = []
    cache_size = 5000
    
    while True:
        try:
            chunk = queue.get(timeout=60)
            if chunk is None:
                break
                
            for seq_id, seq, qua in chunk:
                confi_id, strand = validate_read(
                    seq, 
                    mer_query, 
                    mer_size, 
                    5
                    )
                if confi_id > -1:
                    confi_reads.append((seq_id, strand, confi_id, seq))
                
                if len(confi_reads) >= cache_size:
                    out_write(output_file, confi_reads)
                    confi_reads = []
                    
        except Empty:
            break
            
    if confi_reads:
        out_write(output_file, confi_reads)


def out_write(file, infos):
    """
    Write assembled results to output file.
    
    Args:
        file (str): Output file path
        infos (list): List of results to write
    """
    with open(file, "a") as f:
        f.write("\n".join(["\t".join([str(ele) for ele in info]) for info in infos]))
        f.write("\n")


def progress_monitor(progress_queue, total_depth):
    """
    Monitor and display progress bar.
    
    Args:
        progress_queue (Queue): Queue containing progress updates
        total_depth (int): Total number of reads to process
    """
    pbar = tqdm.tqdm(total=total_depth, desc="Processing reads")
    while True:
        count = progress_queue.get()
        if count is None:
            break
        pbar.update(count)
    pbar.close()


def main():
    parser = argparse.ArgumentParser(
        description="Assemble element based on confident mers."
    )
    parser.add_argument("reads_files", nargs="+", help="Input reads files (gzipped or plain text)")
    parser.add_argument("--mers", type=str, required=True, help="Mers table file")
    parser.add_argument(
        "--depth", type=int, default=200_000_000, help="Number of reads to assemble"
    )
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use")
    parser.add_argument("--output_dir", type=str, default="out", help="Output directory")
    parser.add_argument("--chunk_size", type=int, default=10000, help="Chunk size for reading")
    args = parser.parse_args()

    # clear output_dir
    shutil.rmtree(args.output_dir) if os.path.exists(args.output_dir) else None
    os.makedirs(args.output_dir)
    
    # Initialize queues for data and progress
    data_queue = mp.Queue(maxsize=200)
    progress_queue = mp.Queue()
    
    # Start progress monitor
    monitor = mp.Process(
        target=progress_monitor,
        args=(progress_queue, args.depth)
    )
    monitor.start()
    
    # Start reader process
    reader = mp.Process(
        target=reader_process,
        args=(args.reads_files, data_queue, args.chunk_size, args.depth, progress_queue)
    )
    reader.start()
    
    # Start worker processes
    workers = []
    for i in range(args.threads):
        output_file = os.path.join(args.output_dir, f"Assemble.{i}.reads")
        open(output_file, "w").close()  # Clear output file
        
        p = mp.Process(
            target=worker_process,
            args=(data_queue, args.mers, output_file, i)
        )
        workers.append(p)
        p.start()
    
    # Wait for all processes to complete
    reader.join()
    for w in workers:
        w.join()
    monitor.join()


if __name__ == "__main__":
    main()
