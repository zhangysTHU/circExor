#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
import itertools
import collections
import argparse
import subprocess

# ---------- 工具函数 ----------
def generate_all_kmers(k, ignore_N=True):
    alphabet = "ACGT"
    if not ignore_N:
        alphabet += "N"
    possible_kmers = itertools.product(alphabet, repeat=k)
    return collections.OrderedDict((''.join(kmer), i) for i, kmer in enumerate(possible_kmers))

def sequence_kmer_pileup(seq, query_kmers):
    assert isinstance(query_kmers, list)
    lengths = set([len(kmer) for kmer in query_kmers])
    retval = np.zeros((len(query_kmers), len(seq))).astype(int)
    for length in lengths:
        assert length <= len(seq), "Cannot query a kmer against a seq shorter than that kmer"
        kmers = [seq[i:i+length] for i in range(len(seq) - length + 1)]
        kmer_to_idx = generate_all_kmers(length)
        kmers_int = np.array([kmer_to_idx[k] for k in kmers if "N" not in k], dtype=int)
        query_int = np.atleast_2d(np.array([kmer_to_idx[k] for k in query_kmers if len(k) == length and "N" not in k], dtype=int)).T
        hits = np.where(query_int == kmers_int)
        this_rows = np.zeros((len(query_int), len(seq)))
        for i in range(length):
            this_rows[hits[0], hits[1] + i] += 1
        retval_idx = np.array([i for i, k in enumerate(query_kmers) if len(k) == length], dtype=int)
        retval[retval_idx, ] = this_rows
    return retval

def find_long_runs(num_sequence, l):
    chunked = [(k, list(g)) for k, g in itertools.groupby(num_sequence)]
    return [(i, len(g)) for i, (k, g) in enumerate(chunked) if k and len(g) > l]

def connect_nearby_runs(pileup_flat, allowed_gap_num):
    chunked = [(k, list(g)) for k, g in itertools.groupby(list(pileup_flat))]
    retval = []
    for i, (item, group) in enumerate(chunked):
        if not item and len(group) <= allowed_gap_num and 0 < i < len(chunked) - 1:
            retval.extend([1] * len(group))
        else:
            retval.extend(group)
    return np.array(retval, dtype=int)

def assemble_kmer_motifs(seq, kmers, min_len=10, gap_allowed=2):
    try:
        pileup = sequence_kmer_pileup(seq, kmers)
    except AssertionError:
        return []
    pileup_flat = np.clip(np.sum(pileup, axis=0), 0, 1)
    pileup_flat = connect_nearby_runs(pileup_flat, gap_allowed)
    motif_idx = find_long_runs(pileup_flat, l=min_len)
    retval = [seq[i:i+l] for i, l in motif_idx]
    assert all([len(s) == l for s, (_i, l) in zip(retval, motif_idx)])
    return retval

def read_msa(fname, combine_lines=True):
    retval = collections.defaultdict(list)
    with open(fname) as source:
        for line in source:
            if line.startswith("CLUSTAL") or not line.strip():
                continue
            tokens = line.strip().split()
            seq_id, aln = tokens
            retval[seq_id].append(aln)
    if combine_lines:
        retval = {k: ''.join(v) for k, v in retval.items()}
    return retval

def _fetch_kmer_from_msa_i(i, seed_seq, msa, min_len, min_reps):
    relevant_seqs = [m[i:] for m in msa if m[i:i+min_len] == seed_seq]
    if not relevant_seqs:
        return seed_seq
    extended = []
    for combo in itertools.combinations(relevant_seqs, min_reps):
        this_seq = []
        for j in range(len(combo[0])):
            jth_bases = set([c[j] for c in combo])
            if (len(jth_bases) != 1) or ('-' in jth_bases):
                break
            this_seq.append(jth_bases.pop())
        extended.append(''.join(this_seq))
    extended_properties = [(len(s), len([m for m in relevant_seqs if s in m])) for s in extended]
    extended_sorted = [seq for _prop, seq in sorted(zip(extended_properties, extended))]
    return extended_sorted[0] if extended_sorted else seed_seq

def find_motifs_in_msa(msa, min_len=7, min_reps=3):
    unique_lengths = set([len(m) for m in msa])
    assert len(unique_lengths) == 1, "All MSA sequences must be of the same length"
    msa_len = unique_lengths.pop()

    hits = []
    for i in range(msa_len - min_len):
        block = [m[i:i+min_len] for m in msa]
        block_no_gaps = [m for m in block if "-" not in m]
        if not block_no_gaps:
            continue
        block_counter = collections.Counter(block_no_gaps)
        for kmer, count in block_counter.items():
            if count >= min_reps:
                hits.append((i, kmer))

    hits_extended = [_fetch_kmer_from_msa_i(hit[0], hit[1], msa=msa, min_len=min_len, min_reps=min_reps) for hit in hits]
    hits_extended.sort(key=len, reverse=True)
    hits_extended_dedup = []
    for hit in hits_extended:
        if any([hit in longer for longer in hits_extended_dedup]):
            continue
        hits_extended_dedup.append(hit)
    return hits_extended_dedup

# def run_pipeline(input_fasta, output_dir, threads=8, min_len=7, min_reps=3):
#     os.makedirs(output_dir, exist_ok=True)

#     # Step 1: MSA
#     cluster_file = os.path.join(output_dir, "msa.cluster")
#     clustalo_cmd = [
#         "clustalo",
#         "--threads", str(threads),
#         "-i", input_fasta,
#         "-o", cluster_file,
#         "--outfmt", "clu"
#     ]
#     subprocess.run(clustalo_cmd, check=True)

#     # Step 2: MSA to motifs
#     msa_dict = read_msa(cluster_file, combine_lines=True)
#     msa_list = list(msa_dict.values())
#     motifs = find_motifs_in_msa(msa_list, min_len=min_len, min_reps=min_reps)

#     # 保存为FASTA
#     fasta_file = os.path.join(output_dir, "motifs.fa")
#     with open(fasta_file, "w") as f:
#         for i, motif in enumerate(motifs, start=1):
#             f.write(f">motif{i}\n{motif}\n")

#     print(f"[OK] MSA saved to: {cluster_file}")
#     print(f"[OK] Motifs saved to: {fasta_file}")
    
# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="FASTA -> MSA -> motifs pipeline")
#     parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
#     parser.add_argument("-o", "--output", required=True, help="Output directory")
#     parser.add_argument("--threads", type=int, default=8, help="Number of threads for Clustal Omega")
#     parser.add_argument("--min_len", type=int, default=7, help="Minimum motif length")
#     parser.add_argument("--min_reps", type=int, default=3, help="Minimum repetitions in MSA")
#     args = parser.parse_args()

#     run_pipeline(args.input, args.output, args.threads, args.min_len, args.min_reps)

def run_pipeline(input_fasta, output_dir, threads=8, min_len=7, min_reps=3, force=False):
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: MSA
    cluster_file = os.path.join(output_dir, "msa.cluster")
    
    # 检查文件是否存在，如果存在且未使用force参数则报错
    if os.path.exists(cluster_file) and not force:
        raise FileExistsError(f"FATAL: Cowardly refusing to overwrite already existing file '{cluster_file}'. Use --force to force overwriting.")
    
    clustalo_cmd = [
        "clustalo",
        "--threads", str(threads),
        "-i", input_fasta,
        "-o", cluster_file,
        "--outfmt", "clu",
        "--force"  # 添加clustalo的force参数
    ]
    subprocess.run(clustalo_cmd, check=True)

    # Step 2: MSA to motifs
    msa_dict = read_msa(cluster_file, combine_lines=True)
    msa_list = list(msa_dict.values())
    motifs = find_motifs_in_msa(msa_list, min_len=min_len, min_reps=min_reps)

    # 保存为FASTA
    fasta_file = os.path.join(output_dir, "motifs.fa")
    # 检查文件是否存在，如果存在且未使用force参数则报错
    if os.path.exists(fasta_file) and not force:
        raise FileExistsError(f"FATAL: Cowardly refusing to overwrite already existing file '{fasta_file}'. Use --force to force overwriting.")
    
    with open(fasta_file, "w") as f:
        for i, motif in enumerate(motifs, start=1):
            f.write(f">motif{i}\n{motif}\n")

    print(f"[OK] MSA saved to: {cluster_file}")
    print(f"[OK] Motifs saved to: {fasta_file}")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FASTA -> MSA -> motifs pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads for Clustal Omega")
    parser.add_argument("--min_len", type=int, default=7, help="Minimum motif length")
    parser.add_argument("--min_reps", type=int, default=3, help="Minimum repetitions in MSA")
    parser.add_argument("--force", action="store_true", help="Force overwrite existing files")
    args = parser.parse_args()

    run_pipeline(args.input, args.output, args.threads, args.min_len, args.min_reps, args.force)
