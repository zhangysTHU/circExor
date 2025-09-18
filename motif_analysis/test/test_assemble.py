#!/usr/bin/env python3
import os
import argparse
import collections
import itertools
import numpy as np
import pandas as pd


# ---------- 核心函数 ----------

def generate_all_kmers(k, ignore_N=True):
    alphabet = "ACGT"
    if not ignore_N:
        alphabet += "N"
    kmers = itertools.product(alphabet, repeat=k)
    return {''.join(kmer): i for i, kmer in enumerate(kmers)}


def sequence_kmer_pileup(seq, query_kmers):
    """返回 pileup 矩阵：每个 kmer 在序列中的位置"""
    assert isinstance(query_kmers, list)
    lengths = set([len(k) for k in query_kmers])
    retval = np.zeros((len(query_kmers), len(seq))).astype(int)

    for length in lengths:
        if length > len(seq):
            continue
        kmers = [seq[i:i+length] for i in range(len(seq)-length+1)]
        kmer_to_idx = generate_all_kmers(length)
        kmers_int = np.array([kmer_to_idx[k] for k in kmers if "N" not in k], dtype=int)
        query_int = np.atleast_2d(np.array(
            [kmer_to_idx[k] for k in query_kmers if len(k) == length and "N" not in k],
            dtype=int)).T
        hits = np.where(query_int == kmers_int)
        this_rows = np.zeros((len(query_int), len(seq)))
        for i in range(length):
            this_rows[hits[0], hits[1] + i] += 1
        retval_idx = np.array([i for i, k in enumerate(query_kmers) if len(k) == length], dtype=int)
        retval[retval_idx, ] = this_rows
    return retval


def connect_nearby_runs(pileup_flat, allowed_gap_num):
    """合并相邻的 motif 片段（允许一定gap）"""
    chunked = [(k, list(g)) for k, g in itertools.groupby(list(pileup_flat))]
    retval = []
    for i, (item, group) in enumerate(chunked):
        if not item and len(group) <= allowed_gap_num and 0 < i < len(chunked)-1:
            retval.extend([1]*len(group))
        else:
            retval.extend(group)
    return np.array(retval, dtype=int)


def find_long_runs(num_sequence, l):
    """找到连续片段"""
    chunked = [(k, list(g)) for k, g in itertools.groupby(num_sequence)]
    return [(i, len(g)) for i, (k, g) in enumerate(chunked) if k and len(g) > l]


def assemble_kmer_motifs(seq, kmers, min_len=10, gap_allowed=2):
    """从序列和kmers组装motifs"""
    try:
        pileup = sequence_kmer_pileup(seq, kmers)
    except AssertionError:
        return []
    pileup_flat = np.clip(np.sum(pileup, axis=0), 0, 1)
    pileup_flat = connect_nearby_runs(pileup_flat, gap_allowed)
    motif_idx = find_long_runs(pileup_flat, l=min_len)
    return [seq[i:i+l] for i, l in motif_idx]


# ---------- 主流程 ----------
def run_pipeline(kmer_file, seq_file, output_fasta, min_len=10, gap_allowed=2):
    # 读 kmer
    kmers = pd.read_csv(kmer_file, sep="\t", header=None)[0].tolist()
    # 读序列
    if seq_file.endswith(".fasta") or seq_file.endswith(".fa"):
        from Bio import SeqIO
        seqs = [str(record.seq) for record in SeqIO.parse(seq_file, "fasta")]
    else:  # 兼容 tsv 格式
        df = pd.read_csv(seq_file, sep="\t")
        seqs = df.iloc[:, -1].tolist()

    motifs = []
    for seq in seqs:
        motifs += assemble_kmer_motifs(seq, kmers, min_len=min_len, gap_allowed=gap_allowed)

    # 输出 fasta
    with open(output_fasta, "w") as f:
        for i, m in enumerate(motifs, 1):
            f.write(f">motif_{i}\n{m}\n")

    print(f"[OK] Saved motifs to {output_fasta}")


# ---------- CLI 接口 ----------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assemble kmers into motifs (cyto vs ev)")
    parser.add_argument("--kmer", required=True, help="Input kmer file (txt/tsv)")
    parser.add_argument("--seq", required=True, help="Input sequence file (fasta/tsv)")
    parser.add_argument("--out", required=True, help="Output fasta file")
    parser.add_argument("--min_len", type=int, default=10, help="Minimum motif length (default: 10)")
    parser.add_argument("--gap_allowed", type=int, default=2, help="Allowed gap length (default: 2)")

    args = parser.parse_args()
    run_pipeline(args.kmer, args.seq, args.out, args.min_len, args.gap_allowed)
