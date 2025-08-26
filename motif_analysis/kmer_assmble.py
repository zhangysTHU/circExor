#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
kmer_assemble_to_meme.py
------------------------
输入：
  --kmers  : k-mer 的 FASTA（EV 或 Cyto 均可）
  --circ   : circRNA 的 FASTA（对应 EV 或 Cyto，均可）
输出：
  --out    : 可直接用于 MEME 的 FASTA（默认基于 k-mer 中心 ±10 nt 的窗口片段）

两种“拼接/片段生成”模式（可选）：
  1) window（默认）：每当在 circRNA 中命中某个 k-mer，取其“中心”±flank 的窗口
  2) coverage      ：复用 RNAlight 的 pileup→连通→长 run 思路，输出长度≥min_len 的覆盖片段

"""

import argparse
import os
import sys
import itertools
import collections
from typing import List, Tuple, Dict, Iterable, Iterator, Optional, Set

try:
    import numpy as np
except ImportError:
    np = None  # coverage 模式需要 numpy；window 模式不需要


# -------------------------
#      通用工具函数
# -------------------------

def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """简单 FASTA 读取器：返回 (seq_id, seq)；seq 自动合并多行。"""
    sid, buf = None, []
    with open(path, 'r', encoding='utf-8') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if sid is not None:
                    yield sid, ''.join(buf)
                sid = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if sid is not None:
        yield sid, ''.join(buf)


def write_fasta(records: Iterable[Tuple[str, str]], path: str) -> None:
    """把 (seq_id, seq) 可迭代写成 FASTA。"""
    with open(path, 'w', encoding='utf-8') as out:
        for sid, seq in records:
            out.write(f">{sid}\n")
            # 切行更友好（可选）
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def revcomp(seq: str) -> str:
    """DNA 反向互补（若需要）"""
    table = str.maketrans("ACGTNacgtnUu", "TGCANtgcanTT")
    return seq.translate(table)[::-1]


def normalize_dna(seq: str, u2t: bool = True, upper: bool = True) -> str:
    """U->T（可选）并统一大小写"""
    if u2t:
        seq = seq.replace('U', 'T').replace('u', 't')
    return seq.upper() if upper else seq


# -------------------------
#   组装函数（含 bug 修复）
# -------------------------

def generate_all_kmers(k: int, ignore_N: bool = True) -> "collections.OrderedDict[str,int]":
    """
    生成长度 k 的所有 k-mer 映射（字符串 -> 索引）。
    """
    alphabet = "ACGT"
    if not ignore_N:
        alphabet += "N"
    possible = itertools.product(alphabet, repeat=k)
    ret = collections.OrderedDict()
    for i, kmer in enumerate(possible):
        ret[''.join(kmer)] = i
    return ret


def sequence_kmer_pileup(seq: str, query_kmers: List[str]):
    """
    给定序列和查询 k-mer，返回“覆盖矩阵”：shape = (len(query_kmers), len(seq))
    单元值表示该 k-mer 覆盖该碱基位置的次数（重叠会累计）。
    需要 numpy；在 window 模式不会用到此函数。
    """
    if np is None:
        raise RuntimeError("sequence_kmer_pileup 需要 numpy，请安装 numpy 或改用 --mode window（默认）。")
    assert isinstance(query_kmers, list)
    lengths = set(len(k) for k in query_kmers)
    retval = np.zeros((len(query_kmers), len(seq)), dtype=int)
    for length in lengths:
        assert length <= len(seq), "Cannot query a kmer against a seq shorter than that kmer"
        kmers = [seq[i:i+length] for i in range(len(seq) - length + 1)]
        kmer_to_idx = generate_all_kmers(length)
        kmers_int = np.array([kmer_to_idx[k] for k in kmers if "N" not in k], dtype=int)
        query_int = np.atleast_2d(np.array([kmer_to_idx[k] for k in query_kmers if len(k) == length and "N" not in k], dtype=int)).T
        hits = np.where(query_int == kmers_int)  # 自动广播
        this_rows = np.zeros((len(query_int), len(seq)))
        for i in range(length):
            this_rows[hits[0], hits[1] + i] += 1
        retval_idx = np.array([i for i, k in enumerate(query_kmers) if len(k) == length], dtype=int)
        retval[retval_idx, ] = this_rows
    return retval


def connect_nearby_runs(pileup_flat, allowed_gap_num: int):
    """
    把两段 1-run 之间长度<=allowed_gap_num 的 0 段填成 1，实现连接。
    支持 list/np.array。
    """
    arr = list(pileup_flat)
    chunked = [(k, list(g)) for k, g in itertools.groupby(arr)]
    retval = []
    for i, (item, group) in enumerate(chunked):
        if not item and len(group) <= allowed_gap_num and 0 < i < len(chunked) - 1:
            retval.extend([1] * len(group))
        else:
            retval.extend(group)
    if np is not None and isinstance(pileup_flat, np.ndarray):
        return np.array(retval, dtype=int)
    return retval


def find_long_runs(num_sequence: Iterable[int], l: int) -> List[Tuple[int, int]]:
    """
    返回所有“连续 1-run 且长度 > l”的 (start, length)（碱基级坐标）

    ===== 修复点（Bug #1）=====
    旧逻辑返回的是“run 的块序号”，不是“碱基起点”，会导致后续 seq[i:i+l] 起点错误。
    下面给出【修正版】并保留【原版（BUGGY）】为注释，方便快速切换。
    """
    # ---------- 【修正版：返回碱基级 (start, length)】 ----------
    retval = []
    pos = 0  # 当前遍历到的碱基位置
    for val, group in itertools.groupby(num_sequence):
        g = list(group)
        length = len(g)
        if val and length > l:
            retval.append((pos, length))
        pos += length
    return retval

    # ---------- 【原版（BUGGY）：返回 run 的块序号】 ----------
    # chunked = [(k, list(g)) for k, g in itertools.groupby(num_sequence)]
    # retval = [(i, len(g)) for i, (k, g) in enumerate(chunked) if k and len(g) > l]
    # return retval


def assemble_kmer_motifs_by_coverage(seq: str, kmers: List[str], min_len: int = 10, gap_allowed: int = 2) -> List[str]:
    """
    利用覆盖矩阵→连通→长 run 的思路，从一条序列中输出长度>=min_len 的候选片段。
    这是 RNAlight 思路的“覆盖拼接模式”，与 window 模式互补。
    """
    if np is None:
        raise RuntimeError("assemble_kmer_motifs_by_coverage 需要 numpy，请安装 numpy 或改用 --mode window。")
    try:
        pileup = sequence_kmer_pileup(seq, kmers)
    except AssertionError:
        return []
    pileup_flat = (pileup.sum(axis=0) > 0).astype(int)  # 是否被任意 k-mer 覆盖的 0/1 串
    pileup_flat = connect_nearby_runs(pileup_flat, gap_allowed)
    motif_spans = find_long_runs(pileup_flat, l=min_len)
    ret = [seq[start:start+length] for (start, length) in motif_spans]
    # 保守断言
    assert all(len(s) == length for s, (_, length) in zip(ret, motif_spans))
    return ret


# ====== MSA 辅助：从 MSA 提 motif（含 Bug #2 修正） ======

def _fetch_kmer_from_msa_i(i: int, seed_seq: str, msa: List[str], min_len: int, min_reps: int) -> str:
    """
    给定 MSA 的起点 i 和长度为 min_len 的 seed_seq，尝试在不含 gap 且有 >=min_reps 条序列一致的前提下尽可能向右延伸。
    返回延伸后的最大一致片段。

    ===== 修复点（Bug #2）=====
    旧逻辑 sorted(...)[0] 会拿到“最短/支持少”的片段，违背“尽量延伸”的意图。
    这里改为按 (长度, 支持数) 选最大。
    """
    # 筛出相关序列并去掉 i 之前的前缀
    relevant_seqs = [m[i:] for m in msa if m[i:i+min_len] == seed_seq]
    if not relevant_seqs:
        return seed_seq
    # 尝试所有满足 min_reps 的组合，逐列比较一致性，累积最远一致长度
    extended: List[str] = []
    for combo in itertools.combinations(relevant_seqs, min_reps):
        this_seq: List[str] = []
        L = len(combo[0])
        for j in range(L):
            jth_bases = {c[j] for c in combo}
            if (len(jth_bases) != 1) or ('-' in jth_bases):
                break
            this_seq.append(jth_bases.pop())
        extended.append(''.join(this_seq))
    if not extended:
        return seed_seq

    # ---------- 【修正版：选择“最长、且支持数最多”的延伸】 ----------
    def score(s: str) -> Tuple[int, int]:
        support = sum(1 for m in relevant_seqs if s in m)
        return (len(s), support)
    best = max(extended, key=score)
    return best

    # ---------- 【原版（BUGGY）：会选到最短】 ----------
    # extended_properties = [(len(s), len([m for m in relevant_seqs if s in m])) for s in extended]
    # extended_sorted = [seq for _prop, seq in sorted(zip(extended_properties, extended))]
    # return extended_sorted[0]


def find_motifs_in_msa(msa: List[str], min_len: int = 7, min_reps: int = 3) -> List[str]:
    """
    在 MSA 上用长度为 min_len 的窗口找“至少有 min_reps 条序列一致且无 gap”的种子，再调用 _fetch_kmer_from_msa_i 延伸，并去掉被包含的短序列。
    """
    unique_lens = set(len(m) for m in msa)
    assert len(unique_lens) == 1, "All MSA sequences must be of the same length"
    msa_len = unique_lens.pop()

    hits: List[Tuple[int, str]] = []
    for i in range(msa_len - min_len):
        block = [m[i:i+min_len] for m in msa]
        block_no_gaps = [m for m in block if "-" not in m]
        if not block_no_gaps:
            continue
        counter = collections.Counter(block_no_gaps)
        for kmer, count in counter.items():
            if count >= min_reps:
                hits.append((i, kmer))

    hits_extended = [_fetch_kmer_from_msa_i(i, kmer, msa=msa, min_len=min_len, min_reps=min_reps)
                     for (i, kmer) in hits]
    # 去除被包含的短序列
    hits_extended.sort(key=len, reverse=True)
    dedup: List[str] = []
    for h in hits_extended:
        if any(h in longer for longer in dedup):
            continue
        dedup.append(h)
    return dedup


# -------------------------
#     片段生成两种模式
# -------------------------

def windows_around_kmer_hits(
    circ_records: List[Tuple[str, str]],
    kmers: List[str],
    flank: int = 10,
    both_strands: bool = False,
    pad_with_N: bool = False,
    dedup: bool = True
) -> List[Tuple[str, str]]:
    """
    “窗口模式”：在每条 circRNA 里查找每个 k-mer 的命中，围绕命中 k-mer 的“中心”截取窗口（默认 ±10 nt）。
    - even 长度 k-mer 的“中心”取左中位（floor）。
    - near-end 的窗口若越界：pad_with_N=True 则补 N；否则截短（MEME 可接受不同长度）。
    """
    out: List[Tuple[str, str]] = []
    seen: Set[str] = set()
    rc_cache: Dict[str, str] = {}
    for sid, seq in circ_records:
        L = len(seq)
        seq_upper = seq  # 已在上游 normalize
        for k in kmers:
            klen = len(k)
            # 正向
            start = 0
            while True:
                idx = seq_upper.find(k, start)
                if idx == -1:
                    break
                center = idx + (klen // 2)
                left = center - flank
                right = center + flank + 1  # 包含中心 => 长度=2*flank+1
                if left < 0 or right > L:
                    if pad_with_N:
                        lpad = max(0, -left)
                        rpad = max(0, right - L)
                        sub = seq_upper[max(0, left):min(L, right)]
                        sub = ("N" * lpad) + sub + ("N" * rpad)
                    else:
                        sub = seq_upper[max(0, left):min(L, right)]
                else:
                    sub = seq_upper[left:right]
                header = f"{sid}|pos={idx}|kmer={k}|strand=+"
                if (not dedup) or (sub not in seen):
                    out.append((header, sub))
                    if dedup:
                        seen.add(sub)
                start = idx + 1  # 允许重叠命中

            # 反向互补（可选）
            if both_strands:
                if k in rc_cache:
                    k_rc = rc_cache[k]
                else:
                    k_rc = revcomp(k)
                    rc_cache[k] = k_rc
                start = 0
                while True:
                    idx = seq_upper.find(k_rc, start)
                    if idx == -1:
                        break
                    center = idx + (klen // 2)
                    left = center - flank
                    right = center + flank + 1
                    if left < 0 or right > L:
                        if pad_with_N:
                            lpad = max(0, -left)
                            rpad = max(0, right - L)
                            sub = seq_upper[max(0, left):min(L, right)]
                            sub = ("N" * lpad) + sub + ("N" * rpad)
                        else:
                            sub = seq_upper[max(0, left):min(L, right)]
                    else:
                        sub = seq_upper[left:right]
                    header = f"{sid}|pos={idx}|kmer={k}|strand=-"
                    if (not dedup) or (sub not in seen):
                        out.append((header, sub))
                        if dedup:
                            seen.add(sub)
                    start = idx + 1
    return out


def coverage_assembled_fragments(
    circ_records: List[Tuple[str, str]],
    kmers: List[str],
    min_len: int = 10,
    gap_allowed: int = 2,
    dedup: bool = True
) -> List[Tuple[str, str]]:
    """
    “覆盖模式”：对每条 circRNA，使用 assemble_kmer_motifs_by_coverage 把相邻命中（允许 gap）连成更长片段。
    """
    out: List[Tuple[str, str]] = []
    seen: Set[str] = set()
    for sid, seq in circ_records:
        frags = assemble_kmer_motifs_by_coverage(seq, kmers, min_len=min_len, gap_allowed=gap_allowed)
        for i, frag in enumerate(frags):
            header = f"{sid}|assembled={i}|len={len(frag)}"
            if (not dedup) or (frag not in seen):
                out.append((header, frag))
                if dedup:
                    seen.add(frag)
    return out


# -------------------------
#           主流程
# -------------------------

def main():
    p = argparse.ArgumentParser(
        description="基于 k-mer 生成可用于 MEME 的片段（默认：k-mer 中心 ±flank）。提供可选覆盖拼接模式；MSA 步骤已写好但默认注释。"
    )
    p.add_argument("--kmers", required=True, help="k-mer FASTA 文件（EV/Cyto 均可）")
    p.add_argument("--circ",  required=True, help="circRNA FASTA 文件（EV/Cyto 均可）")
    p.add_argument("--out",   required=True, help="输出 FASTA（可直接用于 MEME）")

    # 窗口模式参数
    p.add_argument("--flank", type=int, default=10, help="窗口模式：k-mer 中心两侧碱基数（默认 10，对应 ±10nt）")
    p.add_argument("--both_strands", action="store_true", help="也在反向互补上搜（默认否）")
    p.add_argument("--pad_with_N", action="store_true", help="越界时用 N 补齐固定长度窗口（默认截短即可）")

    # 覆盖模式参数
    p.add_argument("--mode", choices=["window", "coverage"], default="window", help="片段生成模式（默认 window）")
    p.add_argument("--min_len", type=int, default=10, help="coverage 模式：最短拼接长度（默认 10）")
    p.add_argument("--gap_allowed", type=int, default=2, help="coverage 模式：允许连通的 0-gap 最大长度（默认 2）")

    # 通用
    p.add_argument("--dedup", action="store_true", help="去重（相同序列只保留一次）")
    p.add_argument("--keep_case", action="store_true", help="保留大小写（默认转大写）")
    p.add_argument("--keep_u", action="store_true", help="保留 U（默认 U->T）")

    args = p.parse_args()

    # 读取并规范化 k-mers
    kmers: List[str] = []
    # for kid, kseq in read_fasta(args.kmers):
    #     k = normalize_dna(kseq, u2t=(not args.keep_u), upper=(not args.keep_case))
    #     if "N" in k:
    #         # N 会导致匹配不明确：这里直接跳过，也可改成用正则/模糊匹配
    #         continue
    #     kmers.append(k)
    for kid, kseq in read_fasta(args.kmers):
        k = normalize_dna(kseq, u2t=(not args.keep_u), upper=(not args.keep_case))
        if "N" in k:
            continue
        if any(base not in "ACGT" for base in k):
            continue  # 跳过非 ACGT k-mer
        kmers.append(k)
    if not kmers:
        print("ERROR: 未读取到有效 k-mer（可能都含 N）。", file=sys.stderr)
        sys.exit(1)

    # 读取并规范化 circRNA
    circ_records = []
    for sid, s in read_fasta(args.circ):
        seq = normalize_dna(s, u2t=(not args.keep_u), upper=(not args.keep_case))
        circ_records.append((sid, seq))
    if not circ_records:
        print("ERROR: 未读取到 circRNA 序列。", file=sys.stderr)
        sys.exit(1)

    # 生成片段
    if args.mode == "window":
        frags = windows_around_kmer_hits(
            circ_records, kmers,
            flank=args.flank,
            both_strands=args.both_strands,
            pad_with_N=args.pad_with_N,
            dedup=args.dedup
        )
    else:
        frags = coverage_assembled_fragments(
            circ_records, kmers,
            min_len=args.min_len,
            gap_allowed=args.gap_allowed,
            dedup=args.dedup
        )

    if not frags:
        print("WARNING: 未生成任何片段。请检查 k-mer 与序列是否匹配、或放宽参数。", file=sys.stderr)

    # 写出 FASTA
    write_fasta(frags, args.out)
    print(f"[OK] 片段已写出：{args.out}（{len(frags)} 条）")


if __name__ == "__main__":
    main()
