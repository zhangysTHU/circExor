#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
kmer_assemble_to_meme.py
------------------------
Input:
  --kmers  : k-mer FASTA (EV or Cyto)
  --circ   : circRNA FASTA (corresponding EV or Cyto)
Output:
  --out    : FASTA suitable for MEME (default windows centered on k-mer ±10 nt)

Two fragment generation modes (optional):
  1) window  (default): for each k-mer match in circRNA, take window around its center
  2) coverage: reuse RNAlight's pileup→connect→long-run idea, output covered fragments with length >= min_len
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
    np = None  # coverage mode requires numpy; window mode does not


# -------------------------
#      General utilities
# -------------------------

def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Simple FASTA reader: yields (seq_id, seq) with multi-line sequences concatenated."""
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
    """Write (seq_id, seq) iterable to FASTA."""
    with open(path, 'w', encoding='utf-8') as out:
        for sid, seq in records:
            out.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")


def revcomp(seq: str) -> str:
    """Reverse complement for DNA (if needed)."""
    table = str.maketrans("ACGTNacgtnUu", "TGCANtgcanTT")
    return seq.translate(table)[::-1]


def normalize_dna(seq: str, u2t: bool = True, upper: bool = True) -> str:
    """Optionally convert U->T and unify case."""
    if u2t:
        seq = seq.replace('U', 'T').replace('u', 't')
    return seq.upper() if upper else seq


# -------------------------
#   Assembly functions (with bug fixes)
# -------------------------

def generate_all_kmers(k: int, ignore_N: bool = True) -> "collections.OrderedDict[str,int]":
    """Generate mapping of all k-mers of length k to indices."""
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
    Given a sequence and query k-mers, return coverage matrix with shape (len(query_kmers), len(seq)).
    Each cell counts how many times the k-mer covers that base (overlaps accumulate).
    Requires numpy; not used in window mode.
    """
    if np is None:
        raise RuntimeError("sequence_kmer_pileup requires numpy; install numpy or use --mode window (default).")
    assert isinstance(query_kmers, list)
    lengths = set(len(k) for k in query_kmers)
    retval = np.zeros((len(query_kmers), len(seq)), dtype=int)
    for length in lengths:
        assert length <= len(seq), "Cannot query a kmer against a seq shorter than that kmer"
        kmers = [seq[i:i+length] for i in range(len(seq) - length + 1)]
        kmer_to_idx = generate_all_kmers(length)
        kmers_int = np.array([kmer_to_idx[k] for k in kmers if "N" not in k], dtype=int)
        query_int = np.atleast_2d(np.array([kmer_to_idx[k] for k in query_kmers if len(k) == length and "N" not in k], dtype=int)).T
        hits = np.where(query_int == kmers_int)  # broadcasting
        this_rows = np.zeros((len(query_int), len(seq)))
        for i in range(length):
            this_rows[hits[0], hits[1] + i] += 1
        retval_idx = np.array([i for i, k in enumerate(query_kmers) if len(k) == length], dtype=int)
        retval[retval_idx, ] = this_rows
    return retval


def connect_nearby_runs(pileup_flat, allowed_gap_num: int):
    """
    Fill zero runs of length <= allowed_gap_num between ones to connect nearby runs.
    Accepts list or np.array.
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
    Return all continuous runs of 1 with length > l as (start, length) using base-level coordinates.

    Fix (Bug #1):
    The old logic returned run block indices instead of base positions, causing wrong slice start.
    Below is the corrected version; the original buggy version is retained as comment.
    """
    # ---------- Corrected: return base-level (start, length) ----------
    retval = []
    pos = 0  # current base position
    for val, group in itertools.groupby(num_sequence):
        g = list(group)
        length = len(g)
        if val and length > l:
            retval.append((pos, length))
        pos += length
    return retval

    # ---------- Original (BUGGY): returned run block indices ----------
    # chunked = [(k, list(g)) for k, g in itertools.groupby(num_sequence)]
    # retval = [(i, len(g)) for i, (k, g) in enumerate(chunked) if k and len(g) > l]
    # return retval


def assemble_kmer_motifs_by_coverage(seq: str, kmers: List[str], min_len: int = 10, gap_allowed: int = 2) -> List[str]:
    """
    Using pileup->connect->long-run idea to output fragments with length >= min_len from a sequence.
    This is the 'coverage' assembly mode complementary to the window mode.
    """
    if np is None:
        raise RuntimeError("assemble_kmer_motifs_by_coverage requires numpy; install numpy or use --mode window.")
    try:
        pileup = sequence_kmer_pileup(seq, kmers)
    except AssertionError:
        return []
    pileup_flat = (pileup.sum(axis=0) > 0).astype(int)
    pileup_flat = connect_nearby_runs(pileup_flat, gap_allowed)
    motif_spans = find_long_runs(pileup_flat, l=min_len)
    ret = [seq[start:start+length] for (start, length) in motif_spans]
    # sanity check
    assert all(len(s) == length for s, (_, length) in zip(ret, motif_spans))
    return ret


# MSA helper: extract motifs from MSA (with Bug #2 fix)

def _fetch_kmer_from_msa_i(i: int, seed_seq: str, msa: List[str], min_len: int, min_reps: int) -> str:
    """
    Given MSA start i and seed_seq of length min_len, try to extend to the right as far as possible
    under the constraint of no gaps and at least min_reps sequences agreeing.
    Return the maximal consistent extension.

    Fix (Bug #2):
    The old logic used sorted(...)[0] which could pick the shortest/least-supported fragment,
    contrary to the goal of maximal extension. Now choose by (length, support).
    """
    relevant_seqs = [m[i:] for m in msa if m[i:i+min_len] == seed_seq]
    if not relevant_seqs:
        return seed_seq
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

    # ---------- Choose the extension with maximal length and highest support ----------
    def score(s: str) -> Tuple[int, int]:
        support = sum(1 for m in relevant_seqs if s in m)
        return (len(s), support)
    best = max(extended, key=score)
    return best

    # ---------- Original (BUGGY): could pick the shortest ----------
    # extended_properties = [(len(s), len([m for m in relevant_seqs if s in m])) for s in extended]
    # extended_sorted = [seq for _prop, seq in sorted(zip(extended_properties, extended))]
    # return extended_sorted[0]


def find_motifs_in_msa(msa: List[str], min_len: int = 7, min_reps: int = 3) -> List[str]:
    """
    Scan an MSA for seeds of length min_len that appear in at least min_reps sequences without gaps,
    extend them via _fetch_kmer_from_msa_i, and remove sequences contained within longer ones.
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
    # remove contained shorter sequences
    hits_extended.sort(key=len, reverse=True)
    dedup: List[str] = []
    for h in hits_extended:
        if any(h in longer for longer in dedup):
            continue
        dedup.append(h)
    return dedup


# -------------------------
#     Two fragment generation modes
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
    Window mode: search each circRNA for each k-mer; take a window around the k-mer 'center' (default ±flank).
    - For even-length k-mers, the center is the left-middle (floor).
    - If window exceeds sequence ends: pad_with_N=True will pad with N; otherwise the window is truncated
      (MEME accepts variable-length sequences).
    """
    out: List[Tuple[str, str]] = []
    seen: Set[str] = set()
    rc_cache: Dict[str, str] = {}
    for sid, seq in circ_records:
        L = len(seq)
        seq_upper = seq  # normalized upstream
        for k in kmers:
            klen = len(k)
            # forward strand
            start = 0
            while True:
                idx = seq_upper.find(k, start)
                if idx == -1:
                    break
                center = idx + (klen // 2)
                left = center - flank
                right = center + flank + 1  # inclusive center => length = 2*flank+1
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
                start = idx + 1  # allow overlapping matches

            # reverse complement (optional)
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
    Coverage mode: for each circRNA use assemble_kmer_motifs_by_coverage to connect nearby matches
    (allowing gaps) into longer fragments.
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
#           Main
# -------------------------

def main():
    p = argparse.ArgumentParser(
        description="Generate fragments for MEME from k-mers (default: windows centered on k-mer ±flank). Coverage mode available."
    )
    p.add_argument("--kmers", required=True, help="k-mer FASTA file (EV/Cyto)")
    p.add_argument("--circ",  required=True, help="circRNA FASTA file (EV/Cyto)")
    p.add_argument("--out",   required=True, help="Output FASTA (suitable for MEME)")

    # window mode parameters
    p.add_argument("--flank", type=int, default=10, help="Window flank size on each side of k-mer center (default 10)")
    p.add_argument("--both_strands", action="store_true", help="Also search the reverse complement (default: off)")
    p.add_argument("--pad_with_N", action="store_true", help="Pad out-of-bounds windows with N (default: truncate)")

    # coverage mode parameters
    p.add_argument("--mode", choices=["window", "coverage"], default="window", help="Fragment generation mode (default: window)")
    p.add_argument("--min_len", type=int, default=10, help="Coverage mode: minimum fragment length (default 10)")
    p.add_argument("--gap_allowed", type=int, default=2, help="Coverage mode: max allowed zero-gap to connect runs (default 2)")

    # common
    p.add_argument("--dedup", action="store_true", help="Deduplicate identical sequences (keep one)")
    p.add_argument("--keep_case", action="store_true", help="Preserve case (default: convert to uppercase)")
    p.add_argument("--keep_u", action="store_true", help="Keep U (default: convert U->T)")

    args = p.parse_args()

    # Read and normalize k-mers
    kmers: List[str] = []
    for kid, kseq in read_fasta(args.kmers):
        k = normalize_dna(kseq, u2t=(not args.keep_u), upper=(not args.keep_case))
        if "N" in k:
            continue
        if any(base not in "ACGT" for base in k):
            continue  # skip non-ACGT k-mers
        kmers.append(k)
    if not kmers:
        print("ERROR: No valid k-mers read (they may all contain N).", file=sys.stderr)
        sys.exit(1)

    # Read and normalize circRNA sequences
    circ_records = []
    for sid, s in read_fasta(args.circ):
        seq = normalize_dna(s, u2t=(not args.keep_u), upper=(not args.keep_case))
        circ_records.append((sid, seq))
    if not circ_records:
        print("ERROR: No circRNA sequences read.", file=sys.stderr)
        sys.exit(1)

    # Generate fragments
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
        print("WARNING: No fragments were generated. Check k-mers vs sequences or relax parameters.", file=sys.stderr)

    # Write FASTA
    write_fasta(frags, args.out)
    print(f"[OK] Fragments written: {args.out} ({len(frags)} entries)")


if __name__ == "__main__":
    main()
