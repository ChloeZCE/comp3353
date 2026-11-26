"""COMP3353 A3 Q2 - progressive multiple sequence alignment.

Scores: match +3, mismatch -3, residue vs gap -2, gap vs gap 0.
A single sequence is treated as a profile with one row, so align_profiles()
is used for parts a, b and d.

Run with no arguments to do all parts on the course data, or pass a FASTA
file to get only its progressive MSA (part d).
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

MATCH, MISMATCH, GAP = 3, -3, -2
GAP_CHAR = "-"

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "Q2")
RESULTS = os.path.join(HERE, "..", "results")


def s(a, b):
    """Score for one pair of characters (residue or gap)."""
    if a == GAP_CHAR and b == GAP_CHAR:
        return 0
    if a == GAP_CHAR or b == GAP_CHAR:
        return GAP
    return MATCH if a == b else MISMATCH


def read_fasta(path):
    """Read a FASTA file into a list of (name, sequence).
    For UniProt headers like sp|P69905|HBA_HUMAN the name is HBA_HUMAN."""
    recs, name, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                tok = line[1:].split()[0]
                name = tok.split("|")[2] if tok.count("|") == 2 else tok
                buf = []
            else:
                buf.append(line.upper())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


# alignment DP

def column_score_matrix(A, B):
    """Average score S(C_i, C_j) for every column i of A and column j of B.
    Uses letter counts per column instead of looping over every pair."""
    alphabet = sorted(set("".join(A)) | set("".join(B)) | {GAP_CHAR})
    idx = {c: k for k, c in enumerate(alphabet)}
    M = np.array([[s(a, b) for b in alphabet] for a in alphabet], dtype=float)

    def counts(P):
        C = np.zeros((len(P[0]), len(alphabet)))
        for row in P:
            for i, ch in enumerate(row):
                C[i, idx[ch]] += 1
        return C

    cA, cB = counts(A), counts(B)
    pair = cA @ M @ cB.T / (len(A) * len(B))
    # score of each column against a new all-gap column
    gapA = cA @ M[:, idx[GAP_CHAR]] / len(A)
    gapB = cB @ M[:, idx[GAP_CHAR]] / len(B)
    return pair, gapA, gapB


def align_profiles(A, B):
    """Needleman-Wunsch alignment of two profiles (lists of aligned strings).
    Returns the merged rows (A's rows first, then B's) and the final score.
    A gap is always added to every row of a profile. Ties: diag > up > left.
    """
    pair, gapA, gapB = column_score_matrix(A, B)
    n, m = len(A[0]), len(B[0])
    D = np.zeros((n + 1, m + 1))
    T = np.zeros((n + 1, m + 1), dtype=np.int8)   # traceback: 0 diag, 1 up, 2 left
    D[1:, 0] = np.cumsum(gapA)
    D[0, 1:] = np.cumsum(gapB)
    T[1:, 0] = 1
    T[0, 1:] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = D[i - 1, j - 1] + pair[i - 1, j - 1]
            up = D[i - 1, j] + gapA[i - 1]
            left = D[i, j - 1] + gapB[j - 1]
            best = diag
            t = 0
            if up > best + 1e-9:
                best, t = up, 1
            if left > best + 1e-9:
                best, t = left, 2
            D[i, j] = best
            T[i, j] = t
    # trace back from the bottom-right corner
    i, j = n, m
    colsA, colsB = [], []
    while i > 0 or j > 0:
        t = T[i, j]
        if t == 0:
            colsA.append(i - 1); colsB.append(j - 1); i -= 1; j -= 1
        elif t == 1:
            colsA.append(i - 1); colsB.append(None); i -= 1
        else:
            colsA.append(None); colsB.append(j - 1); j -= 1
    colsA.reverse(); colsB.reverse()
    newA = ["".join(GAP_CHAR if c is None else r[c] for c in colsA) for r in A]
    newB = ["".join(GAP_CHAR if c is None else r[c] for c in colsB) for r in B]
    return newA + newB, D[n, m]


# part a
def needleman_wunsch(x, y):
    (ax, ay), score = align_profiles([x], [y])
    return ax, ay, score


# part b
def align_seq_to_profile(profile, q):
    rows, score = align_profiles(profile, [q])
    return rows, score


# part c
def identity_distance(ax, ay):
    same = sum(1 for a, b in zip(ax, ay) if a == b and a != GAP_CHAR)
    return 1 - same / len(ax)


def distance_matrix(recs):
    n = len(recs)
    Dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            ax, ay, _ = needleman_wunsch(recs[i][1], recs[j][1])
            Dm[i, j] = Dm[j, i] = identity_distance(ax, ay)
    return Dm


def upgma(names, Dm):
    """UPGMA. Returns the merges in order as (cluster1, cluster2, distance),
    where a cluster is a tuple of sequence names."""
    pos = {n: k for k, n in enumerate(names)}
    clusters = [(n,) for n in names]
    merges = []

    def cdist(A, B):
        return np.mean([Dm[pos[a], pos[b]] for a in A for b in B])

    while len(clusters) > 1:
        best = None
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = cdist(clusters[i], clusters[j])
                if best is None or d < best[0] - 1e-12:
                    best = (d, i, j)
        d, i, j = best
        A, B = clusters[i], clusters[j]
        merges.append((A, B, d))
        clusters = [c for k, c in enumerate(clusters) if k not in (i, j)]
        clusters.append(A + B)
    return merges


# part d
def progressive_msa(recs, merges):
    seqs = dict(recs)
    profiles = {(n,): [seqs[n]] for n in seqs}   # each sequence starts as its own profile
    for A, B, _ in merges:
        rows, _ = align_profiles(profiles.pop(A), profiles.pop(B))
        profiles[A + B] = rows
    (cluster, rows), = profiles.items()
    return list(zip(cluster, rows))


# printing and plots
def fmt_block(names, rows, width=None):
    w = max(len(n) for n in names)
    if width is None:
        return "\n".join(f"{n:<{w}}  {r}" for n, r in zip(names, rows))
    out = []
    L = len(rows[0])
    for k in range(0, L, width):
        out.append("\n".join(f"{n:<{w}}  {r[k:k + width]}" for n, r in zip(names, rows)))
        out.append("")
    return "\n".join(out).rstrip()


def fmt_merges(merges):
    return "\n".join(f"  step {k}: ({A!r}, {B!r})   d = {d:.4f}"
                     for k, (A, B, d) in enumerate(merges, 1))


def heatmap(names, Dm, path, title):
    fig, ax = plt.subplots(figsize=(0.55 * len(names) + 2, 0.55 * len(names) + 1.5))
    im = ax.imshow(Dm, cmap="viridis")
    ax.set_xticks(range(len(names)), names, rotation=90 if len(names) > 8 else 0)
    ax.set_yticks(range(len(names)), names)
    if len(names) <= 10:
        for i in range(len(names)):
            for j in range(len(names)):
                ax.text(j, i, f"{Dm[i, j]:.3f}", ha="center", va="center",
                        color="w" if Dm[i, j] < Dm.max() * 0.6 else "k", fontsize=9)
    fig.colorbar(im, ax=ax, label="d = 1 - identity")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def msa_report(path, label):
    recs = read_fasta(path)
    names = [n for n, _ in recs]
    Dm = distance_matrix(recs)
    merges = upgma(names, Dm)
    msa = progressive_msa(recs, merges)
    ok = all(r.replace(GAP_CHAR, "") == dict(recs)[n] for n, r in msa)
    out = [f"Progressive MSA of {label} ({len(recs)} sequences, "
           f"alignment length {len(msa[0][1])}; residues preserved: {ok})",
           "UPGMA merge order:", fmt_merges(merges), "",
           fmt_block([n for n, _ in msa], [r for _, r in msa], width=60)]
    with open(os.path.join(RESULTS, f"q2d_{label}_msa.fasta"), "w") as fh:
        for n, r in msa:
            fh.write(f">{n}\n{r}\n")
    return "\n".join(out), names, Dm


def main():
    os.makedirs(RESULTS, exist_ok=True)
    cyt = read_fasta(os.path.join(DATA, "cytochrome_c_unaligned.fa"))
    seq = dict(cyt)
    rep = []

    # (a)
    rep.append("=" * 70 + "\n(a) Needleman-Wunsch global alignment\n" + "=" * 70)
    for other in ("S4", "S6"):
        ax, ay, sc = needleman_wunsch(seq["S1"], seq[other])
        rep.append(f"S1 vs {other}:  score = {sc:g}\n{fmt_block(['S1', other], [ax, ay])}\n")

    # (b)
    rep.append("=" * 70 + "\n(b) Sequence-to-profile alignment\n" + "=" * 70)
    ex, _ = align_seq_to_profile(["-KTAYIAKQRQISFVKSHFSRQDIL",
                                  "MKTAYIAKQRQISFVKSHFTRQ-IL"],
                                 "MKTAYIAKQREISFVKSHFSRQDIL")
    rep.append("Example from the handout (P1, P2 profile + Q):\n"
               + fmt_block(["P1", "P2", "Q"], ex) + "\n")
    a1, a4, _ = needleman_wunsch(seq["S1"], seq["S4"])
    prof3, sc3 = align_seq_to_profile([a1, a4], seq["S2"])
    rep.append(f"S2 aligned to profile(S1, S4): score = {sc3:g}\n"
               + fmt_block(["S1", "S4", "S2"], prof3) + "\n")
    prof4, sc4 = align_seq_to_profile(prof3, seq["S6"])
    rep.append(f"S6 aligned to profile(S1, S4, S2): score = {sc4:g}\n"
               + fmt_block(["S1", "S4", "S2", "S6"], prof4) + "\n")

    # (c)
    rep.append("=" * 70 + "\n(c) Pairwise distances + UPGMA guide tree\n" + "=" * 70)
    names = [n for n, _ in cyt]
    Dm = distance_matrix(cyt)
    w = 8
    hdr = " " * 4 + "".join(f"{n:>{w}}" for n in names)
    rows = [f"{n:<4}" + "".join(f"{v:>{w}.4f}" for v in Dm[i]) for i, n in enumerate(names)]
    merges = upgma(names, Dm)
    rep.append("Distance matrix d(x,y) = 1 - identical/alignment length:\n"
               + hdr + "\n" + "\n".join(rows) + "\n\nUPGMA merge order:\n" + fmt_merges(merges))
    heatmap(names, Dm, os.path.join(RESULTS, "q2c_cytochrome_c_distance.png"),
            "Cytochrome c pairwise distance")

    # (d)
    rep.append("\n" + "=" * 70 + "\n(d) Progressive MSA\n" + "=" * 70)
    txt, _, _ = msa_report(os.path.join(DATA, "cytochrome_c_unaligned.fa"), "cytochrome_c")
    rep.append(txt + "\n")
    txt, hnames, hD = msa_report(os.path.join(DATA, "hemoglobins_unaligned.fa"), "hemoglobins")
    rep.append(txt)
    heatmap(hnames, hD, os.path.join(RESULTS, "q2d_hemoglobins_distance.png"),
            "Hemoglobin pairwise distance")

    report = "\n".join(rep)
    print(report)
    with open(os.path.join(RESULTS, "q2_output.txt"), "w") as fh:
        fh.write(report + "\n")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        os.makedirs(RESULTS, exist_ok=True)
        for p in sys.argv[1:]:
            print(msa_report(p, os.path.splitext(os.path.basename(p))[0])[0])
    else:
        main()
