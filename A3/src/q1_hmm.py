"""COMP3353 A3 Q1: two-state transmembrane HMM.

Implements (a) Viterbi decoding, (b) Forward, (c) Backward and
(d) Posterior decoding.  All computations are done in log space so the
probabilities do not underflow for long sequences.

Usage:  python q1_hmm.py data/Q1/P02724.fasta [data/Q1/P02945.fasta ...]
"""
import os
import sys

import numpy as np
from scipy.special import logsumexp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STATES = ["M", "S"]                      # index 0 = M (membrane), 1 = S (soluble)
HYDROPHOBIC = set("GAVILMFW")             # group H
POLAR = set("RNDCEQHKPSTY")               # group P

# Model parameters
INIT = np.array([0.1, 0.9])              # P(M), P(S)
TRANS = np.array([[0.95, 0.05],          # a[k' -> k], rows = from, cols = to
                  [0.05, 0.95]])
# Group-level emission probabilities; each amino acid in a group gets an
# equal share (8 hydrophobic, 12 polar amino acids).
GROUP_EMIT = {"M": {"H": 0.85, "P": 0.15}, "S": {"H": 0.30, "P": 0.70}}


def emission_prob(state, aa):
    if aa in HYDROPHOBIC:
        return GROUP_EMIT[state]["H"] / len(HYDROPHOBIC)
    if aa in POLAR:
        return GROUP_EMIT[state]["P"] / len(POLAR)
    raise ValueError(f"unknown amino acid {aa!r}")


def read_single_fasta(path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
            else:
                seq.append(line)
    return header, "".join(seq).upper()


def log_emissions(seq):
    """L x K matrix of log e_k(x_t)."""
    return np.log(np.array([[emission_prob(s, aa) for s in STATES] for aa in seq]))


LOG_INIT = np.log(INIT)
LOG_TRANS = np.log(TRANS)


# ---------------------------------------------------------------- (a) Viterbi
def viterbi(seq):
    E = log_emissions(seq)
    L, K = E.shape
    V = np.zeros((L, K))
    ptr = np.zeros((L, K), dtype=int)
    V[0] = E[0] + LOG_INIT                          # V_k(1) = e_k(x1) P(k)
    for t in range(1, L):
        # cand[k', k] = V_k'(t-1) + log a_{k'->k}
        cand = V[t - 1][:, None] + LOG_TRANS
        ptr[t] = cand.argmax(axis=0)
        V[t] = E[t] + cand.max(axis=0)
    # traceback
    path = np.zeros(L, dtype=int)
    path[-1] = V[-1].argmax()
    for t in range(L - 1, 0, -1):
        path[t - 1] = ptr[t, path[t]]
    return "".join(STATES[k] for k in path), V[-1].max()


# ---------------------------------------------------------------- (b) Forward
def forward(seq):
    E = log_emissions(seq)
    L, K = E.shape
    F = np.zeros((L, K))
    F[0] = E[0] + LOG_INIT
    for t in range(1, L):
        F[t] = E[t] + logsumexp(F[t - 1][:, None] + LOG_TRANS, axis=0)
    return F, logsumexp(F[-1])                      # log P(x) = log sum_k F_k(L)


# --------------------------------------------------------------- (c) Backward
def backward(seq):
    E = log_emissions(seq)
    L, K = E.shape
    B = np.zeros((L, K))                            # B_k(L) = 1  -> log 0
    for t in range(L - 2, -1, -1):
        # B_k(t) = sum_k' a_{k->k'} e_k'(x_{t+1}) B_k'(t+1)
        B[t] = logsumexp(LOG_TRANS + (E[t + 1] + B[t + 1])[None, :], axis=1)
    logp = logsumexp(LOG_INIT + E[0] + B[0])        # sum_k P(k) e_k(x1) B_k(1)
    return B, logp


# ------------------------------------------------------ (d) Posterior decoding
def posterior(seq):
    F, logp_f = forward(seq)
    B, _ = backward(seq)
    log_post = F + B - logsumexp(F + B, axis=1, keepdims=True)
    post = np.exp(log_post)
    path = "".join(STATES[k] for k in post.argmax(axis=1))
    return path, post


# ------------------------------------------------------------------ helpers
def segments(path, state="M"):
    """1-based inclusive (start, end) runs of `state`."""
    out, start = [], None
    for i, s in enumerate(path, 1):
        if s == state and start is None:
            start = i
        if s != state and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(path)))
    return out


def sci(logp):
    """Format exp(logp) in scientific notation without underflow."""
    l10 = logp / np.log(10)
    e = int(np.floor(l10))
    return f"{10 ** (l10 - e):.6f}e{e}"


def wrap(s, n=60):
    return "\n".join(s[i:i + n] for i in range(0, len(s), n))


def plot(name, seq, vit, post_path, post, out_png):
    x = np.arange(1, len(seq) + 1)
    fig, axes = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
    ax = axes[0]
    ax.step(x, [1 if s == "M" else 0 for s in vit], where="mid",
            label="Viterbi", color="C0", lw=2)
    ax.step(x, [1 if s == "M" else 0 for s in post_path], where="mid",
            label="Posterior decoding", color="C3", lw=1.2, ls="--")
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["S", "M"])
    ax.set_ylim(-0.2, 1.3)
    ax.set_ylabel("hidden state")
    ax.legend(loc="upper right", ncol=2)
    ax.set_title(f"{name}: decoded hidden states (M = transmembrane, S = soluble)")
    ax = axes[1]
    ax.plot(x, post[:, 0], color="C2")
    ax.axhline(0.5, color="grey", lw=0.8, ls=":")
    hyd = [1.0 if aa in HYDROPHOBIC else 0.0 for aa in seq]
    ax.bar(x, np.array(hyd) * 0.08, bottom=-0.1, width=1, color="k", alpha=0.4,
           label="hydrophobic residue")
    ax.set_ylim(-0.12, 1.05)
    ax.set_ylabel("P(M | x)")
    ax.set_xlabel("residue position")
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def run(path, outdir):
    header, seq = read_single_fasta(path)
    name = os.path.splitext(os.path.basename(path))[0]
    vit, vit_logp = viterbi(seq)
    F, logp_f = forward(seq)
    B, logp_b = backward(seq)
    post_path, post = posterior(seq)
    diff = [i + 1 for i, (a, b) in enumerate(zip(vit, post_path)) if a != b]

    lines = [
        f"=== {name} ({header.split()[0]}) length = {len(seq)} ===",
        "Sequence:", wrap(seq),
        "",
        "(a) Viterbi path:", wrap(vit),
        f"    ln P(x, best path) = {vit_logp:.6f}",
        f"    M segments (1-based): {segments(vit)}",
        "",
        f"(b) Forward : ln P(x) = {logp_f:.10f}   P(x) = {sci(logp_f)}",
        f"(c) Backward: ln P(x) = {logp_b:.10f}   P(x) = {sci(logp_b)}",
        f"    |difference| = {abs(logp_f - logp_b):.3e}",
        "",
        "(d) Posterior-decoded path:", wrap(post_path),
        f"    M segments (1-based): {segments(post_path)}",
        f"    positions differing from Viterbi: {len(diff)} -> {diff}",
        "",
    ]
    png = os.path.join(outdir, f"q1_{name}_states.png")
    plot(name, seq, vit, post_path, post, png)
    lines.append(f"Plot: results/{os.path.basename(png)}")
    return "\n".join(lines)


if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(here, "..", "results")
    os.makedirs(outdir, exist_ok=True)
    files = sys.argv[1:] or [os.path.join(here, "..", "data", "Q1", f)
                             for f in ("P02724.fasta", "P02945.fasta")]
    report = "\n\n".join(run(f, outdir) for f in files)
    print(report)
    with open(os.path.join(outdir, "q1_output.txt"), "w") as fh:
        fh.write(report + "\n")
