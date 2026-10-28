import sys

from bwt import read_fasta, bwt_encode  # reused from Question 4


def build_fm_index(bwt):
    alphabet = sorted(set(bwt))
    # C[c] = number of symbols smaller than c
    C = {}
    total = 0
    for c in alphabet:
        C[c] = total
        total += bwt.count(c)
    # Occ[c][k] = number of c in the first k symbols of the BWT
    Occ = {c: [0] for c in alphabet}
    for ch in bwt:
        for c in alphabet:
            Occ[c].append(Occ[c][-1] + (1 if ch == c else 0))
    return C, Occ


def backward_search(motif, C, Occ, n):
    lo, hi = 0, n
    for c in reversed(motif):
        lo = C[c] + Occ[c][lo]
        hi = C[c] + Occ[c][hi]
        print(f"  {c}: [{lo}, {hi})")
        if lo >= hi:
            return 0
    return hi - lo


path, motif = sys.argv[1], sys.argv[2]
bwt = bwt_encode(read_fasta(path))
n = len(bwt)
C, Occ = build_fm_index(bwt)
alphabet = sorted(C)

print("BWT:", bwt if n <= 100 else bwt[:100] + "...")
print("C:", ", ".join(f"C[{c}]={C[c]}" for c in alphabet))
print("Occ:")
print("   k  " + "  ".join(f"{c:>5}" for c in alphabet))
rows = list(range(n + 1)) if n <= 20 else list(range(6)) + [n]
for k in rows:
    if k == n and n > 20:
        print("  ...")
    print(f"{k:>5} " + "  ".join(f"{Occ[c][k]:>5}" for c in alphabet))
print(f"Backward search for {motif}:")
print("Occurrences:", backward_search(motif, C, Occ, n))
