"""Q6: Bayesian genotype posteriors for a site with reference base A.

P(G | k T-reads of n) = P(k | n, p_G) P(G) / sum_G' P(k | n, p_G') P(G')
with p_G = probability that one read shows T given genotype G:
    AA: e            (an A misread as T)
    AT: 0.5 * (1 - e) + 0.5 * e = 0.5
    TT: 1 - e
"""
from scipy.stats import binom

ERR = 0.01


def posteriors(n_t, n, prior_het, prior_hom, err=ERR):
    p_t = {"AA": err, "AT": 0.5 * (1 - err) + 0.5 * err, "TT": 1 - err}
    prior = {"AA": 1 - prior_het - prior_hom, "AT": prior_het, "TT": prior_hom}
    joint = {g: binom.pmf(n_t, n, p_t[g]) * prior[g] for g in p_t}
    evidence = sum(joint.values())
    print(f"{n - n_t} A / {n_t} T out of {n} reads, "
          f"priors AA={prior['AA']:.5f} AT={prior_het} TT={prior_hom}")
    print(f"  {'G':<3} {'P(T|G)':>7} {'likelihood':>12} {'prior':>8} "
          f"{'lik x prior':>12} {'posterior':>10}")
    for g in p_t:
        lik = binom.pmf(n_t, n, p_t[g])
        print(f"  {g:<3} {p_t[g]:>7.3f} {lik:>12.4e} {prior[g]:>8.5f} "
              f"{joint[g]:>12.4e} {joint[g] / evidence:>10.6f}")
    mutant = (joint["AT"] + joint["TT"]) / evidence
    print(f"  evidence P(D) = {evidence:.4e}")
    print(f"  P(mutant | D) = P(AT|D) + P(TT|D) = {mutant:.6f}\n")
    return mutant


if __name__ == "__main__":
    print("(a)"); posteriors(1, 5, 0.02, 0.001)
    print("(b)"); posteriors(6, 30, 0.02, 0.001)
    print("(c)"); posteriors(6, 30, 0.0005, 0.00001)
