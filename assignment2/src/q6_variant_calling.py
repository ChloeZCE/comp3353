from scipy.stats import binom

ERROR = 0.01


def mutant_probability(n_T, n, prior_AT, prior_TT):
    prior = {"AA": 1 - prior_AT - prior_TT, "AT": prior_AT, "TT": prior_TT}
    # chance that one read shows T under each genotype
    p_T = {"AA": ERROR, "AT": 0.5, "TT": 1 - ERROR}
    joint = {g: binom.pmf(n_T, n, p_T[g]) * prior[g] for g in prior}
    total = sum(joint.values())
    for g in prior:
        print(f"  P({g} | data) = {joint[g] / total:.4g}")
    return (joint["AT"] + joint["TT"]) / total


print("(a)")
print("  P(mutant) =", round(mutant_probability(1, 5, 0.02, 0.001), 4))
print("(b)")
print("  P(mutant) =", round(mutant_probability(6, 30, 0.02, 0.001), 4))
print("(c)")
print("  P(mutant) =", round(mutant_probability(6, 30, 0.0005, 0.00001), 4))
