"""
Computation & Verification Script (Generalised)
================================================
Given (N, C, n_0) parameters, computes the full breakdown of:
  - Main term coefficient A_k(n)/n
  - Error term E_k(n)/n and its components
  - Margin = Main - Error

Set OMEGA_K_MAX and the optimised parameters from OptimisationComputation.py.
"""

import math
import pandas as pd
import mpmath as mp
from functools import lru_cache
from itertools import combinations
from math import prod

from sympy import factorint, totient
from sympy.ntheory import primerange

# =========================================================
# CONFIGURATION - SET THESE
# =========================================================
OMEGA_K_MAX = 3  # Must match parameter_optimisation.py

# From parameter_optimisation.py results:
N_USE = 60000
C_USE = 0.31
N0_USE = 2e12

# =========================================================
# OTHER INPUTS
# =========================================================
C_THETA_1_HJ = 9.5913e-4  # HJ: c_theta(1) = 9.5913 × 10^-4
C_ARTIN_LOWER = 0.37395   # HJ: c = 0.37395... (safe lower bound)
C_ALL_FILE = "c_all_rounded.txt"  # Bennett table path


# -------------------------
# Compute worst-case k
# -------------------------
def compute_worst_case_k(omega_max: int) -> tuple:
    """
    Compute the worst-case k for ω(k) ≤ omega_max.
    Product of first omega_max odd primes.
    """
    if omega_max < 1:
        raise ValueError("omega_max must be at least 1")
    
    odd_primes = [p for p in primerange(3, 3 + 6 * omega_max) if p >= 3][:omega_max]
    k = prod(odd_primes)
    return k, odd_primes


# -------------------------
# Number theory helpers
# -------------------------
PRIMES = list(primerange(2, 10001))


def is_squarefree(n: int) -> bool:
    """Check if n is squarefree."""
    fac = factorint(n)
    return all(e == 1 for e in fac.values())


@lru_cache(maxsize=None)
def mu2(n: int) -> int:
    """Squarefree indicator function μ²(n)."""
    return 1 if is_squarefree(n) else 0


@lru_cache(maxsize=None)
def phi(n: int) -> int:
    """Euler's totient function φ(n)."""
    return totient(n) if n > 0 else 0


@lru_cache(maxsize=None)
def phi_sq(n: int) -> int:
    """φ(n²)."""
    return totient(n * n)


def primes_of_squarefree(k: int) -> list:
    """Return sorted list of prime factors of squarefree k."""
    return sorted(factorint(k).keys())


def squarefree_divisors_from_primes(prs: list) -> list:
    """Return all squarefree divisors from a list of primes."""
    result = [1]
    for r in range(1, len(prs) + 1):
        result.extend(prod(c) for c in combinations(prs, r))
    return sorted(result)


# -------------------------
# Load c_theta table
# -------------------------
def load_c_theta_map(path: str) -> dict:
    """Load Bennett's c_theta values from TSV file."""
    df = pd.read_csv(path, sep="\t")
    ctheta = {int(r["#q"]): float(r["c_theta"]) for _, r in df.iterrows()}
    ctheta[1] = C_THETA_1_HJ
    return ctheta


# -------------------------
# Tail constant
# -------------------------
mp.mp.dps = 80
S_TOTAL = float(mp.zeta(2) * mp.zeta(3) / mp.zeta(6))


def total_coprime_constant(d: int) -> float:
    """Sum_{b≥1, gcd(b,d)=1} μ²(b)/φ(b²) via Euler product."""
    denom = 1.0
    for p in primes_of_squarefree(d):
        denom *= (1.0 + 1.0 / (p * (p - 1)))
    return S_TOTAL / denom


# -------------------------
# Main term coefficient
# -------------------------
def main_coeff_over_n(k: int, c_artin_lower: float = C_ARTIN_LOWER) -> float:
    """Lower bound for A_k(n)/n."""
    prod_val = 1.0
    for q in primes_of_squarefree(k):
        prod_val *= (1.0 - (q - 1.0) / (q * q - q - 1.0))
    return c_artin_lower * prod_val


# -------------------------
# Error term with detailed breakdown
# -------------------------
def E_k_over_n_detailed(k: int, N: int, C: float, n: float, c_theta_map: dict) -> tuple:
    """
    Compute E_k(n)/n with detailed breakdown.
    
    Returns:
        tuple: (E_over_n, debug_dict)
    """
    prs_k = primes_of_squarefree(k)
    divs_d = squarefree_divisors_from_primes(prs_k)
    sqrtN = math.sqrt(N)
    logn = math.log(n)

    S1 = 0.0
    tail_total = 0.0
    last_total = 0.0
    
    de_contributions = []

    for d in divs_d:
        prs_d = primes_of_squarefree(d)
        divs_e = squarefree_divisors_from_primes(prs_d)
        const_coprime_d = total_coprime_constant(d)

        for e in divs_e:
            A0 = int(math.floor(sqrtN * math.sqrt(e / d)))
            
            s1_de = 0.0
            tail_de = 0.0
            last_de = 0.0

            for a in range(1, A0 + 1):
                if math.gcd(a, d) == e and mu2(a):
                    m = (d * a * a) // e
                    ct = c_theta_map.get(m)
                    if ct is None:
                        raise KeyError(f"Missing c_theta({m})")
                    s1_de += ct
            S1 += s1_de

            b0 = A0 // e
            partial = 0.0
            for b in range(1, b0 + 1):
                if math.gcd(b, d) == 1 and mu2(b):
                    partial += 1.0 / phi_sq(b)
            tail_b = max(0.0, const_coprime_d - partial)
            tail_a = tail_b / phi(e * e)
            tail_de = (1.0 / phi(d // e)) * tail_a
            tail_total += tail_de

            last_de = (
                (n ** (-0.5)) * (1.0 / e - 1.0 / d)
                + (1.0 / math.sqrt(d * e)) * (n ** (-C))
                + (n ** (-2 * C))
            )
            last_total += last_de
            
            de_contributions.append({
                "d": d, "e": e, "A0": A0,
                "s1_de": s1_de, "tail_de": tail_de, "last_de": last_de
            })

    term1 = S1 / logn
    term2 = ((1 + 2 * C) / (1 - 2 * C)) * tail_total
    term3 = logn * last_total
    E_over_n = term1 + term2 + term3

    dbg = {
        "logn": logn,
        "sqrtN": sqrtN,
        "S1_raw": S1,
        "tail_raw": tail_total,
        "last_raw": last_total,
        "term1_cTheta": term1,
        "term2_tail": term2,
        "term3_small": term3,
        "E_over_n": E_over_n,
        "de_contributions": de_contributions,
    }
    return E_over_n, dbg


# -------------------------
# Lemma conditions check
# -------------------------
def check_lemma_conditions(N: int, C: float, n: float) -> dict:
    """Check all conditions for Lemma 3.1."""
    sqrtN = math.sqrt(N)
    nC = n ** C
    
    checks = {
        "n > 0": n > 0,
        "N > 0": N > 0,
        "0 < C < 0.5": 0 < C < 0.5,
        "n^C > sqrt(N)": nC > sqrtN,
    }
    
    return checks


# -------------------------
# MAIN
# -------------------------
if __name__ == "__main__":
    # Compute worst-case k
    K_WORST, primes_used = compute_worst_case_k(OMEGA_K_MAX)

    # Load data
    print(f"\nLoading c_theta table from: {C_ALL_FILE}")
    c_theta_map = load_c_theta_map(C_ALL_FILE)
    print(f"Loaded {len(c_theta_map)} c_theta values.")
    
    # Parameters
    print("\nINPUT PARAMETERS")
    print(f"  OMEGA_K_MAX = {OMEGA_K_MAX}")
    print(f"  k   = {K_WORST} = {' * '.join(map(str, primes_used))}")
    print(f"  N   = {N_USE}")
    print(f"  C   = {C_USE}")
    print(f"  n_0 = {N0_USE:.6e}")
    
    # Check lemma conditions
    print("\nLEMMA CONDITIONS CHECK")
    checks = check_lemma_conditions(N_USE, C_USE, N0_USE)
    all_ok = True
    for condition, passed in checks.items():
        status = "PASS" if passed else "FAIL"
        print(f"  {condition}: {status}")
        if not passed:
            all_ok = False
    
    if not all_ok:
        print("\n  *** WARNING: Not all conditions satisfied! ***")
    
    # Main term
    print("\nMAIN TERM")
    main = main_coeff_over_n(K_WORST)
    print(f"  A_k(n)/n lower bound = {main:.10f}")
    print(f"  (using c_Artin >= {C_ARTIN_LOWER})")
    
    # Error term breakdown
    print("\nERROR TERM BREAKDOWN")
    E_over_n, dbg = E_k_over_n_detailed(K_WORST, N_USE, C_USE, N0_USE, c_theta_map)
    
    print(f"  log(n_0) = {dbg['logn']:.6f}")
    print(f"  sqrt(N)  = {dbg['sqrtN']:.6f}")
    print(f"  Term 1 (c_theta sum / log n):")
    print(f"    S1_raw      = {dbg['S1_raw']:.10f}")
    print(f"    term1       = {dbg['term1_cTheta']:.10f}")
    print(f"  Term 2 (tail sum * (1+2C)/(1-2C)):")
    print(f"    tail_raw    = {dbg['tail_raw']:.10f}")
    print(f"    multiplier  = {(1 + 2*C_USE)/(1 - 2*C_USE):.6f}")
    print(f"    term2       = {dbg['term2_tail']:.10f}")
    print(f"  Term 3 (small terms * log n):")
    print(f"    last_raw    = {dbg['last_raw']:.10e}")
    print(f"    term3       = {dbg['term3_small']:.10e}")
    print(f"  TOTAL E_k(n)/n = {dbg['E_over_n']:.10f}")
    
    print("\nMARGIN (MAIN - ERROR)")
    margin_val = main - E_over_n
    print(f"  Main term  = {main:.10f}")
    print(f"  Error term = {E_over_n:.10f}")
    print(f"  Margin     = {margin_val:.10f}")
    
    if margin_val > 0:
        print(f"\nPOSITIVE MARGIN: Lemma 3.1 applies for n >= {N0_USE:.2e}")
    else:
        print(f"\nNEGATIVE MARGIN: Need larger n_0 or different (N, C)")

"""
Output:

Loading c_theta table from: c_all_rounded.txt
Loaded 99999 c_theta values.

INPUT PARAMETERS
  OMEGA_K_MAX = 3
  k   = 105 = 3 * 5 * 7
  N   = 60000
  C   = 0.31
  n_0 = 2.000000e+12

LEMMA CONDITIONS CHECK
  n > 0: PASS
  N > 0: PASS
  0 < C < 0.5: PASS
  n^C > sqrt(N): PASS

MAIN TERM
  A_k(n)/n lower bound = 0.1512121309
  (using c_Artin >= 0.37395)

ERROR TERM BREAKDOWN
  log(n_0) = 28.324168
  sqrt(N)  = 244.948974
  Term 1 (c_theta sum / log n):
    S1_raw      = 1.3148746300
    term1       = 0.0464223562
  Term 2 (tail sum * (1+2C)/(1-2C)):
    tail_raw    = 0.0195671747
    multiplier  = 4.263158
    term2       = 0.0834179552
  Term 3 (small terms * log n):
    last_raw    = 7.4198949065e-04
    term3       = 2.1016235207e-02
  TOTAL E_k(n)/n = 0.1508565466

MARGIN (MAIN - ERROR)
  Main term  = 0.1512121309
  Error term = 0.1508565466
  Margin     = 0.0003555843

POSITIVE MARGIN: Lemma 3.1 applies for n >= 2.00e+12
"""
