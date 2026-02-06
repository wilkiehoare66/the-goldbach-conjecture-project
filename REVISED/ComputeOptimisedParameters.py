"""
Parameter Optimisation Script (Generalised)
============================================
Finds optimal (N, C, n_0) parameters for HJ Lemma 3.1 error bounds.

Set OMEGA_K_MAX at the top to optimise for any omega(k) <= OMEGA_K_MAX.
The script automatically determines the worst-case k (product of first
OMEGA_K_MAX odd primes).

Examples:
  OMEGA_K_MAX = 1  ->  k = 3
  OMEGA_K_MAX = 2  ->  k = 3 * 5 = 15
  OMEGA_K_MAX = 3  ->  k = 3 * 5 * 7 = 105
  OMEGA_K_MAX = 4  ->  k = 3 * 5 * 7 * 11 = 1155
"""

import math
import pandas as pd
import numpy as np
import mpmath as mp
from functools import lru_cache
from itertools import combinations
from math import prod

from sympy import factorint, totient
from sympy.ntheory import primerange

# =========================================================
# MAIN CONFIGURATION - SET THIS
# =========================================================
OMEGA_K_MAX = 3  # Maximum number of prime factors

# =========================================================
# OTHER INPUTS
# =========================================================
C_THETA_1_HJ = 9.5913e-4  # HJ: c_theta(1) = 9.5913 * 10^-4
C_ARTIN_LOWER = 0.37395   # HJ: c = 0.37395... (safe lower bound)
C_ALL_FILE = "c_all_rounded.txt"  # Bennett table path

# Grid search parameters (can be adjusted for different OMEGA_K_MAX)
COARSE_N_RANGE = (30000, 100001, 5000)   # (start, stop, step)
COARSE_C_RANGE = (0.27, 0.46, 0.01)      # (start, stop, step)
REFINE_N_DELTA = 8000                     # +- around coarse best
REFINE_N_STEP = 500
REFINE_C_DELTA = 0.02                     # +- around coarse best
REFINE_C_STEP = 0.002

# Search bounds for n
N_LOW_DEFAULT = 8e9
N_HIGH_DEFAULT = 5e18


# -------------------------
# Compute worst-case k
# -------------------------
def compute_worst_case_k(omega_max: int) -> tuple:
    """
    Compute the worst-case k for omega(k) <= omega_max.
    
    This is the product of the first omega_max odd primes.
    Smallest primes give the smallest main term coefficient,
    hence the "worst case" for the margin.
    
    Args:
        omega_max: Maximum number of prime factors
        
    Returns:
        Tuple of (k, list of primes used)
    """
    if omega_max < 1:
        raise ValueError("omega_max must be at least 1")
    
    # First omega_max odd primes: 3, 5, 7, 11, 13, ...
    odd_primes = [p for p in primerange(3, 3 + 6 * omega_max) if p >= 3][:omega_max]
    
    k = prod(odd_primes)
    return k, odd_primes


# -------------------------
# Number theory helpers
# -------------------------
PRIMES = list(primerange(2, 10001))  # Increased for larger k


def is_squarefree(n: int) -> bool:
    """Check if n is squarefree."""
    fac = factorint(n)
    return all(e == 1 for e in fac.values())


@lru_cache(maxsize=None)
def mu2(n: int) -> int:
    """Squarefree indicator function mu^2(n)."""
    return 1 if is_squarefree(n) else 0


@lru_cache(maxsize=None)
def phi(n: int) -> int:
    """Euler's totient function phi(n)."""
    return totient(n) if n > 0 else 0


@lru_cache(maxsize=None)
def phi_sq(n: int) -> int:
    """phi(n^2)."""
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
    """Sum_{b>=1, gcd(b,d)=1} mu^2(b)/phi(b^2) via Euler product."""
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
# Error term E_k(n)/n
# -------------------------
def E_k_over_n(k: int, N: int, C: float, n: float, c_theta_map: dict) -> float:
    """Compute E_k(n)/n."""
    prs_k = primes_of_squarefree(k)
    divs_d = squarefree_divisors_from_primes(prs_k)
    sqrtN = math.sqrt(N)
    logn = math.log(n)

    S1 = 0.0
    tail_total = 0.0
    last_total = 0.0

    for d in divs_d:
        prs_d = primes_of_squarefree(d)
        divs_e = squarefree_divisors_from_primes(prs_d)
        const_coprime_d = total_coprime_constant(d)

        for e in divs_e:
            A0 = int(math.floor(sqrtN * math.sqrt(e / d)))

            for a in range(1, A0 + 1):
                if math.gcd(a, d) == e and mu2(a):
                    m = (d * a * a) // e
                    ct = c_theta_map.get(m)
                    if ct is None:
                        raise KeyError(f"Missing c_theta({m}) - increase table or reduce N")
                    S1 += ct

            b0 = A0 // e
            partial = 0.0
            for b in range(1, b0 + 1):
                if math.gcd(b, d) == 1 and mu2(b):
                    partial += 1.0 / phi_sq(b)
            tail_b = max(0.0, const_coprime_d - partial)
            tail_a = tail_b / phi(e * e)
            tail_total += (1.0 / phi(d // e)) * tail_a

            last_total += (
                (n ** (-0.5)) * 0.5 * (1.0 / e - 1.0 / d)
                + (1.0 / math.sqrt(d * e)) * ((n ** (-C)) + (n ** (-2 * C)))
            )

    term1 = S1 / logn
    term2 = ((1 + 2 * C) / (1 - 2 * C)) * tail_total
    term3 = logn * last_total

    return term1 + term2 + term3


def margin(k: int, N: int, C: float, n: float, c_theta_map: dict) -> float:
    """Main - Error, normalised by n."""
    main = main_coeff_over_n(k)
    E_over_n = E_k_over_n(k, N, C, n, c_theta_map)
    return main - E_over_n


# -------------------------
# Threshold finding
# -------------------------
def threshold_n(
    k: int,
    N: int,
    C: float,
    c_theta_map: dict,
    n_low: float = N_LOW_DEFAULT,
    n_high: float = N_HIGH_DEFAULT,
    iters: int = 50,
) -> float | None:
    """Find threshold n_0 via log-bisection."""
    sqrtN = math.sqrt(N)
    if n_low ** C <= sqrtN:
        n_low = (sqrtN + 1.0) ** (1.0 / C)

    try:
        if margin(k, N, C, n_low, c_theta_map) > 0:
            return n_low
        if margin(k, N, C, n_high, c_theta_map) <= 0:
            return None
    except KeyError:
        return None  # c_theta table insufficient for this N

    lo, hi = n_low, n_high
    for _ in range(iters):
        mid = math.sqrt(lo * hi)
        try:
            if margin(k, N, C, mid, c_theta_map) > 0:
                hi = mid
            else:
                lo = mid
        except KeyError:
            return None
    return hi


# -------------------------
# Grid search
# -------------------------
def coarse_optimise(k: int, c_theta_map: dict, verbose: bool = True) -> tuple:
    """Coarse grid search for optimal (N, C)."""
    best = None
    best_params = None

    Ns = list(range(*[int(x) for x in COARSE_N_RANGE]))
    Cs = [round(x, 2) for x in np.arange(*COARSE_C_RANGE)]

    total = len(Ns) * len(Cs)
    count = 0

    for N in Ns:
        for C in Cs:
            count += 1
            if verbose and count % 50 == 0:
                print(f"  Coarse: {count}/{total} ({100*count/total:.0f}%)")
            t = threshold_n(k, N, C, c_theta_map, iters=35)
            if t is None:
                continue
            if best is None or t < best:
                best = t
                best_params = (N, C)

    return best, best_params


def refine_around(k: int, c_theta_map: dict, N0: int, C0: float, verbose: bool = True) -> tuple:
    """Refine search around coarse optimum."""
    best = None
    best_params = None

    N_lo = max(20000, N0 - REFINE_N_DELTA)
    N_hi = N0 + REFINE_N_DELTA + 1
    Ns = list(range(N_lo, N_hi, REFINE_N_STEP))

    C_lo = max(0.26, C0 - REFINE_C_DELTA)
    C_hi = min(0.49, C0 + REFINE_C_DELTA) + 1e-9
    Cs = [round(x, 4) for x in np.arange(C_lo, C_hi, REFINE_C_STEP)]

    total = len(Ns) * len(Cs)
    count = 0

    for N in Ns:
        for C in Cs:
            count += 1
            if verbose and count % 100 == 0:
                print(f"  Refine: {count}/{total} ({100*count/total:.0f}%)")
            t = threshold_n(k, N, C, c_theta_map, iters=45)
            if t is None:
                continue
            if best is None or t < best:
                best = t
                best_params = (N, C)

    return best, best_params


# -------------------------
# MAIN
# -------------------------
if __name__ == "__main__":
    # Compute worst-case k
    K_WORST, primes_used = compute_worst_case_k(OMEGA_K_MAX)

    print(f"  OMEGA_K_MAX = {OMEGA_K_MAX}")
    print(f"  Worst-case k = {K_WORST} = {' * '.join(map(str, primes_used))}")
    print(f"  Number of (d,e) pairs: {2**OMEGA_K_MAX * (OMEGA_K_MAX + 1)}")
    print(f"\n  Main term coefficient A_k(n)/n >= {main_coeff_over_n(K_WORST):.6e}")
    
    print(f"\nLoading c_theta table from: {C_ALL_FILE}")
    c_theta_map = load_c_theta_map(C_ALL_FILE)
    print(f"Loaded {len(c_theta_map)} c_theta values.")

    print("\nSTEP 1: Coarse Grid Search")
    print(f"  N in [{COARSE_N_RANGE[0]}, {COARSE_N_RANGE[1]}) step {COARSE_N_RANGE[2]}")
    print(f"  C in [{COARSE_C_RANGE[0]}, {COARSE_C_RANGE[1]}) step {COARSE_C_RANGE[2]}")
    
    coarse_result = coarse_optimise(K_WORST, c_theta_map)
    
    if coarse_result[0] is None:
        print("\n  No valid parameters found in coarse search!")
        print("    Try increasing N_HIGH_DEFAULT or adjusting grid ranges.")
        exit(1)
    
    best, (N_best, C_best) = coarse_result
    print(f"\n  COARSE RESULT:")
    print(f"    n_0 ≈ {best:.6e}")
    print(f"    N   = {N_best}")
    print(f"    C   = {C_best}")

    # Refine
    print("\nSTEP 2: Refined Grid Search")
    print(f"  N in [{max(20000, N_best - REFINE_N_DELTA)}, {N_best + REFINE_N_DELTA}] step {REFINE_N_STEP}")
    print(f"  C in [{max(0.26, C_best - REFINE_C_DELTA):.2f}, {min(0.49, C_best + REFINE_C_DELTA):.2f}] step {REFINE_C_STEP}")
    
    best2, (N_best2, C_best2) = refine_around(K_WORST, c_theta_map, N_best, C_best)
    
    print(f"\n  REFINED RESULT:")
    print(f"    n_0 ≈ {best2:.6e}")
    print(f"    N   = {N_best2}")
    print(f"    C   = {C_best2}")

"""
Output:

  OMEGA_K_MAX = 3
  Worst-case k = 105 = 3 * 5 * 7
  Number of (d,e) pairs: 32

  Main term coefficient A_k(n)/n >= 1.512121e-01

Loading c_theta table from: c_all_rounded.txt
Loaded 99999 c_theta values.

STEP 1: Coarse Grid Search
  N in [30000, 100001) step 5000
  C in [0.27, 0.46) step 0.01
  Coarse: 50/285 (18%)
  Coarse: 100/285 (35%)
  Coarse: 150/285 (53%)
  Coarse: 200/285 (70%)
  Coarse: 250/285 (88%)

  COARSE RESULT:
    n_0 ≈ 1.767447e+12
    N   = 55000
    C   = 0.31

STEP 2: Refined Grid Search
  N in [47000, 63000] step 500
  C in [0.29, 0.33] step 0.002
  Refine: 100/693 (14%)
  Refine: 200/693 (29%)
  Refine: 300/693 (43%)
  Refine: 400/693 (58%)
  Refine: 500/693 (72%)
  Refine: 600/693 (87%)

  REFINED RESULT:
    n_0 ≈ 1.732573e+12
    N   = 54000
    C   = 0.312
"""
