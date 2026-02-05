import math
import heapq
from typing import List, Tuple, Dict

M = 105  # 3*5*7


def sieve_primes(limit: int) -> List[int]:
    """Simple bytearray sieve; fine up to a few million."""
    if limit < 2:
        return []
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start: limit + 1: step] = b"\x00" * (((limit - start) // step) + 1)
    return [i for i in range(2, limit + 1) if sieve[i]]


def eta1_bound(q: int, n: int) -> float:
    """
    Upper bound from Lemma 5.2:
      - if q|n: fallback
      - q <= 947: explicit AP-theta branch
      - q >= 953 and n > 105q: BT-derived branch
      - else (n <= 105q): fallback 49 log n / n
    """
    ln = math.log(n)

    # Patch: non-primitive class if q|n
    if n % q == 0:
        return 49.0 * ln / n

    if q <= 947:
        return (1.0 / (q - 1.0)) * ((70.0 * ln + 21.0) / (70.0 * ln - 4.0))

    # q >= 953
    if n <= 105.0 * q:
        return 49.0 * ln / n

    L = 105.0 * q
    return (2.0 / (q - 1.0)) * (ln / math.log(n / L)) * ((35.0 * ln) / (35.0 * ln - 2.0))


def eta2_bound(t: int, n: int) -> float:
    """
    Upper bound from Lemma 5.3:
      - if t|n: fallback
      - t <= 29: explicit AP-theta branch
      - t >= 31 and n > 105 t^2: BT-derived branch
      - else (n <= 105 t^2): fallback 49 log n / n
    """
    ln = math.log(n)

    # Patch: non-primitive class if t|n
    if n % t == 0:
        return 49.0 * ln / n

    if t <= 29:
        return (1.0 / (t * (t - 1.0))) * ((70.0 * ln + 21.0) / (70.0 * ln - 4.0))

    if n <= 105.0 * (t * t):
        return 49.0 * ln / n

    return (2.0 / (t * (t - 1.0))) * (
        ln / math.log(n / (105.0 * (t * t)))
    ) * ((35.0 * ln) / (35.0 * ln - 2.0))


def top_m_primes_by_eta1(n: int, m: int, prime_search_limit: int = 200_000) -> Tuple[List[int], float]:
    primes = sieve_primes(prime_search_limit)
    cand = [p for p in primes if p >= 11 and (M % p) != 0]  # exclude 3,5,7

    # keep a min-heap of size m storing (bound, prime)
    top: List[Tuple[float, int]] = []
    for p in cand:
        v = eta1_bound(p, n)
        if len(top) < m:
            heapq.heappush(top, (v, p))
        else:
            if v > top[0][0]:
                heapq.heapreplace(top, (v, p))

    top_sorted = sorted(top, reverse=True)
    sel_primes = sorted([p for _, p in top_sorted])
    total = sum(v for v, _ in top_sorted)
    return sel_primes, total


def eta1_breakdown(n: int, k_primes: List[int]) -> Dict[str, float]:
    """Return eta1 totals split by the Lemma 5.2 cases."""
    ln = math.log(n)
    out = {"case_q_le_947": 0.0, "case_BT": 0.0, "case_fallback": 0.0, "total": 0.0}
    for q in k_primes:
        v = eta1_bound(q, n)
        out["total"] += v

        # classify in the same way as eta1_bound
        if n % q == 0:
            out["case_fallback"] += v
        elif q <= 947:
            out["case_q_le_947"] += v
        elif n <= 105 * q:
            out["case_fallback"] += v
        else:
            out["case_BT"] += v
    return out


def eta2_total_and_breakdown(n: int, k_primes: List[int]) -> Dict[str, float]:
    """
    Compute sum_{t prime, t>=11, t^2<=n, t∤M, t∤k} eta2(t),
    with subtotals by Lemma 5.3 case.
    """
    k_set = set(k_primes)

    sqrt_n = int(math.isqrt(n))
    primes_t = sieve_primes(sqrt_n)

    cut_bt = int(math.floor(math.sqrt(n / 105.0)))  # largest t with n > 105 t^2 (approx)

    out = {"case_t_le_29": 0.0, "case_BT": 0.0, "case_fallback": 0.0, "total": 0.0,
           "t_max": sqrt_n, "t_cut_BT": cut_bt, "count_terms": 0}

    for t in primes_t:
        if t < 11:
            continue
        if M % t == 0:
            continue
        if t in k_set:
            continue  # exclude t|k to avoid double counting with the first family

        v = eta2_bound(t, n)
        out["total"] += v
        out["count_terms"] += 1

        # classify
        if n % t == 0:
            out["case_fallback"] += v
        elif t <= 29:
            out["case_t_le_29"] += v
        elif n <= 105 * (t * t):
            out["case_fallback"] += v
        else:
            out["case_BT"] += v

    return out


if __name__ == "__main__":
    n = 8_000_000_000
    m = 85

    sel, eta1_total = top_m_primes_by_eta1(n, m, prime_search_limit=200_000)

    print(f"n={n}, m={m}")
    print(f"maximising set (for eta1): min={sel[0]}, max={sel[-1]}, size={len(sel)}")

    print("\nFull maximising set:")
    print(sel)

    # eta1 totals (and breakdown)
    e1 = eta1_breakdown(n, sel)
    print("\neta1 totals:")
    print(f"  sum_{'{'}q|k{'}'} eta1(q) = {e1['total']:.12f}")
    print(f"    q<=947   : {e1['case_q_le_947']:.12f}")
    print(f"    BT-case  : {e1['case_BT']:.12f}")
    print(f"    fallback : {e1['case_fallback']:.12f}")

    # eta2 totals (and breakdown)
    e2 = eta2_total_and_breakdown(n, sel)
    print("\neta2 totals (summing over primes t with t^2<=n, t∤M, t∤k):")
    print(f"  t_max = {e2['t_max']}   (i.e. floor(sqrt(n)))")
    print(f"  t_cut_BT ≈ {e2['t_cut_BT']}   (largest t with n > 105 t^2)")
    print(f"  number of t-terms summed = {e2['count_terms']}")
    print(f"  sum eta2(t) = {e2['total']:.12f}")
    print(f"    t<=29     : {e2['case_t_le_29']:.12f}")
    print(f"    BT-case   : {e2['case_BT']:.12f}")
    print(f"    fallback  : {e2['case_fallback']:.12f}")

    print("\nCombined (eta1 + eta2):")
    print(f"  {e1['total'] + e2['total']:.12f}")

"""
Output:

n=8000000000, m=85
maximising set (for eta1): min=11, max=1129, size=85

Full maximising set:
[11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129]

eta1 totals:
  sum_{q|k} eta1(q) = 0.988810144171
    q<=947   : 0.879017658366
    BT-case  : 0.109792485805
    fallback : 0.000000000000

eta2 totals (summing over primes t with t^2<=n, t∤M, t∤k):
  t_max = 89442   (i.e. floor(sqrt(n)))
  t_cut_BT ≈ 8728   (largest t with n > 105 t^2)
  number of t-terms summed = 8571
  sum eta2(t) = 0.007819876971
    t<=29     : 0.000000000000
    BT-case   : 0.006762181940
    fallback  : 0.001057695031

Combined (eta1 + eta2):
  0.996630021142
"""
