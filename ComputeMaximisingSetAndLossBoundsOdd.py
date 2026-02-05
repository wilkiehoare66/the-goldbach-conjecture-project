import math
import heapq
from typing import List, Tuple, Dict

M = 420  # odd-n branch: include mod 4, so M = lcm(105,4) = 420


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
    Odd-branch eta1 bound:
      - if q|n: fallback
      - q <= 233: explicit AP-theta branch (since 420*q <= 1e5 => q <= 238, last prime 233)
      - q >= 239 and n > 420q: BT-derived branch
      - else (n <= 420q): fallback 98 log n / n
    """
    ln = math.log(n)

    # Patch: non-primitive class if q|n
    if n % q == 0:
        return 98.0 * ln / n

    if q <= 233:
        return (1.0 / (q - 1.0)) * ((35.0 * ln + 21.0) / (35.0 * ln - 4.0))

    # q >= 239
    if n <= 420.0 * q:
        return 98.0 * ln / n

    L = 420.0 * q
    return (2.0 / (q - 1.0)) * (ln / math.log(n / L)) * ((35.0 * ln) / (35.0 * ln - 4.0))


def eta2_bound(t: int, n: int) -> float:
    """
    Odd-branch eta2 bound:
      - if t|n: fallback
      - t <= 13: explicit AP-theta branch (since 420*t^2 <= 1e5 => t^2 <= 238, last prime 13)
      - t >= 17 and n > 420 t^2: BT-derived branch
      - else (n <= 420 t^2): fallback 98 log n / n
    """
    ln = math.log(n)

    # Patch: non-primitive class if t|n
    if n % t == 0:
        return 98.0 * ln / n

    if t <= 13:
        return (1.0 / (t * (t - 1.0))) * ((35.0 * ln + 21.0) / (35.0 * ln - 4.0))

    if n <= 420.0 * (t * t):
        return 98.0 * ln / n

    return (2.0 / (t * (t - 1.0))) * (
        ln / math.log(n / (420.0 * (t * t)))
    ) * ((35.0 * ln) / (35.0 * ln - 4.0))


def candidate_primes(prime_search_limit: int) -> List[int]:
    """
    Candidate primes for divisors of k: primes >= 11 excluding primes dividing 420 (2,3,5,7).
    """
    primes = sieve_primes(prime_search_limit)
    return [p for p in primes if p >= 11 and (M % p) != 0]


def top_m_by_key(primes: List[int], m: int, key_fn) -> Tuple[List[int], float]:
    """Generic top-m selector by key_fn; returns selected primes and sum of key values."""
    top: List[Tuple[float, int]] = []
    for p in primes:
        v = key_fn(p)
        if len(top) < m:
            heapq.heappush(top, (v, p))
        else:
            if v > top[0][0]:
                heapq.heapreplace(top, (v, p))

    top_sorted = sorted(top, reverse=True)
    sel_primes = sorted([p for _, p in top_sorted])
    total_key = sum(v for v, _ in top_sorted)
    return sel_primes, total_key


def top_m_primes_by_eta1(n: int, m: int, prime_search_limit: int = 200_000) -> Tuple[List[int], float]:
    cand = candidate_primes(prime_search_limit)
    return top_m_by_key(cand, m, key_fn=lambda q: eta1_bound(q, n))


def top_m_primes_by_joint(n: int, m: int, prime_search_limit: int = 200_000) -> Tuple[List[int], float]:
    """
    Joint maximiser for:
        sum_{q|k} eta1(q) + sum_{t<=sqrt(n), t prime, t∤M, t∤k} eta2(t)
    For fixed |K|=m, this is equivalent to maximising sum_{q in K} Delta(q), where:
        Delta(q) = eta1(q) - eta2(q)   if q <= sqrt(n)
                 = eta1(q)            if q >  sqrt(n)
    (Because including q in k removes the eta2(q) term from the second-family sum when q<=sqrt(n).)
    """
    cand = candidate_primes(prime_search_limit)
    sqrt_n = int(math.isqrt(n))

    def delta(q: int) -> float:
        e1 = eta1_bound(q, n)
        if q <= sqrt_n:
            return e1 - eta2_bound(q, n)
        return e1

    return top_m_by_key(cand, m, key_fn=delta)


def eta1_breakdown(n: int, k_primes: List[int]) -> Dict[str, float]:
    """Return eta1 totals split by the eta1 cases."""
    out = {"case_q_le_233": 0.0, "case_BT": 0.0, "case_fallback": 0.0, "total": 0.0}
    for q in k_primes:
        v = eta1_bound(q, n)
        out["total"] += v

        if n % q == 0:
            out["case_fallback"] += v
        elif q <= 233:
            out["case_q_le_233"] += v
        elif n <= 420 * q:
            out["case_fallback"] += v
        else:
            out["case_BT"] += v
    return out


def eta2_total_and_breakdown(n: int, k_primes: List[int]) -> Dict[str, float]:
    """
    Compute sum_{t prime, 11<=t<=sqrt(n), t∤M, t∤k} eta2(t),
    with subtotals by eta2 cases.
    """
    k_set = set(k_primes)

    sqrt_n = int(math.isqrt(n))
    primes_t = sieve_primes(sqrt_n)

    cut_bt = int(math.floor(math.sqrt(n / 420.0)))  # largest t with n > 420 t^2 (approx)

    out = {
        "case_t_le_13": 0.0,
        "case_BT": 0.0,
        "case_fallback": 0.0,
        "total": 0.0,
        "t_max": sqrt_n,
        "t_cut_BT": cut_bt,
        "count_terms": 0,
    }

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

        if n % t == 0:
            out["case_fallback"] += v
        elif t <= 13:
            out["case_t_le_13"] += v
        elif n <= 420 * (t * t):
            out["case_fallback"] += v
        else:
            out["case_BT"] += v

    return out


def joint_objective(n: int, k_primes: List[int]) -> float:
    """Compute eta1_sum + eta2_sum (with eta2 summed over t not dividing k)."""
    e1 = sum(eta1_bound(q, n) for q in k_primes)
    e2 = eta2_total_and_breakdown(n, k_primes)["total"]
    return e1 + e2


if __name__ == "__main__":
    n = 8_000_000_000
    m = 41
    prime_search_limit = 200_000

    # eta1-only maximiser
    sel_eta1, eta1_total = top_m_primes_by_eta1(n, m, prime_search_limit=prime_search_limit)
    obj_eta1 = joint_objective(n, sel_eta1)

    # joint maximiser
    sel_joint, delta_total = top_m_primes_by_joint(n, m, prime_search_limit=prime_search_limit)
    obj_joint = joint_objective(n, sel_joint)

    print(f"n={n}, m={m}")
    print(f"prime_search_limit={prime_search_limit}\n")

    print("ETA1-only maximiser:")
    print(f"  min={sel_eta1[0]}, max={sel_eta1[-1]}, size={len(sel_eta1)}")
    print(f"  sum_{{q|k}} eta1(q) = {eta1_total:.12f}")
    print(f"  joint objective (eta1+eta2) = {obj_eta1:.12f}\n")

    print("JOINT maximiser:")
    print(f"  min={sel_joint[0]}, max={sel_joint[-1]}, size={len(sel_joint)}")
    print(f"  sum_{{q in K}} Delta(q) (ranking score) = {delta_total:.12f}")
    print(f"  joint objective (eta1+eta2) = {obj_joint:.12f}\n")

    same = (set(sel_eta1) == set(sel_joint))
    print(f"Do the sets match? {same}")

    if not same:
        removed = sorted(set(sel_eta1) - set(sel_joint))
        added = sorted(set(sel_joint) - set(sel_eta1))
        print("\nDifferences:")
        print(f"  removed ({len(removed)}): {removed}")
        print(f"  added   ({len(added)}): {added}")

    # Breakdowns on the joint set
    e1 = eta1_breakdown(n, sel_joint)
    e2 = eta2_total_and_breakdown(n, sel_joint)

    print("\nBreakdown on JOINT set:")
    print(f"  eta1 total = {e1['total']:.12f}  (q<=233: {e1['case_q_le_233']:.12f}, BT: {e1['case_BT']:.12f}, fallback: {e1['case_fallback']:.12f})")
    print(f"  eta2 total = {e2['total']:.12f}  (t<=13: {e2['case_t_le_13']:.12f}, BT: {e2['case_BT']:.12f}, fallback: {e2['case_fallback']:.12f})")
    print(f"  eta2 terms counted = {e2['count_terms']}")
    print(f"  t_max = {e2['t_max']}   (i.e. floor(sqrt(n)))")
    print(f"  t_cut_BT ≈ {e2['t_cut_BT']}   (largest t with n > 420 t^2)")

"""
Output:

n=8000000000, m=41
prime_search_limit=200000

ETA1-only maximiser:
  min=11, max=359, size=41
  sum_{q|k} eta1(q) = 0.964858855916
  joint objective (eta1+eta2) = 0.988858665184

JOINT maximiser:
  min=11, max=367, size=41
  sum_{q in K} Delta(q) (ranking score) = 0.877956271332
  joint objective (eta1+eta2) = 0.989354722719

Do the sets match? False

Differences:
  removed (1): [89]
  added   (1): [367]

Breakdown on JOINT set:
  eta1 total = 0.964672103493  (q<=233: 0.655711347067, BT: 0.308960756426, fallback: 0.000000000000)
  eta2 total = 0.024682619226  (t<=13: 0.000000000000, BT: 0.022430076580, fallback: 0.002252542646)
  eta2 terms counted = 8615
  t_max = 89442   (i.e. floor(sqrt(n)))
  t_cut_BT ≈ 4364   (largest t with n > 420 t^2)
"""