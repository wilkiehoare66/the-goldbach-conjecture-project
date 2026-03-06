import time
import numpy as np
from numba import njit

MAX_N = 8_000_000_000
START_N = 4_810_000_000  # already verified up to and including this
BLOCK = 10_000_000


def _miller_rabin_witness(a, s, d, n):
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return False
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return False
    return True


def is_prime_u64(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n == p: return True
        if n % p == 0: return False
    d, s = n - 1, 0
    while d % 2 == 0:
        s += 1
        d //= 2
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % n == 0: continue
        if _miller_rabin_witness(a, s, d, n): return False
    return True


def squarefree_table(limit):
    sf = np.ones(limit + 1, dtype=np.uint8)
    sf[0] = 0
    p = 2
    while p * p <= limit:
        sf[p*p::p*p] = 0
        p += 1
    return sf


@njit(cache=True)
def gcd_numba(a, b):
    while b:
        a, b = b, a % b
    return a


@njit(cache=True)
def verify_interval_kernel(n_start, n_end, primes_desc, squarefree, all_n, excluded):
    out = np.empty(n_end - n_start + 1, dtype=np.int64)
    bad = 0
    for n in range(n_start, n_end + 1):
        if all_n == 0 and n % 2 == 0: continue
        g = 0
        success = False
        for idx in range(primes_desc.size):
            p = primes_desc[idx]
            d = n - p
            if d <= 0 or d >= squarefree.size: continue
            if squarefree[d] and excluded[d] == 0:
                g = d if g == 0 else gcd_numba(g, d)
                if g <= 2:
                    success = True
                    break
        if not success:
            out[bad] = n
            bad += 1
    return out[:bad]


def largest_primes_below(hi, count=100):
    found = []
    x = hi - 1 if hi % 2 == 0 else hi
    while len(found) < count:
        if is_prime_u64(x):
            found.append(x)
        x -= 2
    return np.array(found, dtype=np.int64)


def main():
    t0 = time.time()
    sf = squarefree_table(2 * BLOCK)
    excluded = np.zeros(2 * BLOCK + 1, dtype=np.uint8)
    exceptions = []

    # START_N = 481 * BLOCK, so the first block has prev_hi = START_N, cur_start = START_N + 1
    a_start = START_N // BLOCK - 1
    max_a = (MAX_N - 1) // BLOCK - 1
    for a in range(a_start, max_a + 1):
        prev_lo = a * BLOCK + 1
        prev_hi = (a + 1) * BLOCK
        cur_start = prev_hi + 1
        cur_end = min((a + 2) * BLOCK, MAX_N)
        if cur_start > cur_end:
            break
        pool = largest_primes_below(prev_hi, 100)
        bad = verify_interval_kernel(cur_start, cur_end, pool, sf, 0, excluded)
        exceptions.extend(int(x) for x in bad)
        print(f"  done block a={a}, up to {cur_end:,}  [{time.time()-t0:.0f}s]")

    exceptions = sorted(set(exceptions))
    print(f"Checked up to {MAX_N:,} in {time.time()-t0:.1f}s")
    if exceptions:
        print(f"Exceptions ({len(exceptions)}): {exceptions}")
    else:
        print("No exceptions found.")


if __name__ == "__main__":
    main()