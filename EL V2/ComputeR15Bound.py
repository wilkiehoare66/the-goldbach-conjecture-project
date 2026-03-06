import math
from pathlib import Path
from math import gcd

JSET = (1, 3, 5, 15)
C_ARTIN = 0.3739558136192023


def load_c_theta(path):
    c_theta = {}
    for line in path.read_text(errors="ignore").splitlines():
        if line.strip() and not line.startswith("#"):
            parts = line.split()
            c_theta[int(parts[0])] = float(parts[5])
    return c_theta


def mobius_sieve(n):
    mu = [1] * (n + 1)
    mu[0] = 0
    primes, is_comp = [], [False] * (n + 1)
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n: break
            is_comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu


def spf_sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n**0.5) + 1):
        if spf[i] == i:
            for j in range(i*i, n+1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def phi_square(a, spf):
    res, x, last = a * a, a, 0
    while x > 1:
        p = spf[x]
        if p != last:
            res -= res // p
            last = p
        while x % p == 0:
            x //= p
    return res


def compute_constants(c, Z, k, c_theta, mu, spf):
    Csum = sum(
        c_theta.get(j * a * a, 0.0)
        for a in range(1, c + 1) if gcd(a, k) == 1 and mu[a] != 0
        for j in JSET if j * a * a > 1
    )
    tail = sum(
        1.0 / phi_square(a, spf)
        for a in range(c + 1, Z + 1) if gcd(a, k) == 1 and mu[a] != 0
    ) + 4.0 / Z
    return Csum, tail


def error_bound(n, A, Csum, tail, k):
    logn = math.log(n)
    return (
        Csum / logn
        + 0.375 / logn**3
        + (0.375 + 2.0 / (1.0 - 2.0 * A)) * tail
        + (k * n**(-2.0 * A) + n**-0.5) * logn
    )


def minimize_A(n, Csum, tail, k):
    best = min(range(500, 4900), key=lambda i: error_bound(n, i/10000, Csum, tail, k))
    bestA = best / 10000.0
    lo, hi = max(1e-6, bestA - 0.001), min(0.499999, bestA + 0.001)
    steps = 4000
    candidates = [lo + t * (hi - lo) / steps for t in range(steps + 1)]
    bestA2 = min(candidates, key=lambda A: error_bound(n, A, Csum, tail, k))
    return bestA2, error_bound(n, bestA2, Csum, tail, k)


N, C, Z, K = 8e9, 81, 100000, 15

def main():
    here = Path(__file__).resolve().parent
    c_theta = load_c_theta(here / "c_all_rounded.txt")
    mu = mobius_sieve(Z)
    spf = spf_sieve(Z)

    Csum, tail = compute_constants(C, Z, K, c_theta, mu, spf)
    A_best, err_best = minimize_A(N, Csum, tail, K)

    print(f"Csum = {Csum:.12f}")
    print(f"tail = {tail:.12f}")
    print(f"A*   = {A_best:.7f}")
    print(f"err  = {err_best:.12f}")
    print(f"R_15(n)/n >= {(9/19)*C_ARTIN - err_best:.12f}")


if __name__ == "__main__":
    main()
