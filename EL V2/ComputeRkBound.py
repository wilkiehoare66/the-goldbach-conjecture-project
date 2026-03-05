import argparse, math
from pathlib import Path
from math import gcd

C_ARTIN = 0.3739558136192023
MAX_TABLE_MOD = 100000
N, Z = 4.81e9, 100000


def load_c_theta(path):
    c_theta = {}
    for line in path.read_text(errors="ignore").splitlines():
        if line.strip() and not line.startswith("#"):
            parts = line.split()
            c_theta[int(parts[0])] = float(parts[5])
    return c_theta


def load_xcols(path):
    xcols = {}
    for line in path.read_text(errors="ignore").splitlines():
        if line.strip() and not line.startswith("#"):
            parts = line.split()
            xcols[int(parts[0])] = (int(parts[2]), int(parts[3]))
    return xcols


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


def factor_squarefree(k):
    if k % 2 == 0:
        raise ValueError(f"k={k} must be odd.")
    primes, x, p = [], k, 3
    while p * p <= x:
        if x % p == 0:
            primes.append(p)
            x //= p
            if x % p == 0:
                raise ValueError(f"k={k} is not squarefree.")
        p += 2
    if x > 1:
        primes.append(x)
    return sorted(primes)


def divisors_from_primes(primes):
    divs = [1]
    for p in primes:
        divs += [d * p for d in divs]
    return sorted(divs)


def alpha_coeff(primes):
    a = 1.0
    for p in primes:
        a *= (p - 2) / (p - 1)
    return a


def beta_coeff(primes):
    b = 1.0
    for p in primes:
        b *= (p * (p - 2)) / (p*p - p - 1)
    return b


def compute_constants(k, c, c_theta, xcols, mu, spf):
    divs = divisors_from_primes(factor_squarefree(k))
    Csum, x0_req = 0.0, 0
    for a in range(1, c + 1):
        if gcd(a, k) != 1 or mu[a] == 0:
            continue
        for d in divs:
            m = d * a * a
            if m == 1:
                continue
            Csum += c_theta.get(m, 0.0)
            x0_req = max(x0_req, xcols[m][1])
    tail = sum(
        1.0 / phi_square(a, spf)
        for a in range(c + 1, Z + 1) if gcd(a, k) == 1 and mu[a] != 0
    ) + 4.0 / Z
    return Csum, tail, x0_req


def err_bound(n, A, k, alpha, Csum, tail):
    logn = math.log(n)
    return (
        Csum / logn
        + 0.375 / logn**3
        + (alpha + 2.0 / (1.0 - 2.0 * A)) * tail
        + (k * n**(-2.0 * A) + n**-0.5) * logn
    )


def minimize_A(k, alpha, Csum, tail):
    best = min(range(500, 4900), key=lambda i: err_bound(N, i/10000, k, alpha, Csum, tail))
    bestA = best / 10000.0
    lo, hi = max(1e-6, bestA - 0.002), min(0.499999, bestA + 0.002)
    steps = 8000
    bestA = min((lo + t*(hi-lo)/steps for t in range(steps+1)),
                key=lambda A: err_bound(N, A, k, alpha, Csum, tail))
    return bestA, err_bound(N, bestA, k, alpha, Csum, tail)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    k = parser.parse_args().k

    c = int(math.floor(math.sqrt(MAX_TABLE_MOD / k)))
    here = Path(__file__).resolve().parent
    c_theta = load_c_theta(here / "c_all_rounded.txt")
    xcols = load_xcols(here / "x0-all-xm-xe.txt")
    mu = mobius_sieve(Z)
    spf = spf_sieve(Z)

    primes = factor_squarefree(k)
    alpha = alpha_coeff(primes)
    beta = beta_coeff(primes)

    Csum, tail, x0_req = compute_constants(k, c, c_theta, xcols, mu, spf)
    A_star, err = minimize_A(k, alpha, Csum, tail)

    print(f"k={k}, primes={primes}, c={c}, x0_req={x0_req}")
    print(f"alpha = {alpha:.12f},  beta = {beta:.12f}")
    print(f"Csum  = {Csum:.12f},  tail = {tail:.12f}")
    print(f"A*    = {A_star:.7f},  err  = {err:.12f}")
    print(f"R_k(n)/n >= {beta*C_ARTIN - err:.12f}")


if __name__ == "__main__":
    main()