import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from math import gcd

C_ARTIN = 0.3739558136192023
BROADBENT_THETA_CONST = 0.375
BROADBENT_LOWER_X0 = math.exp(20.0)
RR_TABLE2_MAX_X = 1.0e10
N = 8e9
Z = 100_000
MAX_TABLE_MOD = 100_000


def bennett_c0(modulus):
    return 1.0 / 840.0 if modulus <= 10_000 else 1.0 / 160.0


def normalize_header(header):
    return "".join(ch.lower() for ch in header.strip() if ch.isalnum())


def read_tsv_rows(path):
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8", errors="ignore") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        if not sample.strip():
            return []
        dialect = csv.excel_tab
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            pass
        reader = csv.DictReader(fh, dialect=dialect)
        rows = []
        for raw in reader:
            row = {normalize_header(k): (v or "").strip() for k, v in raw.items() if k is not None}
            if any(row.values()):
                rows.append(row)
        return rows


def pick_value(row, aliases, required=True):
    for alias in aliases:
        if alias in row and row[alias] != "":
            return row[alias]
    if required:
        raise KeyError(f"Missing required column. Tried: {', '.join(aliases)}")
    return None


def parse_int(text):
    return int(float(text.replace("_", "")))


def parse_float(text):
    return float(text.replace("_", ""))


@dataclass(frozen=True)
class BennettBound:
    c_theta: float
    x_theta: int


@dataclass(frozen=True)
class RRTable2Bound:
    theta_sqrt_const: float


@dataclass(frozen=True)
class CandidateBound:
    source: str
    modulus: int
    contribution: float
    lower_x0: float


def load_bennett_c_theta(path):
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["m", "modulus", "k", "q"]))
        c_theta = parse_float(pick_value(row, ["ctheta", "c_theta", "rawctheta", "raw_c_theta", "thetabound", "thetaepsilon"]))
        out[modulus] = c_theta
    return out


def load_bennett_x0(path):
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["m", "modulus", "k", "q"]))
        x_theta = parse_int(pick_value(row, ["xtheta", "x_theta", "x0theta", "x0_theta", "xthetamin", "xthetah"]))
        out[modulus] = x_theta
    return out


def load_rr_table1(path):
    threshold_aliases = {
        1.0e10: ["eps1e10", "rr1e10", "1e10", "1010", "10pow10", "10to10"],
        1.0e13: ["eps1e13", "rr1e13", "1e13", "1013", "10pow13", "10to13"],
        1.0e30: ["eps1e30", "rr1e30", "1e30", "1030", "10pow30", "10to30"],
        1.0e100: ["eps1e100", "rr1e100", "1e100", "10100", "10pow100", "10to100"],
    }
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["k", "modulus", "m", "q"]))
        entries = sorted(
            [(x0, parse_float(v)) for x0, aliases in threshold_aliases.items()
             if (v := pick_value(row, aliases, required=False)) is not None],
            key=lambda t: t[0],
        )
        if entries:
            out[modulus] = entries
    return out


def load_rr_table2(path):
    out = {}
    for row in read_tsv_rows(path):
        modulus = parse_int(pick_value(row, ["k", "modulus", "m", "q"]))
        theta = parse_float(pick_value(row, ["theta", "rrtheta", "thetasqrt", "thetaconst"]))
        out[modulus] = RRTable2Bound(theta_sqrt_const=theta)
    return out


def mobius_sieve(n):
    mu = [1] * (n + 1)
    mu[0] = 0
    primes, is_comp = [], [False] * (n + 1)
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            ip = i * p
            if ip > n:
                break
            is_comp[ip] = True
            if i % p == 0:
                mu[ip] = 0
                break
            mu[ip] = -mu[i]
    return mu


def spf_sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n**0.5) + 1):
        if spf[i] == i:
            for j in range(i*i, n+1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def phi_sieve(n):
    phi = list(range(n + 1))
    for p in range(2, n + 1):
        if phi[p] == p:
            for j in range(p, n + 1, p):
                phi[j] -= phi[j] // p
    return phi


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


def collect_required_terms(k, c, mu):
    divs = divisors_from_primes(factor_squarefree(k))
    return [
        (a, d, d * a * a)
        for a in range(1, c + 1) if gcd(a, k) == 1 and mu[a] != 0
        for d in divs
    ]


def rr_table1_candidate(modulus, phi_m, rr_table1):
    entries = rr_table1.get(modulus)
    if not entries:
        return None
    valid = [(x0, eps) for x0, eps in entries if N >= x0]
    if not valid:
        return None
    x0, eps = max(valid, key=lambda t: t[0])
    return CandidateBound(source=f"RR-Table1(x0={x0:.0e})", modulus=modulus, contribution=eps / phi_m, lower_x0=x0)


def rr_table2_candidate(modulus, rr_table2):
    entry = rr_table2.get(modulus)
    if entry is None or N > RR_TABLE2_MAX_X:
        return None
    return CandidateBound(source="RR-Table2", modulus=modulus, contribution=entry.theta_sqrt_const / math.sqrt(N), lower_x0=0.0)


def bennett_candidate(modulus, bennett):
    entry = bennett.get(modulus)
    if entry is None or N < entry.x_theta:
        return None
    return CandidateBound(source="Bennett-raw", modulus=modulus, contribution=entry.c_theta / math.log(N), lower_x0=float(entry.x_theta))


def bennett_c0_candidate(modulus):
    if not (3 <= modulus <= MAX_TABLE_MOD) or N < 8.0e9:
        return None
    return CandidateBound(source="Bennett-c0", modulus=modulus, contribution=bennett_c0(modulus) / math.log(N), lower_x0=8.0e9)


def choose_best_bound(modulus, phi_m, bennett, rr_table1, rr_table2):
    candidates = [c for c in [
        bennett_candidate(modulus, bennett),
        rr_table2_candidate(modulus, rr_table2),
        rr_table1_candidate(modulus, phi_m, rr_table1),
        bennett_c0_candidate(modulus),
    ] if c is not None]
    if not candidates:
        raise ValueError(f"No valid bound for modulus m={modulus} at n={N:.0f}.")
    return min(candidates, key=lambda c: c.contribution)


def compute_first_sum(k, c, mu, phi, bennett, rr_table1, rr_table2):
    total = broadbent_total = non_m1_total = max_x0 = 0.0
    for a, d, modulus in collect_required_terms(k, c, mu):
        if modulus == 1:
            contrib = BROADBENT_THETA_CONST / math.log(N)**3
            total += contrib
            broadbent_total += contrib
            max_x0 = max(max_x0, BROADBENT_LOWER_X0)
        else:
            best = choose_best_bound(modulus, phi[modulus], bennett, rr_table1, rr_table2)
            total += best.contribution
            non_m1_total += best.contribution
            max_x0 = max(max_x0, best.lower_x0)
    return total, broadbent_total, non_m1_total, max_x0


def tail_constant(k, c, mu, spf):
    return sum(
        1.0 / phi_square(a, spf)
        for a in range(c + 1, Z + 1) if gcd(a, k) == 1 and mu[a] != 0
    ) + 4.0 / Z


def err_bound(A, k, alpha, first_sum_total, tail):
    logn = math.log(N)
    return (
        first_sum_total
        + (alpha + 2.0 / (1.0 - 2.0 * A)) * tail
        + (k * N**(-2.0 * A) + math.sqrt(k) * N**(-A)) * logn
    )


def minimize_A(k, alpha, first_sum_total, tail):
    best_i = min(range(500, 4900), key=lambda i: err_bound(i/10000, k, alpha, first_sum_total, tail))
    bestA = best_i / 10000.0
    lo, hi = max(1e-6, bestA - 0.002), min(0.499999, bestA + 0.002)
    bestA = min((lo + t*(hi-lo)/8000 for t in range(8001)), key=lambda A: err_bound(A, k, alpha, first_sum_total, tail))
    return bestA, err_bound(bestA, k, alpha, first_sum_total, tail)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, required=True)
    k = parser.parse_args().k

    c = int(math.floor(math.sqrt(MAX_TABLE_MOD / k)))
    here = Path(__file__).resolve().parent

    bennett_c_theta = load_bennett_c_theta(here / "bennett_c_theta.tsv")
    bennett_x0 = load_bennett_x0(here / "bennett_x0.tsv")
    bennett = {m: BennettBound(c_theta=bennett_c_theta[m], x_theta=bennett_x0[m])
               for m in bennett_c_theta if m in bennett_x0}
    rr_table1 = load_rr_table1(here / "rr_theta_table1.tsv")
    rr_table2 = load_rr_table2(here / "rr_theta_table2.tsv")

    mu = mobius_sieve(Z)
    spf = spf_sieve(Z)
    phi = phi_sieve(MAX_TABLE_MOD)

    primes = factor_squarefree(k)
    alpha = alpha_coeff(primes)
    beta = beta_coeff(primes)

    total, broadbent_total, non_m1_total, x0_req = compute_first_sum(k, c, mu, phi, bennett, rr_table1, rr_table2)
    tail = tail_constant(k, c, mu, spf)
    A_star, err = minimize_A(k, alpha, total, tail)
    lower_bound = beta * C_ARTIN - err

    print(f"k={k}, primes={primes}, c={c}, x0_req={x0_req:.0e}")
    print(f"alpha = {alpha:.12f},  beta = {beta:.12f}")
    print(f"first_sum = {total:.12f}  (broadbent = {broadbent_total:.12f},  non_m1 = {non_m1_total:.12f})")
    print(f"tail  = {tail:.12f}")
    print(f"A*    = {A_star:.7f},  err = {err:.12f}")
    print(f"R_k(n)/n >= {lower_bound:.12f}")


if __name__ == "__main__":
    main()
