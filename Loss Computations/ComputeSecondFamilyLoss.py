import sympy
import math

n = 8 * 10**9
omega_k = 105
M = 105

def get_nth_prime_excluding_M_factors(n, M):
    excluded_primes = []
    for p in sympy.primefactors(M):
        excluded_primes.append(p)
    
    count = 0
    for p in sympy.primerange(11, 10**7):
        if p not in excluded_primes:
            count += 1
            if count == n:
                return p
    return None

def eta(t, n, M):
    if t <= 1:
        return 0
    if t < 11:
        return 0

    try:
        log_n = math.log(n)
        L = sympy.lcm(M, t**2)

        # Case 1: L <= 1e5  (Lemma 5.3 Case 1)
        if L <= 10**5:
            numerator = 70 * log_n + 21
            denominator = 70 * log_n - 4
            if abs(denominator) < 1e-10:
                return 0
            return (1 / (t * (t - 1))) * (numerator / denominator)

        # Case 2: L > 1e5  (Lemma 5.3 Case 2)
        if n <= L:
            return 49 * log_n / n

        denom_log_arg = n / L
        if denom_log_arg <= 1:
            return 49 * log_n / n

        term1 = 2 / (t * (t - 1))
        term2 = log_n / math.log(denom_log_arg)

        denominator3 = 35 * log_n - 2
        if abs(denominator3) < 1e-10:
            return 0
        term3 = (35 * log_n) / denominator3

        return term1 * term2 * term3

    except (ValueError, ZeroDivisionError, OverflowError):
        return 0

def sum_eta_over_primes(omega_k, n, M):
    # Find starting prime (the prime after the omega_k-th prime)
    last_first_family_prime = get_nth_prime_excluding_M_factors(omega_k, M)
    starting_prime = sympy.nextprime(last_first_family_prime)
    
    total = 0
    prime_count = 0
    prime_contributions = []
    
    # Detailed sum for primes up to 8719 (where Case 2 applies)
    for p in sympy.primerange(starting_prime, 8720):
        if M % p == 0:
            continue
        eta_value = eta(p, n, M)
        total += eta_value
        prime_contributions.append((p, eta_value))
        prime_count += 1
    
    # Tail sum for primes 8731 to 89431 (constant fallback bound)
    tail_prime_count = sum(1 for _ in sympy.primerange(8731, 89432))
    fallback_constant = 49 * math.log(n) / n
    total += fallback_constant * tail_prime_count
    prime_count += tail_prime_count
    
    return total, prime_count, prime_contributions

# Calculate second family
sum_result, num_primes, contributions = sum_eta_over_primes(omega_k, n, M)

# Display results
print(f"\nResults:")
print(f"Number of primes summed: {num_primes}")
print(f"Sum of η(t) over primes: {sum_result:.4f}")
print(f"Sum (more precision): {sum_result:.10f}")
