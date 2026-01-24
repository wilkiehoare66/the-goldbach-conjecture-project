import sympy
import math

n = 8 * 10**9
omega_k = 105
M = 105
max_prime = 10000

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

def sum_eta_over_primes(omega_k, max_prime, n, M):
    # Find starting prime (the prime after the omega_k-th prime)
    last_first_family_prime = get_nth_prime_excluding_M_factors(omega_k, M)
    starting_prime = sympy.nextprime(last_first_family_prime)
    
    total = 0
    prime_count = 0
    prime_contributions = []
    
    # Generate primes and calculate sum
    for p in sympy.primerange(starting_prime, max_prime + 1):
        if p < 11:
            continue
        if M % p == 0:
            continue

        eta_value = eta(p, n, M)  # Pass M to eta function
        total += eta_value
        prime_contributions.append((p, eta_value))
        prime_count += 1
    
    return total, prime_count, starting_prime, prime_contributions, last_first_family_prime

# Calculate second family
sum_result, num_primes, starting_prime, contributions, last_first = sum_eta_over_primes(omega_k, max_prime, n, M)

# Display results
print(f"\nResults:")
print(f"Number of primes summed: {num_primes}")
print(f"Sum of η(t) over primes: {sum_result:.4f}")
print(f"Sum (more precision): {sum_result:.10f}")

"""
Output:

Results:
Number of primes summed: 1120
Sum of η(t) over primes: 0.0050
Sum (more precision): 0.0050033675
"""
