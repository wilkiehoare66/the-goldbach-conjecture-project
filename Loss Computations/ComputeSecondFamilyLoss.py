import sympy
import math

n = 8 * 10**9
omega_k = 26
M = 900
max_prime = 10000

def get_nth_prime_excluding_M_factors(n, M):
    excluded_primes = []
    for p in sympy.primefactors(M):
        excluded_primes.append(p)
    
    count = 0
    for p in sympy.primerange(2, 10**7):
        if p not in excluded_primes:
            count += 1
            if count == n:
                return p
    return None

def eta(t, n):
    if t <= 1:
        return 0
    
    try:
        log_n = math.log(n)
        
        # Calculate each component
        term1 = 2 / (t * (t - 1))
        
        # Check if n/900t² > 1 to avoid log of negative/zero
        denominator_arg = n / (900 * t**2)
        if denominator_arg <= 1:
            return 0
        term2 = log_n / math.log(denominator_arg)
        
        # Check if 7*log(n) - 2 != 0 to avoid division by zero
        denominator3 = 7 * log_n - 2
        if abs(denominator3) < 1e-10:
            return 0
        term3 = (7 * log_n) / denominator3
        
        result = term1 * term2 * term3
        return result
    
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
        eta_value = eta(p, n)
        total += eta_value
        prime_contributions.append((p, eta_value))
        prime_count += 1
    
    return total, prime_count, starting_prime, prime_contributions

# Calculate second family
sum_result, num_primes, starting_prime, contributions = sum_eta_over_primes(omega_k, max_prime, n, M)

# Display results
print(f"\nResults:")
print(f"Number of primes summed: {num_primes}")
print(f"Sum of η(t) over primes: {sum_result:.4f}")
print(f"Sum (more precision): {sum_result:.10f}")

"""
Output:

Results:
Number of primes summed: 1200
Sum of η(t) over primes: 0.0210
Sum (more precision): 0.0209938476
"""
