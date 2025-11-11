import sympy
import math

n = 8 * 10**9
starting_prime = 17
max_prime = 10000

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

def sum_eta_over_primes(starting_prime, max_prime, n):
    total = 0
    prime_count = 0
    
    # Generate primes and calculate sum
    for p in sympy.primerange(starting_prime, max_prime + 1):
        eta_value = eta(p, n)
        total += eta_value
        prime_count += 1
    
    return total, prime_count

# Run the calculation
sum_result, num_primes = sum_eta_over_primes(starting_prime, max_prime, n)

# Display results
print(f"\nResults:")
print(f"Number of primes summed: {num_primes}")
print(f"Sum of η(t) over primes: {sum_result:.4f}")
print(f"Sum (more precision): {sum_result:.10f}")