import sympy
import math

n = 8 * 10**9
omega_k = 26
M = 900

def eta(q, n, M):
    if q <= 1:
        return 0
    
    L = sympy.lcm(M, q)
    
    try:
        log_n = math.log(n)
        
        # Case 1: L ≤ 10^5 (from Lemma 5.2)
        if L <= 10**5:
            if q <= 109:
                # First subcase: q ≤ 109
                numerator = 14 * log_n + 21
                denominator = 14 * log_n - 14
                
                # Check for division by zero
                if abs(denominator) < 1e-10:
                    return 0
                    
                eta_value = (1 / (q - 1)) * (numerator / denominator)
                return eta_value
            
            else:  # q ≥ 113 and n > 900q
                # Check the condition n > 900q
                if n <= 900 * q:
                    return 0
                
                # Check if n/L > 1 to avoid log of negative/zero
                if n/L <= 1:
                    return 0
                
                # Calculate each component
                term1 = 2 / (q - 1)
                term2 = log_n / math.log(n / L)
                
                # Check if 7*log(n) - 2 != 0 to avoid division by zero
                denominator3 = 7 * log_n - 2
                if abs(denominator3) < 1e-10:
                    return 0
                term3 = (7 * log_n) / denominator3
                
                eta_value = term1 * term2 * term3
                return eta_value
        
        # Case 2: L > 10^5 (from Lemma 5.2)
        else:
            # Check the condition n > 900q^2
            if n <= 900 * q**2:
                return 0
            
            # Check if n/L > 1 to avoid log of negative/zero
            if n/L <= 1:
                return 0
            
            # Calculate each component
            term1 = 2 / (q - 1)
            term2 = log_n / math.log(n / L)
            
            # Check if 7*log(n) - 2 != 0 to avoid division by zero
            denominator3 = 7 * log_n - 2
            if abs(denominator3) < 1e-10:
                return 0
            term3 = (7 * log_n) / denominator3
            
            eta_value = term1 * term2 * term3
            return eta_value
    
    except (ValueError, ZeroDivisionError, OverflowError):
        return 0

def sum_eta_over_primes(omega_k, n, M):
    total = 0
    prime_count = 0
    prime_contributions = []
    
    for q in sympy.primerange(7, 10**7):  # Arbitrary large upper bound
        if q in [2, 3, 5]:
            continue  # Skip primes that divide M
        
        eta_value = eta(q, n, M)
        total += eta_value
        prime_contributions.append((q, eta_value))
        prime_count += 1
        
        # Stop once we've considered omega_k primes
        if prime_count >= omega_k:
            break
    
    return total, prime_count, prime_contributions

# Run the calculation
sum_result, num_primes, contributions = sum_eta_over_primes(omega_k, n, M)

# Display results
print(f"\nResults:")
print(f"Number of primes summed: {num_primes}")
print(f"Sum of η(q) over primes: {sum_result:.4f}")
print(f"Sum (more precision): {sum_result:.10f}")

"""
Output:

Number of primes summed: 26
Sum of η(q) over primes: 0.9616
Sum (more precision): 0.9615567361
"""