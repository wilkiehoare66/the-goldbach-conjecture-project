import sympy
import math

n = 8 * 10**9
omega_k = 91
M = 210

def eta(q, n, M):
    if q <= 1:
        return 0
    
    L = sympy.lcm(M, q)
    
    try:
        log_n = math.log(n)
        
        # Case 1: L ≤ 10^5 (from Lemma 5.2)
        if L <= 10**5:
            # Use the Case 1 bound uniformly for all primes with L ≤ 10^5
            numerator = 70 * log_n + 21
            denominator = 70 * log_n - 4
            
            # Check for division by zero
            if abs(denominator) < 1e-10:
                return 0
            
            eta_value = (1 / (q - 1)) * (numerator / denominator)
            return eta_value

        # Case 2: L > 10^5 (from Lemma 5.2)
        else:
            # Check the condition n > Mq
            if n <= M * q:
                return 0
            
            # Check if n/L > 1 to avoid log of negative/zero
            if n/L <= 1:
                return 0
            
            # Calculate each component
            term1 = 2 / (q - 1)
            term2 = log_n / math.log(n / L)
            
            # Check if 35*log(n) - 2 != 0 to avoid division by zero
            denominator3 = 35 * log_n - 2
            if abs(denominator3) < 1e-10:
                return 0
            term3 = (35 * log_n) / denominator3
            
            eta_value = term1 * term2 * term3
            return eta_value
    
    except (ValueError, ZeroDivisionError, OverflowError):
        return 0

def sum_eta_over_primes(omega_k, n, M):
    total = 0
    prime_count = 0
    prime_contributions = []
    
    # Find primes that do not divide M (since those are the q we care about)
    q = 2
    while True:
        q = sympy.nextprime(q)
        
        # Skip primes dividing M
        if M % q == 0:
            continue
        
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
print(f"Results:")
print(f"Number of primes summed: {num_primes}")
print(f"Sum of η(q) over primes: {sum_result:.4f}")
print(f"Sum (more precision): {sum_result:.10f}")

"""
Output:

Results:
Number of primes summed: 91
Sum of η(q) over primes: 0.9928
Sum (more precision): 0.9928399840
"""
