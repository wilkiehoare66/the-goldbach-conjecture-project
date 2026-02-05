#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <mutex>
#include <omp.h>
#include <atomic>
#include <iomanip>
#include <thread>
#include <numeric>

const int OMEGA_K_MAX = 41;
const long long MAX_PRIME_SIEVE = 100000000;

const int MAX_THREADS = 6;
const int CHUNKS_PER_BATCH = 20;
const int COOLDOWN_SECONDS = 30;
const bool ENABLE_COOLDOWN = true;

// For automatic resume after crash
const std::string CHECKPOINT_FILE = "lemma5_12_checkpoint.dat";

// ========================================================

// Global prime storage
std::vector<bool> is_prime_sieve;
std::vector<long long> prime_list;

// Thread-safe progress tracking
std::atomic<int> completed_chunks(0);
std::mutex output_mutex;
std::map<int, std::vector<long long>> all_results;

// Checkpoint management
struct Checkpoint {
    int last_completed_chunk;
    std::chrono::system_clock::time_point timestamp;
};

void saveCheckpoint(int chunk) {
    std::ofstream ckpt(CHECKPOINT_FILE, std::ios::binary);
    Checkpoint cp;
    cp.last_completed_chunk = chunk;
    cp.timestamp = std::chrono::system_clock::now();
    ckpt.write(reinterpret_cast<const char*>(&cp), sizeof(cp));
    ckpt.close();
}

int loadCheckpoint() {
    std::ifstream ckpt(CHECKPOINT_FILE, std::ios::binary);
    if (!ckpt) return 0;
    
    Checkpoint cp;
    ckpt.read(reinterpret_cast<char*>(&cp), sizeof(cp));
    ckpt.close();
    
    return cp.last_completed_chunk;
}

// Initialize prime sieve using Sieve of Eratosthenes
void initializePrimes() {
    std::cout << "Initializing prime sieve up to " << MAX_PRIME_SIEVE << "..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    is_prime_sieve.resize(MAX_PRIME_SIEVE, true);
    is_prime_sieve[0] = is_prime_sieve[1] = false;
    
    for (long long i = 2; i * i < MAX_PRIME_SIEVE; i++) {
        if (is_prime_sieve[i]) {
            for (long long j = i * i; j < MAX_PRIME_SIEVE; j += i) {
                is_prime_sieve[j] = false;
            }
        }
    }
    
    prime_list.reserve(5761455);
    for (long long i = 2; i < MAX_PRIME_SIEVE; i++) {
        if (is_prime_sieve[i]) {
            prime_list.push_back(i);
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "Prime sieve complete in " << duration.count() << " seconds. Found " 
              << prime_list.size() << " primes." << std::endl;
}

bool isPrime(long long n) {
    if (n < MAX_PRIME_SIEVE) {
        return is_prime_sieve[n];
    }
    
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (long long i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}

long long nextPrime(long long n) {
    if (n < 2) return 2;
    n++;
    
    if (n < MAX_PRIME_SIEVE) {
        auto it = std::lower_bound(prime_list.begin(), prime_list.end(), n);
        if (it != prime_list.end()) {
            return *it;
        }
    }
    
    while (!isPrime(n)) n++;
    return n;
}

long long prevPrime(long long n) {
    if (n <= 2) return 2;
    n--;
    
    if (n < MAX_PRIME_SIEVE) {
        auto it = std::upper_bound(prime_list.begin(), prime_list.end(), n);
        if (it != prime_list.begin()) {
            return *(--it);
        }
    }
    
    while (n > 2 && !isPrime(n)) n--;
    return n;
}

thread_local std::unordered_map<long long, std::map<long long, int>> thread_factorCache;

std::map<long long, int> factorize(long long n) {
    if (thread_factorCache.count(n)) {
        return thread_factorCache[n];
    }
    
    std::map<long long, int> factors;
    if (n <= 1) {
        thread_factorCache[n] = factors;
        return factors;
    }
    
    long long original = n;
    
    for (long long p : prime_list) {
        if (p * p > n) break;
        
        while (n % p == 0) {
            factors[p]++;
            n /= p;
        }
        
        if (n == 1) break;
    }
    
    if (n > 1) {
        factors[n]++;
    }
    
    thread_factorCache[original] = factors;
    
    if (thread_factorCache.size() > 100000) {
        thread_factorCache.clear();
    }
    
    return factors;
}

bool squarefreeOdds(long long a) {
    auto factors = factorize(a);
    for (const auto& [prime, exponent] : factors) {
        if (prime > 2 && exponent > 1) {
            return false;
        }
    }
    return true;
}

static inline long long odd_part(long long x) {
    if (x == 0) return 0;
    unsigned long long ux = static_cast<unsigned long long>(x);
    return x >> __builtin_ctzll(ux);
}

bool coprimeWeaker(long long a, long long b) {
    return std::gcd(odd_part(a), odd_part(b)) == 1;
}

std::set<std::tuple<long long, long long, long long>> threeProgressions() {
    std::set<std::tuple<long long, long long, long long>> options;
    long long p = 3;
    
    std::cout << "\nComputing three progressions (this is single-threaded)..." << std::endl;
    int counter = 0;
    
    while (p < 5000) {
        if (++counter % 100 == 0) {
            std::cout << "Progress: " << p << " / 5000" << std::endl;
        }
        
        long long q = nextPrime(p);
        
        while (q < p + 500) {
            if (squarefreeOdds(p + q)) {
                bool allCoprime = true;
                for (const auto& [x, y, z] : options) {
                    if (!coprimeWeaker(p + q, z)) {
                        allCoprime = false;
                        break;
                    }
                }
                if (allCoprime) {
                    options.insert({p, q, p + q});
                }
            }
            q = nextPrime(q);
        }
        p = nextPrime(p);
    }
    
    std::cout << "Three progressions complete. Found " << options.size() << " progressions." << std::endl;
    return options;
}

bool squarefree(long long a) {
    auto factors = factorize(a);
    for (const auto& [prime, exponent] : factors) {
        if (exponent > 1) {
            return false;
        }
    }
    return true;
}

bool check(long long w) {
    std::vector<long long> reps;
    long long pr = 2;
    
    while (pr < w) {
        if (squarefree(w - pr)) {
            bool allCoprime = true;
            for (long long j : reps) {
                if (!coprimeWeaker(w - pr, j)) {
                    allCoprime = false;
                    break;
                }
            }
            if (allCoprime) {
                reps.push_back(w - pr);
                if (reps.size() == OMEGA_K_MAX) {
                    return true;
                }
            }
        }
        pr = nextPrime(pr);
    }
    return false;
}

std::vector<long long> proc(long long fromme, long long upto, 
                            const std::set<std::tuple<long long, long long, long long>>& preloads) {
    std::map<long long, int> representations;
    long long prime = prevPrime(fromme - 1000);
    
    while (prime < upto) {
        for (const auto& [a, b, c] : preloads) {
            int e = 0;
            long long m = prime + c * (1LL << e);
            while (m < upto) {
                representations[m]++;
                e++;
                m = prime + c * (1LL << e);
            }
        }
        prime = nextPrime(prime + 500);
    }
    
    std::vector<long long> exceptions;
    long long w = fromme;
    
    while (w <= upto) {
        if (representations[w] < OMEGA_K_MAX) {
            if (!check(w)) {
                exceptions.push_back(w);
            }
        }
        w += 2;
    }
    
    return exceptions;
}

void saveResults(const std::string& prefix) {
    std::ofstream csv(prefix + "_results.csv");
    csv << "Chunk,Range_Start,Range_End,Exception_Count,Exceptions\n";
    
    std::ofstream log(prefix + "_results.txt");
    log << "LEMMA 5.12 COMPUTATION RESULTS (LAPTOP-SAFE VERSION)\n";
    log << "OMEGA_K_MAX = " << OMEGA_K_MAX << "\n";
    log << "Max threads used: " << MAX_THREADS << "\n";
    log << "=====================================\n\n";
    
    for (const auto& [ell, exceptions] : all_results) {
        long long range_start = (ell == 1) ? 100001 : (ell - 1) * 1000000 + 1;
        long long range_end = ell * 1000000;
        
        // CSV
        csv << ell << "," 
            << range_start << ","
            << range_end << ","
            << exceptions.size() << ",\"";
        
        for (size_t i = 0; i < exceptions.size(); i++) {
            csv << exceptions[i];
            if (i < exceptions.size() - 1) csv << ", ";
        }
        csv << "\"\n";
        
        // Log
        log << ell << " [";
        for (size_t i = 0; i < exceptions.size(); i++) {
            log << exceptions[i];
            if (i < exceptions.size() - 1) log << ", ";
        }
        log << "]\n";
    }
    
    csv.close();
    log.close();
}

int main() {
    auto total_start = std::chrono::high_resolution_clock::now();
    
    // Set thread limit for laptop safety
    omp_set_num_threads(MAX_THREADS);
    
    std::cout << "========================================\n";
    std::cout << "LAPTOP-SAFE VERSION - THERMAL PROTECTION\n";
    std::cout << "========================================\n";
    std::cout << "Max threads: " << MAX_THREADS << " (reduce if still overheating)\n";
    std::cout << "Batch size: " << CHUNKS_PER_BATCH << " chunks\n";
    std::cout << "Cooldown: " << (ENABLE_COOLDOWN ? std::to_string(COOLDOWN_SECONDS) + " seconds" : "disabled") << "\n\n";
    
    // Initialize prime sieve
    initializePrimes();
    
    // Compute three progressions
    auto tp = threeProgressions();
    
    // Check for previous checkpoint
    int start_chunk = loadCheckpoint();
    if (start_chunk > 0) {
        std::cout << "\nResuming from chunk " << start_chunk << " (found checkpoint)\n";
        completed_chunks = start_chunk;
    }
    
    const long long max_ell = 8000;
    std::cout << "\nProcessing chunks " << (start_chunk + 1) << " to " << max_ell << " in batches...\n\n";
    
    auto processing_start = std::chrono::high_resolution_clock::now();
    
    // Process in batches with cooling periods
    for (int batch_start = start_chunk + 1; batch_start <= max_ell; batch_start += CHUNKS_PER_BATCH) {
        int batch_end = std::min((long long)(batch_start + CHUNKS_PER_BATCH - 1), max_ell);
        
        std::cout << "Processing batch: chunks " << batch_start << "-" << batch_end << std::endl;
        
        #pragma omp parallel for schedule(dynamic)
        for (int ell = batch_start; ell <= batch_end; ell++) {
            long long range_start = (ell == 1) ? 100001 : (ell - 1) * 1000000 + 1;
            long long range_end = ell * 1000000;
            
            auto exceptions = proc(range_start, range_end, tp);
            
            {
                std::lock_guard<std::mutex> lock(output_mutex);
                all_results[ell] = exceptions;
                completed_chunks++;
                
                if (completed_chunks % 5 == 0) {
                    std::cout << "  Completed: " << completed_chunks << "/" << max_ell 
                              << " (" << (completed_chunks * 100 / max_ell) << "%)" << std::endl;
                    
                    // Save checkpoint
                    saveCheckpoint(completed_chunks);
                    
                    // Periodic save of results
                    if (completed_chunks % 50 == 0) {
                        saveResults("lemma5_10_partial");
                        std::cout << "  Partial results saved." << std::endl;
                    }
                }
            }
        }
        
        // Cooling period between batches
        if (ENABLE_COOLDOWN && batch_end < max_ell) {
            std::cout << "Cooling down for " << COOLDOWN_SECONDS << " seconds..." << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(COOLDOWN_SECONDS));
        }
    }
    
    auto processing_end = std::chrono::high_resolution_clock::now();
    auto processing_duration = std::chrono::duration_cast<std::chrono::seconds>(processing_end - processing_start);
    
    std::cout << "\nProcessing complete! Saving final results..." << std::endl;
    
    // Save final results
    saveResults("lemma5_10_final");
    
    // Summary
    int total_exceptions = 0;
    for (const auto& [ell, exceptions] : all_results) {
        total_exceptions += exceptions.size();
    }
    
    auto total_end = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(total_end - total_start);
    
    std::cout << "\n=====================================\n";
    std::cout << "COMPUTATION COMPLETE\n";
    std::cout << "=====================================\n";
    std::cout << "Total chunks processed: " << completed_chunks << "\n";
    std::cout << "Total exceptions found: " << total_exceptions << "\n";
    std::cout << "Processing time: " << processing_duration.count() / 3600 << " hours\n";
    std::cout << "Total time: " << total_duration.count() / 3600 << " hours\n";
    std::cout << "\nResults saved to:\n";
    std::cout << "  - lemma5_10_final_results.csv\n";
    std::cout << "  - lemma5_10_final_results.txt\n";
    std::cout << "\nCheckpoint file can be deleted: " << CHECKPOINT_FILE << "\n";
    
    return 0;
}


