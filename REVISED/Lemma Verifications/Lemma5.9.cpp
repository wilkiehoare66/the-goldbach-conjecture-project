#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <omp.h>

using namespace std;
using i64 = int64_t;
using u64 = uint64_t;

// Configuration
const i64 RANGE_START = 4810000000LL;       // 4.81 * 10^9 (verified up to here)
const i64 RANGE_END = 2000000000000LL;      // 2 * 10^12
const i64 CHUNK_SIZE = 10000000LL;          // 10^7 (same as original)
const i64 SMALL_PRIME_LIMIT = 500;          // Small primes < 500 (same as original)

const string RESULTS_FILE = "lemma59_results.csv";
const string CHECKPOINT_FILE = "lemma59_checkpoint.csv";

const int NUM_THREADS = 6;  // Number of parallel threads to use

// Primality Testing (Miller-Rabin, deterministic for n < 3.3 * 10^24)
u64 mod_pow(u64 base, u64 exp, u64 mod) {
    u64 result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) {
            result = (__uint128_t)result * base % mod;
        }
        exp >>= 1;
        base = (__uint128_t)base * base % mod;
    }
    return result;
}

bool miller_rabin_witness(u64 n, u64 a) {
    if (n % a == 0) return n == a;
    
    u64 d = n - 1;
    int r = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        r++;
    }
    
    u64 x = mod_pow(a, d, n);
    if (x == 1 || x == n - 1) return true;
    
    for (int i = 0; i < r - 1; i++) {
        x = (__uint128_t)x * x % n;
        if (x == n - 1) return true;
    }
    return false;
}

bool is_prime(u64 n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0) return false;
    
    // Deterministic witnesses for n < 3,317,044,064,679,887,385,961,981
    static const u64 witnesses[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    for (u64 a : witnesses) {
        if (n == a) return true;
        if (!miller_rabin_witness(n, a)) return false;
    }
    return true;
}

// Prime Navigation
u64 next_prime(u64 n) {
    if (n < 2) return 2;
    n++;
    if (n == 2) return 2;
    if (n % 2 == 0) n++;
    while (!is_prime(n)) n += 2;
    return n;
}

u64 prev_prime(u64 n) {
    if (n <= 2) return 0;
    n--;
    if (n == 2) return 2;
    if (n % 2 == 0) n--;
    while (!is_prime(n)) n -= 2;
    return n;
}

// Segmented Sieve for Fast Prime Generation

// Simple sieve for small primes (used to initialize segmented sieve)
vector<i64> sieve_small_primes(i64 limit) {
    vector<bool> is_prime_arr(limit + 1, true);
    is_prime_arr[0] = is_prime_arr[1] = false;
    
    for (i64 i = 2; i * i <= limit; i++) {
        if (is_prime_arr[i]) {
            for (i64 j = i * i; j <= limit; j += i) {
                is_prime_arr[j] = false;
            }
        }
    }
    
    vector<i64> primes;
    for (i64 i = 2; i <= limit; i++) {
        if (is_prime_arr[i]) primes.push_back(i);
    }
    return primes;
}

// Segmented sieve: get all primes in [from, to]
vector<i64> segmented_sieve(i64 from, i64 to, const vector<i64>& small_primes) {
    if (from < 2) from = 2;
    i64 size = to - from + 1;
    vector<bool> is_prime_seg(size, true);
    
    // Use small primes to sieve the segment
    for (i64 p : small_primes) {
        if (p * p > to) break;
        
        // Find first multiple of p in [from, to]
        i64 start = ((from + p - 1) / p) * p;
        if (start == p) start += p;  // Don't mark the prime itself
        
        for (i64 j = start; j <= to; j += p) {
            is_prime_seg[j - from] = false;
        }
    }
    
    // Collect primes
    vector<i64> primes;
    primes.reserve(size / 20);  // Approximate density
    for (i64 i = 0; i < size; i++) {
        if (is_prime_seg[i]) {
            primes.push_back(from + i);
        }
    }
    return primes;
}

// Global: small primes for sieving (primes up to sqrt(2*10^12) ≈ 1.5*10^6)
vector<i64> SIEVE_PRIMES;

// Global: primes for fallback iteration (primes from 503 to 10^7)
// Using precomputed array instead of next_prime() gives ~3x speedup on fallback
vector<i64> FALLBACK_PRIMES;

void generate_sieve_primes() {
    // Need primes up to sqrt(RANGE_END) for segmented sieve
    i64 limit = (i64)sqrt((double)RANGE_END) + 1000;
    SIEVE_PRIMES = sieve_small_primes(limit);
}

void generate_fallback_primes() {
    // Generate primes from 503 to 10^7 for fallback iteration
    // This avoids repeated next_prime() calls which are expensive
    const i64 FALLBACK_LIMIT = 10000000;  // 10^7 - more than enough to find 2 reps
    auto all_primes = sieve_small_primes(FALLBACK_LIMIT);
    for (i64 p : all_primes) {
        if (p >= 503) FALLBACK_PRIMES.push_back(p);
    }
}

// Square-free and Coprimality Checks (for is_good verification)
bool is_squarefree(u64 n) {
    if (n <= 1) return true;
    
    int count = 0;
    while (n % 2 == 0) {
        count++;
        if (count > 1) return false;
        n /= 2;
    }
    
    for (u64 i = 3; i * i <= n; i += 2) {
        count = 0;
        while (n % i == 0) {
            count++;
            if (count > 1) return false;
            n /= i;
        }
    }
    return true;
}

// Get prime factors of n (just the primes, not their multiplicities)
vector<i64> prime_factors(i64 n) {
    vector<i64> factors;
    if (n <= 1) return factors;
    
    while (n % 2 == 0) {
        if (factors.empty() || factors.back() != 2) {
            factors.push_back(2);
        }
        n /= 2;
    }
    
    for (i64 i = 3; i * i <= n; i += 2) {
        while (n % i == 0) {
            if (factors.empty() || factors.back() != i) {
                factors.push_back(i);
            }
            n /= i;
        }
    }
    
    if (n > 1) {
        factors.push_back(n);
    }
    return factors;
}

// Check if a and b are coprime (share no prime factors)
bool coprime(i64 a, i64 b) {
    vector<i64> fa = prime_factors(a);
    vector<i64> fb = prime_factors(b);
    
    for (i64 p : fa) {
        for (i64 q : fb) {
            if (p == q) return false;
        }
    }
    return true;
}

// Check if n is coprime to all elements in coprime_to
bool coprime_to_all(i64 n, const vector<i64>& coprime_to) {
    for (i64 c : coprime_to) {
        if (!coprime(n, c)) return false;
    }
    return true;
}

// Find next square-free after y that is coprime to all elements in coprime_to
i64 next_coprime_squarefree(i64 y, const vector<i64>& coprime_to) {
    i64 i = 1;
    while (true) {
        i64 candidate = y + i;
        if (is_squarefree(candidate) && coprime_to_all(candidate, coprime_to)) {
            return candidate;
        }
        i++;
    }
}

// Check if n can be written as p + s in at least 4 ways where:
// - p is prime
// - s is square-free  
// - all the s values are pairwise coprime
// Returns: {is_good, list of square-free witnesses}
pair<bool, vector<i64>> is_good(i64 n) {
    i64 sf = 2;  // First nontrivial square-free
    vector<i64> cpt;  // Coprime witnesses found
    
    while (sf < n) {
        if (is_prime(n - sf)) {
            cpt.push_back(sf);
        }
        if (cpt.size() < 4) {
            sf = next_coprime_squarefree(sf, cpt);
        } else {
            return {true, cpt};
        }
    }
    return {false, cpt};
}

// Global Small Primes (generated once)
vector<i64> SMALL_PRIMES;

void generate_small_primes() {
    i64 p = 2;
    while (p < SMALL_PRIME_LIMIT) {
        SMALL_PRIMES.push_back(p);
        p = next_prime(p);
    }
}

// Core Verification Function (direct translation of verify_ft_semiprimes_quick)
struct ChunkResult {
    i64 chunk_start;
    i64 chunk_end;
    vector<i64> exceptions;
    i64 checked_count;
    double elapsed_seconds;
};

ChunkResult verify_chunk(i64 from_me, i64 to_me) {
    ChunkResult result;
    result.chunk_start = from_me;
    result.chunk_end = to_me;
    result.checked_count = 0;
    
    auto t_start = chrono::high_resolution_clock::now();
    
    // Generate primes in range using segmented sieve (MUCH faster than next_prime iteration)
    i64 prime_start = from_me - SMALL_PRIME_LIMIT;
    if (prime_start < 2) prime_start = 2;
    vector<i64> pc_primes = segmented_sieve(prime_start, to_me, SIEVE_PRIMES);
    
    // Quick sieve: count representations using small primes
    unordered_map<i64, int> num_representations;
    num_representations.reserve(CHUNK_SIZE + 1000);
    
    for (i64 q : pc_primes) {
        for (i64 sp : SMALL_PRIMES) {
            i64 m = sp + q;
            if (m >= from_me && m <= to_me && m % 2 == 0) {
                num_representations[m]++;
            }
        }
    }
    
    // Check even numbers in range
    i64 m = from_me;
    if (m % 2 != 0) m++;
    
    while (m <= to_me) {
        result.checked_count++;
        
        auto it = num_representations.find(m);
        int sieve_reps = (it != num_representations.end()) ? it->second : 0;
        
        if (sieve_reps < 2) {
            // Fallback: use direct is_prime() with precomputed prime list
            // The sieve already checked small primes < 500 with primes in chunk range
            // FALLBACK_PRIMES contains primes from 503 onwards (precomputed at startup)
            // Using array iteration instead of next_prime() gives ~3x speedup
            int reps_needed = 2 - sieve_reps;
            int reps_found = 0;
            
            for (i64 q : FALLBACK_PRIMES) {
                if (2 * q >= m) break;  // No point checking if q > m/2
                if (is_prime(m - q)) {
                    reps_found++;
                    if (reps_found >= reps_needed) break;
                }
            }
            
            if (reps_found < reps_needed) {
                result.exceptions.push_back(m);
            }
        }
        m += 2;
    }
    
    auto t_end = chrono::high_resolution_clock::now();
    result.elapsed_seconds = chrono::duration<double>(t_end - t_start).count();
    
    return result;
}

// CSV I/O
string get_timestamp() {
    auto now = chrono::system_clock::now();
    time_t t = chrono::system_clock::to_time_t(now);
    char buf[64];
    strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", localtime(&t));
    return string(buf);
}

void init_results_file() {
    ifstream check(RESULTS_FILE);
    if (!check.good()) {
        ofstream f(RESULTS_FILE);
        f << "timestamp,chunk_start,chunk_end,status,exceptions,checked_count,elapsed_seconds" << endl;
        f.close();
    }
}

void append_result(const ChunkResult& r) {
    ofstream f(RESULTS_FILE, ios::app);
    
    string status = r.exceptions.empty() ? "OK" : "EXCEPTION";
    string exc_str = "";
    for (size_t i = 0; i < r.exceptions.size(); i++) {
        if (i > 0) exc_str += ";";
        exc_str += to_string(r.exceptions[i]);
    }
    if (exc_str.empty()) exc_str = "none";
    
    f << get_timestamp() << ","
      << r.chunk_start << ","
      << r.chunk_end << ","
      << status << ","
      << exc_str << ","
      << r.checked_count << ","
      << fixed << setprecision(2) << r.elapsed_seconds
      << endl;
    
    f.close();
}

set<pair<i64,i64>> load_completed_chunks() {
    set<pair<i64,i64>> completed;
    
    ifstream f(RESULTS_FILE);
    if (!f.good()) return completed;
    
    string line;
    getline(f, line); // Skip header
    
    while (getline(f, line)) {
        if (line.empty()) continue;
        
        stringstream ss(line);
        string timestamp, chunk_start_str, chunk_end_str;
        
        getline(ss, timestamp, ',');
        getline(ss, chunk_start_str, ',');
        getline(ss, chunk_end_str, ',');
        
        i64 cs = stoll(chunk_start_str);
        i64 ce = stoll(chunk_end_str);
        completed.insert({cs, ce});
    }
    
    return completed;
}

// Progress Display
string format_time(double seconds) {
    int h = (int)(seconds / 3600);
    int m = (int)(fmod(seconds, 3600) / 60);
    int s = (int)(fmod(seconds, 60));
    
    stringstream ss;
    if (h > 0) {
        ss << h << "h " << m << "m";
    } else if (m > 0) {
        ss << m << "m " << s << "s";
    } else {
        ss << fixed << setprecision(1) << seconds << "s";
    }
    return ss.str();
}

string format_number(i64 n) {
    if (n >= 1e12) return to_string(n / (i64)1e9) + "B";
    if (n >= 1e9) return to_string((double)n / 1e9).substr(0, 5) + "B";
    if (n >= 1e6) return to_string((double)n / 1e6).substr(0, 5) + "M";
    return to_string(n);
}

// Main
int main(int argc, char* argv[]) {
    int num_threads = NUM_THREADS;
    bool resume = false;
    bool test_mode = false;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "--resume") {
            resume = true;
        } else if (arg == "--test") {
            test_mode = true;
        } else if (arg == "--help") {
            cout << "Usage: " << argv[0] << " [options]" << endl;
            cout << "  --resume       Resume from checkpoint" << endl;
            cout << "  --test         Run single chunk test" << endl;
            return 0;
        }
    }
    
    omp_set_num_threads(num_threads);
    
    // Header
    cout << "============================================================" << endl;
    cout << "Lemma 5.9 Verification (C++ with OpenMP)" << endl;
    cout << "============================================================" << endl;
    cout << "Range: " << format_number(RANGE_START) << " to " << format_number(RANGE_END) << endl;
    cout << "Chunk size: " << format_number(CHUNK_SIZE) << endl;
    cout << "Threads: " << num_threads << endl;
    cout << "Results file: " << RESULTS_FILE << endl;
    cout << "============================================================" << endl;
    
    // Generate small primes for quick sieve
    cout << "Generating small primes < " << SMALL_PRIME_LIMIT << "..." << flush;
    generate_small_primes();
    cout << " done (" << SMALL_PRIMES.size() << " primes)" << endl;
    
    // Generate sieve primes for segmented sieve
    cout << "Generating sieve primes up to sqrt(" << RANGE_END << ")..." << flush;
    generate_sieve_primes();
    cout << " done (" << SIEVE_PRIMES.size() << " primes)" << endl;
    
    // Generate fallback primes for fast fallback iteration
    cout << "Generating fallback primes (503 to 10^7)..." << flush;
    generate_fallback_primes();
    cout << " done (" << FALLBACK_PRIMES.size() << " primes)" << endl;
    
    // Test mode
    if (test_mode) {
        cout << "\n=== TEST MODE ===" << endl;
        
        // Test is_good() with known case from original Python
        cout << "Testing is_good(740000138)..." << endl;
        auto [good, witnesses] = is_good(740000138);
        cout << "  Result: " << (good ? "True" : "False") << endl;
        cout << "  Witnesses: ";
        for (size_t i = 0; i < witnesses.size(); i++) {
            if (i > 0) cout << ", ";
            cout << witnesses[i];
        }
        cout << endl;
        cout << "  Expected: True, [21, 235, 247, 391]" << endl;
        
        cout << "\nRunning single chunk: " << RANGE_START << " to " << (RANGE_START + CHUNK_SIZE - 1) << endl;
        
        ChunkResult r = verify_chunk(RANGE_START, RANGE_START + CHUNK_SIZE - 1);
        
        cout << "Time: " << r.elapsed_seconds << "s" << endl;
        cout << "Even integers checked: " << r.checked_count << endl;
        cout << "Exceptions: " << (r.exceptions.empty() ? "none" : to_string(r.exceptions.size())) << endl;
        if (!r.exceptions.empty()) {
            cout << "  Values: ";
            for (i64 e : r.exceptions) cout << e << " ";
            cout << endl;
        }
        
        i64 total_chunks = (RANGE_END - RANGE_START + CHUNK_SIZE - 1) / CHUNK_SIZE;
        double est_hours = (total_chunks * r.elapsed_seconds) / 3600.0 / num_threads;
        cout << "\nEstimated total time with " << num_threads << " threads: " 
             << fixed << setprecision(1) << est_hours << " hours" << endl;
        
        return 0;
    }
    
    // Build chunk list
    vector<pair<i64, i64>> all_chunks;
    for (i64 cs = RANGE_START; cs < RANGE_END; cs += CHUNK_SIZE) {
        i64 ce = min(cs + CHUNK_SIZE - 1, RANGE_END);
        all_chunks.push_back({cs, ce});
    }
    
    i64 total_chunks = all_chunks.size();
    cout << "Total chunks: " << total_chunks << endl;
    
    // Load completed chunks if resuming
    set<pair<i64,i64>> completed;
    if (resume) {
        completed = load_completed_chunks();
        cout << "Resuming: " << completed.size() << " chunks already completed" << endl;
    }
    
    // Filter to remaining chunks
    vector<pair<i64, i64>> chunks;
    for (auto& c : all_chunks) {
        if (completed.find(c) == completed.end()) {
            chunks.push_back(c);
        }
    }
    
    i64 remaining_chunks = chunks.size();
    cout << "Remaining: " << remaining_chunks << " chunks" << endl;
    
    if (remaining_chunks == 0) {
        cout << "\nAll chunks already verified!" << endl;
        return 0;
    }
    
    // Initialize results file
    init_results_file();
    
    cout << "\nStarting verification..." << endl;
    cout << "------------------------------------------------------------" << endl;
    
    // Track progress
    i64 completed_count = 0;
    double total_time = 0;
    vector<i64> all_exceptions;
    auto global_start = chrono::high_resolution_clock::now();
    
    // Process chunks - sequential processing but with thread pool
    // Process in batches to avoid memory issues
    const size_t BATCH_SIZE = num_threads * 2;
    
    for (size_t batch_start = 0; batch_start < chunks.size(); batch_start += BATCH_SIZE) {
        size_t batch_end = min(batch_start + BATCH_SIZE, chunks.size());
        size_t batch_count = batch_end - batch_start;
        
        // Store results for this batch
        vector<ChunkResult> batch_results(batch_count);
        
        // Process batch in parallel
        #pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < batch_count; i++) {
            size_t idx = batch_start + i;
            i64 cs = chunks[idx].first;
            i64 ce = chunks[idx].second;
            batch_results[i] = verify_chunk(cs, ce);
        }
        
        // Process results sequentially (thread-safe)
        for (size_t i = 0; i < batch_count; i++) {
            ChunkResult& r = batch_results[i];
            completed_count++;
            total_time += r.elapsed_seconds;
            
            // Append to CSV
            append_result(r);
            
            // Collect exceptions
            for (i64 e : r.exceptions) {
                all_exceptions.push_back(e);
            }
            
            // Progress output
            double avg_time = total_time / completed_count;
            double eta_seconds = (remaining_chunks - completed_count) * avg_time / num_threads;
            
            cout << "[" << completed_count << "/" << remaining_chunks << "] "
                 << format_number(r.chunk_start) << "-" << format_number(r.chunk_end) << ": ";
            
            if (r.exceptions.empty()) {
                cout << "OK";
            } else {
                cout << "EXCEPTION(";
                for (size_t j = 0; j < r.exceptions.size(); j++) {
                    if (j > 0) cout << ",";
                    cout << r.exceptions[j];
                }
                cout << ")";
            }
            
            cout << " [" << r.checked_count << " checked, " 
                 << fixed << setprecision(1) << r.elapsed_seconds << "s]"
                 << " ETA: " << format_time(eta_seconds) << endl;
        }
    }
    
    auto global_end = chrono::high_resolution_clock::now();
    double wall_time = chrono::duration<double>(global_end - global_start).count();
    
    cout << "------------------------------------------------------------" << endl;
    cout << "\nVerification complete!" << endl;
    cout << "Wall time: " << format_time(wall_time) << endl;
    cout << "Total CPU time: " << format_time(total_time) << endl;
    cout << "Results saved to: " << RESULTS_FILE << endl;
    
    if (all_exceptions.empty()) {
        cout << "\n*** NO EXCEPTIONS FOUND ***" << endl;
        cout << "All even integers in [" << format_number(RANGE_START) << ", " 
             << format_number(RANGE_END) << "] have at least 2 Goldbach representations." << endl;
    } else {
        cout << "\n*** EXCEPTIONS FOUND BY QUICK CHECK: " << all_exceptions.size() << " ***" << endl;
        cout << "Verifying with full is_good() check..." << endl;
        cout << "------------------------------------------------------------" << endl;
        
        vector<i64> true_exceptions;
        for (i64 e : all_exceptions) {
            auto [good, witnesses] = is_good(e);
            cout << e << ": ";
            if (good) {
                cout << "VERIFIED OK (witnesses: ";
                for (size_t i = 0; i < witnesses.size(); i++) {
                    if (i > 0) cout << ", ";
                    cout << witnesses[i];
                }
                cout << ")" << endl;
            } else {
                cout << "TRUE EXCEPTION (only " << witnesses.size() << " witnesses found)" << endl;
                true_exceptions.push_back(e);
            }
        }
        
        cout << "------------------------------------------------------------" << endl;
        if (true_exceptions.empty()) {
            cout << "\n*** ALL EXCEPTIONS VERIFIED - NO TRUE EXCEPTIONS ***" << endl;
        } else {
            cout << "\n*** TRUE EXCEPTIONS: " << true_exceptions.size() << " ***" << endl;
            for (i64 e : true_exceptions) {
                cout << "  " << e << endl;
            }
        }
    }
    
    return 0;
}