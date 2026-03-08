# In:
from sympy import prevprime, primerange, nextprime, factorint, isprime

def coprime(a,b):
    fi_A = factorint(a)
    fi_B = factorint(b)
    for p in fi_A:
        if p in fi_B.keys():
            return False
    return True

def nextcoprimesquarefree(y,coprime_to):
    i = 1
    while i > 0:
        fi = factorint(y + i)
        condition_met = [coprime(j,y+i) for j in coprime_to]
        if max(fi.values()) == 1 and False not in condition_met:
            return y + i
        i += 1

def is_good(n):
    sf = 2 #the first nontrivial squarefree
    cpt = []
    while sf < n:
        if isprime(n-sf) == True:
            cpt.append(sf)
        if len(cpt) < 4:
            sf = nextcoprimesquarefree(sf, cpt)
        else:
            return True, cpt
    return False, cpt

def verify_ft_semiprimes_small(from_me, to_me):
    not_verified = []
    m = from_me
    while m <= to_me+1:
        ig = is_good(m)
        if ig[0] == False:
            print(m)
            not_verified.append(m)
        m += 2
    return not_verified

def verify_ft_semiprimes_quick(from_me, to_me):
    small_primes = set()
    p = 2
    while p < 500:
        small_primes.add(p)
        p = nextprime(p)
    
    pc_primes = set()
    p = prevprime(from_me - 10**3)
    while p < to_me:
        pc_primes.add(p)
        p = nextprime(p)
    
    num_representations = {}
    for q in pc_primes:
        for p in small_primes:
            m = p + q
            if m not in num_representations:
                num_representations[m] = 1
            else:
                num_representations[m] += 1
    #print("Generated!")
    
    not_verified = []
    m = from_me
    while m <= to_me:
        if m not in num_representations or num_representations[m] < 2:
            reps = []
            q = 2
            while 2*q < m and len(reps) < 2:
                if m-q in pc_primes:
                    #print(q,m-q)
                    reps.append((q,m-q))
                q = nextprime(q)
            #print(m, reps)
            if len(reps) < 2:
                not_verified.append(m)
        m += 2
    return not_verified

all_exceptions = []

for l in range(481, 800):
    exptons = verify_ft_semiprimes_quick(l*10**7, (l+1)*10**7)
    for e in exptons:
        all_exceptions.append(e)
    print(l, l+1, all_exceptions)

print(all_exceptions)

# Out:
"""
4
6
8
10
12
14
[4, 6, 8, 10, 12, 14]
0 1 [4, 6, 8, 10, 12, 14]
1 2 [4, 6, 8, 10, 12, 14]
2 3 [4, 6, 8, 10, 12, 14]
3 4 [4, 6, 8, 10, 12, 14]
4 5 [4, 6, 8, 10, 12, 14]
5 6 [4, 6, 8, 10, 12, 14]
6 7 [4, 6, 8, 10, 12, 14]
7 8 [4, 6, 8, 10, 12, 14]
8 9 [4, 6, 8, 10, 12, 14]
9 10 [4, 6, 8, 10, 12, 14]
10 11 [4, 6, 8, 10, 12, 14]
11 12 [4, 6, 8, 10, 12, 14]
12 13 [4, 6, 8, 10, 12, 14]
13 14 [4, 6, 8, 10, 12, 14]
14 15 [4, 6, 8, 10, 12, 14]
15 16 [4, 6, 8, 10, 12, 14]
16 17 [4, 6, 8, 10, 12, 14]
17 18 [4, 6, 8, 10, 12, 14]
18 19 [4, 6, 8, 10, 12, 14]
19 20 [4, 6, 8, 10, 12, 14]
20 21 [4, 6, 8, 10, 12, 14]
21 22 [4, 6, 8, 10, 12, 14]
22 23 [4, 6, 8, 10, 12, 14]
23 24 [4, 6, 8, 10, 12, 14]
24 25 [4, 6, 8, 10, 12, 14]
25 26 [4, 6, 8, 10, 12, 14]
26 27 [4, 6, 8, 10, 12, 14]
27 28 [4, 6, 8, 10, 12, 14]
28 29 [4, 6, 8, 10, 12, 14]
29 30 [4, 6, 8, 10, 12, 14]
30 31 [4, 6, 8, 10, 12, 14]
31 32 [4, 6, 8, 10, 12, 14]
...
71 72 [4, 6, 8, 10, 12, 14]
72 73 [4, 6, 8, 10, 12, 14]
73 74 [4, 6, 8, 10, 12, 14]
74 75 [4, 6, 8, 10, 12, 14, 740000138]
75 76 [4, 6, 8, 10, 12, 14, 740000138]
76 77 [4, 6, 8, 10, 12, 14, 740000138]
77 78 [4, 6, 8, 10, 12, 14, 740000138]
78 79 [4, 6, 8, 10, 12, 14, 740000138]
79 80 [4, 6, 8, 10, 12, 14, 740000138]
...
476 477 [4, 6, 8, 10, 12, 14, 740000138]
477 478 [4, 6, 8, 10, 12, 14, 740000138]
478 479 [4, 6, 8, 10, 12, 14, 740000138]
479 480 [4, 6, 8, 10, 12, 14, 740000138]
480 481 [4, 6, 8, 10, 12, 14, 740000138]
[4, 6, 8, 10, 12, 14, 740000138]
"""

# In:
for e in all_exceptions:
    if e >= 40:
        print(e, is_good(e))

# Out:
"""
740000138 (True, [21, 235, 247, 391])
"""
