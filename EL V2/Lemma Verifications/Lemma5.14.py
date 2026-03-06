from sympy import nextprime, isprime, prevprime, factorint, primerange

# Step 1: Precompute a set of pairwise almost-coprime values. Here, A and B are "almost-coprime" if gcd(A,B) = 2^k for any k \geq 0. 

def squarefree_odds(a):
    for p, e in factorint(a).items():
        if p > 2 and e > 1:
            return False
    return True

def coprime_weaker(a,b):
    fa = factorint(a)
    fb = factorint(b)
    for p in fa:
        if p in fb and p > 2:
            return False
    return True

def three_progressions():
    options = set()
    p = 3
    while p < 5000:
        print(p, len(options))
        q = nextprime(p)
        while q < p + 500:
            if squarefree_odds(p+q) == True and False not in [coprime_weaker(p+q,j[-1]) for j in options]:
                options.add((p,q,p+q))
            q = nextprime(q)
        p = nextprime(p)
    return options

tp = three_progressions()
print("Step 1 complete.")

# Step 2: The main computation.. 

def squarefree(a):
    for p, e in factorint(a).items():
        if e > 1:
            return False
    return True

def check(w):
    reps = []
    pr = 2
    while pr < w:
        if squarefree(w - pr) == True and False not in [coprime_weaker(w - pr,j) for j in reps]:
            reps.append(w-pr)
            if len(reps) == 3:
                return True
        pr = nextprime(pr)
    return False

def proc(fromme, upto, preloads):
    representations = {}
    prime = prevprime(fromme - 1000)
    while prime < upto:
        for a,b,c in preloads:
            e = 0
            m = prime + c*(2**e)
            while m < upto:
                if m not in representations:
                    representations[m] = 1
                else:
                    representations[m] += 1
                e += 1
                m = prime + c*(2**e)
        prime = nextprime(prime + 500) # Toggle 500 to increase/decrease speed (500 good for width of 10**6 following trial and error..)

    exceptions = []
    w = fromme
    while w <= upto:
        if w not in representations or representations[w] < 3:
            if check(w) == False:
                exceptions.append(w)
        w += 2
    return exceptions

# Run the following code..
ell = 4810
while ell*10**6 <= 8*(10**9):
    print(ell, proc((ell - 1)*(10**6) + 1, ell*(10**6),tp))
    ell += 1

"""
Output:
Step 1 complete.

1 []
2 []
3 []
4 []
5 []
6 []
7 []
8 []
9 []
10 []
11 []
12 []
13 []
14 []
15 []
16 []
17 []
18 []
19 []
20 []
21 []
22 []
23 []
24 []
25 []
26 []
...
4806 []
4807 []
4808 []
4809 []
4810 []

Summary: No exceptions detected.
"""