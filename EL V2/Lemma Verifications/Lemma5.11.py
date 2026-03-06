from sympy import prevprime, primerange, nextprime, factorint, isprime

def verify(to_me):
    #generate square frees up to N:
    sq_frees = set([1])
    for k in range(2,10**5):
        if max(factorint(k).values()) == 1:
            sq_frees.add(k)
    print("Square-frees generated.")

    #populate a set of verified integers..
    verified = set()
    primes = [q for q in primerange(2,10)]
    for q in primes:
        for s in sq_frees:
            if q + s not in verified:
                verified.add(q + s)
    print(len(verified))

    np = nextprime(10)
    checkpoint = 10**7
    while np < to_me:
        for s in sq_frees:
            if np + s not in verified:
                verified.add(np + s)
        lv = len(verified)
        if lv >= checkpoint:
            print(lv)
            checkpoint += 10**7
        np = nextprime(np + 10**4)
    print(lv)

    #re-generate square frees up to N:
    sq_frees = [1]
    for k in range(2,10**4):
        if max(factorint(k).values()) == 1:
            sq_frees.append(k)
    print("Square-frees re-generated.")

	#create a list of exceptional integers using brute force on integers not in verified..
    exceptions = []
    for m in range(3, to_me + 1):
        if m not in verified:
            found = False
            ind = 0
            while found == False:
                if isprime(m - sq_frees[ind]) == True:
                    found = True
                else:
                    ind += 1
            if found == False:
                exceptions.append(m)
    print(exceptions)

def verify_ft(from_me, to_me):
    #generate square frees up to N:
    sq_frees = set([])
    for k in range(2,10**5):
        if max(factorint(k).values()) == 1:
            sq_frees.add(k)

    verified = set()
    np = prevprime(from_me)
    while np < to_me:
        for s in sq_frees:
            if np + s not in verified:
                verified.add(np + s)
        np = nextprime(np + 10**4)

    #re-generate square frees up to N:
    sq_frees = []
    for k in range(2,10**4):
        if max(factorint(k).values()) == 1:
            sq_frees.append(k)

    exceptions = []
    for m in range(from_me, to_me + 1):
        if m not in verified:
            #print(m)
            found = False
            ind = 0
            while found == False:
                if isprime(m - sq_frees[ind]) == True:
                    found = True
                else:
                    ind += 1
            if found == False:
                exceptions.append(m)
    return exceptions

for l in range(481, 800):
    print(l, l+1, verify_ft(l*10**7, (l+1)*10**7))

"""
Output: 

1 2 []
2 3 []
3 4 []
4 5 []
5 6 []
6 7 []
7 8 []
8 9 []
9 10 []
10 11 []
11 12 []
12 13 []
13 14 []
14 15 []
15 16 []
16 17 []
17 18 []
18 19 []
19 20 []
20 21 []
21 22 []
22 23 []
23 24 []
24 25 []
25 26 []
26 27 []
27 28 []
28 29 []
29 30 []
30 31 []
31 32 []
32 33 []
33 34 []
34 35 []
35 36 []
36 37 []
37 38 []
38 39 []
39 40 []
40 41 []
41 42 []
42 43 []
43 44 []
44 45 []
45 46 []
46 47 []
47 48 []
48 49 []
49 50 []
50 51 []
51 52 []
52 53 []
53 54 []
54 55 []
55 56 []
56 57 []
...
471 472 []
472 473 []
473 474 []
474 475 []
475 476 []
476 477 []
477 478 []
478 479 []
479 480 []
480 481 []
"""
