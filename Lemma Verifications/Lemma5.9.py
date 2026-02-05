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

verify(10**7)

"""
Output:

Square-frees generated.
99290
10004369
10071310
Square-frees re-generated.
[]
"""

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

for l in range(481,800):
    print(l, l+1, verify_ft(l*10**7, (l+1)*10**7))

"""
Output - NB the output from 1 to 481 is taken for granted from Lee-O'Clarey's work:

481 482 []
482 483 []
483 484 []
484 485 []
485 486 []
486 487 []
487 488 []
488 489 []
489 490 []
490 491 []
491 492 []
...
789 790 []
790 791 []
791 792 []
792 793 []
793 794 []
794 795 []
795 796 []
796 797 []
797 798 []
798 799 []
799 800 []
"""
