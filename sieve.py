#!/usr/bin/env python3

import sys
from math import log10, log2, log, floor, ceil, gcd, exp, isqrt
import time
import random
from sparse import *

DEBUG = 0

primeMemo = {}
primeSets = {}
primeCounts = {}

primeNumbers = [3]
maxPrime = 3

def primefacMemo(n, start=3, b=999999999999999999999999, idx=len(primeNumbers)):
    global primeMemo
    if n in primeMemo.keys():
        for i in primeMemo[n]:
            yield i
    else:
        primeMemo[n] = []
        primeSets[n] = set([])
        primeCounts[n] = {}
        for i in primefac(n, start=start, b=b, idx=idx):
            primeMemo[n].append(i)
            primeSets[n].add(i)
            if i in primeCounts[n]:
                primeCounts[n][i] += 1
            else:
                primeCounts[n][i] = 1
            yield i

def pollardRho(n):
    print(n)
    g = lambda x: (x**2+1)%n
    x = random.randint(2, n)
    y = x
    d = 1

    loops = 0

    while d == 1:
        if loops >= 42000:
            return 1
        x = g(x)
        y = g(g(y))
        d = gcd(abs(x-y))
        loops += 1
    return d

def pollardRho_pm1(n, b = 10):
    b = int(b)
    if b > primeNumbers[-1]:
        return 1
    m = 1
    i = 0
    q = primeNumbers[i]
    while q < b:
        print(m)
        m *= q**int(log(b, q))
        q = primeNumbers[i]
        i += 1
    g = gcd(pow(2, m, n)-1, n)
    if g == 1:
        return pollardRho_pm1(n, b*10)
    if g == n:
        return pollardRho_pm1(n, b//2)
    return g


def primefac(n, start=3, b=999999999999999999999999, idx=len(primeNumbers)):
    if n == 1:
        return
    if n <= 0:
        raise ValueError
    if n%2 == 0:
        yield 2
        for p in primefacMemo(n//2, 2, b):
            yield p
        return
    sq = intSqrt(n)

    if idx < len(primeNumbers):
        for i in range(idx, len(primeNumbers)):
            prime = primeNumbers[i]
            if n%prime == 0:
                yield prime
                for p in primefacMemo(n//prime, start=prime, b=b, idx=i):
                    yield p
                return
            if prime > b:
                yield b*2

    # d = pollardRho_pm1(n)
    # if d != n and d != 0:
    #     yield d
    #     # print("solved with pollardRho")
    #     for p in primefacMemo(n//d, start=maxPrime, b=b, idx=len(primeNumbers)):
    #         yield p

    i = max(start, maxPrime)
    while True:
        if i > sq:
            yield n
            return
        if i > b: #stop when reaching bsmooth limit
            yield b*2
            return
        if n%i == 0:
            yield i
            if n//i == 1:
                return
            for p in primefacMemo(n//i, i, b):
                yield p
            return
        i += 2


def bsmooth(n, b): # can this primefac be memoized without bad things happening?
    # print(n)
    for i in primefacMemo(n, idx=0, b=b):
        # print(i)
        if i > b:
            # print(i)
            return False
    return True

def pi(b): #x/(ln(x)-1)
    return ceil(b/(log(b, 2))) # approximation
    ret = 0
    for i in range(2, b):
        ret += millerPrimeTest(i)*1 
    return ret

def pr(mat):
    if DEBUG:
        for row in mat:
            print(list([int(i)%2 for i in row]))
        print()

sqrts = {}

def intSqrt(n):
    global sqrts
    if n in sqrts:
        return sqrts[n]
    ret = isqrt(n)
    sqrts[n] = ret
    return ret

def millerPrimeTest(n): # Uses Miller-Rabin three times to test primality.
    bases = [2, 3, 5] # Bases to use Miller-Rabin with
    if not isinstance(n, int): # Handle edge cases ...
        return False
    if n in [2, 3, 5, 7]:
        return True
    if n <= 10:
        return False
    if n % 2 == 0:
        return False
    m = n - 1
    k = 0
    while m % 2 == 0:
        m = m // 2
        k = k + 1
    # Now note that n - 1 = 2^k * m, with m odd
    for a in bases:
        x = pow(a, m, n)  # Sets x to a^m mod n
        if x in [n - 1, 1]:
            return True
        iterations = 0
        while iterations <= k - 1:
            x = (x * x) % n
            if x == 1:
                return False
            if x == n - 1:
                break
            iterations = iterations + 1
        if x != n - 1:
            return False
    return True  # Passed three rounds of Miller-Rabin

def legendre(a, p):
    ret = pow(a, (p-1)//2, p)
    return -1 if ret == p-1 else ret

def sqrtModP(n, p):
    q = p-1
    s = 0
    while q%2 == 0:
        q //= 2
        s += 1
    z = 1
    # quadres = []
    # for i in range(1, ((p-1)//2+1)):
    #     quadres.append(pow(i, 2, p))
    # while z in quadres:
    while legendre(z, p) != -1:
        z += 1
    assert z != p
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q+1)//2, p)
    while True:
        if t == 0:
            return 0, 0
        if t == 1:
            return r, p-r
        i = 0
        temp = t
        while temp != 1:
            temp = pow(temp, 2, p)
            i += 1
        if m == i:
            return 0,0
        b = pow(c, pow(2, m-i-1), p)
        m = i
        c = pow(b, 2, p)
        t = (t*c) % p
        r = (r*b) % p


def quadsieve(n):
    global primeNumbers
    primetime = time.time()
    print("reading primes from file")
    primeNumbers = eval(open("primes.txt", "r").read())
    maxPrime = primeNumbers[-1]
    print("loaded %d primes"%len(primeNumbers))
    timeElapsed(primetime)
    fac = 1
    while True:
        print("Trying fac = %f"%fac)
        p, q = quadsieveloop(n, fac)
        if p != 1 and q != 1:
            return p, q
        fac += .05

def timeElapsed(t):
    diff = time.time()-t
    h = int(diff)//3600
    m = int(diff-3600*h)//60
    s = round(diff%60, 3)
    print("\tTime elapsed: %dh %dm %0.3fs\n"%(h,m,s))
    return time.time()

smoothnums = []
cong = {}

def quadsieveloop(n, fac):
    global smoothnums, cong
    timevar = time.time()
    looptime = timevar
    for i in range(2, floor(log2(n))):
        if n % i == 0:
            return i, n//i
    eps = 0
    x = 2*int(n**(.5+eps))
    o = .15
    b = int(exp((2**(-.5)+o)*((log(x)*log(log(x)))**.5)))
    print("b =", b)
    t = int(fac*int(pi(b)))
    print("pi(b):", t)
    count = 0
    testnum = intSqrt(n)
    if testnum**2 == n:
        return testnum, testnum
    testnum += 1
    
    print("finding %d bsmooth numbers"%t)
    # thousandtime = time.time()
    # while count <= t:
    #     sq = pow(testnum, 2, n)
    #     if bsmooth(sq, b):
    #         smoothnums.append(sq)
    #         cong[sq] = testnum
    #         smallprimes = smallprimes.union(set(primefacMemo(sq)))
    #         count += 1
    #         if count %100 == 0:
    #             print(count, "bsmooth nums found")
    #             thousandtime = timeElapsed(thousandtime)
    #     testnum += 1
    # print(sq)

    sq = intSqrt(n)+1
    y = lambda x: (x+sq)**2 - n
    vlen = int(t*t*intSqrt(t)*1.2)
    print("v len is", vlen)

    fb = [2]
    idx = 1
    while len(fb) < t:
        p = primeNumbers[idx]
        if legendre(n, p) == 1:
            fb.append(p)
        idx += 1
    # print("factor base", fb)

    V = [y(x) for x in range(vlen)]
    if V[0]%2 == 0:
        for i in range(0, vlen, 2):
            assert V[i] %2 == 0
            V[i] //= 2
    else:
        for i in range(1, vlen, 2):
            assert V[i] %2 == 0
            V[i] //= 2
    for i in range(1, t):
        p = fb[i]
        # print("running sieve mod", p)
        r1, r2 = sqrtModP(n, p)
        if r1 == 0:
            print("sqrt is 0")
            continue
        # print("found sqrt of %d mod %d"%(n, p))
        # print(r1, r2)
        if pow(r1, 2, p)!=n%p or pow(r2, 2, p)!=n%p:
            print("not actually root?", r1, r2)
            continue
        r1 = (r1 - sq) % p
        r2 = (r2 - sq) % p

        for idx in range(r1, vlen, p):
            try:
                assert V[idx] %p == 0
                V[idx] //= p
            except IndexError as e:
                print("bad index")
        if p != 2:
            for idx in range(r2, vlen, p):
                assert V[idx] %p == 0
                V[idx] //= p
        if i % 100 == 0:
            print("through %d bases"%i)
    smallprimes = fb
    for i, test in enumerate(V):
        if test == 1 or test == 0:
            num = y(i)
            if num in smoothnums:
                continue
            smoothnums.append(num)
            cong[num] = sq+i
            list(primefacMemo(num)) # calculating for later
            # smallprimes = smallprimes.union(set(primefacMemo(num)))
    # print(smoothnums)


    print("found bsmooth numbers")
    timevar = timeElapsed(timevar)

    print("indexing small primes")

    smallprimes = sorted(list(smallprimes))
    primeIndices = {}
    for i, p in enumerate(smallprimes):
        primeIndices[p] = i
    rowlen = len(smallprimes)
    numnums = len(smoothnums)
    print((numnums,rowlen))

    print("indexed small primes")
    timevar = timeElapsed(timevar)

    print("factoring numbers")

    mat = SparseMatrix(numnums, rowlen)
    iters = 0
    for i, num in enumerate(smoothnums):
        # we can probably do row generation in the last loop
        # but that might mean we'll be doing some unnecessary computation...

        for f in primeSets[num]:
            try:
                mat.add(i, primeIndices[f], primeCounts[num][f]%2)
                iters += 1
            except KeyError as e:
                print("bad prime", f)
                continue


    print("iters run:", iters)
    print("elements:", numnums*rowlen)
    print("done factoring loop")
    timevar = timeElapsed(timevar)

    print("doing gaussian reduction")
    pr(mat.to_array())
    congruences = mat.get_congruences()
    pr(mat.to_array())
    # print(congruences)

    timevar = timeElapsed(timevar)

    print("evaluating congruences")

    smooths = []
    congs = []
    done = 0
    factor = 1
    for row in congruences:
        primesToUse = [smoothnums[i] for i in row]

        # print(primesToUse)

        smoothsq = 1
        congsq = 1
        for i in primesToUse:
            smoothsq *= i
            congsq *= cong[i]

        # print(smoothsq)
        # print(congsq)
        assert pow(smoothsq, 1, n) == pow(congsq, 2, n)

        factor = gcd(abs(intSqrt(smoothsq)-congsq), n)
        if factor != 1 and factor != n:
            done = 1
            break

    print("done evaluating congruences")
    timevar = timeElapsed(timevar)

    print("loop finished")
    timeElapsed(looptime)

    return factor, n//factor

if __name__ == '__main__':
    args = sys.argv
    try:
        if len(args) == 1:
            n = int(input("Enter a product of two primes: "))
            ret = quadsieve(n)
        else:
            n = int(args[1])
            ret = quadsieve(n)
    except ValueError as e:
        print("Bad int!")
        raise e

    # assert ret[0]*ret[1] == n
    print("The factors are %d and %d."%ret)