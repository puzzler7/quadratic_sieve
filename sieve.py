#!/usr/bin/env python3

import sys
from math import log10, log2, log, floor, ceil, gcd, exp
import time
from itertools import product
from sparse import *

# library stuff
import numpy as np # pip install numpy
from scipy.linalg import lu # pip install scipy

DEBUG = 0

primeMemo = {}
primeSets = {}
primeCounts = {}

primeNumbers = []
maxPrime = 3

def primefacMemo(n, b=99999999999999999999999):
    global primeMemo
    if n in primeMemo.keys():
        for i in primeMemo[n]:
            yield i
    else:
        primeMemo[n] = []
        primeSets[n] = set([])
        primeCounts[n] = {}
        for i in primefac(n, b=b):
            primeMemo[n].append(i)
            primeSets[n].add(i)
            if i in primeCounts[n]:
                primeCounts[n][i] += 1
            else:
                primeCounts[n][i] = 1
            yield i

def primefac(n, start=3, b=999999999999999999999999, idx=len(primeNumbers)):
    if n == 1:
        return
    if n <= 0:
        raise ValueError
    if n%2 == 0:
        yield 2
        for p in primefac(n//2, 2, b):
            yield p
        return

    if idx < len(primeNumbers):
        for i in range(idx, len(primeNumbers)):
            prime = primeNumbers[i]
            if n%prime == 0:
                yield prime
                for p in primefac(n//prime, start=prime, b=b, idx=i):
                    yield p
                return
            if prime > b:
                yield b*2

    i = max(start, maxPrime)
    while True:
        if i > intSqrt(n):
            yield n
            return
        if i > b: #stop when reaching bsmooth limit
            yield b*2
            return
        if n%i == 0:
            yield i
            if n//i == 1:
                return
            for p in primefac(n//i, i, b):
                yield p
            return
        i += 2


def bsmooth(n, b): # can this primefac be memoized without bad things happening?
    # print(n)
    for i in primefac(n, b=b):
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
    bot = 1
    top = n//2
    last = -1
    while True:
        val = (bot+top)//2
        sq = val**2
        if sq == n:
            sqrts[n] = val
            return val
        if sq < n:
            bot = val
        if sq > n:
            top = val
        if last == val:
            if (val+1)**2 == n:
                sqrts[n] = val+1
                return val+1
            sqrts[n] = val
            return val
        last = val

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

def quadsieve(n):
    global primeNumbers
    primetime = time.time()
    print("reading primes from file")
    # primeNumbers = eval(open("primes.txt", "r").read())
    print("loaded %d primes"%len(primeNumbers))
    timeElapsed(primetime)
    fac = 1
    while True:
        print("Trying fac = %f"%fac)
        p, q = quadsieveloop(n, fac)
        if p != 1 and q != 1:
            return p, q
        fac += .5

def timeElapsed(t):
    diff = time.time()-t
    h = int(diff)//3600
    m = int(diff-3600*h)//60
    s = round(diff%60, 3)
    print("\tTime elapsed: %dh %dm %0.3fs\n"%(h,m,s))
    return time.time()

def quadsieveloop(n, fac):
    timevar = time.time()
    looptime = timevar
    for i in range(2, floor(log2(n))):
        if n % i == 0:
            return i, n//i
    eps = .01
    x = 2*int(n**(.5+eps))
    o = .1
    b = int(exp((2**.5+o)*((log(x)*log(log(x)))**.5)))
    print("b =", b)
    t = int(fac*int(pi(b)))
    print("pi(b):", t)
    count = 0
    testnum = intSqrt(n)
    if testnum**2 == n:
        return testnum, testnum
    testnum += 1
    smoothnums = []
    smallprimes = set([])
    cong = {}
    print("finding %d bsmooth numbers"%t)
    while count <= t:
        sq = pow(testnum, 2, n)
        if bsmooth(sq, b):
            smoothnums.append(sq)
            cong[sq] = testnum
            smallprimes = smallprimes.union(set(primefacMemo(sq)))
            count += 1
        testnum += 1

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
            mat.add(i, primeIndices[f], primeCounts[num][f]%2)
            iters += 1


    print("iters run:", iters)
    print("elements:", numnums*rowlen)
    print("done factoring loop")
    timevar = timeElapsed(timevar)

    print("transposing and adding identity")
    height = mat.cols
    # mat = mat.transpose()
    # mat = mat.add_identity()
    # print(smallprimes)
    # print(cong)

    # mat = np.append(mat, np.identity(len(mat[0])), axis=0)

    print("done prepping for reduction")
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