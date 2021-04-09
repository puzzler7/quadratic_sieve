#!/usr/bin/env python3

import sys
from math import log2, log, floor, ceil, gcd
import time

from primefac import primefac # pip install git+git://github.com/elliptic-shiho/primefac-fork@master
from gmpy import mpz # pip install gmpy
import numpy as np # pip install numpy
from scipy.linalg import lu # pip install scipy

DEBUG = 0

primeMemo = {}

def primefacMemo(n):
    global primeMemo
    if n in primeMemo.keys():
        for i in primeMemo[n]:
            yield i
    else:
        primeMemo[n] = []
        for i in primefac(n):
            primeMemo[n].append(i)
            yield i

def bsmooth(n, b): # can this primefac be memoized? don't think so, but it should be fine here
    for i in primefac(n): # FIXME need to write own primefac
        if i > b:
            # print(i)
            return False
    return True

def pi(b): #1.3x/ln(x)
    return ceil(2*b/log(b)) # approximation
    ret = 0
    for i in range(2, b):
        ret += millerPrimeTest(i)*1 # FIXME improve, need own isPrime?
    return ret

def pr(mat):
    if DEBUG:
        for row in mat:
            print(list([int(i)%2 for i in row]))

def intSqrt(n): #uses lib fn rn, should change?
    return int(mpz(n).root(2)[0])

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
    fac = 1
    while True:
        print("Trying fac = %f"%fac)
        p, q = quadsieveloop(n, fac)
        if p != 1 and q != 1:
            return p, q
        if p == 1:
            fac += 1
        if q == 1:
            fac -= .5

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
    b = intSqrt(n)*fac
    t = int(pi(b))
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
        if bsmooth(testnum, b):
            smoothnums.append(testnum)
            cong[testnum] = testnum**2 - n
            smallprimes = smallprimes.union(set(primefacMemo(testnum)))
            count += 1
        testnum += 1

    print("found bsmooth numbers")
    timevar = timeElapsed(timevar)

    smallprimes = sorted(list(smallprimes))
    rowlen = len(smallprimes)
    numnums = len(smoothnums)
    print((numnums,rowlen))

    print("factoring numbers")

    mat = np.array([np.zeros(rowlen)])
    for i in smoothnums:
        # print(i)
        factor = list(primefacMemo(i)) # we can probably do row generation in the last loop
                                       # but that might mean we'll be doing some unnecessary computation...
        newrow = np.zeros(rowlen)

        for f in set(factor): 
            newrow[smallprimes.index(f)] = factor.count(f)%2

        mat = np.append(mat, [newrow], axis=0)

    mat = mat[1:]
    print("done factoring loop")
    timevar = timeElapsed(timevar)
    mat = mat.transpose()
    height = len(mat)
    # print(smallprimes)
    # print(cong)

    mat = np.append(mat, np.identity(len(mat[0])), axis=0)

    print("doing gaussian reduction")
    # pr(mat)
    p, l, u = lu(mat)
    print("solved")
    timevar = timeElapsed(timevar)
    pr(l)
    print()
    pr(u)

    mat = l.transpose()

    print("finding basis")
    basis = []
    for row in mat[::-1]:
        if sum([int(i)%2 for i in row[:height]]) == 0:
            basis.append(row[-numnums:])

    print() 
    pr(basis)
    print("basis found")
    timevar = timeElapsed(timevar)

    # basis is now basis for the kernel (I think?)
    # too tired to keep going, but the hard part should be done
    # just need to actually do the calculation (congruence) on the primes now

    primesToUse = basis[0]
    for row in basis[1:]:
        primesToUse += row

    primesToUse = [int(i)%2 for i in primesToUse]
    primesToUse = [smoothnums[i] for i, test in enumerate(primesToUse) if test]

    # print(primesToUse)

    smoothsq = 1
    congsq = 1
    for i in primesToUse:
        smoothsq *= i
        congsq *= cong[i]

    # print(smoothsq)
    # print(congsq)
    # assert pow(smoothsq, 2, n) == pow(congsq, 2, n)

    factor = gcd(abs(smoothsq-congsq), n)

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