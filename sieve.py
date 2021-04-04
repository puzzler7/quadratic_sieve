#!/usr/bin/env python3

import sys
from math import log2, log, floor, ceil, gcd

from primefac import primefac # pip install git+git://github.com/elliptic-shiho/primefac-fork@master
from gmpy import mpz # pip install gmpy
from Crypto.Util.number import isPrime # pip install pycryptodome
from pyfinite.genericmatrix import GenericMatrix # pip install pyfinite

def bsmooth(n, b):
    for i in primefac(n): # FIXME need to write own primefac
        if i > b:
            print(i)
    return True

def pi(b): #1.3x/ln(x)
    return ceil(2*b/log(b)) # approximation
    ret = 0
    for i in range(2, b):
        ret += isPrime(i)*1 # FIXME improve, need own isPrime?
    return ret

def exp(base, exponent, modulus): # Fast exponentiation function using successive squaring
    toMultiply = base # Number we should multiply by next: Successively square this
    ans = 1 # Return value
    z = exponent # Keeps track of binary representation of exponent
    while (z > 0):
        digit = z % 2 # Last digit in binary representation of exponent
        z = z // 2 # Remove last digit from z
        if (digit == 1):
            ans = (ans * toMultiply) % modulus
        toMultiply = (toMultiply * toMultiply) % modulus # Square part of successive square
    return ans

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
        x = exp(a, m, n)  # Sets x to a^m mod n
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
    for i in range(2, floor(log2(n))):
        if n % i == 0:
            return i, n//i
    b = 10000
    t = pi(b)
    print("pi(b):", t)
    count = 0
    testnum = int(mpz(n).root(2)[0])
    smoothnums = []
    smallprimes = set([])
    while count <= t:
        if bsmooth(testnum, b):
            smoothnums.append(testnum)
            smallprimes = smallprimes.union(set(primefac(testnum)))
            count += 1
        testnum += 1

    smallprimes = sorted(list(smallprimes))
    rowlen = len(smallprimes)
    numnums = len(smoothnums)
    print((numnums,rowlen))

    XOR = lambda x,y: x^y
    AND = lambda x,y: x&y
    DIV = lambda x,y: x
    mat = GenericMatrix(size=(numnums,rowlen),
        zeroElement=0,identityElement=1,add=XOR,mul=AND,sub=XOR,div=DIV)
    rowcount = 0
    for i in smoothnums:
        # print(i)
        factor = list(primefac(i))
        newrow = []
        for p in smallprimes:
            newrow.append(factor.count(p)%2)
        mat.SetRow(rowcount, newrow)
        rowcount += 1
    b = [0]*numnums
    x = mat.LowerGaussianElim()
    print(x)


    return n, n

if __name__ == '__main__':
    args = sys.argv
    try:
        if len(args) == 1:
            n = int(input("Enter a product of two primes: "))
            ret = quadsieve(n)
        else:
            ret = quadsieve(int(args[1]))
    except ValueError as e:
        print("Bad int!")
        print(e)

    assert ret[0]*ret[1] == n
    print("The factors are %d and %d."%ret)