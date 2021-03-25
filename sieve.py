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