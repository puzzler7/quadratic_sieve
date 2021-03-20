#!/usr/bin/env python3

import sys
from math import log2, floor, gcd

from primefac import primefac
from Crypto.Util.number import isPrime

def bsmooth(n, b):
    for i in primefac(n): # FIXME need to write own primefac
        if i > b:
            return False
    return True

def pi(b):
    ret = 0
    for i in range(2, b):
        ret += isPrime(i)*1 # FIXME improve, need own isPrime?
    return ret

def quadsieve(n):
    for i in range(2, floor(log2(n))):
        if n % i == 0:
            return i, n//i
    b = 100000
    t = pi(b)
    for i in range(1, t+2):


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