#!/usr/bin/env python3

from Crypto.Util.number import *

if __name__ == '__main__':
    i = 3
    primes = [2]
    try:
        while True:
            if isPrime(i):
                primes.append(i)
            i += 2
            if (i-1)%10000 == 0:
                print(i)
    except KeyboardInterrupt as e:
        print("Ctrl-C caught")
        print('Num primes detected:', len(primes))
        f = open("primes.txt", 'w')
        f.write(str(primes))
        f.close()
        print("Primes written to primes.txt")

