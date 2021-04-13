#!/usr/bin/env python3

from sparse import *

# based on https://link.springer.com/content/pdf/10.1007%2F3-540-49264-X_9.pdf

def c(i, j, w):
    num = w[j].transpose().multiply(A).multiply(A).multiply(w[i-1]).get(0,0)
    den = w[j].transpose().multiply(A).multiply(w[j]).get(0,0)
    return num/den

def std_lanczos(A, b):
    w = [b]
    w1 = A.multiply(b).sub_mat(w[0].scale(c(1,0,w)))
    w.append(w1)
    i = 2
    print(w[0].to_array())
    print(w[1].to_array())
    assert(len(w) == 2)
    while not w[-1].eq0():
        wnext = A.multiply(w[i-1]).sub_mat(w[i-1].scale(c(i, i-1, w)))
        wnext = wnext.sub_mat(w[i-2].scale(c(i,i-2,w)))
        w.append(wnext)
        i += 1
        print(wnext.to_array())
        # if len(w) == 10:
        #     raise RuntimeError
    x = SparseMatrix(A.cols, 1)
    test = SparseMatrix(A.cols, 1)
    for wv in w[:-1]:
        test = test.add_mat(wv)
        sc = wv.transpose().multiply(b).get(0,0)
        sc /= wv.transpose().multiply(A).multiply(wv).get(0,0)
        add = wv.scale(sc)
        x = x.add_mat(add)
        # print(x.to_array())
    print("sum", test.to_array())
    return x


if __name__ == "__main__":
    A = SparseMatrix(3, 3)
    b = SparseMatrix(3, 1)

    for i in range(3):
        A.add(i,i,1)
    A.add(0,2, 1)

    for i in range(3):
        b.add(i, 0, i+1)

    print("soln", std_lanczos(A, b).to_array())