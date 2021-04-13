#!/usr/bin/env python3

class SparseMatrix:
    def __init__(self, r, c):
        self.pointsx = {}
        self.pointsy = {}
        self.rows = r
        self.cols = c

    def add(self, x, y, val):
        if x  not in self.pointsx:
            self.pointsx[x] = {}
        if y not in self.pointsy:
            self.pointsy[y] = {}
        self.pointsx[x][y] = val
        self.pointsy[y][x] = val

    def get(self, x, y):
        return self.pointsx[x][y]

    def transpose(self):
        ret = SparseMatrix(self.cols, self.rows)
        ret.pointsx = self.pointsy
        ret.pointsy = self.pointsx
        return ret

    def to_array(self):
        import numpy as np
        ret = np.zeros((self.rows, self.cols))
        print("{}, {}".format(self.rows, self.cols))
        for x in self.pointsx:
            for y in self.pointsx[x]:
                ret[x][y] = self.pointsx[x][y]
        return ret

    def add_identity(self):
        start = self.rows
        self.rows += self.cols
        i = 0
        while start+i < self.rows:
            self.add(start+i, i, 1)
            i += 1
        return self

    def multiply(self, other):
        if self.cols != other.rows:
            raise ValueError
        ret = SparseMatrix(self.rows, other.cols)
        for i in range(self.rows):
            for j in range(other.cols):
                if i not in self.pointsx or j not in self.pointsy:
                    continue
                v = 0
                for k in self.cols:
                    v += self.get(i, k)*other.get(k, j)
                ret.add(i, j, v)
        return ret

    def __str__(self):
        return(str(self.points))
