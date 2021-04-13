#!/usr/bin/env python3

class SparseMatrix:
    def __init__(self, r, c):
        self.pointsx = {}
        self.pointsy = {}
        self.rows = r
        self.cols = c

    def add(self, x, y, val):
        if val == 0:
            return
        if x  not in self.pointsx:
            self.pointsx[x] = {}
        if y not in self.pointsy:
            self.pointsy[y] = {}
        self.pointsx[x][y] = val
        self.pointsy[y][x] = val

    def get(self, x, y):
        try:
            return self.pointsx[x][y]
        except KeyError as e:
            return 0

    def transpose(self):
        ret = SparseMatrix(self.cols, self.rows)
        ret.pointsx = self.pointsy
        ret.pointsy = self.pointsx
        return ret

    def to_array(self):
        import numpy as np
        ret = np.zeros((self.rows, self.cols))
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
                for k in range(self.cols):
                    v += self.get(i, k)*other.get(k, j)
                ret.add(i, j, v)
        return ret

    def add_mat(self, other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError
        ret = SparseMatrix(self.rows, self.cols)
        for i in range(self.rows):
            for j in range(self.cols):
                ret.add(i, j, self.get(i, j)+other.get(i, j))
        return ret

    def sub_mat(self, other):
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError
        ret = SparseMatrix(self.rows, self.cols)
        for i in range(self.rows):
            for j in range(self.cols):
                ret.add(i, j, self.get(i, j)-other.get(i, j))
        return ret

    def scale(self, r):
        ret = SparseMatrix(self.rows, self.cols)
        for x in self.pointsx:
            for y in self.pointsx[x]:
                ret.add(x, y, r*self.get(x,y))
        return ret

    def eq0(self):
        for x in self.pointsx:
            for y in self.pointsx[x]:
                if abs(self.get(x, y)) > 1e-20:
                    return False
        return True

    def __str__(self):
        return str(self.pointsx) 

if __name__ == "__main__":
    a = SparseMatrix(3, 3)
    b = SparseMatrix(3, 3)

    for i in range(3):
        a.add(i,i,1)

    for i in range(9):
        b.add(i//3, i%3, i+1)
    print(a.multiply(b).to_array())
    print(b.multiply(a).to_array())

