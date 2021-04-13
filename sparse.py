#!/usr/bin/env python3

class SparseMatrix:
    def __init__(self, r, c):
        self.pointsx = {}
        self.pointsy = {}
        self.rows = r
        self.cols = c

    def add(self, x, y, val):
        if x >= self.rows:
            self.rows = x+1
        if y >= self.cols:
            self.cols = y+1
        if x  not in self.pointsx:
            self.pointsx[x] = {}
        if y not in self.pointsy:
            self.pointsy[y] = {}
        self.pointsx[x][y] = val
        self.pointsy[y][x] = val
        # if val == 0:
        #     del self.pointsx[x][y]
        #     del self.pointsy[y][x]

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
                if i not in self.pointsx or j not in other.pointsy:
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
                if abs(self.get(x, y)) > 1e-10:
                    return False
        return True

    def add_col(self, fr, to):
        for x in self.pointsy[fr]:
            self.add(x, to, (self.get(x, fr)+self.get(x, to))%2)

    def find_pivot(self, j):
        if j not in self.pointsy:
            return -1
        for x in self.pointsy[j]:
            if self.get(x, j) and x not in self.marked:
                return x
        return -1

    def reduce(self):
        if self.rows < self.cols:
            print("Matrix size (%d, %d) bad for reducing!"%(self.rows, self.cols))
        print(self)
        self.marked = {}
        self.marks = {}
        for j in range(self.cols):
            i = self.find_pivot(j)
            if i == -1:
                continue
            self.marked[i] = 1 #dictionary for faster lookup
            self.marks[j] = i
            for k in range(self.cols):
                if j != k and self.get(i, k) == 1:
                    self.add_col(j, k)

    def get_congruences(self):
        self.reduce()
        ret = []
        for i in range(self.rows):
            cong = [i]
            if i in self.marked or i not in self.pointsx:
                continue
            for y in self.pointsx[i]:
                if self.get(i, y):
                    cong.append(self.marks[y])
            if len(cong) > 1:
                ret.append(cong)
        return ret

    def __str__(self):
        return str(self.to_array())

if __name__ == "__main__":
    a = SparseMatrix(3, 3)
    b = SparseMatrix(3, 3)

    for i in range(3):
        a.add(i,i,1)

    for i in range(9):
        b.add(i//3, i%3, i+1)
    # print(a.multiply(b).to_array())
    # print(b.multiply(a).to_array())
    # print(a.add_mat(b).to_array())
    # print(b.add_mat(a).to_array())
    # print(a.sub_mat(b).to_array())
    # print(b.sub_mat(a).to_array())
    # print(b.multiply(b).to_array())
    # print(b.scale(-5).to_array())
    # print(b.transpose().to_array())

    c = SparseMatrix(1, 3)
    d = SparseMatrix(3, 1)

    for i in range(3):
        c.add(0, i, i+1)
        d.add(i, 0, i+1)

    print(c.to_array())
    print(d.to_array())
    print(c.multiply(d).to_array())
    print(c.multiply(c.transpose()).to_array())
    print()
    print(d.multiply(c).to_array())
    print(d.multiply(d.transpose()).to_array())
    print(c.transpose().transpose().to_array())
    print(c.transpose().multiply(c).to_array())



