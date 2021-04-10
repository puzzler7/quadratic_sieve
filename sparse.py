#!/usr/bin/env python3

class Pair:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def swap(self):
        temp = self.x
        self.x = self.y
        self.y = temp

    def __str__(self):
        return "({}, {})".format(self.x, self.y)

class SparseMatrix:
    def __init__(self, r, c):
        self.points = {}
        self.rows = r
        self.cols = c

    def add(self, *args):
        if args[-1] == 0:
            return
        if len(args) == 2:
            self.points[args[0]] = args[1]
        elif len(args) == 3:
            self.points[Pair(args[0], args[1])] = args[2]
        else:
            raise NotImplemented

    def get(self, *args):
        if len(args) == 1:
            return self.points[args[0]]
        elif len(args) == 2:
            return self.points[Pair(args[0], args[1])]

    def transpose(self):
        for p in self.points.keys():
            p.swap()
        temp = self.rows
        self.rows = self.cols
        self.cols = temp
        return self

    def to_array(self):
        import numpy as np
        ret = np.zeros((self.rows, self.cols))
        print("{}, {}".format(self.rows, self.cols))
        for p in self.points.keys():
            ret[p.x][p.y] = self.points[p]
        return ret

    def add_identity(self):
        start = self.rows
        self.rows += self.cols
        i = 0
        while start+i < self.rows:
            self.add(start+i, i, 1)
            i += 1
        return self

    def __str__(self):
        return(str(self.points))
