from LabelledMap import *
from MutableLabelledMap import *
from DynamicShow import *
from MapGenerator import *

from sage.all import *

import igraph as ig
import networkx as nx
import matplotlib.pyplot as plt
import timeit


def m1():
    return MutableLabelledMap(Permutation([(9, 1, 11), (6, 4, 8), (5, 12, 10, 7)]), Permutation(
        [(1, 2), (11, 12), (9, 10), (7, 8), (5, 6), (3, 4)]))


def m2():
    return MutableLabelledMap(Permutation([(9, 1, 8), (3, 2), (4, 10, 5), (6, 7)]), Permutation(
        [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10)]))


def m22():
    return MutableLabelledMap(Permutation([(6, 9, 13), (10, 4, 8), (1, 2, 3), (7, 11, 5), (14,)]), Permutation(
        [(1, 5), (2, 4), (3, 6), (7, 8), (9, 10), (11, 12), (13, 14)]))


def m3():
    return MutableLabelledMap(Permutation(
        [(3, 6, 5, 1), (2, 4)]), Permutation([(1, 2), (3, 4), (5, 6)]))


def star(n):
    return MutableLabelledMap(
        adj=[list(range(2, n + 2))] + [(1,) for i in range(n)])

# m = MutableLabelledMap(adj = ([2,4],[1,3],[2,4],[1,3]))
# m = MutableLabelledMap(sigma = Permutation([(9,1,5,3,6),(2,10),(4,7,8)]), alpha = Permutation([(1,2),(3,4),(5,6),(7,8),(9,10)]))
# m = MutableLabelledMap(adj = [(2,4),(1,4,3),(2,4),(3,2,1)])
# m = MutableLabelledMap(adj = [(2,3,4,5,6,7,8,9), (1,3),(2,1),(1,),(1,),(1,),(1,),(1,),(1,)])
# m = MutableLabelledMap(adj = [(2,3,4,), (1,), (1,), (1,)])
# m = star(15)
# m = MutableLabelledMap(adj = [(4,3,2),(3,4,1),(2,1,5),(2,1),(3,6),(5,7),(6,)])


# m = MutableLabelledMap(sigma = Permutation([(1,3,5,2), (4,6)]), alpha = Permutation([(1,2),(3,4),(5,6)]))
# m = MutableLabelledMap(sigma = Permutation([(1,3),(2,4)]), alpha = Permutation([(1,2),(3,4)]))
mm = star(6)

# plt.ion()

seed = 25

gen = MapGenerator()
m = gen.getRandomPlanarMap(15, seed)
c = gen.cube()

l = gen.getRandomPlanarQuadrangulation(30, seed)

d = DynamicShow(c.derivedMap())
d.start(True)

plt.show()
