#!/bin/python3

from triangulationDev import getRandomRootedTwoLeafTree, rootedTwoLeafTreeFromBit, checkPrefixCondition
from MapPermutation import MapPermutation
from LabelledMap import LabelledMap
from MapGenerator import MapGenerator
from PrimitiveMutableLabelledMap import PrimitiveMutableLabelledMap

mapGenerator = MapGenerator()


def test_0():
    for i in range(40, 61):
        for _ in range(1000):
            try:
                myMap = getRandomRootedTwoLeafTree(i)
            except Exception as e:
                print("Size", i)
                print("Error", e)
                exit()


tricky = [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]


def test_1():
    alpha = [(1, 2), (3, 4), (5, 6), (7, 8),
             (9, 10), (11, 12), (13, 14), (15, 16)]
    sigma = [(1,), (2, 3, 5), (4,), (6, 7, 9, 11),
             (8,), (10,), (12, 13, 15), (14,), (16,)]
    tree = PrimitiveMutableLabelledMap(alpha=MapPermutation(
        alpha), sigma=MapPermutation(sigma))
    tree.pretty_print()
    mapGenerator.randomTreeToTriangulation(tree)


def test_2():
    alpha = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10), (11, 12), (13, 14),
             (15, 16), (17, 18), (19, 20), (21, 22), (23, 24), (25, 26), (27, 28)]
    sigma = [(1,), (2, 3, 9), (4, 5, 7), (6,), (8,), (10, 11, 13, 15, 27), (12,),
             (14,), (16, 17, 23, 25), (18, 19, 21), (20,), (22,), (24,), (26,), (28,)]

    tree = PrimitiveMutableLabelledMap(alpha=MapPermutation(
        alpha), sigma=MapPermutation(sigma))
    tree.pretty_print()
    mapGenerator.randomTreeToTriangulation(tree)


def test_3():
    tree = PrimitiveMutableLabelledMap(lmap=rootedTwoLeafTreeFromBit(tricky))
    out = mapGenerator.randomTreeToTriangulation(tree)
    out.show()


def test_4():
    N = 100000
    tree = PrimitiveMutableLabelledMap(lmap=getRandomRootedTwoLeafTree(N))
    tree.pretty_print()
    out = mapGenerator.randomTreeToTriangulation(tree)


def test_5():
    rootedTwoLeafTreeFromBit(tricky)


def test_6():
    for i in range(1, 2000):
        print(i)
        out = mapGenerator.getRandomTriangulation(i)
        if not out.isTriangulation():
            print("Oh no", i)
            exit()


def test_7():

    out = mapGenerator.getRandomTriangulation(1)
    out.show()


if __name__ == "__main__":
    test_6()
