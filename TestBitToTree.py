#!/bin/python3

from triangulationDev import bitToTree
from triangulationDev import bitToTreeDebug
from LabelledMap import LabelledMap


def test_1():
    b = [1, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    sigma, alpha = bitToTreeDebug(b)
    sigma.pretty_print()
    alpha.pretty_print()
    try:
        myMap = LabelledMap(sigma, alpha)
        myMap.show()
    except Exception as e:
        print("="*50)
        print(e)


def test_0():
    for i in range(40, 61):
        for _ in range(1000):
            try:
                myMap = bitToTree(i)
            except Exception as e:
                print("Size", i)
                print("Error", e)
                exit()


if __name__ == "__main__":
    test_0()
