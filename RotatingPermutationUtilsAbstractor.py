
from PermutationUtilsAbstractor import *


class RotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    def __init__(self, rpermutation) -> None:
        self.rpermutation = rpermutation

    def numberInCycle(self, index):
        return self.rpermutation.numberInCycle(index)

    def sameCycle(self, i, j):
        return self.rpermutation.sameCycle(i, j)

    def numberOfCycles(self):
        return self.rpermutation.numberOfCycles()

    def numberOfFixedPoint(self):
        return self.rpermutation.number_of_fixed_points()

    def checkTwoInTheSameCycle(self, listIndexes):
        return self.rpermutation.checkTwoInTheSameCycle(listIndexes)
