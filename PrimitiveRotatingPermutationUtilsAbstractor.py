from PermutationUtilsAbstractor import *
from MapError import NotImplemented


class PrimitiveRotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    def __init__(self, rpermutation) -> None:
        self.rpermutation = rpermutation

    def numberInCycle(self, index):
        raise NotImplemented(self)

    def sameCycle(self, i, j):
        raise NotImplemented(self)

    def numberOfCycles(self):
        return self.rpermutation.numberOfCycles()

    def numberOfFixedPoint(self):
        return self.rpermutation.number_of_fixed_points()

    def checkTwoInTheSameCycle(self, listIndexes):
        raise NotImplemented(self)
