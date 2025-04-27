from PermutationUtilsAbstractor import *
from MapError import NotImplemented


class PrimitiveRotatingPermutationUtilsAbstractor(PermutationUtilsAbstractor):
    def __init__(self, rpermutation) -> None:
        self.rpermutation = rpermutation

    def numberInCycle(self, index):
        """
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor
        """
        raise NotImplemented(self)

    def sameCycle(self, i, j):
        """
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor
        """
        raise NotImplemented(self)

    def numberOfCycles(self):
        """
        Returns:
            The number of cycles of the permutation
        ---
        O(1)
        """
        return self.rpermutation.numberOfCycles()

    def numberOfFixedPoint(self):
        """
        Returns:
            The number of fixed point of the permutation
        ---
        O(1)
        """
        return self.rpermutation.number_of_fixed_points()

    def checkTwoInTheSameCycle(self, listIndexes):
        """ 
        Not implemented for PrimitiveRotatingPermutationUtilsAbstractor
        """
        raise NotImplemented(self)
