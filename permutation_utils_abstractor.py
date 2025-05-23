"""Define some internal utils tools used in LabelledMap and MutableLabelledMap."""

import numpy as np
from map_permutation import MapPermutation


class PermutationUtilsAbstractor:
    """
    This class abstract some utils use in LabelledMap and MutableLabelledMap so that they can use the same
    apis but with different underlying implementation hence the version for MutableLabelledMap inherit from
    this class and is call RotatingUtilsAbstractor
    """

    def __init__(self, permutation: MapPermutation) -> None:
        """
        Init the PermutationUtilsAbstractor

        INPUT:

        - ``permutation`` -- MapPermutation

        EXAMPLES:

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor
            sage: t = MapPermutation([(1,2,4)])
            sage: tAbstr = PermutationUtilsAbstractor(t)

        NOTE:

            O(n) where n is the size of permutation
        """
        cycles = permutation.to_cycles()

        self._numberOfCycles = len(cycles)

        self._numberOfFixedPoint = sum(len(c) == 1 for c in cycles)

        self._cyclesLength = [len(c) for c in cycles]

        self._cycleIndexes = np.zeros(permutation.size() + 1, dtype=int)

        for j, c in enumerate(cycles):
            for i in c:
                self._cycleIndexes[i] = j

    def numberInCycle(self, index: int) -> int:
        """
        INPUT:

        - ``index`` -- int

        OUTPUT:

            The size of the cycle containing index

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor
            sage: t = MapPermutation([(1,2,4)])
            sage: tAbstr = PermutationUtilsAbstractor(t)
            sage: tAbstr.numberInCycle(2)
            3

        NOTE:

            O(1)
        """
        return self._cyclesLength[self._cycleIndexes[index]]

    def sameCycle(self, i: int, j: int) -> bool:
        """
        INPUT:

        - ``i`` -- int
        - ``j`` -- int

        OUTPUT:

            A boolean indicating if i and j are on the same cycle

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor
            sage: t = MapPermutation([(1,2,4)])
            sage: tAbstr = PermutationUtilsAbstractor(t)
            sage: tAbstr.sameCycle(2,3)
            False
            sage: tAbstr.sameCycle(2,1)
            True

        NOTE:

            O(1)
        """
        return bool(self._cycleIndexes[i] == self._cycleIndexes[j])

    def numberOfCycles(self) -> int:
        """
        OUTPUT:

            The number of cycles of the permutation


        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor
            sage: t = MapPermutation([(1,2,4)])
            sage: tAbstr = PermutationUtilsAbstractor(t)
            sage: tAbstr.numberOfCycles()
            2

        NOTE:

            O(1)
        """
        return self._numberOfCycles

    def numberOfFixedPoint(self) -> int:
        """
        OUTPUT:

            The number of fixed point of the permutation

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor
            sage: t = MapPermutation([(1,2,4)])
            sage: tAbstr = PermutationUtilsAbstractor(t)
            sage: tAbstr.numberOfFixedPoint()
            1

        NOTE:

            O(1)
        """
        return self._numberOfFixedPoint

    def checkTwoInTheSameCycle(self, listIndexes: list[int]) -> bool:
        """
        INPUT:

        - ``listIndexes`` -- List[int]

        OUTPUT:
            A boolean indicating if there are two indices in listIndexes on the sameCycle

        EXAMPLES::

            sage: from sage.graphs.maps.map_permutation import MapPermutation
            sage: from sage.graphs.maps.permutation_utils_abstractor import PermutationUtilsAbstractor
            sage: t = MapPermutation([(1,2,4)])
            sage: tAbstr = PermutationUtilsAbstractor(t)
            sage: tAbstr.checkTwoInTheSameCycle([1,2,3])
            True
            sage: tAbstr.checkTwoInTheSameCycle([1,3])
            False

        NOTE:

            O(len(listIndexes))
        """
        checkSet = set()
        for i in listIndexes:
            if self._cycleIndexes[i] in checkSet:
                return True
            checkSet.add(self._cycleIndexes[i])
        return False
