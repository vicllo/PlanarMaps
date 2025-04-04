from time import time
from CycleUtilsProvider import CycleUtilsProvider
from CyclicChainedList import CyclicChainedList
from sage.all import Permutation
from MapPermutation import MapPermutation


class RotatingPermutation(MapPermutation):

    """
    A class representing permutation where it is fast to:
    - delete (O(log(n))) element,
    - check if two indices are in the same cycle (O(log(n))),
    - add (O(log(n))) element in its cycles representation,
    - and more things useful in MutableLabelledMap.

    Note that compared to simple MapPermutation,
    RotatingPermutation are more heavy objects; hence they are more demanding when initializing.
    If you don't need all the power of RotatingPermutation, consider using the simple MapPermutation.

    Another thing: for compatibility reasons between MutableLabelledMap and LabelledMap,
    every method that returns a permutation must return MapPermutation.
    Hence, don't assume that the permutation you get is a RotatingPermutation;
    you should do it yourself.

    WARNING: We take as a convention for this class that if i is bigger than the size of self,
    then self(i) = i.
      """

    def __init__(self, lst):
        """
        This function initiate the rotating permutation, lst can be  a Permutation or a list of int or list of tuple representing the cycle of
        the permutation or a MapPermutation or an integer representing the size of the permutation(in this case self will represent the identity permutation of size lst).
        """
        if isinstance(lst, Permutation) or isinstance(lst, MapPermutation):
            self.__init__(list(lst))
            return

        # Important will stay
        # Associate to a node its image
        self._permCycle = {}
        # Size of the permutation
        self._n = 0
        # Number of cycles in the permutation
        self._numCycles = 0

        # Number of fixed point in the permutation
        self._numberOfFixedPoint = 0

        try:
            if lst == int(lst) and lst > 0:
                # If lst is an integer we just set our permutation to be the
                # identity
                self._n = lst
                self._numCycles = self._n
                self.provider = CycleUtilsProvider([])
                return
        except BaseException:
            pass

        mx = 0
        seen = []
        try:
            # We're directly using the cycle representation to initialise the
            # permutation
            if isinstance(lst[0], type((42,))):
                for l in lst:
                    for i in l:
                        mx = max(i, mx)
                        if i != int(i) or i <= 0:
                            raise ValueError(
                                f"Invalid argument: {i} isn't a strictly positive integer in the list given")
                seen = [False for i in range(mx + 1)]
                cnt = 0
                for l in lst:
                    k = 0
                    while k < len(l):
                        i = l[k]
                        if seen[i]:
                            raise ValueError(
                                f"Invalid argument: {i} appears at least two times in list given it cannot be a permutation.")
                        seen[i] = True
                        k += 1
                        while k < len(l) and l[k] == i:
                            k += 1
                            continue
                        cnt += 1
                self._numCycles += mx - cnt
                self._numberOfFixedPoint += mx - cnt
                for l in lst:
                    prevNode = None
                    for i in l:
                        newNode = CyclicChainedList(i)
                        self._permCycle[i] = newNode
                        if prevNode is not None:
                            prevNode.insertAfter(newNode)
                        prevNode = newNode
                    self._numberOfFixedPoint += len(l) == 1
                    self._numCycles += 1

            else:
                mx = len(lst)
                for i in lst:
                    if i != int(i) or i <= 0:
                        raise ValueError(
                            f"Invalid argument : {i} isn't a strictly positive integer in the list given")
                    if i > len(lst):
                        raise ValueError(
                            f"{i} is bigger than the size of the given list")
                seen = [False for i in range(mx + 1)]
                for i in lst:
                    if seen[i]:
                        raise ValueError(
                            f"Invalid argument: {i} appears at least two time in the list given it cannot be a permutation..")
                    seen[i] = True

                seen = [False for i in range(mx + 1)]
                for i in range(1, mx + 1):
                    if seen[i]:
                        continue
                    prevNode = None
                    curElement = i
                    cnt = 0
                    while not seen[curElement]:
                        newNode = CyclicChainedList(curElement)
                        self._permCycle[curElement] = newNode
                        if prevNode is not None:
                            prevNode.insertAfter(newNode)
                        prevNode = newNode
                        seen[curElement] = True
                        curElement = lst[curElement - 1]
                        cnt += 1
                    self._numCycles += 1
                    self._numberOfFixedPoint += cnt == 1
        except ValueError as e:
            raise
        except BaseException:
            raise ValueError("Invalid argument: The argument given must be Permutation or MapPermutation or a non empty list of integers representing the permutation or a non empty list of tuples representing the cycles of the permutations or a positive integer.")
        self._n = mx
        self.provider = CycleUtilsProvider(self.to_cycles())

    # OK

    def size(self):
        """
        Returns: the size of the permutation
        -------
        O(1)
        """
        return self._n

    # OK
    def deleteLastKIndex(self, k):
        """
        This function will delete the last k index from self
        ------
        Args:
            - k the number of node to delete
        """
        if k > self.size():
            raise ValueError(
                f"Cannot delete {k} last element in a RotatingPermutation of size {self.size()}")
        for _ in range(k):
            self.delete(self._n)

    # OK

    def delete(self, index):
        """
        This will delete index of the corresponding cycle note that after this operation if we note the original
        size of self as n, the which contained index will count one less element,
        self will be of size n-1 and if n != index the element numbered n will relabeled as index.
        For instance if self is the permutation(1, 2, 3)(4, 5) and we delete 2 it will become(1, 3)(4, 2),
        If n = 1 an error or index is not a strictly positive integer <= n an error will be raised.
        -------
        Args:
            -index: an integer representing the index to delete

        Note that if index must be an strictly positive integer and self.size() >= 2 otherwise an error will be raised
        -------
        O(log(n))
        """
        if self.size() == 1:
            raise ValueError(
                "Cannot delete an element from a Permutation of size 1")
        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                "{index} isn't a strictly positive integer <= self.size()")

        nPrev = self.size()
        node = self.getNode(index)

        node.remove()

        if self.provider.numberInCycle(index) == 2:
            self._numberOfFixedPoint += 1

        if self.provider.numberInCycle(index) == 1:
            self._numberOfFixedPoint -= 1
            self._numCycles -= 1

        self._n -= 1

        self._permCycle.pop(index)

        self.provider.swapIndex(nPrev, index)
        self.provider.detach(nPrev)

        if nPrev != index:
            try:
                self._permCycle[nPrev].val = index

                self._permCycle[index] = self._permCycle[nPrev]

                self._permCycle.pop(nPrev)
            except BaseException:
                pass

    # OK
    def inverseApply(self, i):
        """
        This function apply  the inverse self on i, we take as a convention i if i is an integer > self.size(), self.inverseApply(i) = i
        ------
        Args:
            i an index
        ------
        O(1)
        """
        if i != int(i) or i <= 0:
            raise ValueError("{i} isn't a positive integer")
        try:
            return self._permCycle[i].prev.val
        except BaseException:
            return i
    # OK

    def checkTwoInTheSameCycle(self, listIndexes):
        """
        This function will return a boolean indicating if there is two index in listIndexes in the sameCycle
        ------
        Args:
            listIndexes: A list of indexes
        Returns:
            A boolean indicating if two indexes are in the same cycle
        ------
        O(plog(n)) where p = len(listIndexes)
        """
        return self.provider.checkTwoInTheSameCycle(listIndexes)
    # OK

    def swapIndex(self, index, otherIndex):
        """
        This function swap the index role in the permutation
        ------
        Args:
            index, otherIndex the two indexes <= self.size()
        ------
        O(log(n))
        """
        self.provider.swapIndex(index, otherIndex)
        nodeIndex = self.getNode(index)
        nodeOther = self.getNode(otherIndex)
        self._permCycle[otherIndex] = nodeIndex
        self._permCycle[index] = nodeOther
        self._permCycle[otherIndex].val = otherIndex
        self._permCycle[index].val = index

    # OK
    def cutAdd(self, startIndex, endIndex, newIndexStart, newIndexEnd):
        """
        This implement a special operation.In a nutshell it cut a cycle and add two index in each cycle,
        let denote A = startIndex, B = endIndex, C = newIndexStart, D = newIndexEnd and say the cycle is of the form F -> A -> S -> .. -> T -> B -> R -> ... -> F
        than the situation will be the following after a call to this function, A -> S -> ... -> T -> D -> A and F -> C -> B -> R -> ... -> F
        ------
        Args:
            startIndex, endIndex, newIndexStart, newIndexEnd: 4 indexes, startIndex and endIndex must be on the same cycle
            and {newIndexEnd, newIndexStart} = {n+1, n+2} and should be fixed point
        ------
        O(log(n))
        """
        if newIndexEnd == newIndexStart:
            raise ValueError(
                f"{newIndexEnd} and {newIndexStart} must be different")
        if newIndexStart <= self.size() or newIndexStart > self.size() + 2:
            raise ValueError(
                f"{newIndexStart} must be  >{self.size()} and <= {self.size() + 2}")
        if newIndexEnd <= self.size() or newIndexEnd > self.size() + 2:
            raise ValueError(
                f"{newIndexEnd} must be  >{self.size()} and <= {self.size() + 2}")
        if not self.sameCycle(startIndex, endIndex):
            raise ValueError(
                f"{newIndexEnd} and {newIndexStart} must be in the same cycle to use cutAdd")
        if startIndex == endIndex:
            self.addBefore(startIndex)
            self.addBefore(startIndex)
            return

        nodeStartIndex = self.getNode(startIndex)
        nodeEndIndex = self.getNode(endIndex)

        # Updating scalar attribute accordingly
        self._numCycles += 1
        self._n += 2

        nodeNewIndexStart = self.getNode(newIndexStart)
        nodeNewIndexEnd = self.getNode(newIndexEnd)

        # NodeNewIndexStart processing
        comeBeforeEnd = self.inverseApply(endIndex)
        tmpNode = self.getNode(self.inverseApply(startIndex))
        nodeNewIndexStart.prev = tmpNode
        tmpNode.nxt = nodeNewIndexStart

        nodeNewIndexStart.nxt = nodeEndIndex

        # NodeNewIndexEnd processing
        tmpNode = self.getNode(self.inverseApply(endIndex))
        nodeNewIndexEnd.prev = tmpNode
        tmpNode.nxt = nodeNewIndexEnd

        nodeNewIndexEnd.nxt = nodeStartIndex

        nodeEndIndex.prev = nodeNewIndexStart
        nodeStartIndex.prev = nodeNewIndexEnd

        # Updating the provider
        self.provider.cut(startIndex,
                          comeBeforeEnd)

        self.provider.addBefore(startIndex, newIndexEnd)
        self.provider.addBefore(endIndex, newIndexStart)

    # OK
    def labelToTheEnd(self, listIndexes):
        """
        This is a helper function  it just move all of the element in listIndexes to the last indices
        -----
        Args:
            listIndexes
        -----
        O(len(listIndexes)*log(n))
        """
        for index in listIndexes:
            if index != int(index) or index <= 0 or index > self.size():
                raise ValueError(
                    f"In labelToTheEnd : {index} isn't a strictly positive integer <= {self.size()}")

        indexMap = set()
        for index in listIndexes:
            indexMap.add(index)

        indexCandidate = set()
        for j in range(len(indexMap)):
            indexCandidate.add(self.size() - j)

        for index in list(indexCandidate):
            if index in indexMap:
                indexCandidate.remove(index)
                indexMap.remove(index)

        corresOut = {}
        for index in list(indexMap):
            if index not in indexMap:
                continue
            corresIndex = indexCandidate.pop()
            corresOut[index] = corresIndex
            indexMap.remove(index)
            self.swapIndex(index, corresIndex)
        return corresOut
    # OK

    def bruteAddCycles(self, cycles):
        """
        Another helper function that add cyclein cycles, this one assumed is more dangerous than addCycles
        cause it assumed that the cycles are well formed and not > self.size()
        thus the term brute
        ----
        Args:
            cycles: list of cycles
        ----
        O(len(cycles)*log(n))
        """
        for c in cycles:
            for i in range(len(c) - 1):
                self.addAfterGeneral(c[i], c[i + 1])

    # OK
    def addCycles(self, cycles):
        """
        Another helper function it will raise an error if element of the cycles
        are not > self.size() and <= self.size()+len(cycles), the cycle must be well formed
        ----
        Args:
            cycles: list of cycles
        ----
        O(len(cycles)*log(n))
        """

        testSet = set()
        N = 0
        for c in cycles:
            N += len(c)

        for c in cycles:
            for e in c:
                testSet.add(e)
                if e <= self.size() or e <= 0 or e != int(e) or e > self.size() + N:
                    raise ValueError("{cycles} isn't valid")
        if len(testSet) != N:
            raise ValueError("{cycles} isn't valid")

        self.stretch(N)

        self.bruteAddCycles(cycles)

    # OK
    def isValidIndex(self, index):
        """
        Check if index is a integer > 0 and <=self.size()
        otherwise raise an Error
        ----
        Args:
            index
        ----
        O(1)
        """
        if index <= 0 or index != int(index) or index > self.size():
            raise ValueError(f"{index} isn't valid")
    # OK

    def addAfterGeneral(self, index, otherIndex):
        """
        This is a more general version of addAfter it only assumed that otherIndex is a fixed point
        and will add it after index in its cycle
        -----
        Args:
            index, otherIndex
        -----
        O(log(n))
        """
        self.isValidIndex(index)
        self.isValidIndex(otherIndex)
        if index == otherIndex:
            return
        if not self.provider.isFixedPoint(otherIndex):
            raise ValueError(
                f"Can only add after fixed point {otherIndex} isn't one")

        self._numberOfFixedPoint -= self.provider.isFixedPoint(index)
        self.provider.addAfter(index, otherIndex)

        node = self.getNode(index)

        newNode = self.getNode(otherIndex)

        self._numberOfFixedPoint -= 1
        self._numCycles -= 1

        node.insertAfter(newNode)

    # OK
    def addBeforeGeneral(self, index, otherIndex):
        """
        More general version of addBeforeit only assumed that otherIndex is a fixed point
        and will add it before index in its cycle
        ----
        Args:
            index, otherIndex
        ----
        O(log(n))
        """
        self.isValidIndex(index)
        indexPrev = self.inverseApply(index)
        self.addAfterGeneral(indexPrev, otherIndex)

    # OK
    def mergeDelete(self, index, otherIndex):
        """
        Assuming that index and otherIndex are not in the same cycle it will do the
        following first index and otherIndex will be sent to self.size() self.size()-1 they will be deleted and given
        that before we add: U -> ... -> V -> index -> R -> U and F -> ... -> T -> otherIndex -> Q -> F, we will have after
        U -> ... -> V -> Q -> F -> ... -> T -> R -> U
        ----
        index, otherIndex two node not on the same cycle
        ----
        O(log(n))
        """
        if self.sameCycle(index, otherIndex):
            raise ValueError("Cannot merge delete two index on the sameCycle")

        backUpNumberOfFixedPoint = self.number_of_fixed_points()

        self.labelToTheEnd([index, otherIndex])

        if self.provider.isFixedPoint(
                self._n) or self.provider.isFixedPoint(self._n - 1):
            self.deleteLastKIndex(2)
            return

        beforeIndex = self.inverseApply(self._n)
        afterIndex = self.apply(self._n - 1)

        self.deleteLastKIndex(2)

        nodeBefore = self.getNode(beforeIndex)
        nodeAfter = self.getNode(afterIndex)

        if nodeBefore.nxt == nodeBefore:
            if nodeAfter.nxt == nodeAfter:
                nodeBefore.nxt = nodeAfter
                nodeBefore.prev = nodeAfter
                nodeAfter.nxt = nodeBefore
                nodeAfter.prev = nodeBefore
            else:
                tmpNode = nodeAfter.prev
                tmpNode.nxt = nodeBefore
                nodeBefore.prev = tmpNode
                nodeAfter.prev = nodeBefore
                nodeBefore.nxt = nodeAfter

        else:
            if nodeAfter.nxt == nodeAfter:
                tmpNode = nodeBefore.nxt
                tmpNode.prev = nodeAfter
                nodeAfter.nxt = tmpNode
                nodeAfter.prev = nodeBefore
                nodeBefore.nxt = nodeAfter
            else:
                tmpNodeBefore = nodeBefore.nxt
                tmpNodeAfter = nodeAfter.prev
                tmpNodeBefore.prev = tmpNodeAfter
                tmpNodeAfter.nxt = tmpNodeBefore

                nodeBefore.nxt = nodeAfter
                nodeAfter.prev = nodeBefore

        self._numCycles -= 1
        self._numberOfFixedPoint = backUpNumberOfFixedPoint
        self.provider.merge(beforeIndex, afterIndex)

    # OK
    def getNode(self, index):
        """
        This function will return the node associated to index
        and if it doesn't exit it will create one note that if index > self.size()
        it will raise an error.
        -----
        O(1)
        -----
        """
        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                "{index} isn't a strictly positive integer <= self.size()")

        try:
            node = self._permCycle[index]
        except BaseException:
            node = CyclicChainedList(index)
            self._permCycle[index] = node
        return node

    # OK

    def stretch(self, m):
        """
        This function will increase the size of the permutation by m,all the new index will
        be fixed point
        -----
        O(1)
        -----
        """
        self._n += m
        self._numberOfFixedPoint += m
        self._numCycles += m
    # OK

    def addAfter(self, index):
        """
        Let denote n=self.size() given that  n>=index>=1, this will increase the size of self by one and add
        the new element n+1 on the cycle of index after index.You should note that if index>self.size() this will raise an error.
        -----
        O(log(n))
        """
        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                f"{index} isn't a strictly positive integer <= {self.size()}")

        self._numberOfFixedPoint -= self.provider.isFixedPoint(index)
        nPrev = self.size()

        self.stretch(1)

        self.provider.addAfter(index, nPrev + 1)

        node = self.getNode(index)

        newNode = self.getNode(nPrev + 1)

        self._numberOfFixedPoint -= 1
        self._numCycles -= 1

        node.insertAfter(newNode)
    # OK

    def addBefore(self, index):
        """
        Let denote n=self.size() given that  n>=index>=1, this will increase the size of self by one and add
        the new element n+1 on the cycle of index before index.You should note that if index>self.size() this will raise an error.
        -----
        O(log(n))
        """
        if index != int(index) or index <= 0 or index > self.size():
            raise ValueError(
                "{index} isn't a strictly positive integer <= self.size()")

        node = self.getNode(index)
        prevIndex = node.prev.val

        self.addAfter(prevIndex)

    # OK

    def numberInCycle(self, index):
        """
        Args:
            -index : A strictly positive integer
        Returns:
            -A integer representing the number of element in the same cycle as index note that
            if index > self.size() it will return 1(which is coherent with the convention that self(i) = i)
        -------------
        O(log(n))
        """
        return self.provider.numberInCycle(index)

    # OK
    def numberOfCycles(self):
        """
        Returns: the number of cycle of self
        -------
        O(1)
        """
        return self._numCycles

    # OK
    def sameCycle(self, i, j):
        """
        Args:
            -i an strictly positive integer
            -j an strictly positive integer

        Returns:
            A boolean indicating whether of not i and j are on the same cycle of self
        -------
        O(log(n))
        """
        if i <= 0 or j <= 0 or i != int(i) or j != int(j):
            raise ValueError("{i} or {j} isn't a strictly positive integer")

        return self.provider.sameCycle(i, j)

    # OK
    def __repr__(self):
        return str(list(self))

    # OK

    def pretty_repr(self):
        """
        Return a string of self in the form of his cycle decomposition
        """
        return f"Rotating permutation: {self.to_cycles()}"

    # OK
    def pretty_print(self):
        """
        Print self in a more pretty form
        """
        print(self.pretty_repr())

    # OK
    def to_cycles(self):
        """
        This method calculate a list of tuple representing the cycle of self
        -------
        Returns:
            - lst a list of tuples representing the cycles of self given in increasing order of their minimum elements
        -------
        O(n)
        where n is the number of element of self
        """
        seen = [False for i in range(self.size() + 1)]
        cycles = []
        for i in range(1, self.size() + 1):
            if seen[i]:
                continue
            try:
                node = self._permCycle[i]
                cycle = node.getValList()
                cycles.append(tuple(cycle))
                for j in cycle:
                    seen[j] = True
            except BaseException:
                cycles.append((i,))

        return cycles

    # OK
    def inverse(self):
        """
        This function calculate  the inverse of self
        -------
        Returns:
            - The inverse of self
        -------
        O(n)
        where n is the number of element of the permutation
        -------
        """
        cycles = self.to_cycles()
        return MapPermutation([tuple(reversed(e)) for e in cycles])

    # OK
    def apply(self, i):
        """
        This function apply self on i , we take as a convention i if i is an integer > self.size() , self.apply(i) = i
        """
        if i != int(i) or i <= 0:
            raise ValueError("{i} isn't a positive integer")
        try:
            return self._permCycle[i].nxt.val
        except BaseException:
            return i

    # OK

    def number_of_fixed_points(self):
        """
        Returns: the number of fixed point ( we only consider i such that i<=self.size())
        """
        return self._numberOfFixedPoint

    def __eq__(self, other):
        if isinstance(other, MapPermutation):
            return list(other) == list(self)
        return False
