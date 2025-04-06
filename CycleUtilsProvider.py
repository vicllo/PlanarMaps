from SplayTree import SplayTree, Node


class CycleUtilsProvider:
    """
    This class is an abstraction of some utils used by RotatingPermutation
    It is mainly use for abstracting some operations
    for now I didn't add an exhaustive documentation
    """

    def __init__(self, cycles):
        self.nodeMap = {}
        for c in cycles:
            if len(c) == 1:
                continue
            splayTree = SplayTree(range(len(c)))
            for i in range(len(c)):
                e = c[i]
                self.nodeMap[e] = splayTree.getNode(i)
                self.nodeMap[e].index = e
    # OK
    # Return the number of element in the same cycle as index
    # O(log(n))

    def numberInCycle(self, index):
        if self.isFixedPoint(index):
            return 1
        return self.getSplayTree(index).size()
    # OK
    # Return whether i and j are in the same cycle
    # O(log(n))

    def sameCycle(self, i, j):
        if self.isFixedPoint(i) or self.isFixedPoint(j):
            return False
        rootI = self.nodeMap[i].getRoot()
        rootJ = self.nodeMap[j].getRoot()
        self.nodeMap[i].SafeSplay()
        return rootI == rootJ

    # If otherIndex is a fixed node it will add it after index
    # Otherwise it will raise an Error
    # O(log(n))
    def addAfter(self, index, otherIndex):
        assert index != otherIndex
        if not self.isFixedPoint(otherIndex):
            raise ValueError(f"{otherIndex} isn't a fixed point ")
        splayTree = None
        value = 0
        if index not in self.nodeMap:
            splayTree = SplayTree()
            splayTree.insert(index)
            self.nodeMap[index] = splayTree.getNode(index)
            value = index
        else:
            splayTree = self.nodeMap[index].getSplayTree()
            value = splayTree.root.value + splayTree.root.offset

        leftSplayTree, rightSplayTree = splayTree.split(value + 1 / 2)
        leftSplayTree.insert(value + 1)

        self.nodeMap[otherIndex] = leftSplayTree.getNode(value + 1)

        leftSplayTree.shift(-1)
        rightSplayTree.merge(leftSplayTree)

        self.nodeMap[otherIndex].index = otherIndex
        self.nodeMap[index].index = index

    # Return a boolean indicating if two index in  listIndexes
    # are in the sameCycle efficiently (here O( len(listIndexes)log(n)))
    def checkTwoInTheSameCycle(self, listIndexes):
        mapIndexes = set(listIndexes)
        for e in mapIndexes:
            if self.isFixedPoint(e):
                mapIndexes.remove(e)

        splayTreeMap = set()
        for e in mapIndexes:

            splayTree = self.getSplayTree(e)

            if splayTree in splayTreeMap:
                return True

            splayTreeMap.add(splayTree)

        return False
    # If otherIndex is a fixed node it will add  it before index
    # O(log(n))

    def addBefore(self, index, otherIndex):
        assert index != otherIndex
        if not self.isFixedPoint(otherIndex):
            raise ValueError(f"{otherIndex} isn't a fixed point ")

        self._safeIndex(index)

        self.nodeMap[otherIndex] = self.nodeMap[index]

        self.nodeMap[otherIndex].index = otherIndex

        self.nodeMap.pop(index)

        self.addAfter(otherIndex, index)

    # O(log(n))
    def _safeIndex(self, index):
        if not self.isFixedPoint(index):
            return
        splayTree = SplayTree()
        splayTree.insert(index)
        self.nodeMap[index] = splayTree.getNode(index)
        self.nodeMap[index].index = index

    # Def swapLabel of index and otherIndex but keep their topology
    # OK
    # O(log(n))
    def swapIndex(self, index, otherIndex):
        if index not in self.nodeMap and otherIndex not in self.nodeMap:
            return
        if index not in self.nodeMap and otherIndex in self.nodeMap:
            self.nodeMap[index] = self.nodeMap[otherIndex]
            self.nodeMap[index].index = index
            self.nodeMap.pop(otherIndex)
            return
        if index in self.nodeMap and otherIndex not in self.nodeMap:
            self.nodeMap[otherIndex] = self.nodeMap[index]
            self.nodeMap[otherIndex].index = otherIndex
            self.nodeMap.pop(index)
            return
        tmp = self.nodeMap[index]
        self.nodeMap[index] = self.nodeMap[otherIndex]
        self.nodeMap[otherIndex] = tmp

        self.nodeMap[index].index = index
        self.nodeMap[otherIndex].index = otherIndex

    # OK
    # O(sizeOfCycle)
    def getCycleList(self, index):
        if self.isFixedPoint(index):
            return [index]
        splayTree = self.getSplayTree(index)
        return splayTree.indexList()

    # Detach the node if it isn't a fixed point and make it a fixed point
    # O(log(n))
    # OK
    def detach(self, index):
        if self.isFixedPoint(index):
            return
        splayTree = self.getSplayTree(index)
        value = splayTree.root.value + splayTree.root.offset
        splayTree.delete(value)

        self.nodeMap.pop(index)

        if splayTree.size() == 1:
            node = splayTree.root
            self.detach(node.index)

    # O(log(n))
    # OK
    def isFixedPoint(self, index):
        return index not in self.nodeMap

    # OK
    # O(log(n))
    def getValue(self, index):
        if self.isFixedPoint(index):
            raise ValueError(
                f"A fixed point doesn't have a value: {index} is a fixed point")
        splayTree = self.getSplayTree(index)
        return splayTree.root.value + splayTree.root.offset

    # OK
    # O(log(n))
    def getSplayTree(self, index):
        if self.isFixedPoint(index):
            raise ValueError(
                f"A fixed point doesn't have a splay tree : {index} is a fixed point")
        return self.nodeMap[index].getSplayTree()
    # OK
    # O(log(n))

    def makeMin(self, index):
        if self.isFixedPoint(index):
            return
        value = self.getValue(index)
        splayTree = self.getSplayTree(index)

        if splayTree.min() == value:
            return

        left, right = splayTree.split(value - 1 / 2)
        left.shift(right.max() + 1 - left.min())
        left.merge(right)
    # OK
    # O(log(n))

    def makeMax(self, index):
        if self.isFixedPoint(index):
            return
        value = self.getValue(index)
        splayTree = self.getSplayTree(index)

        if splayTree.max() == value:
            return
        left, right = splayTree.split(value + 1 / 2)
        left.shift(right.max() + 1 - left.min())
        left.merge(right)

    # OK
    # O(log(n))
    def merge(self, beforeIndex, afterIndex):
        if self.sameCycle(beforeIndex, afterIndex):
            raise ValueError(
                f"Cannot merge detache two indexes in the same cycle  : {beforeIndex} and {afterIndex}")

        self._safeIndex(beforeIndex)
        self._safeIndex(afterIndex)

        self.makeMax(beforeIndex)
        self.makeMin(afterIndex)

        leftTree = self.getSplayTree(beforeIndex)
        rightSplayTree = self.getSplayTree(afterIndex)

        rightSplayTree.shift(leftTree.max() + 1 - rightSplayTree.min())
        rightSplayTree.merge(leftTree)

    # OK
    # O(log(n))
    def cut(self, startIndex, endIndex):
        if not self.sameCycle(startIndex, endIndex):
            raise ValueError(
                f"Cut can only be called for index in the same cycle {startIndex} and {endIndex} are not")
        if startIndex == endIndex:
            self.detach(startIndex)
            return
        startValue = self.getValue(startIndex)
        endValue = self.getValue(endIndex)
        splayTree = self.getSplayTree(startIndex)

        if startValue < endValue:
            left, otherRight = splayTree.split(endValue + 1 / 2)
            otherLeft, _ = left.split(startValue - 1 / 2)
            otherRight.merge(otherLeft)
        else:
            otherLeft, right = splayTree.split(startValue - 1 / 2)
            left, otherRight = otherLeft.split(endValue + 1 / 2)
            left.shift(right.max() + 1 - left.min())
            left.merge(right)
