
# This file contain an implementation of a custom
# Variant of SplayTree used in CycleUtilsProviderMain
# Mainly because it isn't the main focus I didn't add exhaustive
# documentation for now


# OK
def swapRight(node, otherNode):
    tmp = node.right
    node.right = otherNode.right
    otherNode.right = tmp

    if node.right is not None:
        node.right.parent = node

    if otherNode.right is not None:
        otherNode.right.parent = otherNode


# OK
def swapLeft(node, otherNode):
    tmp = node.left
    node.left = otherNode.left
    otherNode.left = tmp

    if node.left is not None:
        node.left.parent = node

    if otherNode.left is not None:
        otherNode.left.parent = otherNode


# OK
def SwapNonTopologicalExceptIndex(node, otherNode):
    # Value
    tmp = node.value
    node.value = otherNode.value
    otherNode.value = tmp

    # Offset
    tmp = node.offset
    node.offset = otherNode.offset
    otherNode.offset = tmp

    # cnt
    tmp = node.cnt
    node.cnt = otherNode.cnt
    otherNode.cnt = tmp

    # SplayTree
    tmp = node.splayTree
    node.splayTree = otherNode.splayTree
    otherNode.splayTree = tmp


# OK
def makeParentKnow(node):
    if node.parent is None:
        return
    parent = node.parent
    if isLeftChild(parent, node):
        parent.left = node
    else:
        parent.right = node


# OK
def swapNodeButNotIndex(node, otherNode):

    # OK
    if node == otherNode:
        return

    if node.parent == otherNode.parent:
        swapLeft(node, otherNode)
        swapRight(node, otherNode)
        SwapNonTopologicalExceptIndex(node, otherNode)
        makeParentKnow(node)
        makeParentKnow(otherNode)
        return

    # OK
    if node.parent == otherNode:
        swapNodeButNotIndex(otherNode, node)
        return

    # OK
    if otherNode.parent == node:
        swapLeft(node, otherNode)
        swapRight(node, otherNode)

        otherNode.parent = node.parent

        SwapNonTopologicalExceptIndex(node, otherNode)

        makeParentKnow(otherNode)

        node.parent = otherNode
        return

    # Left
    swapLeft(node, otherNode)
    # Right
    swapRight(node, otherNode)
    # Non topological field(i.e != left,right,parent) and != index
    SwapNonTopologicalExceptIndex(node, otherNode)

    # Parent
    tmp = node.parent
    node.parent = otherNode.parent
    otherNode.parent = tmp

    makeParentKnow(node)
    makeParentKnow(otherNode)


def numberOfElement(node):
    return 0 if node is None else node.cnt


def valueToTheLeft(parentValue, value):
    return parentValue > value


def valueToTheRight(parentValue, value):
    return parentValue < value


def isLeftChild(parentNode, node):
    return valueToTheLeft(parentNode.value, node.value + node.offset)


def isRightChild(parentNode, node):
    return valueToTheRight(parentNode.value, node.value + node.offset)


class Node:
    def __init__(self, value=None, parent=None) -> None:
        self.value = value

        self.left = None
        self.right = None
        self.parent = parent

        self.cnt = 0
        self.offset = 0
        self.splayTree = None
        self.index = None

    def SafeSplay(self):
        oldRoot = self.getRoot()
        if oldRoot != self:
            self.splayTree = oldRoot.splayTree
            oldRoot.splayTree = None
        self.splay()

    def getSplayTree(self):
        # It is guaranteed after colling this function that self is the root
        root = self.getRoot()
        self.splay()
        splayTree = root.splayTree
        if splayTree is not None:
            splayTree.changeRoot(self)
        return splayTree

    def getRoot(self):
        node = self
        while not node.isRoot():
            node = node.parent

        return node

    def height(self):
        if self.isEmpty():
            return -1
        maxChildHeight = -1
        if self.right is not None:
            maxChildHeight = self.right.height()
        if self.left is not None:
            maxChildHeight = max(self.left.height(), maxChildHeight)
        return 1 + maxChildHeight

    def sortedList(self, offset=0):
        if self.isEmpty():
            return []
        out = [self.value + offset + self.offset]
        if self.right is not None:
            out = out + self.right.sortedList(self.offset + offset)
        if self.left is not None:
            out = self.left.sortedList(offset + self.offset) + out
        return out

    def indexList(self):
        if self.isEmpty():
            return []
        out = [self.index]
        if self.right is not None:
            out = out + self.right.indexList()
        if self.left is not None:
            out = self.left.indexList() + out
        return out

    def isEmpty(self):
        return self.value is None

    def findSmallestGreater(self, value):
        if self.isEmpty():
            return self, 0

        node = self
        offset = self.offset
        bestNode = node
        bestOffset = self.offset
        while True:
            if value == node.value + offset:
                return node, offset
            if valueToTheLeft(node.value + offset, value):
                bestNode = node
                bestOffset = offset
                if node.left is None:
                    return bestNode, bestOffset
                node = node.left
            else:
                if node.right is None:
                    return bestNode, bestOffset
                node = node.right
            offset += node.offset

    def findBiggestSmaller(self, value):
        if self.isEmpty():
            return self, 0

        node = self
        offset = self.offset
        bestNode = node
        bestOffset = self.offset
        while True:
            if value == node.value + offset:
                return node, offset
            if valueToTheLeft(node.value + offset, value):
                if node.left is None:
                    return bestNode, bestOffset
                node = node.left
            else:
                bestNode = node
                bestOffset = offset
                if node.right is None:
                    return bestNode, bestOffset
                node = node.right
            offset += node.offset

    def insert(self, newValue):
        if self.isEmpty():
            self.value = newValue
            self.addCountUpward(1)
            return True, self, 0
        node = self
        offset = self.offset
        while (True):
            if node.value + offset == newValue:
                return False, node, offset
            if valueToTheLeft(node.value + offset, newValue):
                if node.left is None:
                    node.left = Node(newValue - offset, node)
                    node.left.addCountUpward(1)
                    return True, node.left, offset
                else:
                    node = node.left
            else:
                if node.right is None:
                    node.right = Node(newValue - offset, node)
                    node.right.addCountUpward(1)
                    return True, node.right, offset
                else:
                    node = node.right
            offset += node.offset

    def addCountUpward(self, toAdd):
        node = self
        node.cnt += toAdd
        while not node.isRoot():
            node = node.parent
            node.cnt += toAdd

    def find(self, value):
        if self.isEmpty():
            return self, 0
        node = self
        offset = self.offset
        while True:
            if value == node.value + offset:
                return node, offset
            if valueToTheLeft(node.value + offset, value):
                if node.left is None:
                    return node, offset
                node = node.left
            else:
                if node.right is None:
                    return node, offset
                node = node.right
            offset += node.offset

    def min(self):
        minimum = self
        offset = self.offset
        while minimum.left is not None:
            minimum = minimum.left
            offset += minimum.offset
        return minimum, offset

    def delete(self, value):
        if self.isEmpty():
            return self

        node, offset = self.find(value)
        if node.value + offset != value:
            return node

        if node.left is None and node.right is None:
            if node.isRoot():
                node.value = None
                node.addCountUpward(-1)
                return node
            node.parent.addCountUpward(-1)
            if isLeftChild(node.parent, node):
                node.parent.left = None
                return node.parent
            node.parent.right = None
            return node.parent

        if node.left is None and node.right is not None:
            node.right.parent = node.parent
            if node.parent is None:
                return node.right

            if (isLeftChild(node.parent, node)):
                node.parent.left = node.right
            else:
                node.parent.right = node.right

             # Offset
            node.right.offset += node.offset

            node.parent.addCountUpward(-1)
            return node.right

        if node.right is None and node.left is not None:
            node.left.parent = node.parent
            if node.parent is None:
                return node.left
            if isLeftChild(node.parent, node):
                node.parent.left = node.left
            else:
                node.parent.right = node.left

            # Offset
            node.left.offset += node.offset

            node.parent.addCountUpward(-1)
            return node.left

        rep, repOffset = node.right.min()

        if isLeftChild(rep.parent, rep):
            rep.parent.left = rep.right
        else:
            rep.parent.right = rep.right

        if rep.right is not None:
            rep.right.parent = rep.parent
            rep.right.offset += rep.offset

        node.value = rep.value + repOffset
        rep.parent.addCountUpward(-1)

        # This line is mainly here to not break the association
        # in CycleUtilsProviderMain between index <-> node
        swapNodeButNotIndex(node, rep)

        return node.parent

    def max(self):
        maximum = self
        offset = self.offset
        while maximum.right is not None:
            maximum = maximum.right
            offset += maximum.offset
        return maximum, offset

    def _isBst(self, offset=0):
        if self.isEmpty():
            return None, None
        min = self.value + self.offset + offset
        max = self.value + self.offset + offset
        if self.left is not None:
            minLeft, maxLeft = self.left._isBst(self.offset + offset)
            assert valueToTheLeft(self.value + offset + self.offset, maxLeft)
            min = minLeft
        if self.right is not None:
            minRight, maxRight = self.right._isBst()
            assert valueToTheRight(self.value + offset + self.offset, minRight)
            max = maxRight
        return min, max

    def isBst(self):
        self._isBst()

    def rightRotation(self):
        if self.isEmpty():
            return self

        oldLeft = self.left

        if oldLeft is None:
            return self

        # Parent
        if self.parent is not None:
            parent = self.parent
            if isLeftChild(parent, self):
                parent.left = oldLeft
            else:
                parent.right = oldLeft

        if oldLeft.right is not None:
            oldLeft.right.parent = self

        oldLeft.parent = self.parent
        self.parent = oldLeft

        # Offset
        x = oldLeft.offset
        y = self.offset

        oldLeft.offset = x + y
        self.offset = -x
        if oldLeft.right is not None:
            oldLeft.right.offset += x

        # Cnt
        totCnt = self.cnt

        self.cnt = totCnt - oldLeft.cnt + numberOfElement(oldLeft.right)

        oldLeft.cnt = totCnt

        # Child

        self.left = oldLeft.right
        oldLeft.right = self

        return oldLeft

    def leftRotation(self):
        if self.isEmpty():
            return self

        oldRight = self.right
        if oldRight is None:
            return self

        # Parent
        if self.parent is not None:
            parent = self.parent
            if isLeftChild(parent, self):
                parent.left = oldRight
            else:
                parent.right = oldRight

        if oldRight.left is not None:
            oldRight.left.parent = self

        oldRight.parent = self.parent
        self.parent = oldRight

        # Offset
        y = self.offset
        x = oldRight.offset
        oldRight.offset = x + y
        self.offset = -x
        if oldRight.left is not None:
            oldRight.left.offset += x

        # Cnt
        totCnt = self.cnt

        self.cnt = totCnt - oldRight.cnt + numberOfElement(oldRight.left)

        oldRight.cnt = totCnt

        # Child
        self.right = oldRight.left
        oldRight.left = self

        return oldRight

    def isRoot(self):
        return self.parent is None

    def getOffset(self):
        offset = self.offset
        node = self
        while not node.isRoot():
            node = node.parent
            offset += node.offset
        return offset

    def splay(self):
        while not self.isRoot():
            parent = self.parent

            if parent.isRoot():
                if isLeftChild(parent, self):
                    parent.rightRotation()
                else:
                    parent.leftRotation()
                continue

            grandparent = parent.parent

            if isLeftChild(grandparent, parent) and isLeftChild(parent, self):
                grandparent.rightRotation()
                parent.rightRotation()

            elif isRightChild(grandparent, parent) and isRightChild(parent, self):
                grandparent.leftRotation()
                parent.leftRotation()

            elif isRightChild(grandparent, parent) and isLeftChild(parent, self):
                parent.rightRotation()
                grandparent.leftRotation()

            else:
                parent.leftRotation()
                grandparent.rightRotation()
        return self


class SplayTree():
    def __init__(self, lst=[], root=None):
        """
        This class implement a special version of splay tree where it is possible to schift all the
        value contained in the splay tree by a constant.
        -------
        Initialises the splay tree
        -------
        Args:
            -lst: A list of element to add too self
            -root: A possible Node on which to root the splayTree if set to None it will be set to
                Node()
        -------
        O(n * log(max(n,m)) ) where n is the length of lst and m is the number of element in the subtree
            of root
        """

        self.root = root
        if self.root is None:
            self.root = Node()
        self.root.splayTree = self
        self.valid = True
        self.insertList(lst)

    def changeRoot(self, root):
        if not root.isRoot():
            raise ValueError("The argument must be the root of his tree")
        self.root.splayTree = None
        self.root = root
        self.root.splayTree = self

    def shift(self, toAdd):
        """
        This will shift all the value contained in the tree by toAdd
        -------
        Args:
            toAdd
        -------
        O(1)
        """
        self.checkValid()
        self.root.offset += toAdd

    def checkValid(self):
        """
        Check whether self is a valid SplayTree if not raise an Error
        """
        if self.valid:
            return
        raise ValueError("This isn't a valid instance of a splayTree anymore")

    def split(self, value):
        """
        This will split self into two splayTree a,b such that a contains all the element <= value
        and b all the element >value.You should note that a == self after the operation.
        -------
        Args:
            value : The value on which to split
        Returns:
            (a,b) such that a contains all the element <= value and b all the element > value
        -------
        O(log(n)) where n is the size of self
        """
        self.checkValid()
        if self.isEmpty():
            return self, SplayTree()

        self.findBiggestSmaller(value)
        right = self.root.right
        self.root.right = None
        self.root.cnt -= numberOfElement(right)

        if right is not None:
            right.parent = None
            right.offset += self.root.offset
        return self, SplayTree(root=right)

    def indexList(self):
        self.checkValid()
        return self.root.indexList()

    def toList(self):
        """
        This will return a sorted list of all the element contains in self
        -------
        Returns:  A sorted list of all the element of self
        -------
        O(n)
        """
        self.checkValid()
        return self.root.sortedList()

    def merge(self, otherSplayTree):
        """
        This function return a splay tree containing all the element of self and otherSplayTree ,
        you need to have otherSplayTree.min()<self.max() or it will raise an Error, moreover after
        applied self isn't anymore a valid SplayTree and otherSplayTree point to the merged tree.
        -------
        Args:
            otherSplayTree : The otherSplayTree such that for otherSplayTree.max()<self.min()
        Returns:
            The merge splay tree containing all the element of self and otherSplayTree
        -------
        O(log(n)+log(m)) where n is the size of self and m is the size of otherSplayTree
        """
        self.checkValid()

        if self.isEmpty():
            self.valid = False
            return otherSplayTree

        if self.min() <= otherSplayTree.max():
            raise ValueError(
                "All the element of the left splay tree must be < than the element of self")
        self.valid = False
        assert otherSplayTree.root.right is None

        otherSplayTree.root.right = self.root
        self.root.parent = otherSplayTree.root
        otherSplayTree.root.cnt += self.root.cnt
        self.root.offset -= otherSplayTree.root.offset

        return otherSplayTree

    def height(self):
        """
        This function returns the height of self
        -------
        Returns: the height of self
        -------
        O(n) where n = self.size()
        """
        self.checkValid()
        return self.root.height()

    def isEmpty(self):
        """
        Returns: A boolean indicating if self is Empty
        -------
        O(1)
        """
        self.checkValid()
        return self.root.isEmpty()

    def _isBst(self):
        """
        Returns: A boolean indicating if self is a valid binary search tree
        -------
        O(n) where n =  self.size()
        """

        self.checkValid()
        try:
            self.root.isBst()
            return True
        except BaseException:
            return False

    def _isRoot(self):
        """
        Returns : A boolean indicating if self.root is really the root
        -------
        O(1)
        """
        self.checkValid()
        return self.root.isRoot()

    def getNode(self, value):
        """
        The node corresponding to value in the tree if value is in the tree otherwise it return None,
        note that node.value+node.getOffSet() == value and not node.value == value.
        ------
        Args: value
        Returns:  The node corresponding to value in the tree if value is in the tree otherwise
                it return None
        ------
        O(log(n))
        """
        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.find(value)
        self.changeRoot(node.splay())
        if self.root.value + self.root.offset != value:
            return None
        return self.root

    def insert(self, newValue):
        """
        Insert newValue inside a tree, return a boolean True if newValue was inserted
        and False if newValue was already present in self
        -------
        Returns: True if newValue was inserted and False if newValue was already present in the tree
        -------
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        isNew, node, _ = self.root.insert(newValue)
        self.changeRoot(node.splay())
        return isNew

    def delete(self, value):
        """
        Delete value inside self
        -------
        Args : value
        -------
        O(log(n)) where n = self.size()
        """
        self.checkValid()
        self.changeRoot(self.root.delete(value).splay())

    def min(self):
        """
        This method return the minimum of self
        -------
        Returns: the minimum of self if self isn't empty and None otherwise
        -------
        O(log(n)) where n = self.size()
        """
        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.min()
        self.changeRoot(node.splay())

        return self.root.value + self.root.offset

    def max(self):
        """
        This method return the maximum of self
        -------
        Returns: the maximum of self if self isn't empty and None otherwise
        -------
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.max()
        self.changeRoot(node.splay())
        return self.root.value + self.root.offset

    def find(self, value):
        """
        This return a boolean indicating if value is in self
        -------
        Args: value
        Returns: A boolean indicating if value is in self
        -------
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return False
        node, _ = self.root.find(value)
        self.changeRoot(node.splay())
        return self.root.value + self.root.offset == value

    def insertList(self, list):
        """
        This will insert the element of list inside self
        -------
        Args: list element to add
        -------
        O(m*log(n+m)) where n = self.size() and m = len(list)
        """

        self.checkValid()
        for e in list:
            self.insert(e)

    def findSmallestGreater(self, value):
        """
        This will return the smallest element of self >= value or None if it doesn't exist
        -------
        Args: value
        Returns: the smallest element of self >= value or None if it  doesn't exist
        -------
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.findSmallestGreater(value)
        self.changeRoot(node.splay())
        if self.root.value + self.root.offset < value:
            return None
        return self.root.value + self.root.offset

    def size(self):
        """
        Returns: Size of self
        -------
        O(1)
        """
        self.checkValid()
        return self.root.cnt

    def findBiggestSmaller(self, value):
        """
        This will return the smallest element of self <= value or None if it doesn't exist
        -------
        Args: value
        Returns: the biggest element of self <= value or None if doesn't exist
        -------
        O(log(n)) where n = self.size()
        """

        self.checkValid()
        if self.isEmpty():
            return None
        node, _ = self.root.findBiggestSmaller(value)
        self.changeRoot(node.splay())
        if self.root.value + self.root.offset > value:
            return None
        return self.root.value + self.root.offset
