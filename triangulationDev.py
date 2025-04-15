import random
import numpy as np
from CustomSwap import CustomSwap
from MapPermutation import MapPermutation
from LabelledMap import LabelledMap
from queue import deque
from MutableLabelledMap import MutableLabelledMap
from MutableTopologicalDemiEdge import MutableTopologicalDemiEdge
from RootedMap import RootedMap


def generateRandomBitstring(n):
    L = 4 * n - 2
    weight = n - 1
    bits = [0] * L
    one_positions = random.sample(range(L), weight)
    for pos in one_positions:
        bits[pos] = 1
    return bits


def checkPrefixCondition(bits):
    current_sum = 0
    for j in range(len(bits)-1):
        bit = bits[j]
        current_sum += 3 if bit == 1 else -1
        if current_sum <= -2:
            return False
    return True


def cyclicShift(bits, shift):
    return bits[shift:] + bits[:shift]


def generateValidCodeword(n):
    L = 4 * n - 2
    while True:
        b = generateRandomBitstring(n)
        candidate_shifts = [i for i in range(L)]
        valid_shifts = []
        for shift in candidate_shifts:
            if checkPrefixCondition(cyclicShift(b, shift)):
                valid_shifts.append(shift)
        if valid_shifts:
            chosen = random.choice(valid_shifts)
            return cyclicShift(b, chosen)


def fastGenerateValidCodeword(n):
    L = 4 * n - 2
    b = generateRandomBitstring(n)
    elems = [4 * bit - 1 for bit in b]          # maps 0 to -1 and 1 to 3

    q = deque()
    s = 0
    # our goal is to find the valid shifts of b, i. e. such that the value of every prefix is > -2 (and = -2 for the whole string)
    # note that here, the value of a string is 3 * (number of bits equal ro 1) - (number of bits equal to 0)

    # at each iteration i, s is the value of b from index 0 to i excluded
    # q is a deque of pairs (index, value) such that, at each iteration i:
    #   - if (j, v) is in q, the value of the bitstring from i to j included is v - s
    #   - the deque is ordered both by increasing indices and values (ie if (j, v) is before (j', v') in q, it holds that j<j' and v<=v')

    # note that we aren't interested in (j, v), (j', v') such that j < j' and v > v' because
    # the sum from i to j' will always be greater than the sum from i to j

    # at each iteration i, we just have to check that the value of the minimum, ie. the first element of the deque - s, is equal to -2
    # and that this minimum is reached at the end of the deque

    def add(j, v):
        while q and q[-1][1] > v:
            q.pop()
        q.append((j, v))

    # initialize the deque
    for i in range(L):
        s += elems[i]
        add(i, s)

    s = 0

    valid_shifts = []

    for i in range(L):
        if q[0][1] - s == -2 and q[0][0] == i + L - 1:
            valid_shifts.append(i)

        if q[0][0] == i:
            q.popleft()
        s += elems[i]

        add(i + L, -2 + s)

    assert len(valid_shifts) == 2

    return cyclicShift(b, random.choice(valid_shifts))


def getRandomRootedTwoLeafTree(n):
    b = fastGenerateValidCodeword(n)

    return rootedTwoLeafTreeFromBit(b)


def rootedTwoLeafTreeFromBit(b):
    n = (len(b)+2)//4

    alphaCycle = [(i, i+1) for i in range(1, 6*n-1, 2)]
    sigmaCycle = []

    p = [1]  # pile contenant les points fixes non couplés
    h = 0  # hauteur dans l'arbre

    H = [0]  # pile des hauteurs associés aux elements de p
    D = [0]*(6*n-1)  # hauteur de toutes les demi-edges
    Q = [[] for _ in range(n)]
    isFull = [False for _ in range(n)]  # Is the node full (in term of leaf)
    Q[0].append(2)
    D[1] = 0
    D[2] = 0
    Q.append([1])
    j = 4
    i = 0
    while i < len(b) and j <= 6 * n - 2:
        if b[i] == 1:
            D[j] = h+1
            D[j-1] = h
            Q[h].append(j-1)
            Q[h+1].append(j)
            j += 2
            h += 1
            i += 1
            continue
        # b[i] == 0
        if H and H[-1] == h:
            D[j] = h
            D[j-1] = h
            Q[h].append(j-1)
            Q.append([j])
            p.pop()
            H.pop()

            isFull[h] = True

            j += 2
            i += 1
            continue

        if isFull[h]:
            if Q[h]:
                sigmaCycle.append(tuple(Q[h]))
                Q[h] = []
            isFull[h] = False
            i += 1
            h -= 1
            continue

        D[j] = h
        D[j-1] = h

        Q[h].append(j-1)
        p.append(j)
        Q.append([j])
        H.append(h)

        j += 2
        i += 1

    for level in range(len(Q)):
        if Q[level]:
            sigmaCycle.append(tuple(Q[level]))

    sigma = MapPermutation(sigmaCycle)
    alpha = MapPermutation(alphaCycle)

    tree = RootedMap(sigma=sigma, alpha=alpha)

    return tree


def treeToTriangulation(tree):

    def isOnInnerEdge(Z):
        return Z.n != Z and (Z.c).n != Z.c

    def isOnLeafEdge(Z):
        return not isOnInnerEdge(Z)

    def extremeOnEdgeAfter(Z):
        return ((Z.c).n).c

    def isSpecial(Z):
        return isOnLeafEdge(extremeOnEdgeAfter(Z)) and isOnLeafEdge(Z)

    triangulation = MutableLabelledMap(lmap=tree)

    root = triangulation.X(1)

    outerList = []
    isLeaf = np.zeros(2*triangulation.m+1)

    for A in root.face():
        if A.n == A:
            outerList.append(A)
            isLeaf[A.raw] = True

    cntSpecial = 0
    for A in outerList:
        cntSpecial += isSpecial(A)

    A = triangulation.X(1)
    B, C = A.f, (A.f).f

    while cntSpecial > 2:
        if isOnInnerEdge(A) and isOnInnerEdge(B) and isOnLeafEdge(C):

            cntSpecial -= isSpecial(C.c)

            N, _ = A.link(C.c)
            N.contract()
            A = (C.c).pf
            B, C = A.f, (A.f).f
        else:
            A = A.f
            B = B.f
            C = C.f

    specialList = []
    for A in outerList:
        if isSpecial(A):
            specialList.append((A, extremeOnEdgeAfter(A)))

    A, B = specialList[0]
    C, D = specialList[1]

    def closeOp(U, V):
        W = U.addEdgeAfter()
        W.link(V)[0].contract()
        W.contract()
        return V

    def enclose(A, D):
        V = A.c.pf
        U = A
        while True:
            nxt = V.pf
            if V.n == V:
                U = closeOp(U, V)
            if V == D:
                break
            V = nxt
    enclose(A, D)
    enclose(C, B)

    X, Y = A.link(C)

    rng = random.Random()
    tau = CustomSwap([(X.raw, 1)])

    if rng.random() < 0.5:
        tau = CustomSwap([(Y.raw, 1)])

    triangulation = triangulation.relabel(tau)
    return RootedMap(triangulation)


def getRandomTriangulation(n):
    return treeToTriangulation(getRandomRootedTwoLeafTree(n))
