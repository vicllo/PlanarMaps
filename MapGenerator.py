from collections import deque
import random
from MapPermutation import MapPermutation

from RootedMap import RootedMap


class MapGenerator:
    """
    This class represents an abstraction containing
    methods to generate a Map.
    """

    def __init__(self):
        # Set it to true when in production
        # during debugging to False
        self._production = True

    def cube(self):
        """Returns the standard cube map."""
        return RootedMap(
            adj=[
                (5, 4, 2),
                (1, 3, 6),
                (4, 7, 2),
                (8, 3, 1),
                (8, 1, 6),
                (5, 2, 7),
                (3, 8, 6),
                (7, 4, 5),
            ], trust=self._production,
        )

    def complete_map(self, n):
        """
        Returns an arbitrary rooted map corresponding to the complete
        graph with n nodes. The genus is guaranteed to be zero if the
        graph is planar (i.e., n <= 4).
        """
        adj = list(
            tuple((j + i) % n + 1 for j in range(1, n)) for i in range(n)
        )
        m = RootedMap(adj=adj, trust=self._production,)
        if n <= 4:
            m = m.force_planar()
        return m

    def getRandomDyckPath(self, n, seed=None):
        """
        Returns a random Dyck path of size n (uniform random generation).

        INPUT:
            - ``n`` -- int; size of the path
            - ``seed`` -- int; A random seed;
            if None is used, no random seed will be set.

        OUTPUT:
            A list of size 2*n with +1 for up and
            -1 for down steps in the Dyck path.

        EXAMPLE::
            sage: dyckPath = MapGenerator().getRandomDyckPath(10)
            sage: dyckPath
            [1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1,
            1, 1, -1, -1, 1, -1, 1, -1]

        TESTS::
            sage: dyckPath = MapGenerator().getRandomDyckPath(50)
            sage: level = 0
            sage: for step in dyckPath:
            ....:     level += step
            ....:     assert level >= 0
            ....:
            sage: assert level == 0
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))
        N = 2 * n + 1
        dyck = [1] * n + [-1] * (n + 1)
        rng.shuffle(dyck)
        level = 0
        minlevel = 0
        posmin = 0
        for i in range(N):
            level += dyck[i]
            if level < minlevel:
                posmin = i + 1
                minlevel = level
        Dyckfinal = dyck[posmin:] + dyck[:posmin]
        return Dyckfinal[:-1]

    def getRandomPermutation(self, n, seed=None):
        """
        Returns a random permutation of size n.

        Args:
            n : The size of the permutation.
            seed : A random seed; if None is used,
                   no random seed will be set.

        Returns:
            A random permutation of size n.
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))
        lst = [i + 1 for i in range(n)]
        rng.shuffle(lst)
        return MapPermutation(lst)

    def isValidDyckPath(self, dyckPathCandidate):
        """
        Checks whether the given Dyck path candidate is valid.

        Args:
            dyckPathCandidate : A list representing a potential Dyck path.

        Returns:
            A boolean indicating whether or not dyckPathCandidate is a
            correct Dyck path.
        """
        if len(dyckPathCandidate) == 0 or len(dyckPathCandidate) % 2 == 1:
            return False

        for e in dyckPathCandidate:
            if e != -1 and e != 1:
                return False

        S = 0
        for e in dyckPathCandidate:
            S += e
            if S < 0:
                return False
        return S == 0

    def getTreeFromDyckPath(self, dyckPath, trust=False):
        """
        Given a Dyck path, this function returns the associated rooted tree.

        Args:
            dyckPath : A list representing a Dyck path, with +1 for up and
                       -1 for down.
            trust: A boolean indicating whether to trust that we have a dyckPath 
        Returns:
            The corresponding rooted plane tree if dyckPath is valid;
            otherwise, raises an error.

        Complexity:
            O(k), where k = len(dyckPath)
        """
        if not trust and not self.isValidDyckPath(dyckPath):
            raise ValueError("The given list isn't a Dyck path")

        phiCycle = []
        alphaCycle = []
        p = []

        for i in range(len(dyckPath)):
            phiCycle.append(i + 1)
            if dyckPath[i] < 0:
                otherDemiEdge = p.pop()
                alphaCycle.append((i + 1, otherDemiEdge))
            else:
                p.append(i + 1)

        phi = MapPermutation([tuple(phiCycle)])
        alpha = MapPermutation(alphaCycle)
        sigma = phi.left_action_product(alpha)

        return RootedMap(alpha=alpha, sigma=sigma, trust=self._production, )

    def getRandomLabellingTree(self, tree, seed=None):
        """
        Generates a uniformly random labelling of a tree.

        Args:
            tree : The input rooted tree.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A list of size 2*tree.m + 1 where labelling[i] (for i >= 1)
            represents the label of demi-edge i. The first value
            (labelling[0]) is set to -1 but has no meaning.

        Complexity:
            O(m), where m is the number of edges in the tree.
        """
        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        sigma = tree.sigma
        alpha = tree.alpha

        p = deque()
        nodes = tree.nodes()

        demiEdgeToNodeId = [-1 for _ in range(2 * tree.m + 1)]

        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                demiEdgeToNodeId[nodes[i][j]] = i

        labelling = [-1 for _ in range(2 * tree.m + 1)]
        startNodeId = demiEdgeToNodeId[1]

        p.append(startNodeId)
        labellingNodes = [-1 for _ in range(len(nodes))]
        seen = [False for _ in range(len(nodes))]

        seen[startNodeId] = True
        transition = [-1, 1, 0]
        labellingNodes[startNodeId] = 0

        while len(p) > 0:
            nodeId = p.popleft()

            for demiEdge in nodes[nodeId]:
                labelling[demiEdge] = labellingNodes[nodeId]

                alphaDemiEdge = alpha(demiEdge)
                alphaNodeId = demiEdgeToNodeId[alphaDemiEdge]
                if not seen[alphaNodeId]:
                    dLabel = rng.sample(transition, 1)[0]
                    labellingNodes[alphaNodeId] = (
                        dLabel + labellingNodes[nodeId]
                    )
                    seen[alphaNodeId] = True
                    p.append(alphaNodeId)

        return labelling

    def getRandomTree(self, numberOfEdge, seed=None):
        """
        Generates a uniformly random rooted tree.

        Args:
            numberOfEdge : The number of edges in the tree.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A randomly selected rooted tree with numberOfEdge edges.

        Complexity:
            O(numberOfEdge)
        """
        return self.getTreeFromDyckPath(
            self.getRandomDyckPath(numberOfEdge, seed=seed), trust=self._production
        )

    def getRandomLabelledTree(self, numberOfEdge, seed=None):
        """
        Generates a uniformly random rooted tree along with a labelling.

        Args:
            numberOfEdge : The number of edges in the tree.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A tuple (tree, labelling) where:
            - tree : A randomly selected rooted tree with numberOfEdge edges.
            - labelling : A list of labels for the tree’s demi-edges.

        Complexity:
            O(numberOfEdge)
        """
        tree = self.getRandomTree(numberOfEdge, seed=seed)
        return tree, self.getRandomLabellingTree(tree, seed=seed)

    def getRandomPlanarQuadrangulation(self, numberOfFace, seed=None):
        """
        Generates a uniformly random rooted planar quadrangulation with a
        specified number of faces.

        Args:
            numberOfFace : The number of faces in the quadrangulation.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A randomly selected rooted planar quadrangulation with
            numberOfFace faces.

        Complexity:
            O(numberOfFace)
        """
        tree, labelling = self.getRandomLabelledTree(numberOfFace, seed=seed)
        quadA, quadB = tree.inverseShaefferTree(
            returnMarkedDemiEdge=False, labelled=labelling
        )

        rng = random.Random()
        if seed is not None:
            rng.seed(int(seed))

        if rng.random() < 0.5:
            return quadB
        return quadA

    def getRandomPlanarMap(self, numberOfEdge, seed=None):
        """
        Generates a uniformly random rooted planar map with a specified
        number of edges.

        Args:
            numberOfEdge : The number of edges in the rooted map.
            seed : A random seed; if None is used, no random seed will be set.

        Returns:
            A randomly selected rooted planar map with numberOfEdge edges.

        Complexity:
            O(numberOfEdge)
        """
        quad = self.getRandomPlanarQuadrangulation(numberOfEdge, seed=seed)

        return quad.inverseQuadrangulation()
