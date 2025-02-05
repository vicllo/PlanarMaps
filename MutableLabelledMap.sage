load("LabelledMap.sage")

class MutableLabelledMap(LabelledMap):
    def _updateAttributes(self):
        """
        Recompute the attributes m, n, f and phi from sigma and alpha.
        """

        self.phi = self.alpha.right_action_product(self.sigma)
        self.size = self.sigma.size()
        self.m = self.size // 2 
        self.f = len(self.phi.to_cycles())
        self.n = len(self.sigma.to_cycles())

    def contractEdge(self, iEdge):
        """
        Contracts the given half-edge (i. e. merge the two nodes that it links, and removes the edge itself).

        INPUT:

        - ``iEdge`` integer; the index of the half-edge to contract

        TESTS::
            sage: adj = [(5,4,2),(1,3,6),(4,7,2),(8,3,1),(8,1,6),(5,2,7),(3,8,6),(7,4,5)]
            sage: cube = MutableLabelledMap(adj = adj)
            sage: cube.contractEdge(1)
            sage: cube.contractEdge(3)
            sage: cube.numberOfEdges()
            10
            sage: cube.numberOfNodes()
            6
        """

        if iEdge < 1 or iEdge > 2 * self.m:
            raise ValueError("Invalid half-edge number.")

        def buildTransp(l):            # build a permutation from a list of possibly null transpositions
            return Permutation(list(filter(lambda t: t[0] != t[1], l)))        # ie. buildTransp([(1,1),(2,4),(3,3)]) = (1,4,3,2)
        
        # swaps the iEdge and its dual with the half-edges 2*self.m-1 and 2*self.m to allow easily removing them

        swap = buildTransp([(iEdge, 2*self.m-1), (self.alpha(iEdge), 2*self.m)])

        self.alpha = swap * self.alpha * swap
        self.sigma = swap * self.sigma * swap

        # merge the two neighbors lists

        swp = buildTransp([(2*self.m-1, self.sigma(2*self.m)), (2*self.m, self.sigma(2*self.m-1))])
        
        self.sigma = self.sigma * swp

        # now sigma is correct; we just have to delete the leftover (2*m-1, 2*m) transposition

        self.sigma = Permutation(self.sigma.to_cycles()[:-1])
        self.alpha = Permutation(self.alpha.to_cycles()[:-1])

        self._updateAttributes()

    def addEdge(self, iEdge1, iEdge2, keepGenus = True):
        """
        Add an edge between the two nodes from whom start iEdge1 and iEdge2.
        The resulting half-edges will be immediately before iEdge1 and iEdge2, respectively, in the new permutation sigma.
        If keepGenus is True and this addition would change the map's genus (i.e. the two half-edges are part of the same face), an error will be raised.
                Initializes the labelled map from either the permutations alpha and sigma, or an adjacency list (giving for
                vertex a list of its neighbors in order; vertices must be numbered from 1 to n).

        INPUT:

        - ``iEdge1`` -- int; Index of the first half-edge before which the edge should be added
        - ``iEdge2`` -- int; Index of the second half-edge before which the edge should be added
        - ``keepGenus`` -- bool; indicates whether an error should be raised if the new edge increases the map's genus.

        EXAMPLES:
        sage: m = MutableLabelledMap(adj = ([2,4],[1,3],[2,4],[1,3]))
        sage: m.addEdge(1, 6)
        sage: m.addEdge(2, 5)
        sage: m
        Labelled map | Sigma : [11, 9, 4, 3, 10, 12, 8, 7, 1, 6, 2, 5], Alpha : [3, 7, 1, 5, 4, 8, 2, 6, 10, 9, 12, 11]

        TESTS::
            sage: m = MutableLabelledMap(adj = ([2,4],[1,3],[2,4],[1,3]))
            sage: m.addEdge(1, 6)
            sage: m.addEdge(4, 7)
            Traceback (most recent call last):
            ...
            ValueError: Adding an edge between 4 and 7 would increase the map's genus
            sage: m.addEdge(4, 7, keepGenus = False)
            sage: m.genus()
            1
        """

        if keepGenus and iEdge1 != iEdge2:
            i = self.sigma(self.alpha(iEdge1))
            while i != iEdge1 and i != iEdge2:
                i = self.sigma(self.alpha(i))

            if i == iEdge1:        # we went through the whole face without encountering iEdge2
                raise ValueError(f"Adding an edge between {iEdge1} and {iEdge2} would increase the map's genus")
        
        h1, h2 = 2*self.m+1, 2*self.m+2

        self.alpha = Permutation([(h1, h2)]) * self.alpha        # the new half-edges are part of the same edge
        if iEdge1 != iEdge2:
            self.sigma = Permutation([(self.sigma.inverse()(iEdge1), h1), (self.sigma.inverse()(iEdge2), h2)]) * self.sigma
            # ...->iEdge1->j->... is now ...->h1->iEdge1->j->... in sigma and same for h2/iEdge2
        else:
            self.sigma = Permutation([(self.sigma.inverse()(iEdge1), h1, h2)]) * self.sigma
        
        self._updateAttributes()