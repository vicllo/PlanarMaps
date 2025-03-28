from sage.all import Permutation

from LabelledMap import LabelledMap


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

    def copy(self):
        return MutableLabelledMap(self.sigma, self.alpha)

    def _buildTransp(self, swaps):
        """
        Builds a transposition from a list of couples representing swaps (possibly null).
        For example, _buildTransp([(1,1),(2,4),(3,3)]) returns (1,4,3,2).
        """

        return Permutation(list(filter(lambda t: t[0] != t[1], swaps)))

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

        # swaps the iEdge and its dual with the half-edges 2*self.m-1 and 2*self.m to allow easily removing them

        swap = self._buildTransp([(iEdge, 2*self.m-1), (self.alpha(iEdge), 2*self.m)])

        self.alpha = swap * self.alpha * swap
        self.sigma = swap * self.sigma * swap

        # merge the two neighbors lists

        swp = self._buildTransp([(2*self.m-1, self.sigma(2*self.m)), (2*self.m, self.sigma(2*self.m-1))])

        self.sigma = self.sigma * swp

        # now sigma is correct; we just have to delete the leftover (2*m-1, 2*m) transposition

        self.sigma = Permutation(self.sigma.to_cycles()[:-1])
        self.alpha = Permutation(self.alpha.to_cycles()[:-1])

        self._updateAttributes()

    def addEdge(self, iEdge1, iEdge2, keepGenus=True):
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

    def deleteEdge(self, iEdge):
        """ Delete the given half-edge

        INPUT:

        - ``iEdge`` integer; the index of the half-edge to delete

        TEST::
            sage: sigma = Permutation([(8,1),(2,3),(4,5),(6,7)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6),(7,8)])
            sage: P = MutableLabelledMap(sigma,alpha)
            sage: P.deleteEdge(8)
            sage: P.sigma.to_cycles()
            [(1,), (2, 3), (4, 5), (6,)]
            sage: P.alpha.to_cycles()
            [(1, 2), (3, 4), (5, 6)]
        """
        if iEdge < 1 or iEdge > 2 * self.m:
            raise ValueError("Invalid half-edge number.")
        sigma = self.sigma
        alpha = self.alpha
        # Swap iEdge with 2*m if it is not already 2*m
        if iEdge != 2*self.m :
            swap1 = Permutation([(iEdge,2*self.m)])
            sigma = swap1 * self.sigma * swap1
            alpha = swap1 * self.alpha * swap1
        # Swap alpha(iEdge) with 2*m-1 if it is not already 2*m-1
        if self.alpha(iEdge) != 2*self.m-1:
            if self.alpha(iEdge) != 2*self.m :
                swap2 = Permutation([(self.alpha(iEdge),2*self.m-1),(2*self.m,)])
            else:
                swap2 = Permutation([(self.alpha(iEdge),2*self.m-1)])
            sigma = swap2 * sigma * swap2
            alpha = swap2 * alpha * swap2

        inverseSigma = sigma.inverse()
        #Apply two permutations to sigma to update its structure
        if 2*self.m-1 != sigma(2*self.m-1):
            t1 = Permutation ([(2*self.m-1 , inverseSigma(2*self.m-1) ),(2*self.m,)])
            sigma = sigma.left_action_product(t1)
        if 2*self.m!=sigma(2*self.m):
            t2 = Permutation ([(2*self.m , inverseSigma(2*self.m) )])
            sigma = sigma.left_action_product(t2)

        # Construct the new alpha and sigma
        new_domain = list(range(1, 2*self.m-1))
        new_sigma = Permutation([sigma(i) for i in new_domain])
        new_alpha = Permutation(alpha.to_cycles()[:-1])

        #test if the new Graph is connected
        seen = [False] * (self.size - 1)
        seen[0] = seen[1] = True            # half-edges are numbered from 1 to size, included
        todo = [1]
        while todo:
            i = todo.pop()
            if not seen[new_alpha(i)]:
                todo.append(new_alpha(i))
                seen[new_alpha(i)] = True
            if not seen[new_sigma(i)]:
                todo.append(new_sigma(i))
                seen[new_sigma(i)] = True

        if False in seen:
            raise ValueError("The graph is not connected")
        self.alpha = new_alpha
        self.sigma = new_sigma
        self._updateAttributes()

    def deleteVertex(self,iEdge):
        sigma = self.sigma
        alpha = self.alpha
        m=self.m

        def transposition(i,j):
            if i!=2*m and j!=2*m:
                return Permutation([(i,j),(2*m,)])
            else :
                return Permutation([(min(i,j),2*m)])
        for cycle in sigma.to_cycles():
            if iEdge in cycle:
                Cycle = cycle
                break
        for i in range(len(Cycle)):
            # swap (Cycle[i],alpha(Cycle[i)) --> (2m-2i,2m-2i-1)
            swap1=Permutation([i for i in range(1, 2*m+1)])
            swap2=Permutation([i for i in range(1, 2*m+1)])
            if Cycle[i]!=2*m-2*i:
                swap1 = transposition(Cycle[i],2*m-2*i)
            if alpha(Cycle[i])!=2*m-2*i-1:
                swap2 = transposition(alpha(Cycle[i]),2*m-2*i-1)
            sigma = swap1 * sigma * swap1
            alpha = swap1 * alpha * swap1
            sigma = swap2 * sigma * swap2
            alpha = swap2 * alpha * swap2

        self.sigma=sigma
        self.alpha=alpha

        for i in range(len(Cycle)):
            self.deleteEdge(2*m-2*i)
        self._updateAttributes()

    def contractFace(self, iEdge):
        """
        Contracts the face corresponding to the given half-edge into a single node.
        Note that half-edges go clockwise around a face.
        """

        # first, we loop through the face, and swap all its half-edges (and their counterparts) with the last half-edges
        # this way, after the contraction, we will just have to remove the last cycles of sigma

        faceEdge = iEdge
        index = 2*self.m

        sigma = self.sigma
        alpha = self.alpha

        while faceEdge != 2*self.m or index == 2*self.m:
            swp = Permutation([(faceEdge, index)]) if faceEdge != index else Permutation([])
            print ("swapping", faceEdge, index)
            if index % 2 == 0:
                faceEdge = alpha(faceEdge)
            else:
                faceEdge = sigma(faceEdge)

            sigma = swp * sigma * swp
            alpha = swp * alpha * swp

            index -= 1

        iEdge = 2 * self.m

        self.sigma = sigma
        self.alpha = alpha

        #return

        sigmaInv = sigma.inverse()
        swaps = []

        faceEdge = sigmaInv(iEdge)     # since sigma goes anticlockwise around nodes, we need to go through the face in the same order
        outgoingEdge = -1              # we don't know any outgoing edge yet
        firstOutgoing = -1

        # we loop through the face and maps each outgoing edge to the following one in the new sigma

        firstIter = True        # we need to actually enter the while loop... and stop when we have been through the whole face

        while firstIter or sigma(faceEdge) != iEdge:
            if sigma(sigma(faceEdge)) != faceEdge:          # there is at least one outgoing edge
                currentEdge = sigma(sigma(faceEdge))
                if outgoingEdge != -1:                                  # if we've seen another outgoing edge, we map it to the current one
                    swaps.append((outgoingEdge, sigma(faceEdge)))  # after the swap, sigma will be outgoingEdge (-> sigma(faceEdge)) -> currentEdge
                else:
                    firstOutgoing = currentEdge
                outgoingEdge = currentEdge

            faceEdge = sigmaInv(alpha(faceEdge))

            print ("faceEdge", faceEdge, "outgoingEdge", outgoingEdge, "firstOutgoing", firstOutgoing)
            firstIter = False

        # we also need to map the last outgoing edge to the first one we have seen
        if firstOutgoing != -1:
            swaps.append((outgoingEdge, sigmaInv(firstOutgoing)))

        sigma = self._buildTransp(swaps) * sigma

        # then, the half-edges we need to remove are from index + 1 to 2*m

        #def rem_

        self.sigma = sigma
        self.alpha = alpha

        self._updateAttributes()
