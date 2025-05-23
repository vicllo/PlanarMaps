
from sage.all import Permutation  # Import sage library
from labelled_map import LabelledMap
from map_error import NotImplementedError
from primitive_mutable_topological_demi_edge import (
    PrimitiveMutableTopologicalDemiEdge,
)
from primitive_rotating_permutation import PrimitiveRotatingPermutation
from primitive_rotating_permutation_utils_abstractor import (
    PrimitiveRotatingPermutationUtilsAbstractor,
)


class PrimitiveMutableLabelledMap(LabelledMap):
    """
    This class represent a more primitive version of MutableLabelledMap
    It implements only the following operations for now ( addAfter, addBefore,deleteEdge,contractEdge)
    and place more responsibility on the user and remove some of the other methods( such as checking if two demi edge are on the same node or face),
    the only advantage is the fact that the operations listed are in O(1) instead of O(log(m)).

    Attributes
    ----------
    sigma :  Permutation or MapPermutation
        Permutation that maps a half-edge to the half-edge incident to
        it in a anti-clockwise direction around the vertex it belongs to.

    alpha : Permutation or MapPermutation
        Fixed-point free involution whose cycles are given by the edges.

    phi: Permutation or MapPermutation
        Permutation that maps a half-edges to the half-edge next to it
        in his face in the clockwise orde.

    m: The number of edge of the map

    g: The genus of the map

    f: The number of face of the map

    q: The number of demi edges of the map

    """

    def __init__(
        self,
        sigma: Permutation | None = None,
        alpha: Permutation | None = None,
        adj: list[tuple[int, ...]] | None = None,
        lmap: LabelledMap | None = None
    ):
        r"""
        Init the PrimitiveMutableLabelledMap

        INPUT:

        - ``sigma`` -- Permutation ; Permutation that maps a half-edge
          to the half-edge incident to it in anti-clockwise direction around
          the vertex it belongs to.
        - ``alpha`` -- Permutation ; Permutation that maps a half-edge
            Fixed-point free involution whose cycles are given by the edges.
        - ``ajd``-- and adjacency list be careful the order of the
            node in your adjaceny will be used to choose the embedding
        - ``lmap`` a LabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)

        NOTE:

            O(m) where m is the size of the map
        """

        if isinstance(lmap, LabelledMap):
            self.__init__(lmap.sigma, lmap.alpha)
            return

        super().__init__(sigma, alpha, adj)
        self.sigma = PrimitiveRotatingPermutation(self.sigma)
        self.alpha = PrimitiveRotatingPermutation(self.alpha)
        self.phi = PrimitiveRotatingPermutation(self.phi)
        self.phiUtilsAbstractor = PrimitiveRotatingPermutationUtilsAbstractor(
            self.phi)
        self.sigmaUtilsAbstractor = PrimitiveRotatingPermutationUtilsAbstractor(
            self.sigma)

        for e in range(1, self.q + 1):
            self.topologicalMap[e] = PrimitiveMutableTopologicalDemiEdge(
                self, e)

    def X(self, demiEdge: int) -> PrimitiveMutableTopologicalDemiEdge:
        """
        Return the PrimitiveMutableTopologicalDemiEdge associated to demiEdge.

        INPUT:

            - ``demiEdge`` -- int ; an index associated to a demiEdge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = MutableLabelledMap(alpha = alpha,sigma=sigma)
            sage: m.X(1)
            X(1)

        NOTE:

            Complexity is O(1)
        """
        return self.getTopologicalDemiEdge(demiEdge)

    def getTopologicalDemiEdge(self, demiEdge: int) -> PrimitiveMutableTopologicalDemiEdge:
        """
        The PrimitiveMutableTopologicalDemiEdge associated to demiEdge

        INPUT:

            - ``demiEdge`` -- int ; An index associated to a demiEdge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: m = MutableLabelledMap(alpha = alpha,sigma=sigma)
            sage: m.getTopologicalDemiEdge(1)
            X(1)

        NOTE:

            Complexity is O(1)
        """

        return self.topologicalMap[demiEdge]

    def addEdge(self, startDemiEdge: int, endDemiEdge: int) -> tuple[PrimitiveMutableTopologicalDemiEdge, PrimitiveMutableTopologicalDemiEdge]:
        """
        This will add an edge between the node of startDemiEdge to endDemiEdge, the edge will be added as follow ,
        let denote (A,B) the demi edges composing this new edge A will be on the same node as startDemiEdge but before it
        and B on the same node as endDemiEdge but before it.
        It will return two MutableLabelledMap topoDemiEdgeA,topoDemiEdgeB corresponding to the new demi edge A and B

        WARNING: Nothing is guaranteed if startDemiEdge and endDemiEdge are not on the same face

        INPUT:

        - ``startDemiEdge`` -- int ; A demi edge of self
        - ``endDemiEdge`` -- int : A demi edge of self

        OUTPUT:

            topoDemiEdgeA,topoDemiEdgeB as described above

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.addEdge(3,7)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 4, 21, 3)]

        NOTE:

            O(1)
        """

        # This will go before endDemiEdge
        newIndexEnd = self.q + 1

        # This will go before startDemiEdge
        newIndexStart = self.q + 2

        # We add the new edge to alpha
        self._addEdgeToAlpha()

        # We add the new edge to sigma
        # The order must be coherent with the choice made
        # newIndexStart,newIndexEnd
        self.sigma.addBefore(endDemiEdge)
        self.sigma.addBefore(startDemiEdge)

        # We add  the new demi edge to phi
        self.phi.cutAdd(startDemiEdge, endDemiEdge, newIndexStart, newIndexEnd)

        self.topologicalMap[newIndexStart] = PrimitiveMutableTopologicalDemiEdge(
            self, newIndexStart)
        self.topologicalMap[newIndexEnd] = PrimitiveMutableTopologicalDemiEdge(
            self, newIndexEnd)

        return self.X(newIndexStart), self.X(newIndexEnd)

    def _swapTopologicalDemiEdgeValue(self, demiEdge: int, otherDemiEdge: int) -> None:
        """
        This will change the index of the topologica demi edge associate to the demi edge of label demi edge to otherDemiEdge
        and same thing but swapping the role or demiEdge and otherDemiEdge

        INPUT:

        - ``demiEdge``  -- int
        - ``otherDemiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(1)
            sage: B = mm.X(2)
            sage: A,B
            (X(1), X(2))
            sage: mm._swapTopologicalDemiEdgeValue(1,2)
            sage: A,B
            (X(2), X(1))

        NOTE:

            O(1)
        """

        topoDemiEdge = self.X(demiEdge)
        otherTopoDemiEdge = self.X(otherDemiEdge)
        topoDemiEdge._swapIndex(otherTopoDemiEdge)
        self.topologicalMap[demiEdge] = otherTopoDemiEdge
        self.topologicalMap[otherDemiEdge] = topoDemiEdge

    def simpleSwap(self, demiEdge: int, otherDemiEdge: int) -> None:
        """
        This function relabel demiEdge into otherDemiEdge and otherDemiEdge into demiEdge

        INPUT:

        - ``demiEdge`` -- int
        - ``otherDemiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(1)
            sage: B = mm.X(2)
            sage: A,B
            (X(1), X(2))
            sage: mm.simpleSwap(2,1)
            sage: A,B
            (X(2), X(1))

        NOTE:

            O(1))
        """

        self._swapTopologicalDemiEdgeValue(demiEdge, otherDemiEdge)

        self.phi.swapIndex(demiEdge, otherDemiEdge)
        self.alpha.swapIndex(demiEdge, otherDemiEdge)
        self.sigma.swapIndex(demiEdge, otherDemiEdge)

    def labelToTheEnd(self, listIndexes: list[int]) -> None:
        """
        This function relabel the indexes in listIndexes to take
        if there is say k indexes in listIndexes (and they are valid) it
        will relabel them into value in m,...,m-k+1

        INPUT:

        - ``listIndexes`` -- List[int] ;A list of valid demi edge in self( otherwise it will raise an error)

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.labelToTheEnd([7,5,3])
            sage: mm.faces()
            [(1, 18, 2, 19, 4, 20, 11, 16, 3, 13, 12, 15, 14, 5, 7, 17, 9, 8, 10, 6)]

        NOTE:

            O(k) where k = len(listIndexes)
        """

        for i in listIndexes:
            if i != int(i) or i > self.q or i <= 0:
                raise ValueError("The list of indexes isn't correct")

        corres = self.phi.labelToTheEnd(listIndexes)
        self.alpha.labelToTheEnd(listIndexes)
        self.sigma.labelToTheEnd(listIndexes)
        for index in listIndexes:
            if index in corres:
                self._swapTopologicalDemiEdgeValue(index, corres[index])

    def _BruteDeleteEdge(self, demiEdge: int, sameFace: int) -> None:
        """
        INPUT:

        - ``demiEdge`` -- int; demi edge on  the edge to delete

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.addEdge(3,11)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 4, 7, 21, 3)]
            sage: mm.deleteEdge(22,sameFace = False)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]

        NOTE:

            O(1)
        """

        otherDemiEdge = self.alpha(demiEdge)

        if demiEdge != self.q or otherDemiEdge != self.q-1:
            self.labelToTheEnd([demiEdge, otherDemiEdge])
            self._BruteDeleteEdge(self.q, sameFace=sameFace)
            return

        # TopologicalDemiEdge thing
        topoDemiEdge = self.X(demiEdge)
        otherTopoDemiEdge = topoDemiEdge.c
        self._removeTopologicalDemiEdge(topoDemiEdge)
        self._removeTopologicalDemiEdge(otherTopoDemiEdge)

        # Actual removing
        self.sigma.deleteLastKIndex(2)
        self.alpha.deleteLastKIndex(2)
        if not sameFace:
            self.phi.mergeDelete(demiEdge, otherDemiEdge)
        else:
            # Because demiEdge doesn't break the connectivity
            # demiEdge and otherDemiEdge are on the same face
            # And the face will be cut
            # Look Strange  cause it mean that deleting a edge will increase the
            # Number of face but not so much because it won't increase the number of
            # On the current surface just the map obteined has more face on
            # The new lower genus surface
            # Note that this case only happen in genus > 0
            self.phi.cutDelete(demiEdge, otherDemiEdge)

    def deleteEdge(self, demiEdge: int, sameFace: bool) -> None:
        """
        INPUT:

        - ``demiEdge`` -- int; a demi edge corresponding to the edge to delete
        - ``sameFace``-- bool ;  a boolean indicating if demiEdge and alpha(demiEdge) are on the same face or not

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.addEdge(3,11)
            (X(22), X(21))
            sage: mm.faces()
            [(1, 22, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6),
             (2, 5, 4, 7, 21, 3)]
            sage: mm.deleteEdge(22, False)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]

        NOTE:

            O(1)
        """

        self._BruteDeleteEdge(demiEdge, sameFace=sameFace)

    def _addEdgeToAlpha(self) -> None:
        """
        Add one edge compose of demi edges self.size() and self.size()+1 to alpha toward the end

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.alpha.pretty_print()
            Primitive Rotating permutation: [(1, 3), (2, 5), (4, 6), (7, 9), (8, 10), (11, 13), (12, 15), (14, 17), (16, 18), (19, 20)]
            sage: mm._addEdgeToAlpha()
            sage: mm.alpha.pretty_print()
            Primitive Rotating permutation: [(1, 3), (2, 5), (4, 6), (7, 9), (8, 10), (11, 13), (12, 15), (14, 17), (16, 18), (19, 20), (21, 22)]

        NOTE:

            O(1)
        """
        self.alpha.stretch(2)
        self.alpha.addAfterGeneral(self.alpha.size(), self.alpha.size() - 1)

    def _addTopologicalDemiEdge(self, demiEdge: int) -> None:
        """
        Add a new MutableTopologicalDemiEdge to self with index associated to demiEdge

        INPUT:

        - ``demiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.X(42)
            ....: except:
            ....:    print("NOT TOPO")
            ....:
            NOT TOPO
            sage: mm._addTopologicalDemiEdge(42)
            sage: mm.X(42)
            X(42)

        NOTE:

            O(1)
        """
        self.topologicalMap[demiEdge] = PrimitiveMutableTopologicalDemiEdge(
            self, demiEdge)

    def _removeTopologicalDemiEdge(self, topoDemiEdge: PrimitiveMutableTopologicalDemiEdge) -> None:
        """
        Remove and make invalid topoDemiEdge

        INPUT:

        - ``topoDemiEdge`` -- PrimitiveMutableTopologicalDemiEdge

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(10)
            sage: mm._removeTopologicalDemiEdge(A)
            sage: try:
            ....:    mm.X(10)
            ....: except:
            ....:    print("NOT TOPO")
            ....:
            NOT TOPO

        NOTE:

            O(1)
        """

        self.topologicalMap.pop(topoDemiEdge.raw)
        topoDemiEdge._invalidate()

    def addEdgeAfter(self, demiEdge: int) -> PrimitiveMutableTopologicalDemiEdge:
        """
        This function will create an new edge attached to the same node as demi edge such that it is after
        demiEdge in the trigonometric order.

        INPUT:

        - ``demiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(10)
            sage: A.n,A.pn
            (X(10), X(10))
            sage: mm.addEdgeAfter(10)
            X(22)
            sage: mm.addEdgeBefore(10)
            X(24)
            sage: A.n,A.pn
            (X(21), X(23))

        NOTE:

            O(1)
        """

        newDemiEdge = self.q + 1
        otherNewDemiEdge = self.q + 2

        # TopologicalDemiEdge things
        self._addTopologicalDemiEdge(newDemiEdge)
        self._addTopologicalDemiEdge(otherNewDemiEdge)

        # Updating sigma
        self.sigma.addAfter(demiEdge)
        self.sigma.stretch(1)

        # Updating alpha
        self._addEdgeToAlpha()

        # Updating phi
        self.phi.addAfter(self.alpha(demiEdge))
        self.phi.addAfter(newDemiEdge)

        return self.X(otherNewDemiEdge)

    def addEdgeBefore(self, demiEdge: int) -> PrimitiveMutableTopologicalDemiEdge:
        """
        This function will create an new edge attached to the same node as demi edge such that it is before
        demiEdge in the trigonometric order.

        INPUT:

        -``demiEdge`` -- int

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: A = mm.X(10)
            sage: A.n,A.pn
            (X(10), X(10))
            sage: mm.addEdgeAfter(10)
            X(22)
            sage: mm.addEdgeBefore(10)
            X(24)
            sage: A.n,A.pn
            (X(21), X(23))

        NOTE:

            O(1)
        """

        return self.addEdgeAfter(self.sigma.inverseApply(demiEdge))

    def copy(self) -> "PrimitiveMutableLabelledMap":
        """
        OUTPUT:

        Returns A copy of self

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.copy() == mm
            True

        NOTE:

            O(m)
        """
        return PrimitiveMutableLabelledMap(sigma=Permutation(self.sigma.to_list()), alpha=Permutation(self.alpha.to_list()))

    def contractEdge(self, demiEdge: int) -> None:
        """
        Contract in self the edge corresponding to demiEdge.
        WARNING: For this version you demiEdge cannot be on a loop anyway because
        contracting an loop is defined(in this library) as deleting an edge just
        call deleteEdge.

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: mm.faces()
            [(1, 3, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14, 19, 20, 17, 9, 8, 10, 6)]
            sage: mm.contractEdge(3)
            sage: mm.faces()
            [(1, 3, 17, 9, 8, 10, 6, 2, 5, 4, 7, 11, 16, 18, 13, 12, 15, 14)]

        NOTE:

            O(1)
        """

        if self.m == 1:
            raise ValueError(
                "Cannot contract an edge  in a map with only one edge")

        if demiEdge != self.q:
            self.labelToTheEnd([demiEdge, self.alpha(demiEdge)])
            self.contractEdge(self.q)
            return

        otherDemiEdge = self.alpha(demiEdge)

        self.sigma.mergeDelete(demiEdge, otherDemiEdge)
        self.phi.deleteLastKIndex(2)
        self.alpha.deleteLastKIndex(2)

    def areOnTheSameNode(self, demiEdgeA: int, demiEdgeB: int) -> bool:
        """
        Not implemented for PrimitiveMutableLabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.areOnTheSameNode(1,2)
            ....: except:
            ....:    print("OK")
            ....:
            OK
        """
        raise NotImplementedError(self)

    def areOnTheSameFace(self, demiEdgeA: int, demiEdgeB: int) -> bool:
        """
        Not implemented for PrimitiveMutableLabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.areOnTheSameFace(1,2)
            ....: except:
            ....:    print("OK")
            ....:
            OK
        """
        raise NotImplementedError(self)

    def numberInTheSameFace(self, demiEdge: int) -> int:
        """
        Not implemented for PrimitiveMutableLabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.numberInTheSameFace(3)
            ....: except:
            ....:    print("OK")
            ....:
            OK
        """
        raise NotImplementedError(self)

    def numberInTheSameNode(self, demiEdge: int) -> int:
        """
        Not implemented for PrimitiveMutableLabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.numberInTheSameNode(3)
            ....: except:
            ....:    print("OK")
            ....:
            OK


        """
        raise NotImplementedError(self)

    def checkTwoInTheSameFace(self, listDemiEdges: list[int]) -> bool:
        """
        Not implemented for PrimitiveMutableLabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.checkTwoInTheSameFace([3,4])
            ....: except:
            ....:    print("OK")
            ....:
            OK
        """
        raise NotImplementedError(self)

    def checkTwoInTheSameNode(self, listDemiEdges: list[int]) -> bool:
        """
        Not implemented for PrimitiveMutableLabelledMap

        EXAMPLES::

            sage: alpha = Permutation([3, 5, 1, 6, 2, 4, 9, 10, 7, 8, 13, 15, 11, 17, 12, 18, 14, 16, 20, 19])
            sage: sigma = Permutation([2, 4, 3, 1, 5, 7, 8, 6, 11, 10, 12, 14, 16, 9, 15, 13, 19, 18, 17, 20])
            sage: mm = PrimitiveMutableLabelledMap(alpha=alpha,sigma=sigma)
            sage: try:
            ....:    mm.checkTwoInTheSameNode([3,4])
            ....: except:
            ....:    print("OK")
            ....:
            OK
        """
        raise NotImplementedError(self)
