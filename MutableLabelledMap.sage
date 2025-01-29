load("LabelledMap.sage")

class MutableLabelledMap(LabelledMap):
    def contractEdge(self, iEdge):
        """
        Contracts the given half-edge (i. e. merge the two nodes that it links, and removes the edge itself).
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

        self.m -= 1
