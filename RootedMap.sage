from collections import deque
load("LabelledMap.sage")
class RootedMap(LabelledMap):
    """
    This class represent rooted map
    """
    def __init__(self,labelledMap = None,sigma = None,alpha = None,adj = None,isAlreadyCanonical = False):
        if labelledMap == None: 
            labelledMap = LabelledMap(sigma = sigma,alpha = alpha,adj = adj)
        canonicalRepresentant = labelledMap
        if not isAlreadyCanonical:
            canonicalRepresentant  = labelledMap.canonicalRepresentant()
        super().__init__(canonicalRepresentant.sigma,canonicalRepresentant.alpha)

    def tetravalance(self):
        """ 
        This method is a bijection between rooted map with m edge of genus g and face-bicolorable tetravalant rooted map of genus g 
        with m vertices,this function  return a the rooted face-bicolorable tetravalance associated to self.
        -------
        Returns:
            A tetravalent face-bicolorable rooted map associated to self
        -------
        O(m)
        where m is the number of edges
        """
        return RootedMap(labelledMap=super().tetravalance(),isAlreadyCanonical = True)
    
    def edgeMap(self):
        """ 
        A method that return the edge map of the map  
        -------
        Returns:
            The edge map of self
        -------
        O(m)
        where m is the number of edges
        """
        return RootedMap(labelledMap = super().edgeMap(),isAlreadyCanonical = True)

    def incidenceMap(self):
        """ 
        A method that return the incidence map of the map  
        -------
        Returns:
            The incidence map of self
        -------
        O(m)
        where m is the number of edges
        """
        return RootedMap(labelledMap = super().incidenceMap(),isAlreadyCanonical = True)

    def quadrangulation(self):
        """ 
        This function is a bijection between rooted map of genus g with m edge and bipartie rooted quandrangulation og genus g with m faces 
        -------
        Returns:
            A tetravalent bi-colorable rooted map associated to self
        -------
        O(m)
        where m is the number of edges
        """
        return RootedMap(labelledMap = super().quadrangulation(),isAlreadyCanonical = True)
    
    def derivedMap(self):
        """
        This function will return the derived map of self
        -------
        Returns:
             The derived map of self
        -------
        O(m)
        """
        return RootedMap(labelledMap = super().derivedMap(),isAlreadyCanonical = True)

    def dual(self):
        """
        A method that return the dual of the map
        -------
        Returns:
             The dual of self
        -------
        O(m)
        where m is the number of edges
        """
        return RootedMap(labelledMap = super().dual())

    def __repr__(self):
        return "Rooted map | Sigma : " + str(self.sigma) +" Alpha : "+str(self.alpha)

    def inverseQuadrangulation(self):
        """
        This function is the inverse of quadrangulation given that self is a bipartite rooted quadrangulation,
        it will return the only rooted map M such that M.quadrangulation() = self.If self isn't a bipartite 
        quadrangulation it will raise an error.
        -------
        Returns:
            The inverse of self from quandrangulation if self is a rooted bipartite quadrangulation otherwise None
        -------
        O(m)
        where m is the number of edges
        """
        return RootedMap(labelledMap=super().inverseQuadrangulation(),isAlreadyCanonical = True)

    def relabel(self,tau):
        """
        This method herited from LabelledMap isn't applicable to RootedMap thus it will just return a copy of self
        """
        return RootedMap( labelledMap = self,isAlreadyCanonical = True)

    def schaefferTree(self,markedDemiEdge):
        """
        The Schaeffer surjection from rooted bipartite quadrangulation of genus g with k face and a marked node to
        rooted one face map (tree in the case g=0) of genus g with k edges and a labelling of its nodes (i.e a function on the nodes of the tree considered up to translation 
        such that if u and v are adjacent f(u) and f(v) differs by atmost one) such that for every rooted one face map T only two rooted marked bipartite quadrangulation give T.
        Given a markDemiEdge which is the corresponding marked node(a node is just a cycle of self.sigma) , this method will return the rooted one face map associated to self 
        and a labelling on its demi edge such that f(node) is the common value of all its demi edge(note that labelling[0] is present but it deosn't have any meaning).
        If self isn't a bipartite quandrangulation this function will raise an error.
        -------
        Args:
            -markedDemiEdge a demi edge on the node which is marked
        Returns:
            - tree: A rooted tree corresponding to the above description
            - labelling: A list of labelling on the demi edge of tree cooresponding to the above description
        -------
        O(m)
        where m is the number of edges
        """
        tree,labelled = super().schaefferTree(markedDemiEdge = markedDemiEdge)
        return RootedMap(labelledMap = tree,isAlreadyCanonical = True),labelled
    
    
    def inverseShaefferTree(self,labelled,returnMarkedDemiEdge = True):
        """
        This method is the inverse of the schaefferTree method given that self is a one face map it will return a quadruple
        (quadA,quadB,markedDemiEdgeA,markedDemiEdgeB) where quadA and quaB are rooted quadrangulations form and we have 
        the following( we will note nodeA and nodeB the node on which markedDemiEdgeA(resp markedDemiEdgeB) is attached in A(resp B)) 
        (quadA,nodeA) are (quadB,nodeB) are the only marked rooted quandrangulation such that calling schaefferTree with quadA (resp quadB) with
        any demi edge attached to nodeA (resp nodeB) give self in particular quadA.schaefferTree(markedDemiEdgeA) = self same for (quadB,markedDemiEdgeB)
        Note that if returnMarkedDemiEdge = False it will only return (quadA,quadB)
        -------
        Args:
            -labelled a list of size 2*m+1 such that for the demiEdge i labelled[i] is the labelled of its attached node,
             0 isn't a valid demiEdge so labelled[0] can take any value it will be ignored. 
            -returnMarkedDemiEdge : a parameter indicating whether or not to return the markedDemiEdge default to true
        Returns:
            -(quadA,quadB,markedDemiEdgeA,markedDemiEdgeB) as in the above description if returnMarkedDemiEdge = True otherwise (quadA,quadB) corresponding to the above description
            ,if self isn't a one face map it will raise an error	
        -------
        O(m)
        where m is the number of edges
        """
        if returnMarkedDemiEdge:
            quadA,quadB,markedDemiEdgeA,markedDemiEdgeB = super().inverseShaefferTree(labelled,returnMarkedDemiEdge = returnMarkedDemiEdge)

            return RootedMap(labelledMap = quadA,isAlreadyCanonical = True),RootedMap(labelledMap = quadB,isAlreadyCanonical = True),markedDemiEdgeA,markedDemiEdgeB
        quadA,quadB = super().inverseShaefferTree(labelled,returnMarkedDemiEdge = returnMarkedDemiEdge)
        return RootedMap(labelledMap = quadA,isAlreadyCanonical = True),RootedMap(labelledMap = quadB,isAlreadyCanonical = True)