from sage.all_cmdline import *   # import sage library
import warnings
from collections import deque
try:
    import igraph
except:
    igraph = None

try:
    import matplotlib.pyplot as plt
except:
    plt = None

class LabelledMap:
    """
    A class to represent a  labelled map.

    Attributes
    ----------
    sigma : Permutation
        Fixed-point free involution whose cycles are given by the edges
    alpha : Permutation
        Permutation that maps a half-edge to the half-edge incident to it in clockwise direction, around the vertex it belongs to.
    
    Methods:
    -------
    """
    
    def __init__(self, sigma:Permutation = None, alpha:Permutation = None, adj = None):
        r"""
        Initializes the labelled map from either the permutations alpha and sigma, or an adjacency list (giving for
                vertex a list of its neighbors in order; vertices must be numbered from 1 to n).

        INPUT:

        - ``sigma`` -- Permutation; Fixed-point free involution whose cycles are given by the edges

        - ``alpha`` -- Permutation; Permutation that maps a half-edge to the half-edge incident to it in clockwise direction, 
          around the vertex it belongs to.

        EXAMPLES:
        sage: sigma = Permutation( [1,3,2,5,4,6])
        sage: alpha = Permutation([(1,2),(3,4),(5,6)])
        sage: LabelledMap(sigma, alpha)
        Labelled map | Sigma : [1, 3, 2, 5, 4, 6], Alpha : [2, 1, 4, 3, 6, 5]

        TESTS::
            sage: sigma = Permutation( [3,4,1,2,6,5])
            sage: alpha = Permutation( [(1,2),(3,4)])
            sage: map = LabelledMap(sigma,alpha)
            Traceback (most recent call last):
            ...
            ValueError: The two permutations do not have the same size
        
            sage: sigma = Permutation([3, 4, 1, 2, 5])
            sage: alpha = Permutation([(1,2),(3,4,5)])
            sage: LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The permutation alpha is not an involution

            sage: sigma = Permutation([3, 4, 1, 2, 5])
            sage: alpha = Permutation([2, 1, 3, 5, 4])
            sage: LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The permutation alpha should not have fixed points

            sage: sigma = Permutation([1,2,3,4,5,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The graph is not connected

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha)
            Labelled map | Sigma : [1, 3, 2, 5, 4, 6], Alpha : [2, 1, 4, 3, 6, 5]

            sage: adj = [(3,),(1,3),(2,)]
            sage: LabelledMap(adj = adj)
            Traceback (most recent call last):
            ...
            ValueError: Invalid adjacency list

        """
        
        if sigma == None and alpha == None and adj == None:
            sigma = Permutation()
            alpha = Permutation()

        if adj != None and (sigma != None or alpha != None):
            raise ValueError("Cannot build the map from both an adjacency list and permutations")

        if adj == None:
            self._build_from_permutations(sigma, alpha)
        else:
            self._build_from_adj(adj)
    
    
    def _build_from_permutations(self, sigma, alpha):
        r"""
        Initializes the labelled  map from the underlying permutations.
        """
        self.sigma = sigma
        self.alpha = alpha
        self.phi = self.alpha.right_action_product(self.sigma)
        self.size = self.sigma.size()
        self.m = self.size // 2 

        if self.sigma.size() != self.alpha.size():
            raise ValueError("The two permutations do not have the same size")

        if self.alpha.right_action_product(self.alpha) != Permutations(self.size).identity():
            raise ValueError("The permutation alpha is not an involution")
        
        if self.alpha.number_of_fixed_points() != 0:
            raise ValueError("The permutation alpha should not have fixed points")

        seen = [False] * (self.size + 1)
        seen[0] = seen[1] = True            # half-edges are numbered from 1 to size, included

        todo = [1]
        while todo:
            i = todo.pop()
            if not seen[self.alpha(i)]:
                todo.append(self.alpha(i))
                seen[self.alpha(i)] = True
            if not seen[self.sigma(i)]:
                todo.append(self.sigma(i))
                seen[self.sigma(i)] = True

        if False in seen:
            raise ValueError("The graph is not connected")

    def _build_from_adj(self, adj):
        n = len(adj)
        
        if sum(map(len, adj)) % 2 != 0:
            raise ValueError("Invalid adjacency list")

        self.m = sum(map(len, adj)) // 2
        pairs = []        # pairs of half-edges corresponding to a single edge (ie. the transpositions of alpha)
        cycles = []        # lists of outgoing half-edges for each vertex (ie. the cycles of sigma)

        edges = {}
        iEdge = 1
        
        for u in range(1, n+1):
            c = []

            for v in adj[u-1]:
                other = None
                if (v, u) in edges:                # the test must be done before setting (u, v) to account for loops
                    other = edges[(v, u)]
                    edges.pop((v, u))

                if other:
                    pairs.append((iEdge, other))
                else:
                    edges[(u, v)] = iEdge

                c.append(iEdge)
                iEdge += 1
            
            cycles.append(tuple(c))
        
        if edges != {}:
            raise ValueError("Invalid adjacency list")

        self._build_from_permutations(Permutation(cycles), Permutation(pairs))

    def buildGraph(self):
        """
        A method that build the multigraph corresponding to the  map.
        Vertices are numbered from 1 to n.
        -------
        Returns:
            A multigraph corresponding to self
        -------
        O(m)
        where m is the number of edges
        """
        vertices = self.sigma.to_cycles()
        corres = [0] * int(2 * self.m + 1)            # associe à une demi-arête le sommet correspondant
        for i in range(1, len(vertices)+1):
            for k in vertices[i-1]:
                corres[k] = i
        
        edges = []

        for i in range(1, 2*self.m+1):                # pour chaque demi-arête, on ajoute une arête entre corres[i] et corres[alpha(i)]
            if i < self.alpha(i):                    # on évite d'ajouter les arêtes en double
                edges.append((corres[i], corres[self.alpha(i)]))

        return Graph(edges, loops = True, multiedges = True)

    def show(self, use_sage_viewer = True, show_vertices = True, show_labels = True):
        """
        Show the planar map, using the default sage viewer (unless use_sage_viewer is set to False, in which case igraph is used).

        Red nodes are actual graph nodes; white nodes are artificial nodes added to draw graphs with multiedges or loops.

        Note that the fancy igraph viewer might display crossing edges for large graphes, and sometimes will not
            display the correct map (ie. change the order of edges around a node).
        The default sage viewer is guaranteed to be correct, but will often display ugly graphs.
        """

        vertices = self.sigma.to_cycles()

        real_n_vertices = len(vertices)                # all the new vertices are added to remove multiedges and loops
                                                    # and thus should not be drawn

        alpha = self.alpha
        sigma = self.sigma
        m = self.m

        corres = [0] * int(2 * self.m + 1)            # corres[i] is the vertex corresponding to the half-edge i 
        for i in range(1, len(vertices)+1):
            for k in vertices[i-1]:
                corres[k] = i

        def break_down(i):
            nonlocal alpha, sigma, corres, vertices, m

            # add a new vertex v, and break down the edge whose half-edges are i & alpha(i) by 2 edges (i, 2*m+1) and (2*m+2, alpha(i)) 
            alpha *= Permutation((int(alpha(i)), int(2*m+1), int(i), int(2*m+2)))
            sigma *= Permutation((int(2*m+1), int(2*m+2)))

            corres.append(len(vertices) + 1)        # the two new half-edges 2*m+1 and 2*m+2 are linked to the new vertex
            corres.append(len(vertices) + 1)

            vertices.append((2*m+1, 2*m+2))

            m += 1

        # for each loop a-a, add a new vertex v and replace the edge a-a with two edges a-v, v-a

        for i in range(1, 2*m+1):
            if corres[i] == corres[alpha(i)]:
                break_down(i)

        # for each vertex v, for each half-edge of v, break down the corresponding edge if there is already an edge between these two vertices
        for v in range(1, len(vertices) + 1):
            seen_vertices = set()
            for i in vertices[v-1]:
                if corres[alpha(int(i))] in seen_vertices:
                    break_down(i)
                else:
                    seen_vertices.add(corres[alpha(int(i))])

        embedding = {}

        for i in range(1, len(vertices)+1):            # build the embedding (list of the neighbors of each edge, in clockwise order)
            embedding[i] = list(corres[alpha(int(k))] for k in vertices[i-1])
            embedding[i].reverse()                    # clockwise order!

        edges = []

        for i in range(1, 2*m+1):                # for each half-edge i, add an edge between corres[i] and corres[alpha(i)]
            if i < alpha(i):                    # should not add the same edge twice
                edges.append((corres[i], corres[alpha(i)]))

        g = Graph(edges, loops = False, multiedges = False)
        g.set_embedding(embedding)

        if not use_sage_viewer and igraph is None:
            warnings.warn("Package igraph not found; falling back to default sage viewer. Consider installing igraph using sage --pip install igraph")
            use_sage_viewer = True
        if not use_sage_viewer and plt is None:
            warnings.warn("Package matplotlib not found; falling back to default sage viewer. Consider installing matplotlib using sage --pip install matplotlib")
            use_sage_viewer = True

        if show_vertices:
            vertex_size = 140 / max(1, len(vertices))**.5
        else:
            vertex_size = 0

        if use_sage_viewer:
            if self.genus() == 0:
                layout = "planar"
            else:
                layout = "spring"
            g.show(layout = layout, vertex_size = vertex_size * 8, vertex_labels = {i: str(i) if i <= real_n_vertices and show_labels else "" for i in range(1, len(vertices)+1)}, vertex_colors = {"red": list(range(1, real_n_vertices+1)), "white": list(range(real_n_vertices+1, len(vertices)+1))}, figsize = (8, 8))

        else:
            if self.genus() == 0:
                layout_dict = g.layout_planar()
            else:
                layout_dict = g.layout()
            
            layout_seed = [layout_dict[i] for i in range(1, len(layout_dict)+1)]

            gg = g.igraph_graph()
            layout = gg.layout_davidson_harel(seed = layout_seed, weight_edge_crossings = float(1e30) * len(vertices)**3)
            layout.fit_into((0,0,1,1))
            
            fig, ax = plt.subplots()
            ax.set_xlim(-0.1,1.1)
            ax.set_ylim(-0.1,1.1)

            igraph.plot(gg, layout = layout, target = ax, vertex_size = vertex_size, vertex_label = list(range(1, len(vertices)+1)), vertex_colors = ["red"] * (real_n_vertices+1) + ["black"] * (len(vertices) - real_n_vertices))
            fig.tight_layout()
            plt.show()


    def __repr__(self):
        return "Labelled map | Sigma : " + str(self.sigma) + ", Alpha : " + str(self.alpha)


    def numberOfFaces(self):
        """
        A method that return the number of faces of the  map
        -------
        Returns:
             The number of faces of self
        -------
        O(m)
        where m is the number of edges
        """
        return len(self.phi.to_cycles())
    

    def numberOfNodes(self):
        """
        A method that returns the number of vertices of the map
        -------
        Returns:
             The number of nodes of self
        -------
        O(m)
        where m is the number of edges
        """
        return len(self.sigma.to_cycles())
    


    def numberOfEdges(self):
        """
        A method that returns the number of edges of the map
        -------
        Returns:
             The number of edge of self
        -------
        O(1)
        """
        return self.m
    
    def genus(self):
        """
        A method that returns the genus of a map
        -------
        Returns:
             The genus of self
        -------
        O(m)
        where m is the number of edges
        """

        return (self.numberOfEdges() + 2 - self.numberOfFaces() - self.numberOfNodes()) // 2

    
    def force_planar(self):
        """
        If the underlying graph is planar, returns a map of genus 0 with the same underlying graph.
        For example, the adjacency list [(2,3,4),(3,4,1),(4,1,2),(1,2,3)] is a valid adjacency list
            for the complete graph with 4 nodes, but its genus is 1; this method could, for example,
            build a map from the list [(2,4,3),(3,4,1),(4,2,1),(1,2,3)] which is the same graph, of genus 0.
        Raises an error if the underlying graph is not planar.
        -------
        Returns:
             Another map corresponding to the above description.
        -------
        """

        g = self.buildGraph()

        if not g.is_planar(set_embedding = True):
            raise ValueError("The force_planar method can be used on maps whose underlying graph is planar.")
        
        e = g.get_embedding()
        adj = [tuple(reversed(e[i])) for i in range(1, len(e)+1)]

        return LabelledMap(adj = adj)

    def getSpanningTree(self):
        """
        A method that returns any spanning tree of the  map
        -------
        Returns:
             A spanning tree of self
        -------
        O(m)
        """

        g = self.buildGraph()
        n = g.order()

        tree = Graph()

        seen = [False] * (n+1)
        seen[0] = True
        
        def dfs(u):
            seen[u] = True
            for v in g.neighbor_iterator(u):
                if not seen[v]:
                    tree.add_edge(u, v)
                    dfs(v)
        
        dfs(1)

        assert False not in seen

        return tree


    def dual(self):
        """
        A method that return the dual of the  map
        -------
        Returns:
             The dual of self
        -------
        O(m)
        where m is the number of edges
        """
        return LabelledMap(self.phi.inverse(),self.alpha)
    

    def diameter(self):
        """
        A method that return the diameter of the  map
        -------
        Returns:
             The diameter of self
        -------
        O(m*n)
        where m is the number of edges and n is the number of nodes
        """
        graph = self.buildGraph()
        return Graph.diameter(graph)
    

    def derivedMap(self):

        """ 
        A method that return the derived Map of the  map
        -------
        Returns:
             The canonical representant of the derived map of self
        -------
        O(m)
        where m is the number of edges
        """
        K = 8*self.m+1
        
        derivedAlphaList = list(range(1,K))
        derivedSigmaList = list(range(1,K))

        invPhi = self.phi.inverse()

        m = int(self.m)


        for i in list(range(1,K)):
            if i<=2*m:
                derivedAlphaList[i-1] = i+2*m
                derivedSigmaList[i-1] = self.sigma(i)
            elif i>2*m and i<=4*m:
                derivedAlphaList[i-1] = i-2*m
                derivedSigmaList[i-1] = i+4*m
            elif i>4*m and i<=6*m:
                derivedAlphaList[i-1] = i+2*m
                derivedSigmaList[i-1] = invPhi(i-4*m)+4*m
            else:
                derivedAlphaList[i-1] = i-2*m
                derivedSigmaList[i-1] = self.alpha(i-6*m)+2*m

        derivedSigma = Permutation(derivedSigmaList)
        derivedAlpha = Permutation(derivedAlphaList)
        return LabelledMap(derivedSigma,derivedAlpha).canonicalRepresentant()
    

    def quadrangulation(self):
        """ 
        There is bijection between rooted map with m edge of genus g and bipartite quadrangulation rooted map
        of genus g with m faces ,
        this function  return the canonical representant of the rooted bipartite quadrangulation associated 
        to self if rooted.
        -------
        Returns:
             The canonical representant of the bipartite rooted quadrangulation associated to rooted(self)
        -------
        O(m)
        where m is the number of edges
        """
        return self.incidenceMap()

    def inverseQuadrangulation(self):
        """
        This function is the inverse of quadrangulation give that self is a bipartite quadrangulation,
        it will return a M a labelled  map such that M.quadrangulation() = self.canonicalRepresentant() 
        and M = M.canonicalRepresentant(), if self isn't a bipartite quadrangulation it will raise an error.
        -------
        Returns:
             The canonical representant of inverse of rooted(self) by quadrangulation if self is a bipartite quadrangulation
            otherwise it will raise an error
        -------
        O(m)
        where m is the number of edges
        """
        bipartition = self.getBipartition()
        
        if bipartition == None or not self.isQuandrangulation():
            raise ValueError("Self isn't a bipartite quadrangulation")
        
        alpha = self.alpha
        sigma = self.sigma
        phi = self.phi

        alphaInvList = [-1 for i in range(self.m)]
        sigmaInvList = [-1 for i in range(self.m)]
        
        colorFace = bipartition[1]
        
        corres = [-1 for i in range(self.m+1)]
        invCorres = [-1 for i in range(2*self.m+1)] 
        
        corres[1] = 1
        invCorres[1] = 1
        
        cnt = 2
        for i in range(2,2*self.m+1):
            if bipartition[i] == colorFace:
                corres[cnt] = i
                invCorres[i] = cnt
                cnt+=1
        
        for invDemiEdge in range(1,self.m+1):
            alphaInvList[invDemiEdge-1] = invCorres[sigma(alpha(sigma(alpha(corres[invDemiEdge]))))]
            sigmaInvList[invDemiEdge-1] = invCorres[alpha(sigma(alpha(corres[invDemiEdge])))]
        alphaInv = Permutation(alphaInvList)
        sigmaInv = Permutation(sigmaInvList)
        return LabelledMap(sigma = sigmaInv,alpha = alphaInv).canonicalRepresentant()

    def incidenceMap(self):
        """ 
        A method that return the incidence map of the  map as its canonical representant
        -------
        Returns:
             Incidence map of self 
        -------
        O(m)
        where m is the number of edges
        """
        
        invPhi = self.phi.inverse()
        invPhiCycles = invPhi.to_cycles()

        quadDemiEdge = 1

        corres = [-1]
        invCorres = list(range(2*self.m+1))

        sigmaQuadList = []

        for k in range(len(invPhiCycles)):
            startQuadDemiEdge = quadDemiEdge
            for demiEdge in invPhiCycles[k]:
                if quadDemiEdge!=startQuadDemiEdge:
                    sigmaQuadList.append(quadDemiEdge)

                corres.append(demiEdge)

                invCorres[demiEdge] = quadDemiEdge
                quadDemiEdge+=1

            sigmaQuadList.append(startQuadDemiEdge)             
        
        numberOfQuadEdge = quadDemiEdge-1

        alphaQuadList = list(range(2*numberOfQuadEdge))

        for quadDemiEdge in range(1,numberOfQuadEdge+1):
            demiEdge = corres[quadDemiEdge]
            turnedDemiEdge = self.sigma(demiEdge)

            quadDemiEdgePrime = invCorres[turnedDemiEdge]

            sigmaQuadList.append(quadDemiEdgePrime+numberOfQuadEdge)

            alphaQuadList[quadDemiEdge-1] = quadDemiEdge+numberOfQuadEdge
            alphaQuadList[quadDemiEdge+numberOfQuadEdge-1] = quadDemiEdge

        alphaQuad = Permutation(alphaQuadList)
        sigmaQuad = Permutation(sigmaQuadList)

        relabelList = [i+1 for i in range(2*numberOfQuadEdge)]

        for quadDemiEdge in range(1,numberOfQuadEdge+1):
            relabelList[quadDemiEdge-1] = corres[quadDemiEdge]
            relabelList[quadDemiEdge+numberOfQuadEdge-1] = relabelList[quadDemiEdge-1]+numberOfQuadEdge

        relabelPerm = Permutation(relabelList)
        return LabelledMap(sigmaQuad,alphaQuad).relabel(relabelPerm).canonicalRepresentant()

    def getRootedMapCorrespondance(self,otherMap,rootDemiEdge):
        """ 
        A method that return a labelling of the demi-edge of self giving otherMap while letting rootDemiEdge 
        invariant if self and otherMap represent the same rooted map at rootDemiEdge otherwise None
        -------
        Args:
              otherMap: The other  map
            rootDemiEdge: The edge on which to root
        Returns:
             t where t is None if they don't represent the same rooted map at rootDemiEdge otherwise 
            t is a permutaion mapping the demi-edge of self to the one of otherMap 
        -------
        O(m)
        where m is the number of edges
        """
        if otherMap.numberOfEdges() != self.numberOfEdges():
            return None
        
        m = self.numberOfEdges()

        tList = [-1 for k in range(2*m)]
        seen = [ False for k in range(2*m)]


        alpha = self.alpha
        sigma = self.sigma 

        sigmaOther = otherMap.sigma
        alphaOther = otherMap.alpha

        tList[rootDemiEdge-1] = rootDemiEdge

        p = []

        p.append(rootDemiEdge)

        seen[rootDemiEdge-1] = True

        while len(p)>0:
            u = p.pop()
            if not seen[alpha(u)-1]:
                seen[alpha(u)-1] = True
                tList[alpha(u)-1] = alphaOther(tList[u-1])
                p.append(alpha(u))

            if not seen[sigma(u)-1]:
                seen[sigma(u)-1] = True
                tList[sigma(u)-1] = sigmaOther(tList[u-1])
                p.append(sigma(u))

        try:
            t = Permutation(tList)
        except:
            return None

        if self.relabel(t) != otherMap: 
            return None

        return t
    
    def relabel(self,tau):
        """ 
        A method that return a relabel LabelledMap , relabelling the demi-edge i by tau(i)
        -------
        Args:
              tau:  A permutation on the demi-edges representing the relabelling
        Returns:
             The relabeled map
        -------
        O(m)
        where m is the number of edges
        """
        
        invTau = tau.inverse()

        relabeledSigma = tau.left_action_product(invTau.right_action_product(self.sigma))

        relabeledAlpha = tau.left_action_product(invTau.right_action_product(self.alpha))

        return LabelledMap(relabeledSigma,relabeledAlpha)

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return self.sigma == other.sigma and self.alpha == other.alpha
        return False

    def tetravalance(self):
        """ 
        There is bijection between rooted map with m edge of genus g and face-bicolorable tetravalant rooted map of genus g 
        with m vertices ,
        this function  return a the canonical representant of the rooted face-bicolorable tetravalance 
        associated to rooted(self).
        -------
        Returns:
             The canonical representant of a tetravalent bicolorable rooted map associated to rooted(self)
        -------
        O(m)
        where m is the number of edges
        """
        return self.edgeMap()


    def edgeMap(self):
        """ 
        A method that return the edge map of the map as its canonical representant 
        -------
        Returns:
             A canonical representant of the edge map of self
        -------
        O(m)
        where m is the number of edges
        """

        invSigma = self.sigma.inverse()
        alpha = self.alpha
        sigma = self.sigma
        m = self.m


        #The number of edge in the edge map
        L = int(2*m)
        
        alphaListEdgeMap = [-1 for k in range(2*L) ]
        sigmaListEdgeMap = [-1 for k in range(2*L) ]

        
        #Construction of alpha and sigma for the edge map
        for k in range(1,L+1):
            alphaListEdgeMap[k-1] = k+L
            alphaListEdgeMap[k+L-1] = k

            t = invSigma(k)
            sigmaListEdgeMap[k-1] = L+t

            j = sigma(k)

            sigmaListEdgeMap[k+L-1] = alpha(j)
        

        alphaEdgeMap = Permutation(alphaListEdgeMap)
        sigmaEdgeMap = Permutation(sigmaListEdgeMap)

        return LabelledMap(sigmaEdgeMap,alphaEdgeMap).canonicalRepresentant()

    def isQuandrangulation(self):
        """
        A function to test wheter or not self is a quandrangulation
        ---
        Returns:
            A boolean indicating if self is a quandrangulation or not
        """
        phi_cycles = self.phi.to_cycles()

        for i in range(len(phi_cycles)):
            if len(phi_cycles[i])!=4:
                return False

        return True
    
    def isBipartite(self):
        """
        Returns : A boolean indicating whether or not self is bipartite
        -------
        O(m)
        where m is the number of edges
        """
        return not (self.getBipartition() is None)

    def getBipartition(self):
        """
        If self isn't bipartite this method will return none otherwise
        it will return a tab clr such that clr[i](=0,1) for a demi edge i give the color of 
        the node on which it is attached.So a node (i.e a  cycle of sigma) will be white(resp black) if 
        all of his element are of color 0(resp 1).clr[0] = -1 cause 0 isn't a valid demi-edge.
        -------
        Returns:
             clr None if self isn't bipartite otherwise give by the above description.
        -------
        O(m)
        where m is the number of edges
        """
        clr = [ -1 for i in range(2*self.m+1)]
        clr[1] = 0
        alpha = self.alpha
        sigma = self.sigma
        phi = self.phi
        p = []
        p.append(1)

        seen = [False for i in range(2*self.m)]

        seen[0] = True
        cnt = 2
        while len(p)>0:
            u = p.pop()
            if not seen[sigma(u)-1]:
                seen[sigma(u)-1] = True
                p.append(sigma(u))
                clr[sigma(u)] = clr[u]

            if not seen[alpha(u)-1]:
                seen[alpha(u)-1] = True
                p.append(alpha(u))
                clr[alpha(u)] = (1+clr[u])%2
        sigma_cycles = sigma.to_cycles()
        for i in range(len(sigma_cycles)):
            r = clr[sigma_cycles[i][0]]
            for j in range(len(sigma_cycles[i])):
                if r!= clr[sigma_cycles[i][j]]:
                    return None
        for i in range(1,2*self.m+1):
            if clr[i] == clr[alpha(i)]:
                return None

        return clr
    
    

    def canonicalRepresentant(self):
        """
        This function return the canonical representant of rooted(self),i.e a
        labelled  map such that M and self are representant of the same rooted map
        and M is the canonical representant.
        -------
        Returns:
             The canonical representant of rooted(self)
        -------
        O(m)
        where m is the number of edges
        """
        relabelList =[-1 for i in range(2*self.m)] 
        alphaCycles = self.alpha.to_cycles()

        sigma = self.sigma
        alpha = self.alpha
        rootDemiEdge = 1
        relabelList[rootDemiEdge-1] = rootDemiEdge

        p = deque()

        p.append(rootDemiEdge)

        seen = [False for i in range(2*self.m)]

        seen[rootDemiEdge-1] = True
        cnt = 2
        while len(p)>0:
            u = p.popleft()
            if not seen[sigma(u)-1]:
                seen[sigma(u)-1] = True
                p.append(sigma(u))        
                relabelList[sigma(u)-1] = cnt
                cnt+=1

            if not seen[alpha(u)-1]:
                seen[alpha(u)-1] = True
                p.append(alpha(u))
                relabelList[alpha(u)-1] = cnt
                cnt+=1
            
        relabel = Permutation(relabelList)
        
        return self.relabel(relabel)

    def isPlaneTree(self):
        """
        A method return a boolean indicating if self is a plane Tree or not
        -------
        Returns:
             A boolean indicating if self is plane tree or not
        -------
        O(m)
        where m is the number of edges
        """
        return self.numberOfFaces()==1 and self.numberOfEdges() == self.numberOfNodes()-1


    def schaefferTree(self,markedDemiEdge):
        """
        The Schaeffer surjection from rooted bipartite quadrangulation of genus g with k face and a marked node to
        rooted one face map (tree in the case g=0) of genus g with k edges and a labelling of its nodes (i.e a function on the nodes of the tree considered up to translation 
        such that if u and v are adjacent f(u) and f(v) differs by atmost one) such that for every rooted one face map T only two rooted marked bipartite quadrangulation give T.
        Given a markDemiEdge which is the corresponding marked node(a node is just a cycle of self.sigma) , this method will return the canonical representant of the 
        rooted one face map associated to rooted(self) and a labelling on its demi edge such that f(node) is the common value of all its demi edge(note that labelling[0] is present but it deosn't have any meaning).
        If self isn't a bipartite quandrangulation this function will raise an error.
        -------
        Args:
            -markedDemiEdge a demi edge on the node which is marked
        Returns:
            - tree: The canonical representant of the rooted one face map corresponding to the above description
            - labelling: A list of labelling on the demi edge of tree cooresponding to the above description
        -------
        O(m)
        where m is the number of edges
        """
        if not self.isBipartite() or not self.isQuandrangulation():
            raise ValueError("Self isn't a bipartite quadrangulation")
        
        sigma = self.sigma
        alpha = self.alpha
        phi = self.phi
        labellingQuad = [-1 for i in range(2*self.m+1)]

        p = deque()

        nodes = sigma.to_cycles()

        demiEdgeToNodeId = [-1 for i in range(2*self.m+1)]

        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                demiEdgeToNodeId[nodes[i][j]] = i

        startNodeId = demiEdgeToNodeId[markedDemiEdge]

        p.append(startNodeId)

        distNodes = [-1 for i in range(len(nodes))]

        distNodes[startNodeId] = 0

        seen = [False for i in range(len(nodes))]

        seen[startNodeId] = True
        
        while len(p)>0:
            nodeId = p.popleft()

            for demiEdge in nodes[nodeId]:
                labellingQuad[demiEdge] = distNodes[nodeId]

                alphaDemiEdge = alpha(demiEdge)
                alphaNodeId = demiEdgeToNodeId[alphaDemiEdge]
                if not seen[alphaNodeId]:
                    distNodes[alphaNodeId] = 1+distNodes[nodeId]
                    seen[alphaNodeId] = True
                    p.append(alphaNodeId)
        
        
        phi_cycles = phi.to_cycles()
        
        corres = [-1 for i in range(2*len(phi_cycles)+1)]
        invCorres = [-1 for i in range(2*self.m+1)]
        
        cnt = 1
        for i in range(len(phi_cycles)):
            A,D,C,B = phi_cycles[i][0],phi_cycles[i][1],phi_cycles[i][2],phi_cycles[i][3]
            link = None
            if labellingQuad[A] == labellingQuad[C]:
                if labellingQuad[B] == labellingQuad[D]:
                    if labellingQuad[A]>labellingQuad[B]:
                        link = (A,C)
                    else:
                        link = (D,B)
                else:
                    if labellingQuad[B]>labellingQuad[D]:
                        link = (C,B)
                    else:
                        link = (D,A)
            else:
                if labellingQuad[A]>labellingQuad[C]:
                    link = (A,B)
                else:
                    link = (C,D)
        
            
            corres[cnt] = link[0]
            invCorres[link[0]] = cnt
            cnt+=1
            corres[cnt] = link[1]
            invCorres[link[1]] = cnt
            cnt+=1

        numberOfTreeDemiEdge = cnt-1

        alphaTreeList = [-1 for i in range(numberOfTreeDemiEdge)]
        sigmaTreeList = [-1 for i in range(numberOfTreeDemiEdge)]

        prelabelling = [-1 for i in range(numberOfTreeDemiEdge+1)]

        for treeDemiEdge in range(1,numberOfTreeDemiEdge+1):
            if treeDemiEdge%2 == 0:
                alphaTreeList[treeDemiEdge-1] = treeDemiEdge-1 
            else:
                alphaTreeList[treeDemiEdge-1] = treeDemiEdge+1
            
            U = corres[treeDemiEdge]
            prelabelling[treeDemiEdge] = labellingQuad[U]
            
            turnU = sigma(U)
            while invCorres[turnU] == -1:
                turnU = sigma(turnU)
            
            sigmaTreeList[treeDemiEdge-1] = invCorres[turnU]
            


        A,D,C,B = 1,phi(1),phi(phi(1)),phi(phi(phi(1)))


        treeRoot = None
        
        if invCorres[A]!=-1 and invCorres[B]!=-1:
            if labellingQuad[A]>labellingQuad[B]:
                treeRoot = invCorres[A]
            else:
                treeRoot = invCorres[B]
        
        elif invCorres[D]!=-1 and invCorres[A]!=-1: 
            if labellingQuad[D]>labellingQuad[A]:
                treeRoot = invCorres[D]
            else:
                treeRoot = invCorres[A]
        
        elif invCorres[B]!=-1 and invCorres[C]!=-1:
            if labellingQuad[B]>labellingQuad[C]:
                treeRoot = invCorres[C]
            else:
                treeRoot = invCorres[B]
        
        elif invCorres[C]!=-1 and invCorres[D]!=-1:
            if labellingQuad[C]>labellingQuad[D]:
                treeRoot = invCorres[D]
            else:
                treeRoot = invCorres[C]

        elif invCorres[A]!=-1 and invCorres[C]!=-1:
            treeRoot = invCorres[A]
            
        else:
            treeRoot = invCorres[B]


        tau = Permutations(len(alphaTreeList)).reflection((1,treeRoot))
        alphaTree = Permutation(alphaTreeList)
        sigmaTree = Permutation(sigmaTreeList)
        
        tree = LabelledMap(alpha=alphaTree,sigma=sigmaTree).relabel(tau)
        
        
        canonicalTree = tree.canonicalRepresentant()
        tauCanonical = tree.getRootedMapCorrespondance(canonicalTree,rootDemiEdge = 1)

        labelling = [-1 for i in range(numberOfTreeDemiEdge+1)]

        for i in range(1,numberOfTreeDemiEdge+1):
            labelling[tauCanonical(tau(i))] = prelabelling[i]

        return canonicalTree,labelling

    def inverseShaefferTree(self,labelled,returnMarkedDemiEdge = True):
        """
        This method is the inverse of the schaefferTree method given that self is a one face map it will return a quadruple
        (quadA,quadB,markedDemiEdgeA,markedDemiEdgeB) where quadA and quaB are in canonical form and we have 
        the following( we will note nodeA and nodeB the node on which markedDemiEdgeA(resp markedDemiEdgeB) is attached in A(resp B)) 
        (rooted(quadA),nodeA) are (rooted(quadB),nodeB) are the only marked rooted quandrangulation such that calling schaefferTree with quadA (resp quadB) with
        any demi edge attached to nodeA (resp nodeB) give rooted(self) in particular quadA.schaefferTree(markedDemiEdgeA) = self.canonicalForm() same for B.
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
        alpha = self.alpha
        sigma = self.sigma

        phi = self.phi
        nextAction = alpha.right_action_product(sigma.inverse())
        maxLabel = labelled[1]
        minLabel = labelled[1]
        for i in range(1,len(labelled)):
            if labelled[i]<minLabel: 
                minLabel = labelled[i]
            if labelled[i]>maxLabel:
                maxLabel = labelled[i]

        for i in range(1,len(labelled)):
            labelled[i] -= minLabel-1

        maxLabel -= minLabel-1
        minLabel = 1

        alphaQuadList = [-1 for i in range(4*self.m)]

        p = [[] for i in range(maxLabel+1)]

        nodes = sigma.to_cycles()
        
        nodesId = [ -1 for i in range(2*self.m+1)]
        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                nodesId[nodes[i][j]] = i+1
                
        sigmaQuadCycleDemiEdge = [deque() for i in range(2*self.m+1)]

        corres = [-1 for i in range(2*self.m+1)]
        cnt = 1
        
        sigmaQuadCycle = [[]]

        root = 1

        partner = [-1 for i in range(2*self.m+1)]

        curDemiEdge = root
        
        oneOrdered = []
        N = 2*self.m
        while N!=0:
            curLabel = labelled[curDemiEdge]
            curNode = -1
            otherNodeId = -1
            if curLabel>1 and len(p[curLabel-1])>0 and partner[curDemiEdge] == -1:
                partner[curDemiEdge] =p[curLabel-1][-1]
                N-=1

            if labelled[curDemiEdge] == 1 and partner[curDemiEdge] == -1:
                partner[curDemiEdge] == -2
                oneOrdered.append(curDemiEdge)
                N-=1
            p[curLabel].append(curDemiEdge)
            curDemiEdge = nextAction(curDemiEdge)
        

        p = [[] for i in range(maxLabel+1)]
        backPartner = [-1 for i in range(2*self.m+1)]
        backAction = nextAction.inverse()
        curDemiEdge = root
        N = 2*self.m
        while N!=0:
            curLabel = labelled[curDemiEdge]
            curNode = -1
            otherNodeId = -1
            if curLabel+1<len(p) and len(p[curLabel+1])>0 and backPartner[curDemiEdge] == -1:
                backPartner[curDemiEdge] = p[curLabel+1][-1]
                N-=1
            if curLabel+1>=len(p) and backPartner[curDemiEdge] == -1:
                backPartner[curDemiEdge] = -2
                N-=1
            p[curLabel].append(curDemiEdge)
            curDemiEdge = backAction(curDemiEdge)
    
        curDemiEdge = root
        
        for i in range(2*self.m):
            curLabel = labelled[curDemiEdge]
            if curLabel>1:
                otherDemiEdge = partner[curDemiEdge]    
                corres[curDemiEdge] = cnt

                sigmaQuadCycleDemiEdge[curDemiEdge].appendleft(cnt)
                sigmaQuadCycleDemiEdge[otherDemiEdge].append(cnt+1)

                alphaQuadList[cnt-1] = cnt+1
                alphaQuadList[cnt] = cnt

                cnt+=2
            
            curDemiEdge = nextAction(curDemiEdge)
        
        for curDemiEdge in oneOrdered:
            otherNodeId = 0
            sigmaQuadCycleDemiEdge[curDemiEdge] = list(sigmaQuadCycleDemiEdge[curDemiEdge])
            corres[curDemiEdge] = cnt

            lst = []
            lst.append(cnt)
            lst += sigmaQuadCycleDemiEdge[curDemiEdge]

            sigmaQuadCycleDemiEdge[curDemiEdge] = lst

            sigmaQuadCycle[otherNodeId].append(cnt+1)

            alphaQuadList[cnt-1] = cnt+1
            alphaQuadList[cnt] = cnt
            cnt+=2

        for curDemiEdge in range(1,2*self.m+1):

            sigmaQuadCycleDemiEdge[curDemiEdge] = list(sigmaQuadCycleDemiEdge[curDemiEdge])
            if len(sigmaQuadCycleDemiEdge[curDemiEdge])==1:
                continue
            curLabel = labelled[curDemiEdge]
            
            firstDemiEdgeQuad = sigmaQuadCycleDemiEdge[curDemiEdge][0]
            restDemiEdgeQuad = sigmaQuadCycleDemiEdge[curDemiEdge][1:]
            

            startId =  -1

            otherDemiEdge = backPartner[curDemiEdge]
            
            P = len(restDemiEdgeQuad)                    

            for i in range(P):
                if alphaQuadList[restDemiEdgeQuad[i]-1] == corres[otherDemiEdge]:
                    startId = i

            newDemiEdgeQuadCycle = []
            newDemiEdgeQuadCycle.append(firstDemiEdgeQuad)

            for i in range(P):
                newDemiEdgeQuadCycle.append(restDemiEdgeQuad[(i+startId)%P])

            sigmaQuadCycleDemiEdge[curDemiEdge] = newDemiEdgeQuadCycle
            

        for node in nodes:
            accum = []
            for demiEdge in node:
                for e in sigmaQuadCycleDemiEdge[demiEdge][1:]:
                    accum.append(e)
                accum.append(sigmaQuadCycleDemiEdge[demiEdge][0])
                
            sigmaQuadCycle.append(accum)

        for i in range(len(sigmaQuadCycle)):
            sigmaQuadCycle[i] = tuple(sigmaQuadCycle[i])
        
        
        sigmaQuad = Permutation(sigmaQuadCycle)
        alphaQuad = Permutation(alphaQuadList)        
        quad = LabelledMap(sigma = sigmaQuad,alpha = alphaQuad)
        numberOfQuadDemiEdge = len(alphaQuadList)

        
        phiQuad = quad.phi
        alphaQuad = quad.alpha
        U = corres[root]
        X,W,V = phiQuad(U),phiQuad(phiQuad(U)),phiQuad(phiQuad(phiQuad(U))) 

        tauA = None
        tauB = None
        if labelled[root] == labelled[alpha(root)]:
            tauA = Permutations(numberOfQuadDemiEdge).reflection((U,root))
            tauB = Permutations(numberOfQuadDemiEdge).reflection((X,root))
            
        elif labelled[root] > labelled[alpha(root)]:
            tauA = Permutations(numberOfQuadDemiEdge).reflection((U,root))
            tauB = Permutations(numberOfQuadDemiEdge).reflection((V,root))
            
        else:
            U = corres[alpha(root)]
            X,W,V = phiQuad(U),phiQuad(phiQuad(U)),phiQuad(phiQuad(phiQuad(U)))
            tauA = Permutations(numberOfQuadDemiEdge).reflection((W,root))
            tauB = Permutations(numberOfQuadDemiEdge).reflection((X,root))        
            


        quadA = quad.relabel(tauA)
        quadB = quad.relabel(tauB)


        if returnMarkedDemiEdge == True:
            quadACanonical = quadA.canonicalRepresentant()
            quadBCanonical = quadB.canonicalRepresentant()
            
            canonicalTauA = quadA.getRootedMapCorrespondance(otherMap = quadACanonical,rootDemiEdge = root)
            canonicalTauB = quadB.getRootedMapCorrespondance(otherMap = quadBCanonical,rootDemiEdge = root)
            
            markedDemiEdge = sigmaQuadCycle[0][0]

            markedDemiEdgeA = tauA(markedDemiEdge)
            markedDemiEdgeB = tauB(markedDemiEdge)

            markedDemiEdgeA = canonicalTauA(markedDemiEdgeA)
            markedDemiEdgeB = canonicalTauB(markedDemiEdgeB)

            return quadACanonical,quadBCanonical,markedDemiEdgeA,markedDemiEdgeB

        return quadA.canonicalRepresentant(),quadB.canonicalRepresentant()
    
    def nodes(self):
        """    
        This function return the nodes of self as cycle of self.sigma
        ----
        Returns :
            A list of cycle of sigma representing nodes of self
        ----
        O(m)
        """
        return self.sigma.to_cycles()
    
    def faces(self):
        """    
        This function return the faces of self as cycle of self.phi
        ----
        Returns : 
            A list of cycle of phi representing faces of self
        ----
        O(m)
        """
        return self.phi.to_cycles()

    def getDyckPath(self,isCanonical=False):
        """
        There is a canonical bijection between rooted planer tree with m edges and dyck path of size m, this method return the associated dyck path
        to rooted(self) if self is a tree if self isn't a tree it will raise an error.
        -----
        Args:
            isCanonical: A boolean indicating if self is already in canonical form
        Returns : 
            A list of size 2*m representing the dyck path associated to rooted(self) +1 for step up and  -1 for step down in the dyck path
            if self is a plane tree otherwise it will raise an error
        -----
        O(m)
        """
        
        if not self.isPlaneTree():
            raise ValueError("Self isn't a plane tree.")
        
        canonicalTree = self.canonicalRepresentant()
        canonicalAlpha = canonicalTree.alpha
        canonicalPhi = canonicalTree.phi

        seen = [False for i in range(2*self.m+1)]

        curDemiEdge = 1

        dyckPath = []

        while not seen[curDemiEdge]:
            print(curDemiEdge,canonicalAlpha(curDemiEdge))
            dyckPath.append(-1 if seen[canonicalAlpha(curDemiEdge)] else 1)
            seen[curDemiEdge] = True
            curDemiEdge = canonicalPhi(curDemiEdge)
        return dyckPath
    