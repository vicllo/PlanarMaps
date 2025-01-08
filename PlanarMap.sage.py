

# This file was *autogenerated* from the file PlanarMap.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_90 = Integer(90); _sage_const_p5 = RealNumber('.5'); _sage_const_8 = Integer(8); _sage_const_10 = Integer(10); _sage_const_1e40 = RealNumber('1e40'); _sage_const_3 = Integer(3); _sage_const_0p1 = RealNumber('0.1'); _sage_const_1p1 = RealNumber('1.1'); _sage_const_4 = Integer(4); _sage_const_6 = Integer(6)
from sage.all_cmdline import *   # import sage library
import warnings

try:
	import igraph
except:
	igraph = None

try:
	import matplotlib.pyplot as plt
except:
	plt = None

class PlanarMap:
	"""
	A class to represent a planar map.

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
		Initializes the planar map from either the permutations alpha and sigma, or an adjacency list (giving for
				vertex a list of its neighbors in order; vertices must be numbered from 1 to n).

		INPUT:

		- ``sigma`` -- Permutation; Fixed-point free involution whose cycles are given by the edges

		- ``alpha`` -- Permutation; Permutation that maps a half-edge to the half-edge incident to it in clockwise direction, 
		  around the vertex it belongs to.

		EXAMPLES:
		sage: sigma = Permutation( [1,3,2,5,4,6])
		sage: alpha = Permutation([(1,2),(3,4),(5,6)])
		sage: p = PlanarMap(sigma, alpha)
		sage: p
		Sigma : [1, 3, 2, 5, 4, 6], Alpha : [2, 1, 4, 3, 6, 5]
		
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
		Initializes the planar map from the underlying permutations.
		"""
		self.sigma = sigma
		self.alpha = alpha
		self.phi = self.alpha.right_action_product(self.sigma)
		self.size = self.sigma.size()
		self.m = self.size / _sage_const_2  

		if self.sigma.size() != self.alpha.size():
			raise ValueError("The two permutations does not have the same size")

		if self.alpha.right_action_product(self.alpha) != Permutations(self.size).identity():
			raise ValueError("The permutation alpha is not an involution")
		
		if self.alpha.number_of_fixed_points() != _sage_const_0 :
			raise ValueError("The permutation alpha should not have fixed points")

		seen = [False] * (self.size + _sage_const_1  )
		seen[_sage_const_0 ] = True  # On s'évite les décalages d'indices de la sorte

		def dfs(i):
			seen[i] = True
			if not seen[self.alpha(i)]:
				dfs(self.alpha(i))
			if not seen[self.sigma(i)]:
				dfs(self.sigma(i))
		dfs(_sage_const_1 )

		if False in seen:
			raise ValueError("The graph is not connected")


	def _build_from_adj(self, adj):
		n = len(adj)
		
		if sum(map(len, adj)) % _sage_const_2  != _sage_const_0 :
			raise ValueError("Invalid adjacency list")

		self.m = sum(map(len, adj)) // _sage_const_2 
		pairs = []		# pairs of half-edges corresponding to a single edge (ie. the transpositions of alpha)
		cycles = []		# lists of outgoing half-edges for each vertex (ie. the cycles of sigma)

		edges = {}
		iEdge = _sage_const_1 
		
		for u in range(_sage_const_1 , n+_sage_const_1 ):
			c = []

			for v in adj[u-_sage_const_1 ]:
				other = None
				if (v, u) in edges:				# the test must be done before setting (u, v) to account for loops
					other = edges[(v, u)]
					edges.pop((v, u))

				if other:
					pairs.append((iEdge, other))
				else:
					edges[(u, v)] = iEdge

				c.append(iEdge)
				iEdge += _sage_const_1 
			
			cycles.append(tuple(c))
		
		if edges != {}:
			raise ValueError("Invalid adjacency list")

		self._build_from_permutations(Permutation(cycles), Permutation(pairs))

	def buildGraph(self):
		"""
		A method that build the multigraph corresponding to the planar map.
		Vertices are numbered from 1 to n.
		-------
		Returns:
    		A multigraph corresponding to self
		-------
		O(m)
		where m is the number of edges
		"""
		vertices = self.sigma.to_cycles()
		corres = [_sage_const_0 ] * int(_sage_const_2  * self.m + _sage_const_1 )			# associe à une demi-arête le sommet correspondant
		for i in range(_sage_const_1 , len(vertices)+_sage_const_1 ):
			for k in vertices[i-_sage_const_1 ]:
				corres[k] = i
		
		edges = []

		for i in range(_sage_const_1 , _sage_const_2 *self.m+_sage_const_1 ):				# pour chaque demi-arête, on ajoute une arête entre corres[i] et corres[alpha(i)]
			if i < self.alpha(i):					# on évite d'ajouter les arêtes en double
				edges.append((corres[i], corres[self.alpha(i)]))

		return Graph(edges, loops = True, multiedges = True)

	def show(self, show_vertices = True, use_sage_viewer = False, weight = _sage_const_1 ):
		"""
		Show the planar map, using igraph if possible and the default sage viewer otherwise.
		Does not work (yet) if the map contains loops or multiedges.

		Note that the fancy igraph viewer might display crossing edges for large graphes.
		In this case, increase the weight parameter (a bigger value will decrease the odds of introducing crossing edges,
				but will also reduce the overall display quality); you may also re-run the method until the result is satisfying.
		"""

		# todo: allow multiedges & loops by adding imaginary nodes wherever necessary

		vertices = self.sigma.to_cycles()
		embedding = {}
		corres = [_sage_const_0 ] * int(_sage_const_2  * self.m + _sage_const_1 )			# associe à une demi-arête le sommet correspondant
		for i in range(_sage_const_1 , len(vertices)+_sage_const_1 ):
			for k in vertices[i-_sage_const_1 ]:
				corres[k] = i

		for i in range(_sage_const_1 , len(vertices)+_sage_const_1 ):			# pour chaque sommet, on construit som embedding (liste des voisins dans le sens horaire)
			embedding[i] = list(corres[self.alpha(k)] for k in vertices[i-_sage_const_1 ])
			embedding[i].reverse()					# sens horaire

		edges = []

		for i in range(_sage_const_1 , _sage_const_2 *self.m+_sage_const_1 ):				# pour chaque demi-arête, on ajoute une arête entre corres[i] et corres[alpha(i)]
			if i < self.alpha(i):					# on évite d'ajouter les arêtes en double
				edges.append((corres[i], corres[self.alpha(i)]))

		g = Graph(edges, loops = False, multiedges = False)
		g.set_embedding(embedding)

		if not use_sage_viewer and igraph is None:
			warnings.warn("Package igraph not found; falling back to default sage viewer. Consider installing igraph using sage --pip install igraph")
			use_sage_viewer = True
		if not use_sage_viewer and plt is None:
			warnings.warn("Package matplotlib not found; falling back to default sage viewer. Consider installing matplotlib using sage --pip install matplotlib")
			use_sage_viewer = True

		if show_vertices:
			vertex_size = _sage_const_90  / max(_sage_const_1 , len(vertices))**_sage_const_p5 
		else:
			vertex_size = _sage_const_0 

		if use_sage_viewer:
			if self.genus() == _sage_const_0 :
				layout = "planar"
			else:
				layout = "spring"
			g.show(layout = layout, vertex_size = vertex_size * _sage_const_8 , vertex_labels = False, vertex_color = "red")

		else:
			if weight is None:
				weight = _sage_const_10 **len(vertices)
			if self.genus() == _sage_const_0 :
				layout_dict = g.layout_planar()
			else:
				layout_dict = g.layout()
			
			layout_seed = [layout_dict[i] for i in range(_sage_const_1 , len(layout_dict)+_sage_const_1 )]

			gg = g.igraph_graph()
			layout = gg.layout_davidson_harel(seed = layout_seed, weight_edge_crossings = float(_sage_const_1e40 ) * len(vertices)**_sage_const_3  * weight)
			layout.fit_into((_sage_const_0 ,_sage_const_0 ,_sage_const_1 ,_sage_const_1 ))
			
			fig, ax = plt.subplots()
			ax.set_xlim(-_sage_const_0p1 ,_sage_const_1p1 )
			ax.set_ylim(-_sage_const_0p1 ,_sage_const_1p1 )

			ig.plot(gg, layout = layout, target = ax, vertex_size = vertex_size)
			fig.tight_layout()
			plt.show()


	def __repr__(self):
		return "Sigma : " + str(self.sigma) + ", Alpha : " + str(self.alpha)


	def numberOfFaces(self):
		"""
		A method that return the number of faces of the planar map
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
		A method that returns the number of vertices of the planar map
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
		A method that returns the number of edges of the planar map
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

		return (self.numberOfEdges() + _sage_const_2  - self.numberOfFaces() - self.numberOfNodes()) // _sage_const_2 


	def getSpanningTree(self):
		"""
		A method that returns any spanning tree of the planar map
		-------
		Returns:
     		A spanning tree of self
		-------
		O(m)
		"""

		g = self.buildGraph()
		n = g.order()

		tree = Graph()

		seen = [False] * (n+_sage_const_1 )
		seen[_sage_const_0 ] = True
		
		def dfs(u):
			seen[u] = True
			for v in g.neighbor_iterator(u):
				if not seen[v]:
					tree.add_edge(u, v)
					dfs(v)
		
		dfs(_sage_const_1 )

		assert False not in seen

		return tree


	def dual(self):
		"""
		A method that return the dual of the planar map
		-------
		Returns:
     		The dual of self
		-------
		O(m)
		where m is the number of edges
		"""
		return PlanarMap(self.phi.inverse(),self.alpha)
	

	def diameter(self):
		"""
		A method that return the diameter of the planar map
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
		A method that return the derived Map of the planar map
		-------
		Returns:
     		The derived map of self
		-------
		O(m)
		where m is the number of edges
		"""
		K = _sage_const_8 *self.m+_sage_const_1 
		
		derivedAlphaList = list(range(_sage_const_1 ,K))
		derivedSigmaList = list(range(_sage_const_1 ,K))

		invPhi = self.phi.inverse()

		m = int(self.m)


		for i in list(range(_sage_const_1 ,K)):
			if i<=_sage_const_2 *m:
				derivedAlphaList[i-_sage_const_1 ] = i+_sage_const_2 *m
				derivedSigmaList[i-_sage_const_1 ] = self.sigma(i)
			elif i>_sage_const_2 *m and i<=_sage_const_4 *m:
				derivedAlphaList[i-_sage_const_1 ] = i-_sage_const_2 *m
				derivedSigmaList[i-_sage_const_1 ] = i+_sage_const_4 *m
			elif i>_sage_const_4 *m and i<=_sage_const_6 *m:
				derivedAlphaList[i-_sage_const_1 ] = i+_sage_const_2 *m
				derivedSigmaList[i-_sage_const_1 ] = invPhi(i-_sage_const_4 *m)+_sage_const_4 *m
			else:
				derivedAlphaList[i-_sage_const_1 ] = i-_sage_const_2 *m
				derivedSigmaList[i-_sage_const_1 ] = self.alpha(i-_sage_const_6 *m)+_sage_const_2 *m

		derivedSigma = Permutation(derivedSigmaList)
		derivedAlpha = Permutation(derivedAlphaList)
		return PlanarMap(derivedSigma,derivedAlpha)
    

	def quadrangulation(self):
		""" 
		There is bijection between rooted map with m edge and bipartite quadrangulation rooted map with m vertices ,
		this function  return a labelled map(say Q) representant of a rooted quadrangulation associated to self rooted,
		the coloration is given as follow,a node is black ( i.e a cycle of Q.sigma ) if every demi-edge 
		inside it have label <=2*m otherwise it is white.
		-------
		Returns:
     		A quadrangulation of self
		-------
		O(m)
		where m is the number of edges
		"""
		return self.incidenceMap()
    

	def incidenceMap(self):
		""" 
		A method that return the incidence map of the planar map
		-------
		Returns:
     		Incidence map of self
		-------
		O(m)
		where m is the number of edges
		"""
		
		invPhi = self.phi.inverse()
		invPhiCycles = invPhi.to_cycles()

		quadDemiEdge = _sage_const_1 

		corres = [-_sage_const_1 ]
		invCorres = list(range(_sage_const_2 *self.m+_sage_const_1 ))

		sigmaQuadList = []

		for k in range(len(invPhiCycles)):
			startQuadDemiEdge = quadDemiEdge
			for demiEdge in invPhiCycles[k]:
				if quadDemiEdge!=startQuadDemiEdge:
					sigmaQuadList.append(quadDemiEdge)

				corres.append(demiEdge)

				invCorres[demiEdge] = quadDemiEdge
				quadDemiEdge+=_sage_const_1 

			sigmaQuadList.append(startQuadDemiEdge) 			
		
		numberOfQuadEdge = quadDemiEdge-_sage_const_1 

		alphaQuadList = list(range(_sage_const_2 *numberOfQuadEdge))

		for quadDemiEdge in range(_sage_const_1 ,numberOfQuadEdge+_sage_const_1 ):
			demiEdge = corres[quadDemiEdge]
			turnedDemiEdge = self.sigma(demiEdge)

			quadDemiEdgePrime = invCorres[turnedDemiEdge]

			sigmaQuadList.append(quadDemiEdgePrime+numberOfQuadEdge)

			alphaQuadList[quadDemiEdge-_sage_const_1 ] = quadDemiEdge+numberOfQuadEdge
			alphaQuadList[quadDemiEdge+numberOfQuadEdge-_sage_const_1 ] = quadDemiEdge

		alphaQuad = Permutation(alphaQuadList)
		sigmaQuad = Permutation(sigmaQuadList)

		relabelList = [i+_sage_const_1  for i in range(_sage_const_2 *numberOfQuadEdge)]

		for quadDemiEdge in range(_sage_const_1 ,numberOfQuadEdge+_sage_const_1 ):
			relabelList[quadDemiEdge-_sage_const_1 ] = corres[quadDemiEdge]
		
		relabelPerm = Permutation(relabelList)
		return PlanarMap(sigmaQuad,alphaQuad).relabel(relabelPerm)
    

	def getRootedMapCorrespondance(self,otherMap,rootDemiEdge):
		""" 
		A method that return a labelling of the demi-edge of self giving otherMap while letting rootDemiEdge 
		invariant if self and otherMap represent the same rooted map at rootDemiEdge otherwise None
		-------
		Args:
      		otherMap: The other planar map
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

		tList = [-_sage_const_1  for k in range(_sage_const_2 *m)]
		seen = [ False for k in range(_sage_const_2 *m)]


		alpha = self.alpha
		sigma = self.sigma 

		sigmaOther = otherMap.sigma
		alphaOther = otherMap.alpha

		tList[rootDemiEdge-_sage_const_1 ] = rootDemiEdge

		p = []

		p.append(rootDemiEdge)

		seen[rootDemiEdge-_sage_const_1 ] = True

		while len(p)>_sage_const_0 :
			u = p.pop()
			if not seen[alpha(u)-_sage_const_1 ]:
				seen[alpha(u)-_sage_const_1 ] = True
				tList[alpha(u)-_sage_const_1 ] = alphaOther(tList[u-_sage_const_1 ])
				p.append(alpha(u))

			if not seen[sigma(u)-_sage_const_1 ]:
				seen[sigma(u)-_sage_const_1 ] = True
				tList[sigma(u)-_sage_const_1 ] = sigmaOther(tList[u-_sage_const_1 ])
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
		A method that return a relabel PlanarMap , relabelling the demi-edge i by tau(i)
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

		return PlanarMap(relabeledSigma,relabeledAlpha)
    

	def tetravalance(self):
		""" 
		There is bijection between rooted map with m edge and face-bicolored tetravalant rooted map with m vertices ,
		this function  return a labelled map(say T) representant of a face-bicolored rooted tetravalance associated to self rooted,
		the coloration is given as follow,a face ( i.e a cycle of T.phi ) is black if every demi-edge 
		inside it have label <=2*m otherwise it is white.
		-------
		Returns:
     		A tetravalent map representant of a bi-colored tetravalant map associated to self rooted
		-------
		O(m)
		where m is the number of edges
		"""
		return self.edgeMap()

    
	

	def edgeMap(self):
		""" 
		A method that return the edge Map of the planar map 
		-------
		Returns:
     		The edge map of self
		-------
		O(m)
		where m is the number of edges
		"""

		invSigma = self.sigma.inverse()
		alpha = self.alpha
		sigma = self.sigma
		m = self.m


		#The number of edge in the edge map
		L = int(_sage_const_2 *m)
		
		alphaListEdgeMap = [-_sage_const_1  for k in range(_sage_const_2 *L) ]
		sigmaListEdgeMap = [-_sage_const_1  for k in range(_sage_const_2 *L) ]

		
		#Construction of alpha and sigma for the edge map
		for k in range(_sage_const_1 ,L+_sage_const_1 ):
			alphaListEdgeMap[k-_sage_const_1 ] = k+L
			alphaListEdgeMap[k+L-_sage_const_1 ] = k

			t = invSigma(k)
			sigmaListEdgeMap[k-_sage_const_1 ] = L+t

			j = sigma(k)

			sigmaListEdgeMap[k+L-_sage_const_1 ] = alpha(j)
		

		alphaEdgeMap = Permutation(alphaListEdgeMap)
		sigmaEdgeMap = Permutation(sigmaListEdgeMap)

		return PlanarMap(sigmaEdgeMap,alphaEdgeMap)
    
	
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
		return self.numberOfFaces()==_sage_const_1  and self.numberOfEdges() == self.numberOfNodes()-_sage_const_1 
	

