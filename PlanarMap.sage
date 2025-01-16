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
		Note that it is not possible to build a map with multiedges from an adjacency list.

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
	
	def __eq__(self,other):
		if isinstance(other,self.__class__):
			return self.sigma == other.sigma and self.alpha == other.alpha
		return False

	def _build_from_permutations(self, sigma, alpha):
		r"""
		Initializes the planar map from the underlying permutations.
		"""
		self.sigma = sigma
		self.alpha = alpha
		self.phi = self.alpha.right_action_product(self.sigma)
		self.size = self.sigma.size()
		self.m = self.size // 2 

		if self.sigma.size() != self.alpha.size():
			raise ValueError("The two permutations does not have the same size")

		if self.alpha.right_action_product(self.alpha) != Permutations(self.size).identity():
			raise ValueError("The permutation alpha is not an involution")
		
		if self.alpha.number_of_fixed_points() != 0:
			raise ValueError("The permutation alpha should not have fixed points")

		seen = [False] * (self.size + 1)
		seen[0] = True  # On s'évite les décalages d'indices de la sorte

		def dfs(i):
			seen[i] = True
			if not seen[self.alpha(i)]:
				dfs(self.alpha(i))
			if not seen[self.sigma(i)]:
				dfs(self.sigma(i))
		dfs(1)

		if False in seen:
			raise ValueError("The graph is not connected")


	def _build_from_adj(self, adj):
		n = len(adj)
		
		if sum(map(len, adj)) % 2 != 0:
			raise ValueError("Invalid adjacency list")

		self.m = sum(map(len, adj)) // 2
		pairs = []		# pairs of half-edges corresponding to a single edge (ie. the transpositions of alpha)
		cycles = []		# lists of outgoing half-edges for each vertex (ie. the cycles of sigma)

		edges = {}
		iEdge = 1
		
		for u in range(1, n+1):
			c = []

			for v in adj[u-1]:
				other = None
				if (v, u) in edges:				# the test must be done before setting (u, v) to account for loops
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
		corres = [0] * int(2 * self.m + 1)			# associe à une demi-arête le sommet correspondant
		for i in range(1, len(vertices)+1):
			for k in vertices[i-1]:
				corres[k] = i
		
		edges = []

		for i in range(1, 2*self.m+1):				# pour chaque demi-arête, on ajoute une arête entre corres[i] et corres[alpha(i)]
			if i < self.alpha(i):					# on évite d'ajouter les arêtes en double
				edges.append((corres[i], corres[self.alpha(i)]))

		return Graph(edges, loops = True, multiedges = True)

	def show(self, use_sage_viewer = True, show_vertices = True, show_labels = True):
		"""
		Show the planar map, using the default sage viewer (unless use_sage_viewer is set to False, in which case igraph is used).

		Red nodes are actual graph nodes; black nodes are artificial nodes added to draw graphs with multiedges or loops.

		Note that the fancy igraph viewer might display crossing edges for large graphes, and sometimes will not
			display the correct map (ie. change the order of edges around a node).
		The default sage viewer is guaranteed to be correct, but will often display ugly graphs.
		"""

		vertices = self.sigma.to_cycles()

		real_n_vertices = len(vertices)				# all the new vertices are added to remove multiedges and loops
													# and thus should not be drawn

		alpha = self.alpha
		sigma = self.sigma
		m = self.m

		corres = [0] * int(2 * self.m + 1)			# corres[i] is the vertex corresponding to the half-edge i 
		for i in range(1, len(vertices)+1):
			for k in vertices[i-1]:
				corres[k] = i

		def break_down(i):
			nonlocal alpha, sigma, corres, vertices, m

			# add a new vertex v, and break down the edge whose half-edges are i & alpha(i) by 2 edges (i, 2*m+1) and (2*m+2, alpha(i)) 
			alpha *= Permutation((alpha(i), 2*m+1, i, 2*m+2))
			sigma *= Permutation((2*m+1, 2*m+2))

			corres.append(len(vertices) + 1)		# the two new half-edges 2*m+1 and 2*m+2 are linked to the new vertex
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
				if corres[alpha(i)] in seen_vertices:
					break_down(i)
				else:
					seen_vertices.add(corres[alpha(i)])

		embedding = {}

		for i in range(1, len(vertices)+1):			# build the embedding (list of the neighbors of each edge, in clockwise order)
			embedding[i] = list(corres[alpha(k)] for k in vertices[i-1])
			embedding[i].reverse()					# clockwise order!

		edges = []

		for i in range(1, 2*m+1):				# for each half-edge i, add an edge between corres[i] and corres[alpha(i)]
			if i < alpha(i):					# should not add the same edge twice
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
			vertex_size = 90 / max(1, len(vertices))**.5
		else:
			vertex_size = 0

		if use_sage_viewer:
			if self.genus() == 0:
				layout = "planar"
			else:
				layout = "spring"
			g.show(layout = layout, vertex_size = vertex_size * 8, vertex_labels = show_labels, vertex_colors = {"red": list(range(1, real_n_vertices+1)), "black": list(range(real_n_vertices+1, len(vertices)+1))}, figsize = (8, 8))

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

		return PlanarMap(adj = adj)


	def contractEdge(self, iEdge):
		"""
		Contracts the given half-edge (i. e. merge the two nodes that it links, and removes the edge itself).
		"""

		if iEdge < 1 or iEdge > 2 * self.m:
			raise ValueError("Invalid half-edge number.")

		def buildTransp(l):			# build a permutation from a list of possibly null transpositions
			return Permutation(list(filter(lambda t: t[0] != t[1], l)))		# ie. buildTransp([(1,1),(2,4),(3,3)]) = (1,4,3,2)
		
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
		return PlanarMap(derivedSigma,derivedAlpha)
	


	def quadrangulation(self):
		""" 
		There is bijection between rooted map with m edge and bipartite quadrangulation rooted map with m vertices ,
		this function  return a labelled map(say Q) representant of a rooted quadrangulation associated to self if rooted,
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
		this function  return a labelled map(say T) representant of a face-bicolored rooted tetravalance associated to self if rooted,
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
		return self.numberOfFaces()==1 and self.numberOfEdges() == self.numberOfNodes()-1
