from sage.all_cmdline import *   # import sage library


class PlanarMap:
	def __init__(self, sigma:Permutation, alpha:Permutation):
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
		self.sigma = sigma
		self.alpha = alpha
		self.phi = self.alpha.right_action_product(self.sigma)
		self.size = self.sigma.size()
		self.m = self.size / 2 

		if self.sigma.size() != self.alpha.size():
			raise ValueError("The two permutations does not have the same size")

		if self.alpha.right_action_product(self.alpha) != Permutations(self.size).identity():
			raise ValueError("The permutation alpha is not an involution")
		
		if self.alpha.number_of_fixed_points() != 0:
			raise ValueError("The permutation alpha should not have fixed points")

		seen = [False] * (self.size + 1 )
		seen[0] = True  # On s'évite les décalages d'indices de la sorte

		def dfs(i):
			seen[i] = True
			if not seen[self.alpha(i)]:
				dfs(self.alpha(i))
			if not seen[self.sigma(i)]:
				dfs(self.sigma(i))
		dfs(1)

		if False in seen:
			raise ValueError("The graph isn't connected")

	def buildGraph(self):
		"""
		A method that build the multigraph corresponding to the planar map
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

		for i in range(1, 2*self.m+1):			# pour chaque demi-arête, on ajoute une arête entre corres[i] et corres[alpha(i)]
			if i < self.alpha(i):					# on évite d'ajouter les arêtes en double
				edges.append((corres[i], corres[self.alpha(i)]))

		return Graph(edges, loops = True, multiedges = True)
			
	def __repr__(self):
		return "Sigma : " + str(self.sigma) + ", Alpha : " + str(self.alpha)

	def numberOfFaces(self):
		"""
		A method that return the number of faces of the planar map
		-------
		O(m)
		where m is the number of edges
		"""
		return len(self.phi.to_cycles())

	def numberOfNodes(self):
		"""
		A method that return the number of vertices of the planar map
		-------
		O(m)
		where m is the number of edges
		"""
		return len(self.sigma.to_cycles())
	
	def numberOfEdges(self):
		"""
		A method that return the number of edges of the planar map
		-------
		O(1)
		"""
		return self.m

	def dual(self):
		"""
		A method that return the dual of the planar map
		-------
		O(m)
		where m is the number of edges
		"""
		return  PlanarMap(self.phi.inverse(),self.alpha)
	
	def diameter(self):
		"""
		A method that return the diameter of the planar map
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
		O(n)
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

	def incidenceMap(self):
		""" 
		A method that return the incidence Map of the planar map
		-------
		O(n)
		"""

		face = list(range(2*self.m+1))

		phiCycles = self.phi.to_cycles()
		invPhi = self.phi.inverse()

		for k in range(len(phiCycles)):
			for demiEdge in phiCycles[k]:
				face[demiEdge] = k
		

		t = 1
		sigmaInciList = []
		corres = [-1]
		invCorres  = list(range(2*self.m+1))

		for k in range(len(phiCycles)):
			l = 0
			t0 = t
			while l<len(phiCycles[k]):
				demiEdge = phiCycles[k][l]
				
				if face[demiEdge] == face[self.alpha(demiEdge)]:
					cnt = 1
					l+=1

					while l<len(phiCycles[k]) and face[phiCycles[k][l]] == face[self.alpha(phiCycles[k][l])]:
						cnt+=1
						l+=1

					for j in range(cnt//2):
						t+=1
						corres.append(-1)
						sigmaInciList.append(-1)

				else:
					sigmaInciList.append(-1)
					corres.append(-1)
					corres[t] = demiEdge
					invCorres[demiEdge] = t
					t+=1
					l+=1	

			if len(phiCycles) == 1:
				sigmaInciList.append(-1)
				corres.append(-1)
				t+=1

			for x in range(t0+1,t):
				sigmaInciList[x-1] = x-1
					
			sigmaInciList[t0-1	] = t-1

		N = len(sigmaInciList)


		for j in range(N):
			sigmaInciList.append(-1)

		alphaInciList = list(range(2*N))

		for j in range(1,N+1):
			alphaInciList[j-1] = j+N 
			alphaInciList[j+N-1] = j
			
			if corres[j] == -1:
				sigmaInciList[j+N-1] = j+N
			else:
				demiEdge = corres[j]
				sigmaInciList[j+N-1] = invCorres[self.sigma(demiEdge)]+N

		alphaInci = Permutation(alphaInciList)
		sigmaInci = Permutation(sigmaInciList)

		return PlanarMap(sigmaInci,alphaInci)
		