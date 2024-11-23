

# This file was *autogenerated* from the file PlanarMap.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_8 = Integer(8); _sage_const_10 = Integer(10); _sage_const_5 = Integer(5); _sage_const_3 = Integer(3); _sage_const_9 = Integer(9); _sage_const_6 = Integer(6); _sage_const_7 = Integer(7); _sage_const_4 = Integer(4)# This file was *autogenerated* from the file PlanarMap.sage
from sage.all_cmdline import *   # import sage library

_sage_const_2 = Integer(_sage_const_2 ); _sage_const_0 = Integer(_sage_const_0 ); _sage_const_1 = Integer(_sage_const_1 ); _sage_const_8 = Integer(_sage_const_8 ); _sage_const_10 = Integer(_sage_const_10 ); _sage_const_5 = Integer(_sage_const_5 ); _sage_const_3 = Integer(_sage_const_3 ); _sage_const_9 = Integer(_sage_const_9 ); _sage_const_6 = Integer(_sage_const_6 ); _sage_const_7 = Integer(_sage_const_7 ); _sage_const_4 = Integer(_sage_const_4 )
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
		
		Methods
		-------

		"""
		self.sigma = sigma
		self.alpha = alpha
		self.phi = self.sigma.right_action_product(self.alpha)
		self.size = self.sigma.size()
		self.n = self.size / _sage_const_2 

		if self.sigma.size() != self.alpha.size():
			raise ValueError("The two permutations does not have the same size")

		if self.alpha.right_action_product(self.alpha) != Permutations(self.size).identity():
			raise ValueError("The permutation alpha is not an involution")
		
		if self.alpha.number_of_fixed_points() != _sage_const_0 :
			raise ValueError("The permutation alpha should not have fixed points")

		seen = [False] * (self.size + _sage_const_1 )
		seen[_sage_const_0 ] = True  # On s'évite les décalages d'indices de la sorte
		def dfs(i):
			seen[i] = True
			if not seen[self.alpha(i)]:
				dfs(self.alpha(i))
			if not seen[self.sigma(i)]:
				dfs(self.sigma(i))
		dfs(_sage_const_1 )

		if false in seen:
			raise ValueError("Le graphe représenté n'est pas connexe")

	def buildGraph(self):
		"""
		A method that build the multigraph corresponding to the planar map
		-------
		O(n)
		"""
		vertices = self.sigma.to_cycles()
		corres = [_sage_const_0 ] * int(_sage_const_2  * self.n + _sage_const_1 )			# associe à une demi-arête le sommet correspondant
		for i in range(_sage_const_1 , len(vertices)+_sage_const_1 ):
			for k in vertices[i-_sage_const_1 ]:
				corres[k] = i
		
		edges = []

		for i in range(_sage_const_1 , _sage_const_2 *self.n+_sage_const_1 ):			# pour chaque demi-arête, on ajoute une arête entre corres[i] et corres[alpha(i)]
			if i < self.alpha(i):					# on évite d'ajouter les arêtes en double
				edges.append((corres[i], corres[self.alpha(i)]))

		return Graph(edges, loops = True, multiedges = True)
			

	def __repr__(self):
		return "Sigma : " + str(self.sigma) + ", Alpha : " + str(self.alpha)

	def numberOfFace(self):
		"""
		A method that return the number of faces of the planar map
		Methods
		-------
		O(n)
		"""
		return len(self.phi.to_cycles())

	def numberOfNodes(self):
		"""
		A method that return the number of vertices of the planar map
		Methods
		-------
		O(n)
		"""
		return len(self.sigma.to_cycles())



