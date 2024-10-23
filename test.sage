class PlanarMap:
	def __init__(self, sigma:Permutation, alpha:Permutation):
		self.sigma = sigma
		self.alpha = alpha
		self.phi = self.sigma.right_action_product(self.alpha)
		self.size = self.sigma.size()
		self.n = self.size / 2

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

		if false in seen:
			raise ValueError("Le graphe représenté n'est pas connexe")

	def __repr__(self):
		return "Sigma : " + str(self.sigma) + ", Alpha : " + str(self.alpha)



sigmaP = Permutation([(1, 8, 10, 5, 2), (3, 9, 6, 7, 4)])
alphaP = Permutation([(1, 10), (2, 6), (3, 7), (4, 8), (5, 9)])

P = PlanarMap(sigmaP, alphaP)

print(P.sigma.to_cycles())
print(P.alpha.to_cycles())
print(P.phi.to_cycles())
print(P)

sigmaP = Permutation([1, 2, 3, 4])
alphaP = Permutation([2, 1, 4, 3])

P = PlanarMap(sigmaP, alphaP)

print(P.sigma.to_cycles())
print(P.alpha.to_cycles())
print(P.phi.to_cycles())
print(P)
