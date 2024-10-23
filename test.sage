class PlanarMap:
	def __init__(self, sigma:Permutation, alpha:Permutation):
		self.sigma = sigma
		self.alpha = alpha
		self.phi = self.sigma.right_action_product(self.alpha)

	def __repr__(self):
		return "Sigma : " + str(self.sigma) + ", Alpha : " + str(self.alpha)



sigmaP = Permutation([(1, 8, 10, 5, 2), (3, 9, 6, 7, 4)])
alphaP = Permutation([(1, 10), (2, 6), (3, 7), (4, 8), (5, 9)])
P = PlanarMap(sigmaP, alphaP)


print(P.sigma.to_cycles())
print(P.alpha.to_cycles())
print(P.phi.to_cycles())
print(P)
