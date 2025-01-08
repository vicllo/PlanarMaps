load("PlanarMap.sage")

def cube():
	"""Returns the standard cube map."""
	return PlanarMap(adj = [(5,4,2),(1,3,6),(4,7,2),(8,3,1),(8,1,6),(5,2,7),(3,8,6),(7,4,5)])

def complete_map(n):
	"""
	Returns an arbitrary map corresponding to the complete graph with n nodes.
	The genus is guaranteed to be zero if the graph is planar (i.e. n <= 4).
	"""

	adj = list(tuple((j+i)%n + 1 for j in range(1,n)) for i in range(n))
	m = PlanarMap(adj = adj)

	if n <= 4:
		m = m.force_planar()
	
	return m