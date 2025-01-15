load("PlanarMap.sage")
import random

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

def getRandomDyckPath(n):
    """
    Returns a random dyck path of size n (uniform random generation)

    Args: 
		int n : size of path
    """
    N = 2 * n + 1
    dyck = [1] * n + [-1] * (n + 1) 
    random.shuffle(dyck)

    level = 0
    minlevel = 0
    posmin = 0

    for i in range(N):
        level += dyck[i]
        if level < minlevel:
            posmin = i + 1
            minlevel = level

    Dyckfinal = [dyck[(posmin + i) % N] for i in range(N - 1)]
    return Dyckfinal

