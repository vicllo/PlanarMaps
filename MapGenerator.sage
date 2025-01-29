load("RootedMap.sage")
from collections import deque

import random

class MapGenerator():
    """
    This class representant an abstraction containing methods to generate Map
    """
    def __init__(self):
        pass
        
    def cube():
        """Returns the standard cube map."""
        return RootedMap(adj = [(5,4,2),(1,3,6),(4,7,2),(8,3,1),(8,1,6),(5,2,7),(3,8,6),(7,4,5)])

    def complete_map(n):
        """
        Returns an arbitrary rooted map corresponding to the complete graph with n nodes.
        The genus is guaranteed to be zero if the graph is planar (i.e. n <= 4).
        """
        adj = list(tuple((j+i)%n + 1 for j in range(1,n)) for i in range(n))
        m = RootedMap(adj = adj)
        if n <= 4:
            m = m.force_planar()
        return m

    def getRandomDyckPath(self,n,seed = None):
        """
        Returns a random dyck path of size n (uniform random generation)
        Args: 
            n : size of path
            seed : a random seed if None is used no random seed will be set
        Returns:
            A list of size 2*n +1 for up and  -1 for down in the dyck path
        """
        rng = random.Random()
        if not seed is None:
            rng.seed(int(seed))
        N = 2 * n + 1
        dyck = [1] * n + [-1] * (n + 1)
        rng.shuffle(dyck)
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
    def getRandomPermutation(self,n,seed = None):
        """
        Args:
            -n the size of the permutation
            -seed a random seed if None is used no random seed will be set
        Returns:
            A random permutation of size n
        """
        rng = random.Random()
        if not seed is None:
            rng.seed(int(seed))
        
        lst = [i+1 for i in range(n)]
        rng.shuffle(lst)
        return Permutation(lst)

    def isValidDyckPath(self,dyckPathCandidate):
        """
        Args:
            dyckPathCandidate: A list representing a potential dyckPath
        Returns: 
            A boolean indicating whether or not dyckPathCandidate is a correct dyck path 
        """

        if len(dyckPathCandidate) == 0 or len(dyckPathCandidate)%2 == 1:
            return False

        for e in dyckPathCandidate:
            if e!=-1 and e!=1:
                return False
        
        S = 0

        for e in dyckPathCandidate:
            S+=e 
            if S<0:
                return False
        return S == 0

    def getTreeFromDyckPath(self,dyckPath):
        """
        Given a dyckPath this function will return the associated
        ----
        Args:
            dyckPath: A list representing a dyckPath +1 for up and  -1 for down in the dyck path 
        Returns:
            The corresponding rooted plane tree if dyckPath is a valid dyckPath otherwise it will raise an error
        ----
        O(k) where k = len(dyckPath) 
        """
        if not self.isValidDyckPath(dyckPath):
            raise ValueError("The given list isn't a dyck path")


        phiCycle= []
        alphaCycle = []

        p = []

        for i in range(len(dyckPath)):
            
            phiCycle.append(i+1)

            if dyckPath[i]<0:
                otherDemiEdge = p.pop()
                alphaCycle.append((i+1,otherDemiEdge))
            else:
                p.append(i+1)


        
        phi = Permutation([tuple(phiCycle)])
        alpha = Permutation(alphaCycle)
        
        sigma = phi.left_action_product(alpha)
        
        return RootedMap(alpha = alpha,sigma = sigma)
    
    def getRandomLabellingTree(self,tree,seed = None):
        """
        This generate a  uniformly random labelling( i.e function with value in Z on the demi edge constant taking the same value for each 
        demi edge on the same node and such that if u and v are adjacents nodes abs(f(u)-f(v))<=1 considered up to common translation) of tree.
        It is guaranted that the root demiEdge(i.e 1) has value 0.
        ----
        Args:
            dyckPath: A list representing a dyckPath +1 for up and  -1 for down in the dyck path 
        Returns:
            labelling : A list of size 2*tree.m+1 such that labelling[i] for i>=1 represent the label of the demi edge i , labelling[0] = -1
            but it doesn't have any meaning
        ----
        O(m) where m is the number of edge of tree
        """
        rng = random.Random()
        if not seed is None:
            rng.seed(int(seed))

        sigma = tree.sigma
        alpha = tree.alpha
        
        p = deque()

        nodes = tree.nodes()

        demiEdgeToNodeId = [-1 for i in range(2*tree.m+1)]

        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                demiEdgeToNodeId[nodes[i][j]] = i

        labelling = [-1 for i in range(2*tree.m+1)] 

        startNodeId = demiEdgeToNodeId[1]

        p.append(startNodeId)

        labellingNodes = [-1 for i in range(len(nodes))]


        seen = [False for i in range(len(nodes))]

        seen[startNodeId] = True
        
        transition = [-1,1,0]

        labellingNodes[startNodeId] = 0

        while len(p)>0:
            nodeId = p.popleft()

            for demiEdge in nodes[nodeId]:
                labelling[demiEdge] = labellingNodes[nodeId]

                alphaDemiEdge = alpha(demiEdge)
                alphaNodeId = demiEdgeToNodeId[alphaDemiEdge]
                if not seen[alphaNodeId]:
                    dLabel = rng.sample(transition,1)[0]
                    labellingNodes[alphaNodeId] = dLabel+labellingNodes[nodeId]
                    seen[alphaNodeId] = True
                    p.append(alphaNodeId)
        
        return labelling
    
    def getRandomTree(self,numberOfEdge,seed = None):
        """
        This generate a random random rooted tree uniformly
        ----
        Args:
            numberOfEdges: The number of edge the tree
            seed: a random seed if None no random seed will be set
        Returns:
            A uniformly selected rooted tree with numberOfEdges edges
        ----
        O(numberOfEdge) 
        """

        return self.getTreeFromDyckPath(self.getRandomDyckPath(numberOfEdge,seed = seed))

    def getRandomLabelledTree(self,numberOfEdge,seed =None):
        """
        This generate a random rooted tree and a labelling(a function with value in Z on the demi edge constant taking the same value for each 
        demi edge on the same node and such that if u and v are adjacents nodes abs(f(u)-f(v))<=1 considered up to common translation) uniformly on
        rooted tree of size numberOfEdge
        ----
        Args:
            numberOfEdges: The number of edge the tree
            seed: a random seed if None no random seed will be set
        Returns:
            (tree,labelling) a uniformly selected tuple of tree and labelling on tree demi edges such that tree has numberOfEdge edges
            -tree :  A rooted tree
            -labelling : A list of size 2*numberOfEdge+1 such that labelling[i] for i>=1 represent the label of the demi edge i , labelling[0] = -1
            but it doesn't have any meaning
        ----
        O(numberOfEdge) 
        """
        tree = self.getRandomTree(numberOfEdge,seed = seed)
        return tree,self.getRandomLabellingTree(tree,seed = seed)

    def getRandomPlanarQuadrangulation(self,numberOfFace,seed = None):
        """
        This generate a uniformly random planar rooted quadrangulation  with numberOfFace faces
        ----
        Args:
            numberOfFace: The number of the quadrangulation
            seed: a random seed if None no random seed will be set
        Returns:
            A uniformly selected rooted planar quadrangulation with numberOfFace faces
        ----
        O(numberOfFace) 
        """
        tree,labelling = self.getRandomLabelledTree(numberOfFace,seed = seed)

        quadA,quadB = tree.inverseShaefferTree(returnMarkedDemiEdge = False,labelled = labelling)

        rng = random.Random()

        if not seed is None:
            rng.seed(int(seed))
        
        if rng.random()<0.5:
            return quadB
        return quadA
    
    def getRandomPlanarMap(self,numberOfEdge,seed = None):
        """
        This generate a uniformly random rooted planar map  with numberOfEdge edges
        ----
        Args:
            numberOfEdge: The number of edge of the rooted map
            seed: a random seed if None no random seed will be set
        Returns:
            A uniformly selected rooted planar map with numberOfEdge edge
        ----
        O(numberOfEdge) 
        """
        
        quad = self.getRandomPlanarQuadrangulation(numberOfEdge,seed = seed)

        return quad.inverseQuadrangulation()