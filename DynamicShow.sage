load("LabelledMap.sage")

class Vector2D:
    def __init__(self, x = 0, y = 0):
        self.x = x
        self.y = y

    def sqnorm(self):
        return self.x * self.x + self.y * self.y

    def norm(self):
        return sqrt(self.sqnorm())s
    
    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y)

    def __neg__(self):
        return Vector(-self.x, -self.y)

    def __mul__(self, l):
        return Vector(self.x * l, self.y * l)

    def __rmul__(self, l):
        return Vector(self.x * l, self.y * l)

class DynamicShow:
    repulsionVerticesCoef = 1.0			# controls the vertex-vertex repulsion force
    torsionCoef = 1.0					# controls the torsion force between the consecutive edges around a vertex
    springCoef = 1.0					# controls the strength of the spring force for each edge
    repulsionVertexEdgeCoef = 1.0       # controls the strength of the repulsive force between nodes & edges
    weakSpringCoef = 0.01               # controls the strength of the force pulling all nodes towards (0, 0)

    def __init__(self, map: LabelledMap)
        # only works on simple maps so far!
        
        self.map = map

        self.nVertices = map.numberOfNodes()
        self.nEdges = map.numberOfEdges()

        # initialize positions to a valid (but ugly) layout

        self.vertices = map.sigma.to_cycles()

        corres = [0] * (2 * self.m + 1)  # Map half-edge i to its corresponding vertex
        for i in range(1, len(vertices) + 1):
            for k in vertices[i - 1]:
                corres[k] = i

        embedding = {
            i: list(corres[map.alpha(k)] for k in vertices[i - 1])[::-1]
            for i in range(1, self.nVertices + 1)
        }

        edges = [
            (corres[i], corres[alpha(i)]) for i in range(1, 2 * m + 1) if i < alpha(i)
        ]

        g = Graph(edges, loops=False, multiedges=False)
        g.set_embedding(embedding)

        layout = g.layout_planar() if map.genus() == 0 else g.layout()

        self.pos = [layout[i+1] for i in range(self.nVertices)]

    