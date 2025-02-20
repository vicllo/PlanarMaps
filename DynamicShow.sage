load("LabelledMap.sage")

class Vector2D:
    def __init__(self, x = 0, y = 0):
        self.x = x
        self.y = y

    def sqnorm(self):
        return self.x * self.x + self.y * self.y

    def norm(self):
        return sqrt(self.sqnorm())

    def normalized(self):
        return self / self.norm()

    def __add__(self, other):
        return Vector2D(self.x + other.x, self.y + other.y)

    def __neg__(self):
        return Vector2D(-self.x, -self.y)

    def __mul__(self, l):
        return Vector2D(self.x * l, self.y * l)

    def __rmul__(self, l):
        return Vector2D(self.x * l, self.y * l)

    def __truediv__(self, l):
        return Vector2D(self.x / l, self.y / l)

    def __str__(self):
        return "(" + str(self.x) + "; " + str(self.y) + ")"

    def __repr__(self):
        return "Vector2D" + str(self)

class DynamicShow:
    repulsionVerticesCoef = 1.0			# controls the vertex-vertex repulsion force
    torsionCoef = 1.0					# controls the torsion force between the consecutive edges around a vertex
    springCoef = 1.0					# controls the strength of the spring force for each edge
    springLength = 1.0                  # controls the default length of an edge
    repulsionVertexEdgeCoef = 1.0       # controls the strength of the repulsive force between nodes & edges
    weakSpringCoef = 0.01               # controls the strength of the force pulling all nodes towards (0, 0)

    def __init__(self, map: LabelledMap):
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

        self.edges = [
            (corres[i], corres[alpha(i)]) for i in range(1, 2 * map.m + 1) if i < alpha(i)
        ]

        g = Graph(self.edges, loops=False, multiedges=False)
        g.set_embedding(embedding)

        layout = g.layout_planar() if map.genus() == 0 else g.layout()

        self.pos = [layout[i+1] for i in range(self.nVertices)]
        
    def computeRepulsionVVForces(self, forces):
        for i in range(self.nVertices):
            for j in range(i+1, self.nVertices):
                r = self.pos[j] - self.pos[i]
                force = -self.repulsionVerticesCoef * r.normalized() / r.sqnorm()     # force applied to i (going away from j)
                forces[i] += force
                forces[j] -= force

    def computeSpringForces(self, forces):
        for (i, j) in self.edges:
            r = self.pos[j] - self.pos[i]
            force = self.springCoef * (r.norm() - self.springLength) * r.normalized()

            forces[i] += force
            forces[j] -= force

    def computeWeakSpringForces(self, forces):
        for i in range(self.nVertices):
            r = self.pos[i]
            force = -self.weakSpringCoef * r.norm() * r.normalized()

            forces[i] += force

    def update(self, delta_t):
        forces = [Vector2D() for i in range(self.nVertices)]

        self.computeRepulsionVVForces(forces)
        self.computeSpringForces(forces)
        self.computeWeakSpringForces(forces)

        for i in range(nVertices):
            pos[i] += 0.5 * forces[i] * delta_t ** 2