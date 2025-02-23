load("LabelledMap.sage")


def check_segments_intersecting(p1, q1, p2, q2):
    def onSegment(p, q, r):
        "Check if q lies on [p, r] assuming p,q,r are colinear"
        return q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)
    
    def orient(p, q, r):
        "Returns the orientation of (p, q, r)."
        # 0 if colinear, 1 if clockwise, 2 if counterclockwise

        val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

        if val == 0:
            return 0
        elif val > 0:
            return 1
        return 2
    
    o1, o2, o3, o4 = orient(p1,q1,p2), orient(p1,q1,q2), orient(p2,q2,p1), orient(p2,q2,q1)

    return (
        (o1 != o2 and o3 != o4) or
        (o1 == 0 and onSegment(p1, p2, q1)) or
        (o2 == 0 and onSegment(p1, q2, q1)) or
        (o3 == 0 and onSegment(p2, p1, q2)) or
        (o4 == 0 and onSegment(p2, q1, q2)))

def check_polygon_intersecting(segments):
    "Check whether the given segments intersect each other."
    # TODO: do this on O(n log n) instead of O(N^2)!

    for (p1, q1) in segments:
        for (p2, q2) in segments:
            if p1 == p2 or p1 == q2 or q1 == p2 or q1 == q2:        # this case is already handled by the embedding check
                continue
            
            if check_segments_intersecting(p1, q1, p2, q2):
                return True
    
    return False


class Vector2D:
    def __init__(self, x = 0, y = 0):
        self.x = float(x)
        self.y = float(y)

    def normSq(self):
        return self.x * self.x + self.y * self.y

    def norm(self):
        return sqrt(self.normSq())

    def normalized(self):
        if self.norm() == 0:
            return Vector2D(1, 0)           # hopefully, this should cancel out if trying to normalize a null vector
        return self / self.norm()

    def force_float(self):
        self.x = float(self.x)
        self.y = float(self.y)

    def rotate90(self):
        "Returns a Vector2D corresponding to a rotation of 90 degrees (pi/2 rad) counterclockwise."
        return Vector2D(-self.y, self.x)

    def dot(self, other):
        return self.x * other.x + self.y * other.y

    def angle_towards(self, other):
        "Returns the angle between self and other (counterclockwise) in the range [0, 2pi[."
        a = atan2(self.x * other.y - self.y * other.x, self.x * other.x + self.y * other.y)
        if a < 0:
            return float(a + 2 * pi)
        return float(a)

    def __add__(self, other):
        return Vector2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return Vector2D(self.x - other.x, self.y - other.y)

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

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __lt__(self, other):        # lexicographical order
        return self.x < other.x or (self.x == other.x and self.y < other.y)

    def __hash__(self):
        return hash((self.x, self.y))

class DynamicShow:
    # force expressions are heavily inspired by planarmap.js (https://github.com/tgbudd/planarmap.js)
    
    repulsionVerticesCoef = 3.0			# controls the vertex-vertex repulsion force

    torsionCoef = 2.0					# controls the torsion force between the consecutive edges around a vertex
    torsionPower = 2.0                  # exponent of the angle in the torsion force
    torsionScale = pi / 6               # angles are divided by this quantity in the torsion force
    
    springCoef = 2  					# controls the strength of the spring force for each edge
    springLength = 1.0                  # controls the default length of an edge
    
    repulsionVertexEdgeCoef = 1.0       # controls the strength of the repulsive force between nodes & edges
    repulsionVertexEdgePower = 2.0      # exponent of the node-edge distance   

    weakSpringCoef = 0.01               # controls the strength of the force pulling all nodes towards (0, 0)
    
    maxDispl = 0.5                        # maximum distance a node is allowed to travel during a single iteration
    fast_delta_t = 1                  # delta_t of the first iterations
    slow_delta_t = 0.01                 # delta_t of the last iterations
    max_iter_tick = 5                   # max number of iterations during a single tick

    def __init__(self, map: LabelledMap):
        # only works on simple maps so far!

        # initialize the map
        
        self.map = map

        self.nVertices = map.numberOfNodes()
        self.nEdges = map.numberOfEdges()

        # initialize positions to a valid (but ugly) layout

        self.vertices = map.sigma.to_cycles()

        corres = [0] * (2 * map.m + 1)  # Map half-edge i to its corresponding vertex
        for i in range(1, self.nVertices + 1):
            for k in self.vertices[i - 1]:
                corres[k] = i

        embedding = {
            i: list(corres[map.alpha(k)] for k in self.vertices[i - 1])[::-1]
            for i in range(1, self.nVertices + 1)
        }

        edges_num1 = [          # vertices must be numbered from 1 to nVertices included in networkx
            (corres[i], corres[map.alpha(i)]) for i in range(1, 2 * map.m + 1) if i < map.alpha(i)
        ]

        self.edges = [(i-1, j-1) for (i, j) in edges_num1]      # but internally, it's more convenient to use zero-indexing
        self.faces = [[corres[j]-1 for j in face] for face in map.faces()]
        self.embedding = [[j-1 for j in embedding[i]] for i in range(1, self.nVertices+1)]
        g = Graph(edges_num1, loops=False, multiedges=False)
        g.set_embedding(embedding)

        def minmax(i, j):
            return min(i, j), max(i, j)
        
        self.G = nx.DiGraph()
        self.G.add_edges_from(minmax(i-1, j-1) for (i, j, _) in g.edges())

        layout = g.layout_planar() if map.genus() == 0 else g.layout()

        self.pos = [Vector2D(*layout[i+1]) for i in range(self.nVertices)]
        # self.speed = [Vector2D(0, 0) for i in range(self.nVertices)]

    def centerPos(self):
        "Returns a list of positions (as tuple) centered in [0,1] (with a .05 margin on each side)."

        xs, ys = [v.x for v in self.pos], [v.y for v in self.pos]

        minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)

        for i in range(self.nVertices):
            if minx < maxx:
                xs[i] = (xs[i] - minx) / (maxx - minx) * 0.9 + 0.05
            else:
                xs[i] = 0.5

            if miny < maxy:
                ys[i] = (ys[i] - miny) / (maxy - miny) * 0.9 + 0.05
            else:
                ys[i] = 0.5

        return [(xs[i], ys[i]) for i in range(self.nVertices)]

    def start(self):
        # initialize the matplotlib figure
        size = 7

        self.fig, (self.top_ax, self.ax) = plt.subplots(2, 1, figsize = (size, size+1), height_ratios = [1, size])
        
        pos_centered = self.centerPos()

        self.nodes_plt = nx.draw_networkx_nodes(self.G, pos_centered, ax = self.ax, nodelist = list(range(self.nVertices)), node_color="red", node_size = min(300, 1000 / self.nVertices**.5))
        self.edges_plt = nx.draw_networkx_edges(self.G, pos_centered, ax = self.ax, arrows = False)

        self.ax.axis("off")
        self.ax.set_xlim(left=0, right=1)
        self.ax.set_ylim(bottom=0, top=1)

        self.top_ax.axis("off")
        self.top_ax.set_xlim(left=0, right=1)
        self.top_ax.set_ylim(bottom=0, top=1)
        
        self.text = self.top_ax.text(0, 0.5, "Frame: 0")

        self.fig.tight_layout()

        self.anim = FuncAnimation(self.fig, self.update_fig, cache_frame_data = False, blit = True, interval = 1)

    def computeRepulsionVVForces(self):
        for i in range(self.nVertices):
            for j in range(i+1, self.nVertices):
                r = self.pos[j] - self.pos[i]
                force = -self.repulsionVerticesCoef * r.normalized() / r.normSq()     # force applied to i (going away from j)
                self.forces[i] += force
                self.forces[j] -= force

    def computeSpringForces(self):
        for (i, j) in self.edges:
            r = self.pos[j] - self.pos[i]
            force = self.springCoef * (r.norm() - self.springLength) * r.normalized()

            self.forces[i] += force
            self.forces[j] -= force

    def computeWeakSpringForces(self):
        for i in range(self.nVertices):
            r = self.pos[i]
            force = -self.weakSpringCoef * r.norm() * r.normalized()

            self.forces[i] += force

    def computeTorsionForces(self):
        for i in range(self.nVertices):
            if len(self.embedding[i]) > 1:
                for k in range(len(self.embedding[i])):
                    u1 = self.embedding[i][k]                                   # simulates an angular spring force between u1 and u2
                    u2 = self.embedding[i][(k+1) % len(self.embedding[i])]
                    
                    v1, v2 = self.pos[u1] - self.pos[i], self.pos[u2] - self.pos[i]

                    angle = abs(v1.angle_towards(v2))
                    
                    energy = self.torsionCoef * float(tanh(angle / self.torsionScale)) ** self.torsionPower

                    scale = -2 * self.torsionPower * energy / self.torsionScale / sinh(2 * angle / self.torsionScale)

                    force1 = v1.rotate90() * (scale / v1.normSq())
                    force2 = v2.rotate90() * (scale / v2.normSq())

                    self.forces[u1] += force1
                    self.forces[i] += force2 - force1
                    self.forces[u2] -= force2

    def computeNodeEdgeForce(self, node, e1, e2):   # compute the repulsion force between node and edge (e1, e2)
        if node == e1 or node == e2:
            return

        distdiff = (self.pos[node] - self.pos[e1]).norm() + (self.pos[node] - self.pos[e2]).norm() - (self.pos[e2] - self.pos[e1]).norm()

        energy = (distdiff / self.springLength) ** (-self.repulsionVertexEdgePower)

        scale = self.repulsionVertexEdgeCoef * energy / distdiff

        e1ton = (self.pos[node] - self.pos[e1]).normalized()
        e2ton = (self.pos[node] - self.pos[e2]).normalized()
        e1toe2 = (self.pos[e2] - self.pos[e1]).normalized()

        self.forces[node] += (e1ton + e2ton) * scale
        self.forces[e1] -= (e1ton - e1toe2) * scale
        self.forces[e2] -= (e2ton + e1toe2) * scale

    def computeRepulsionEEForces(self):
        for face in self.faces:
            for i in range(len(face)):      # iterate over every pair of edges
                for j in range(len(face)):
                    self.computeNodeEdgeForce(face[i], face[j], face[(j+1) % len(face)])

    def check_pos_correct(self, pos):
        # check if the embedding is still respected
        for i in range(self.nVertices):
            if len(self.embedding) > 2:
                angle_sum = 0
                for k in range(len(self.embedding[i])):
                    u1 = self.embedding[i][k]                                   # simulates an angular spring force between u1 and u2
                    u2 = self.embedding[i][(k+1) % len(self.embedding[i])]
                    
                    angle_sum += float(2*pi - (pos[u1] - pos[i]).angle_towards(pos[u2] - pos[i]))   # embedding is given clockwise

                # if the embedding is correct, angle_sum = 2pi
                # otherwise, it's way bigger (>= 4pi?)

                if angle_sum >= 2.1 * pi:  # account for float inaccuracies
                    print ("bad!! bad angle (sum", angle_sum, "!)")
                    print (pos[i])
                    for k in self.embedding[i]:
                        print (pos[k])
                    return False
        
        # check that there is no edge crossing
        # it suffices to check that for each face, none of the edges cross each other
        
        for face in self.faces:
            segments = []
            for i in range(len(face)):
                segments.append(tuple(sorted((pos[face[i]], pos[face[(i+1) % len(face)]]))))

            # remove all identical segments
            segments = list(set(segments))

            if check_polygon_intersecting(segments):
                print ("bad!! intersecting")
                print (segments)
                return False

        return True


    def update_forces(self):
        self.forces = [Vector2D() for i in range(self.nVertices)]

        self.computeRepulsionVVForces()
        self.computeSpringForces()
        self.computeWeakSpringForces()
        self.computeTorsionForces()
        self.computeRepulsionEEForces()

    def tick(self, frame):
        full_force_frames = self.nVertices * 2     # after this number of iterations, we should be closer to the optimum; start going slower to refine the positions
        power = 1.5                                   # delta_t will be proportional to 1 / frame ^ power when this number is reached

        if frame < full_force_frames:
            delta_t = self.fast_delta_t
            iters = 1
        else:
            delta_t = self.slow_delta_t * full_force_frames**power / frame**power
            iters = self.max_iter_tick      # we do several iterations with geometrically decreasing delta_t to reach the optimum without flickering

        while iters > 0:
            self.update_forces()
            
            maxForce = max(map(Vector2D.normSq, self.forces))
            #if maxForce > 0:
            #    delta_t = min(delta_t, 2.0 / maxForce)

            newpos = [Vector2D() for i in range(self.nVertices)]

            for i in range(self.nVertices):
                displ = self.forces[i] * delta_t
                if displ.normSq() > self.maxDispl ** 2:
                    displ = displ.normalized() * self.maxDispl
                
                newpos[i] = self.pos[i] + displ

                # considering speed (ie. assuming nodes have inertia) seems to make the convergence slower
                # self.pos[i] += self.speed[i] * self.maxDelta + 0.5 * self.forces[i] * self.maxDelta * self.maxDelta / self.vertexInertia
                # self.speed[i] /= 2              # friction forces (we need energy loss somewhere)
                # self.speed[i] += self.forces[i] * self.maxDelta / self.vertexInertia
                
            if self.check_pos_correct(newpos):
                self.pos = newpos
                iters -= 1
                delta_t /= 1.5
            else:
                delta_t /= 2
        
    def update_fig(self, frame):
        if frame > 5:
            self.tick(frame)
            #self.pos[0].x += (random.random() - 0.5) / 5

        self.text.set_text("Frame: " + str(frame))

        #plt.axis("on")
        #plt.cla()

        # self.pos[0] += Vector2D((random.random()-.5), (random.random()-.5))
        
        # nx.draw(
        #     self.G,
        #     {i: (float(self.pos[i-1].x), float(self.pos[i-1].y)) for i in range(1, self.nVertices + 1)},
        #     ax=self.ax,
        #     labels={
        #         i: str(i)
        #         for i in range(1, self.nVertices + 1)
        #     },
        #     node_size=[300] * self.nVertices,
        #     nodelist=list(range(1, self.nVertices + 1)),
        #     arrows=False,
        #     with_labels=True,
        #     node_color="red",
        # )

        pos_centered = self.centerPos()

        edges_pos = [(pos_centered[i], pos_centered[j]) for (i, j) in self.edges]

        self.nodes_plt.set_offsets(pos_centered)
        self.edges_plt.set_segments(edges_pos)

        return self.nodes_plt, self.edges_plt, self.text