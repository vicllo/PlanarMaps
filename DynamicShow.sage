load("LabelledMap.sage")

import matplotlib.animation as animation
import matplotlib.widgets
import matplotlib.gridspec

import math         # we use python's float type and math functions to avoid using the more accurate but much slower sage types
                    # eg. we want pi - 1 to be stored as its floating-point value instead of this formal expression
                    # this is why all numbers in this file are written with "r"s: type(1.0r+math.pi) = float, but type(1.0 + pi) = sage.symbolic.expression.Expression
import time

def check_segments_intersecting(p1, q1, p2, q2):
    def onSegment(p, q, r):
        "Check if q lies on [p, r] assuming p,q,r are colinear"
        return q.x <= max(p.x, r.x) and q.x >= min(p.x, r.x) and q.y <= max(p.y, r.y) and q.y >= min(p.y, r.y)
    
    def orient(p, q, r):
        "Returns the orientation of (p, q, r)."
        # 0 if colinear, 1 if clockwise, 2 if counterclockwise

        val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y)

        if val == 0r:
            return 0r
        elif val > 0r:
            return 1r
        return 2r
    
    o1, o2, o3, o4 = orient(p1,q1,p2), orient(p1,q1,q2), orient(p2,q2,p1), orient(p2,q2,q1)

    return (
        (o1 != o2 and o3 != o4) or
        (o1 == 0r and onSegment(p1, p2, q1)) or
        (o2 == 0r and onSegment(p1, q2, q1)) or
        (o3 == 0r and onSegment(p2, p1, q2)) or
        (o4 == 0r and onSegment(p2, q1, q2)))

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
    def __init__(self, x = 0.0r, y = 0.0r):
        self.x = float(x)
        self.y = float(y)

    def normSq(self):
        return self.x * self.x + self.y * self.y

    def norm(self):
        return sqrt(self.normSq())

    def normalized(self):
        if self.norm() == 0.0r:
            return Vector2D(1.0r, 0.0r)           # hopefully, this should cancel out if trying to use a normalized null vector
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
        a = math.atan2(self.x * other.y - self.y * other.x, self.x * other.x + self.y * other.y)
        if a < 0.0r:
            return a + 2.0r * math.pi
        return a

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
    
    repulsionVVRRCoef = 3.0r			# controls the vertex-vertex repulsion force between red nodes (real nodes)
    repulsionVVRWCoef = 2.0r           # controls the vertex-vertex repulsion force between a red node and a white node (artificial node to split multiedges & loops)
    repulsionVVWWCoef = 2.0r            # controls the vertex-vertex repulsion force between white nodes

    torsionCoef = 3.0r					# controls the torsion force between the consecutive edges around a vertex
    torsionPower = 3.0r                  # exponent of the angle in the torsion force
    torsionScale = math.pi / 6r               # angles are divided by this quantity in the torsion force
    
    springCoef = 3.0r  					# controls the strength of the spring force for each edge
    springLength = 1.0r                 # controls the default length of an edge
    useLogSpring = True                 # if true, spring force is proportional to log(r / springLength); if false, proportional to r - springLength
    
    repulsionVertexEdgeCoef = 3.0r       # controls the strength of the repulsive force between nodes & edges
    repulsionVertexEdgePower = 2.0r      # exponent of the node-edge distance   

    weakSpringCoef = 0.01r               # controls the strength of the force pulling all nodes towards (0, 0)
    
    maxDispl = springLength * 2.0r       # maximum distance a node is allowed to travel during a single iteration
    fast_delta_t = 1.0r                  # delta_t of the first iterations
    slow_delta_t = 0.01r                 # delta_t of the last iterations

    convergence_limit = 0.01r              # consider that the optimum is found when all nodes move less than convergenceLimit * max_extent in a single slow_delta_t frame
    
    max_iter_tick = 4r                   # max number of iterations during a single tick

    time_profile = True                  # print informations about time use of different functions
    break_down_num = 5r                  # break each duplicate edge into this number of small edges

    useVVForce = True
    useTorsionForce = True
    useSpringForce = True
    useWeakSpringForce = True
    useEEForce = True

    validFrames = 0
    validFramesThreshold = 20           # start increasing current_delta_t -> delta_t after this number of valid frames

    def __init__(self, map: LabelledMap):
        # initialize the map
        
        self.map = map
        
        vertices = map.sigma.to_cycles()

        real_n_vertices = len(vertices)  # Remove multiedges and loops
        real_n_halfedges = map.m         # These half-edges should not be drawn

        self.real_n_vertices = real_n_vertices

        alpha = map.alpha
        sigma = map.sigma
        m = map.m

        def minmax(i, j):
            """Ensure edges always go from lowest to highest vertex id."""
            return min(i, j), max(i, j)

        corres = [0] * (2 * m + 1)  # Map half-edge i to its corresponding vertex
        for i in range(1, len(vertices) + 1):
            for k in vertices[i - 1]:
                corres[k] = i

        def break_down(i, break_down_num):
            nonlocal alpha, sigma, corres, vertices, m
            # Add a new vertex v, and break down the edge whose half-edges 
            # are i & alpha(i) into break_down_num edges (i, 2*m+1), (2*m+2, 2*m+3), .., (2*m+2*(break_down_num-1), alpha(i)).

            # print ("breaking down", i)

            alpha_cycles = [(alpha(i), 2*m + 1, i, 2*m + 2*(break_down_num-1))] + [(2*k, 2*k + 1) for k in range(m + 1, m + break_down_num - 1)]
            alpha *= Permutation(alpha_cycles)

            sigma_cycles = [(2*k - 1, 2 * k) for k in range(m+1, m + break_down_num)]
            sigma *= Permutation(sigma_cycles)

            for k in range(break_down_num - 1):
                corres.append(len(vertices) + 1)
                corres.append(len(vertices) + 1)

                vertices.append((2 * m + 1, 2 * m + 2))
                m += 1

        def break_loop(i):
            j = alpha(i)
            vertex = len(vertices)

            break_down(i, max(self.break_down_num, 3))

        # For each loop a-a, add a new vertex v and replace the edge a-a
        # with two edges a-v, v-a.
        for i in range(1, 2 * m + 1):
            if corres[i] == corres[alpha(i)]:
                break_loop(i)

        break_all_down = True           # if True, break down each multiedge; if False, break down all but one

        # Handle each vertex and break down edges if needed.
        for v in range(1, len(vertices) + 1):
            seen_vertices = {}
            for i in vertices[v - 1]:
                if corres[alpha(i)] in seen_vertices:
                    if seen_vertices[corres[alpha(i)]] != -1:
                        duplicate_he = seen_vertices[corres[alpha(i)]]
                        seen_vertices[corres[alpha(i)]] = -1
                        break_down(duplicate_he, self.break_down_num)
                    break_down(i, self.break_down_num)
                else:
                    seen_vertices[corres[alpha(i)]] = i

        # Build the graph embedding
        embedding = {
            i: list(corres[alpha(k)] for k in vertices[i - 1])[::-1]
            for i in range(1, len(vertices) + 1)
        }

        edges_num1 = [       # vertices must be numbered from 1 to nVertices included in networkx
            (corres[i], corres[alpha(i)]) for i in range(1, 2 * m + 1) if i < alpha(i)
        ]

        self.alpha, self.sigma = alpha, sigma

        self.nVertices = len(vertices)
        self.nEdges = m

        # initialize positions to a valid (but ugly) layout

        self.edges = [(i-1, j-1) for (i, j) in edges_num1]      # but internally, it's more convenient to use zero-indexing
        self.faces = [[corres[j]-1 for j in face] for face in alpha.right_action_product(sigma).to_cycles()]
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

        self.max_extent = 1

        self.anim_running = True
        self.done = False

    def centerPos(self):
        "Returns a list of positions (as tuple) centered in [0,1] (with a .05 margin on each side)."

        xs, ys = [v.x for v in self.pos], [v.y for v in self.pos]

        minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)

        for i in range(self.nVertices):
            if minx < maxx:
                xs[i] = (xs[i] - minx) / (maxx - minx) * 0.9r + 0.05r
            else:
                xs[i] = 0.5r

            if miny < maxy:
                ys[i] = (ys[i] - miny) / (maxy - miny) * 0.9r + 0.05r
            else:
                ys[i] = 0.5r

        self.max_extent = max(maxx - minx, maxy - miny)

        return [(xs[i], ys[i]) for i in range(self.nVertices)]

    def start(self, frameByFrame = False):
        # initialize the matplotlib figure
        size = 7

        self.fig = plt.figure(figsize = (size, size+1), layout = "constrained")
        gs = matplotlib.gridspec.GridSpec(2, 2, figure = self.fig, width_ratios = [1, 4], height_ratios = [1, size])

        self.ax = self.fig.add_subplot(gs[1, :])
        txt_ax = self.fig.add_subplot(gs[0, 0])
        slider_ax = self.fig.add_subplot(gs[0, 1])

        pos_centered = self.centerPos()

        self.nodes_plt = nx.draw_networkx_nodes(self.G, pos_centered, ax = self.ax, nodelist = list(range(self.nVertices)),
                    node_color="red", node_size=[min(300, 1000 / self.nVertices**.5)] * self.real_n_vertices + [0] * (self.nVertices - self.real_n_vertices))
        self.edges_plt = nx.draw_networkx_edges(self.G, pos_centered, ax = self.ax, arrows = False)

        self.ax.axis("off")
        self.ax.set_xlim(left=0, right=1)
        self.ax.set_ylim(bottom=0, top=1)

        txt_ax.axis("off")
        txt_ax.set_xlim(left=0, right=1)
        txt_ax.set_ylim(bottom=0, top=1)
        
        self.text = txt_ax.text(0, 0.5, "Frame: 0")

        self.slider = matplotlib.widgets.Slider(ax = slider_ax, valmin = -3, valmax = 1, valinit = 0, label = "delta_t")
        self.slider.on_changed(self.update_slider_val)
        self.slider.valtext.set_text("1.000")
        self.force_delta_t = False

        #self.fig.tight_layout()

        self.fig.canvas.mpl_connect('button_press_event', self.onClick)
        self.fig.canvas.mpl_connect('key_press_event', self.onKey)

        self.currentMaxDispl = self.maxDispl
        
        self.doFrame = False
        self.frame = 0

        self.anim = FuncAnimation(self.fig, self.update_fig, cache_frame_data = False, blit = True, interval = 1)
        
        self.ax.callbacks.connect('xlim_changed', lambda event: self.anim._blit_cache.clear())
        self.ax.callbacks.connect('ylim_changed', lambda event: self.anim._blit_cache.clear())

        if frameByFrame:
            self.anim_running = False
            #self.anim.event_source.stop()
        else:
            self.anim_running = True

        self.reset_forces()
            
        if self.time_profile:
            self.forces_times = {}
            for force in self.forces_to_compute:
                self.forces_times[force.__name__] = 0.0r

            self.update_forces_time = 0.0r
            self.check_pos_correct_time = 0.0r

        if self.nVertices > 30:
            self.useLogSpring = True
        else:
            self.useLogSpring = False

        self.validFrames = 0


    def reset_forces(self):
        self.forces_to_compute = []

        if self.useVVForce:
            self.forces_to_compute.append(self.computeRepulsionVVForces)
        if self.useSpringForce:
            self.forces_to_compute.append(self.computeSpringForces)
        if self.useWeakSpringForce:
            self.forces_to_compute.append(self.computeWeakSpringForces)
        if self.useTorsionForce:
            self.forces_to_compute.append(self.computeTorsionForces)
        if self.useEEForce:
            self.forces_to_compute.append(self.computeRepulsionEEForces)

    def print_timers(self):
        if self.time_profile and self.frame > 0:
            print ("Average correctness check time:", int(self.check_pos_correct_time * 1000 / self.frame))
            print ("Average forces computing time:", int(self.update_forces_time * 1000 / self.frame))

            for force in self.forces_to_compute:
                print ('\t', force.__name__, '\t', int(self.forces_times[force.__name__] * 1000 / self.frame))
            print()

    def onClick(self, event):
        return

        if self.anim_running:
            self.anim.event_source.stop()
            self.anim_running = False
        else:
            self.anim.event_source.start()
            self.anim_running = True

    def onKey(self, event):
        if event.key == ' ':
            self.doFrame = True

        if event.key == "enter":
            if self.anim_running:
                #self.anim.event_source.stop()
                self.anim_running = False
            else:
                #self.anim.event_source.start()
                self.anim_running = True

        if event.key == "p":
            self.print_timers()

        if event.key == "escape":
            plt.close()

    def update_slider_val(self, val):
        self.fixed_delta_t = 10.0r ** val
        self.slider.valtext.set_text(round(self.fixed_delta_t, 3))
        self.force_delta_t = True

    def computeRepulsionVVForces(self):
        for i in range(self.nVertices):
            for j in range(i+1, self.nVertices):
                r = self.pos[j] - self.pos[i]
                if i < self.real_n_vertices and j < self.real_n_vertices:
                    coef = self.repulsionVVRRCoef
                elif i >= self.real_n_vertices and j >= self.real_n_vertices:
                    coef = self.repulsionVVWWCoef
                else:
                    coef = self.repulsionVVRWCoef

                force = coef * r.normalized() / r.normSq()     # force applied to i (going away from j)
                self.forces[i] -= force
                self.forces[j] += force

    def computeSpringForces(self):
        for (i, j) in self.edges:
            r = self.pos[j] - self.pos[i]
            
            if i < self.real_n_vertices and j < self.real_n_vertices:
                length = self.springLength
            else:
                length = self.springLength / self.break_down_num

            if self.useLogSpring:
                force = self.springCoef * math.log(r.norm() / length) * r.normalized()
            else:
                force = self.springCoef * (r.norm() - length) * r.normalized()

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
                    goal = 2.0r * math.pi / len(self.embedding[i])
                    
                    energy = self.torsionCoef * math.tanh(angle / self.torsionScale) ** self.torsionPower

                    #scale = -2.0r * self.torsionPower * 10 * energy / self.torsionScale / math.sinh(2.0r * angle / self.torsionScale)

                    scale = self.torsionCoef * abs(angle - goal) ** self.torsionPower
                    if angle < goal:
                        scale = -scale

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
            for i in range(len(face)):      # iterate over every pair (node, edge)
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
                    
                    angle_sum += 2.0r * math.pi - (pos[u1] - pos[i]).angle_towards(pos[u2] - pos[i])   # embedding is given clockwise

                # if the embedding is correct, angle_sum = 2pi
                # otherwise, it's way bigger (>= 4pi?)

                if angle_sum >= 2.1r * math.pi:  # account for float inaccuracies
                    print ("bad!! bad angle (sum", angle_sum, "!)")
                    # print (pos[i])
                    # for k in self.embedding[i]:
                    #    print (pos[k])
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
                #print (segments)
                return False

        return True

    def update_forces(self):
        self.forces = [Vector2D() for i in range(self.nVertices)]
        
        for force in self.forces_to_compute:
            if self.time_profile:
                begin = time.perf_counter()
            force()
            if self.time_profile:
                end = time.perf_counter()
                self.forces_times[force.__name__] += end - begin
        """
        self.computeRepulsionVVForces()
        self.computeSpringForces()
        self.computeWeakSpringForces()
        self.computeTorsionForces()
        self.computeRepulsionEEForces()
        """
        
    def tick(self):
        full_force_frames = int(self.nVertices * 1.5r)   # after this number of iterations, we should be closer to the optimum; start going slower to refine the positions
        power = 1r                                     # delta_t will be proportional to 1 / frame ^ power when this number is reached

        if self.frame < full_force_frames:
            delta_t = self.fast_delta_t
            iters = 1
        else:
            delta_t = self.slow_delta_t * full_force_frames**power / self.frame**power
            iters = self.max_iter_tick      # we do several iterations with geometrically decreasing delta_t to reach the optimum without flickering

        if self.force_delta_t:
            delta_t = self.fixed_delta_t
        else:
            self.slider.valtext.set_text(round(delta_t, 3))
        
        while iters > 0:
            if self.time_profile:
                debut = time.perf_counter()
            self.update_forces() # LONG ? Et si oui, quelle force ?
            if self.time_profile:
                fin = time.perf_counter() 
                self.update_forces_time += float(fin - debut)

            maxForce = max(map(Vector2D.normSq, self.forces))
            #if maxForce > 0:
            #    delta_t = min(delta_t, 2.0 / maxForce)

            if maxForce * self.slow_delta_t < self.max_extent * self.convergence_limit and not self.done:
                print ("STOP")
                #self.anim.event_source.stop()
                self.anim_running = False
                self.done = True

            newpos = [Vector2D() for i in range(self.nVertices)]

            for i in range(self.nVertices):
                displ = self.forces[i] * delta_t
                if displ.normSq() > self.currentMaxDispl ** 2.0r:
                    displ = displ.normalized() * self.currentMaxDispl
                
                newpos[i] = self.pos[i] + displ

                # considering speed (ie. assuming nodes have inertia) seems to make the convergence slower
                # self.pos[i] += self.speed[i] * self.maxDelta + 0.5 * self.forces[i] * self.maxDelta * self.maxDelta / self.vertexInertia
                # self.speed[i] /= 2              # friction forces (we need energy loss somewhere)
                # self.speed[i] += self.forces[i] * self.maxDelta / self.vertexInertia
            
            if self.time_profile:
                debut = time.perf_counter()

            if self.check_pos_correct(newpos): # LONG ?
                self.pos = newpos
                iters -= 1r
                delta_t /= 2.0r
            
                self.frame += 1r
                self.validFrames += 1r
            else:
                # we're going too fast... reduce delta_t & maxDispl
                delta_t /= 2.0r
                self.currentMaxDispl /= 2.0r

                self.validFrames = 0

            if self.validFrames > self.validFramesThreshold:
                self.currentMaxDispl = 0.2 * self.maxDispl + 0.8 * self.currentMaxDispl
            
            if self.time_profile:
                fin = time.perf_counter()
                self.check_pos_correct_time += float(fin - debut)

    def update_fig(self, frame):
        if frame > 5 and (self.anim_running or self.doFrame):
            self.tick()
            self.doFrame = False
            #self.pos[0].x += (random.random() - 0.5) / 5

        self.text.set_text("Frame: " + str(self.frame) + ("; done" if self.done else ""))

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