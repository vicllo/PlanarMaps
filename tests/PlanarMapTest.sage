from sage.all_cmdline import *   # import sage library


load("LabelledMap.sage")
load("MutableLabelledMap.sage")
load("MapGenerator.sage")


#Testing the init method of the LabelledMap class
def test_init_planar_map():

    #Test case 1
    #The case where sigma and alpha don't have the same size    
    sigma_1 = Permutation( [3,4,1,2,6,5])
    alpha_1 = Permutation( [(1,2),(3,4)])
   
    try: 
        LabelledMap(sigma_1,alpha_1)
    except ValueError as e:
        pass
    else:
        assert False

    #Test case 2
    #The case where alpha isn't an involution
    sigma_2 = Permutation( [3,4,1,2,5])
    alpha_2 = Permutation( [(1,2),(3,4,5)])
    

    involution_fail = False

    try: 
        map_2 = LabelledMap(sigma_2,alpha_2)
    except ValueError as e:
        if str(e) == "The permutation alpha is not an involution" :
            involution_fail = True
    assert involution_fail is True

    #Test case 3
    #The case where alpha has a fixed point
    sigma_3 = Permutation( [3,4,1,2,5])
    alpha_3 = Permutation( [2,1,3,5,4])
    
    fixed_point_fail = False

    try: 
        map_3 = LabelledMap(sigma_3,alpha_3)
    except ValueError as e:
        if str(e) == "The permutation alpha should not have fixed points" :
            fixed_point_fail = True
    assert fixed_point_fail is True

    #Test case 4
    #Correspond to a linear tree with 4 nodes
    #The case where the graph isn't connected
    sigma_4 = Permutation( [1,2,3,5,4,6])
    alpha_4 = Permutation( [(1,2),(3,4),(5,6)])
    
    connected_fail = False

    try: 
        map_4 = LabelledMap(sigma_4,alpha_4)
    except ValueError:
        connected_fail = True
    assert connected_fail is True

    #Test case 5
    #The case where the graph is correctly defined
    sigma_5 = Permutation( [1,3,2,5,4,6])
    alpha_5 = Permutation( [(1,2),(3,4),(5,6)])
    
    correct_graph = True

    try: 
        map_5 = LabelledMap(sigma_5,alpha_5)   
    except:
        correct_graph = False
    assert correct_graph is True

    # Test case 6
    # Invalid adjacency list
    adj_1 = [(3,), (1,3), (2,)]
    invalid_adj_fail = False

    try:
        map_6 = LabelledMap(adj = adj_1)
    except ValueError:
        invalid_adj_fail = True
    assert invalid_adj_fail

    # Test case 7
    # Too much information
    invalid_args_fail = False

    try:
        map_7 = LabelledMap(sigma_5, adj = adj_1)
    except ValueError:
        invalid_args_fail = True
    assert invalid_args_fail

#Test the repr method of planar map
def test_repr_map():
    #Test case 1
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)
    correct_repr = "Labelled map | Sigma : [1, 3, 2, 5, 4, 6], Alpha : [2, 1, 4, 3, 6, 5]"
    print(str(map_1))
    print(correct_repr)
    assert str(map_1) == correct_repr

#Test the numberOfFaces method
def test_number_of_faces():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)
    assert map_1.numberOfFaces()==1

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = LabelledMap(sigma_2,alpha_2)
    assert map_2.numberOfFaces()==2

    #Test case 3
    #Correspond to a segment
    sigma_3 = Permutation( [1,2])
    alpha_3 = Permutation( [(1,2)])
    map_3 = LabelledMap(sigma_3,alpha_3)
    
    assert map_3.numberOfFaces() == 1

    #Test case 4
    #Correspond to a linear tree of 3 nodes
    sigma_4 = Permutation( [1,3,2,4])
    alpha_4 = Permutation( [(1,2),(3,4)])
    map_4 = LabelledMap(sigma_4,alpha_4)

    assert map_4.numberOfFaces() == 1

    #Test case 5
    #Correspond to a loop on a node
    sigma_5 = Permutation( [2,1])
    alpha_5 =Permutation( [(1,2)])
    map_5 = LabelledMap(sigma_5,alpha_5)
    
    assert map_5.numberOfFaces() == 2

    #Test case 6
    #Correspond to two nodes link by 3 edges
    sigma_6 = Permutation( [(1,3,5),(2,6,4)])
    alpha_6 =Permutation( [(1,2),(3,4),(5,6)])
    map_6 = LabelledMap(sigma_6,alpha_6)
    
    assert map_6.numberOfFaces() == 3

    
    #Test case 7
    #Correspond to two nodes link with 4 edges
    sigma_7 = Permutation( [(1,7,3,5),(2,6,4,8)])
    alpha_7 = Permutation( [(1,2),(3,4),(5,6),(7,8)])
    map_7 = LabelledMap(sigma_7,alpha_7)
    
    assert map_7.numberOfFaces() == 4
    

#Test the numberOfNodes method
def test_number_of_nodes():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)
    assert map_1.numberOfNodes()==4

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = LabelledMap(sigma_2,alpha_2)
    assert map_2.numberOfNodes()==3

    #Test case 3
    #Correspond to a segment
    sigma_3 = Permutation( [1,2])
    alpha_3 = Permutation( [(1,2)])
    map_3 = LabelledMap(sigma_3,alpha_3)
    
    assert map_3.numberOfNodes() == 2

    #Test case 4
    #Correspond to a linear tree of 3 nodes
    sigma_4 = Permutation( [1,3,2,4])
    alpha_4 = Permutation( [(1,2),(3,4)])
    map_4 = LabelledMap(sigma_4,alpha_4)

    assert map_4.numberOfNodes() == 3

    #Test case 5
    #Correspond to a loop on a node
    sigma_5 = Permutation( [2,1])
    alpha_5 =Permutation( [(1,2)])
    map_5 = LabelledMap(sigma_5,alpha_5)
    
    assert map_5.numberOfNodes() == 1

    #Test case 6
    #Correspond to two nodes link by 3 edges
    sigma_6 = Permutation( [(1,3,5),(2,6,4)])
    alpha_6 =Permutation( [(1,2),(3,4),(5,6)])
    map_6 = LabelledMap(sigma_6,alpha_6)
    
    assert map_6.numberOfNodes() == 2

    
    #Test case 7
    #Correspond to two nodes link with 4 edges
    sigma_7 = Permutation( [(1,7,3,5),(2,6,4,8)])
    alpha_7 = Permutation( [(1,2),(3,4),(5,6),(7,8)])
    map_7 = LabelledMap(sigma_7,alpha_7)
    
    assert map_7.numberOfNodes() == 2

#Test the numberOfEdges method
def test_number_of_edges():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)
    assert map_1.numberOfEdges()==3

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = LabelledMap(sigma_2,alpha_2)
    
    assert map_2.numberOfEdges()==3

    #Test case 3
    #Correspond to a segment
    sigma_3 = Permutation( [1,2])
    alpha_3 = Permutation( [(1,2)])
    map_3 = LabelledMap(sigma_3,alpha_3)
    
    assert map_3.numberOfEdges() == 1

    #Test case 4
    #Correspond to a linear tree of 3 nodes
    sigma_4 = Permutation( [1,3,2,4])
    alpha_4 = Permutation( [(1,2),(3,4)])
    map_4 = LabelledMap(sigma_4,alpha_4)

    assert map_4.numberOfEdges() == 2

    #Test case 5
    #Correspond to a loop on a node
    sigma_5 = Permutation( [2,1])
    alpha_5 =Permutation( [(1,2)])
    map_5 = LabelledMap(sigma_5,alpha_5)
    
    assert map_5.numberOfEdges() == 1

    #Test case 6
    #Correspond to two nodes link by 3 edges
    sigma_6 = Permutation( [(1,3,5),(2,6,4)])
    alpha_6 =Permutation( [(1,2),(3,4),(5,6)])
    map_6 = LabelledMap(sigma_6,alpha_6)
    
    assert map_6.numberOfEdges() == 3

    
    #Test case 7
    #Correspond to two nodes link with 4 edges
    sigma_7 = Permutation( [(1,7,3,5),(2,6,4,8)])
    alpha_7 = Permutation( [(1,2),(3,4),(5,6),(7,8)])
    map_7 = LabelledMap(sigma_7,alpha_7)
    
    assert map_7.numberOfEdges() == 4

#Test the buildGraph method
def test_build_graph():

    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_1 = LabelledMap(sigma_1,alpha_1)
    
    graph_1= map_1.buildGraph()
    edges_1 = graph_1.edges(labels = False)
    
    correctEdges_1  = [(1,2),(2,3),(3,4)]

    passed_test_1 = True
    for k in range(map_1.numberOfEdges()):
        passed_test_1 &= edges_1[k]==correctEdges_1[k]

    assert passed_test_1 is True

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_2 = LabelledMap(sigma_2,alpha_2)
    
    graph_2 = map_2.buildGraph()
    edges_2 = graph_2.edges(labels = False)

    correctEdges_2 = [(1,2),(1,3),(2,3)]
    
    passed_test_2 = True

    for k in range(map_2.numberOfEdges()):
        passed_test_2 &= edges_2[k]==correctEdges_2[k]

    assert passed_test_2 is True

#Test the genus method
def test_genus():
    # Test case 1: basic planar graph
    assert LabelledMap(Permutation([(2,3,1,5),(4,),(6,)]),Permutation([(1,2),(3,4),(5,6)])).genus() == 0

    # Test case 2: cube graph
    assert LabelledMap(adj = [(5,4,2),(1,3,6),(4,7,2),(8,3,1),(8,1,6),(5,2,7),(3,8,6),(7,4,5)]).genus() == 0

    # Test case 2: complete bipartite graph with 2 * 3 nodes
    assert LabelledMap(adj = [(4,5,6),(4,5,6),(4,5,6),(1,2,3),(1,2,3),(1,2,3)]).genus() == 1

#Test the contract_edge method
def test_contract_edge():
    cube = MutableLabelledMap(adj = [(5,4,2),(1,3,6),(4,7,2),(8,3,1),(8,1,6),(5,2,7),(3,8,6),(7,4,5)])
    
    cube.contractEdge(1)
    cube.contractEdge(3)
    
    assert cube.numberOfEdges() == 10
    assert cube.numberOfNodes() == 6

#Test the dual method
def test_dual():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_1 = LabelledMap(sigma_1,alpha_1)
    
    dualMap_1 = map_1.dual()


    graphDual_1 = dualMap_1.buildGraph()

    correctDualEdges_1 = [(1,1),(1,1),(1,1)]

    dualEdges_1 = graphDual_1.edges(labels = False)

    passed_test_1 = True

    for k in range(map_1.numberOfEdges()):
        passed_test_1 &= dualEdges_1[k]==correctDualEdges_1[k]
    assert passed_test_1 is True

    #Test case 2    
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_2 = LabelledMap(sigma_2,alpha_2)

    dualMap_2 = map_2.dual()

    graphDual_2  = dualMap_2.buildGraph()

    correctDualEdges_2 = [(1,2),(1,2),(1,2)]

    dualEdges_2 = graphDual_2.edges(labels = False)

    passed_test_2 = True

    for k in range(map_2.numberOfEdges()):
        passed_test_2 &= correctDualEdges_2[k]==dualEdges_2[k]
    assert passed_test_2 is True

#Test the diameter method
def test_diameter():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_1 = LabelledMap(sigma_1,alpha_1)

    assert map_1.diameter() == 3

    #Test case 2    
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_2 = LabelledMap(sigma_2,alpha_2)

    assert map_2.diameter() == 1



def printVEF(pMap):
    print("V: ",pMap.numberOfNodes(),"E: ",pMap.numberOfEdges(),"F: ",pMap.numberOfFaces())

def printSigmaAlphaPhi(pMap):
    print("Alpha:",pMap.alpha.to_cycles())
    print("Sigma:",pMap.alpha.to_cycles())
    print("Phi",pMap.phi.to_cycles())

#Test the derivedMap method
def test_derivedMap():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)
    derived_1 = map_1.derivedMap()
    assert derived_1.numberOfFaces() == 6
    assert derived_1.numberOfEdges() == 12
    assert derived_1.numberOfNodes() == 8

    #Test case 2    
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = LabelledMap(sigma_2,alpha_2)
    derived_2 = map_2.derivedMap()
    assert derived_2.numberOfFaces() == 6
    assert derived_2.numberOfEdges() == 12
    assert derived_2.numberOfNodes() == 8

    #Test case 3
    #Correspond to a segment
    sigma_3 = Permutation( [1,2])
    alpha_3 = Permutation( [(1,2)])
    map_3 = LabelledMap(sigma_3,alpha_3)
    derived_3 = map_3.derivedMap()
    assert derived_3.numberOfFaces() == 2
    assert derived_3.numberOfEdges() == 4
    assert derived_3.numberOfNodes() == 4

    #Test case 4
    #Correspond to a linear tree of 3 nodes
    sigma_4 = Permutation( [1,3,2,4])
    alpha_4 = Permutation( [(1,2),(3,4)])
    map_4 = LabelledMap(sigma_4,alpha_4)
    derived_4 = map_4.derivedMap()

    assert derived_4.numberOfNodes() == 6
    assert derived_4.numberOfEdges() == 8
    assert derived_4.numberOfFaces() == 4

    #Test case 5
    #Correspond to a loop on a node
    sigma_5 = Permutation( [2,1])
    alpha_5 =Permutation( [(1,2)])
    map_5 = LabelledMap(sigma_5,alpha_5)
    derived_5 = map_5.derivedMap()
    
    assert derived_5.numberOfNodes() == 4
    assert derived_5.numberOfEdges() == 4
    assert derived_5.numberOfFaces() == 2

    #Test case 6
    #Correspond to two nodes link by 3 edges
    sigma_6 = Permutation( [(1,3,5),(2,6,4)])
    alpha_6 =Permutation( [(1,2),(3,4),(5,6)])
    map_6 = LabelledMap(sigma_6,alpha_6)
    derived_6 = map_6.derivedMap()
    assert derived_6.numberOfNodes() == 8
    assert derived_6.numberOfEdges() == 12
    assert derived_6.numberOfFaces() == 6

    #Test case 7
    #Correspond to two nodes link with 4 edges
    sigma_7 = Permutation( [(1,7,3,5),(2,6,4,8)])
    alpha_7 = Permutation( [(1,2),(3,4),(5,6),(7,8)])
    map_7 = LabelledMap(sigma_7,alpha_7)
    derived_7 = map_7.derivedMap()

    assert derived_7.numberOfFaces() == 8
    assert derived_7.numberOfEdges() == 16
    assert derived_7.numberOfNodes() == 10

#Test the incidence map
def test_incidenceMap():

    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [_sage_const_1 ,_sage_const_3 ,_sage_const_2 ,_sage_const_5 ,_sage_const_4 ,_sage_const_6 ])
    alpha_1 = Permutation( [(_sage_const_1 ,_sage_const_2 ),(_sage_const_3 ,_sage_const_4 ),(_sage_const_5 ,_sage_const_6 )])
    map_1 = LabelledMap(sigma_1,alpha_1)
    
    incidence_1 = map_1.incidenceMap()


    assert incidence_1.numberOfNodes() == 5
    assert incidence_1.numberOfFaces() == 3
    assert incidence_1.numberOfEdges() == 6

    #Test case 2    
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = LabelledMap(sigma_2,alpha_2)

    incidence_2 = map_2.incidenceMap()


    assert incidence_2.numberOfNodes() == 5
    assert incidence_2.numberOfEdges() == 6
    assert incidence_2.numberOfFaces() == 3


    #Test case 3
    #Correspond to a segment
    sigma_3 = Permutation( [1,2])
    alpha_3 = Permutation( [(1,2)])
    map_3 = LabelledMap(sigma_3,alpha_3)

    incidence_3 = map_3.incidenceMap()

    assert incidence_3.numberOfNodes() == 3
    assert incidence_3.numberOfEdges() == 2
    assert incidence_3.numberOfFaces() == 1

    #Test case 4
    #Correspond to a linear tree of 3 nodes
    sigma_4 = Permutation( [1,3,2,4])
    alpha_4 = Permutation( [(1,2),(3,4)])
    map_4 = LabelledMap(sigma_4,alpha_4)
    incidence_4 = map_4.incidenceMap()


    assert incidence_4.numberOfNodes() == 4
    assert incidence_4.numberOfEdges() == 4
    assert incidence_4.numberOfFaces() == 2


    #Test case 5
    #Correspond to a loop on a node
    sigma_5 = Permutation( [2,1])
    alpha_5 =Permutation( [(1,2)])
    map_5 = LabelledMap(sigma_5,alpha_5)
    
    incidence_5 = map_5.incidenceMap()

    assert incidence_5.numberOfNodes() == 3
    assert incidence_5.numberOfEdges() == 2
    assert incidence_5.numberOfFaces() == 1

    #Test case 6
    #Correspond to two nodes link by 3 edges
    sigma_6 = Permutation( [(1,3,5),(2,6,4)])
    alpha_6 =Permutation( [(1,2),(3,4),(5,6)])
    map_6 = LabelledMap(sigma_6,alpha_6)
    
    incidence_6 = map_6.incidenceMap()

    assert incidence_6.numberOfNodes() == 5
    assert incidence_6.numberOfEdges() == 6
    assert incidence_6.numberOfFaces() == 3

    #Test case 7
    #Correspond to two nodes link with 4 edges
    sigma_7 = Permutation( [(1,7,3,5),(2,6,4,8)])
    alpha_7 = Permutation( [(1,2),(3,4),(5,6),(7,8)])
    map_7 = LabelledMap(sigma_7,alpha_7)
    
    incidence_7 = map_7.incidenceMap()

    assert incidence_7.numberOfNodes() == 6
    assert incidence_7.numberOfEdges() == 8
    assert incidence_7.numberOfFaces() == 4

#Test the edge map
def test_edgeMap():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)

    edgeMap_1 = map_1.edgeMap()

    assert edgeMap_1.numberOfNodes() == 3
    assert edgeMap_1.numberOfFaces() == 5
    assert edgeMap_1.numberOfEdges() == 6

    #Test case 2    
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = LabelledMap(sigma_2,alpha_2)

    edgeMap_2 = map_2.edgeMap()

    assert edgeMap_2.numberOfNodes() == 3
    assert edgeMap_2.numberOfEdges() == 6
    assert edgeMap_2.numberOfFaces() == 5


    #Test case 3
    #Correspond to a segment
    sigma_3 = Permutation( [1,2])
    alpha_3 = Permutation( [(1,2)])
    map_3 = LabelledMap(sigma_3,alpha_3)
    
    edgeMap_3 = map_3.edgeMap()

    assert edgeMap_3.numberOfNodes() == 1
    assert edgeMap_3.numberOfEdges() == 2
    assert edgeMap_3.numberOfFaces() == 3

    #Test case 4
    #Correspond to a linear tree of 3 nodes
    sigma_4 = Permutation( [1,3,2,4])
    alpha_4 = Permutation( [(1,2),(3,4)])
    map_4 = LabelledMap(sigma_4,alpha_4)
    edgeMap_4 = map_4.edgeMap()

    assert edgeMap_4.numberOfNodes() == 2
    assert edgeMap_4.numberOfEdges() == 4
    assert edgeMap_4.numberOfFaces() == 4


    #Test case 5
    #Correspond to a loop on a node
    sigma_5 = Permutation( [2,1])
    alpha_5 =Permutation( [(1,2)])
    map_5 = LabelledMap(sigma_5,alpha_5)
    
    edgeMap_5 = map_5.edgeMap()

    assert edgeMap_5.numberOfNodes() == 1
    assert edgeMap_5.numberOfEdges() == 2
    assert edgeMap_5.numberOfFaces() == 3

    #Test case 6
    #Correspond to two nodes link by 3 edges
    sigma_6 = Permutation( [(1,3,5),(2,6,4)])
    alpha_6 =Permutation( [(1,2),(3,4),(5,6)])
    map_6 = LabelledMap(sigma_6,alpha_6)
    
    edgeMap_6 = map_6.edgeMap()

    assert edgeMap_6.numberOfNodes() == 3
    assert edgeMap_6.numberOfEdges() == 6
    assert edgeMap_6.numberOfFaces() == 5

    #Test case 7
    #Correspond to two nodes link with 4 edges
    sigma_7 = Permutation( [(1,7,3,5),(2,6,4,8)])
    alpha_7 = Permutation( [(1,2),(3,4),(5,6),(7,8)])
    map_7 = LabelledMap(sigma_7,alpha_7)
    
    edgeMap_7 = map_7.edgeMap()

    assert edgeMap_7.numberOfNodes() == 4
    assert edgeMap_7.numberOfEdges() == 8
    assert edgeMap_7.numberOfFaces() == 6

#Test the getRootedMapCorrespondance method
def test_getRootedMapCorrespondace():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)

    tau_1 = Permutation((1,3))
    relabelMap_1 = map_1.relabel(tau_1)

    map_1.getRootedMapCorrespondance(relabelMap_1,2) == tau_1

#Test the relabel method
def test_relabel():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = LabelledMap(sigma_1,alpha_1)

    tau_1 = Permutation((1,3))

    relabelMap_1 = map_1.relabel(tau_1)

    assert relabelMap_1.sigma == Permutation([2,1,3,5,4,6])
    assert relabelMap_1.alpha == Permutation([4,3,2,1,6,5])

#Test the getRandomDyckPath method
def test_getRandomDyckPath():
    dyckPath = MapGenerator().getRandomDyckPath(50)
    negative = False 
    level = 0
    for step in dyckPath:
        level += step
        if level < 0: 
            negative = True
    assert (negative is False) and (level == 0)  

if __name__ == "__main__":
    #Test init method
    test_init_planar_map()
    print("test_init_planar_map passed")
    
    #Test repr method
    test_repr_map()
    print("test_repr_map passed")

    #Test the numberOfFaces method
    test_number_of_faces()
    print("test_number_of_faces passed")

    #Test the numberOfNodes method
    test_number_of_nodes()
    print("test_number_of_nodes passed")

    #Test the numberOfEdges method
    test_number_of_edges()
    print("test_number_of_edges passed")

    #Test the buildGraph method
    test_build_graph()
    print("test_build_graph passed")

    #Test the genus method
    test_genus()
    print("test_genus passed")

    #Test the contractEdge methid
    test_contract_edge()
    print("test_contract_edge passed")

    #Test the dual method
    test_dual()
    print("test_dual_graph passed")

    #Test the diameter method
    test_diameter()
    print("test_diameter passed")

    #Test the derivedMap method
    test_derivedMap()
    print("test_derivedMap passed")

    #Test the incidenceMap method
    test_incidenceMap()
    print("test_incidence_map passed")

    #Test the edgeMap method
    test_edgeMap()
    print("test_edgeMap passed")

    #Test the relabel method
    test_relabel()
    print("test_relabel passed")

    #Test the getRootedMapCorrespondance method
    test_getRootedMapCorrespondace()
    print("test_getRootedMapCorrespondance passed")

    #Test the getRandomDyckPath method
    test_getRandomDyckPath()
    print("test_getRandomDyckPath passed")