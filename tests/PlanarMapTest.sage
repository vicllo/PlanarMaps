from sage.all_cmdline import *   # import sage library


load("../PlanarMap.sage")


#Testing the init function of the PlanarMap class
def test_init_planar_map():

    #Test case 1
    #The case where sigma and alpha doesn't have the same size    
    sigma_1 = Permutation( [3,4,1,2,6,5])
    alpha_1 = Permutation( [(1,2),(3,4)])
   
    same_size_fail = False

    try: 
        map_1 = PlanarMap(sigma_1,alpha_1)
    except ValueError as e:
        if str(e) == "The two permutations does not have the same size" :
            same_size_fail = True
    assert same_size_fail is True

    #Test case 2
    #The case where alpha isn't an involution
    sigma_1 = Permutation( [3,4,1,2,6,5])
    alpha_1 = Permutation( [(1,2),(3,4)])
   
    same_size_fail = False

    try: 
        map_1 = PlanarMap(sigma_1,alpha_1)
    except ValueError as e:
        if str(e) == "The two permutations does not have the same size" :
            same_size_fail = True
    assert same_size_fail is True

    #Test case 3
    #The case where alpha isn't an involution
    sigma_2 = Permutation( [3,4,1,2,5])
    alpha_2 = Permutation( [(1,2),(3,4,5)])
    

    involution_fail = False

    try: 
        map_2 = PlanarMap(sigma_2,alpha_2)
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
        map_3 = PlanarMap(sigma_3,alpha_3)
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
        map_4 = PlanarMap(sigma_4,alpha_4)
    except ValueError as e:
        if str(e) == "The graph isn't connected" :
            connected_fail = True
    assert connected_fail is True

    #Test case 5
    #The case where the graph is correctly defined
    sigma_5 = Permutation( [1,3,2,5,4,6])
    alpha_5 = Permutation( [(1,2),(3,4),(5,6)])
    
    correct_graph = True

    try: 
        map_5 = PlanarMap(sigma_5,alpha_5)   
    except:
        correct_graph = False
    assert correct_graph is True


#Test the repr function of planar map
def test_repr_map():
    #Test case 1
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = PlanarMap(sigma_1,alpha_1)
    correct_repr = "Sigma : [1, 3, 2, 5, 4, 6], Alpha : [2, 1, 4, 3, 6, 5]"
    assert str(map_1) == correct_repr

#Test the numberOfFaces function
def test_number_of_faces():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = PlanarMap(sigma_1,alpha_1)
    assert map_1.numberOfFaces()==1

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = PlanarMap(sigma_2,alpha_2)
    assert map_2.numberOfFaces()==2

#Test the numberOfNodes function
def test_number_of_nodes():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = PlanarMap(sigma_1,alpha_1)
    assert map_1.numberOfNodes()==4

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = PlanarMap(sigma_2,alpha_2)
    assert map_2.numberOfNodes()==3

#Test the numberOfEdges function
def test_number_of_edges():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    map_1 = PlanarMap(sigma_1,alpha_1)
    assert map_1.numberOfEdges()==3

    #Test case 2
    #Correspond to a triangle
    sigma_2 = Permutation( [(1,6),(2,3),(4,5)])
    alpha_2 = Permutation( [(1,2),(3,4),(5,6)])
    map_2 = PlanarMap(sigma_2,alpha_2)
    assert map_2.numberOfEdges()==3

#Test the buildGraph function
def test_build_graph():

    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_1 = PlanarMap(sigma_1,alpha_1)
    
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
    
    map_2 = PlanarMap(sigma_2,alpha_2)
    
    graph_2 = map_2.buildGraph()
    edges_2 = graph_2.edges(labels = False)

    correctEdges_2 = [(1,2),(1,3),(2,3)]
    
    passed_test_2 = True

    for k in range(map_2.numberOfEdges()):
        passed_test_2 &= edges_2[k]==correctEdges_2[k]

    assert passed_test_2 is True

#Return the dual of the Planar map
def test_dual():
    #Test case 1
    #Correspond to a linear tree with 4 nodes
    sigma_1 = Permutation( [1,3,2,5,4,6])
    alpha_1 = Permutation( [(1,2),(3,4),(5,6)])
    
    map_1 = PlanarMap(sigma_1,alpha_1)
    
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
    
    map_2 = PlanarMap(sigma_2,alpha_2)

    dualMap_2 = map_2.dual()

    graphDual_2  = dualMap_2.buildGraph()

    correctDualEdges_2 = [(1,2),(1,2),(1,2)]

    dualEdges_2 = graphDual_2.edges(labels = False)

    passed_test_2 = True

    for k in range(map_2.numberOfEdges()):
        passed_test_2 &= correctDualEdges_2[k]==dualEdges_2[k]
    assert passed_test_2 is True

if __name__ == "__main__":
    #Test init function
    test_init_planar_map()
    print("test_init_planar_map passed")
    
    #Test repr function
    test_repr_map()
    print("test_repr_map passed")

    #Test the numberOfFaces function
    test_number_of_faces()
    print("test_number_of_faces passed")

    #Test the numberOfNodes function
    test_number_of_nodes()
    print("test_number_of_nodes passed")

    #Test the numberOfEdges function
    test_number_of_edges()
    print("test_number_of_edges passed")

    #Test the build graph function
    test_build_graph()
    print("test_build_graph passed")

    #Test the dual graph
    test_dual()
    print("test_dual_graph passed")