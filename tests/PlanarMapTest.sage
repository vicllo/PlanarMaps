from sage.all_cmdline import *   # import sage library


load("LabelledMap.sage")
load("MutableLabelledMap.sage")
load("MapGenerator.sage")


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