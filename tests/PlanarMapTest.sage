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