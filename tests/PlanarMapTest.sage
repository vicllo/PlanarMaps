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