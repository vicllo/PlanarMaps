from sage.all_cmdline import *   # import sage library

load("LabelledMap.sage")
load("MutableLabelledMap.sage")
load("MapGenerator.sage")


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
    #Test the getRandomDyckPath method
    test_getRandomDyckPath()
    print("test_getRandomDyckPath passed")