
import os
import sys

relPath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1,relPath)
sys.path.insert(1,relPath + "/modules/")
sys.path.insert(1,relPath + "/modules/nelmin")
sys.path.insert(1,relPath + "/modules/bulletprobe")
sys.path.insert(1,relPath + "/modules/medialaxis")
sys.path.insert(1,relPath + "/modules/alphashape")
sys.path.insert(1,relPath + "/modules/ogre")



from opaque.maps.Paths import getTangentIntersections


inputVals = []
#f = open("testData/tangentData/input.txt", 'r')
f = open("testData/tangentData/input2.txt", 'r')
str_f = f.readline()

while str_f != '':
    inputVals.append(eval(str_f))
    print inputVals[-1][-1]
    str_f = f.readline()
    
plotCount = 0
for entry in inputVals:
    getTangentIntersections(entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], plotCount)
    plotCount += 1















