
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

from opaque.maps.MapState import getBranchPoint


class DataObj:
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/mapStateSave_%04u_%04u.txt" % (0, numNodes+1)
		f = open(dirName + "/mapStateSave_%04u_%04u.txt" % (0, numNodes+1), 'r')		 
		saveStr = f.read()
		#print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		exec(saveStr)


def recoverData(startIndex, endIndex):

	foo = DataObj()

	dataSet = []
	for k in range(startIndex, endIndex):
		foo.restoreState(".", k)
		dataSet.append((foo.pathClasses, foo.paths))

	saveFile = ""
	saveFile = "dataSet = " + repr(dataSet) + "\n"

	f = open("intersectData1.dat", 'w')
	f.write(saveFile)
	f.close()

def loadData(fileName):
	f = open(fileName, 'r')
	saveStr = f.read()
	f.close()
	saveStr = saveStr.replace('\r\n','\n')
	exec(saveStr)
	
	return dataSet

#recoverData(14, 47)


dataSet = loadData("intersectData1.dat")


for dataEntry in dataSet:

	pathClasses, paths = dataEntry

	for pathID in pathClasses.keys():
		if pathID != 0:

			globalJunctionPose = pathClasses[pathID]["globalJunctionPose"]
			parentID = pathClasses[pathID]["parentID"]
			branchNodeID = pathClasses[pathID]["branchNodeID"]

			parentPath = paths[parentID]
			childPath = paths[pathID]

			newGlobJuncPose, controlPoint = getBranchPoint(globalJunctionPose, parentID, pathID, parentPath, childPath, plotIter = True, hypothesisID = 0, nodeID = branchNodeID)
			newGlobJuncPose2, controlPoint2 = getBranchPoint(globalJunctionPose, pathID, parentID, childPath, parentPath, plotIter = True, hypothesisID = 0, nodeID = branchNodeID)
			
			exit()
		

	#self.pathClasses[newPathID] = {"parentID" : parentID, "branchNodeID" : branchNodeID, "localJunctionPose" : localJunctionPose, "sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : newGlobJuncPose, "controlPoint" : controlPoint }		

#newGlobJuncPose, controlPoint = getBranchPoint(globalJunctionPose, parentID, pathID, parentPath, childPath, plotIter = True, hypothesisID = 0, nodeID = branchNodeID)
#newGlobJuncPose2, controlPoint2 = getBranchPoint(globalJunctionPose, newPathID, parentID, globalMedial0, self.trimmedPaths[parentID], plotIter = True, hypothesisID = self.hypothesisID, nodeID = branchNodeID)



