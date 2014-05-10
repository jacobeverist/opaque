#!/usr/bin/python
import getopt, sys, os, shutil, csv
from pylab import *
#from Numeric import *
from math import *
from copy import *
import gen_icp

import toromod

def writeConstrainedToroGraph(testDir, filename, v_list, e_list):

	filename = testDir + "/" + filename
	sf = open(filename, 'w')


	nodeCount = 0
	poses = {}


	for vertex in v_list:

		nodeID = vertex[0]
		pose = vertex[1]
		
		poses[nodeID] = pose
		
		if nodeID > nodeCount:
			nodeCount = nodeID

		sf.write("VERTEX " + str(nodeID) + " " + str(pose[0]) + " " + str(pose[1]) + " " + str(pose[2]))
		sf.write("\n")
		
	edgeHash = {}
	regularEdge = []
	extraEdge = []

	for edge in e_list:
		node1 = edge[0]
		node2 = edge[1]
		
		try:
			if edgeHash[node1] == node2:
				extraEdge.append(edge)
			else:
				regularEdge.append(edge)
		except:
			regularEdge.append(edge)
			
		edgeHash[node1] = node2
	
	#print "regular edges:"
	#print regularEdge
	##print
	#print "extra edges:"
	#print extraEdge


	" for each extra edge, change the node ID and add an EQUIV statement "
	equivs = []
	finalExtra = []
	extraNodes = []
	for edge in extraEdge:
		node1 = edge[0]
		node2 = edge[1]
		
		nodeCount += 1
		newNodeID = nodeCount
		
		equivs.append([node2, newNodeID])		
		equiv_edge = [node2, newNodeID, [0.0,0.0,0.0], [10,0,10,10,0,0] ]

		newEdge = deepcopy(edge)		
		newEdge[1] = newNodeID
		
		finalExtra.append(newEdge)
		finalExtra.append(equiv_edge)
		
		extraNodes.append([newNodeID, poses[node2]])

	for vertex in extraNodes:

		nodeID = vertex[0]
		pose = vertex[1]

		sf.write("VERTEX " + str(nodeID) + " " + str(pose[0]) + " " + str(pose[1]) + " " + str(pose[2]))
		sf.write("\n")
		


	for edge in regularEdge:
		sf.write("EDGE ")

		node1 = edge[0]
		node2 = edge[1]
		offset = edge[2]
		covar = edge[3]

		sf.write(str(node1) + " ")
		sf.write(str(node2) + " ")
		sf.write(str(offset[0]) + " ")
		sf.write(str(offset[1]) + " ")
		sf.write(str(offset[2]) + " ")
		sf.write(str(covar[0]) + " ")
		sf.write(str(covar[1]) + " ")
		sf.write(str(covar[2]) + " ")
		sf.write(str(covar[3]) + " ")
		sf.write(str(covar[4]) + " ")
		sf.write(str(covar[5]))
		sf.write("\n")

	for edge in finalExtra:
		sf.write("EDGE ")

		node1 = edge[0]
		node2 = edge[1]
		offset = edge[2]
		covar = edge[3]

		sf.write(str(node1) + " ")
		sf.write(str(node2) + " ")
		sf.write(str(offset[0]) + " ")
		sf.write(str(offset[1]) + " ")
		sf.write(str(offset[2]) + " ")
		sf.write(str(covar[0]) + " ")
		sf.write(str(covar[1]) + " ")
		sf.write(str(covar[2]) + " ")
		sf.write(str(covar[3]) + " ")
		sf.write(str(covar[4]) + " ")
		sf.write(str(covar[5]))
		sf.write("\n")

	for equiv in equivs:
		sf.write("EQUIV ")
		
		node1 = equiv[0]
		node2 = equiv[1]
		
		sf.write(str(node1) + " ")
		sf.write(str(node2) + "\n")


	sf.close()

def writeToroGraph(testDir, filename, v_list, e_list, edgeConstraints):

	filename = testDir + "/" + filename
	sf = open(filename, 'w')

	for vertex in v_list:

		nodeID = vertex[0]
		pose = vertex[1]

		sf.write("VERTEX " + str(nodeID) + " " + str(pose[0]) + " " + str(pose[1]) + " " + str(pose[2]))
		sf.write("\n")

	for edge in e_list:
		sf.write("EDGE ")

		node1 = edge[0]
		node2 = edge[1]
		offset = edge[2]
		covar = edge[3]

		sf.write(str(node1) + " ")
		sf.write(str(node2) + " ")
		sf.write(str(offset[0]) + " ")
		sf.write(str(offset[1]) + " ")
		sf.write(str(offset[2]) + " ")
		sf.write(str(covar[0]) + " ")
		sf.write(str(covar[1]) + " ")
		sf.write(str(covar[2]) + " ")
		sf.write(str(covar[3]) + " ")
		sf.write(str(covar[4]) + " ")
		sf.write(str(covar[5]))
		sf.write("\n")

	#constraint = [nodeID1, nodeID2, xDiff, yDiff, rDiff, nomVar, nomVar, nomVar, nomVar, nomVar, nomVar]

	for constraint in edgeConstraints:
		sf.write("EDGE ")
		for item in constraint:
			sf.write(str(item) + " ")
		sf.write("\n")

	sf.close()

def readToroGraph(fileName):

	pointReader = csv.reader(file(fileName), delimiter=' ')

	# read data into lists
	vertexList = []
	edgeList = []

	for row in pointReader:

		#print row
		# add to the vertex list
		if row[0] == "VERTEX" or row[0] == "VERTEX2":

			#vertexList.append([float(row[2]),float(row[3]),float(row[4]),int(row[1])])	

			nodeID = int(row[1])

			# x, y, theta
			pose = [float(row[2]),float(row[3]),float(row[4])]	

			vertex = [nodeID, pose]

			vertexList.append(vertex)


		# add to the edge list
		if row[0] == "EDGE" or row[0] == "EDGE2":
			#edge = [int(row[1]),int(row[2]),float(row[3]),float(row[4]),float(row[5]),
			#	float(row[6]),float(row[7]),float(row[8]),float(row[9]),float(row[10]),float(row[11])]	

			# node 1, node 2
			node1 = int(row[1])
			node2 = int(row[2])

			# forward, slide, rotate
			offset = [float(row[3]), float(row[4]), float(row[5])]
			covar = [float(row[6]),float(row[7]),float(row[8]),float(row[9]),float(row[10]),float(row[11])]	

			edge = [ node1, node2, offset, covar]
			edgeList.append(edge)	

	#for i in range(0,len(vertexList)):
	#	if i != vertexList[i][0]:
	#		print "Nodes in TORO graph out of order"
	#		raise

	return vertexList, edgeList

def plotToro(vertexList, edgeList = [], hulls = [], drawedge=False, clr='0.8'):

	vx = []
	vy = []
	for item in vertexList:
		pose1 = item[1]
		vx.append(pose1[0])
		vy.append(pose1[1])

	if drawedge:
		ex1 = []
		ey1 = []
		ex2 = []
		ey2 = []
		for item in edgeList:
			n1 = item[0]
			n2 = item[1]
			
			#print "%d -> %d" % (n1, n2)
			pose1 = vertexList[n1][1]
			pose2 = vertexList[n2][1]
			ex1.append(pose1[0])
			ey1.append(pose1[1])
			ex2.append(pose2[0])
			ey2.append(pose2[1])

		for index in range(len(ex1)):
			xSeg = [ ex1[index], ex2[index] ]
			ySeg = [ ey1[index], ey2[index] ]
			plot(xSeg,ySeg, linewidth=0.2, color='0.5')

			#plot(ySeg,xSeg, linewidth=0.2, color='0.5')

	scatter(vx,vy,color=clr,faceted=False)
	#scatter(vy,vx,color=clr,faceted=False)

	
	if len(hulls) > 0:
		for i in range(len(vertexList)):
			pose1 = vertexList[i][1]
	
	
			hull1 = hulls[i]
			
			hull_trans = []
			for p in hull1:
				hull_trans.append(gen_icp.dispPoint(p, pose1))	
							
			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			xP.append(hull_trans[0][0])
			yP.append(hull_trans[0][1])
			plot(xP,yP)		

	#xlim(-16,16)
	#ylim(-16,16)
	#xlim(-5,5)
	#ylim(-2,7)
	xlim(-8,8)
	ylim(-5,10)

def runTORO(v_list, e_list):

	v_list2, e_list2 = toromod.doToro(v_list,e_list)

	return v_list2, e_list2


def executeToro(fileName):
	
	# 1. make temporary work directory
	# 2. ch working directory
	# 3. execute TORO on filename
	# 4. cleanup

	fileRoot = copy(fileName)
	fileRoot = fileRoot[:-6]

	if sys.platform == "win32":
		os.system("mkdir toro_tmp")
		os.system("cp " + fileName + " toro_tmp")
		os.system("cp toro.exe toro_tmp")
		os.chdir("toro_tmp")
		result = os.system("toro -nde " + fileName + " 2> toro.out.txt")
		
	else:
		os.system("mkdir toro_tmp")
		os.system("cp " + fileName + " toro_tmp")
		os.system("cp toro toro_tmp")
		#os.system("cp toro.exe toro_tmp")
		os.chdir("toro_tmp")
		result = os.system("./toro -nde " + fileName + " 2> toro.out.txt")

	if result != 0 and result != 1:
		print "returned", result
		raise
	os.chdir("..")
	os.system("rm " + fileRoot + "-treeopt-final.graph")
	os.system("mv toro_tmp/" + fileRoot + "-treeopt-final.graph .")
	#os.system("rm toro_tmp/*")
	#os.system("rmdir toro_tmp")

	#os.system("./toro ../toro_tmp/" + fileName)
	#os.chdir("toro_tmp")
	#os.system("../toro ../" + fileName)
	#os.chdir("..")
	#os.system("mv toro_tmp/" + fileRoot + "-treeopt-final.graph .")
	#os.system("rm toro_tmp/*")
	#os.system("rmdir toro_tmp")


if __name__ == '__main__':

	#testDir = "test_0001"
	fileName = "test.graph"
	#fileName = "w10000-odom.graph"
	vlist, elist = readToroGraph("./" + fileName)
	print len(vlist), len(elist)
	plotToro(vlist, elist, drawedge=True)
	savefig("w10000_uncorrect.png")
	clf()

	executeToro("./" + fileName)

	finalFileName = "test-treeopt-final.graph"
	#finalFileName = "w10000-odom-treeopt-final.graph"
	vlist, elist = readToroGraph("./" + finalFileName)
	#print len(vlist), len(elist)
	plotToro(vlist, elist, drawedge=True)
	savefig("w10000_correct.png")

	#executeToro(fileName)

	#vlist, elist = readToroGraph(testDir + "/" + fileName)
	#writeToroGraph(".", vlist, elist, [])
	#show()

