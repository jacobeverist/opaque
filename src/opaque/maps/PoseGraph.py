import sys

import scipy.linalg
import random

import numpy
import math
from LocalNode import getLongestPath
from StableCurve import StableCurve
from SplineFit import SplineFit
from Pose import Pose
import gen_icp
from functions import *
import toro
import scsgp
import bayes
from operator import itemgetter

import pylab
import Image

from subprocess import Popen, PIPE
from medialaxis import computeMedialAxis
import graph

renderCount = 0
topCount = 0

GND_PRIORITY = 10
PATH_PRIORITY = 6
INPLACE_PRIORITY = 5
CORNER_PRIORITY = 3
SHAPE_PRIORITY = 2
OVERLAP_PRIORITY = 1
MOTION_PRIORITY = 0
INPLACE_PRIORITY2 = -1
INPLACE_PRIORITY3 = -2


def computeBareHull(node1, sweep = False, static = False):
	
	if static:
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		a_data = decimatePoints(a_data)

		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
	else:
				
		" Read in data of Alpha-Shapes without their associated covariances "
		node1.computeAlphaBoundary(sweep = sweep)
		a_data = node1.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
	
	return a_data_GPAC

def computeHull(node1, sweep = False, static = False):
	
	if static:
		
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		a_data = decimatePoints(a_data)

		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		" treat the points with the point-to-line constraint "
		gen_icp.addPointToLineCovariance(a_data_GPAC, high_var=1.0, low_var=0.001)
	
		" treat the points with the distance-from-origin increasing error constraint "
		gen_icp.addDistanceFromOriginCovariance(a_data_GPAC, tan_var=0.1, perp_var=0.01)

	else:			
		" Read in data of Alpha-Shapes and add their associated covariances "
		node1.computeAlphaBoundary(sweep = sweep)
		a_data = node1.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		" treat the points with the point-to-line constraint "
		gen_icp.addPointToLineCovariance(a_data_GPAC, high_var=1.0, low_var=0.001)
	
		" treat the points with the distance-from-origin increasing error constraint "
		gen_icp.addDistanceFromOriginCovariance(a_data_GPAC, tan_var=0.1, perp_var=0.01)
	
	return a_data_GPAC

def computeBoundary(node1, sweep = False):
	
	" Read in data of Alpha-Shapes and add their associated covariances "
	node1.computeAlphaBoundary(sweep)
	#node1.boundaryMap.getBoundaryPoints()
	#computeAlphaBoundary(sweep = sweep)
	a_data = node1.getAlphaBoundary(sweep = sweep)
	a_data = decimatePoints(a_data)
	
	" convert hull points to GPAC coordinates before adding covariances "
	localGPACPose = node1.getLocalGPACPose()
	localGPACProfile = Pose(localGPACPose)
	
	a_data_GPAC = []
	for pnt in a_data:
		a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
	
	" treat the points with the point-to-line constraint "
	gen_icp.addPointToLineCovariance(a_data_GPAC, high_var=1.0, low_var=0.001)

	" treat the points with the distance-from-origin increasing error constraint "
	gen_icp.addDistanceFromOriginCovariance(a_data_GPAC, tan_var=0.1, perp_var=0.01)
	
	return a_data_GPAC

def computeHullAxis(nodeID, node2, tailCutOff = False):

	medial2 = node2.getBestMedialAxis()

	if tailCutOff:
		medial2 = node2.medialTailCuts[0]
	else:
		medial2 = node2.medialLongPaths[0]


	if node2.isBowtie:			
		hull2 = computeBareHull(node2, sweep = False, static = True)
		hull2.append(hull2[0])

	else:
		hull2 = computeBareHull(node2, sweep = False)
		hull2.append(hull2[0])
	
	
	return hull2, medial2

	"""
	if node2.isBowtie:			
		hull2 = computeBareHull(node2, sweep = False, static = True)
		hull2.append(hull2[0])
		medial2 = node2.getStaticMedialAxis()

	else:
		hull2 = computeBareHull(node2, sweep = False)
		hull2.append(hull2[0])
		medial2 = node2.getMedialAxis(sweep = False)
	"""
	
	#print "len(medial2) =", len(medial2)
	#print "medial2:", medial2


	" take the long length segments at tips of medial axis"
	edge1 = medial2[0:2]
	edge2 = medial2[-2:]
	
	frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
	backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
	
	frontVec[0] /= frontMag
	frontVec[1] /= frontMag
	backVec[0] /= backMag
	backVec[1] /= backMag
	
	" make a smaller version of these edges "
	newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
	newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

	edge1 = [newP1, edge1[1]]
	edge2 = [edge2[0], newP2]

	
	" find the intersection points with the hull "
	interPoints = []
	for k in range(len(hull2)-1):
		hullEdge = [hull2[k],hull2[k+1]]
		isIntersect1, point1 = Intersect(edge1, hullEdge)
		if isIntersect1:
			interPoints.append(point1)
			break

	for k in range(len(hull2)-1):
		hullEdge = [hull2[k],hull2[k+1]]
		isIntersect2, point2 = Intersect(edge2, hullEdge)
		if isIntersect2:
			interPoints.append(point2)
			break
	
	" replace the extended edges with a termination point at the hull edge "			
	medial2 = medial2[1:-2]
	if isIntersect1:
		medial2.insert(0, point1)
	if isIntersect2:
		medial2.append(point2)
	

	" cut off the tail of the non-sweeping side "
	TAILDIST = 0.5

	if tailCutOff:
		
		if nodeID % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]
	
		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
			
	return hull2, medial2

class PoseGraph:

	def __init__(self, probe, contacts):
		
		self.probe = probe
		self.contacts = contacts
		
		self.initPose = self.probe.getActualJointPose(19)
		self.nodeHash = {}
		self.edgeHash = {}
		
		self.numNodes = 0
		self.currNode = 0

		self.inplace1 = []
		self.inplace2 = []
		self.inplace3 = []

		self.gnd_constraints = []		
		self.merged_constraints = []
		
		self.sweepHulls = []
		self.nonSweepHulls = []
		
		self.edgePriorityHash = {}
		self.cornerBins = []
		self.a_medial = []

		self.E_gnd = matrix([[ 0.00001, 0.0, 0.0 ],
							[ 0.0, 0.00001, 0.0],
							[ 0.0, 0.0, 0.0001 ]])

		self.E_junction = matrix([[ 0.0001, 0.0, 0.0 ],
							[ 0.0, 0.0001, 0.0],
							[ 0.0, 0.0, 0.02 ]])

		self.E_corner = matrix([[ 0.01, 0.0, 0.0 ],
							[ 0.0, 0.01, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

		self.E_inplace = matrix([[ 0.05, 0.0, 0.0 ],
							[ 0.0, 0.05, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

		self.E_overlap = matrix([[ 0.1,  0.0, 0.0],
							[ 0.0,  0.01, 0.0],
							[0.0, 0.0,  0.02]])

		self.E_featureless = matrix([[ 1.0,  0.0, 0.0],
							[ 0.0,  0.01, 0.0],
							[0.0, 0.0,  0.02]])
		#elf.E_featureless = matrix([[ 0.1,  0.0, 0.0],
		#					[ 0.0,  0.0001, 0.0],
		#					[0.0, 0.0,  0.02]])

		#self.E_overlap = matrix([[ 0.2,  0.0, 0.0],
		#					[ 0.0,  0.02, 0.0],
		#					[0.0, 0.0,  0.1]])

		self.E_motion = matrix([[ 0.06, 0.0, 0.0],
							[0.0,  5.0, 0.0],
							[0.0, 0.0,  0.1 ]])
		
		self.E_sensor = matrix([[0.1,0.0,0.0],
							[0.0,0.05,0.0],
							[0.0,0.0,0.02]])
		
		
		self.pathPlotCount = 0
		self.overlapPlotCount = 0
		self.pathDrawCount = 0
		#self.numPaths = 2
		#self.paths = {0 : [], 1 : []}
		#self.hulls = {0 : [], 1 : []}
		#self.nodeSets = {0 : [], 1 : []}
		self.numPaths = 1
		self.paths = {0 : []}
		self.hulls = {0 : []}
		self.nodeSets = {0 : []}
		self.currPath = 0
		self.junctionPoint = 0
		self.pathParents = [[None,None,None,None]]

		self.trimCount = 0
		self.spliceCount = 0
		self.orderCount = 0		
		self.stateCount = 0
		
		self.colors = []
		#for i in range(1000):
		#	self.colors.append((random.random(),random.random(),random.random()))

		self.colors = ['b','r','g','k','y']

	def saveState(self):
		
		saveFile = ""
		
		saveFile += "self.initPose = " + repr(self.initPose) + "\n"
		saveFile += "self.numNodes = " + repr(self.numNodes) + "\n"

		saveFile += "self.cornerBins = " + repr(self.cornerBins) + "\n"

		#self.cornerHash = {}
		"TODO: save localNode corner candidates "
		
		saveFile += "self.edgePriorityHash = " + repr(self.edgePriorityHash) + "\n"
		#saveFile += "self.edgeHash = " + repr(self.edgeHash) + "\n"
		
		#self.edgePriorityHash = {}
		
		saveFile += "self.numPaths = " + repr(self.numPaths) + "\n"
		saveFile += "self.nodeSets = " + repr(self.nodeSets) + "\n"
		saveFile += "self.junctionPoint = " + repr(self.junctionPoint) + "\n"

		saveFile += "self.currPath = " + repr(self.currPath) + "\n"

		saveFile += "self.pathParents = " + repr(self.pathParents) + "\n"
		

		f = open("stateSave_%04u.txt" % (self.numNodes-1), 'w')
		f.write(saveFile)
		f.close()		
		
		
		for nodeID in range(self.numNodes):
			
			f = open("cornerCandidates_%04u.txt" % (nodeID), 'w')
			f.write(repr(self.nodeHash[nodeID].cornerCandidates))
			f.close()
		
		
		
		
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/stateSave_%04u.txt" % (numNodes-1)
		f = open(dirName + "/stateSave_%04u.txt" % (numNodes-1), 'r')		
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		exec(saveStr)
		
		print self.numNodes
		print self.edgePriorityHash

		
	def splicePathIDs(self, pathIDs):
		
		if len(pathIDs) == 0:
			return []
		
		if len(pathIDs) == 1:
			return [self.trimmedPaths[pathIDs[0]]]

		" Assumption:  the pathIDs are connected paths with parent-child relations "


		" find the root "
		" pick any node, and go up the tree until we hit the root "
		currPathID = pathIDs[0]
		
		while True:
			parent = self.pathParents[currPathID]
			thisParentID = parent[0]
			
			if thisParentID == None:
				break
			
			if pathIDs.count(thisParentID) == 0:
				break
			
			currPathID = thisParentID

		" currPathID is the root path "
		rootPathID = currPathID
		
		" for each path, attempt to join with its parent path "
		joins = []
		for pathID in pathIDs:

			if self.pathParents[pathID] == None:
				continue
						
			parentPathID = self.pathParents[pathID][0]
			
			" parent does not concern us "
			if pathIDs.count(parentPathID) < 1:
				continue
			
			junctionNodeID = self.pathParents[pathID][1]
			localJunctionPoint = self.pathParents[pathID][2]
			poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
			junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)

			path1 = self.trimmedPaths[pathID]
			path2 = self.trimmedPaths[parentPathID]

			
			#print "junctionPoint =", junctionPoint
			#print "path1:", self.paths[0]
			#print "path2:", self.paths[1]
	
			minDist1 = 1e100
			minI1 = 0		
			for i in range(len(path1)):
				pnt = path1[i]
				dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
			
				if dist < minDist1:
					minDist1 = dist
					minI1 = i
	
			minDist2 = 1e100
			minI2 = 0		
			for i in range(len(path2)):
				pnt = path2[i]
				dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
			
				if dist < minDist2:
					minDist2 = dist
					minI2 = i
			
			joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])

		print "joins:", joins

		" create a tree with the node IDs and then stitch them together with joins "
		pathGraph = graph.graph()
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for k in range(len(path)):
				pathGraph.add_node((pathID, k), path[k])
			for k in range(len(path)-1):
				pathGraph.add_edge((pathID, k), (pathID, k+1))

		" join with the junction in between the join points "
		for k in range(len(joins)):
			join = joins[k]
			pathGraph.add_node(k, join[2])
			pathGraph.add_edge(join[0], k)
			pathGraph.add_edge(k, join[1])
			
		" if both terminal paths are child paths, 1 resultant spliced path "
		" if one terminal path is the root, then 2 resultant spliced paths "
		rightPathID = pathIDs[0]
		leftPathID = pathIDs[-1]

		newPaths = []		
		if rightPathID == rootPathID or leftPathID == rootPathID:
			
			print "one of paths is root make 2 resultant spliced paths"
			
			if rightPathID != rootPathID:
				childPathID = rightPathID
			else:
				childPathID = leftPathID
			
			" find path 1 from rootPathID to childPathID "
			" find path 2 from rootPathID to childPathID "
			rootPath = self.trimmedPaths[rootPathID]
			childPath = self.trimmedPaths[childPathID]


			for join in joins:
				if join[0][0] == childPathID:
					childJunctionPoint = join[2]
					childI = join[0][1]

			if childI < abs(childI-len(childPath)-1):
				childTermI = len(childPath)-1
			else:
				childTermI = 0

			"find path between: (childPathID, childTermI) to (rootPathID, 0) and (rootPathID, len(rootPath)-1)"
			shortestPathSpanTree, shortestDist = pathGraph.shortest_path((childPathID, childTermI))
			startNode1 = shortestPathSpanTree[(rootPathID, 0)]
			startNode2 = shortestPathSpanTree[(rootPathID, len(rootPath)-1)]
			
			print "len(shortestPathSpanTree) =", len(shortestPathSpanTree)
			#print "shorestPathSpanTree:", shortestPathSpanTree

			currNode = startNode1
			splicedPath1 = []
			while currNode != (childPathID, childTermI):
				splicedPath1.append(pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath1.append(pathGraph.get_node_attributes(currNode))

			currNode = startNode2
			splicedPath2 = []
			while currNode != (childPathID, childTermI):
				splicedPath2.append(pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath2.append(pathGraph.get_node_attributes(currNode))

			#print "splicedPath1:", splicedPath1
			#print "splicedPath2:", splicedPath2
			
			newPaths.append(splicedPath1)
			newPaths.append(splicedPath2)
			
		else:
			"find entire path from rightPathID to leftPathID "
			" start from point farthest from child's junction point "

			print "both paths are childs,  make 1 resultant spliced paths"

			rightPath = self.trimmedPaths[rightPathID]
			leftPath = self.trimmedPaths[leftPathID]
		
			" find the junction of this child's to its parent "	
			#joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])
			
			for join in joins:
				if join[0][0] == rightPathID:
					rightJunctionPoint = join[2]
					rightI = join[0][1]

				if join[0][0] == leftPathID:
					leftJunctionPoint = join[2]
					leftI = join[0][1]
					
			
			if leftI < abs(leftI-len(leftPath)-1):
				leftTermI = len(leftPath)-1
			else:
				leftTermI = 0

			if rightI < abs(rightI-len(rightPath)-1):
				rightTermI = len(rightPath)-1
			else:
				rightTermI = 0
				
			"find path between: (rightPathID, rightTermI) to (leftPathID, rightTermI) "
			shortestPathSpanTree, shortestDist = pathGraph.shortest_path((rightPathID, rightTermI))
			currNode = shortestPathSpanTree[(leftPathID, leftTermI)]

			splicedPath = []
			while currNode != (rightPathID, rightTermI):
				splicedPath.append(pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath.append(pathGraph.get_node_attributes(currNode))

			print "len(shortestPathSpanTree) =", len(shortestPathSpanTree)
			#print "splicedPath:", splicedPath
			#print "shorestPathSpanTree:", shortestPathSpanTree

			newPaths.append(splicedPath)


	
			pylab.clf()
	
			print "spliced path has", len(splicedPath), "points"
			xP = []
			yP = []
			for p in splicedPath:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP)
	
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
	
			pylab.title("Spliced Paths 1, numNodes = %d" % self.numNodes)
			pylab.savefig("splicedPath_%04u.png" % self.spliceCount)
			self.spliceCount += 1



		return newPaths
			


	def getTopology(self, nodes = 0):
		global topCount
		
		def convertAlphaUniform(a_vert, max_spacing = 0.04):
			
			" make the vertices uniformly distributed "
			
			new_vert = []
		
			for i in range(len(a_vert)):
				p0 = a_vert[i]
				p1 = a_vert[(i+1) % len(a_vert)]
				dist = math.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
	
				vec = [p1[0]-p0[0], p1[1]-p0[1]]
				vec[0] /= dist
				vec[1] /= dist
				
				new_vert.append(copy(p0))
				
				if dist > max_spacing:
					" cut into pieces max_spacing length or less "
					numCount = int(math.floor(dist / max_spacing))
					
					for j in range(1, numCount+1):
						newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
						new_vert.append(newP)
			
			return new_vert			

		medialPointSoup = []

		if nodes == 0:
			nodes = range(self.numNodes)
			
		if nodes == []:
			return [], []


		#pylab.clf()	

		for nodeID in nodes:
			estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()		
	
			if self.nodeHash[nodeID].isBowtie:			
				hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
			else:
				hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
	
			" set the origin of pose 1 "
			poseOrigin = Pose(estPose1)
	
			#xP = []
			#yP = []	
			for p in hull1:
				p1 = poseOrigin.convertLocalToGlobal(p)
				medialPointSoup.append(p1)
				#xP.append(p1[0])
				#yP.append(p1[1])

			#pylab.scatter(xP,yP)
		"""
		#pylab.xlim(-4.4)
		#pylab.ylim(-3,3)
		pylab.xlim(-4,4)
		pylab.ylim(-4,4)
		pylab.savefig("plotMedialSoup_%04u.png" % (topCount))
		topCount += 1
		"""
		


		nodeID = self.numNodes - 1
			
		radius = 0.2

		numPoints = len(medialPointSoup)
		inputStr = str(numPoints) + " "

		" alpha shape circle radius "
		inputStr += str(radius) + " "

		#for p in medialPointSoup:
		#	p2 = copy(p)
		#	" add a little bit of noise to avoid degenerate conditions in CGAL "
		#	inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		#
		#inputStr += "\n"


		isDone = False
		
		while not isDone:

			inputStr = str(numPoints) + " "
	
			" alpha shape circle radius "
			inputStr += str(radius) + " "
			
			for p in medialPointSoup:
				p2 = copy(p)
				" add a little bit of noise to avoid degenerate conditions in CGAL "
				p2[0] += random.gauss(0.0,0.000001)
				p2[1] += random.gauss(0.0,0.000001)
	
				inputStr += str(p2[0]) + " " + str(p2[1]) + " "
			
			inputStr += "\n"
			
			try:			
				" start the subprocess "
				if sys.platform == "win32":
					subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
				else:
					subProc = Popen(["./alpha2"], stdin=PIPE, stdout=PIPE)
					
				
				" send input and receive output "
				sout, serr = subProc.communicate(inputStr)
		
				" convert string output to typed data "
				sArr = sout.split(" ")
				
				numVert = int(sArr[0])
				
				sArr = sArr[1:]
				
				
				vertices = []
				for i in range(len(sArr)/2):
					vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
				isDone = True
			except:
				print "hull has holes!  retrying..."
				#print sArr	

		
		" cut out the repeat vertex "
		vertices = vertices[:-1]
		
		vertices = convertAlphaUniform(vertices)

		vertices.append(vertices[0])
		
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in vertices:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]
	
	
		PIXELSIZE = 0.05
		mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
		pixelSize = PIXELSIZE
		numPixel = int(2.0*mapSize / pixelSize + 1.0)
		divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)
		
		def realToGrid(point):
			indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
			indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
			return indexX, indexY
	
		def gridToReal(indices):
			i = indices[0]
			j = indices[1]
			point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
			return point
	
		gridHull = []
		for i in range(len(vertices)):
			p = vertices[i]
			gridHull.append(realToGrid(p))
	
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in gridHull:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]

		resultImg = Image.new('L', (numPixel,numPixel))
		resultImg = computeMedialAxis(nodeID, numPixel,numPixel, resultImg, len(gridHull[:-2]), gridHull[:-2])

		imgA = resultImg.load()
	
		points = []
		for i in range(1, numPixel-1):
			for j in range(1, numPixel-1):
				if imgA[i,j] == 0:
					points.append((i,j))
	
		medialGraph = graph.graph()
		for p in points:
			medialGraph.add_node(p, [])
	
		builtGraph = {}
		for i in range(2, numPixel-2):
			for j in range(2, numPixel-2):
				if imgA[i,j] == 0:
					builtGraph[(i,j)] = []
					for k in range(i-1, i+2):
						for l in range(j-1, j+2):
							if imgA[k,l] == 0:
								builtGraph[(i,j)].append((k,l))
								medialGraph.add_edge((i,j), (k,l))
								
		mst = medialGraph.minimal_spanning_tree()
		
		uni_mst = {}
		isVisited = {}
		nodeSum = {}
		for k, v in mst.items():
			uni_mst[k] = []
			isVisited[k] = 0
			nodeSum[k] = 0
		
		for k, v in mst.items():
			if v != None:
				uni_mst[k].append(v)
				uni_mst[v].append(k)
		
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
	
	
		" find the longest path from each leaf"
		maxPairDist = 0
		maxPair = None
		maxPath = []
		for leaf in leaves:
	
			isVisited = {}
			nodeSum = {}
			nodePath = {}
			for k, v in uni_mst.items():
				isVisited[k] = 0
				nodeSum[k] = 0
				nodePath[k] = []
	
			getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
	
			maxDist = 0
			maxNode = None
			for k, v in nodeSum.items():
				if v > maxDist:
					maxNode = k
					maxDist = v
	
			if maxDist > maxPairDist:
				maxPairDist = maxDist
				maxPair = [leaf, maxNode]
				maxPath = nodePath[maxNode]

		" convert path points to real "
		realPath = []
		for p in maxPath:
			realPath.append(gridToReal(p))
		
		return realPath, vertices

	def checkSupport(self, nodeID1, nodeID2, offset, supportLine):

		node1 = self.nodeHash[nodeID1]
		node2 = self.nodeHash[nodeID2]

		#hull1, medial1 = computeHullAxis(nodeID1, node1, tailCutOff = True)
		#hull2, medial2 = computeHullAxis(nodeID2, node2, tailCutOff = True)
		hull2, medial2 = computeHullAxis(nodeID2, node2, tailCutOff = False)
		
		estPose1 = node1.getGlobalGPACPose()		
			
		#minMatchDist2 = 0.5
		minMatchDist2 = 2.0
			
		" set the initial guess "
		poseOrigin = Pose(estPose1)
		
		localSupport = []
		for pnt in supportLine:
			localSupport.append(poseOrigin.convertGlobalToLocal(pnt))
	
		supportSpline = SplineFit(localSupport, smooth=0.1)		
		supportPoints = gen_icp.addGPACVectorCovariance(supportSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		points2 = gen_icp.addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)

	
		" transform pose 2 by initial offset guess "	
		" transform the new pose "
		points2_offset = []
		for p in points2:
			result = gen_icp.dispPoint(p, offset)		
			points2_offset.append(result)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])			
		
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
				
		support_pairs = []
		for i in range(len(poly2)):
			
			p_2 = poly2[i]
	
			" for every transformed point of A, find it's closest neighbor in B "
			p_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

			if gen_icp.isInCircle(p_1, radius2, center2):

				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = p_1[2]

					" we store the untransformed point, but the transformed covariance of the A point "
					support_pairs.append([points2[i],p_1,C2,C1])

		vals = []
		sum1 = 0.0
		for pair in support_pairs:
	
			a = pair[0]
			b = pair[1]
			Ca = pair[2]
			Cb = pair[3]
	
			ax = a[0]
			ay = a[1]		
			bx = b[0]
			by = b[1]
	
			c11 = Ca[0][0]
			c12 = Ca[0][1]
			c21 = Ca[1][0]
			c22 = Ca[1][1]
					
			b11 = Cb[0][0]
			b12 = Cb[0][1]
			b21 = Cb[1][0]
			b22 = Cb[1][1]	
		
			val = gen_icp.computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
			
			vals.append(val)
			sum1 += val
			
		return sum1


	def addPriorityEdge(self, edge, priority):

		try:
			self.edgePriorityHash[(edge[0],edge[1])]
		except:
			self.edgePriorityHash[(edge[0],edge[1])] = []
		
		" corner constraints with priority of 3 "
		self.edgePriorityHash[(edge[0],edge[1])].append([edge[2], edge[3], priority])

	def getEdges(self, nodeID1, nodeID2):

		resultEdges = []
		
		for k, v in self.edgePriorityHash.items():

			if k[0] == nodeID1 and k[1] == nodeID2 or k[0] == nodeID2 and k[1] == nodeID1:
				
				for edge in v:
					resultEdges.append(edge)
				
		return resultEdges

	def getPriorityEdges(self, priorityLevel = -1):

		priorityEdges = {}
		
		if priorityLevel == -1:
			for k, v in self.edgePriorityHash.items():
				
				if len(v) > 0:
			
					id1 = k[0]
					id2 = k[1]
					
					" take the highest priority constraints only "
					maxPriority = -1
					for const in v:
						thisPriority = const[2]
						if thisPriority > maxPriority:
							maxPriority = thisPriority
			
					if maxPriority == -1:
						raise
					
					priorityEdges[k] = []
					for const in v:
						if const[2] == maxPriority:
							priorityEdges[k].append(const)

		else:
			
			for k, v in self.edgePriorityHash.items():
				
				if len(v) > 0:
			
					id1 = k[0]
					id2 = k[1]
					
					" take the highest priority constraints only "
					priorityEdges[k] = []
					for const in v:
						if priorityLevel == const[2]:
							
							priorityEdges[k].append(const)


		return priorityEdges

	def deleteAllPriority(self, priorityLevel):
		
		for k, v in self.edgePriorityHash.items():
			newV = []
			
			for const in v:
				if const[2] != priorityLevel:
					newV.append(const)
			self.edgePriorityHash[k] = newV


	def deleteAllEdges(self):

		self.edgePriorityHash = {}

	def getPathTerms(self):

		terms = []

		numPaths = len(self.paths)
		for k in range(numPaths):

			path = self.paths[k]
			
			pathSpline = SplineFit(path, smooth=0.1)
			vecPoints = pathSpline.getUniformSamples()
			
			angleSum = 0.0
			for p in vecPoints[-11:]:
				angleSum += p[2]
				
			angleAvg = angleSum / 10.
			
			terms.append([vecPoints[-1][0],vecPoints[-1][1], angleAvg])

		return terms

	def trimPaths(self, paths):

		trimmedPaths = []
		
		if len(paths) <= 1:
			for k, v in paths.iteritems():
				trimmedPaths.append(v)
				
			return trimmedPaths

		" path0 is root, and path1 is a child of path0 "
		#parents = [None, [0, self.junctionNodeID, self.localJunctionPoint]]
		parents = self.pathParents

		" information on parentage is required to determine which path belongs to which "

		" global junction point for parent 0 and child 1 "

		print "pathParents:", parents
		
		for pathID in range(len(paths)):
			
			if parents[pathID][0] == None:
				trimmedPaths.append(deepcopy(paths[pathID]))

			else:

				path1 = paths[parents[pathID][0]]
				path2 = paths[pathID]
				
				parentPathID = parents[pathID][0]
				childPathID = pathID

				#secP1, secP2 = self.getOverlapDeparture(path1, path2)				
				secP1, secP2 = self.getOverlapDeparture(parentPathID, childPathID, paths)				

				minI_1 = 0		
				minI_2 = 0
				minDist_1 = 1e100
				minDist_2 = 1e100		
				for i in range(len(path2)):
					pnt = path2[i]
					dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
					dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
				
					if dist1 < minDist_1:
						minDist_1 = dist1
						minI_1 = i

					if dist2 < minDist_2:
						minDist_2 = dist2
						minI_2 = i

				term0 = path2[0]
				termN = path2[-1]

				" smallest distance is the terminal point "
				dist0_1 = sqrt((term0[0]-secP1[0])**2 + (term0[1]-secP1[1])**2)
				distN_1 = sqrt((termN[0]-secP1[0])**2 + (termN[1]-secP1[1])**2)
				dist0_2 = sqrt((term0[0]-secP2[0])**2 + (term0[1]-secP2[1])**2)
				distN_2 = sqrt((termN[0]-secP2[0])**2 + (termN[1]-secP2[1])**2)
				print "terminal distances:", dist0_1, distN_1, dist0_2, distN_2
				
				distList = [dist0_1, distN_1, dist0_2, distN_2]
				
				minI = -1
				minDist = 1e100
				for i in range(len(distList)):
					if distList[i] < minDist:
						minDist = distList[i]
						minI = i

				junctionNodeID = parents[pathID][1]
				localJunctionPoint = parents[pathID][2]

				poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
				junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)
				
				
				if minDist_1 < minDist_2:
					junctionPoint_K = secP1
					juncI_K = minI_1
				else:
					junctionPoint_K = secP2
					juncI_K = minI_2
						
				#junctionPoint = junctionPoint_K
				#juncI = juncI_K
				
				
				minDist2 = 1e100
				juncI = 0		
				for i in range(len(path2)):
					pnt = path2[i]
					dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
				
					if dist < minDist2:
						minDist2 = dist
						juncI = i
 				
 				
 				
				print "len(path1):", len(path1)
				print "len(path2):", len(path2)
				print "juncI:", juncI
				print "minDist:", minDist_1, minDist_2
				
				" now we have closest point to departure point. "
				" Which side is the departing side? "	
				#sum1 = self.getOverlapSum(path1, path2[0:minI2+1] + [junctionPoint])
				#sum2 = self.getOverlapSum(path1, [junctionPoint] + path2[minI2:])
				

						
				if minI == 0:
					"secP1 is terminal 0"
					index = juncI-10
					if index < 1:
						index = 1
					#newPath2 = path2[:index+1] + [junctionPoint]
					newPath2 = path2[:index+1]
					newPath2.reverse()
				
				elif minI == 1:
					"secP1 is terminal N"
					index = juncI+10
					if index >= len(path2):
						index = len(path2)-1
					#newPath2 = [junctionPoint] + path2[index:]
					newPath2 = path2[index:]
					
				elif minI == 2:
					"secP2 is terminal 0"
					index = juncI-10
					if index < 1:
						index = 1
					#newPath2 = path2[:index+1] + [junctionPoint]
					newPath2 = path2[:index+1]
					newPath2.reverse()
					
				elif minI == 3:
					"secP2 is terminal N"
					index = juncI+10
					if index >= len(path2):
						index = len(path2)-1
					#newPath2 = [junctionPoint] + path2[index:]
					newPath2 = path2[index:]
				
				else:
					print "no terminal found"
					raise
								
				#if minI_1 > minI_2:
				#	path2[minI_2:minI_1+1]
				#else:
				#	path2[minI_1:minI_2+1]
				
				
				#if startI == 0:
				#	newPath2 = path2[startI:endI+1] + [junctionPoint]
				#	newPath2.reverse()
				#else:
				#	newPath2 = [junctionPoint] + path2[startI:endI+1]
					

				#print "Sums:", sum1, sum2
				
				#if sum1 > sum2:
				#	index = minI2+1-10
				#	if index < 1:
				#		index = 1
				#	newPath2 = path2[0:index] + [junctionPoint]
				#	newPath2.reverse()
				#else:
				#	index = minI2+10
				#	if index >= len(path2):
				#		index = len(path2)-1
				#	newPath2 = [junctionPoint] + path2[index:]

				" convert path so that the points are uniformly distributed "
				
				
				max_spacing = 0.08
				#max_spacing = 0.04
				#max_spacing = 0.1
				newPath3 = []

				" make sure path is greater than 5 points "
				while len(newPath3) <= 5:
	
					max_spacing /= 2
					
					newPath3 = [copy(newPath2[0])]
										
					for i in range(len(newPath2)-1):
						p0 = newPath2[i]
						p1 = newPath2[(i+1)]
						dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
			
						vec = [p1[0]-p0[0], p1[1]-p0[1]]
						vec[0] /= dist
						vec[1] /= dist
						
						
						if dist > max_spacing:
							" cut into pieces max_spacing length or less "
							numCount = int(floor(dist / max_spacing))
							
							for j in range(1, numCount+1):
								newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
								newPath3.append(newP)
	
						newPath3.append(copy(p1))			
				
				trimmedPaths.append(deepcopy(newPath3))

		"""
		self.drawTrimmedPaths(self.trimmedPaths)
		"""
		
		
		return trimmedPaths


	def getOrderedOverlappingPaths(self, nodeID):

		overlappedPaths = []
		
		"FIXME:  add departure detection and junction continuity checking "
		"        to reject paths that we are only incident to "
	
		print "getOrderedOverlappingPaths(", nodeID, ")"
	
		overlapSums = {}
	
		for pathID in range(len(self.trimmedPaths)):
			sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
			
			print "overlap sum", pathID, "=", sum1
			
			if sum1 <= 1e10:
				overlappedPaths.append(pathID)
				overlapSums[pathID] = sum1

		print "overlapSums:", overlapSums

		orderedPathIDs = self.getPathOrdering(nodeID, overlappedPaths)

		print "orderedPathIDs:", orderedPathIDs

		departures = []
		interiors = []
		depPoints = []

		for pathID in orderedPathIDs:
			departurePoint1, depAngle2, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID)
			departures.append([isExist1,isExist2])
			interiors.append([isInterior1, isInterior2])
			depPoints.append([departurePoint1, departurePoint2])
		
		print "departures:", departures
		print "interiors:", interiors
			
		" determine which paths are leaves "
		isAParent = [False for k in range(self.numPaths)]
		for k in orderedPathIDs:
			currPath = self.pathParents[k]
			currParent = currPath[0]
			if currParent != None:
				isAParent[currParent] = True

		print "isAParent:", isAParent
		
		toBeRemoved = []
		seqSize = len(orderedPathIDs)
		
		" TODO: is this necessary still? "
		" check for proper parent-child relationship and "
		" that there is interior departure between paths "

		
		
		for k in range(seqSize-1):
			pathID1 = orderedPathIDs[k]
			pathID2 = orderedPathIDs[k+1]

			" there should be a interior departure between these paths on the parent path "
			if self.pathParents[pathID1][0] == pathID2:
				pass
				#parentPath = pathID2
				
			elif self.pathParents[pathID2][0] == pathID1:
				pass
				#parentPath = pathID1

					
			else:
				print "ordered path IDs could not determine relationship", pathID1, pathID2, orderedPathIDs
				raise
		
		"""		
		for k in range(seqSize-1):
			pathID1 = orderedPathIDs[k]
			pathID2 = orderedPathIDs[k+1]

			" there should be a interior departure between these paths on the parent path "
			if self.pathParents[pathID1][0] == pathID2:
				parentPath = pathID2

				" if the departure doesn't show, remove the least overlapped path section "
				if not interiors[k+1][0]:
					if overlapSums[pathID1] > overlapSums[pathID2]:
						print "remove pathID1", pathID1
						#toBeRemoved.append(k)
					else:
						print "remove pathID2", pathID2
						#toBeRemoved.append(k+1)
					
					if k != 0 and k != seqSize-1:
						print "paths in the middle violate a departure constraint!"
						print pathID1, pathID2, orderedPathIDs
						print departures, interiors
						print overlapSums
						print toBeRemoved
						raise
				
				
			elif self.pathParents[pathID2][0] == pathID1:
				parentPath = pathID1
				if not interiors[k][1]:
					if overlapSums[pathID1] > overlapSums[pathID2]:
						print "remove pathID1", pathID1
						#toBeRemoved.append(k)
					else:
						print "remove pathID2", pathID2
						#toBeRemoved.append(k+1)

					if k != 0 and k != seqSize-1:
						print "paths in the middle violate a departure constraint!"
						print pathID1, pathID2, orderedPathIDs
						print departures, interiors
						print overlapSums
						print toBeRemoved
						raise
					
			else:
				print "ordered path IDs could not determine relationship", pathID1, pathID2, orderedPathIDs
				raise
		"""
		
		print "toBeRemoved:", toBeRemoved

		finalPathIDs = []
		for k in range(seqSize):
			
			if toBeRemoved.count(k) == 0:
				finalPathIDs.append(orderedPathIDs[k])
			
		print "finalPathIDs:", finalPathIDs	
		return finalPathIDs

	def getPathOrdering(self, nodeID, pathIDs):

		plotIter = False

		"""
		for pathID in pathIDs:
			
			path = self.paths[pathID]
			hull = self.hulls[pathID]
			overlapCost1 = self.getOverlapCondition(path, nodeID1)
			overlapCost2 = self.getOverlapCondition(path, nodeID2)


			departurePoint1, isInterior1, isExist1 = self.getDeparturePoint(path, nodeID1, hull)
			departurePoint2, isInterior2, isExist2 = self.getDeparturePoint(path, nodeID2, hull)
		"""

		node2 = self.nodeHash[nodeID]

		hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)

		estPose2 = node2.getGlobalGPACPose()		
			
		#minMatchDist2 = 0.2
		#minMatchDist2 = 0.5
		#minMatchDist2 = 2.0
			
		" set the initial guess "
		poseOrigin = Pose(estPose2)
						
		localizedPaths = []
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			
			localPath = []
			for pnt in path:
				localPath.append(poseOrigin.convertGlobalToLocal(pnt))
			
			localizedPaths.append(localPath)
			
		" find minimum distance and its associated path "					
		#supportSpline = SplineFit(localSupport, smooth=0.1)		
		#vecPoints1 = supportSpline.getUniformSamples()
		#supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			
		
		pathSelect = []
		#minMatchDist2 = 0.5
		minMatchDist2 = 1.0

		for i in range(len(vecPoints2)):

			minDist = 1e100
			minPathID = -1
			
			p_2 = vecPoints2[i]

			for j in range(len(localizedPaths)):


				path = localizedPaths[j]
				medialSpline1 = SplineFit(path, smooth=0.1)
				vecPoints1 = medialSpline1.getUniformSamples()
				supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
	
				" for every transformed point of A, find it's closest neighbor in B "
				try:
					p_1, minI, minDist_I = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
		
					if minDist_I <= minMatchDist2:
						C2 = points2[i][2]
						C1 = supportPoints[minI][2]

						ax = supportPoints[minI][0]
						ay = supportPoints[minI][1]		
						bx = points2[i][0]
						by = points2[i][1]
				
						c11 = C1[0][0]
						c12 = C1[0][1]
						c21 = C1[1][0]
						c22 = C1[1][1]
								
						b11 = C2[0][0]
						b12 = C2[0][1]
						b21 = C2[1][0]
						b22 = C2[1][1]	
						
						dist = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])

						if dist < minDist:
							minDist = dist
							minPathID = pathIDs[j]
						#" we store the untransformed point, but the transformed covariance of the A point "
						#support_pairs.append([points2[i],supportPoints[minI],C2,C1])

				
				except:
					#raise
					pass

				#print "minPathID:", minPathID

			if minPathID != -1:
				pathSelect.append(minPathID)
		
		print "pathSelect:", pathSelect
		
		" find the average position of each path ID"
		avgIDPosition = {}
		avgIDCount = {}
		for pathID in pathIDs:
			avgIDPosition[pathID] = 0
			avgIDCount[pathID] = 0
			
		for i in range(len(pathSelect)):
			avgIDPosition[pathSelect[i]] += i
			avgIDCount[pathSelect[i]] += 1


		for pathID in pathIDs:
			if avgIDCount[pathID] != 0:
				avgIDPosition[pathID] /= avgIDCount[pathID]
			else:
				del avgIDPosition[pathID]
			
		" now sort out the order "
		orderedPaths = sorted(avgIDPosition, key=avgIDPosition.__getitem__)		
		
		#print "orderedPaths:", orderedPaths
		
		
		if plotIter:
			pylab.clf()
			
			xP = []
			yP = []
			
			for p in poly2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP)
			
			
			for i in range(len(localizedPaths)):
				xP = []
				yP = []
				path = localizedPaths[i]
				for p in path:
					xP.append(p[0])
					yP.append(p[1])
				
				pylab.plot(xP,yP)
	
			pylab.title("%d: %s %s" % (nodeID, repr(pathIDs), repr(orderedPaths)))
			pylab.savefig("orderedPath_%04u.png" % self.orderCount)
			
			self.orderCount += 1
		
		return orderedPaths

	def getDeparturePoint(self, currPath, nodeID):
		
		isExist1 = False
		isInterior1 = False
		departurePoint1 = 0
		angle1 = 0.0
		isExist2 = False
		isInterior2 = False
		departurePoint2 = 0
		angle2 = 0.0
		
		if len(currPath) == 0:
			return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
		
		node2 = self.nodeHash[nodeID]

		hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		
		estPose2 = node2.getGlobalGPACPose()		
		
		"Assumption:  one section of the medial axis is closely aligned with the path "
		poseOrigin = Pose(estPose2)
		
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		points2 = medialSpline2.getUniformSamples()

		pathSpline = SplineFit(currPath, smooth=0.1)
		pathPoints = pathSpline.getUniformSamples()

		
		points2_offset = []
		for p in points2:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points2_offset.append(result)


		" tip angles "
		angSum1 = 0.0
		angSum2 = 0.0
		angs1 = []
		angs2 = []
		phi1 = normalizeAngle(points2_offset[0][2])
		phi2 = normalizeAngle(points2_offset[-1][2])
		for i in range(10):
			ang1 = normalizeAngle(points2_offset[i][2]-phi1)
			ang2 = normalizeAngle(points2_offset[-i-1][2]-phi2)
			angSum1 += ang1
			angSum2 += ang2
			
			angs1.append(ang1+phi1)
			angs2.append(ang2+phi2)

		angle1 = angSum1 / 10.0 + phi1
		angle2 = angSum2 / 10.0 + phi2

		" invert one angle so opposite tips have opposite angles "
		angle1 = normalizeAngle(angle1 + pi)

		print "ang1:", angle1, angs1
		print "ang2:", angle2, angs2
		print "diff:", diffAngle(angle1, angle2)

		distances = []
		indices = []
		for i in range(0,len(points2_offset)):
			p_2 = points2_offset[i]
			#p_1, minDist = gen_icp.findClosestPointInB(pathPoints, p_2, [0.0,0.0,0.0])
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints, p_2)
			distances.append(minDist)
			indices.append(i_1)
		
		
		" Compute the front and back departure points by finding the inflection point on the distance curve "
		" these indices become frontDepI and backDepI respectively "
		maxFront = distances[0]
		maxBack = distances[-1]

		currI = 1
		try:
			while distances[currI+3] < maxFront:
				maxFront = distances[currI]
				currI += 1
		except:
			pass
		
		frontDepI = currI
		frontPoint = [frontDepI, distances[frontDepI]]

		" FIXME:  index out of bounds case "
		currI = 2
		try:
			while distances[-currI-3] < maxBack:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass

		backDepI = len(distances) - currI
		backPoint = [backDepI, distances[backDepI]]

		"reset to the maximum distances "
		maxFront = distances[0]
		maxBack = distances[-1]
		
		" perform a line fit on the distance curves from the departure point to the tip"	
		#(ar1,br1)= scipy.polyfit(range(0,frontDepI+1),distances[0:frontDepI+1],1)
		#(ar2,br2)= scipy.polyfit(range(backDepI,len(distances)),distances[backDepI:],1)

		#t1 = range(0,frontDepI+1)
		#xr1 = scipy.polyval([ar1,br1],t1)			
		
		#t2 = range(backDepI,len(distances))
		#xr2 = scipy.polyval([ar2,br2],t2)			
		
		" compute the average of distances between departure and termination "
		#frontAvg = 0.0
		#for i in range(0,frontDepI+1):
		#	frontAvg += distances[i]
		#frontAvg /= (frontDepI+1)

		#backAvg = 0.0
		#for i in range(backDepI,len(distances)):
		#	backAvg += distances[i]
		#backAvg /= len(distances)-backDepI

		" for all the points matched from the local curve to the path curve "
		" count the number of times they are matched "
		foo = indices[0:frontDepI+1]
		d1 = {}
		for i in set(foo):
			d1[i] = foo.count(i)

		foo = indices[backDepI:]
		d2 = {}
		for i in set(foo):
			d2[i] = foo.count(i)

		" find the point that has the most matches "
		max1 = max(d1, key=d1.get)
		max2 = max(d2, key=d2.get)



		arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
		arr2 = numpy.array(deepcopy(indices[backDepI:]))
		matchVar1 = arr1.var()
		matchVar2 = arr2.var()

		" compute the selected point index by average instead of by maximum"
		#frontI = 0.0
		#for i in range(0,frontDepI+1):
		#	frontI += indices[i]
		#frontI /= (frontDepI+1)

		#backI = 0.0
		#for i in range(backDepI,len(distances)):
		#	backI += indices[i]
		#backI /= len(distances)-backDepI

		tipMatch1 = pathPoints[indices[0]]

		" discrepancy distance between tip closest point and average departure point "
		dist1 = sqrt((tipMatch1[0]-pathPoints[max1][0])**2 + (tipMatch1[1] - pathPoints[max1][1])**2)

		tipMatch2 = pathPoints[indices[-1]]

		" discrepancy distance between tip closest point and average departure point "
		dist2 = sqrt((tipMatch2[0]-pathPoints[max2][0])**2 + (tipMatch2[1] - pathPoints[max2][1])**2)


		DEP_THRESH = 0.3

		if maxFront > DEP_THRESH:
			departurePoint1 = pathPoints[max1]
			isExist1 = True



			if max1 == 0 or max1 == len(pathPoints)-1:
				isInterior1 = False
			else:
				isInterior1 = True

			# Test purposes only
			#if nodeID == 4 or nodeID == 5:
			#	isInterior1 = True


		if maxBack > DEP_THRESH:
			departurePoint2 = pathPoints[max2]
			isExist2 = True


			if max2 == 0 or max2 == len(pathPoints)-1:
				isInterior2 = False
			else:
				isInterior2 = True



			
		"""
		if maxFront > 0.5 and maxFront > maxBack:
			departurePoint = pathPoints[max1]
			isExist = True

			if max1 == 0 or max1 == len(pathPoints)-1:
				isInterior = False
			else:
				isInterior = True

		if maxBack > 0.5 and maxBack > maxFront:
			departurePoint = pathPoints[max2]
			isExist = True

			if max2 == 0 or max2 == len(pathPoints)-1:
				isInterior = False
			else:
				isInterior = True
		"""	
			
		" sum of closest points on front and back "
		" select the one with minimal cost "
		
		if False:
			pylab.clf()
			xP = range(len(points2_offset))
			yP = distances
			pylab.plot(xP,yP, color ='b')
			
			if maxFront > 0.5:
				xP = [frontPoint[0]]
				yP = [frontPoint[1]]
				pylab.scatter(xP,yP, color='b')
	
			if maxBack > 0.5:
				xP = [backPoint[0]]
				yP = [backPoint[1]]
				pylab.scatter(xP,yP, color='b')
			
			pylab.xlim(0,200)
			pylab.ylim(0,2)
			#pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
			pylab.title("nodeID %d: %1.2f %1.2f %d %d %d" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2))
			pylab.savefig("distances_%04u.png" % self.pathPlotCount)

		if True:
			pylab.clf()
			xP = []
			yP = []
			for p in points2_offset:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')
	
			if True:	
				xP = [pathPoints[max1][0]]
				yP = [pathPoints[max1][1]]
				pylab.scatter(xP,yP, color='b')		
	
				xP = [tipMatch1[0]]
				yP = [tipMatch1[1]]
				pylab.scatter(xP,yP, color='r')		
	
	
			if True:
				xP = [pathPoints[max2][0]]
				yP = [pathPoints[max2][1]]
				pylab.scatter(xP,yP, color='b')		
	
				xP = [tipMatch2[0]]
				yP = [tipMatch2[1]]
				pylab.scatter(xP,yP, color='r')		
	
	
			xP = []
			yP = []
			for p in currPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			#pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
			pylab.title("%d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2))
			pylab.savefig("departure_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1
		
		
		#return departurePoint, isInterior, isExist
		return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2
		

	" returns the endpoints of a subpath of path 2 that does not overlap path 1 "

	" get the trimmed version of child and parent paths that are overlapping in some fashion "
	def getOverlapDeparture(self, parentPathID, childPathID, paths):

		"Assumption:  one section of the medial axis is closely aligned with the path "		
			
		plotIter = True
		print "getOverlapDeparture():"
		
		path1 = paths[parentPathID]
		path2 = paths[childPathID]
		
		
		isExist1 = False
		isInterior1 = False
		departurePoint1 = 0
		isExist2 = False
		isInterior2 = False
		departurePoint2 = 0


		branchNodeID = self.pathParents[childPathID][1]
		poseOrigin = Pose(self.nodeHash[branchNodeID].getGlobalGPACPose())
		localJunctionPoint = self.pathParents[childPathID][2]
		
		
		poseOrigin2 = Pose(self.nodeHash[branchNodeID].getEstPose())
		globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)
		
		
		" orienting the medial axis of the branch node correctly "
		branchDir = self.pathParents[childPathID][3]
		
		branchHull, branchMedial = computeHullAxis(branchNodeID, self.nodeHash[branchNodeID], tailCutOff = False)

		branchMedialSpline = SplineFit(branchMedial, smooth=0.1)
		branchPoints = branchMedialSpline.getUniformSamples()
		
		globalBranchPoints = []
		for p in branchPoints:
			result = poseOrigin.convertLocalToGlobal(p)
			globalBranchPoints.append(result)
		
		if not branchDir:
			" backward "
			globalBranchPoints.reverse()
			
		globalBranchSpline = SplineFit(globalBranchPoints, smooth = 0.1)

			

		" return exception if we receive an invalid path "		
		if len(path1) == 0:
			print "path1 has zero length"
			raise
			#return departurePoint1, isInterior1, isExist1, departurePoint2, isInterior2, isExist2


		" make sure the overlap of both paths are oriented the same way "
		path1Spline = SplineFit(path1, smooth=0.1)
		path2Spline = SplineFit(path2, smooth=0.1)

		path2Reverse = deepcopy(path2)
		path2Reverse.reverse()
		path2SplineReverse = SplineFit(path2Reverse, smooth=0.1)
		
		overlapMatch = []
		angleSum1 = 0.0
		angleSum2 = 0.0
		for i in range(0,len(path1)):
			p_1 = path1[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(path2, p_1)
			if minDist < 0.1:
				overlapMatch.append((i,j,minDist))

				pathU1 = path1Spline.findU(p_1)	
				pathU2 = path2Spline.findU(p_2)	
				pathU2_R = path2SplineReverse.findU(p_2)	

				pathVec1 = path1Spline.getUVector(pathU1)
				pathVec2 = path2Spline.getUVector(pathU2)
				pathVec2_R = path2SplineReverse.getUVector(pathU2_R)

				val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
				if val1 > 1.0:
					val1 = 1.0
				elif val1 < -1.0:
					val1 = -1.0
				ang1 = acos(val1)
				
				val2 = pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1]
				if val2 > 1.0:
					val2 = 1.0
				elif val2 < -1.0:
					val2 = -1.0
				ang2 = acos(val2)
				
				angleSum1 += ang1
				angleSum2 += ang2
	
		" select global path orientation based on which has the smallest angle between tangent vectors "
		print "angle sums for orienting child path with parent path:", angleSum1, angleSum2

		if angleSum1 > angleSum2:
			orientedPath2 = path2Reverse
		else:
			orientedPath2 = path2
			
		path2Spline = SplineFit(orientedPath2, smooth=0.1)



		" for each point on the child path, find its closest pair on the parent path "
		
		pathPoints1 = path1Spline.getUniformSamples()
		pathPoints2 = path2Spline.getUniformSamples()

		distances = []
		indices = []
		for i in range(0,len(pathPoints2)):
			p_2 = pathPoints2[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
			
			" keep the distance information "
			distances.append(minDist)
			" and the associated index of the point on the parent path "
			indices.append(i_1)

		" match distances of the tip points of child path "		
		maxFront = distances[0]
		maxBack = distances[-1]
		print "match distances of tip points:", maxFront, maxBack

		" walk back from tip point until we have a non-monotic increase in match distance "
		" this becomes our departure point "
		
		"TODO:  why is there a 3-offset in this comparison? "
		currI = 1
		try:
			while distances[currI+3] < maxFront:
				maxFront = distances[currI]
				currI += 1
		except:
			pass
		
		" departure index on child path "
		frontDepI = currI
		
		" departure index of child path and match distance "
		frontPoint = [frontDepI, distances[frontDepI]]

		" FIXME:  index out of bounds case "
		currI = 2
		try:
			while distances[-currI-3] < maxBack:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass

		" departure index on child path "
		backDepI = len(distances) - currI

		" departure index of child path and match distance "
		backPoint = [backDepI, distances[backDepI]]

		print "lengths of parent and child paths:", len(pathPoints1), len(pathPoints2)
		print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint


		"reset to the tip match distance "
		maxFront = distances[0]
		maxBack = distances[-1]
		
		
		" count the number of times a matched point on the parent path is used "
		foo = indices[0:frontDepI+1]
		d1 = {}
		for i in set(foo):
			d1[i] = foo.count(i)

		foo = indices[backDepI:]
		d2 = {}
		for i in set(foo):
			d2[i] = foo.count(i)

		" select the match point that is used the most on the parent path "
		max1 = max(d1, key=d1.get)
		max2 = max(d2, key=d2.get)


		" NOTE:  This guarantees detection of a departure.  We've assumed that there exists a"
		" departure between these two paths "
		
		" compute two candidate departure points "
		
		" NOTE:  departure point on parent path is not used here "
		if True:

			#departurePoint1 = pathPoints1[max1]
			isExist1 = True

			" if the highest match point is either of the tips of the parent path, then this is an external departure "
			if max1 == 0 or max1 == len(pathPoints1)-1:
				isInterior1 = False
			else:
				isInterior1 = True

		if True:
			
			#departurePoint2 = pathPoints1[max2]
			isExist2 = True

			" if the highest match point is either of the tips of the parent path, then this is an external departure "
			if max2 == 0 or max2 == len(pathPoints1)-1:
				isInterior2 = False
			else:
				isInterior2 = True

		print "isExist1 =", isExist1, "isInterior1 =", isInterior1
		print "isExist2 =", isExist2, "isInterior2 =", isInterior2

		" sum of closest points on front and back "
		" select the one with minimal cost "		
		
		" now we compare our candidates to the known direction and location of the branching point "
		angleSum1 = 0.0
		overlapSum1 = 0.0
		matchCount1 = 0
		
		" path section for our front departure hypothesis "
		pathSec1 = pathPoints2[:frontDepI+1]
		pathSec1.reverse()
		
		if len(pathSec1) > 5:
			pathSec1Spline = SplineFit(pathSec1, smooth = 0.1)		
			
			
			for i in range(0,len(pathSec1)):
				p_1 = pathSec1[i]
				p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
				if minDist < 0.1:
					overlapMatch.append((i,j,minDist))
					matchCount1 += 1
					
					overlapSum1 += minDist
	
					pathU1 = pathSec1Spline.findU(p_1)	
					pathUB = globalBranchSpline.findU(p_2)	
	
					pathVec1 = pathSec1Spline.getUVector(pathU1)
					pathVecB = globalBranchSpline.getUVector(pathUB)
	
					val1 = pathVec1[0]*pathVecB[0] + pathVec1[1]*pathVecB[1]
					if val1 > 1.0:
						val1 = 1.0
					elif val1 < -1.0:
						val1 = -1.0
					ang1 = acos(val1)
					
					angleSum1 += ang1

			if matchCount1 > 0:
				angleSum1 /= matchCount1
				overlapSum1 /= matchCount1
		

		angleSum2 = 0.0
		overlapSum2 = 0.0
		matchCount2 = 0
		
		" path section for our back departure hypothesis "
		pathSec2 = pathPoints2[backDepI:]

		if len(pathSec2) > 5:
			pathSec2Spline = SplineFit(pathSec2, smooth = 0.1)		
	
	
			for i in range(0,len(pathSec2)):
				p_1 = pathSec2[i]
				p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
				if minDist < 0.1:
					overlapMatch.append((i,j,minDist))
					matchCount2 += 1
					
					overlapSum2 += minDist
	
					pathU2 = pathSec2Spline.findU(p_1)	
					pathUB = globalBranchSpline.findU(p_2)	
	
					pathVec2 = pathSec2Spline.getUVector(pathU2)
					pathVecB = globalBranchSpline.getUVector(pathUB)
	
					val2 = pathVec2[0]*pathVecB[0] + pathVec2[1]*pathVecB[1]
					if val2 > 1.0:
						val2 = 1.0
					elif val2 < -1.0:
						val2 = -1.0
					ang2 = acos(val2)
					
					angleSum2 += ang2

			if matchCount2 > 0:
				angleSum2 /= matchCount2
				overlapSum2 /= matchCount2

		print "pathSec1 hypothesis angle and overlap sum and match count:", angleSum1, overlapSum1, matchCount1
		print "pathSec2 hypothesis angle and overlap sum and match count:", angleSum2, overlapSum2, matchCount2

		" distance of departure point from known junction point "
		juncDist1 = -1
		juncDist2 = -1
		if len(pathSec1) > 0:
			p0 = pathSec1[0]
			juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

		if len(pathSec2) > 0:
			p0 = pathSec2[0]
			juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

		print "pathSec1 hypothesis discrepancy distance:", juncDist1
		print "pathSec2 hypothesis discrepancy distance:", juncDist2


		"""
		if isInterior1 and isInterior2:
			if juncDist1 < juncDist2:
				secP1 = pathPoints2[0]
				secP2 = pathPoints2[frontDepI]
			else:
				secP1 = pathPoints2[backDepI]
				secP2 = pathPoints2[len(distances)-1]
		
		elif isExist1 and isInterior1:
			secP1 = pathPoints2[0]
			secP2 = pathPoints2[frontDepI]
			
		elif isExist2 and isInterior2:
			secP1 = pathPoints2[backDepI]
			secP2 = pathPoints2[len(distances)-1]
		
		else:
			print "Warning: external departure found"
			if isExist1:
				secP1 = pathPoints2[0]
				secP2 = pathPoints2[frontDepI]
				
			elif isExist2:
				secP1 = pathPoints2[backDepI]
				secP2 = pathPoints2[len(distances)-1]
		"""

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in path2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(0.5,0.5,1.0))
	
			if True:	
				P1 = pathPoints2[0]
				P2 = pathPoints2[frontDepI]
				
				P3 = pathPoints2[backDepI]
				P4 = pathPoints2[-1]

				" 0, frontDepI "
				xP = [P1[0], P2[0]]
				yP = [P1[1], P2[1]]
				pylab.scatter(xP,yP, color='b')		

				" backDepI, -1 "
				xP = [P3[0], P4[0]]
				yP = [P3[1], P4[1]]
				pylab.scatter(xP,yP, color='g')		
	
				
				#xP = [P2[0], P3[0]]
				#yP = [P2[1], P3[1]]
				#xP = [P1[0], P2[0], P3[0], P4[0]]
				#yP = [P1[1], P2[1], P3[1], P4[1]]
	
			pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')		
	
			xP = []
			yP = []
			for p in path1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(1.0,0.5,0.5))
	
			#xP = []
			#yP = []
			#for p in globalBranchPoints:
			#	xP.append(p[0])
			#	yP.append(p[1])
			#pylab.plot(xP,yP, color='g')
	
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d: %d %d %d %d %d %3.2f %3.2f %3.2f %d %3.2f %3.2f %3.2f" % (self.numNodes, isExist1, isExist2, isInterior1, isInterior2, matchCount1, overlapSum1, angleSum1, juncDist1, matchCount2, overlapSum2, angleSum2, juncDist2))
			pylab.savefig("trimDeparture_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1

		secP1 = []
		secP2 = []

		if matchCount1 > 0 and matchCount2 > 0:
			if juncDist1 < juncDist2:
				secP1 = pathPoints2[0]
				secP2 = pathPoints2[frontDepI]
			else:
				secP1 = pathPoints2[backDepI]
				secP2 = pathPoints2[len(distances)-1]

		elif matchCount1 > 0:
			secP1 = pathPoints2[0]
			secP2 = pathPoints2[frontDepI]
			
		elif matchCount2 > 0:
			secP1 = pathPoints2[backDepI]
			secP2 = pathPoints2[len(distances)-1]
		
		else:
			print "no overlap detected!"
			raise
				
		

		#if isInterior1 and isInterior2:
		#	print "2 departures found"
		#	raise
		
		if len(secP1) == 0:
			print "no departures found"
			raise
		
			
		
		return secP1, secP2


	def getOverlapCondition(self, supportLine, nodeID):

		plotIter = False

		if len(supportLine) == 0:
			return 1e100

		node2 = self.nodeHash[nodeID]

		hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
	
		estPose2 = node2.getGlobalGPACPose()		
			
		#minMatchDist2 = 0.2
		#minMatchDist2 = 0.5
		minMatchDist2 = 1.0
		#minMatchDist2 = 2.0
			
		" set the initial guess "
		poseOrigin = Pose(estPose2)
		
		localSupport = []
		for pnt in supportLine:
			localSupport.append(poseOrigin.convertGlobalToLocal(pnt))
	
		supportSpline = SplineFit(localSupport, smooth=0.1)		
		vecPoints1 = supportSpline.getUniformSamples()
		supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			
		
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2)
				
		support_pairs = []
		#for i in range(len(poly2)):
		for i in range(len(vecPoints2)):
			
			p_2 = vecPoints2[i]
	
			" for every transformed point of A, find it's closest neighbor in B "
			#_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

			try:
				p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
	
				#if gen_icp.isInCircle(p_1, radius2, center2):
	
				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = supportPoints[minI][2]
					#C1 = p_1[2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					#support_pairs.append([points2[i],p_1,C2,C1])
					support_pairs.append([points2[i],supportPoints[minI],C2,C1])
			except:
				pass

		cost = 0.0
		if len(support_pairs) == 0:
			cost = 1e100
		else:
			vals = []
			sum1 = 0.0
			for pair in support_pairs:
		
				a = pair[0]
				b = pair[1]
				Ca = pair[2]
				Cb = pair[3]
		
				ax = a[0]
				ay = a[1]		
				bx = b[0]
				by = b[1]
		
				c11 = Ca[0][0]
				c12 = Ca[0][1]
				c21 = Ca[1][0]
				c22 = Ca[1][1]
						
				b11 = Cb[0][0]
				b12 = Cb[0][1]
				b21 = Cb[1][0]
				b22 = Cb[1][1]	
			
				val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
				
				vals.append(val)
				sum1 += val
				
			cost = sum1 / len(support_pairs)
			

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in supportPoints:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')
	
			xP = []
			yP = []
			for p in poly2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			for pair in support_pairs:
				p1 = pair[0]
				p2 = pair[1]
				xP = [p1[0],p2[0]]
				yP = [p1[1],p2[1]]
				pylab.plot(xP,yP)
					
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("nodeID %d, cost = %f, count = %d" % (nodeID, cost, len(support_pairs)))
			pylab.savefig("overlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 1e100

		return cost

	def addNode(self, newNode):
		
		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1

	def correctNode2(self, nodeID):

		"FIXME:  add latest changes from loadNewNode() "

		newNode = self.nodeHash[nodeID]		

		return self.integrateNode(newNode, nodeID)

						
	def loadNewNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		return self.integrateNode(newNode, nodeID)



	def integrateNode(self, newNode, nodeID):

		" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
		direction = newNode.travelDir


		
		if nodeID > 0:
			
			" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
			if nodeID % 2 == 0:

				if self.nodeHash[nodeID-2].getNumLeafs() > 2 or self.nodeHash[nodeID].getNumLeafs() > 2:
					#transform1, covE1, hist1 = self.makeMultiJunctionMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )
					transform, covE, hist1 = self.makeMultiJunctionMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )

				else:
				
					transform1, covE1, hist1 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )
					if hist1[2] > 0 or hist1[1] > 10 and hist1[2] == 0:
					
						transform2, covE2, hist2 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = not direction )
	
	
						if hist1[2] < hist2[2]:
							transform = transform1
							covE = covE1
						else:
							if hist1[1] <= hist2[1]:
								transform = transform1
								covE = covE1
							else:
								transform = transform2
								covE = covE2
					else:
						transform = transform1
						covE = covE1
					
				self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)
	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" ODD NUMBER POSES RECEIVE INPLACE CONSTRAINT WITH EVEN PAIR "
				" PERFORM MEDIAL OVERLAP WITH PREVIOUS ODD NUMBER POSE "
			else:
			
				" FIXME:  apply this support check to spliced paths "
				supportLine = self.paths[self.currPath]
				
				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				#self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY2)
				self.inplace1.append([nodeID-1,nodeID,transform,covE])
				transform1 = transform
				offset1 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE1 = covE
				
				if len(supportLine) == 0:
					resultSum1 = 1e100
				else:
					resultSum1 = self.checkSupport(nodeID-1, nodeID, offset1, supportLine)

				node1 = self.nodeHash[nodeID-1]
				node2 = self.nodeHash[nodeID]
				
				" create the ground constraints "
				gndGPAC1Pose = node1.getGndGlobalGPACPose()
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = node2.getGndGlobalGPACPose()
				offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
				transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
				offset2 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE2 = covE
				#resultSum2 = self.checkSupport(nodeID-1, nodeID, offset2, supportLine)

	
				#self.addPriorityEdge([nodeID-1,nodeID,transform,self.E_inplace], INPLACE_PRIORITY3)				
				self.inplace2.append([nodeID-1,nodeID,transform,covE])

				
				if self.nodeHash[nodeID-1].getNumLeafs() > 2 or self.nodeHash[nodeID].getNumLeafs() > 2:
					transform, covE, overHist = self.makeMultiJunctionMedialOverlapConstraint(nodeID-1, nodeID, isMove = False, inPlace = False, isForward = direction )
				else:
					transform, covE, overHist = self.makeMedialOverlapConstraint(nodeID-1, nodeID, isMove = False, inPlace = True, isForward = direction)
				
				
				transform3 = transform
				offset3 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE3 = covE
				if len(supportLine) == 0:
					resultSum3 = 1e100
				else:
					resultSum3 = self.checkSupport(nodeID-1, nodeID, offset3, supportLine)
				#self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY)
				self.inplace3.append([nodeID-1,nodeID,transform,covE])

				print "INPLACE sums:", resultSum1, resultSum3

				if len(supportLine) > 0 and resultSum1 < resultSum3 and resultSum3 - resultSum1 > 100.0:
					self.addPriorityEdge([nodeID-1,nodeID,transform1,covE1], INPLACE_PRIORITY)
				else:
					self.addPriorityEdge([nodeID-1,nodeID,transform3,covE3], INPLACE_PRIORITY)
					


				if nodeID > 2:
					if self.nodeHash[nodeID-2].getNumLeafs() > 2 or self.nodeHash[nodeID].getNumLeafs() > 2:
						transform, covE, hist1 = self.makeMultiJunctionMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )
					else:
						transform1, covE1, hist1 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction)
	
						if hist1[2] > 0 or hist1[1] > 10 and hist1[2] == 0:
							transform2, covE2, hist2 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = not direction )
	
							if hist1[2] < hist2[2]:
								transform = transform1
								covE = covE1
							else:
								if hist1[1] <= hist2[1]:
									transform = transform1
									covE = covE1
								else:
									transform = transform2
									covE = covE2
						else:
							transform = transform1
							covE = covE1
						
					
					self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)

	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()


		" CHECK FOR A BRANCHING EVENT "
		
		if self.numNodes >= 4 and self.numNodes % 2 == 0:


			" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
			for k in range(len(self.nodeSets)):
				print "computing path for node set", k, ":", self.nodeSets[k]
				self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])
					
			self.drawPathAndHull()
			
			" DETECT BRANCHING EVENTS FOR THE 2 NODES OF LAST STOP "
			" AFTER ODD-NUMBER OVERLAP OF CURRENT STOP HAS IRONED OUT ERRORS "
			nodeID1 = self.numNodes-4
			nodeID2 = self.numNodes-3
		
			" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
			if len(self.paths[0]) == 0:
				" first nodes in path 0 "					
				self.nodeSets[0].append(nodeID1)
				self.nodeSets[0].append(nodeID2)

				" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
			else:
				
				" departure events for node 1 "
				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
	
				" departure events for node 2 "
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []

				
				" raw medial axis of union of path nodes "
				paths = {}
				for k in range(len(self.paths)):
					path = self.paths[k]
					if len(path) > 0:
						paths[k] = path
	
				print "paths:", len(paths)
				for k in range(len(paths)):
					print k, len(paths[k])
	
	
				" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
				self.trimmedPaths = self.trimPaths(paths)
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)
	
				
				" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
					
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
				
				print "orderedPathIDs1:", orderedPathIDs1
				print "orderedPathIDs2:", orderedPathIDs2
				
				
				" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
				
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					
				print "node departures", nodeID1, ":", departures1
				print "node  interiors", nodeID1, ":", interiors1
				print "node departures", nodeID2, ":", departures2
				print "node  interiors", nodeID2, ":", interiors2
	
	
				" new junction finding logic "
				" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
				" if a terminal departure exists that is internal, than we have a new junction "
				DISC_THRESH = 0.2

				" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
				frontExist1 = departures1[0][0]
				backExist1 =  departures1[-1][1]
				frontInterior1 = interiors1[0][0]
				backInterior1 = interiors1[-1][1]

				foreTerm1 = frontInterior1 and frontExist1
				backTerm1 = backInterior1 and backExist1
				
				" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
				discForeTerm1 = distances1[0][0] > DISC_THRESH
				discBackTerm1 = distances1[-1][1] > DISC_THRESH
				
				foreAngle1 = depAngles1[0][0]
				backAngle1 = depAngles1[-1][1]
				
				" NODE2: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
				frontExist2 = departures2[0][0]
				backExist2 =  departures2[-1][1]
				frontInterior2 = interiors2[0][0]
				backInterior2 = interiors2[-1][1]

				foreTerm2 = frontInterior2 and frontExist2
				backTerm2 = backInterior2 and backExist2


				discForeTerm2 = distances2[0][0] > DISC_THRESH
				discBackTerm2 = distances2[-1][1] > DISC_THRESH

				foreAngle2 = depAngles2[0][0]
				backAngle2 = depAngles2[-1][1]

				frontAngDiff = diffAngle(foreAngle1, foreAngle2)
				backAngDiff = diffAngle(backAngle1, backAngle2)


				print "pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
				print "pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2

				" check for the cases that internal departure may have a departure point discrepancy, "
				" then we should perform a pose adjustment and recompute departures "
				print "checking discrepancy in departure:", nodeID1, nodeID2, foreTerm1, discForeTerm1, backTerm1, discBackTerm1, foreTerm2, discForeTerm2, backTerm2, discBackTerm2

				if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

					print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2

					#if discForeTerm1 or discBackTerm1:
					if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1:


						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]

						splicedPaths1 = self.splicePathIDs(orderedPathIDs1)

						totalGuesses = []
						for path in splicedPaths1:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
							angDiff = abs(estPose[2]-guessPose1[2])
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID1].setGPACPose(guessPose)

					elif foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]


						splicedPaths2 = self.splicePathIDs(orderedPathIDs2)

						totalGuesses = []
						for path in splicedPaths2:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
							angDiff = abs(estPose[2]-guessPose1[2])
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID2].setGPACPose(guessPose)

						#if foreTerm2 and discForeTerm2:
						#	pathID = orderedPathIDs2[0]
						#elif backTerm2 and discBackTerm2:
						#	pathID = orderedPathIDs2[-1]
							
						#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
						#self.nodeHash[nodeID2].setGPACPose(guessPose1)

					departures1 = []
					interiors1 = []
					depPoints1 = []
					distances1 = []
					depAngles1 = []
		
					departures2 = []
					interiors2 = []
					depPoints2 = []
					distances2 = []
					depAngles2 = []
					
					" the overlapping paths are computed from the initial guess of position "
					orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
					orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
					
					print "orderedPathIDs1:", orderedPathIDs1
					print "orderedPathIDs2:", orderedPathIDs2
					
					" now compute whether there are departure points after we have guessed a better position in synch with the paths "
					for pathID in orderedPathIDs1:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
						departures1.append([isExist1,isExist2])
						interiors1.append([isInterior1, isInterior2])
						depPoints1.append([departurePoint1, departurePoint2])
						distances1.append([discDist1, discDist2])
						depAngles1.append([depAngle1, depAngle2])
					
					for pathID in orderedPathIDs2:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
						departures2.append([isExist1,isExist2])
						interiors2.append([isInterior1, isInterior2])
						depPoints2.append([departurePoint1, departurePoint2])
						distances2.append([discDist1, discDist2])
						depAngles2.append([depAngle1, depAngle2])
						
					print "node", nodeID1, ":", departures1, interiors1
					print "node", nodeID2, ":", departures2, interiors2
		
					" new junction finding logic "
					" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
					" if a terminal departure exists that is internal, than we have a new junction "
					
					frontExist1 = departures1[0][0]
					backExist1 =  departures1[-1][1]
					frontInterior1 = interiors1[0][0]
					backInterior1 = interiors1[-1][1]

					foreTerm1 = frontInterior1 and frontExist1
					backTerm1 = backInterior1 and backExist1
					
					
					discForeTerm1 = distances1[0][0] > DISC_THRESH
					discBackTerm1 = distances1[-1][1] > DISC_THRESH
					
					foreAngle1 = depAngles1[0][0]
					backAngle1 = depAngles1[-1][1]
					
					#foreTerm2 = departures2[0][0] and interiors2[0][0]
					#backTerm2 = departures2[-1][1] and interiors2[-1][1]
					frontExist2 = departures2[0][0]
					backExist2 =  departures2[-1][1]
					frontInterior2 = interiors2[0][0]
					backInterior2 = interiors2[-1][1]

					foreTerm2 = frontInterior2 and frontExist2
					backTerm2 = backInterior2 and backExist2

	
					discForeTerm2 = distances2[0][0] > DISC_THRESH
					discBackTerm2 = distances2[-1][1] > DISC_THRESH

					foreAngle2 = depAngles2[0][0]
					backAngle2 = depAngles2[-1][1]
					

					frontAngDiff = diffAngle(foreAngle1, foreAngle2)
					backAngDiff = diffAngle(backAngle1, backAngle2)

					print "pass2:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
					print "pass2:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2
					
				print "frontAngDiff:", frontAngDiff
				print "backAngDiff:", backAngDiff



				newPaths = []
	
				frontExist1 = departures1[0][0]
				frontInterior1 = interiors1[0][0]
		
				frontExist2 = departures2[0][0]
				frontInterior1 = interiors2[0][0]
				
				depAngle1 = depAngles1[0][0]
				depAngle2 = depAngles2[0][0]
				
				depPoint1 = depPoints1[0][0]
				depPoint2 = depPoints2[0][0]
				
				parentPathID1 = orderedPathIDs1[0]
				parentPathID2 = orderedPathIDs2[0]
	
	
				isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)

				print "determineBranch:"
				print isBranch
				print pathBranchIDs
	
				if isBranch[0]:
					orderedPathIDs1.insert(0,pathBranchIDs[0])
					departures1.insert(0, [False, False])
					interiors1.insert(0, [False, False])
					depPoints1.insert(0, [None, None])
	
				if isBranch[1]:
					orderedPathIDs2.insert(0,pathBranchIDs[1])
					departures2.insert(0, [False, False])
					interiors2.insert(0, [False, False])
					depPoints2.insert(0, [None, None])
					
				for pathID in pathBranchIDs:
					if pathID != -1 and newPaths.count(pathID) == 0:
						newPaths.append(pathID)




				backExist1 = departures1[-1][1]
				backInterior1 = interiors1[-1][1]
		
				backExist2 = departures2[-1][1]
				backInterior1 = interiors2[-1][1]
				
				depAngle1 = depAngles1[-1][1]
				depAngle2 = depAngles2[-1][1]
				
				depPoint1 = depPoints1[-1][1]
				depPoint2 = depPoints2[-1][1]
	
				parentPathID1 = orderedPathIDs1[-1]
				parentPathID2 = orderedPathIDs2[-1]
				print "self.pathParents:", self.pathParents
				print "self.nodeSets:", self.nodeSets
				print "newPaths:", newPaths
				print "self.numPaths:", self.numPaths
	
	
				isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)
				print "determineBranch:"
				print isBranch
				print pathBranchIDs
	
				if isBranch[0]:
					orderedPathIDs1.append(0,pathBranchIDs[0])
					departures1.append([False, False])
					interiors1.append([False, False])
					depPoints1.append([None, None])
	
				if isBranch[1]:
					orderedPathIDs2.append(pathBranchIDs[1])
					departures2.append([False, False])
					interiors2.append([False, False])
					depPoints2.append([None, None])
	
				for pathID in pathBranchIDs:
					if pathID != -1 and newPaths.count(pathID) == 0:
						newPaths.append(pathID)
				print "self.pathParents:", self.pathParents
				print "self.nodeSets:", self.nodeSets
				print "newPaths:", newPaths
				print "self.numPaths:", self.numPaths
				


				" determine which paths are leaves "
				print "self.pathParents:", self.pathParents
				isAParent = [False for k in range(self.numPaths)]
				for k in orderedPathIDs1:
					print "index:", k
					currPath = self.pathParents[k]
					currParent = currPath[0]
					if currParent != None:
						isAParent[currParent] = True
				
				"add nodes to paths that are the leaves "
				for pathID in orderedPathIDs1:
					if not isAParent[pathID]:				
						self.nodeSets[pathID].append(nodeID1)

				isAParent = [False for k in range(self.numPaths)]
				for k in orderedPathIDs2:
					print "index:", k
					currPath = self.pathParents[k]
					currParent = currPath[0]
					if currParent != None:
						isAParent[currParent] = True

				for pathID in orderedPathIDs2:
					if not isAParent[pathID]:				
						self.nodeSets[pathID].append(nodeID2)

				paths = {}
				for k in range(len(self.paths)):
					path = self.paths[k]
					if len(path) > 0:
						paths[k] = path	
				self.trimmedPaths = self.trimPaths(paths)

				" perform path constraints only if we have not created a new path "
				isNew1 = False
				isNew2 = False
				for pathID in newPaths:
					if orderedPathIDs1.count(pathID) > 0:
						isNew1 = True
					if orderedPathIDs2.count(pathID) > 0:
						isNew2 = True
				
	
				
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)



				"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
				"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
				"3)  select choice with the lowest cost "
				"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
				
				" nodeID1:  is it in a junction or a single path? "
				
				
				if not isNew1:
					self.constrainToPaths(nodeID1, orderedPathIDs1, departures1, interiors1, depPoints1)
					
				if not isNew2:
					self.constrainToPaths(nodeID2, orderedPathIDs2, departures2, interiors2, depPoints2)

				self.mergePriorityConstraints()

			for k in range(len(self.nodeSets)):
				print "computing path for node set", k, ":", self.nodeSets[k]
				self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])

			self.drawPathAndHull()
			
			return


	def addNewPath(self, parentID, branchNodeID, localJunctionPose):
				
		" FIXME:  Branch direction does not seem to be used properly yet.  Only used in getOverlapDeparture() "
		self.pathParents.append([parentID, branchNodeID, localJunctionPose, True])

		newPathID = self.numPaths
		self.nodeSets[newPathID] = []
		self.numPaths += 1

		return newPathID
			
	def constrainToPaths(self, nodeID, orderedPathIDs, departures, interiors, depPoints):
			
		"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
		"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
		"3)  select choice with the lowest cost "
		"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
		
		" nodeID1:  is it in a junction or a single path? "
		
		splicedPaths1 = self.splicePathIDs(orderedPathIDs)

		if len(orderedPathIDs) == 1:

			" in a single path "

			pathID = orderedPathIDs[0]
			nodeSet = self.nodeSets[pathID]
			
			" check that we're not the only pose in the path "
			" does it contain at least one node that is not nodeID1 or nodeID2 "
			doConstraint = False
			
			" companion node iD "
			if nodeID % 2 == 0:
				pairNodeID = nodeID + 1
			else:
				pairNodeID = nodeID - 1
				 
			if len(nodeSet) == 1:						
				if nodeSet.count(nodeID) == 0:
					doConstraint = True
			
			elif len(nodeSet) == 2:						
				if nodeSet.count(nodeID) == 0 or nodeSet.count(pairNodeID) == 0:
					doConstraint = True
			
			elif len(nodeSet) > 2:
				doConstraint = True

			if doConstraint:
				
				" control point nearest the GPAC origin "
				globalPoint1 = self.nodeHash[nodeID].getGlobalGPACPose()[:2]
					
				guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
				self.nodeHash[nodeID].setGPACPose(guessPose1)

				" make path constraint "
				print "pathID: addPathConstraints(", pathID, nodeID
				self.addPathConstraints(self.nodeSets[pathID], nodeID)
			
		else:
			" in at least one junction "
			
			" we don't intersect the junction sufficiently "
			" add constraints to most overlapped path"
			#if interiors.count(True) == 0:
			" no recognized internal departure point "
			if [item for inner_list in interiors for item in inner_list].count(True) == 0:


				" determine which paths are leaves "
				isAChild = [True for k in range(self.numPaths)]
				for k in orderedPathIDs:
					currPath = self.pathParents[k]
					currParent = currPath[0]
					if currParent != None and orderedPathIDs.count(currParent) > 0:
						isAChild[currParent] = False

				print "isAChild:", isAChild

				minPathID = 0
				minSum = 1e100
				for pathID in orderedPathIDs:
					
					if isAChild[pathID]:
						sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
						
						if sum1 < minSum:
							minPathID = pathID
							minSum = sum1
						print "overlap sum", pathID, "=", sum1						

				print "maximum overlap path is", minPathID

				nodeSet = self.nodeSets[minPathID]
				
				" check that we're not the only pose in the path "
				" does it contain at least one node that is not nodeID or nodeID2 "
				doConstraint = False

				" companion node iD "
				if nodeID % 2 == 0:
					pairNodeID = nodeID + 1
				else:
					pairNodeID = nodeID - 1

				if len(nodeSet) == 1:						
					if nodeSet.count(nodeID) == 0:
						doConstraint = True
				
				elif len(nodeSet) == 2:						
					if nodeSet.count(nodeID) == 0 or nodeSet.count(pairNodeID) == 0:
						doConstraint = True
				
				elif len(nodeSet) > 2:
					doConstraint = True

				if doConstraint:
					
					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID].getGlobalGPACPose()[:2]

					totalGuesses = []
					
					for path in splicedPaths1:

						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					
					guessPose1 = totalGuesses[0][2]
					
						
					#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, self.trimmedPaths[minPathID], localPoint1, localPoint1)
					self.nodeHash[nodeID].setGPACPose(guessPose1)

					" make path constraint "
					print "minPathID: addPathConstraints(", minPathID, nodeID
					self.addPathConstraints(self.nodeSets[minPathID], nodeID)






			
			elif len(orderedPathIDs) == 2:
				
				" one junction point "
				pathID1 = orderedPathIDs[0]
				pathID2 = orderedPathIDs[1]

				#[pathID, branchNodeID, poseOrigin.convertGlobalToLocal(globalJunctionPoint)]
				
				" find parent-child relationship "
				val1 = self.pathParents[pathID1]
				val2 = self.pathParents[pathID2]


				if val1[0] == None:
					parentPath = pathID1
					childPath = pathID2
					
				elif val2[0] == None:					
					parentPath = pathID2
					childPath = pathID1
					
				else:
					parentID1 = val1[0]
					parentID2 = val2[0]
					
					if parentID1 == pathID2:
						parentPath = pathID2
						childPath = pathID1
					elif parentID2 == pathID1:
						parentPath = pathID1
						childPath = pathID2
					else:
						raise

				" global path point "
				localJunctionNodeID = self.pathParents[childPath][1]
				localPathJunctionPoint = self.pathParents[childPath][2]						
				
				localJunctionNode = self.nodeHash[localJunctionNodeID]
				poseOrigin = Pose(localJunctionNode.getEstPose())
				globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
			
				print departures, interiors
				print pathID1, childPath, parentPath
			
				" new node departure point "
				" try childPath first "
				if childPath == pathID1:
					index1 = orderedPathIDs.index(pathID1)			
					isDepExists1 = departures[index1][1] and interiors[index1][1]
					depPoint1 = depPoints[index1][1]
					print "caseA:", index1, isDepExists1
				
				else:
					index2 = orderedPathIDs.index(pathID2)
					isDepExists1 = departures[index2][0] and interiors[index2][0]
					depPoint1 = depPoints[index2][0]
					print "caseB:", index2, isDepExists1

				" try parentPath first "
				if parentPath == pathID1:
					index1 = orderedPathIDs.index(pathID1)
					isDepExists2 = departures[index1][1] and interiors[index1][1]
					depPoint2 = depPoints[index1][1]
					print "caseC:", index1, isDepExists2
				else:
					index2 = orderedPathIDs.index(pathID2)
					isDepExists2 = departures[index2][0] and interiors[index2][0]
					depPoint2 = depPoints[index2][0]
					print "caseD:", index2, isDepExists2
				
				
				if isDepExists1:
					
					totalGuesses = []
					
					for path in splicedPaths1:

						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPathJunctionPoint, globalPathJunctionPoint)
						
						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					
					guessPose = totalGuesses[0][2]
					
					self.nodeHash[nodeID].setGPACPose(guessPose)

					" make path constraint "
					print "childPath: addPathConstraints(", childPath, nodeID
					self.addPathConstraints(self.nodeSets[childPath], nodeID)

				
				elif isDepExists2:
					totalGuesses = []
					
					for path in splicedPaths1:

						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPathJunctionPoint, depPoint2)

						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					
					guessPose = totalGuesses[0][2]

					self.nodeHash[nodeID].setGPACPose(guessPose)

					" make path constraint "
					print "childPath: addPathConstraints(", childPath, nodeID
					self.addPathConstraints(self.nodeSets[childPath], nodeID)

				
				else:
					print "no local junction point "
					raise
		
			else:
				" two or more junction points "
				pathID1_1 = orderedPathIDs[0]
				pathID1_2 = orderedPathIDs[1]
				pathID2_1 = orderedPathIDs[-1]
				pathID2_2 = orderedPathIDs[-2]

				" pick the newest path from where to sample the local junction "
				if pathID1_1 > pathID2_2:
					pathID1 = pathID1_1
					pathID2 = pathID1_2
				else:
					pathID1 = pathID2_1
					pathID2 = pathID2_2
					
				" find parent-child relationship "
				val1 = self.pathParents[pathID1]
				val2 = self.pathParents[pathID2]


				if val1[0] == None:
					parentPath = pathID1
					childPath = pathID2

				elif val2[0] == None:					
					parentPath = pathID2
					childPath = pathID1

				else:
					parentID1 = val1[0]
					parentID2 = val2[0]
					
					if parentID1 == pathID2:
						parentPath = pathID2
						childPath = pathID1
					elif parentID2 == pathID1:
						parentPath = pathID1
						childPath = pathID2
					else:
						raise

				if self.pathParents[childPath][0] != parentPath:
					print "WARNING: selected childPath", childPath, "does not have parent", parentPath, "but has instead", self.pathParents[childPath][0]

				" global path point "
				" we don't check, but we assume the parent here is the parent in our orderedPaths"
				localJunctionNodeID = self.pathParents[childPath][1]
				localPathJunctionPoint = self.pathParents[childPath][2]						
				
				localJunctionNode = self.nodeHash[localJunctionNodeID]
				poseOrigin = Pose(localJunctionNode.getEstPose())
				globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
			
				totalGuesses = []
				
				for path in splicedPaths1:

					" make the path constraints "								
					guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPathJunctionPoint, globalPathJunctionPoint)
					estPose = self.nodeHash[nodeID].getGlobalGPACPose()
					angDiff = abs(estPose[2]-guessPose1[2])
					totalGuesses.append((angDiff, cost1, guessPose1))
				
				totalGuesses.sort()
				
				guessPose = totalGuesses[0][2]
				
				self.nodeHash[nodeID].setGPACPose(guessPose)

				" make path constraint "
				print "childPath: addPathConstraints(", childPath, nodeID
				self.addPathConstraints(self.nodeSets[childPath], nodeID)

	def determineBranch(self, nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2):
		
		"""
		1)  both nodes departing
			- are same departure  ( equal departing angle and are on top of each other)  1 new path
			- are different departures ( different departing angle and are not necessarily on top of each other) 2 new paths
			
		2) neither node is departing
			- are on same path
			
		3) one node is departing only
			- is only one departure ( other node has no departure, or has external departure but not same angle )
			- is a false departure  ( falseness happens if bad map data, ignore )
			- are both departing ( should have at least external departure on both, compare departure angle )
		
		Now, we have the departing angles, but we need to compare them with each other.	
		"""

		
		
		foreTerm1 = frontInterior1 and frontExist1
		foreTerm2 = frontInterior2 and frontExist2

		frontAngDiff = diffAngle(depAngle1, depAngle2)

		pathBranchIDs = [-1, -1]
		isBranch = [False, False]
		
		" 60 degree threshold "
		ANG_THRESH = 1.047
		#ANG_THRESH = 0.3
		
		if foreTerm1 or foreTerm2:

			if foreTerm1 and foreTerm2:
				" both departing case: check if same or different "
				if fabs(frontAngDiff) < ANG_THRESH:
					" are on top of each other and are branching to the same path "

					if self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1) and self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2):

						
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						pathBranchIDs[1] = newPathID

						isBranch[0] = True
						isBranch[1] = True

						print "foreA"
				else:
					" these are branching to separate new paths "

					if self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1) and self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2):


						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID1 = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID1
						isBranch[0] = True



						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID2 = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID2
						isBranch[1] = True
						print "foreB"
			
			elif foreTerm1 and not foreTerm2:
				
				if self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1):
					
					if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
						" both have same departure "
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						pathBranchIDs[1] = newPathID
						isBranch[0] = True
						isBranch[1] = True

						print "foreC"
						
					else:
						" foreTerm1 has unique departure "	
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True

						print "foreD"

			elif foreTerm2 and not foreTerm1:

				if self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2):

					if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

						" both have same departure "
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						pathBranchIDs[1] = newPathID
						isBranch[0] = True
						isBranch[1] = True

						print "foreE"

					else:
						" foreTerm2 has unique departure "	
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True

						print "foreF"
			else:
				print "foreG"
				" no new departure, just add nodes to the leaf paths "
				pathID = parentPathID1
				branchNodeID = nodeID1
				globalJunctionPoint = depPoint1
				
				if parentPathID1 != parentPathID2:
					print "departing path IDs of two colocated nodes are not the same"
					raise
		
		
		return isBranch, pathBranchIDs
		
		"""
		pathID = parentPathID1
		branchNodeID = nodeID1
		globalJunctionPoint = depPoint1
		depAng = depAngle1
		junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
		orderedPathIDs1.insert(0,self.numPaths)
		departures1.insert(0, [False, False])
		interiors1.insert(0, [False, False])
		depPoints1.insert(0, [None, None])
	
		orderedPathIDs2.insert(0,self.numPaths)
		departures2.insert(0, [False, False])
		interiors2.insert(0, [False, False])
		depPoints2.insert(0, [None, None])
	
		poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
		self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
		self.nodeSets[self.numPaths] = []
		newPaths.append(self.numPaths)
		self.numPaths += 1
		"""



	def checkForBranches(self, nodeID1, nodeID2):

		" CHECK FOR A BRANCHING EVENT FROM LEAF "

		" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
		if len(self.paths[0]) == 0:
			" first nodes in path 0 "					
			self.nodeSets[0].append(nodeID1)
			self.nodeSets[0].append(nodeID2)

			" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
		else:
			
			" departure events for node 1 "
			departures1 = []
			interiors1 = []
			depPoints1 = []
			distances1 = []
			depAngles1 = []

			" departure events for node 2 "
			departures2 = []
			interiors2 = []
			depPoints2 = []
			distances2 = []
			depAngles2 = []

			
			" raw medial axis of union of path nodes "
			paths = {}
			for k in range(len(self.paths)):
				path = self.paths[k]
				if len(path) > 0:
					paths[k] = path

			print "paths:", len(paths)
			for k in range(len(paths)):
				print k, len(paths[k])


			" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
			self.trimmedPaths = self.trimPaths(paths)
			self.drawTrimmedPaths(self.trimmedPaths)
			print "trimmed paths:", len(self.trimmedPaths)

			
			" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
				
			" the overlapping paths are computed from the initial guess of position "
			orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
			orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
			
			print "orderedPathIDs1:", orderedPathIDs1
			print "orderedPathIDs2:", orderedPathIDs2
			
			
			" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
			for pathID in orderedPathIDs1:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
				departures1.append([isExist1,isExist2])
				interiors1.append([isInterior1, isInterior2])
				depPoints1.append([departurePoint1, departurePoint2])
				distances1.append([discDist1, discDist2])
				depAngles1.append([depAngle1, depAngle2])
			
			for pathID in orderedPathIDs2:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
				departures2.append([isExist1,isExist2])
				interiors2.append([isInterior1, isInterior2])
				depPoints2.append([departurePoint1, departurePoint2])
				distances2.append([discDist1, discDist2])
				depAngles2.append([depAngle1, depAngle2])
				
			print "node departures", nodeID1, ":", departures1
			print "node  interiors", nodeID1, ":", interiors1
			print "node departures", nodeID2, ":", departures2
			print "node  interiors", nodeID2, ":", interiors2


			" new junction finding logic "
			" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
			" if a terminal departure exists that is internal, than we have a new junction "
			DISC_THRESH = 0.2

			" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
			frontExist1 = departures1[0][0]
			backExist1 =  departures1[-1][1]
			frontInterior1 = interiors1[0][0]
			backInterior1 = interiors1[-1][1]

			foreTerm1 = frontInterior1 and frontExist1
			backTerm1 = backInterior1 and backExist1
			
			" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
			discForeTerm1 = distances1[0][0] > DISC_THRESH
			discBackTerm1 = distances1[-1][1] > DISC_THRESH
			
			foreAngle1 = depAngles1[0][0]
			backAngle1 = depAngles1[-1][1]
			
			" NODE2: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
			frontExist2 = departures2[0][0]
			backExist2 =  departures2[-1][1]
			frontInterior2 = interiors2[0][0]
			backInterior2 = interiors2[-1][1]

			foreTerm2 = frontInterior2 and frontExist2
			backTerm2 = backInterior2 and backExist2


			discForeTerm2 = distances2[0][0] > DISC_THRESH
			discBackTerm2 = distances2[-1][1] > DISC_THRESH

			foreAngle2 = depAngles2[0][0]
			backAngle2 = depAngles2[-1][1]

			frontAngDiff = diffAngle(foreAngle1, foreAngle2)
			backAngDiff = diffAngle(backAngle1, backAngle2)


			print "pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
			print "pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2

			" check for the cases that internal departure may have a departure point discrepancy, "
			" then we should perform a pose adjustment and recompute departures "
			print "checking discrepancy in departure:", nodeID1, nodeID2, foreTerm1, discForeTerm1, backTerm1, discBackTerm1, foreTerm2, discForeTerm2, backTerm2, discBackTerm2

			#if discForeTerm1 or discBackTerm1 or discForeTerm2 or discBackTerm2:
			if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

				print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2

				#if discForeTerm1 or discBackTerm1:
				if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1:


					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]

					splicedPaths1 = self.splicePathIDs(orderedPathIDs1)

					totalGuesses = []
					for path in splicedPaths1:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					guessPose = totalGuesses[0][2]
					self.nodeHash[nodeID1].setGPACPose(guessPose)


					#if foreTerm1 and discForeTerm1:
					#	pathID = orderedPathIDs1[0]
					#elif backTerm1 and discBackTerm1:
					#	pathID = orderedPathIDs1[-1]
					
						
					#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
					#self.nodeHash[nodeID1].setGPACPose(guessPose1)
						
				#elif discForeTerm2 or discBackTerm2:
				elif foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]


					splicedPaths2 = self.splicePathIDs(orderedPathIDs2)

					totalGuesses = []
					for path in splicedPaths2:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					guessPose = totalGuesses[0][2]
					self.nodeHash[nodeID2].setGPACPose(guessPose)

					#if foreTerm2 and discForeTerm2:
					#	pathID = orderedPathIDs2[0]
					#elif backTerm2 and discBackTerm2:
					#	pathID = orderedPathIDs2[-1]
						
					#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
					#self.nodeHash[nodeID2].setGPACPose(guessPose1)

				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
	
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []
				
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
				
				print "orderedPathIDs1:", orderedPathIDs1
				print "orderedPathIDs2:", orderedPathIDs2
				
				" now compute whether there are departure points after we have guessed a better position in synch with the paths "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
				
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					
				print "node", nodeID1, ":", departures1, interiors1
				print "node", nodeID2, ":", departures2, interiors2
	
				" new junction finding logic "
				" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
				" if a terminal departure exists that is internal, than we have a new junction "
				
				frontExist1 = departures1[0][0]
				backExist1 =  departures1[-1][1]
				frontInterior1 = interiors1[0][0]
				backInterior1 = interiors1[-1][1]

				foreTerm1 = frontInterior1 and frontExist1
				backTerm1 = backInterior1 and backExist1
				
				
				discForeTerm1 = distances1[0][0] > DISC_THRESH
				discBackTerm1 = distances1[-1][1] > DISC_THRESH
				
				foreAngle1 = depAngles1[0][0]
				backAngle1 = depAngles1[-1][1]
				
				#foreTerm2 = departures2[0][0] and interiors2[0][0]
				#backTerm2 = departures2[-1][1] and interiors2[-1][1]
				frontExist2 = departures2[0][0]
				backExist2 =  departures2[-1][1]
				frontInterior2 = interiors2[0][0]
				backInterior2 = interiors2[-1][1]

				foreTerm2 = frontInterior2 and frontExist2
				backTerm2 = backInterior2 and backExist2


				discForeTerm2 = distances2[0][0] > DISC_THRESH
				discBackTerm2 = distances2[-1][1] > DISC_THRESH

				foreAngle2 = depAngles2[0][0]
				backAngle2 = depAngles2[-1][1]
				

				frontAngDiff = diffAngle(foreAngle1, foreAngle2)
				backAngDiff = diffAngle(backAngle1, backAngle2)

				print "pass2:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
				print "pass2:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2
				
			print "frontAngDiff:", frontAngDiff
			print "backAngDiff:", backAngDiff


			newPaths = []

			frontExist1 = departures1[0][0]
			frontInterior1 = interiors1[0][0]
	
			frontExist2 = departures2[0][0]
			frontInterior1 = interiors2[0][0]
			
			depAngle1 = depAngles1[0][0]
			depAngle2 = depAngles2[0][0]
			
			depPoint1 = depPoints1[0][0]
			depPoint2 = depPoints2[0][0]
			
			parentPathID1 = orderedPathIDs1[0]
			parentPathID2 = orderedPathIDs2[0]


			isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)

			if isBranch[0]:
				orderedPathIDs1.insert(0,pathBranchIDs[0])
				departures1.insert(0, [False, False])
				interiors1.insert(0, [False, False])
				depPoints1.insert(0, [None, None])

			if isBranch[1]:
				orderedPathIDs2.insert(0,pathBranchIDs[1])
				departures2.insert(0, [False, False])
				interiors2.insert(0, [False, False])
				depPoints2.insert(0, [None, None])
				
			for pathID in pathBranchIDs:
				if pathID != -1 and newPaths.count(pathID) == 0:
					newPaths.append(pathID)
			
			#for newPath in addThesePaths:
			#	self.pathParents.append(newPath)
			#	self.nodeSets[newPath[0]] = []
			#	newPaths.append(newPath[0])
			#	self.numPaths += 1


			backExist1 = departures1[-1][1]
			backInterior1 = interiors1[-1][1]
	
			backExist2 = departures2[-1][1]
			backInterior1 = interiors2[-1][1]
			
			depAngle1 = depAngles1[-1][1]
			depAngle2 = depAngles2[-1][1]
			
			depPoint1 = depPoints1[-1][1]
			depPoint2 = depPoints2[-1][1]

			parentPathID1 = orderedPathIDs1[-1]
			parentPathID2 = orderedPathIDs2[-1]


			isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)

			if isBranch[0]:
				orderedPathIDs1.append(0,pathBranchIDs[0])
				departures1.append([False, False])
				interiors1.append([False, False])
				depPoints1.append([None, None])

			if isBranch[1]:
				orderedPathIDs2.append(pathBranchIDs[1])
				departures2.append([False, False])
				interiors2.append([False, False])
				depPoints2.append([None, None])

			for pathID in pathBranchIDs:
				if pathID != -1 and newPaths.count(pathID) == 0:
					newPaths.append(pathID)
			
			#for newPath in addThesePaths:
			#	self.pathParents.append(newPath)
			#	self.nodeSets[newPath[0]] = []
			#	newPaths.append(newPath[0])
			#	self.numPaths += 1

	



			"""
			newPaths = []
			
			" 60 degree threshold "
			ANG_THRESH = 1.047
			#ANG_THRESH = 0.3
			
			if foreTerm1 or foreTerm2:

				if foreTerm1 and foreTerm2:
					" both departing case: check if same or different "
					if fabs(frontAngDiff) < ANG_THRESH:
						" are on top of each other and are branching to the same path "

						if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]) and self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):

							
							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1

							print "foreA"
					else:
						" these are branching to separate new paths "

						if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]) and self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):


							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1



							pathID = orderedPathIDs2[0]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[0][0]
							depAng = depAngles2[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])


							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreB"
				
				elif foreTerm1 and not foreTerm2:
					
					if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]):
						
						if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
							" both have same departure "
							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
							#depAng1 = depAngles1[0][0]
							#depAng2 = depAngles2[0][0]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreC"
							
						else:
							" foreTerm1 has unique departure "	
							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreD"

				elif foreTerm2 and not foreTerm1:

					if self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):

						if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

							" both have same departure "
							pathID = orderedPathIDs2[0]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[0][0]
							depAng = depAngles2[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreE"

						else:
							" foreTerm2 has unique departure "	
							pathID = orderedPathIDs2[0]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[0][0]
							depAng = depAngles2[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1

							print "foreF"
				else:
					print "foreG"
					" no new departure, just add nodes to the leaf paths "
					pathID = orderedPathIDs1[0]
					branchNodeID = nodeID1
					globalJunctionPoint = depPoints1[0][0]
					
					if orderedPathIDs1[0] != orderedPathIDs2[0]:
						print "departing path IDs of two colocated nodes are not the same"
						raise

			if backTerm1 or backTerm2:

				if backTerm1 and backTerm2:
					" both departing case: check if same or different "
					if fabs(backAngDiff) < ANG_THRESH:
						" are on top of each other and are branching to the same path "

						if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]) and self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):

							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backA"

					else:
						" these are branching to separate new paths "
						if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]) and self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):

							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1

							pathID = orderedPathIDs2[-1]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[-1][1]
							depAng = depAngles2[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backB"					
				elif backTerm1 and not backTerm2:
					
					if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]):
						
						if backExist2 and fabs(backAngDiff) < ANG_THRESH:
							" both have same departure "
							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backC"
						else:
							" backTerm1 has unique departure "	
							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backD"

				elif backTerm2 and not backTerm1:

					if self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):

						if backExist1 and fabs(backAngDiff) < ANG_THRESH:

							" both have same departure "
							pathID = orderedPathIDs2[-1]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[-1][1]
							depAng = depAngles2[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backE"
						else:
							" backTerm2 has unique departure "	
							self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1])

							pathID = orderedPathIDs2[-1]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[-1][1]
							depAng = depAngles2[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backF"
				else:
					print "backG"
					" no new departure, just add nodes to the leaf paths "
					pathID = orderedPathIDs1[-1]
					branchNodeID = nodeID1
					globalJunctionPoint = depPoints1[-1][1]
					if orderedPathIDs1[-1] != orderedPathIDs2[-1]:
						print "departing path IDs of two colocated nodes are not the same"
						raise
			"""





			" determine which paths are leaves "
			print "self.pathParents:", self.pathParents
			isAParent = [False for k in range(self.numPaths)]
			for k in orderedPathIDs1:
				print "index:", k
				currPath = self.pathParents[k]
				currParent = currPath[0]
				if currParent != None:
					isAParent[currParent] = True
			
			"add nodes to paths that are the leaves "
			for pathID in orderedPathIDs1:
				if not isAParent[pathID]:				
					self.nodeSets[pathID].append(nodeID1)

			isAParent = [False for k in range(self.numPaths)]
			for k in orderedPathIDs2:
				print "index:", k
				currPath = self.pathParents[k]
				currParent = currPath[0]
				if currParent != None:
					isAParent[currParent] = True

			for pathID in orderedPathIDs2:
				if not isAParent[pathID]:				
					self.nodeSets[pathID].append(nodeID2)

			paths = {}
			for k in range(len(self.paths)):
				path = self.paths[k]
				if len(path) > 0:
					paths[k] = path	
			self.trimmedPaths = self.trimPaths(paths)

			" perform path constraints only if we have not created a new path "
			isNew1 = False
			isNew2 = False
			for pathID in newPaths:
				if orderedPathIDs1.count(pathID) > 0:
					isNew1 = True
				if orderedPathIDs2.count(pathID) > 0:
					isNew2 = True
			


			self.drawTrimmedPaths(self.trimmedPaths)
			print "trimmed paths:", len(self.trimmedPaths)



			"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
			"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
			"3)  select choice with the lowest cost "
			"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
			
			" nodeID1:  is it in a junction or a single path? "

			if not isNew1:
				self.constrainToPaths(nodeID1, orderedPathIDs1, departures1, interiors1, depPoints1)
				
			if not isNew2:
				self.constrainToPaths(nodeID2, orderedPathIDs2, departures2, interiors2, depPoints2)

			self.mergePriorityConstraints()


		for k in range(len(self.nodeSets)):
			print "computing path for node set", k, ":", self.nodeSets[k]
			self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])

		self.drawPathAndHull()
		
		return
	
		


	def checkUniqueBranch(self, parentPathID, nodeID1, depAngle, depPoint):

		" check cartesian distance to similar junction points from parent path "

		" cartesian distance "
		DISC_THRESH = 1.0

		" 60 degree threshold "
		ANG_THRESH = 1.047
		
		for path in self.pathParents:
			if path[0] == parentPathID:
				
				junctionNodeID = path[1]
				localJunctionPoint = path[2]
				" check dist "

				poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
				junctionPoint = poseOrigin.convertLocalOffsetToGlobal(localJunctionPoint)
				
				dist = sqrt((depPoint[0]-junctionPoint[0])**2 + (depPoint[1]-junctionPoint[1])**2 )

				" check difference of tangential angle "
				angDiff = diffAngle(depAngle, junctionPoint[2])

				print "checkUniqueBranch(", parentPathID, ",", nodeID1, ",", depAngle, ",", depPoint, ") = ", junctionNodeID, dist, angDiff
		
				if dist < DISC_THRESH:
					if fabs(angDiff) < ANG_THRESH:
						
						print "DUPLICATE junction of", junctionPoint, "rejecting with differences of", dist, angDiff
						return False
		
		return True
				

	def mergePriorityConstraints(self):
		
		totalConstraints = []
		
		#print "merging:"
		
		" merge the constraint "
		for k, v in self.edgePriorityHash.items():
			
			if len(v) > 0:
		
				id1 = k[0]
				id2 = k[1]
				
				" take the highest priority constraints only "
				maxPriority = -1
				for const in v:
					thisPriority = const[2]
					if thisPriority > maxPriority:
						maxPriority = thisPriority
		
				if maxPriority == -1:
					raise
				
				for const in v:
					if const[2] == maxPriority:
						transform = const[0]
						covE = const[1]
						#print id1, id2, maxPriority
						totalConstraints.append([id1,id2,transform,covE])	
		
		#print "merging totalConstraints:", totalConstraints
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")
		
		self.edgeHash = {}
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1, node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]
		
		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])

	def computeMedialError(self, i, j, offset, minMatchDist = 2.0, tail1=0, tail2=0):


		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		#posture1 = node1.getStableGPACPosture()
		#posture2 = node2.getStableGPACPosture()		
		#hull1, medial1 = computeHullAxis(i, node1, tailCutOff = True)
		#hull2, medial2 = computeHullAxis(j, node2, tailCutOff = True)
		
		medial1 = node1.medialLongPaths[tail1]
		medial2 = node2.medialLongPaths[tail2]
		
		#estPose1 = node1.getGlobalGPACPose()		
		#estPose2 = node2.getGlobalGPACPose()
		
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medialSpline2 = SplineFit(medial2, smooth=0.1)

		points1 = gen_icp.addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
		points2 = gen_icp.addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)
	
		" transform pose 2 by initial offset guess "	
		" transform the new pose "
		points2_offset = []
		for p in points2:
			result = gen_icp.dispPoint(p, offset)		
			points2_offset.append(result)


		#costThresh = 0.004
		#minMatchDist = 2.0
		#lastCost = 1e100
		
		poly1 = []		
		for p in points1:
			poly1.append([p[0],p[1]])	
		
		" find the matching pairs "
		match_pairs = []
		
		" transformed points without associated covariance "
		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])	
		
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
		radius1, center1 = gen_icp.computeEnclosingCircle(points1)
	
		for i in range(len(points2_offset)):
			p_2 = poly2[i]

			if gen_icp.isInCircle(p_2, radius1, center1):

				" for every transformed point of A, find it's closest neighbor in B "
				p_1, minDist = gen_icp.findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2_offset[i][2]
					C1 = p_1[2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i],p_1,C2,C1])

		for i in range(len(points1)):
			p_1 = poly1[i]
	
			if gen_icp.isInCircle(p_1, radius2, center2):
		
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2_offset[i_2][2]
					
					C1 = points1[i][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i_2],points1[i],C2,C1])

		vals = []
		sum1 = 0.0
		for pair in match_pairs:
	
			a = pair[0]
			b = pair[1]
			Ca = pair[2]
			Cb = pair[3]
	
			ax = a[0]
			ay = a[1]		
			bx = b[0]
			by = b[1]
	
			c11 = Ca[0][0]
			c12 = Ca[0][1]
			c21 = Ca[1][0]
			c22 = Ca[1][1]
					
			b11 = Cb[0][0]
			b12 = Cb[0][1]
			b21 = Cb[1][0]
			b22 = Cb[1][1]	
		
			val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
						
			vals.append(val)
			sum1 += val


		matchCount = len(vals)

		return sum1, matchCount


	def makeGlobalMedialOverlapConstraint(self, nodeID, globalPath, globalJunctionPoint, globalDeparturePoint):

		" compute the medial axis for each pose "
		
		node1 = self.nodeHash[nodeID]
		posture1 = node1.getStableGPACPosture()
		#hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = True)
		hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = False)
		
		globalPathReverse = deepcopy(globalPath)
		globalPathReverse.reverse()
		

		
		estPose1 = node1.getGlobalGPACPose()		
		poseOrigin = Pose(estPose1)
		
		globalMedial = []
		for p in medial1:
			globalMedial.append(poseOrigin.convertLocalToGlobal(p))
		
		medialSpline1 = SplineFit(globalMedial, smooth=0.1)
		globalSpline = SplineFit(globalPath, smooth=0.1)
		globalSplineReverse = SplineFit(globalPathReverse, smooth=0.1)


		overlapMatch = []
		angleSum1 = 0.0
		angleSum2 = 0.0
		for i in range(0,len(globalMedial)):
			p_1 = globalMedial[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(globalPath, p_1)
			if minDist < 0.5:
				overlapMatch.append((i,j,minDist))

				pathU1 = medialSpline1.findU(p_1)	
				pathU2 = globalSpline.findU(p_2)	
				pathU2_R = globalSplineReverse.findU(p_2)	

				pathVec1 = medialSpline1.getUVector(pathU1)
				pathVec2 = globalSpline.getUVector(pathU2)
				pathVec2_R = globalSplineReverse.getUVector(pathU2_R)

				#path1Mag = sqrt(pathVec1[0]*pathVec1[0] + pathVec1[1]*pathVec1[1])
				#path2Mag = sqrt(pathVec2[0]*pathVec2[0] + pathVec2[1]*pathVec2[1])
				#path2MagR = sqrt(pathVec2_R[0]*pathVec2_R[0] + pathVec2_R[1]*pathVec2_R[1])

				val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
				if val1 > 1.0:
					val1 = 1.0
				elif val1 < -1.0:
					val1 = -1.0
				ang1 = acos(val1)
				
				val2 = pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1]
				if val2 > 1.0:
					val2 = 1.0
				elif val2 < -1.0:
					val2 = -1.0
				ang2 = acos(val2)
				

				#ang1 = acos(pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1])
				#ang2 = acos(pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1])
		
				angleSum1 += ang1
				angleSum2 += ang2
		
		" select global path orientation based on which has the smallest angle between tangent vectors "
		print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
		if angleSum1 > angleSum2:
			orientedGlobalPath = globalPathReverse
		else:
			orientedGlobalPath = globalPath
			
		globalSpline = SplineFit(orientedGlobalPath, smooth=0.1)
		medialSpline1 = SplineFit(medial1, smooth=0.1)


		globalSamples = globalSpline.getUniformSamples(spacing = 0.04)
		medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
		
		globalMedialSamples = []
		for p in medialSamples:
			result = poseOrigin.convertLocalOffsetToGlobal(p)	
			globalMedialSamples.append(result)
		
		
		globalVar = []
		medialVar = []
		
		" compute the local variance of the angle "
		VAR_WIDTH = 40
		for i in range(len(globalSamples)):
			
			lowK = i - VAR_WIDTH/2
			if lowK < 0:
				lowK = 0
				
			highK = i + VAR_WIDTH/2
			if highK >= len(globalSamples):
				highK = len(globalSamples)-1
			
			localSamp = []
			for k in range(lowK, highK+1):
				localSamp.append(globalSamples[k][2])
			
			sum = 0.0
			for val in localSamp:
				sum += val
				
			meanSamp = sum / float(len(localSamp))
			

			sum = 0
			for k in range(len(localSamp)):
				sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
		
			varSamp = sum / float(len(localSamp))
			
			globalVar.append((meanSamp, varSamp))		

		for i in range(len(globalMedialSamples)):
			
			lowK = i - VAR_WIDTH/2
			if lowK < 0:
				lowK = 0
				
			highK = i + VAR_WIDTH/2
			if highK >= len(globalMedialSamples):
				highK = len(globalMedialSamples)-1
			
			localSamp = []
			for k in range(lowK, highK+1):
				localSamp.append(globalMedialSamples[k][2])
			
			sum = 0.0
			for val in localSamp:
				sum += val
				
			meanSamp = sum / float(len(localSamp))
			

			sum = 0
			for k in range(len(localSamp)):
				sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
		
			varSamp = sum / float(len(localSamp))
			
			medialVar.append((meanSamp, varSamp))		


		" now lets find closest points and save their local variances "			
		closestPairs = []
		TERM_DIST = 20
		
		for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
			pG = globalSamples[i]
			minDist = 1e100
			minJ = -1
			for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
				pM = globalMedialSamples[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minJ = j
					
			if minDist < 0.1:
				closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

		for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
			pM = globalMedialSamples[j]
			minDist = 1e100
			minI = -1
			for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
				pG = globalSamples[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minI = i
					
			if minDist < 0.1:
				closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))

		" remove duplicates "
		closestPairs = list(set(closestPairs))
		#closestPairs = list(closestPairs)
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(5,6))
		#s = sorted(student_objects, key=attrgetter('age'))     # sort on secondary key
		#sorted(s, key=attrgetter('grade'), reverse=True)  
		#closestPairs.sort()

		print len(closestPairs), "closest pairs"
		#for val in closestPairs:
		#	print val


		"""
		overlapMatch = []
		for i in range(0,len(globalMedial)):
			p_1 = globalMedial[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(orientedGlobalPath, p_1)
			if minDist < 0.5:
				overlapMatch.append((i,j,minDist))

				pathU1 = medialSpline1.findU(p_1)	
				pathU2 = globalSpline.findU(p_2)	
				
				pathVec1 = medialSpline1.getUVector(pathU1)
				pathVec2 = globalSpline.getUVector(pathU2)

				val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
				if val1 > 1.0:
					val1 = 1.0
				elif val1 < -1.0:
					val1 = -1.0
				ang1 = acos(val1)
		"""
		if len(closestPairs) > 0:
			originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
			originU1 = globalSpline.findU(globalSamples[closestPairs[0][0]])

		else:
			localDeparturePoint = poseOrigin.convertGlobalToLocal(globalDeparturePoint)		
			originU2 = medialSpline1.findU(localDeparturePoint)	
			originU1 = globalSpline.findU(globalJunctionPoint)
		
		u2 = originU2
		u1 = originU1
		angGuess = 0.0
		
		resultPose, lastCost = gen_icp.globalOverlapICP([u1, u2, angGuess], orientedGlobalPath, hull1, medial1,  plotIter = False, n1 = nodeID, n2 = 0)

		print "estimating branch pose", nodeID, "at",  resultPose[0], resultPose[1], resultPose[2]

		return resultPose, lastCost


	def makeMultiJunctionMedialOverlapConstraint(self, i, j, isMove = True, isForward = True, inPlace = False, uRange = 1.5 ):

		#isMove = False
		#inPlace = True

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]


		hull1 = node1.getBareHull()
		hull1.append(hull1[0])

		hull2 = node2.getBareHull()
		hull2.append(hull2[0])

		estPose1 = node1.getGlobalGPACPose()		
		estPose2 = node2.getGlobalGPACPose()

		originProfile = Pose(estPose1)
		diffOffset = originProfile.convertGlobalPoseToLocal(estPose2)
		
		initDiff1 = diffAngle(estPose2[2], estPose1[2])
		print initDiff1, diffOffset[2]


		results = []

		for k in range(len(node1.medialLongPaths)):
			medial1 = node1.medialLongPaths[k]
			for l in range(len(node2.medialLongPaths)):
				medial2 = node2.medialLongPaths[l]
				


				for medialDir in [True, False]:
					for travelDir in [True, False]:
		
						orientedMedial2 = deepcopy(medial2)
						
						if medialDir:
							orientedMedial2.reverse()
		
						medialSpline1 = SplineFit(medial1, smooth=0.1)
						medialSpline2 = SplineFit(orientedMedial2, smooth=0.1)
				
						originU1 = medialSpline1.findU([0.0,0.0])	
						originU2 = medialSpline2.findU([0.0,0.0])	
				
						if inPlace:
							originU1 = 0.6
							originU2 = 0.4
							u2 = originU2
							print "computed u2 =", u2, "from originU2 =", originU2
				
						elif isMove:
							
							#if isForward:
							if travelDir:
								
								if len(node2.frontProbeError) > 0:
					
									frontSum = 0.0
									frontProbeError = node2.frontProbeError
									for n in frontProbeError:
										frontSum += n
									foreAvg = frontSum / len(frontProbeError)
													
									if foreAvg >= 1.4:
										u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
									else:	
										u2 = medialSpline2.getUOfDist(originU2, 0.3, distIter = 0.001)
								else:	
									u2 = medialSpline2.getUOfDist(originU2, 0.3, distIter = 0.001)
									
							else:
								u2 = medialSpline2.getUOfDist(originU2, -0.3, distIter = 0.001)
								#u2 = 0.4
							print "computed u2 =", u2, "from originU2 =", originU2
							 
						else:
							poseOrigin = Pose(estPose1)
							offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
							
							points2 = medialSpline2.getUniformSamples()
							p_1 = medialSpline1.getU(0.5)
							
							points2_offset = []
							for p in points2:
								result = gen_icp.dispOffset(p, offset)		
								points2_offset.append(result)
					
							p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
					
							u2 = medialSpline2.findU(points2[i_2])	
							
							#if u2 > 0.9 or u2 < 0.1:
							#	raise
				
						u1 = originU1
						angGuess = 0.0
					
						" create the ground constraints "
						gndGPAC1Pose = node1.getGndGlobalGPACPose()
						currProfile = Pose(gndGPAC1Pose)
						gndGPAC2Pose = node2.getGndGlobalGPACPose()
						gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
		
						
						#result, hist = gen_icp.overlapICP(estPose1, diffOffset, [u1, u2, angGuess], hull1, hull2, orientedMedial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (i,k), n2 = (j,l), uRange = uRange)
						result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (i,k), n2 = (j,l), uRange = uRange)
				
						transform = matrix([[result[0]], [result[1]], [result[2]]])
						covE =  self.E_overlap
						
						print "making overlap constraint:", result[0], result[1], result[2]
		
						angDiff = abs(diffAngle(diffOffset[2], transform[2,0]))			
						#totalGuesses.append((angDiff, result[2], result[3], result[4]))
		
						points1 = medialSpline1.getUniformSamples()
						points2 = medialSpline2.getUniformSamples()
						p_1 = medialSpline1.getU(0.5)
		
		
						offset = [transform[0,0], transform[1,0], transform[2,0]]
						fitProfile = Pose(offset)
						
						points2_offset = []
						for p in points2:
							result = fitProfile.convertLocalOffsetToGlobal(p)
							points2_offset.append(result)		
						
						overlapMatch = []
						matchCount = 0
						overlapSum = 0.0
						angleSum = 0.0
						for m in range(0,len(points2_offset)):
							p_2 = points2_offset[m]
							p_1, n, minDist = gen_icp.findClosestPointInA(points1, p_2)
			
							if minDist < 0.1:
								overlapMatch.append((m,n,minDist))
								matchCount += 1
								
								overlapSum += minDist
								
								ang1 = p_1[2]
								ang2 = p_2[2]
		
								angleSum += abs(diffAngle(ang1,ang2))
									
						if matchCount > 0:
							angleSum /= matchCount
							overlapSum /= matchCount
				
				
						medialError, matchCount = self.computeMedialError(i, j, offset, minMatchDist = 0.5, tail1=k, tail2=l)
		
				
						results.append((hist[0], angleSum, medialError, matchCount, angDiff, k, l, transform, covE, hist, matchCount, overlapSum, angleSum))
		
		results.sort(reverse=True)
		#results.sort(reverse=False)
		
		print "Multi-Junction Overlap of", i, "and", j
		selectedIndex = 0
		#for result in results:
		for k in range(len(results)):
			print results[k]
			

		for k in range(len(results)):
			if results[k][9][1] < 10 and results[k][9][2] == 0:
				selectedIndex = k
				break

		"""
		totalGuesses = []
		for result in results:

			transform = result[2]

			angDiff = abs(diffAngle(diffOffset[2], transform[2,0]))			
			totalGuesses.append((angDiff, result[2], result[3], result[4]))
		
		totalGuesses.sort()
		"""
		#print "Sorted:"
		#for guess in totalGuesses:
		#	print guess

		transform = results[selectedIndex][7]
		covE = results[selectedIndex][8]
		hist = results[selectedIndex][9]
		
		return transform, covE, hist


	def makeMedialOverlapConstraint(self, i, j, isMove = True, isForward = True, inPlace = False, uRange = 1.5 ):

		#print "recomputing hulls and medial axis"
		" compute the medial axis for each pose "
		
			

		#self.nodeHash[i].computeStaticAlphaBoundary()
		#print "static alpha computed"

		#hull2 = computeBareHull(self.nodeHash[j], sweep = False)
		#hull2 = computeBareHull(self.nodeHash[j], sweep = False, static = True)
		#hull2.append(hull2[0])
		#medial2 = self.nodeHash[j].getMedialAxis(sweep = False)
		#medial2 = self.nodeHash[j].getStaticMedialAxis()

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()
		#hull1, medial1 = computeHullAxis(i, node1, tailCutOff = True)
		#hull2, medial2 = computeHullAxis(j, node2, tailCutOff = True)
		hull1, medial1 = computeHullAxis(i, node1, tailCutOff = False)
		hull2, medial2 = computeHullAxis(j, node2, tailCutOff = False)

		estPose1 = node1.getGlobalGPACPose()		
		estPose2 = node2.getGlobalGPACPose()
		
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medialSpline2 = SplineFit(medial2, smooth=0.1)


		#samples = scipy.arange(0.0,1.0,0.01)

		#originU1 = medialSpline2.findU(node1.rootPose)	
		#originU2 = medialSpline2.findU(node2.rootPose)	
		originU1 = medialSpline1.findU([0.0,0.0])	
		originU2 = medialSpline2.findU([0.0,0.0])	

		#medialSpline1.getUOfDist(originU1, dist)

		if inPlace:
			" FULL LENGTH MEDIAL AXIS "
			originU1 = 0.5
			originU2 = 0.5
			
			" TAIL CUT OFF MEDIAL AXIS "
			#originU1 = 0.6
			#originU2 = 0.4
			
			u2 = originU2
			print "computed u2 =", u2, "from originU2 =", originU2

		elif isMove:
			
			if isForward:
				
				if len(node2.frontProbeError) > 0:
	
					frontSum = 0.0
					frontProbeError = node2.frontProbeError
					for n in frontProbeError:
						frontSum += n
					foreAvg = frontSum / len(frontProbeError)
									
					if foreAvg >= 1.4:
						u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
					else:	
						u2 = medialSpline2.getUOfDist(originU2, 0.3, distIter = 0.001)
				else:	
					u2 = medialSpline2.getUOfDist(originU2, 0.3, distIter = 0.001)
					
			else:
				u2 = medialSpline2.getUOfDist(originU2, -0.3, distIter = 0.001)
				#u2 = 0.4
			print "computed u2 =", u2, "from originU2 =", originU2
			 
		else:
			poseOrigin = Pose(estPose1)
			offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
			
			points2 = medialSpline2.getUniformSamples()
			p_1 = medialSpline1.getU(0.5)
			
			points2_offset = []
			for p in points2:
				result = gen_icp.dispOffset(p, offset)		
				points2_offset.append(result)
	
			p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
	
			u2 = medialSpline2.findU(points2[i_2])	
			
			if u2 > 0.9 or u2 < 0.1:
				raise

		#print "u2 =", u2


		#motionT, motionCov = self.makeMotionConstraint(i,j)
		#travelDelta = motionT[0,0]
		
		u1 = originU1
		#u1 = 0.5
		#u2 = 0.6
		angGuess = 0.0
	
		" create the ground constraints "
		gndGPAC1Pose = node1.getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)
		gndGPAC2Pose = node2.getGndGlobalGPACPose()
		gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
		
		#supportLine = self.getTopology(nodes = range(self.numNodes-2))

		#result = gen_icp.overlapICP(estPose1, [u1, u2, angGuess], hull1, hull2, medial1, medial2, node1.rootPose, node2.rootPose, plotIter = True, n1 = i, n2 = j)
		result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)
		#result, hist = gen_icp.overlapICP2(estPose1, [u1, u2, angGuess], hull1, hull2, medial1, medial2, supportLine, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = True, n1 = i, n2 = j)

		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE =  self.E_overlap
		
		print "making overlap constraint:", result[0], result[1], result[2]

		return transform, covE, hist
	
	def makeInPlaceConstraint(self, nodeID1, nodeID2):

		" compute the new posture "
		" 1. compute the GPAC pose "
		" 2. perform correction of GPAC pose "
		" 3. compute location of rootPose from corrected GPAC pose "
		
		" node1 is the front poke node "
		" nodes is the back poke node "

		node1 = self.nodeHash[nodeID1]
		node2 = self.nodeHash[nodeID2]
		
		#print "rootPose1", node1.rootPose
		#print "rootPose2", node2.rootPose
		#print "estPose1", node1.estPose
		#print "estPose2", node2.estPose

		#originPosture = node1.localPosture
		originPosture = node1.correctedPosture
		originCurve = StableCurve(originPosture)
		originForePose, originBackPose = originCurve.getPoses()
		
		originForeProfile = Pose(originForePose)
		originBackProfile = Pose(originBackPose)
		
		originForePosture = []
		originBackPosture = []
		for i in range(len(originPosture)):
			originForePosture.append(originForeProfile.convertGlobalPoseToLocal(originPosture[i]))
		for i in range(len(originPosture)):
			originBackPosture.append(originBackProfile.convertGlobalPoseToLocal(originPosture[i]))

		#newPosture = node2.localPosture
		newPosture = node2.correctedPosture
		#newPosture = []
		#for j in range(self.probe.numSegs-1):
		#	newPosture.append(self.probe.getJointPose([0.0,0.0,0.0], node1.rootNode, j))


		localPosture = newPosture
		localCurve = StableCurve(localPosture)

		localForePose, localBackPose = localCurve.getPoses()
		 
		localForeProfile = Pose(localForePose)
		localBackProfile = Pose(localBackPose)
		
		localForePosture = []
		localBackPosture = []
		for i in range(len(localPosture)):
			localForePosture.append(localForeProfile.convertGlobalPoseToLocal(localPosture[i]))
		for i in range(len(localPosture)):
			localBackPosture.append(localBackProfile.convertGlobalPoseToLocal(localPosture[i]))
	
		foreCost = 0.0
		for j in range(0,20):
			p1 = originForePosture[j]
			p2 = localForePosture[j]
			
			foreCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)		   
		
		backCost = 0.0
		for j in range(20,39):
			p1 = originBackPosture[j]
			p2 = localBackPosture[j]
			
			backCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)		   

		angle = 0.0

		correctedGPACPose1 = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		correctedGPACPose2 = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		correctedProfile1 = Pose(correctedGPACPose1)
		correctedProfile2 = Pose(correctedGPACPose2)

		localRootOffset3_1 = correctedProfile1.convertGlobalPoseToLocal([0.0,0.0,0.0])
		localRootOffset3_2 = correctedProfile2.convertGlobalPoseToLocal([0.0,0.0,0.0])

		offset1 = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3_1)
		offset2 = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3_2)
		
		" except the case where one node has not updated "
		if self.nodeHash[nodeID1].occMap.mapImage == 0 or self.nodeHash[nodeID2].occMap.mapImage == 0:
			cost1 = 1.0
			cost2 = 0.0
		else:
			cost1, matchCount1 = self.computeMedialError(nodeID1, nodeID2, offset1)
			cost2, matchCount2 = self.computeMedialError(nodeID1, nodeID2, offset2)
					
		#if foreCost > backCost:
		if cost1 < cost2:
			#plotPosture(originBackPosture, localBackPosture)
			correctedGPACPose = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		else:
			#plotPosture(originForePosture, localForePosture)					
			correctedGPACPose = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])

		correctedProfile = Pose(correctedGPACPose)

		" offset from 0,0 to newGPACPose "
		" rootPose from neg. offset from correctedGPACPose "
		localRootOffset3 = correctedProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		#localRootOffset3 = correctedProfile.convertGlobalPoseToLocal(node2.rootPose)
		
		
		#if foreCost > backCost:
		#if foreCost > backCost:
		if cost1 < cost2:
			offset = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3)
		else:
			offset = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3)
		
		transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
		covE = self.E_inplace
		
		return transform, covE

		

	def addPathConstraints(self, pathNodes, targetNodeID):

		print "addPathConstraints:"
		print "pathNodes:", pathNodes
		print "targetNodeID:", targetNodeID

		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		#cart_distances = []
	#	targetNode = self.nodeHash[targetNodeID]
		#targetPose = targetNode.getGlobalGPACPose()
		#for i in range(len(self.nodeSets[0])):
		#	
		#	candNode = self.nodeHash[self.nodeSets[0][i]]
		#	candPose = candNode.getGlobalGPACPose()
		#	
		#	dist = sqrt((candPose[0]-targetPose[0])**2 + (candPose[1]-targetPose[1])**2)
		#	cart_distances.append(dist)

		#print "non-candidate node distances:"
		#for i in range(len(self.nodeSets[0])):
		#	print self.nodeSets[0][i], cart_distances[i]

		print "candidate node distances", targetNodeID, ":"

		" compute cartesian distances "
		cart_distances = []
		angle_diffs = []
		targetNode = self.nodeHash[targetNodeID]
		targetPose = targetNode.getGlobalGPACPose()
		for i in range(len(pathNodes)):
			
			candNode = self.nodeHash[pathNodes[i]]
			candPose = candNode.getGlobalGPACPose()
			
			dist = sqrt((candPose[0]-targetPose[0])**2 + (candPose[1]-targetPose[1])**2)
			angDelta = diffAngle(candPose[2],targetPose[2])
			print pathNodes[i], dist, angDelta, targetPose, candPose
			cart_distances.append(dist)

			angle_diffs.append(angDelta)
		#for i in range(len(pathNodes)):
		#	print pathNodes[i], cart_distances[i]

		constraintNodes = deepcopy(pathNodes)

		" remove the node that is the target from consideration "
		delIndex = pathNodes.index(targetNodeID)
		constraintNodes.pop(delIndex)
		cart_distances.pop(delIndex)
		angle_diffs.pop(delIndex)

		" remove the paired inplace node if it exists "
		if targetNodeID % 2 == 0:
			pairNodeID = targetNodeID + 1
		else:
			pairNodeID = targetNodeID - 1

		try:
			delIndex = pathNodes.index(pairNodeID)
			constraintNodes.pop(delIndex)
			cart_distances.pop(delIndex)
			angle_diffs.pop(delIndex)
		except:
			pass

		" remove the neighbor node that has no constraint if it exists "
		if targetNodeID % 2 == 1:
			pairNodeID = targetNodeID + 1
		else:
			pairNodeID = targetNodeID - 1

		try:
			delIndex = pathNodes.index(pairNodeID)
			constraintNodes.pop(delIndex)
			cart_distances.pop(delIndex)
			angle_diffs.pop(delIndex)
		except:
			pass


		
		
		
		#cart_distances.remove(pairNodeID)



		#constraintNodes = pathNodes[:-2]
		#constraintNodes = pNodes
		constraintPairs = []
		
		PATH_MATCH_DIST = 1.0
		#PATH_MATCH_DIST = 0.5
		

		for i in range(len(constraintNodes)):
			nodeID = constraintNodes[i]
			dist = cart_distances[i]
			angDelta = angle_diffs[i]

			if dist < PATH_MATCH_DIST:
				state1 = self.nodeHash[targetNodeID].getDeparture()
				state2 = self.nodeHash[nodeID].getDeparture()
				
				
				"FIXME:  Will always return true, this state is no longer recorded "
				if state1 == state2:
					print targetNodeID, "and", nodeID, "have same departure state", state1, state2
					
					" check if a constraint already exists between these two nodes "
					oldConstraints = self.getEdges(targetNodeID, nodeID)
					isExist = False
					for constraint in oldConstraints:
						isExist = True

					if not isExist:
						constraintPairs.append([targetNodeID,nodeID,angDelta])
					else:
						print "constraint already exists between", targetNodeID, "and", nodeID
						
				else:
					print targetNodeID, "and", nodeID, "have different departure state", state1, state2 
			else:	
				print targetNodeID, "and", nodeID, "are too far away", dist 

		print "received", len(constraintPairs), "possible pairs"

		" now, take only the 5 oldest nodes to constrain with "
		
		"FIXED:  Sorting on targetNodeID instead of nodeID "
		#constraintPairs.sort()
		constraintPairs = sorted(constraintPairs, key = itemgetter(1))
		if len(constraintPairs) > 3:
			constraintPairs = constraintPairs[:3]
			print "reducing to only", len(constraintPairs)

		PAIR_ANG_DELTA = 0.3

		" make hypothesized constraints out of the pairs "
		constraintResults = []
		for i in range(len(constraintPairs)):
			p = constraintPairs[i]
			n1 = p[0]
			n2 = p[1]
			angDiff = p[2]

			try:

				if self.nodeHash[n1].getIsFeatureless():	
					
					if self.nodeHash[n2].getIsFeatureless():	
						if self.nodeHash[n2].getNumLeafs() > 2 or self.nodeHash[n1].getNumLeafs() > 2:
							transform, covE, hist = self.makeMultiJunctionMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)
						else:
							transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)				
	
						if hist[1] > 0 or hist[2] > 0:
							" disregard if the fit is poor "
							pass
						else:
							
							" compare this new orientation to the old, is it too different?  If it is, throw it out."
							angDelta = diffAngle(transform[2,0], angDiff)
							
							if fabs(angDelta) < PAIR_ANG_DELTA:
								#constraintResults.append([n1, n2, transform, deepcopy(self.E_featureless), angDiff])
								constraintResults.append([n1, n2, transform, deepcopy(self.E_featureless), angDelta])
							else:
								print "rejecting", (n1,n2), "constraint because angDelta is", angDelta, "with transform angle", transform[2,0]

				else:
					if not self.nodeHash[n2].getIsFeatureless():

						if self.nodeHash[n2].getNumLeafs() > 2 or self.nodeHash[n1].getNumLeafs() > 2:
							transform, covE, hist = self.makeMultiJunctionMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)
						else:
							transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)				
	
						if hist[1] > 0 or hist[2] > 0:
							" disregard if the fit is poor "
							pass
						else:
	
							" compare this new orientation to the old, is it too different?  If it is, throw it out."
							angDelta = diffAngle(transform[2,0], angDiff)
							
							if fabs(angDelta) < PAIR_ANG_DELTA:								
								constraintResults.append([n1, n2, transform, deepcopy(self.E_junction), angDelta])
							else:
								print "rejecting", (n1,n2), "constraint because angDelta is", angDelta, "with transform angle", transform[2,0]
							
						" TODO: "
						" 1) find the overlapping region curves "
						" 2) determine if the curves are straight and featureless "
						" 3) if so, add high variance along the tangent "
				
			except:
				pass
		
		print "received", len(constraintResults), "constraint results"
		for result in constraintResults:
			print result
		#if len(constraintResults) > 0:
		#	print constraintResults[0]	


		totalHypotheses = constraintResults
		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) == 0:
			print "no hypotheses!  returning "
			newConstraints = []
			#return
		
		elif len(totalHypotheses) == 1:
			newConstraints = totalHypotheses

		else: 
					
			" set all proscribed hypothesis sets to maximum error "	
			maxError = 0.0
			for result in results:
				if result[0] != None and result[0] > maxError:
					maxError = result[0]
	
			maxError = maxError*2
	
			for result in results:
				if result[0] == None:
					result[0] = maxError
		
			" create the consistency matrix "
			A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
			A.resize(len(totalHypotheses), len(totalHypotheses))
				
			" populate consistency matrix "					
			for result in results:
				i = result[1]
				j = result[2]
				A[i,j] = maxError - result[0]
				A[j,i] = maxError - result[0]
	
			" do graph clustering of consistency matrix "
			w = []
			e = []
			for i in range(100):
				e, lmbda = scsgp.dominantEigenvectors(A)
				w = scsgp.getIndicatorVector(e[0])
				w2 = scsgp.getIndicatorVector(e[1])
				if len(e) <= 1:
					break
	
				" threshold test l1 / l2"
				ratio = lmbda[0][0,0] / lmbda[1][0,0]
				if ratio >= 2.0:
					break
	
			selected = []	
			for i in range(len(totalHypotheses)):
				if w[i,0] >= 1.0:
					selected.append(totalHypotheses[i])
	
			selected2 = []	
			for i in range(len(totalHypotheses)):
				if w2[i,0] >= 1.0:
					selected2.append(totalHypotheses[i])
	
			newConstraints = selected2

			newConstraints = totalHypotheses
		print "adding", len(newConstraints), "medial overlap constraints on path"

		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		#self.mergePriorityConstraints()
		
		#return selected, selected2

		return
		



	def doToro(self, constraints, fileName = "probe"):
		
		
		
		" vertices "
		v_list = []

		for k, v in self.nodeHash.items():
			node1 = v
			estPose1 = node1.getGlobalGPACPose()
			v_list.append([k, [estPose1[0], estPose1[1], estPose1[2]]])

		if len(constraints) == 0:
			return v_list, []
		
		e_list = []
		for const in constraints:
				transform, covE = const[2], const[3]
				node1 = const[0]
				node2 = const[1]
				offset = [transform[0,0],transform[1,0],transform[2,0]]

				invMat = scipy.linalg.inv(covE)				
				prec = [invMat[0,0], invMat[0,1], invMat[1,1], invMat[2,2], invMat[0,2], invMat[1,2]]
				
				# inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr 
				e_list.append([node1,node2,offset,prec])			
		
		#print "Running TORO with", len(e_list), "constraints"
		" add hypotheses to the TORO file "
		
		#print "e_list:", e_list
		#print "v_list:", v_list

		toro.writeConstrainedToroGraph(".", fileName + ".graph", v_list, e_list)
		toro.executeToro("./" + fileName + ".graph")
	
		finalFileName = fileName + "-treeopt-final.graph"
		v_list2, e_list2 = toro.readToroGraph("./" + finalFileName)		
	
		final_constraints = []
		
		#print "TORO RESULTS:"
		for edge in e_list2:
			node1 = edge[0]
			node2 = edge[1]
			offset = edge[2]
			prec = edge[3]
			#print node1, node2, offset, prec

			invMat = matrix([[0.,0.,0.],
							[0.,0.,0.],
							[0.,0.,0.]])
			invMat[0,0] = prec[0]
			invMat[0,1] = invMat[1,0] = prec[1]
			invMat[1,1] = prec[2]
			invMat[2,2] = prec[3]
			invMat[0,2] = invMat[2,0] = prec[4]
			invMat[1,2] = invMat[2,1] = prec[5]
			
			transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
			covE = scipy.linalg.inv(invMat)				
			
			final_constraints.append([node1, node2, transform, covE])			

		return v_list2, final_constraints

	def plotEnv(self):
		
		walls = self.probe.getWalls()
	
		for wall in walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = copy(wall[i])
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g')
			


	def drawConstraints(self, id = []):
		
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())

		pylab.clf()
		for i in range(self.numNodes):

			hull = computeBareHull(self.nodeHash[i], sweep = False)
			hull.append(hull[0])

			node1 = self.nodeHash[i]
			currPose = currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	

			hull_trans = []
			for p in hull:
				hull_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	

		for k, v in self.edgeHash.items():	
			node1 = self.nodeHash[k[0]]
			node2 = self.nodeHash[k[1]]
			pose1 = node1.getGlobalGPACPose()
			pose2 = node2.getGlobalGPACPose()

			xP = [pose1[0], pose2[0]]
			yP = [pose1[1], pose2[1]]

			#age = self.edgeAgeHash[k]
			age = 0

			fval = float(age) / 20.0
			if fval > 0.4:
				fval = 0.4
			color = (fval, fval, fval)
			pylab.plot(xP, yP, color=color, linewidth=1)


		xP = []
		yP = []
		for pose in poses:
			xP.append(pose[0])
			yP.append(pose[1])		
		pylab.scatter(xP,yP, color='k', linewidth=1)

		#pylab.title("Corner Constraint %d ---> %d" % (id1, id2))
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		pylab.title("%d Poses" % self.numNodes)
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)

			
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)

	def drawTrimmedPaths(self, trimmedPaths):
		
		pylab.clf()
		
		for k in range(len(trimmedPaths)):
			path = trimmedPaths[k]
			print "path has", len(path), "points"
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color = self.colors[k])

		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		pylab.title("Trimmed Paths, numNodes = %d" % self.numNodes)
		pylab.savefig("trimmedPath_%04u.png" % self.trimCount)
		self.trimCount += 1



	def drawPathAndHull(self):
		
		
		pylab.clf()
		for k in range(len(self.nodeSets)):
			xP = []
			yP = []
			
			for p in self.paths[k]:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=self.colors[k], linewidth=4)


			for nodeID in self.nodeSets[k]:
				xP = []
				yP = []

				estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()		
		
				if self.nodeHash[nodeID].isBowtie:			
					hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
				else:
					hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
		
				" set the origin of pose 1 "
				poseOrigin = Pose(estPose1)
		
				points = []
				for p in hull1:
					p1 = poseOrigin.convertLocalToGlobal(p)
					points.append(p1)
				
				for p in points:
					xP.append(p[0])
					yP.append(p[1])
					
				pylab.plot(xP,yP, color=self.colors[k])

			xP = []
			yP = []
			for p in self.hulls[k]:
				xP.append(p[0])
				yP.append(p[1])
				
			pylab.plot(xP,yP, '--', color=self.colors[k], linewidth=4)
			
		print "pathAndHull:", self.pathDrawCount, self.nodeSets.values()	
		pylab.title("%s" % repr(self.nodeSets.values()))
		pylab.savefig("pathAndHull_%04u.png" % self.pathDrawCount)

		self.pathDrawCount += 1
			
		
				
