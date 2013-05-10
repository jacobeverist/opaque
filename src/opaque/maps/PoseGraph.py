import sys

import hashlib
import scipy.linalg
from scipy.sparse.linalg import eigsh, eigs
import random

import numpy
from numpy import array
import math
from LocalNode import getLongestPath
from StableCurve import StableCurve
from SplineFit import SplineFit
from ParticleFilter import ParticleFilter
from Pose import Pose
from Paths import Paths
import gen_icp
from functions import *
import toro
import scsgp
import bayes
from operator import itemgetter
import cProfile
import time
import traceback

import pylab
#import Image

#from subprocess import Popen, PIPE
#from medialaxis import computeMedialAxis
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


globalFunc = 0
globalArg = 0

def computeBareHull(node1, sweep = False, static = False):
	
	if static:
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		
		m = hashlib.md5()
		m.update(repr(a_data))
		#print "PoseGraph: computeBareHull():", int(m.digest().encode('hex'),16)   
				
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

		m = hashlib.md5()
		m.update(repr(a_data))
		#print "PoseGraph: computeBareHull():", int(m.digest().encode('hex'),16)   

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

def printStack():

	#traceback.print_stack()
	flist = traceback.format_stack()
	flist = flist[:-1]
	
	printStr = ""
	for line in flist:
		printStr += line
		
	print printStr

class PoseGraph:

	def __init__(self, probe, contacts):
		
		self.probe = probe
		self.contacts = contacts
		
		self.walls = self.probe.getWalls()
		
		self.initPose = self.probe.getActualJointPose(19)
		self.nodeHash = {}
		self.edgeHash = {}
		
		self.numNodes = 0
		self.currNode = 0
		
		self.edgePriorityHash = {}
		self.cornerBins = []

		self.paths = Paths(self.nodeHash)
		
		self.pathPlotCount = 0
		self.overlapPlotCount = 0
		self.pathDrawCount = 0
		self.candidateCount = 0
		self.statePlotCount = 0

		#self.numPaths = 1
		#self.paths = {0 : []}
		#self.hulls = {0 : []}
		#self.nodeSets = {0 : []}
		
		
		" FIXME:  self.currPath is used but never modified  "
		self.currPath = 0
		#self.pathParents = [[None,None,None,None]]
		#self.pathTermsVisited = {0: False}

		self.trimCount = 0
		self.spliceCount = 0
		self.orderCount = 0		
		
		self.colors = []
		self.colors.append([215, 48, 39])
		self.colors.append([252, 141, 89])
		self.colors.append([254, 224, 144])
		self.colors.append([224, 243, 248])
		self.colors.append([145, 191, 219])
		self.colors.append([69, 117, 180])

		for color in self.colors:
			color[0] = float(color[0])/256.0
			color[1] = float(color[1])/256.0
			color[2] = float(color[2])/256.0
		
		for i in range(1000):
			self.colors.append((random.random(),random.random(),random.random()))

		#self.colors = ['b','r','g','k','y']

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
		self.E_motion = matrix([[ 0.06, 0.0, 0.0],
							[0.0,  5.0, 0.0],
							[0.0, 0.0,  0.1 ]])
		
		self.E_sensor = matrix([[0.1,0.0,0.0],
							[0.0,0.05,0.0],
							[0.0,0.0,0.02]])
		
		self.E_relaxed = matrix([[ 1.0, 0.0, 0.0],
							[ 0.0, 0.1, 0.0],
							[0.0, 0.0, pi/8.0]])


		self.allPathConstraints = []
		self.hypPlotCount = 0


		self.nodePoses = []

	def doNothing1(self):
		pass
	
	def doNothing2(self):
		
		pass
	
	def doNothing3(self):
		pass

	def saveState(self):

		self.nodePoses = []
		for k in range(self.numNodes):
			self.nodePoses.append(self.nodeHash[k].getGlobalGPACPose())
					
		saveFile = ""

		saveFile += "self.walls = " + repr(self.walls) + "\n"
		saveFile += "self.initPose = " + repr(self.initPose) + "\n"
		#saveFile += "self.nodeHash = " + repr(self.nodeHash) + "\n"
		saveFile += "self.edgePriorityHash = " + repr(self.edgePriorityHash) + "\n"
		saveFile += "self.edgeHash = " + repr(self.edgeHash) + "\n"
		saveFile += "self.numNodes = " + repr(self.numNodes) + "\n"
		saveFile += "self.cornerBins = " + repr(self.cornerBins) + "\n"
		#saveFile += "self.numPaths = " + repr(self.numPaths) + "\n"
		#saveFile += "self.nodeSets = " + repr(self.nodeSets) + "\n"
		saveFile += "self.currPath = " + repr(self.currPath) + "\n"
		#saveFile += "self.pathParents = " + repr(self.pathParents) + "\n"
		saveFile += "self.allPathConstraints = " + repr(self.allPathConstraints) + "\n"

		saveFile += "self.nodePoses = " + repr(self.nodePoses) + "\n"



		f = open("stateSave_%04u.txt" % (self.numNodes-1), 'w')
		f.write(saveFile)
		f.close()		

		" SAVE STATE "
		self.paths.saveState(self.numNodes-1)
		#self.paths = Paths(self.nodeHash)
		self.currPath = 0
		
		
		#self.currNode = 0
					





		
				
		#self.pathPlotCount = 0
		#self.overlapPlotCount = 0
		#self.pathDrawCount = 0
		
		#self.paths = {0 : []}
		#self.hulls = {0 : []}
		#self.pathTermsVisited = {0: False}

		#self.trimCount = 0
		#self.spliceCount = 0
		#self.orderCount = 0		
	
		
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
		
		self.paths.restoreState(dirName, numNodes)
		



	def pairLastTwo(self):
	
		node1 = self.nodeHash[self.numNodes-1]
		node2 = self.nodeHash[self.numNodes-2]

		node1.setPartnerNodeID(self.numNodes-2)
		node2.setPartnerNodeID(self.numNodes-1)



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

			#if k[0] == nodeID1 and k[1] == nodeID2 or k[0] == nodeID2 and k[1] == nodeID1:
			if k[0] == nodeID1 and k[1] == nodeID2:
				
				for edge in v:
					resultEdges.append(edge)
				
		return resultEdges

	def relaxPriorityEdge(self, nodeID1, nodeID2, priorityLevel = -1):

		# self.addPriorityEdge([nodeID1,nodeID2,transform1,covE1], INPLACE_PRIORITY)
		# self.edgePriorityHash[(edge[0],edge[1])].append([edge[2], edge[3], priority])
		
		for k, v in self.edgePriorityHash.items():

			if k[0] == nodeID1 and k[1] == nodeID2 or k[0] == nodeID2 and k[1] == nodeID1:

				for const in v:
					if priorityLevel == -1:
						const[1] = self.E_relaxed
					elif const[2] == priorityLevel:
						const[1] = self.E_relaxed
				
		return

	def delPriorityEdge(self, nodeID1, nodeID2, priorityLevel = -1):
		
		for k, v in self.edgePriorityHash.items():

			if k[0] == nodeID1 and k[1] == nodeID2 or k[0] == nodeID2 and k[1] == nodeID1:
				
				if priorityLevel == -1:
					self.edgePriorityHash[k] = []

				else:							
					keepList = []
					for const in v:
						if const[2] != priorityLevel:
							keepList.append(const)
					
					self.edgePriorityHash[k] = keepList
		return
				
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
					
					" take the target priority constraints only "
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

	" functions added to change the code size but not impact function "
	" testing for variation in function "


	def addNode(self, newNode):
		
		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1

	def correctNode2(self, nodeID):

		newNode = self.nodeHash[nodeID]		

		return self.integrateNode(newNode, nodeID)


	def loadNewNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		return self.integrateNode(newNode, nodeID)


	def insertPose(self, foreNode, backNode, initLocation = [0.0,0.0,0.0]):
		
		" CHECK FOR A BRANCHING EVENT "

			

		for abc in range(0,2):
			
			if abc == 0:
				nodeID1 = self.numNodes-2
				nodeID2 = self.numNodes-1
			else:
				self.currNode = foreNode
				self.currNode.setEstPose(initLocation)
				nodeID = self.numNodes
				self.nodeHash[nodeID] = self.currNode
				self.currNode.nodeID = nodeID
				self.numNodes += 1
				
				#self.insertNode(nodeFore, nodeID, initLocation)
				
				self.currNode = backNode
				self.currNode.setEstPose(initLocation)
				nodeID = self.numNodes
				self.nodeHash[nodeID] = self.currNode
				self.currNode.nodeID = nodeID
				self.numNodes += 1
		
				foreNode.setPartnerNodeID(backNode.nodeID)
				backNode.setPartnerNodeID(foreNode.nodeID)
		
				" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
				direction = foreNode.travelDir
		
				nodeID1 = foreNode.nodeID
				nodeID2 = backNode.nodeID
		
				" ensure the axes are computed before this check "
				computeHullAxis(nodeID1, foreNode, tailCutOff = False)
				computeHullAxis(nodeID2, backNode, tailCutOff = False)
		
				" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
				self.paths.generatePaths()
		
				" PERFORM INPLACE CONSTRAINT BETWEEN PAIR "
				supportLine = self.paths.paths[self.currPath]
				
				transform, covE = self.makeInPlaceConstraint(nodeID1, nodeID2)
				transform1 = transform
				offset1 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE1 = covE
				
				if len(supportLine) == 0:
					resultSum1 = 1e100
				else:
					resultSum1 = self.checkSupport(nodeID1, nodeID2, offset1, supportLine)
				
				if self.nodeHash[nodeID1].getNumLeafs() > 2 or self.nodeHash[nodeID2].getNumLeafs() > 2:
					transform, covE, overHist = self.makeMultiJunctionMedialOverlapConstraint(nodeID1, nodeID2, isMove = False, inPlace = False, isForward = direction )
				else:
					transform, covE, overHist = self.makeMedialOverlapConstraint(nodeID1, nodeID2, isMove = False, inPlace = True, isForward = direction)
				
				
				transform3 = transform
				offset3 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE3 = covE
				if len(supportLine) == 0:
					resultSum3 = 1e100
				else:
					resultSum3 = self.checkSupport(nodeID1, nodeID2, offset3, supportLine)
		
				print "INPLACE sums:", resultSum1, resultSum3
		
				if len(supportLine) > 0 and resultSum1 < resultSum3 and resultSum3 - resultSum1 > 100.0:
					self.addPriorityEdge([nodeID1,nodeID2,transform1,covE1], INPLACE_PRIORITY)
				else:
					self.addPriorityEdge([nodeID1,nodeID2,transform3,covE3], INPLACE_PRIORITY)
				
				
				nodeID1 = foreNode.nodeID
				nodeID2 = backNode.nodeID

			
			" if these nodes are already path-classified, return"
			isContained1 = False
			isContained2 = False
			
			pathIDs = self.paths.getPathIDs()
			for k in pathIDs:
				if self.paths.getNodes(k).count(nodeID1) > 0:
					isContained1 = True
				if self.paths.getNodes(k).count(nodeID2) > 0:
					isContained2 = True
					
					
			if isContained1 or isContained2:
				return
			
			" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
			self.paths.generatePaths()

			self.drawPathAndHull()
		
			" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
			if len(self.paths.paths[0]) == 0:
				" first nodes in path 0 "
				self.paths.addNode(nodeID1,0)
				self.paths.addNode(nodeID2,0)
	
				" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
			else:
				
				" departure events for node 1 "
				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
				contig1 = []
	
				" departure events for node 2 "
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []
				contig2 = []
	
				
				" raw medial axis of union of path nodes "
				paths = {}
				pathIDs = self.paths.getPathIDs()
				for k in pathIDs:
					path = self.paths.paths[k]
					if len(path) > 0:
						paths[k] = path
	
				print "paths:", len(paths)
				for k in pathIDs:
					print k, len(paths[k])
	
	
				" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
				self.trimmedPaths = self.paths.trimPaths(paths)
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)
	
				
				" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
					
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
				
				print nodeID1, "orderedPathIDs1:", orderedPathIDs1
				print nodeID2, "orderedPathIDs2:", orderedPathIDs2
				
				
				" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
					contig1.append((contigFrac, overlapSum))
				
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					contig2.append((contigFrac, overlapSum))
					
				print "node departures", nodeID1, ":", departures1
				print "node  interiors", nodeID1, ":", interiors1
				print "node departures", nodeID2, ":", departures2
				print "node  interiors", nodeID2, ":", interiors2
	
				print "node contiguity", nodeID1, ":", contig1
				print "node contiguity", nodeID2, ":", contig2
					
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
	
				if contig1[0][0] < 0.4 or contig1[-1][0] < 0.4:
					adjustPose1 = True
				else:
					adjustPose1 = False

				if contig2[0][0] < 0.4 or contig2[-1][0] < 0.4:
					adjustPose2 = True
				else:
					adjustPose2 = False
						
				#print "pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
				#print "pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2
				print "insert pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1, contig1[0][0], contig1[-1][0]
				print "insert pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2, contig2[0][0], contig2[-1][0]
					
				" check for the cases that internal departure may have a departure point discrepancy, "
				" then we should perform a pose adjustment and recompute departures "
				print "checking discrepancy in departure:", nodeID1, nodeID2, foreTerm1, discForeTerm1, backTerm1, discBackTerm1, foreTerm2, discForeTerm2, backTerm2, discBackTerm2
	
				if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2 or adjustPose1 or adjustPose2:
	
					print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2
	
					#if discForeTerm1 or discBackTerm1:
					if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or adjustPose1:
	
	
						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]
						print "adjusting pose", nodeID1, "from", self.nodeHash[nodeID1].getGlobalGPACPose()
	
						splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs1)
	
						totalGuesses = []
						for path in splicedPaths1:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
							angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID1].setGPACPose(guessPose)
						print "set to pose", nodeID1, "to", guessPose
	
					if foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2 or adjustPose2:
	
						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]
	
						print "adjusting pose", nodeID2, "from", self.nodeHash[nodeID2].getGlobalGPACPose()
	
						splicedPaths2 = self.paths.splicePathIDs(orderedPathIDs2)
	
						totalGuesses = []
						for path in splicedPaths2:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
							angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
							
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID2].setGPACPose(guessPose)

						print "set to pose", nodeID2, "to", guessPose
	
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
					contig1 = []
		
					departures2 = []
					interiors2 = []
					depPoints2 = []
					distances2 = []
					depAngles2 = []
					contig2 = []
					
					" the overlapping paths are computed from the initial guess of position "
					orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
					orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
					
					print nodeID1, "orderedPathIDs1:", orderedPathIDs1
					print nodeID2, "orderedPathIDs2:", orderedPathIDs2
					
					" now compute whether there are departure points after we have guessed a better position in synch with the paths "
					for pathID in orderedPathIDs1:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
						departures1.append([isExist1,isExist2])
						interiors1.append([isInterior1, isInterior2])
						depPoints1.append([departurePoint1, departurePoint2])
						distances1.append([discDist1, discDist2])
						depAngles1.append([depAngle1, depAngle2])
						contig1.append((contigFrac, overlapSum))
					
					for pathID in orderedPathIDs2:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
						departures2.append([isExist1,isExist2])
						interiors2.append([isInterior1, isInterior2])
						depPoints2.append([departurePoint1, departurePoint2])
						distances2.append([discDist1, discDist2])
						depAngles2.append([depAngle1, depAngle2])
						contig2.append((contigFrac, overlapSum))
					
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
	
	
	
				#newPaths = []
	
				frontExist1 = departures1[0][0]
				frontInterior1 = interiors1[0][0]
		
				frontExist2 = departures2[0][0]
				frontInterior2 = interiors2[0][0]
				
				depAngle1 = depAngles1[0][0]
				depAngle2 = depAngles2[0][0]
				
				depPoint1 = depPoints1[0][0]
				depPoint2 = depPoints2[0][0]
				
				parentPathID1 = orderedPathIDs1[0]
				parentPathID2 = orderedPathIDs2[0]
	

				isFront1 = self.nodeHash[nodeID1].faceDir
				isFront2 = self.nodeHash[nodeID2].faceDir
				if isFront1 and not isFront2:
					dirFlag = 0
				elif not isFront1 and isFront2:
					dirFlag = 1
				else:
					print isFront1, isFront2
					raise
				
				isBranch, pathBranchIDs, isNew = self.paths.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)
	
				print "determineBranch:"
				print isBranch
				print pathBranchIDs
	
				isNew1 = isNew[0]
				isNew2 = isNew[1]

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
					
				#for pathID in pathBranchIDs:
				#	if pathID != -1 and newPaths.count(pathID) == 0:
				#		newPaths.append(pathID)
	
	
	
	
				backExist1 = departures1[-1][1]
				backInterior1 = interiors1[-1][1]
		
				backExist2 = departures2[-1][1]
				backInterior2 = interiors2[-1][1]
				
				depAngle1 = depAngles1[-1][1]
				depAngle2 = depAngles2[-1][1]
				
				depPoint1 = depPoints1[-1][1]
				depPoint2 = depPoints2[-1][1]
	
				parentPathID1 = orderedPathIDs1[-1]
				parentPathID2 = orderedPathIDs2[-1]
				#print "self.pathParents:", self.pathParents
				#print "self.nodeSets:", self.nodeSets
				#print "newPaths:", newPaths
				#print "self.numPaths:", self.numPaths
	

				isFront1 = self.nodeHash[nodeID1].faceDir
				isFront2 = self.nodeHash[nodeID2].faceDir
				if isFront1 and not isFront2:
					dirFlag = 1
				elif not isFront1 and isFront2:
					dirFlag = 0
				else:
					print isFront1, isFront2
					raise
				
				isBranch, pathBranchIDs, isNew = self.paths.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)
				print "determineBranch:"
				print isBranch
				print pathBranchIDs


				isNew1 = isNew1 or isNew[0]
				isNew2 = isNew2 or isNew[1]
	
				if isBranch[0]:
					orderedPathIDs1.append(pathBranchIDs[0])
					departures1.append([False, False])
					interiors1.append([False, False])
					depPoints1.append([None, None])
	
				if isBranch[1]:
					orderedPathIDs2.append(pathBranchIDs[1])
					departures2.append([False, False])
					interiors2.append([False, False])
					depPoints2.append([None, None])
	
				#for pathID in pathBranchIDs:
				#	if pathID != -1 and newPaths.count(pathID) == 0:
				#		newPaths.append(pathID)
				#print "self.pathParents:", self.pathParents
				#print "self.nodeSets:", self.nodeSets
				#print "newPaths:", newPaths
				#print "self.numPaths:", self.numPaths
				
	
	
				" determine which paths are leaves "
				pathIDs = self.paths.getPathIDs()
				isAParent = {}
				for k in pathIDs:
					isAParent[k] = False
				for k in orderedPathIDs1:
					print "index:", k
					currPath = self.paths.getPath(k)
					currParent = currPath["parentID"]
					if currParent != None:
						isAParent[currParent] = True
				
				"add nodes to paths that are the leaves "
				for pathID in orderedPathIDs1:
					if not isAParent[pathID]:				
						self.paths.addNode(nodeID1,pathID)

				pathIDs = self.paths.getPathIDs()
				isAParent = {}
				for k in pathIDs:
					isAParent[k] = False
				for k in orderedPathIDs2:
					print "index:", k
					currPath = self.paths.getPath(k)
					currParent = currPath["parentID"]
					if currParent != None:
						isAParent[currParent] = True

				for pathID in orderedPathIDs2:
					if not isAParent[pathID]:				
						self.paths.addNode(nodeID2, pathID)


				self.paths.generatePaths()
				self.trimmedPaths = self.paths.trimPaths(self.paths.paths)		
				#self.paths.comparePaths()
				self.mergePaths()

				paths = {}
				pathIDs = self.paths.getPathIDs()

				" remove pathIDs that have been merged "
				toBeRemoved = []
				for k in range(len(orderedPathIDs1)):
					currID = orderedPathIDs1[k]
					if not currID in pathIDs:
						toBeRemoved.append(k)

				" reverse so we don't change the indexing when deleting entries "
				toBeRemoved.reverse()	
				for k in toBeRemoved:
					del orderedPathIDs1[k]
					del departures1[k]
					del interiors1[k]
					del depPoints1[k]

				toBeRemoved = []
				for k in range(len(orderedPathIDs2)):
					currID = orderedPathIDs2[k]
					if not currID in pathIDs:
						toBeRemoved.append(k)

				" reverse so we don't change the indexing when deleting entries "
				toBeRemoved.reverse()	
				for k in toBeRemoved:
					del orderedPathIDs2[k]
					del departures2[k]
					del interiors2[k]
					del depPoints2[k]
	
				#paths = {}
				#pathIDs = self.paths.getPathIDs()
				for k in pathIDs:
					path = self.paths.paths[k]
					if len(path) > 0:
						paths[k] = path	
				self.trimmedPaths = self.paths.trimPaths(paths)
	
				" perform path constraints only if we have not created a new path "
				#isNew1 = False
				#isNew2 = False
				#for pathID in newPaths:
				#	if orderedPathIDs1.count(pathID) > 0:
				#		isNew1 = True
				#	if orderedPathIDs2.count(pathID) > 0:
				#		isNew2 = True
				
	
				
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)
	
	
	
				"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
				"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
				"3)  select choice with the lowest cost "
				"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
				
				" nodeID1:  is it in a junction or a single path? "
				
				
				if abc > 0:
					if not isNew1:
						self.constrainToPaths(nodeID1, orderedPathIDs1, departures1, interiors1, depPoints1, insertNode = True)
						
					if not isNew2:
						self.constrainToPaths(nodeID2, orderedPathIDs2, departures2, interiors2, depPoints2, insertNode = True)
				else:
					if not isNew1:
						self.constrainToPaths(nodeID1, orderedPathIDs1, departures1, interiors1, depPoints1, insertNode = False)
						
					if not isNew2:
						self.constrainToPaths(nodeID2, orderedPathIDs2, departures2, interiors2, depPoints2, insertNode = False)
					
		
		self.mergePriorityConstraints()
		
		self.paths.generatePaths()
			
		self.drawPathAndHull()
		
		return




	def integrateNode(self, newNode, nodeID):

		global globalFunc
		global globalArg



		" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
		direction = newNode.travelDir

		" ensure the axes are computed before this check "
		computeHullAxis(nodeID, newNode, tailCutOff = False)
		
		
		if nodeID > 0:
		
			#self.addBatchConstraints(nodeID)

			
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
				#self.mergePriorityConstraints()
								
				" ODD NUMBER POSES RECEIVE INPLACE CONSTRAINT WITH EVEN PAIR "
				" PERFORM MEDIAL OVERLAP WITH PREVIOUS ODD NUMBER POSE "
			else:
			
				" FIXME:  apply this support check to spliced paths "
				supportLine = self.paths.paths[self.currPath]
				
				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				#self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY2)
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

				
				m = hashlib.md5()
				m.update(repr(self.nodeHash[nodeID-1].medialPathA))
				print nodeID-1, "medialPathA", int(m.digest().encode('hex'),16)

				m = hashlib.md5()
				m.update(repr(self.nodeHash[nodeID].medialPathA))
				print nodeID, "medialPathA", int(m.digest().encode('hex'),16)

				m = hashlib.md5()
				m.update(repr(supportLine))
				print "supportLine", int(m.digest().encode('hex'),16)


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
				#self.mergePriorityConstraints()
		
		self.updateLastNode(nodeID)


		" CHECK FOR A BRANCHING EVENT "
		
		if self.numNodes >= 4 and self.numNodes % 2 == 0:

			" DETECT BRANCHING EVENTS FOR THE 2 NODES OF LAST STOP "
			" AFTER ODD-NUMBER OVERLAP OF CURRENT STOP HAS IRONED OUT ERRORS "
			nodeID1 = self.numNodes-4
			nodeID2 = self.numNodes-3
			
			" if these nodes are already path-classified, return"
			isContained1 = False
			isContained2 = False
			
			pathIDs = self.paths.getPathIDs()
			for k in pathIDs:
				if self.paths.getNodes(k).count(nodeID1) > 0:
					isContained1 = True
				if self.paths.getNodes(k).count(nodeID2) > 0:
					isContained2 = True
					
				
			if isContained1 or isContained2:
				return

			" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
			self.paths.generatePaths()

				
			self.drawPathAndHull()


			" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
			if len(self.paths.paths[0]) == 0:
				" first nodes in path 0 "								
				self.paths.addNode(nodeID1,0)
				self.paths.addNode(nodeID2,0)

				" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
			else:
				
				" departure events for node 1 "
				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
				contig1 = []
	
				" departure events for node 2 "
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []
				contig2 = []

				
				" raw medial axis of union of path nodes "
				paths = {}
				pathIDs = self.paths.getPathIDs()

				for k in pathIDs:
					path = self.paths.paths[k]
					if len(path) > 0:
						paths[k] = path
	
				print "paths:", len(paths)
				for k in pathIDs:
					print k, len(paths[k])
	
	
				" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
				self.trimmedPaths = self.paths.trimPaths(paths)
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)
	
				
				" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
				
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
				
				print nodeID1, "orderedPathIDs1:", orderedPathIDs1
				print nodeID2, "orderedPathIDs2:", orderedPathIDs2
				
				
				" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
					contig1.append((contigFrac, overlapSum))
			
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					contig2.append((contigFrac, overlapSum))
					
				print "node departures", nodeID1, ":", departures1
				print "node  interiors", nodeID1, ":", interiors1
				print "node departures", nodeID2, ":", departures2
				print "node  interiors", nodeID2, ":", interiors2

				print "node contiguity", nodeID1, ":", contig1
				print "node contiguity", nodeID2, ":", contig2
	
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

				if contig1[0][0] < 0.4 or contig1[-1][0] < 0.4:
					adjustPose1 = True
				else:
					adjustPose1 = False

				if contig2[0][0] < 0.4 or contig2[-1][0] < 0.4:
					adjustPose2 = True
				else:
					adjustPose2 = False
					
				print "integrate pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1, contig1[0][0], contig1[-1][0]
				print "integrate pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2, contig2[0][0], contig2[-1][0]
				
				" check for the cases that internal departure may have a departure point discrepancy, "
				" then we should perform a pose adjustment and recompute departures "
				print "checking discrepancy in departure:", nodeID1, nodeID2, foreTerm1, discForeTerm1, backTerm1, discBackTerm1, foreTerm2, discForeTerm2, backTerm2, discBackTerm2

				if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2 or adjustPose1 or adjustPose2:

					print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2

					#if discForeTerm1 or discBackTerm1:
					if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or adjustPose1:


						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]
						
						print "adjusting pose", nodeID1, "from", self.nodeHash[nodeID1].getGlobalGPACPose()

						splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs1)
						
						self.fitToSplices(nodeID1, splicedPaths1, globalPoint1, globalPoint1)

						"""
						totalGuesses = []
						for path in splicedPaths1:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
							angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID1].setGPACPose(guessPose)
						
						print "set to pose", nodeID1, "to", guessPose
						"""
						
					if foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2 or adjustPose2:


						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]

						print "adjusting pose", nodeID2, "from", self.nodeHash[nodeID2].getGlobalGPACPose()
						
						splicedPaths2 = self.paths.splicePathIDs(orderedPathIDs2)

						self.fitToSplices(nodeID2, splicedPaths2, globalPoint1, globalPoint1)

						"""
						totalGuesses = []
						for path in splicedPaths2:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
							angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID2].setGPACPose(guessPose)

						print "set to pose", nodeID2, "to", guessPose
						#if foreTerm2 and discForeTerm2:
						#	pathID = orderedPathIDs2[0]
						#elif backTerm2 and discBackTerm2:
						#	pathID = orderedPathIDs2[-1]
							
						#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
						#self.nodeHash[nodeID2].setGPACPose(guessPose1)
						
						"""

					departures1 = []
					interiors1 = []
					depPoints1 = []
					distances1 = []
					depAngles1 = []
					contig1 = []
					
		
					departures2 = []
					interiors2 = []
					depPoints2 = []
					distances2 = []
					depAngles2 = []
					contig2 = []
					
					" the overlapping paths are computed from the initial guess of position "
					orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
					orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
					
					print nodeID1, "orderedPathIDs1:", orderedPathIDs1
					print nodeID2, "orderedPathIDs2:", orderedPathIDs2
					
					" now compute whether there are departure points after we have guessed a better position in synch with the paths "
					for pathID in orderedPathIDs1:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
						departures1.append([isExist1,isExist2])
						interiors1.append([isInterior1, isInterior2])
						depPoints1.append([departurePoint1, departurePoint2])
						distances1.append([discDist1, discDist2])
						depAngles1.append([depAngle1, depAngle2])
						contig1.append((contigFrac, overlapSum))
	
					
					for pathID in orderedPathIDs2:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
						departures2.append([isExist1,isExist2])
						interiors2.append([isInterior1, isInterior2])
						depPoints2.append([departurePoint1, departurePoint2])
						distances2.append([discDist1, discDist2])
						depAngles2.append([depAngle1, depAngle2])
						contig2.append((contigFrac, overlapSum))
						
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

	
				frontExist1 = departures1[0][0]
				frontInterior1 = interiors1[0][0]
		
				frontExist2 = departures2[0][0]
				frontInterior2 = interiors2[0][0]
				
				depAngle1 = depAngles1[0][0]
				depAngle2 = depAngles2[0][0]
				
				depPoint1 = depPoints1[0][0]
				depPoint2 = depPoints2[0][0]
				
				parentPathID1 = orderedPathIDs1[0]
				parentPathID2 = orderedPathIDs2[0]

				isFront1 = self.nodeHash[nodeID1].faceDir
				isFront2 = self.nodeHash[nodeID2].faceDir
				if isFront1 and not isFront2:
					dirFlag = 0
				elif not isFront1 and isFront2:
					dirFlag = 1
				else:
					print isFront1, isFront2
					raise	
	
				isBranch, pathBranchIDs, isNew = self.paths.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)

				print "determineBranch:"
				print isBranch
				print pathBranchIDs

				isNew1 = isNew[0]
				isNew2 = isNew[1]
	
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
					
				#for pathID in pathBranchIDs:
				#	if pathID != -1 and newPaths.count(pathID) == 0:
				#		newPaths.append(pathID)




				backExist1 = departures1[-1][1]
				backInterior1 = interiors1[-1][1]
		
				backExist2 = departures2[-1][1]
				backInterior2 = interiors2[-1][1]
				
				depAngle1 = depAngles1[-1][1]
				depAngle2 = depAngles2[-1][1]
				
				depPoint1 = depPoints1[-1][1]
				depPoint2 = depPoints2[-1][1]
	
				parentPathID1 = orderedPathIDs1[-1]
				parentPathID2 = orderedPathIDs2[-1]
				#print "self.nodeSets:", self.nodeSets
				#print "newPaths:", newPaths
				#print "self.numPaths:", self.numPaths

				isFront1 = self.nodeHash[nodeID1].faceDir
				isFront2 = self.nodeHash[nodeID2].faceDir
				if isFront1 and not isFront2:
					dirFlag = 1
				elif not isFront1 and isFront2:
					dirFlag = 0
				else:
					print isFront1, isFront2
					raise	
	
				isBranch, pathBranchIDs, isNew = self.paths.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)
				print "determineBranch:"
				print isBranch
				print pathBranchIDs

				isNew1 = isNew1 or isNew[0]
				isNew2 = isNew2 or isNew[1]
	
				if isBranch[0]:
					orderedPathIDs1.append(pathBranchIDs[0])
					departures1.append([False, False])
					interiors1.append([False, False])
					depPoints1.append([None, None])
	
				if isBranch[1]:
					orderedPathIDs2.append(pathBranchIDs[1])
					departures2.append([False, False])
					interiors2.append([False, False])
					depPoints2.append([None, None])
	
				#for pathID in pathBranchIDs:
				#	if pathID != -1 and newPaths.count(pathID) == 0:
				#		newPaths.append(pathID)
				#print "self.nodeSets:", self.nodeSets
				#print "newPaths:", newPaths
				#print "self.numPaths:", self.numPaths
				


				" determine which paths are leaves "
				pathIDs = self.paths.getPathIDs()
				isAParent = {}
				for k in pathIDs:
					isAParent[k] = False
				for k in orderedPathIDs1:
					print "index:", k
					currPath = self.paths.getPath(k)
					currParent = currPath["parentID"]
					if currParent != None:
						isAParent[currParent] = True

				
				"add nodes to paths that are the leaves "
				for pathID in orderedPathIDs1:
					if not isAParent[pathID]:				
						self.paths.addNode(nodeID1,pathID)

				pathIDs = self.paths.getPathIDs()
				isAParent = {}
				for k in pathIDs:
					isAParent[k] = False
				for k in orderedPathIDs2:
					print "index:", k
					currPath = self.paths.getPath(k)
					currParent = currPath["parentID"]
					if currParent != None:
						isAParent[currParent] = True

				for pathID in orderedPathIDs2:
					if not isAParent[pathID]:				
						self.paths.addNode(nodeID2,pathID)


				self.paths.generatePaths()
				self.trimmedPaths = self.paths.trimPaths(self.paths.paths)		
				self.drawTrimmedPaths(self.trimmedPaths)
				#self.paths.comparePaths()
				self.mergePaths()

				paths = {}
				pathIDs = self.paths.getPathIDs()

				" remove pathIDs that have been merged "
				toBeRemoved = []
				for k in range(len(orderedPathIDs1)):
					currID = orderedPathIDs1[k]
					if not currID in pathIDs:
						toBeRemoved.append(k)

				" reverse so we don't change the indexing when deleting entries "
				toBeRemoved.reverse()	
				for k in toBeRemoved:
					del orderedPathIDs1[k]
					del departures1[k]
					del interiors1[k]
					del depPoints1[k]

				toBeRemoved = []
				for k in range(len(orderedPathIDs2)):
					currID = orderedPathIDs2[k]
					if not currID in pathIDs:
						toBeRemoved.append(k)

				" reverse so we don't change the indexing when deleting entries "
				toBeRemoved.reverse()	
				for k in toBeRemoved:
					del orderedPathIDs2[k]
					del departures2[k]
					del interiors2[k]
					del depPoints2[k]


				for k in pathIDs:
					path = self.paths.paths[k]
					if len(path) > 0:
						paths[k] = path	
				self.trimmedPaths = self.paths.trimPaths(paths)

				" perform path constraints only if we have not created a new path "
				#for pathID in newPaths:
				#	if orderedPathIDs1.count(pathID) > 0:
				#		isNew1 = True
				#	if orderedPathIDs2.count(pathID) > 0:
				#		isNew2 = True
				
	
				
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

				self.paths.generatePaths()
				self.drawPathAndHull()
	
				self.updatePathNode(nodeID1, nodeID2, orderedPathIDs1, orderedPathIDs2)
	
				self.updateLastNode(nodeID1+2)
				self.updateLastNode(nodeID2+2)
				
			self.paths.generatePaths()
			self.drawPathAndHull()

			
			"""
			try:
				self.particleFilter
			except:
				self.particleFilter = ParticleFilter(self.paths.paths[0], self.nodeHash)
			else:
				#globalFunc = self.particleFilter.update
				#globalArg = self.paths.paths[0]
				#cProfile.run('globalFunc(globalArg)', 'prof_sim')
				#cProfile.runctx('self.particleFilter.update(self.paths.paths[0],nodeID1)', None, locals(), 'icp_prof')
				#exit()
				self.particleFilter.update(self.paths.paths[0], nodeID1, isForward = direction)
			"""
			if nodeID1 >= 2:
				newPath = deepcopy(self.paths.paths[0])
				p0 = newPath[0]
				pN = newPath[-1]
				
				rootPose = [-2.0,0.0]
				
				dist0 = sqrt((rootPose[0]-p0[0])*(rootPose[0]-p0[0]) + (rootPose[1]-p0[1])*(rootPose[1]-p0[1]))
				distN = sqrt((rootPose[0]-pN[0])*(rootPose[0]-pN[0]) + (rootPose[1]-pN[1])*(rootPose[1]-pN[1]))
		
				if dist0 > distN:
					newPath.reverse()
							
				" new path "
				newSpline = SplineFit(newPath, smooth = 0.1)
	
	
				posNew = self.nodeHash[nodeID1].getGlobalGPACPose()
				posOld = self.nodeHash[nodeID1-2].getGlobalGPACPose()
				minDist, minU, closestPoint = newSpline.findClosestPoint(posNew)
				arcDistNew = newSpline.dist(0.0, minU)
	
				minDist, minU, closestPoint = newSpline.findClosestPoint(posOld)
				arcDistOld = newSpline.dist(0.0, minU)
	
				print "arcDistNew, arcDistOld, diff =", arcDistNew, arcDistOld, arcDistNew-arcDistOld
					
			return

	def mergePaths(self):


		def dispOffset(p, offset):
			xd = offset[0]
			yd = offset[1]
			theta = offset[2]
		
			px = p[0]
			py = p[1]
			pa = p[2]
			
			newAngle = normalizeAngle(pa + theta)
			
			p_off = [px*math.cos(theta) - py*math.sin(theta) + xd, px*math.sin(theta) + py*math.cos(theta) + yd, newAngle]
			
			return p_off

		

		toBeMerged = self.paths.comparePaths()
							
		for mergeThis in toBeMerged:

			pathID1 = mergeThis[0]
			pathID2 = mergeThis[1]
			offset = mergeThis[2]

			" verify that this path still exists before we try to merge it "
			allPathIDs = self.paths.getPathIDs()
			if pathID1 in allPathIDs and pathID2 in allPathIDs:
	
				
				" capture transform between pathID1 and pathID2 under correction "

				cOffset = offset
				rootID = pathID1
				lessID = pathID2
				
				
				mergeNodes = self.paths.getNodes(lessID)
				
				
				" capture transform between pathID2 and member nodes "
				" add compute new node locations under corrected transform "
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				
				nodeHasMerged = {}
				
				" offset between paths "
				poseOrigin = Pose(cOffset)
				for nodeID in mergeNodes:

					" adjust node poses to relative path positions "
					nodePose = self.nodeHash[nodeID].getGlobalGPACPose()
					guessPose = poseOrigin.convertLocalOffsetToGlobal(nodePose)
					self.nodeHash[nodeID].setGPACPose(guessPose)
					
					self.drawConstraints(self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull()
					
					pathNodes = self.paths.getNodes(rootID)
					targetNodeID = nodeID
					
					print "addPathConstraints:"
					print "pathNodes:", pathNodes
					print "targetNodeID:", targetNodeID
					
					#newConstraints = self.addTargetBatchConstraints(targetNodeID, pathNodes)
					newConstraints = []
					
					for const in newConstraints:
						n1 = const[0]
						n2 = const[1]
						transform = const[2]
						covE = const[3]
			
						self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

					if len(newConstraints) > 0:
						nodeHasMerged[nodeID] = True
					else:
						nodeHasMerged[nodeID] = False
			
					
				

					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID].getGlobalGPACPose()[:2]					
					print "adjusting pose", nodeID, "from", self.nodeHash[nodeID].getGlobalGPACPose()
					splicedPaths1 = self.paths.splicePathIDs([rootID])

					self.fitToSplices(nodeID, splicedPaths1, globalPoint1, globalPoint1)

					"""
					totalGuesses = []
					for path in splicedPaths1:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					guessPose = totalGuesses[0][2]
					self.nodeHash[nodeID].setGPACPose(guessPose)
					print "set to pose", nodeID, "to", guessPose
					"""

				if self.numNodes-4 in mergeNodes:
					self.updateLastNode(self.numNodes-2)
				if self.numNodes-3 in mergeNodes:
					self.updateLastNode(self.numNodes-1)

				" merge the constraints "
				self.mergePriorityConstraints()

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()


				" relax non-path constraints so that rest of graph can adjust smoothly "
				inplace_edges = self.getPriorityEdges(INPLACE_PRIORITY)
				overlap_edges = self.getPriorityEdges(OVERLAP_PRIORITY)

				" between mergedNodes and their constrained partners, "
				" loosen the coaxial covariance all along the chain "
				" Only disconnect edges if a connection has been made during the merge "
				
				for nodeID in mergeNodes:
					
					if nodeHasMerged[nodeID]:
						for k, v in overlap_edges.items():
							if len(v) > 0:
								#print "edges:", v
								covE = v[0][1]
								if k[0] == nodeID or k[1] == nodeID:
									self.relaxPriorityEdge(k[0],k[1],priorityLevel = OVERLAP_PRIORITY)
									#covE[0,0] *= 100
									#covE[1,1] *= 100
						for k, v in inplace_edges.items():
							if len(v) > 0:
								#print "edges:", v
								covE = v[0][1]
								if k[0] == nodeID or k[1] == nodeID:
									self.relaxPriorityEdge(k[0],k[1],priorityLevel = INPLACE_PRIORITY)


				" k = (nodeID1,nodeID2)  v = [transform1,covE1]"
				"""
				for nodeID in mergeNodes:
					for k, v in inplace_edges.items():
						if len(v) > 0:
							print "edges:", v
							covE = v[0][1]
							if k[0] == nodeID or k[1] == nodeID:
								covE[0,0] *= 100
								covE[1,1] *= 100

					for k, v in overlap_edges.items():
						if len(v) > 0:
							print "edges:", v
							covE = v[0][1]
							if k[0] == nodeID or k[1] == nodeID:
								covE[0,0] *= 100
								covE[1,1] *= 100
				"""
				
				" re-evaluate constraints with relaxed overlap and inplace priorities"
				self.mergePriorityConstraints()

				print "relaxed constraints",  self.statePlotCount
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				
				" constrain nodes to merge to pathID1 "
				
				
				" move nodes from lessID into rootID and delete pathID "
				for nodeID in mergeNodes:
					self.paths.addNode(nodeID, rootID)
	
				" move nodes to pathID1, delete pathID2 "
				self.paths.delPath(lessID, rootID)
	
				self.paths.generatePaths()
				
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				

	def constrainToPaths(self, nodeID, orderedPathIDs, departures, interiors, depPoints, insertNode = False):
			
		"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
		"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
		"3)  select choice with the lowest cost "
		"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
		
		" nodeID1:  is it in a junction or a single path? "
		
		splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs)

		if len(orderedPathIDs) == 1:

			" in a single path "

			pathID = orderedPathIDs[0]
			#nodeSet = self.nodeSets[pathID]
			nodeSet = self.paths.getNodes(pathID)
			
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

				self.fitToSplices(nodeID, [self.trimmedPaths[pathID]], globalPoint1, globalPoint1)
					
				#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
				#self.nodeHash[nodeID].setGPACPose(guessPose1)

				#print "set to pose", nodeID, "to", guessPose1
				
				" make path constraint "
				print "pathID: addPathConstraints3(", pathID, nodeID
				#self.addPathConstraints3(self.paths.getNodes(pathID), nodeID, insertNode = insertNode, orderedPathIDs = orderedPathIDs)
			
		else:
			" in at least one junction "
			
			" we don't intersect the junction sufficiently "
			" add constraints to most overlapped path"
			#if interiors.count(True) == 0:
			" no recognized internal departure point "
			if [item for inner_list in interiors for item in inner_list].count(True) == 0:


				" determine which paths are leaves "
				#isAChild = [True for k in range(self.numPaths)]
				#for k in orderedPathIDs:
				#	currPath = self.pathParents[k]
				#	currParent = currPath[0]
				#	if currParent != None and orderedPathIDs.count(currParent) > 0:
				#		isAChild[currParent] = False


				pathIDs = self.paths.getPathIDs()
				isAChild = {}
				for k in pathIDs:
					isAChild[k] = True
				for k in orderedPathIDs:
					currPath = self.paths.getPath(k)
					currParent = currPath["parentID"]
					if currParent != None and currParent in orderedPathIDs:
						isAChild[currParent] = False




				print "isAChild:", isAChild


				" of all the paths that are children, pick the lowest overlap sum one "
				pathCandidates = []
				for pathID in orderedPathIDs:
					if isAChild[pathID]:
						pathCandidates.append(pathID)
				
				" initialize to the first child path "
				minPathID = pathCandidates[0]
				minSum = 1e100
				
				for pathID in pathCandidates:
					
					if isAChild[pathID]:
						sum1 = self.paths.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
						
						if sum1 < minSum:
							minPathID = pathID
							minSum = sum1
						print "overlap sum", pathID, "=", sum1						

				print "maximum overlap path is", minPathID

				nodeSet = self.paths.getNodes(minPathID)
				
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

					self.fitToSplices(nodeID, splicedPaths1, globalPoint1, globalPoint1)

					" make path constraint "
					print "minPathID: addPathConstraints3(", minPathID, nodeID
					#self.addPathConstraints3(self.paths.getNodes(minPathID), nodeID, insertNode = insertNode)






			
			elif len(orderedPathIDs) == 2:
				
				" one junction point "
				pathID1 = orderedPathIDs[0]
				pathID2 = orderedPathIDs[1]

				#[pathID, branchNodeID, poseOrigin.convertGlobalToLocal(globalJunctionPoint)]
				
				" find parent-child relationship "
				#val1 = self.pathParents[pathID1]
				#val2 = self.pathParents[pathID2]
				val1 = self.paths.getPath(pathID1)
				val2 = self.paths.getPath(pathID2)


				if val1['parentID'] == None:
					parentPath = pathID1
					childPath = pathID2
					
				elif val2['parentID'] == None:					
					parentPath = pathID2
					childPath = pathID1
					
				else:
					parentID1 = val1['parentID']
					parentID2 = val2['parentID']
					
					if parentID1 == pathID2:
						parentPath = pathID2
						childPath = pathID1
					elif parentID2 == pathID1:
						parentPath = pathID1
						childPath = pathID2
					else:
						print "child-parent status not established for paths:", orderedPathIDs
						raise

				" global path point "
				localJunctionNodeID = self.paths.getPath(childPath)["branchNodeID"]
				localPathJunctionPoint = self.paths.getPath(childPath)["localJunctionPose"]						
				
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

					self.fitToSplices(nodeID, splicedPaths1, globalPathJunctionPoint, depPoint1)

					" make path constraint "
					print "childPath: addPathConstraints3(", childPath, nodeID
					#self.addPathConstraints3(self.paths.getNodes(childPath), nodeID, insertNode = insertNode)

				
				elif isDepExists2:

					self.fitToSplices(nodeID, splicedPaths1, globalPathJunctionPoint, depPoint2)

					" make path constraint "
					print "childPath: addPathConstraints3(", childPath, nodeID
					#self.addPathConstraints3(self.paths.getNodes(childPath), nodeID, insertNode = insertNode)

				
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
				#val1 = self.pathParents[pathID1]
				#val2 = self.pathParents[pathID2]
				val1 = self.paths.getPath(pathID1)
				val2 = self.paths.getPath(pathID2)


				if val1['parentID'] == None:
					parentPath = pathID1
					childPath = pathID2

				elif val2['parentID'] == None:					
					parentPath = pathID2
					childPath = pathID1

				else:
					parentID1 = val1['parentID']
					parentID2 = val2['parentID']
					
					if parentID1 == pathID2:
						parentPath = pathID2
						childPath = pathID1
					elif parentID2 == pathID1:
						parentPath = pathID1
						childPath = pathID2
					else:
						raise

				cPath = self.paths.getPath(childPath)
				if cPath['parentID'] != parentPath:
					print "WARNING: selected childPath", childPath, "does not have parent", parentPath, "but has instead", self.pathParents[childPath][0]

				" global path point "
				" we don't check, but we assume the parent here is the parent in our orderedPaths"
				localJunctionNodeID = cPath["branchNodeID"]
				localPathJunctionPoint = cPath["localJunctionPose"]						
				
				localJunctionNode = self.nodeHash[localJunctionNodeID]
				poseOrigin = Pose(localJunctionNode.getEstPose())
				globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)

				
				self.fitToSplices(nodeID, splicedPaths1, globalPathJunctionPoint, globalPathJunctionPoint)
			
				" make path constraint "
				print "childPath: addPathConstraints3(", childPath, nodeID
				#self.addPathConstraints3(self.paths.getNodes(childPath), nodeID, insertNode = insertNode)



	def checkForBranches(self, nodeID1, nodeID2):

		" CHECK FOR A BRANCHING EVENT FROM LEAF "

		" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
		if len(self.paths.paths[0]) == 0:
			" first nodes in path 0 "					
			self.paths.addNode(nodeID1,0)
			self.paths.addNode(nodeID2,0)

			" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
		else:
			
			" departure events for node 1 "
			departures1 = []
			interiors1 = []
			depPoints1 = []
			distances1 = []
			depAngles1 = []
			contig1 = []

			" departure events for node 2 "
			departures2 = []
			interiors2 = []
			depPoints2 = []
			distances2 = []
			depAngles2 = []
			contig2 = []

			
			" raw medial axis of union of path nodes "
			paths = {}
			pathIDs = self.paths.getPathIDs()
			for k in pathIDs:
				path = self.paths.paths[k]
				if len(path) > 0:
					paths[k] = path

			print "paths:", len(paths)
			for k in pathIDs:
				print k, len(paths[k])


			" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
			self.trimmedPaths = self.paths.trimPaths(paths)
			self.drawTrimmedPaths(self.trimmedPaths)
			print "trimmed paths:", len(self.trimmedPaths)

			
			" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
				
			" the overlapping paths are computed from the initial guess of position "
			orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
			orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
			
			print "orderedPathIDs1:", orderedPathIDs1
			print "orderedPathIDs2:", orderedPathIDs2
			
			
			" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
			for pathID in orderedPathIDs1:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
				departures1.append([isExist1,isExist2])
				interiors1.append([isInterior1, isInterior2])
				depPoints1.append([departurePoint1, departurePoint2])
				distances1.append([discDist1, discDist2])
				depAngles1.append([depAngle1, depAngle2])
				contig1.append((contigFrac, overlapSum))
			
			for pathID in orderedPathIDs2:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
				departures2.append([isExist1,isExist2])
				interiors2.append([isInterior1, isInterior2])
				depPoints2.append([departurePoint1, departurePoint2])
				distances2.append([discDist1, discDist2])
				depAngles2.append([depAngle1, depAngle2])
				contig2.append((contigFrac, overlapSum))
				
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


			print "branchCheck:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1, contig1[0][0], contig1[-1][0]
			print "branchCheck:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2, contig2[0][0], contig2[-1][0]

			
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

					splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs1)

					totalGuesses = []
					for path in splicedPaths1:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
						angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
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


					splicedPaths2 = self.paths.splicePathIDs(orderedPathIDs2)

					totalGuesses = []
					for path in splicedPaths2:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
						angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
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
				contig1 = []
	
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []
				contig2 = []
				
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
				
				print "orderedPathIDs1:", orderedPathIDs1
				print "orderedPathIDs2:", orderedPathIDs2
				
				" now compute whether there are departure points after we have guessed a better position in synch with the paths "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
					contig1.append((contigFrac, overlapSum))
				
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					contig2.append((contigFrac, overlapSum))
					
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
			frontInterior2 = interiors2[0][0]
			
			depAngle1 = depAngles1[0][0]
			depAngle2 = depAngles2[0][0]
			
			depPoint1 = depPoints1[0][0]
			depPoint2 = depPoints2[0][0]
			
			parentPathID1 = orderedPathIDs1[0]
			parentPathID2 = orderedPathIDs2[0]

			
			isFront1 = self.nodeHash[nodeID1].faceDir
			isFront2 = self.nodeHash[nodeID2].faceDir
			if isFront1 and not isFront2:
				dirFlag = 0
			elif not isFront1 and isFront2:
				dirFlag = 1
			else:
				print isFront1, isFront2
				raise

			isBranch, pathBranchIDs, isNew = self.paths.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)

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
			backInterior2 = interiors2[-1][1]
			
			depAngle1 = depAngles1[-1][1]
			depAngle2 = depAngles2[-1][1]
			
			depPoint1 = depPoints1[-1][1]
			depPoint2 = depPoints2[-1][1]

			parentPathID1 = orderedPathIDs1[-1]
			parentPathID2 = orderedPathIDs2[-1]

			isFront1 = self.nodeHash[nodeID1].faceDir
			isFront2 = self.nodeHash[nodeID2].faceDir
			if isFront1 and not isFront2:
				dirFlag = 1
			elif not isFront1 and isFront2:
				dirFlag = 0
			else:
				print isFront1, isFront2
				raise

			isBranch, pathBranchIDs, isNew = self.paths.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)

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
			
			
			pathIDs = self.paths.getPathIDs()
			isAParent = {}
			for k in pathIDs:
				isAParent[k] = False
			for k in orderedPathIDs1:
				print "index:", k
				currPath = self.paths.getPath(k)
				currParent = currPath["parentID"]
				if currParent != None:
					isAParent[currParent] = True

		
			"add nodes to paths that are the leaves "
			for pathID in orderedPathIDs1:
				if not isAParent[pathID]:				
					self.paths.addNode(nodeID1,pathID)

			pathIDs = self.paths.getPathIDs()
			isAParent = {}
			for k in pathIDs:
				isAParent[k] = False
			for k in orderedPathIDs2:
				print "index:", k
				currPath = self.paths.getPath(k)
				currParent = currPath["parentID"]
				if currParent != None:
					isAParent[currParent] = True



			for pathID in orderedPathIDs2:
				if not isAParent[pathID]:				
					self.paths.addNode(nodeID2,pathID)

			paths = {}
			pathIDs = self.paths.getPathIDs()
			
			for k in pathIDs:
				path = self.paths.paths[k]
				if len(path) > 0:
					paths[k] = path	
			self.trimmedPaths = self.paths.trimPaths(paths)

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

			
		self.paths.generatePaths()

		self.drawPathAndHull()
		
		return
	
	
	def updateLastNode(self, nodeID):
		
		" find constraints of nodeID and update it according to constraint from nodeID-2"
		
		
		if nodeID >= 2:
			
			resultEdges = self.getEdges(nodeID-2, nodeID) 
			
			resultEdges.sort(key=itemgetter(2), reverse=True)

			edge = resultEdges[0]
			
			transform = edge[0]
			offset = [transform[0,0],transform[1,0],transform[2,0]]

			"compute pose nodeID wrt node-2 given offset "
			
			origPose = self.nodeHash[nodeID-2].getGlobalGPACPose()
			
			originPose = Pose(origPose)
			nodePose = originPose.convertLocalOffsetToGlobal(offset)
			
			self.nodeHash[nodeID].setGPACPose(nodePose)

			print "update position of node", nodeID, "from", origPose, "to", nodePose

	def updatePathNode(self, nodeID1, nodeID2, orderedPathIDs1, orderedPathIDs2):

		"""
		GOALS:
		
		1) node medial axes should be 100% contiguously overlapping path
		2) pair nodes should be approximately in the same location and orientation
		
		"""
		
		print "updatePathNode(", nodeID1, nodeID2, orderedPathIDs1, orderedPathIDs2, ")"


		"""
		departures1 = []
		interiors1 = []
		depPoints1 = []
		distances1 = []
		depAngles1 = []
		isDep = []
		contig1 = []
		sums = []
		counts = []
			
		for nodeID in nodeSet:

			departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1, plotIter = False)
			departures1.append([isExist1,isExist2])
			interiors1.append([isInterior1, isInterior2])
			depPoints1.append([departurePoint1, departurePoint2])
			distances1.append([discDist1, discDist2])
			depAngles1.append([depAngle1, depAngle2])
			contig1.append((contigFrac, overlapSum))
			isDep.append(isExist1 and isInterior1 or isExist2 and isInterior2)

			sum1 = self.paths.getOverlapCondition(self.trimmedPaths[pathID], nodeID, plotIter = False)
			sums.append(sum1)
		
			count = self.paths.getOverlapCondition2(self.trimmedPaths[pathID], nodeID, plotIter = False)
			counts.append(count)
		"""


		if nodeID1 > nodeID2:
			print "nodes called out of order", nodeID1, nodeID2
			raise

	

		" INPLACE CONSTRAINT BETWEEN PAIRED NODES "
		resultEdges = self.getEdges(nodeID1, nodeID2)
		resultEdges.sort(key=itemgetter(2), reverse=True)
		edge = resultEdges[0]
		
		transform = edge[0]
		offset = [transform[0,0],transform[1,0],transform[2,0]]
		
		nomAngDiff = offset[2]

		origPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()
		origPose2 = self.nodeHash[nodeID2].getGlobalGPACPose()


		#if edge[0] == nodeID1 and edge[1] == nodeID2:
		if True:
			offset12 = offset
			offsetPose12 = Pose(offset12)
						
			offset21 = offsetPose12.convertGlobalPoseToLocal([0.0,0.0,0.0])
			
		#elif edge[0] == nodeID2 and edge[1] == nodeID1:
		elif False:
			offset21 = offset
			offsetPose21 = Pose(offset21)
						
			offset12 = offsetPose21.convertGlobalPoseToLocal([0.0,0.0,0.0])

		else:
			print "using edge:", edge
			print "all edges:"
			print resultEdges
			raise


		print "update position of path node", nodeID1, nodeID2
	
		" control point nearest the GPAC origin "
		globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]

		print "adjusting pose", nodeID1, "from", origPose1

		splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs1)
		print "received", len(splicedPaths1), "spliced paths from path IDs", orderedPathIDs1

		#self.fitToSplices(nodeID1, splicedPaths1, globalPoint1, globalPoint1)
		
		
		"""
		totalGuesses1 = []
		for path in splicedPaths1:

			" make the path constraints "								
			guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)

			self.nodeHash[nodeID1].setGPACPose(guessPose1)

			resultArgs = self.paths.getDeparturePoint(path, nodeID1, plotIter = True)

			angDiff1 = abs(diffAngle(origPose1[2], guessPose1[2]))
			totalGuesses1.append((angDiff1, cost1, guessPose1, resultArgs, path))
		
		totalGuesses1.sort()

		print "totalGuesses1:", totalGuesses1
		"""



		totalGuesses1 = []
		
		for k in range(len(splicedPaths1)):

			path = splicedPaths1[k]

			" make the path constraints "								
			guessPoseA, costA = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)

			self.nodeHash[nodeID1].setGPACPose(guessPoseA)
			resultArgs = self.paths.getDeparturePoint(path, nodeID1, plotIter = True)

			angDiffA = abs(diffAngle(origPose1[2], guessPoseA[2]))

			totalGuesses1.append((angDiffA, costA, guessPoseA, resultArgs, path, k))

			self.nodeHash[nodeID1].setGPACPose(origPose1)
		
		totalGuesses1.sort()
		print "totalGuesses1:", totalGuesses1

		angDiff1 = totalGuesses1[0][0]
		guessPose1 = totalGuesses1[0][2]
		resultArgs1 = totalGuesses1[0][3]
		path1 = totalGuesses1[0][4]
		
		if len(splicedPaths1) > 1:
			angDiffA = totalGuesses1[0][0]
			angDiffB = totalGuesses1[1][0]
			
			print "angDiffs:", angDiffA, angDiffB
			if fabs(angDiffA-angDiffB) < 0.1:
				resultArgsA = totalGuesses1[0][3]
				resultArgsB = totalGuesses1[1][3]
				#resultArgsA = self.paths.getDeparturePoint(splicedPaths1[totalGuesses1[0][5]], nodeID1, plotIter = True)
				#resultArgsB = self.paths.getDeparturePoint(splicedPaths1[totalGuesses1[1][5]], nodeID1, plotIter = True)
				
				"select the most contiguous overlap "
				contigFracA = resultArgsA[10]
				contigFracB = resultArgsB[10]
				
				print "contigFracs:", contigFracA, contigFracB
				
				if contigFracA < contigFracB:
					guessPose1 = totalGuesses1[1][2]
					angDiff1 = totalGuesses1[1][0]
					resultArgs1 = totalGuesses1[1][3]
					path1 = totalGuesses1[1][4]




		





		globalPoint2 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]
		
		print "adjusting pose", nodeID2, "from", origPose2

		splicedPaths2 = self.paths.splicePathIDs(orderedPathIDs2)
		print "received", len(splicedPaths2), "spliced paths from path IDs", orderedPathIDs2


		totalGuesses2 = []
		
		for k in range(len(splicedPaths2)):

			path = splicedPaths2[k]


			" make the path constraints "								
			guessPoseA, costA = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint2, globalPoint2)

			self.nodeHash[nodeID2].setGPACPose(guessPoseA)
			resultArgs = self.paths.getDeparturePoint(path, nodeID2, plotIter = True)

			angDiffA = abs(diffAngle(origPose2[2], guessPoseA[2]))

			totalGuesses2.append((angDiffA, costA, guessPoseA, resultArgs, path, k))

			self.nodeHash[nodeID2].setGPACPose(origPose2)

		
		totalGuesses2.sort()
		print "totalGuesses2:", totalGuesses2

		angDiff2 = totalGuesses2[0][0]
		guessPose2 = totalGuesses2[0][2]
		resultArgs2 = totalGuesses2[0][3]
		path2 = totalGuesses2[0][4]
		
		if len(splicedPaths2) > 1:
			angDiffA = totalGuesses2[0][0]
			angDiffB = totalGuesses2[1][0]
			
			print "angDiffs:", angDiffA, angDiffB
			if fabs(angDiffA-angDiffB) < 0.1:
				resultArgsA = totalGuesses2[0][3]
				resultArgsB = totalGuesses2[1][3]
				#resultArgsA = self.paths.getDeparturePoint(splicedPaths2[totalGuesses2[0][5]], nodeID2, plotIter = True)
				#resultArgsB = self.paths.getDeparturePoint(splicedPaths2[totalGuesses2[1][5]], nodeID2, plotIter = True)
				
				"select the most contiguous overlap "
				contigFracA = resultArgsA[10]
				contigFracB = resultArgsB[10]
				
				print "contigFracs:", contigFracA, contigFracB
				
				if contigFracA < contigFracB:
					guessPose2 = totalGuesses2[1][2]
					angDiff2 = totalGuesses2[1][0]
					resultArgs2 = totalGuesses2[1][3]
					path2 = totalGuesses2[1][4]








		"""

		totalGuesses2 = []
		for path in splicedPaths2:

			" make the path constraints "								
			guessPose2, cost2 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint2, globalPoint2)

			self.nodeHash[nodeID2].setGPACPose(guessPose2)
			resultArgs = self.paths.getDeparturePoint(path, nodeID2, plotIter = True)
			
			angDiff2 = abs(diffAngle(origPose2[2], guessPose2[2]))
			totalGuesses2.append((angDiff2, cost2, guessPose2, resultArgs, path))
		
		totalGuesses2.sort()
		print "totalGuesses2:", totalGuesses2

		angDiff1 = totalGuesses1[0][0]
		angDiff2 = totalGuesses2[0][0]

		resultArgs1 = totalGuesses1[0][3]
		resultArgs2 = totalGuesses2[0][3]

		path1 = totalGuesses1[0][4]
		path2 = totalGuesses2[0][4]

		guessPose1 = totalGuesses1[0][2]
		guessPose2 = totalGuesses2[0][2]
		"""

		originPose1 = Pose(origPose1)
		originPose2 = Pose(origPose2)
		
		if angDiff1 < angDiff2:
			
			nodePose2 = originPose1.convertLocalOffsetToGlobal(offset12)	
					
			inPlaceDiff2 = abs(guessPose2[2]-nodePose2[2])
			print "inPlaceDiff2 = |", guessPose2[2], "-", nodePose2[2], "| =", inPlaceDiff2
			
			if inPlaceDiff2 > pi/4.0:
				self.nodeHash[nodeID2].setGPACPose(nodePose2)
				print "set to pose", nodeID2, "to", nodePose2
				
			else:
				self.nodeHash[nodeID2].setGPACPose(guessPose2)
				print "set to pose", nodeID2, "to", guessPose2
						

			self.nodeHash[nodeID1].setGPACPose(guessPose1)
			print "set to pose", nodeID1, "to", guessPose1

		
		else:
			
			nodePose1 = originPose2.convertLocalOffsetToGlobal(offset21)			

			inPlaceDiff1 = abs(guessPose1[2]-nodePose1[2])
			print "inPlaceDiff1 = |", guessPose1[2], "-", nodePose1[2], "| =", inPlaceDiff1
			
			if inPlaceDiff1 > pi/4.0:
				self.nodeHash[nodeID1].setGPACPose(nodePose1)
				print "set to pose", nodeID1, "to", nodePose1
				
			else:
				self.nodeHash[nodeID1].setGPACPose(guessPose1)
				print "set to pose", nodeID1, "to", guessPose1

			self.nodeHash[nodeID2].setGPACPose(guessPose2)
			print "set to pose", nodeID2, "to", guessPose2
		

		"""
		departurePoint1 = resultArgs[0]
		depAngle1 = resultArgs[1]
		isInterior1 = resultArgs[2]
		isExist1 = resultArgs[3]
		discDist1 = resultArgs[4]
		departurePoint2 = resultArgs[5]
		depAngle2 = resultArgs[6]
		isInterior2 = resultArgs[7]
		isExist2 = resultArgs[8]
		discDist2 = resultArgs[9]
		contigFrac = resultArgs[10]
		overlapSum = resultArgs[11]
		"""


		resultArgs1 = self.paths.getDeparturePoint(path, nodeID1, plotIter = True)
		resultArgs2 = self.paths.getDeparturePoint(path, nodeID2, plotIter = True)


		contigFrac1 = resultArgs1[10]
		contigFrac2 = resultArgs2[10]
		
		print "resultArgs1:" , resultArgs1
		print "resultArgs2:" , resultArgs2
		
		if contigFrac1 < 0.98:
			self.findPathLocation(nodeID1, path1)
		if contigFrac2 < 0.98:
			self.findPathLocation(nodeID2, path2)
		


	def findPathLocation(self, nodeID, path):


		print "findPathLocation(", nodeID
		
		node1 = self.nodeHash[nodeID]
		hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = False)
		medialSpline1 = SplineFit(medial1, smooth=0.1)

		estPose1 = node1.getGlobalGPACPose()		
		poseOrigin = Pose(estPose1)
		
		globalMedial = []
		for p in medial1:
			globalMedial.append(poseOrigin.convertLocalToGlobal(p))

		globalMedialSpline1 = SplineFit(globalMedial, smooth=0.1)
					
		originU1 = medialSpline1.findU([0.0,0.0])	



		pathSpline = SplineFit(path, smooth=0.1)
		pathU1 = pathSpline.findU(estPose1[:2])
		pathForeU = pathSpline.getUOfDist(originU1, 1.0, distIter = 0.001)
		pathBackU = pathSpline.getUOfDist(originU1, -1.0, distIter = 0.001)
		
		
		" orient the path correctly with the orientation of the node medial axis "

		pathReverse = deepcopy(path)
		pathReverse.reverse()
		pathSplineReverse = SplineFit(pathReverse, smooth=0.1)

		" FIXME: switches between global and local coordinates for medial axis incorrectly "

		overlapMatch = []
		angleSum1 = 0.0
		angleSum2 = 0.0
		for i in range(0,len(globalMedial)):
			p_1 = globalMedial[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(path, p_1)
			if minDist < 0.5:
				overlapMatch.append((i,j,minDist))

				pathU1 = globalMedialSpline1.findU(p_1)	
				pathU2 = pathSpline.findU(p_2)	
				pathU2_R = pathSplineReverse.findU(p_2)	

				pathVec1 = globalMedialSpline1.getUVector(pathU1)
				pathVec2 = pathSpline.getUVector(pathU2)
				pathVec2_R = pathSplineReverse.getUVector(pathU2_R)

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
		print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
		if angleSum1 > angleSum2:
			orientedPath = pathReverse
		else:
			orientedPath = path
			
		orientedPathSpline = SplineFit(orientedPath, smooth=0.1)

		#medialSpline1 = SplineFit(medial1, smooth=0.1)


		globalSamples = orientedPathSpline.getUniformSamples(spacing = 0.04)
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
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(5,6))
		print len(closestPairs), "closest pairs"


		if len(closestPairs) > 0:
			originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
			pathU1 = orientedPathSpline.findU(globalSamples[closestPairs[0][0]])

		else:
		
			originU2 = medialSpline1.findU([0.0,0.0])
			pathU1 = orientedPathSpline.findU(estPose1[:2])
		
		



		resultArgs = self.paths.getDeparturePoint(orientedPath, nodeID, plotIter = True)
		"""
		departurePoint1 = resultArgs[0]
		depAngle1 = resultArgs[1]
		isInterior1 = resultArgs[2]
		isExist1 = resultArgs[3]
		discDist1 = resultArgs[4]
		departurePoint2 = resultArgs[5]
		depAngle2 = resultArgs[6]
		isInterior2 = resultArgs[7]
		isExist2 = resultArgs[8]
		discDist2 = resultArgs[9]
		contigFrac = resultArgs[10]
		overlapSum = resultArgs[11]
		"""
		
		isExist1 = resultArgs[3]
		isExist2 = resultArgs[8]

		originForeU = globalMedialSpline1.getUOfDist(originU1, -0.2, distIter = 0.001)
		originBackU = globalMedialSpline1.getUOfDist(originU1, 0.2, distIter = 0.001)
							
		forePoint = globalMedialSpline1.getU(originForeU)
		backPoint = globalMedialSpline1.getU(originBackU)

		pathForeU = orientedPathSpline.findU(forePoint)
		pathBackU = orientedPathSpline.findU(backPoint)
		
		if pathForeU > pathBackU:
			deltaForeU = 0.02
			deltaBackU = -0.02
		else:
			deltaForeU = -0.02
			deltaBackU = 0.02
				
		" if departing forward, find the first location going backward that it no longer departs "
		if isExist1 and not isExist2:
			
			currPathU = pathU1
			isDone = False

			print "searching pathOverlapUpdate:", isExist1, isExist2, originForeU, originBackU, pathForeU, pathBackU
			
			while not isDone:

				#currPathU += -0.02
				currPathU += deltaBackU

				u2 = originU2
				u1 = currPathU
				angGuess = 0.0

				resultPose, lastCost, matchCount = gen_icp.globalOverlapICP_GPU2([u1,u2,angGuess], orientedPath, medial1, plotIter = True, n1 = nodeID, n2 = -1)
			 	
				self.nodeHash[nodeID].setGPACPose(resultPose)

				resultArgs = self.paths.getDeparturePoint(orientedPath, nodeID, plotIter = True)
				isExist1 = resultArgs[3]
				
				print "pathOverlapUpdate:", currPathU, isExist1, resultPose
				
				if not isExist1:
					isDone = True
				
				" reset to original position "	
				if deltaBackU > 0.0:
					if currPathU > 1.0:
						self.nodeHash[nodeID].setGPACPose(estPose1)
						isDone = True
				else:					
					if currPathU < 0.0:
						self.nodeHash[nodeID].setGPACPose(estPose1)
						isDone = True
			
		
			" if departing backward, find the first location going forward that it no longer departs "
		elif not isExist1 and isExist2:
			currPathU = pathU1
			isDone = False

			print "searching pathOverlapUpdate:", isExist1, isExist2, originForeU, originBackU, pathForeU, pathBackU
			
			while not isDone:

				currPathU += deltaForeU

				u2 = originU2
				u1 = currPathU
				angGuess = 0.0

				resultPose, lastCost, matchCount = gen_icp.globalOverlapICP_GPU2([u1,u2,angGuess], orientedPath, medial1, plotIter = True, n1 = nodeID, n2 = -1)
				
				self.nodeHash[nodeID].setGPACPose(resultPose)

				resultArgs = self.paths.getDeparturePoint(orientedPath, nodeID, plotIter = True)
				isExist2 = resultArgs[8]

				print "pathOverlapUpdate:", currPathU, isExist2, resultPose
				
				if not isExist2:
					isDone = True
				
				" reset to original position "					
				if deltaForeU > 0.0:
					if currPathU > 1.0:
						self.nodeHash[nodeID].setGPACPose(estPose1)
						isDone = True
				else:					
					if currPathU < 0.0:
						self.nodeHash[nodeID].setGPACPose(estPose1)
						isDone = True
			

		
		else:
			print "node", nodeID, "departs BOTH SIDES, not changing anything!"
			return


	def fitToSplices(self, nodeID, splicedPaths1, globalPoint1, globalPoint2):
		
		totalGuesses = []
		
		for k in range(len(splicedPaths1)):

			path = splicedPaths1[k]

			" make the path constraints "								
			guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPoint1, globalPoint2)

			estPose = self.nodeHash[nodeID].getGlobalGPACPose()
			angDiff = abs(diffAngle(estPose[2], guessPose1[2]))
			totalGuesses.append((angDiff, cost1, guessPose1, k))
		
		totalGuesses.sort()
		print "totalGuesses:", totalGuesses

		guessPose = totalGuesses[0][2]
		
		if len(splicedPaths1) > 1:
			angDiff1 = totalGuesses[0][0]
			angDiff2 = totalGuesses[1][0]
			
			print "angDiffs:", angDiff1, angDiff2
			if fabs(angDiff1-angDiff2) < 0.1:
				resultArgs1 = self.paths.getDeparturePoint(splicedPaths1[totalGuesses[0][3]], nodeID, plotIter = True)
				resultArgs2 = self.paths.getDeparturePoint(splicedPaths1[totalGuesses[1][3]], nodeID, plotIter = True)
				
				"select the most contiguous overlap "
				contigFrac1 = resultArgs1[10]
				contigFrac2 = resultArgs2[10]
				
				print "contigFracs:", contigFrac1, contigFrac2
				
				if contigFrac1 < contigFrac2:
					guessPose = totalGuesses[1][2]

		self.nodeHash[nodeID].setGPACPose(guessPose)
		print "set to pose", nodeID, "to", guessPose

	def mergePriorityConstraints(self):
		
		totalConstraints = []
		
		print "merging: do nothing"

		return
		
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
		v_list, merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")

		#self.edgeHash = {}
		#for i in range(len(merged_constraints)):
		#	node1 = merged_constraints[i][0]
		#	node2 = merged_constraints[i][1]
		#	self.edgeHash[(node1, node2)] = [merged_constraints[i][2], merged_constraints[i][3]]
		
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


	def makeGlobalMedialOverlapConstraint(self, nodeID, globalPath, globalJunctionPoint, globalDeparturePoint, fullOverlap = False ):

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
		
		#resultPose, lastCost = gen_icp.globalOverlapICP([u1, u2, angGuess], orientedGlobalPath, hull1, medial1,  plotIter = False, n1 = nodeID, n2 = 0)
		#resultPose2, lastCost2 = gen_icp.globalOverlapTransformICP([u1, u2, angGuess], orientedGlobalPath, hull1, medial1,  plotIter = True, n1 = nodeID, n2 = 0)

		#uSet = [i*0.01 for i in range(100)]
		#globalSpline = SplineFit(orientedGlobalPath, smooth=0.1)
		#medialSpline = SplineFit(medial1, smooth=0.1)
	
		#poses_1 = globalSpline.getUVecSet(uSet)
		#poses_2 = medialSpline.getUVecSet(uSet)


		#resultPose, lastCost = gen_icp.globalOverlapICP([u1,u2,angGuess], orientedGlobalPath, medial1, poses_1, poses_2)

		
		if fullOverlap:
			resultPose, lastCost, matchCount = gen_icp.globalOverlapICP_GPU2([u1,u2,angGuess], orientedGlobalPath, medial1, plotIter = True, n1 = nodeID, n2 = -1, minMatchDist = 2.0)
		else:
			resultPose, lastCost, matchCount = gen_icp.globalOverlapICP_GPU2([u1,u2,angGuess], orientedGlobalPath, medial1, plotIter = True, n1 = nodeID, n2 = -1)
		


		print "estimating branch pose", nodeID, "at",  resultPose[0], resultPose[1], resultPose[2]

		return resultPose, lastCost


	def makeMultiJunctionMedialOverlapConstraint(self, nodeID1, nodeID2, isMove = True, isForward = True, inPlace = False, uRange = 1.5):

		#isMove = False
		#inPlace = True

		def computeOffset(point1, point2, ang1, ang2):
		
			" corner points and orientations "
			corner1Pose = [point1[0], point1[1], ang1]
			corner2Pose = [point2[0], point2[1], ang2]
			
			" convert the desired intersection point on curve 1 into global coordinates "
			poseProfile1 = Pose([0.0,0.0,0.0])
				
			" now convert this point into a pose, and perform the inverse transform using corner2Pose "
			desGlobalPose2 = Pose(corner1Pose)
			
			" perform inverse offset from the destination pose "
			negCurve2Pose = desGlobalPose2.doInverse(corner2Pose)
			
			" relative pose between pose 1 and pose 2 to make corners coincide and same angle "
			resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
			localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
			
			return [localOffset[0], localOffset[1], localOffset[2]]



		node1 = self.nodeHash[nodeID1]
		node2 = self.nodeHash[nodeID2]


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
				

				if isMove:
					
					medialSpline1 = SplineFit(medial1, smooth=0.1)
					medialSpline2 = SplineFit(medial2, smooth=0.1)

					pose1 = medialSpline1.getUVecSet([0.5, 0.5+0.02])[0]	
					pose2 = medialSpline2.getUVecSet([0.5, 0.5+0.02])[0]
				
					point1 = [pose1[0],pose1[1]]
					point2 = [pose2[0],pose2[1]]
					ang1 = pose1[2]
					ang2 = pose2[2]
				
					offset = computeOffset(point1, point2, ang1, ang2)

					points1 = medialSpline1.getUniformSamples()
					points2 = medialSpline2.getUniformSamples()
				
				
					" transform the past poses "
					points2_offset = []
					medialProfile = Pose(offset)
					for p in points2:
						result = medialProfile.convertLocalOffsetToGlobal(p)
						points2_offset.append(result)

					" transform the past poses "
					medial2_offset = []
					for p in medial2:
						result = medialProfile.convertLocalToGlobal(p)
						medial2_offset.append(result)

				
					" if isMove is True, then we align medial2 with medial1 "
					points2R_offset = deepcopy(points2_offset)
					points2R_offset.reverse()
							
					overlapMatch = []
					angleSum1 = 0.0
					angleSum2 = 0.0
					for i in range(0,len(points1)):
						p_1 = points1[i]
						p_2, j, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
						if minDist < 0.5:
							p_2R, j, minDist = gen_icp.findClosestPointInA(points2R_offset, p_2)
							
							ang1 = abs(diffAngle(p_1[2], p_2[2]))
							ang2 = abs(diffAngle(p_1[2], p_2R[2]))
							
							angleSum1 += ang1
							angleSum2 += ang2
					
					" select global path orientation based on which has the smallest angle between tangent vectors "
					print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
					if angleSum1 > angleSum2:
						#orientedMedial2 = deepcopy(medial2Reverse)
						medial2Reverse = deepcopy(medial2)
						medial2Reverse.reverse()
					else:
						orientedMedial2 = deepcopy(medial2)

					medialSpline1 = SplineFit(medial1, smooth=0.1)
					medialSpline2 = SplineFit(orientedMedial2, smooth=0.1)
			
					originU1 = medialSpline1.findU([0.0,0.0])	
					originU2 = medialSpline2.findU([0.0,0.0])	
					
					if isForward:
						
						moveDist = 0.3
							
						if len(node2.frontProbeError) > 0:
			
							frontSum = 0.0
							frontProbeError = node2.frontProbeError
							for n in frontProbeError:
								frontSum += n
							foreAvg = frontSum / len(frontProbeError)
											
							if foreAvg >= 1.4:
								u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
							else:	
								u2 = medialSpline2.getUOfDist(originU2, moveDist, distIter = 0.001)
						else:	
							u2 = medialSpline2.getUOfDist(originU2, moveDist, distIter = 0.001)
							
					else:
						moveDist = -0.3
						u2 = medialSpline2.getUOfDist(originU2, moveDist, distIter = 0.001)
					
					print "computed u2 =", u2, "from originU2 =", originU2

					u1 = originU1
					angGuess = 0.0
				
					" create the ground constraints "
					gndGPAC1Pose = node1.getGndGlobalGPACPose()
					currProfile = Pose(gndGPAC1Pose)
					gndGPAC2Pose = node2.getGndGlobalGPACPose()
					gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
					
					#result, hist = gen_icp.overlapICP(estPose1, diffOffset, [u1, u2, angGuess], hull1, hull2, orientedMedial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (i,k), n2 = (j,l), uRange = uRange)
					#result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)
					
					result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)
						
					
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
			
			
					medialError, matchCount = self.computeMedialError(nodeID1, nodeID2, offset, minMatchDist = 0.5, tail1=k, tail2=l)
	
			
					results.append((hist[0], angleSum, medialError, matchCount, angDiff, k, l, transform, covE, hist, matchCount, overlapSum, angleSum))


				else:
					
					medialSpline1 = SplineFit(medial1, smooth=0.1)
					medialSpline2 = SplineFit(medial2, smooth=0.1)

					pose1 = medialSpline1.getUVecSet([0.5, 0.5+0.02])[0]	
					pose2 = medialSpline2.getUVecSet([0.5, 0.5+0.02])[0]
				
					point1 = [pose1[0],pose1[1]]
					point2 = [pose2[0],pose2[1]]
					ang1 = pose1[2]
					ang2 = pose2[2]
				
					offset = computeOffset(point1, point2, ang1, ang2)

					points1 = medialSpline1.getUniformSamples()
					points2 = medialSpline2.getUniformSamples()
				
				
					" transform the past poses "
					points2_offset = []
					medialProfile = Pose(offset)
					for p in points2:
						result = medialProfile.convertLocalOffsetToGlobal(p)
						points2_offset.append(result)

					" transform the past poses "
					medial2_offset = []
					for p in medial2:
						result = medialProfile.convertLocalToGlobal(p)
						medial2_offset.append(result)

				
					" if isMove is True, then we align medial2 with medial1 "
					points2R_offset = deepcopy(points2_offset)
					points2R_offset.reverse()
							
					overlapMatch = []
					angleSum1 = 0.0
					angleSum2 = 0.0
					for i in range(0,len(points1)):
						p_1 = points1[i]
						p_2, j, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
						if minDist < 0.5:
							p_2R, j, minDist = gen_icp.findClosestPointInA(points2R_offset, p_2)
							
							ang1 = abs(diffAngle(p_1[2], p_2[2]))
							ang2 = abs(diffAngle(p_1[2], p_2R[2]))
							
							angleSum1 += ang1
							angleSum2 += ang2
					
					" select global path orientation based on which has the smallest angle between tangent vectors "
					print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
					if angleSum1 > angleSum2:
						#orientedMedial2 = deepcopy(medial2Reverse)
						medial2Reverse = deepcopy(medial2)
						medial2Reverse.reverse()
					else:
						orientedMedial2 = deepcopy(medial2)

	
					#medialSpline1 = SplineFit(medial1, smooth=0.1)
					medialSpline2 = SplineFit(orientedMedial2, smooth=0.1)
			
					originU1 = medialSpline1.findU([0.0,0.0])	
					originU2 = medialSpline2.findU([0.0,0.0])	
			
					if inPlace:
						originU1 = 0.6
						originU2 = 0.4
						u2 = originU2
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
					result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)
					
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
			
			
					medialError, matchCount = self.computeMedialError(nodeID1, nodeID2, offset, minMatchDist = 0.5, tail1=k, tail2=l)
	
			
					results.append((hist[0], angleSum, medialError, matchCount, angDiff, k, l, transform, covE, hist, matchCount, overlapSum, angleSum))
		
		results.sort(reverse=True)
		#results.sort(reverse=False)
		
		print "Multi-Junction Overlap of", nodeID1, "and", nodeID2
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


		#result = gen_icp.overlapICP(estPose1, [u1, u2, angGuess], hull1, hull2, medial1, medial2, node1.rootPose, node2.rootPose, plotIter = True, n1 = i, n2 = j)
		#result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)
		result, hist = gen_icp.overlapICP_GPU2(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)
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

	def getNearestPathPoint(self, originPoint):

		newPoint = self.paths.getNearestPathPoint(originPoint)

		return newPoint

	def computeNavigationPath(self, startPose, endPose):
		return self.paths.computeNavigationPath(startPose, endPose)

	def addTargetBatchConstraints(self, targetNodeID, pathNodes):

		allNodes = pathNodes
		
		#paths = []
		#for i in range(self.numNodes):
		#	paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
			#print "paths", i, ":", paths[i]




		constraintPairs = []
		for nodeID in allNodes:
			constraintPairs.append([targetNodeID,nodeID])


		medials = {}
		splines = {}
		originUs = {}
		originForeUs = {}
		originBackUs = {}
		for k in allNodes + [targetNodeID]:
	
			node1 = self.nodeHash[k]
			hull1, medial1 = computeHullAxis(k, node1, tailCutOff = False)
			medials[k] = medial1
			#estPose1 = node1.getGlobalGPACPose()		
			medialSpline1 = SplineFit(medial1, smooth=0.1)
			splines[k] = medialSpline1
			
			originU1 = medialSpline1.findU([0.0,0.0])	
			originUs[k] = originU1
			originForeUs[k] = splines[k].getUOfDist(originUs[k], 1.0, distIter = 0.001)
			originBackUs[k] = splines[k].getUOfDist(originUs[k], -1.0, distIter = 0.001)


		targetNode = self.nodeHash[targetNodeID]
		targetPose = targetNode.getGlobalGPACPose()
		
		isTargFeatureless = targetNode.getIsFeatureless()
		
		candidates = []
		
		u2_1 = splines[targetNodeID].getUOfDist(originUs[targetNodeID], 0.0, distIter = 0.001)
		u2_2 = splines[targetNodeID].getUOfDist(originUs[targetNodeID], 1.0, distIter = 0.001)
		u2_3 = splines[targetNodeID].getUOfDist(originUs[targetNodeID], -1.0, distIter = 0.001)
		for k in allNodes:

			candNode = self.nodeHash[k]
			candPose = candNode.getGlobalGPACPose()
			
			cartDist1 = sqrt((candPose[0]-targetPose[0])**2 + (candPose[1]-targetPose[1])**2)
			angDiff1 = diffAngle(candPose[2],targetPose[2])
			
			isFeatureless = candNode.getIsFeatureless()
			
			if k == targetNodeID-1 or k == targetNodeID-2:
				isNeigh = True
			else:
				isNeigh = False
			
			pathIDs = []
			allPaths = self.paths.getPathIDs()
			for j in allPaths:
				if self.paths.isNodeExist(k,j):
					pathIDs.append(j)


			"  select closest point to origin of candidate pose as u2 "
			estPose1 = candPose
			estPose2 = targetPose
			poseOrigin = Pose(estPose1)
			offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
			
			points2 = splines[k].getUniformSamples()
			p_1 = splines[targetNodeID].getU(originUs[targetNodeID])
			
			points2_offset = []
			for p in points2:
				result = gen_icp.dispOffset(p, offset)		
				points2_offset.append(result)
	
			p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
	
			u2 = splines[k].findU(points2[i_2])	
			
			args = [k, targetNodeID, medials[k], medials[targetNodeID], originUs[k], u2, 0.0]
			candidates.append({ "nodeID1" : k, "nodeID2" : targetNodeID,
							"medial1": medials[k], "medial2": medials[targetNodeID],
							"u1": args[4], "u2": args[5], "initAng": args[6],
							"isFeatureless1": isFeatureless, "isFeatureless2": isTargFeatureless,
							"isNeigh": isNeigh, "cartDist1": cartDist1, "angDiff1": angDiff1,
							"pathIDs": pathIDs, "isAlternate": False})

		PATH_MATCH_DIST = 1.0
		PAIR_ANG_DELTA = 0.3


		candidates2 = []
		argSets = []
		print
		print "CANDIDATES:"
		for cand in candidates:
			nodeID1 = cand['nodeID1']
			nodeID2 = cand['nodeID2']
			medial1 = cand['medial1']
			medial2 = cand['medial2']
			u1 = cand['u1']
			u2 = cand['u2']
			initAng = cand['initAng']
			angDiff1 = cand['angDiff1']
			cartDist1 = cand['cartDist1']
			isFeatureless = cand['isFeatureless1']
			isNeigh = cand['isNeigh']
			pathIDs = cand['pathIDs']

			isInPath = False
			commonPathID = -1
			for j in pathIDs:
				if self.paths.isNodeExist(nodeID2,j):
					isInPath = True
					commonPathID = j

			oldConstraints = self.getEdges(nodeID1, nodeID2) + self.getEdges(nodeID2, nodeID1)
			isExist = False
			for constraint in oldConstraints:
				isExist = True

			print nodeID1, nodeID2, u1, u2, initAng, pathIDs
			
			
			if True and not isNeigh:
				if cartDist1 < PATH_MATCH_DIST and fabs(diffAngle(angDiff1,initAng)) < PAIR_ANG_DELTA:
					if (isTargFeatureless and isFeatureless) or (not isTargFeatureless and not isFeatureless):
						#if isInPath and not isExist:
						if not isExist:

							points1 = splines[nodeID1].getUniformSamples()
							points2 = splines[nodeID2].getUniformSamples()

							frontDepI1 = 0
							backDepI1 = len(points1)-1
							frontDepI2 = 0
							backDepI2 = len(points2)-1
							args = [nodeID1, nodeID2, medial1, medial2, u1, u2, initAng, frontDepI1, backDepI1, frontDepI2, backDepI2, False]
							argSets.append(args)
							candidates2.append(cand)							
							
							
							
			if False and fabs(diffAngle(angDiff1,initAng)) < PAIR_ANG_DELTA:
				args = [nodeID1, nodeID2, medial1, medial2, u1, u2, initAng]
				argSets.append(args)
				candidates2.append(cand)


		print "FILTERED:"
		for cand in candidates2:
			nodeID1 = cand['nodeID1']
			nodeID2 = cand['nodeID2']
			u1 = cand['u1']
			u2 = cand['u2']
			initAng = cand['initAng']
			print nodeID1, nodeID2, u1, u2, initAng

		t1 = time.time()
		results = gen_icp.serialOverlapICP(argSets)
		#results = gen_icp.batchOverlapICP(argSets)
		t2 = time.time()
		print "batchOverlapICP:", t2-t1, "seconds", len(argSets), "pairs"

		"""
		METRICS:
		
		1) initial angle difference
		2) constraint angle difference
		3) initial cartesian difference
		4) constraint cartestian difference
		5) featurelessness state
		6) paired inplace or not status
		7) path classification
		8) ICP histogram
		"""

		for k in range(len(candidates2)):
			cand = candidates2[k]
			nodeID1 = cand['nodeID1']
			nodeID2 = cand['nodeID2']
			medial1 = cand['medial1']
			medial2 = cand['medial2']
			u1 = cand['u1']
			u2 = cand['u2']
			initAng = cand['initAng']
			angDiff1 = cand['angDiff1']
			cartDist1 = cand['cartDist1']
			isFeatureless = cand['isFeatureless1']
			isNeigh = cand['isNeigh']
			pathIDs = cand['pathIDs']
						
			
			res = results[k]
			offset = res[0]
			histogram = res[1]
						
			transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
			covE =  self.E_overlap
			
			candNode = self.nodeHash[nodeID1]
			candPose = candNode.getGlobalGPACPose()
			
			cartDist2 = sqrt((candPose[0]-targetPose[0]-offset[0])**2 + (candPose[0]-targetPose[0]-offset[0])**2)
			angDiff2 = diffAngle(offset[2], angDiff1)

			print "cartDist:", candPose[0]-targetPose[0], offset[0], candPose[1]-targetPose[1], offset[1], angDiff1, angDiff2

			cand["offset"] = offset
			cand["transform"] = transform
			cand["histogram"] = histogram
			cand["angDiff2"] = angDiff2
			cand["cartDist2"] = cartDist2

			if not isFeatureless:
				cand["covE"] = deepcopy(self.E_junction)
			else:
				cand["covE"] = deepcopy(self.E_overlap)

			cand["xiError"] = 0.0
			cand["xiErrorBase"] = 0.0


		candidates2.sort(key=lambda cand: cand["cartDist2"])

		candidates3 = []
		isPicked = {}
		for cand in candidates2:
			
			nodeID1 = cand["nodeID1"]
			nodeID2 = cand["nodeID2"]
			try:
				isPicked[(nodeID1,nodeID2)]
			except:
				isPicked[(nodeID1,nodeID2)] = True
				candidates3.append(cand)
			else:
				pass
				

 
		print "RESULTS:"
		finalCandidates = []
		for cand in candidates3:
			nodeID1 = cand["nodeID1"]
			nodeID2 = cand["nodeID2"]
			transform = cand["transform"]
			covE = cand["covE"]
			offset = cand["offset"]

			angDiff2 = cand["angDiff2"]
			cartDist2 = cand["cartDist2"]
			
			xiError = cand["xiError"]
			xiBaseError = cand["xiErrorBase"]
						
			#print nodeID1, nodeID2, offset, cartDist2, angDiff2, xiError, xiError-xiBaseError
			print nodeID1, nodeID2, cartDist2, angDiff2, xiError, xiError-xiBaseError
			
			self.allPathConstraints.append([nodeID1,nodeID2,transform,covE])
			finalCandidates.append([nodeID1,nodeID2,transform,covE])

		print

		" draw new hypotheses "
		#for cand in candidates3:
		#	self.drawCandidate(cand)
		
		return finalCandidates[:1]
		
		
		"""
		METRICS:
		
		1) initial angle difference
		2) constraint angle difference
		3) initial cartesian difference
		4) constraint cartestian difference
		5) featurelessness state
		6) paired inplace or not status
		7) path classification
		8) ICP histogram
		"""
			

		"""
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
		"""


		#return transform, covE, hist		

		# result, hist = gen_icp.overlapICP_GPU2(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)
		

		#for pair in constraintPairs:
		#	n1 = pair[0]
		#	n2 = pair[1]
		#	transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)				
	
	
	def addBatchConstraints(self, targetNodeID):

		allNodes = range(targetNodeID)

		#paths = []
		#for i in range(self.numNodes):
		#	paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		constraintPairs = []
		for nodeID in allNodes:
			constraintPairs.append([targetNodeID,nodeID])

		medials = []
		splines = []
		originUs = []
		originForeUs = []
		originBackUs = []
		for k in range(targetNodeID+1):
	
			node1 = self.nodeHash[k]
			hull1, medial1 = computeHullAxis(k, node1, tailCutOff = False)
			medials.append(medial1)
			#estPose1 = node1.getGlobalGPACPose()		
			medialSpline1 = SplineFit(medial1, smooth=0.1)
			splines.append(medialSpline1)
			
			originU1 = medialSpline1.findU([0.0,0.0])	
			originUs.append(originU1)
			originForeUs.append(splines[k].getUOfDist(originUs[k], 1.0, distIter = 0.001))
			originBackUs.append(splines[k].getUOfDist(originUs[k], -1.0, distIter = 0.001))


		targetNode = self.nodeHash[targetNodeID]
		targetPose = targetNode.getGlobalGPACPose()
		
		isTargFeatureless = targetNode.getIsFeatureless()
		
		candidates = []
		
		u2_1 = splines[targetNodeID].getUOfDist(originUs[targetNodeID], 0.0, distIter = 0.001)
		u2_2 = splines[targetNodeID].getUOfDist(originUs[targetNodeID], 1.0, distIter = 0.001)
		u2_3 = splines[targetNodeID].getUOfDist(originUs[targetNodeID], -1.0, distIter = 0.001)
		for k in range(targetNodeID-2):

			candNode = self.nodeHash[k]
			candPose = candNode.getGlobalGPACPose()
			
			cartDist1 = sqrt((candPose[0]-targetPose[0])**2 + (candPose[1]-targetPose[1])**2)
			angDiff1 = diffAngle(candPose[2],targetPose[2])
			
			isFeatureless = candNode.getIsFeatureless()
			
			if k == targetNodeID-1 or k == targetNodeID-2:
				isNeigh = True
			else:
				isNeigh = False
			
			pathIDs = []
			allPaths = self.paths.getPathIDs()
			for j in allPaths:
				if self.paths.isNodeExist(k,j):
					pathIDs.append(j)


			"  select closest point to origin of candidate pose as u2 "
			estPose1 = candPose
			estPose2 = targetPose
			poseOrigin = Pose(estPose1)
			offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
			
			points2 = splines[k].getUniformSamples()
			p_1 = splines[targetNodeID].getU(originUs[targetNodeID])
			
			points2_offset = []
			for p in points2:
				result = gen_icp.dispOffset(p, offset)		
				points2_offset.append(result)
	
			p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
	
			u2 = splines[k].findU(points2[i_2])	
			
			

			
			
			args = [k, targetNodeID, medials[k], medials[targetNodeID], originUs[k], u2, 0.0]
			candidates.append({ "nodeID1" : k, "nodeID2" : targetNodeID,
							"medial1": medials[k], "medial2": medials[targetNodeID],
							"u1": args[4], "u2": args[5], "initAng": args[6],
							"isFeatureless1": isFeatureless, "isFeatureless2": isTargFeatureless,
							"isNeigh": isNeigh, "cartDist1": cartDist1, "angDiff1": angDiff1,
							"pathIDs": pathIDs})

		allSplices, terminals, junctions = self.paths.getAllSplices()

		if False:
			print "foo"
			print "foo"


		PATH_MATCH_DIST = 1.0
		PAIR_ANG_DELTA = 0.3


		candidates2 = []
		argSets = []
		icpArgs = []
		print
		print "CANDIDATES:"
		for cand in candidates:
			nodeID1 = cand['nodeID1']
			nodeID2 = cand['nodeID2']
			medial1 = cand['medial1']
			medial2 = cand['medial2']
			u1 = cand['u1']
			u2 = cand['u2']
			initAng = cand['initAng']
			#transform = cand[2]
			#covE = cand[3]
			angDiff1 = cand['angDiff1']
			#angDiff2 = cand[5]
			cartDist1 = cand['cartDist1']
			#cartDist2 = cand[7]
			isFeatureless = cand['isFeatureless1']
			isNeigh = cand['isNeigh']
			pathIDs = cand['pathIDs']
			#histogram = cand[11]

			#isInPath = False
			#for j in pathIDs:
			#	if self.nodeSets[j].count(nodeID1) > 0:
			#		isInPath = True

			isInPath = False
			isAlternate = False
			#allPaths = self.paths.getPathIDs()
			commonPathID = -1
			for j in pathIDs:
				if self.paths.isNodeExist(nodeID2,j):
					isInPath = True
					commonPathID = j

			"  ----------------------- "
			if isInPath:
				path1 = self.paths.getPath(commonPathID)
				parentID1 = path1["parentID"]
				if parentID1 != None:
				
					sPaths = allSplices[(parentID1, commonPathID)]
					termPath1 = sPaths[0]['termPath']
					termPath2 = sPaths[1]['termPath']
					
					startKey = termPath1[0]
					endKey = termPath1[-1]
		
					shortestPathSpanTree, shortestDist = self.paths.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					
					path1 = []
					while currNode != endKey:
						path1.append(self.paths.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path1.append(self.paths.pathGraph.get_node_attributes(currNode))
		
					startKey = termPath2[0]
					endKey = termPath2[-1]
		
					shortestPathSpanTree, shortestDist = self.paths.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					
					path2 = []
					while currNode != endKey:
						path2.append(self.paths.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path2.append(self.paths.pathGraph.get_node_attributes(currNode))
					
					resultSum1_1 = self.paths.getOverlapCondition2(path1, nodeID1)
					resultSum1_2 = self.paths.getOverlapCondition2(path2, nodeID1)
					resultSum2_1 = self.paths.getOverlapCondition2(path1, nodeID2)
					resultSum2_2 = self.paths.getOverlapCondition2(path2, nodeID2)
					
					print "intra-path overlaps:[", parentID1, ",", commonPathID, "]"
					print nodeID1, resultSum1_1, resultSum1_2
					print nodeID2, resultSum2_1, resultSum2_2
					
					if resultSum1_2 == 0:
						ratio1 = 1e100
					else:
						ratio1 = float(resultSum1_1) / resultSum1_2
					
					if resultSum2_1 == 0:
						ratio2 = 1e100
					else:
						ratio2 = float(resultSum2_2) / resultSum2_1

					if ratio1 > 2.0 and ratio2 > 2.0:
						print "alternate junction state!"
						isAlternate = True

					if resultSum1_1 == 0:
						ratio1 = 1e100
					else:
						ratio1 = float(resultSum1_2) / resultSum1_1
					
					if resultSum2_2 == 0:
						ratio2 = 1e100
					else:
						ratio2 = float(resultSum2_1) / resultSum2_2
					
					if ratio1 > 2.0 and ratio2 > 2.0:
						print "alternate junction state!"
						isAlternate = True
						
				elif False:
					sPaths = allSplices[(commonPathID, commonPathID)]
					termPath1 = sPaths[0]['termPath']
					
					startKey = termPath1[0]
					endKey = termPath1[-1]
		
					shortestPathSpanTree, shortestDist = self.paths.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					
					path1 = []
					while currNode != endKey:
						path1.append(self.paths.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path1.append(self.paths.pathGraph.get_node_attributes(currNode))
							
					resultSum1_1 = self.paths.getOverlapCondition2(path1, nodeID1)
					resultSum2_1 = self.paths.getOverlapCondition2(path1, nodeID2)
					
					print "intra-path overlaps:[", commonPathID, ",", commonPathID, "]"
					print nodeID1, resultSum1_1
					print nodeID2, resultSum2_1

		
			"  ----------------------- "


			cand['isAlternate'] = isAlternate

			oldConstraints = self.getEdges(nodeID1, nodeID2) + self.getEdges(nodeID2, nodeID1)
			isExist = False
			for constraint in oldConstraints:
				isExist = True

			print nodeID1, nodeID2, u1, u2, initAng, cartDist1, angDiff1, pathIDs, isAlternate
			

			#sum1 = self.paths.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
			
			
			
			
			if True and not isNeigh:
				if cartDist1 < PATH_MATCH_DIST and fabs(diffAngle(angDiff1,initAng)) < PAIR_ANG_DELTA or isAlternate:
					if (isTargFeatureless and isFeatureless) or (not isTargFeatureless and not isFeatureless):
						if isInPath and not isExist:

							if isAlternate:
								resultArgs1 = self.paths.getDeparturePoint(self.trimmedPaths[commonPathID], nodeID1)
								resultArgs2 = self.paths.getDeparturePoint(self.trimmedPaths[commonPathID], nodeID2)							
				
								frontDeparturePoint1 = resultArgs1[0]
								backDeparturePoint1 = resultArgs1[5]
								frontDeparturePoint2 = resultArgs2[0]
								backDeparturePoint2 = resultArgs2[5]
	
								#frontDeparturePoint1 = 0
								#backDeparturePoint1 = 0
								#frontDeparturePoint2 = 0
								#backDeparturePoint2 = 0
	
								node1 = self.nodeHash[nodeID1]
								estPose1 = node1.getGlobalGPACPose()
								node2 = self.nodeHash[nodeID2]
								estPose2 = node2.getGlobalGPACPose()
	
								points1 = splines[nodeID1].getUniformSamples()
								points2 = splines[nodeID2].getUniformSamples()
	
								points1_offset = []
								for p in points1:
									result = gen_icp.dispOffset(p, estPose1)		
									points1_offset.append(result)
	
								points2_offset = []
								for p in points2:
									result = gen_icp.dispOffset(p, estPose2)		
									points2_offset.append(result)
									
	
	
								if frontDeparturePoint1 != 0:
									p_1, depI, pDist1 = gen_icp.findClosestPointInA(points1_offset, frontDeparturePoint1)
									p_1, frontDepI1, pDist1 = gen_icp.findClosestPointInA(medial1, points1[depI])
								else:
									frontDepI1 = 0
	
								if backDeparturePoint1 != 0:
									p_1, depI, pDist1 = gen_icp.findClosestPointInA(points1_offset, backDeparturePoint1)
									p_1, backDepI1, pDist1 = gen_icp.findClosestPointInA(medial1, points1[depI])
								else:
									backDepI1 = len(points1)-1
	
								if frontDeparturePoint2 != 0:
									p_2, depI, pDist2 = gen_icp.findClosestPointInA(points2_offset, frontDeparturePoint2)
									p_2, frontDepI2, pDist2 = gen_icp.findClosestPointInA(medial2, points2[depI])
								else:
									frontDepI2 = 0
	
								if backDeparturePoint2 != 0:
									p_2, depI, pDist2 = gen_icp.findClosestPointInA(points2_offset, backDeparturePoint2)
									p_2, backDepI2, pDist2 = gen_icp.findClosestPointInA(medial2, points2[depI])
								else:
									backDepI2 = len(points2)-1
								
								print nodeID1, nodeID2, "indices:", frontDepI1, backDepI1, frontDepI2, backDepI2
								
								" prevent case where little or no overlap of parent path prevents constraining "
								if abs(frontDepI1-backDepI1) > 5 and abs(frontDepI2-backDepI2) > 5:
									icpArgs = [nodeID1, nodeID2, medial1, medial2, u1, u2, initAng, frontDepI1, backDepI1, frontDepI2, backDepI2, isAlternate]
							
							else:
								points1 = splines[nodeID1].getUniformSamples()
								points2 = splines[nodeID2].getUniformSamples()

								frontDepI1 = 0
								backDepI1 = len(points1)-1
								frontDepI2 = 0
								backDepI2 = len(points2)-1

								icpArgs = [nodeID1, nodeID2, medial1, medial2, u1, u2, initAng, frontDepI1, backDepI1, frontDepI2, backDepI2, isAlternate]
							
							if len(icpArgs) > 0:
								argSets.append(icpArgs)
								candidates2.append(cand)
								icpArgs = []
								
			if False and fabs(diffAngle(angDiff1,initAng)) < PAIR_ANG_DELTA:
				icpArgs = [nodeID1, nodeID2, medial1, medial2, u1, u2, initAng]
				argSets.append(icpArgs)
				candidates2.append(cand)


		print "FILTERED:"
		for cand in candidates2:
			nodeID1 = cand['nodeID1']
			nodeID2 = cand['nodeID2']
			u1 = cand['u1']
			u2 = cand['u2']
			initAng = cand['initAng']
			print nodeID1, nodeID2, u1, u2, initAng
			
		print "ARGS:"
		for icpArgs in argSets:
			print icpArgs

		t1 = time.time()
		results = gen_icp.serialOverlapICP(argSets)
		#results = gen_icp.batchOverlapICP(argSets)
		t2 = time.time()
		print "batchOverlapICP:", t2-t1, "seconds", len(argSets), "pairs"

		"""
		METRICS:
		
		1) initial angle difference
		2) constraint angle difference
		3) initial cartesian difference
		4) constraint cartestian difference
		5) featurelessness state
		6) paired inplace or not status
		7) path classification
		8) ICP histogram
		"""


		#self.deleteAllPriority(CORNER_PRIORITY)
		#self.mergePriorityConstraints()
		#err, edgeTotal = self.getXiError()
		#errAvg = err/edgeTotal
		#baseLineError = err
		
		#print "PATH RESULT:"		
		for k in range(len(candidates2)):
			cand = candidates2[k]
			nodeID1 = cand['nodeID1']
			nodeID2 = cand['nodeID2']
			medial1 = cand['medial1']
			medial2 = cand['medial2']
			u1 = cand['u1']
			u2 = cand['u2']
			initAng = cand['initAng']
			#transform = cand[2]
			#covE = cand[3]
			angDiff1 = cand['angDiff1']
			#angDiff2 = cand[5]
			cartDist1 = cand['cartDist1']
			#cartDist2 = cand[7]
			isFeatureless = cand['isFeatureless1']
			isNeigh = cand['isNeigh']
			pathIDs = cand['pathIDs']
						
			
			res = results[k]
			offset = res[0]
			histogram = res[1]
						
			transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
			covE =  self.E_overlap
			
			candNode = self.nodeHash[nodeID1]
			candPose = candNode.getGlobalGPACPose()
			
			cartDist2 = sqrt((candPose[0]-targetPose[0]-offset[0])**2 + (candPose[0]-targetPose[0]-offset[0])**2)
			#cartDist2 = sqrt(offset[0]**2 + offset[1]**2)
			angDiff2 = diffAngle(offset[2], angDiff1)

			print "cartDist:", cartDist2, angDiff2, candPose[0]-targetPose[0], offset[0], candPose[1]-targetPose[1], offset[1], angDiff1

			cand["offset"] = offset
			cand["transform"] = transform
			cand["histogram"] = histogram
			cand["angDiff2"] = angDiff2
			cand["cartDist2"] = cartDist2

			if not isFeatureless:
				cand["covE"] = deepcopy(self.E_junction)
			else:
				cand["covE"] = deepcopy(self.E_overlap)

			cand["xiError"] = 0.0
			cand["xiErrorBase"] = 0.0
				
			"""
			self.deleteAllPriority(CORNER_PRIORITY)
			self.addPriorityEdge([nodeID1,nodeID2,transform,cand["covE"]], CORNER_PRIORITY)
			self.mergePriorityConstraints()
			err, edgeTotal = self.getXiError()
			errAvg = err/edgeTotal
			
			cand["xiError"] = err
			cand["xiErrorAvg"] = errAvg
			cand["xiErrorBase"] = baseLineError
			"""

		
			#candidates.append([nodeID1, nodeID2, transform, deepcopy(self.E_junction), angDiff1, angDiff2, cartDist1, cartDist2, isFeatureless, isNeigh, pathIDs, histogram])


		candidates2.sort(key=lambda cand: cand["xiError"])

		candidates3 = []
		isPicked = {}
		for cand in candidates2:
			
			nodeID1 = cand["nodeID1"]
			nodeID2 = cand["nodeID2"]
			try:
				isPicked[(nodeID1,nodeID2)]
			except:
				isPicked[(nodeID1,nodeID2)] = True
				candidates3.append(cand)
			else:
				pass
				


		print "RESULTS:"
		finalCandidates = []
		for cand in candidates3:
			nodeID1 = cand["nodeID1"]
			nodeID2 = cand["nodeID2"]
			transform = cand["transform"]
			covE = cand["covE"]
			offset = cand["offset"]

			angDiff2 = cand["angDiff2"]
			cartDist2 = cand["cartDist2"]
			
			xiError = cand["xiError"]
			xiBaseError = cand["xiErrorBase"]

			isAlternate = cand['isAlternate']
			
			if cartDist2 < 1.0 or isAlternate:
	
				#print nodeID1, nodeID2, offset, cartDist2, angDiff2, xiError, xiError-xiBaseError
				print nodeID1, nodeID2, cartDist2, angDiff2, xiError, xiError-xiBaseError
				
				
				self.allPathConstraints.append([nodeID1,nodeID2,transform,covE])
				finalCandidates.append([nodeID1,nodeID2,transform,covE])

		print

		" draw new hypotheses "
		#for cand in candidates3:
		#	self.drawCandidate(cand)
		
		return finalCandidates
		
		
		"""
		METRICS:
		
		1) initial angle difference
		2) constraint angle difference
		3) initial cartesian difference
		4) constraint cartestian difference
		5) featurelessness state
		6) paired inplace or not status
		7) path classification
		8) ICP histogram
		"""
			

		"""
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
		"""


		#return transform, covE, hist		

		# result, hist = gen_icp.overlapICP_GPU2(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)
		

		#for pair in constraintPairs:
		#	n1 = pair[0]
		#	n2 = pair[1]
		#	transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)				



	def computePathConsistency(self, pathID):
	
		nodeSet = self.paths.getNodes(pathID)
		nodeSet.sort()

		departures1 = []
		interiors1 = []
		depPoints1 = []
		distances1 = []
		depAngles1 = []
		isDep = []
		contig1 = []
		sums = []
		counts = []
			
		for nodeID in nodeSet:

			departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID, plotIter = False)
			departures1.append([isExist1,isExist2])
			interiors1.append([isInterior1, isInterior2])
			depPoints1.append([departurePoint1, departurePoint2])
			distances1.append([discDist1, discDist2])
			depAngles1.append([depAngle1, depAngle2])
			contig1.append((contigFrac, overlapSum))
			isDep.append(isExist1 and isInterior1 or isExist2 and isInterior2)

			sum1 = self.paths.getOverlapCondition(self.trimmedPaths[pathID], nodeID, plotIter = False)
			sums.append(sum1)
		
			count = self.paths.getOverlapCondition2(self.trimmedPaths[pathID], nodeID, plotIter = False)
			counts.append(count)

		print "pathID", pathID, "consistencyCheck for nodeSet:", nodeSet
		print "nodeID, foreDep, backDep, foreDist, backDist, contigFrac, overlapSum, overlapSum, matchCount"
		for k in range(len(nodeSet)):
			print nodeSet[k], isDep[k], departures1[k][0] and interiors1[k][0], departures1[k][1] and interiors1[k][1], distances1[k][0], distances1[k][1], contig1[k][0], contig1[k][1], sums[k], counts[k]
		
		if isDep.count(True) == 0:
			return True
		else:
			return False
		#print pathID, "nodeSet:", nodeSet
		#print "overlapMatchCounts:", overlapSums


	def addPathConstraints3(self, pathNodes, targetNodeID, insertNode = False, orderedPathIDs = []):

		print "addPathConstraints:"
		print "pathNodes:", pathNodes
		print "targetNodeID:", targetNodeID
		print "orderedPathIDs:", orderedPathIDs
		
		newConstraints = self.addBatchConstraints(targetNodeID)
						
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()

		isConsistent = True
		for pathID in orderedPathIDs:
			resultVal = self.computePathConsistency(pathID)
			isConsistent &= resultVal

		if not isConsistent:
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.delPriorityEdge( n1, n2, CORNER_PRIORITY)					
				
		"""
		" do the new constraints cause a branching?  Then they suck!  Revert. "
		if len(orderedPathIDs) == 1:
			print "checking new constraint on-path consistency"
			departures1 = []
			interiors1 = []
			depPoints1 = []
			distances1 = []
			depAngles1 = []
			contig1 = []
			
			sums = []
			for pathID in orderedPathIDs:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], targetNodeID)
				departures1.append([isExist1,isExist2])
				interiors1.append([isInterior1, isInterior2])
				depPoints1.append([departurePoint1, departurePoint2])
				distances1.append([discDist1, discDist2])
				depAngles1.append([depAngle1, depAngle2])
				contig1.append((contigFrac, overlapSum))

				sum1 = self.paths.getOverlapCondition(self.trimmedPaths[pathID], targetNodeID, plotIter = True)
				sums.append(sum1)

			
			print "node departures", targetNodeID, ":", departures1
			print "node  interiors", targetNodeID, ":", interiors1
			print "node contiguity", targetNodeID, ":", contig1
			print "node overlap sums", targetNodeID, ":", sums

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
			
			if contig1[0][0] < 0.4 or contig1[-1][0] < 0.4:
				adjustPose1 = True
			else:
				adjustPose1 = False
				
			print "mergeTest results:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1, contig1[0][0], contig1[-1][0], sums[0], sums[-1]
			

			if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or adjustPose1:
				for const in newConstraints:
					n1 = const[0]
					n2 = const[1]
					transform = const[2]
					covE = const[3]
		
					self.delPriorityEdge( n1, n2, CORNER_PRIORITY)			
		"""
		" merge the constraints "
		self.mergePriorityConstraints()
		
		return
		

		

		

	def addPathConstraints2(self, pathNodes, targetNodeID, insertNode = False):

		print "addPathConstraints:"
		print "pathNodes:", pathNodes
		print "targetNodeID:", targetNodeID

		pathNodes = []
		for k in range(self.numPaths):
			pathNodes += self.nodeSets[k]

		print "pathNodes:", pathNodes
		
		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
				
		self.addBatchConstraints(targetNodeID)

		" analyze over all possible constraints "
		totalHypotheses = self.allPathConstraints

		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
		
		self.mergePriorityConstraints()
		
		" recompute dijkstra projection "
		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
			#print "paths", i, ":", paths[i]
		
		
		
		#totalHypotheses = constraintResults
		
		estPaths = []
		
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

				#Th1 = paths[m1][m2][0]
				#Ch1 = paths[m1][m2][1]
				
				#Th2 = paths[n1][n2][0]
				#Ch2 = paths[n1][n2][1]
				
				#poseOrigin = Pose(self.nodeHash[m1].getGlobalGPACPose())
				#offset = poseOrigin.convertGlobalPoseToLocal(self.nodeHash[m2].getGlobalGPACPose())
				#Th1 = matrix([[offset[0]],[offset[1]],[offset[2]]])
				#Ch1 = matrix([[0.2, 0.0, 0.0],
				#				[0.0, 0.2, 0.0],
				#				[0.0, 0.0, 0.04]])

				#poseOrigin = Pose(self.nodeHash[n1].getGlobalGPACPose())
				#offset = poseOrigin.convertGlobalPoseToLocal(self.nodeHash[n2].getGlobalGPACPose())
				#Th2 = matrix([[offset[0]],[offset[1]],[offset[2]]])
				#Ch2 = matrix([[0.2, 0.0, 0.0],
				#				[0.0, 0.2, 0.0],
				#				[0.0, 0.0, 0.04]])
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				#print
				#print "hypothesis:", i, j
				#print "pairwise hypotheses:", m1, m2, n1, n2
				#print "Th1,Ch1:", Th1, Ch1
				#print "Tp1,Ep1:", Tp1, Ep1
				#print "Th2, Ch2:", Th2, Ch2
				#print "Tp2, Ep2:", Tp2, Ep2
				
				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				#print "results:"
				#print "result1, covE:", result1, covE
				#print "result2, cov2:", result2, cov2
				#print "result3, cov3:", result3, cov3
				#print "result4, cov4:", result4, cov4
				
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])
				#print "invMat:", invMat
				#if [0,1,2].count(i) > 0 and [0,1,2].count(j) > 0:
				#	results.append([0.0, i, j])
				#	print "err:", 0.0
				#	#elif [3,4,5].count(i) > 0 and [3,4,5].count(j) > 0:
				#	#	results.append([0.0, i, j])
				#	#	print "err:", 0.0
				#else:
				#	results.append([1.0, i, j])
				#	print "err:", 1.0

					
		results.sort()
		selected = []	
		selected2 = []
		selected3 = []
		selected4 = []
		selected5 = []
		selected6 = []
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) == 0:
			print "no hypotheses!  returning "
			newConstraints = []
			#return
		
			#elif len(totalHypotheses) == 1:
		elif len(totalHypotheses) <= 2:
			newConstraints = totalHypotheses

		else: 
					
			" set all proscribed hypothesis sets to maximum error "	
			maxError = 0.0
			for result in results:
				if result[0] != None and result[0] > maxError:
					maxError = result[0]
	
			#maxError = maxError*2
	
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
				#A[i,j] = (maxError - result[0])/(2*maxError)
				#A[j,i] = (maxError - result[0])/(2*maxError)
				A[i,j] = 10*(maxError - result[0])
				A[j,i] = 10*(maxError - result[0])
	
			#print A
			#print "A:"
			#for i in range(len(totalHypotheses)):
			#	printStr = ""
			#	for j in range(len(totalHypotheses)):
			#		printStr += "%1.2f" % (A[i,j]) + " "
			#	
			#	print printStr
	
			" do graph clustering of consistency matrix "
			w = []
			e = []
			for i in range(100):

				#print "A:", A

				e, lmbda = scsgp.dominantEigenvectors2(A)
				#print "e:", e
				"""
				if len(totalHypotheses) > 3:
					kVal = 3
					lmbda2, e2 = eigsh(A, k=kVal)
				elif len(totalHypotheses) > 2:
					kVal = 2
					lmbda2, e2 = eigsh(A, k=kVal)
				else:
					kVal = 1
					lmbda2, e2 = eigsh(A, k=kVal)
				kVal = min(3, len(totalHypotheses)-1)
					
				eNew = [matrix([[e2[i,j]] for i in range(len(totalHypotheses))], dtype=float) for j in range(kVal)]
				"""
				#print "eNew:", eNew

				#print "lmbda:", lmbda, lmbda2
				#print "eigen:", e, eNew
				
				
				#exit()
				#print "lambda:", lmbda
				eigVec0 = [e[0][k,0] for k in range(len(totalHypotheses))]
				eigVec1 = [e[1][k,0] for k in range(len(totalHypotheses))]
				#print "eigVec0:", eigVec0
				#print "eigVec1:", eigVec1


				w1 = scsgp.getIndicatorVector(e[0])
				w2 = scsgp.getIndicatorVector(e[1])
				#w3 = scsgp.getIndicatorVector(e[2])
				if len(e) <= 1:
					break
	
				" threshold test l1 / l2"
				ratio = fabs(lmbda[0][0,0]) / fabs(lmbda[1][0,0])
				if ratio >= 2.0:
					break

			#if len(totalHypotheses) >= 2:
			#	kVal = max(3, len(totalHypotheses)-2)			
			#	lmbda2, e2 = eigs(A, k=kVal)
			#	print "eigenvalues:"
			#	print lmbda2
			#	print e2[0]
			#	print e2[1]
			

			if False and 3 < len(totalHypotheses):

				#print "A:", A
				
				kVal = min(3, len(totalHypotheses)-1)
				#print "kVal:", kVal
				
				lmbda2, e2 = eigsh(A, k=kVal)
				#print "e2:", e2

				eNew = [matrix([[e2[i,j]] for i in range(len(totalHypotheses))], dtype=float) for j in range(kVal)]
				#print "eNew:", eNew
				
				w1_2 = scsgp.getIndicatorVector(eNew[0])
				w2_2 = scsgp.getIndicatorVector(eNew[1])
				w3_2 = scsgp.getIndicatorVector(eNew[2])				

				#print "w1_2:", w1_2

				selected4 = []	
				for i in range(len(totalHypotheses)):
					if w1_2[i,0] >= 1.0:
						selected4.append(totalHypotheses[i])
		
				selected5 = []	
				for i in range(len(totalHypotheses)):
					if w2_2[i,0] >= 1.0:
						selected5.append(totalHypotheses[i])
	
				selected6 = []	
				for i in range(len(totalHypotheses)):
					if w3_2[i,0] >= 1.0:
						selected6.append(totalHypotheses[i])
		
			selected = []	
			for i in range(len(totalHypotheses)):
				if w1[i,0] >= 1.0:
					selected.append(totalHypotheses[i])
	
			selected2 = []	
			for i in range(len(totalHypotheses)):
				if w2[i,0] >= 1.0:
					selected2.append(totalHypotheses[i])

			#selected3 = []	
			#for i in range(len(totalHypotheses)):
			#	if w3[i,0] >= 1.0:
			#		selected3.append(totalHypotheses[i])


			newConstraints = selected

			#newConstraints = totalHypotheses


		print "total:", len(totalHypotheses), " = ", len(selected), "vs.", len(selected2), "hypotheses"
			
		print "adding", len(selected2), "medial overlap constraints on path"

		"""
		self.deleteAllPriority(CORNER_PRIORITY)
		self.mergePriorityConstraints()
		self.drawHyp(tit = "Empty")

		
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		self.drawHyp(tit = "Power eig 0")

		"""
		newConstraints = selected2
		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
				
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		#self.drawHyp(tit = "Power eig 1")

		if False and len(totalHypotheses) > 3:
			newConstraints = selected4
			" delete old edges "
			self.deleteAllPriority(CORNER_PRIORITY)
					
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()
			self.drawHyp(tit = "PCA 0")
	
	
			newConstraints = selected5
			" delete old edges "
			self.deleteAllPriority(CORNER_PRIORITY)
					
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()
			self.drawHyp(tit = "PCA 1" )
	
			newConstraints = selected6
			" delete old edges "
			self.deleteAllPriority(CORNER_PRIORITY)
					
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()
			self.drawHyp(tit = "PCA 2")

		"""
		newConstraints = selected2
		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
				
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		self.drawHyp(tit = "Power eig 1")
		"""
		
		#print "xi Error"
		#self.getXiError()
		
		#return selected, selected2
		#exit()

		return
		

		

	def addPathConstraints(self, pathNodes, targetNodeID, insertNode = False):

		print "addPathConstraints:"
		print "pathNodes:", pathNodes
		print "targetNodeID:", targetNodeID

		#pathNodes = []
		#for k in range(self.numPaths):
		#	pathNodes += self.nodeSets[k]

		print "pathNodes:", pathNodes
			
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
					oldConstraints = self.getEdges(targetNodeID, nodeID) + self.getEdges(nodeID, targetNodeID)
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
		#if len(constraintPairs) > 3:
		#	constraintPairs = constraintPairs[:3]
		#	print "reducing to only", len(constraintPairs)

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

		self.allPathConstraints += constraintResults

		if len(constraintResults) > 0 and insertNode:
			const = constraintResults[0]
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
			
			return


		" analyze over all possible constraints "
		totalHypotheses = self.allPathConstraints

		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
		
		" recompute dijkstra projection "
		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
		
		
		
		#totalHypotheses = constraintResults
		
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
		selected = []	
		selected2 = []
		selected3 = []
		selected4 = []
		selected5 = []
		selected6 = []
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) == 0:
			print "no hypotheses!  returning "
			newConstraints = []
			#return
		
			#elif len(totalHypotheses) == 1:
		elif len(totalHypotheses) <= 2:
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

				#print "A:", A

				e, lmbda = scsgp.dominantEigenvectors(A)
				#print "e:", e
				"""
				if len(totalHypotheses) > 3:
					kVal = 3
					lmbda2, e2 = eigsh(A, k=kVal)
				elif len(totalHypotheses) > 2:
					kVal = 2
					lmbda2, e2 = eigsh(A, k=kVal)
				else:
					kVal = 1
					lmbda2, e2 = eigsh(A, k=kVal)
				kVal = min(3, len(totalHypotheses)-1)
					
				eNew = [matrix([[e2[i,j]] for i in range(len(totalHypotheses))], dtype=float) for j in range(kVal)]
				"""
				#print "eNew:", eNew

				#print "lmbda:", lmbda, lmbda2
				#print "eigen:", e, eNew
				
				
				#exit()
				w1 = scsgp.getIndicatorVector(e[0])
				w2 = scsgp.getIndicatorVector(e[1])
				#w3 = scsgp.getIndicatorVector(e[2])
				if len(e) <= 1:
					break
	
				" threshold test l1 / l2"
				ratio = lmbda[0][0,0] / lmbda[1][0,0]
				if ratio >= 2.0:
					break

			if False and 3 < len(totalHypotheses):

				print "A:", A
				
				kVal = min(3, len(totalHypotheses)-1)
				print "kVal:", kVal
				
				lmbda2, e2 = eigsh(A, k=kVal)
				print "e2:", e2

				eNew = [matrix([[e2[i,j]] for i in range(len(totalHypotheses))], dtype=float) for j in range(kVal)]
				print "eNew:", eNew
				
				w1_2 = scsgp.getIndicatorVector(eNew[0])
				w2_2 = scsgp.getIndicatorVector(eNew[1])
				w3_2 = scsgp.getIndicatorVector(eNew[2])				

				print "w1_2:", w1_2

				selected4 = []	
				for i in range(len(totalHypotheses)):
					if w1_2[i,0] >= 1.0:
						selected4.append(totalHypotheses[i])
		
				selected5 = []	
				for i in range(len(totalHypotheses)):
					if w2_2[i,0] >= 1.0:
						selected5.append(totalHypotheses[i])
	
				selected6 = []	
				for i in range(len(totalHypotheses)):
					if w3_2[i,0] >= 1.0:
						selected6.append(totalHypotheses[i])
		
			selected = []	
			for i in range(len(totalHypotheses)):
				if w1[i,0] >= 1.0:
					selected.append(totalHypotheses[i])
	
			selected2 = []	
			for i in range(len(totalHypotheses)):
				if w2[i,0] >= 1.0:
					selected2.append(totalHypotheses[i])

			#selected3 = []	
			#for i in range(len(totalHypotheses)):
			#	if w3[i,0] >= 1.0:
			#		selected3.append(totalHypotheses[i])


			newConstraints = selected

			#newConstraints = totalHypotheses


		print "total:", len(totalHypotheses), " = ", len(selected), "vs.", len(selected2), "vs.", len(selected3), "hypotheses"
			
		print "adding", len(selected2), "medial overlap constraints on path"

		#self.deleteAllPriority(CORNER_PRIORITY)
		#self.mergePriorityConstraints()
		#self.drawHyp(tit = "Empty")

		"""
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		self.drawHyp(tit = "Power eig 0")
		newConstraints = selected3
		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
				
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		self.drawHyp(tit = "Power eig 2")
		"""

		if False and len(totalHypotheses) > 3:
			newConstraints = selected4
			" delete old edges "
			self.deleteAllPriority(CORNER_PRIORITY)
					
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()
			self.drawHyp(tit = "PCA 0")
	
	
			newConstraints = selected5
			" delete old edges "
			self.deleteAllPriority(CORNER_PRIORITY)
					
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()
			self.drawHyp(tit = "PCA 1" )
	
			newConstraints = selected6
			" delete old edges "
			self.deleteAllPriority(CORNER_PRIORITY)
					
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()
			self.drawHyp(tit = "PCA 2")


		newConstraints = selected2
		" delete old edges "
		self.deleteAllPriority(CORNER_PRIORITY)
				
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		#self.drawHyp(tit = "Power eig 1")

		
		
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
		
		v_list.sort()
		e_list.sort()
		
		
		
		for mn in range(10):
			print "v_list:", repr(v_list)
			print "e_list:", repr(e_list)

		#toro.writeConstrainedToroGraph(".", fileName + ".graph", v_list, e_list)
		#toro.executeToro("./" + fileName + ".graph")
	
		#finalFileName = fileName + "-treeopt-final.graph"
		
		v_list2, e_list2 = toro.runTORO(v_list,e_list)
		
		#v_list2, e_list2 = toro.readToroGraph("./" + finalFileName)		

		v_list2.sort()
		e_list2.sort()
		
		print "v_list2:", repr(v_list2)
		print "e_list2:", repr(e_list2)
			
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
			

	def drawHyp(self, tit = ""):
		
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

		err, edgeTotal = self.getXiError()
		errAvg = err/edgeTotal

		
		pylab.title(tit + " %d Poses and Xi2 = %f %f" % (self.numNodes,err, errAvg))
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)

			
		pylab.savefig("plotHypotheses%04u.png" % self.hypPlotCount)
		self.hypPlotCount += 1

	def drawCandidate(self, cand):


		const_count = 0
		nodeID1 = cand["nodeID1"]
		nodeID2 = cand["nodeID2"]
		transform = cand["transform"]
		covE = cand["covE"]
		offset = cand["offset"]
		medial1 = cand["medial1"]
		medial2 = cand["medial2"]

		angDiff1 = cand["angDiff1"]
		cartDist1 = cand["cartDist1"]
		angDiff2 = cand["angDiff2"]
		cartDist2 = cand["cartDist2"]
		histogram = cand["histogram"]

		u1 = cand["u1"]
		u2 = cand["u2"]
		initAng = cand["initAng"]
		
		estPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()		
		
		hull1 = computeBareHull(self.nodeHash[nodeID1], sweep = False)
		hull1.append(hull1[0])

		hull2 = computeBareHull(self.nodeHash[nodeID2], sweep = False)
		hull2.append(hull2[0])

		#offset = [transform[0,0], transform[1,0], transform[2,0]]

		#hull2_trans = []
		#for p in hull2:
		#	hull2_trans.append(gen_icp.dispOffset(p, offset))	

		" set the origin of pose 1 "
		poseOrigin = Pose(estPose1)

		pylab.clf()
		pylab.axes()
		
		xP = []
		yP = []
		for b in hull1:
			p1 = poseOrigin.convertLocalToGlobal(b)

			xP.append(p1[0])	
			yP.append(p1[1])

		pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))


		xP = []
		yP = []
		for p in hull2:
			p1 = gen_icp.dispOffset(p,offset)
			
			p2 = poseOrigin.convertLocalToGlobal(p1)
			xP.append(p2[0])	
			yP.append(p2[1])

		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))


		" create the ground constraints "
		gndGPAC1Pose = self.nodeHash[nodeID1].getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)
		gndGPAC2Pose = self.nodeHash[nodeID2].getGndGlobalGPACPose()
		gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)

		xP = []
		yP = []
		for p in hull2:
			p1 = gen_icp.dispOffset(p,gndOffset)
			
			p2 = poseOrigin.convertLocalToGlobal(p1)
			xP.append(p2[0])	
			yP.append(p2[1])

		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))


		medialPath1 = medial1
		medialPath2 = medial2
		medial_trans = []
		for p in medialPath2:
			medial_trans.append(gen_icp.dispOffset(p, offset))	

		xP = []
		yP = []
		for p in medialPath1:
			p1 = poseOrigin.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

		xP = []
		yP = []
		for p in medial_trans:
			p1 = poseOrigin.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

		xP = []
		yP = []
		for p in medialPath2:
			p1 = gen_icp.dispOffset(p, gndOffset)
			p2 = poseOrigin.convertLocalToGlobal(p1)
			xP.append(p2[0])
			yP.append(p2[1])
		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

		self.drawWalls()

		pylab.title("%u -> %u, [%1.1f %1.1f %1.1f] %1.2f %1.2f %1.2f %1.2f %d %d %d" % (nodeID1, nodeID2, u1,u2,initAng, angDiff1, cartDist1, cartDist2, angDiff2, histogram[0], histogram[1], histogram[2]))


		pylab.xlim(estPose1[0]-4, estPose1[0]+4)					
		pylab.ylim(estPose1[1]-3, estPose1[1]+3)

		pylab.savefig("candidate_%04u.png" % self.candidateCount)
			
		#pylab.savefig("constraints_%04u.png" % const_count)
		pylab.clf()			

		self.candidateCount += 1

	
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
		
		#for k in range(len(trimmedPaths)):
		for k,path in trimmedPaths.iteritems():
			#path = trimmedPaths[k]
			print "path has", len(path), "points"
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color = self.colors[k])

		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		self.drawWalls()

		pylab.title("Trimmed Paths, numNodes = %d" % self.numNodes)
		pylab.savefig("trimmedPath_%04u.png" % self.trimCount)
		self.trimCount += 1

	def drawWalls(self):
		for wall in self.walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = wall[i]
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g')



	def drawPathAndHull(self):
		
		
		pylab.clf()
		pathIDs = self.paths.getPathIDs()

		allNodes = []
		for k in pathIDs:
			nodeSet = self.paths.getNodes(k)
			allNodes += copy(nodeSet)
		
		allNodes.sort()	
		
		if len(allNodes) > 0:
			highestNodeID = allNodes[-1]
		else:
			highestNodeID = 1e100		
			
			
		for k in pathIDs:
			xP = []
			yP = []
			
			for p in self.paths.paths[k]:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=self.colors[k], linewidth=4)


			nodeSet = self.paths.getNodes(k)

			print "drawing pathID", k, "for nodes:", nodeSet
			for nodeID in nodeSet:
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
				
				if nodeID == highestNodeID:
					pylab.plot(xP,yP, color=(0,0,0))
				elif nodeID == highestNodeID-1:
					pylab.plot(xP,yP, color=(0.5,0.5,0.5))
				else:
					pylab.plot(xP,yP, color=self.colors[k])
					
			
			xP = []
			yP = []
			for p in self.paths.hulls[k]:
				xP.append(p[0])
				yP.append(p[1])
				
			pylab.plot(xP,yP, '--', color=self.colors[k], linewidth=4)
			
		print "pathAndHull:", self.pathDrawCount
		printStack()

		pylab.title("paths: %s numNodes: %d %d" % (repr(self.paths.getPathIDs()), self.numNodes, highestNodeID))
		pylab.savefig("pathAndHull_%04u.png" % self.pathDrawCount)

		self.pathDrawCount += 1
			
	def getXiError(self):
		
		edgeHash = self.getPriorityEdges()
		errSum = 0.0
		edgeTotal = 0
		
		for k, v in edgeHash.iteritems():
			
			nodeID1 = k[0]
			nodeID2 = k[1]
			
			
			transform1 = v[0][0]
			covE = v[0][1]
			
			offset1 = [transform1[0,0],transform1[1,0],transform1[2,0]]
			
			#print
			#print offset1
			
			#print transform1
			#print covE
			

			estPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()
			estPose2 = self.nodeHash[nodeID2].getGlobalGPACPose()

			poseOrigin = Pose(estPose1)
			offset2 = poseOrigin.convertGlobalPoseToLocal(estPose2)
			
			transform2 = matrix([[offset2[0]],[offset2[1]],[offset2[2]]])
			#print offset2


			invMat = scipy.linalg.inv(covE)
			#print "invMat:", invMat
			
			diff = transform2 - transform1
			#print "diff:", diff
			
			err = numpy.transpose(diff) * invMat * diff

			#print nodeID1, nodeID2, "Xi^2 err:", err[0,0]

			errSum += err[0,0]
			edgeTotal += 1
			
		return errSum, edgeTotal



