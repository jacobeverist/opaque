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
from Paths import Paths, computePathAngleVariance
import gen_icp
from functions import *
import scsgp
import bayes
from operator import itemgetter
import cProfile
import time
import traceback

from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath
import pylab
import graph

renderCount = 0
topCount = 0

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
		
		
		self.numNodes = 0
		self.currNode = 0

		self.paths = Paths(self.nodeHash)
		
		self.pathPlotCount = 0
		self.overlapPlotCount = 0
		self.pathDrawCount = 0
		self.candidateCount = 0
		self.statePlotCount = 0
		self.commonOriginCount = 0
		self.multiDepCount = 0
		

		#self.numPaths = 1
		#self.paths = {0 : []}
		#self.hulls = {0 : []}
		#self.nodeSets = {0 : []}
		
		
		" FIXME:  self.currPath is used but never modified	"
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
							[ 0.0,	0.01, 0.0],
							[0.0, 0.0,	0.02]])

		self.E_featureless = matrix([[ 1.0,  0.0, 0.0],
							[ 0.0,	0.01, 0.0],
							[0.0, 0.0,	0.02]])
		self.E_motion = matrix([[ 0.06, 0.0, 0.0],
							[0.0,  5.0, 0.0],
							[0.0, 0.0,	0.1 ]])
		
		self.E_sensor = matrix([[0.1,0.0,0.0],
							[0.0,0.05,0.0],
							[0.0,0.0,0.02]])
		
		self.E_relaxed = matrix([[ 1.0, 0.0, 0.0],
							[ 0.0, 0.1, 0.0],
							[0.0, 0.0, pi/8.0]])


		self.allPathConstraints = []
		self.hypPlotCount = 0


		self.nodePoses = {}

	@logFunction
	def saveState(self):

		#self.nodePoses = {}
		for k in range(self.numNodes):
			self.nodePoses[k] = self.nodeHash[k].getGlobalGPACPose()
					
		saveFile = ""

		saveFile += "self.walls = " + repr(self.walls) + "\n"
		saveFile += "self.initPose = " + repr(self.initPose) + "\n"
		#saveFile += "self.nodeHash = " + repr(self.nodeHash) + "\n"
		saveFile += "self.numNodes = " + repr(self.numNodes) + "\n"
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
	
		
	@logFunction
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/stateSave_%04u.txt" % (numNodes-1)
		f = open(dirName + "/stateSave_%04u.txt" % (numNodes-1), 'r')		
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		exec(saveStr)
		
		
		tempDict = {}
		for k in range(len(self.nodePoses)):
			tempDict[k] = self.nodePoses[k]
			
		self.nodePoses = tempDict
		
		print self.numNodes
		
		self.paths.restoreState(dirName, numNodes)
		



	@logFunction
	def pairLastTwo(self):
	
		node1 = self.nodeHash[self.numNodes-1]
		node2 = self.nodeHash[self.numNodes-2]

		node1.setPartnerNodeID(self.numNodes-2)
		node2.setPartnerNodeID(self.numNodes-1)



	@logFunction
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


	@logFunction
	def addNode(self, newNode):
		
		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1


	@logFunction
	def loadNewNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		
		
		return self.integrateNode(newNode, nodeID)



	@logFunction
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

				" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
				if nodeID1 % 2 == 1:		
					" Move node along path "
					self.movePath(nodeID1, direction, distEst = 0.0)
				elif nodeID2 % 2 == 1:		
					" Move node along path "
					self.movePath(nodeID2, direction, distEst = 0.0)



			" if these nodes are already path-classified, return"
			isContained1 = False
			isContained2 = False
			
			pathIDs = self.paths.getPathIDs()
			for k in pathIDs:
				if self.paths.getNodes(k).count(nodeID1) > 0:
					isContained1 = True
				if self.paths.getNodes(k).count(nodeID2) > 0:
					isContained2 = True
					
					
			if not isContained1 and not isContained2:
				
				" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
				self.paths.generatePaths()
	
				self.drawPathAndHull()
				
	
				self.addToPaths(nodeID1, nodeID2)
	
				" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
				if len(self.paths.paths[0]) > 0:
					
					#self.paths.comparePaths()
					self.mergePaths()
	
					paths = {}
					pathIDs = self.paths.getPathIDs()
	
	
	
					orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
					print nodeID1, "recompute orderedPathIDs1:", orderedPathIDs1
	
				
					orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
					print nodeID2, "recompute orderedPathIDs2:", orderedPathIDs2
	
	
					for k in pathIDs:
						path = self.paths.paths[k]
						if len(path) > 0:
							paths[k] = path 
					self.trimmedPaths = self.paths.trimPaths(paths)
	
					
					self.drawTrimmedPaths(self.trimmedPaths)
					print "trimmed paths:", len(self.trimmedPaths)
	
	
	
					"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
					"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
					"3)  select choice with the lowest cost "
					"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
					
					" nodeID1:	is it in a junction or a single path? "
	
					hull1, medial1 = computeHullAxis(nodeID1, self.nodeHash[nodeID1], tailCutOff = False)
					hull2, medial2 = computeHullAxis(nodeID2, self.nodeHash[nodeID2], tailCutOff = False)
					estPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()
					estPose2 = self.nodeHash[nodeID2].getGlobalGPACPose()
					self.currSplicePath = self.selectSplice(nodeID1, nodeID2, medial1, medial2, estPose1, estPose2, orderedPathIDs1, orderedPathIDs2)
					
					self.paths.generatePaths()
					self.drawPathAndHull()
					self.drawTrimmedPaths(self.trimmedPaths)

			
			self.paths.generatePaths()
			self.drawPathAndHull()


		" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
		if self.numNodes >= 4:
			pass
		
		
				

		for k in range(self.numNodes):
			self.nodePoses[k] = self.nodeHash[k].getGlobalGPACPose()


	@logFunction
	def integrateNode(self, newNode, nodeID):


		" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
		direction = newNode.travelDir

		" ensure the axes are computed before this check "
		computeHullAxis(nodeID, newNode, tailCutOff = False)
		
		
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

				
			self.drawConstraints(self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull()


			self.addToPaths(nodeID1, nodeID2)

			" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
			if len(self.paths.paths[0]) > 0:
				
				self.mergePaths()

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				paths = {}
				pathIDs = self.paths.getPathIDs()



				orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
				print nodeID1, "recompute orderedPathIDs1:", orderedPathIDs1

			
				orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
				print nodeID2, "recompute orderedPathIDs2:", orderedPathIDs2


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
				
				" nodeID1:	is it in a junction or a single path? "
				
				
				hull1, medial1 = computeHullAxis(nodeID1, self.nodeHash[nodeID1], tailCutOff = False)
				hull2, medial2 = computeHullAxis(nodeID2, self.nodeHash[nodeID2], tailCutOff = False)
				estPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()
				estPose2 = self.nodeHash[nodeID2].getGlobalGPACPose()
				self.currSplicePath = self.selectSplice(nodeID1, nodeID2, medial1, medial2, estPose1, estPose2, orderedPathIDs1, orderedPathIDs2)
				
				self.paths.generatePaths()
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				self.drawTrimmedPaths(self.trimmedPaths)

				try:
					self.consistentFit(nodeID1,  self.nodeHash[nodeID1].getGlobalGPACPose(), numGuesses = 11)
				except IndexError:
					print "failed to consistentFit node", nodeID1
					pass
					
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				try:
					self.consistentFit(nodeID2,  self.nodeHash[nodeID2].getGlobalGPACPose(), numGuesses = 11)
				except IndexError:
					print "failed to consistentFit node", nodeID2
					pass

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				
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
		
				
		if nodeID > 0:
			
			" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
			if self.numNodes >= 4:

				" Move node along path "
				self.movePath(nodeID, direction)
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				

		for k in range(self.numNodes):
			self.nodePoses[k] = self.nodeHash[k].getGlobalGPACPose()


	@logFunction
	def getInPlaceGuess(self, nodeID1, nodeID2, direction):
		
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

		poseOrigin = Pose(self.nodeHash[nodeID1].getGlobalGPACPose())

		if len(supportLine) > 0 and resultSum1 < resultSum3 and resultSum3 - resultSum1 > 100.0:
			estPose2 = poseOrigin.convertLocalOffsetToGlobal(offset1)
			self.nodeHash[nodeID2].setGPACPose(estPose2)
		else:
			estPose2 = poseOrigin.convertLocalOffsetToGlobal(offset3)
			self.nodeHash[nodeID2].setGPACPose(estPose2)
	

	@logFunction
	def getStepGuess(self, nodeID, direction):
		
		" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
		if nodeID >= 2:
			
			if self.nodeHash[nodeID-2].getNumLeafs() > 2 or self.nodeHash[nodeID].getNumLeafs() > 2:
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

			offset = [transform[0,0], transform[1,0], transform[2,0]]			
			poseOrigin = Pose(self.nodeHash[nodeID-2].getGlobalGPACPose())
			estPose2 = poseOrigin.convertLocalOffsetToGlobal(offset)
			self.nodeHash[nodeID].setGPACPose(estPose2)


	@logFunction
	def addToPaths(self, nodeID1, nodeID2):

		" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
		if len(self.paths.paths[0]) == 0:
			" first nodes in path 0 "								
			self.paths.addNode(nodeID1,0)
			self.paths.addNode(nodeID2,0)
	
			self.paths.generatePaths()
			self.currSplicePath = self.paths.paths[0]
	
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
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1, plotIter = True)
				departures1.append([isExist1,isExist2])
				interiors1.append([isInterior1, isInterior2])
				depPoints1.append([departurePoint1, departurePoint2])
				distances1.append([discDist1, discDist2])
				depAngles1.append([depAngle1, depAngle2])
				contig1.append((contigFrac, overlapSum))
		
			for pathID in orderedPathIDs2:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID2, plotIter = True)
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
				print "reporting discrepancy"
				print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2
	
			
			print "frontAngDiff:", frontAngDiff
			print "backAngDiff:", backAngDiff
			
			
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
	
			isBranch, pathBranchIDs, isNew = self.paths.determineBranchPair(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)
	
			print "determineBranchPair:"
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
	
			isBranch, pathBranchIDs, isNew = self.paths.determineBranchPair(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag)
			print "determineBranchPair:"
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
						
		


	@logFunction
	def movePath(self, nodeID, direction, distEst = 1.0):
		
		print "movePath(", nodeID, ",", direction, ",", distEst, ")"
		
		if nodeID > 0:
			
			if nodeID % 2 == 1:
				self.getStepGuess(nodeID-1, direction)
				self.getInPlaceGuess(nodeID-1, nodeID, direction)
			else:
				print "movePath:  DO NOTHING"
			
			" guess path state from these guesses "
			if False: 
			
				" 1) splice is signified by the paired terminals (t0,t1) "
				" 2) splice can only change by one terminal.  (t0,t1) -> (t0,t2) "
				" assume at most 3 paths for possible splice states "
	
				" junctions indexed by their pathID "
				" terminals indexed by pathID+1,  root is 0 and 1 "
				" t%u % (pathID+1) dictionary returns pathID and index "
				" j%u % (pathID) dictionary returns pathID and index "
				
				" 3) splice state for node0 and node1 are same or separate by one terminal delta "
				" 4) splice state for node2 and node3 are at most one delta away from node0 and node1 "
				" 5) generate possible splice states for node2 given splice state of node0 "
				" 6) apply displacement of node2 from node0's position for every given splice "
				" 7) compute its fit according to getMultiDeparture(), classifying by matchCount, cost, and contiguity (or ICP?) (see selectSplice function)"
				" 8) select its most likely splice state "
				" 9) repeat for node1 and node3 "
				" 10) difference of splice state between node2 and node3 is at most 1, but usually zero "
				
				
				
				"""		
				self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]
				self.terminals[pathID+1] = [(pathID,0), path[0]]
				self.topDict["t%u" % (pathID+1)] = (pathID,0)
				termIDs = [self.topDict["j%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
				"""
				
			" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
			if self.numNodes >= 4 and nodeID % 2 == 1:
				
				" 1) get spliced path that the previous node is on "
				#self.currSplicePath = self.paths.paths[0]
				hull2, medial2 = computeHullAxis(nodeID-1, self.nodeHash[nodeID-1], tailCutOff = False)
				hull3, medial3 = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)
				
				allSplices, terminals, junctions = self.paths.getAllSplices(plotIter = True)

				initPose2 = self.nodeHash[nodeID-1].getGlobalGPACPose()
				initPose3 = self.nodeHash[nodeID].getGlobalGPACPose()
				
				print "junctions:", junctions
				print "initPose2:", initPose2
				print "initPose3:", initPose3
				
				junctions2 = []
				junctions3 = []
				for pathID, params in junctions.iteritems():
					junctionPose = params[1]
					
					print "checking", pathID, params
					
					dist2 = sqrt((junctionPose[0]-initPose2[0])**2 + (junctionPose[1]-initPose2[1])**2)
					dist3 = sqrt((junctionPose[0]-initPose3[0])**2 + (junctionPose[1]-initPose3[1])**2)
					
					print "dist2 =", dist2
					print "dist3 =", dist3
					
					if dist2 < 3.0:
						junctions2.append((pathID,params[2][0]))

					if dist3 < 3.0:
						junctions3.append((pathID,params[2][0]))
							
					
				"self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]"
				
				closePathID2 = self.paths.getClosestPath(initPose2)

				print "closePathID2:", closePathID2
				
				pathSet2 = [closePathID2]

				junctionSet = junctions2 + junctions3
				junctionSet = list(set(junctionSet))

				print "junctions2:", junctions2
				print "junctions3:", junctions3

				for junc in junctionSet:
					pathSet2.append(junc[0])
					pathSet2.append(junc[1])

				print "pathSet2:", pathSet2

				pathSet2 = list(set(pathSet2))
				
				pathSet2.sort()

				print "pathSet2:", pathSet2
				
				termSet2 = []
				for pathID in pathSet2:
					if pathID == 0:
						termSet2.append(terminals[0][0])
						termSet2.append(terminals[1][0])
					else:
						termSet2.append(terminals[pathID+1][0])
						
				for pathID in pathSet2:
					parentID = self.paths.getParentPathID(pathID)

					if parentID != None and not parentID in pathSet2:
						termSet2.append(junctions[pathID][2])
					

				

				splicePaths = []

				print "termSet2:", termSet2
				print "allSplices:"

				for k, result in allSplices.iteritems():
					print "id:", k
					if k[0] in pathSet2 and k[1] in pathSet2:
						
						for sPath in result:
							termPath = sPath['termPath']
							print "termPath:", termPath
														
							if termPath[0] in termSet2 and termPath[-1] in termSet2:
								splicePaths.append(sPath['path'])
								
										
				resultMoves2 = []
				resultMoves3 = []

				for path in splicePaths:		
	
					" 2) get pose of previous node, get point on path curve "
					pose0 = self.nodeHash[nodeID-3].getGlobalGPACPose()
					pose1 = self.nodeHash[nodeID-2].getGlobalGPACPose()
	
					" FIXME:  make sure path is oriented correctly wrt node medial axis "
					hull0, medial0 = computeHullAxis(nodeID-3, self.nodeHash[nodeID-3], tailCutOff = False)
					
					poseOrigin0 = Pose(pose0)
					globalMedial0 = []
					for p in medial0:
						globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))
		
					orientedSplicePath = orientPath(path, globalMedial0)				
					#orientedSplicePath = self.orientPath(self.currSplicePath, medial0, pose0)				
					currPathSpline = SplineFit(orientedSplicePath, smooth=0.1)
	
					
					print "pose0,pose1:", pose0, pose1
					
					minDist0, u0, p0 = currPathSpline.findClosestPoint(pose0[0:2])
					minDist1, u1, p1 = currPathSpline.findClosestPoint(pose1[0:2])
					
					print "u0,u1:", u0, u1
					print "len(distPoints):", len(currPathSpline.distPoints)
					
					arcDist0 = currPathSpline.distPoints[int(u0*1000)]
					arcDist1 = currPathSpline.distPoints[int(u1*1000)]
	
					
					" 3) step distance along the path in direction of travel "
					#distEst = 0.6			
					#distEst = 1.0			
					
					if direction:
						arcDist2 = arcDist0 - distEst
						arcDist3 = arcDist1 - distEst
					else:
						arcDist2 = arcDist0 + distEst
						arcDist3 = arcDist1 + distEst
					
					print "arcDist0, arcDist1, arcDist2, arcDist3:", arcDist0, arcDist1, arcDist2, arcDist3
					
					"""
					pose2 = currPathSpline.getPointOfDist(arcDist2)
					pose3 = currPathSpline.getPointOfDist(arcDist3)
					"""
					
					
					pose2 = initPose2
					pose3 = initPose3
					
					
					print "pose2,pose3:", pose2, pose3
					
					" 4) set as pose of new node "
					self.nodeHash[nodeID-1].setGPACPose(pose2)
					self.nodeHash[nodeID].setGPACPose(pose3)
	
					uPath2, uMedialOrigin2 = self.selectLocalCommonOrigin(orientedSplicePath, medial2, pose2)
	
					uPath3, uMedialOrigin3 = self.selectLocalCommonOrigin(orientedSplicePath, medial3, pose3)
					
					u2 = uPath2
					u3 = uPath3
	
					print "input: uMedialOrigin2, u2, pose2:", uMedialOrigin2, u2, pose2
					print "input: uMedialOrigin3, u3, pose3:", uMedialOrigin3, u3, pose3
	
					
					resultPose2, lastCost2, matchCount2, currAng2 = gen_icp.globalPathToNodeOverlapICP2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial2, plotIter = True, n1 = nodeID-1, n2 = -1, arcLimit = 0.1)
					resultPose3, lastCost3, matchCount3, currAng3 = gen_icp.globalPathToNodeOverlapICP2([u3, uMedialOrigin3, 0.0], orientedSplicePath, medial3, plotIter = True, n1 = nodeID, n2 = -1, arcLimit = 0.1)
					
					print "resultPoses:", resultPose2, resultPose3

					result2 = getMultiDeparturePoint(orientedSplicePath, medial2, pose2, resultPose2, [], nodeID-1, pathPlotCount = self.multiDepCount, plotIter = True)
					self.multiDepCount += 1
					result3 = getMultiDeparturePoint(orientedSplicePath, medial3, pose3, resultPose3, [], nodeID, pathPlotCount = self.multiDepCount, plotIter = True)
					self.multiDepCount += 1

					#results1.append(result+(k,))
					" (departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 )"
					
					" (resultPose2,lastCost2,matchCount2,fabs(currAng2)) "
					
					angDiff2 = abs(diffAngle(pose2[2],resultPose2[2]))
					angDiff3 = abs(diffAngle(pose3[2],resultPose3[2]))
					
					print "angDiff:", angDiff2, angDiff3
					
					" NOTE:  added minimum threshold to angle difference "
					" NOTE:  guess pose is now the inter-nodal estimate instead of path-based estimate "
					if angDiff2 < 0.5:
						#resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(currAng2)) + result2)
						resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(0.0)) + result2)
					if angDiff3 < 0.5:
						#resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(currAng3)) + result3)				
						resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(0.0)) + result3)				
				
				
				
				resultMoves2 = sorted(resultMoves2, key=itemgetter(18))
				resultMoves2 = sorted(resultMoves2, key=itemgetter(16), reverse=True)

				resultMoves3 = sorted(resultMoves3, key=itemgetter(18))
				resultMoves3 = sorted(resultMoves3, key=itemgetter(16), reverse=True)
				
				##resultMoves = sorted(resultMoves, key=itemgetter(3))

				print "resultMoves2:"
				for res in resultMoves2:
					print res
				print "resultMoves3:"
				for res in resultMoves3:
					print res

				if len(resultMoves2) > 0:
					self.nodeHash[nodeID-1].setGPACPose(resultMoves2[0][0])
				else:
					print "node", nodeID-1, "not movePathed because no valid pose"
				if len(resultMoves3) > 0:
					self.nodeHash[nodeID].setGPACPose(resultMoves3[0][0])
				else:
					print "node", nodeID, "not movePathed because no valid pose"


	

	@logFunction
	def selectSplice(self, nodeID1, nodeID2, medial1, medial2, estPose1, estPose2, orderedPathIDs1, orderedPathIDs2):
		
		splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs1)
		splicedPaths2 = self.paths.splicePathIDs(orderedPathIDs2)

		print "received", len(splicedPaths1), "spliced paths from path IDs", orderedPathIDs1
		print "received", len(splicedPaths2), "spliced paths from path IDs", orderedPathIDs2

		" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 |"
		" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 "

		results1 = []		
		for k in range(len(splicedPaths1)):
			path = splicedPaths1[k]			
			result = getMultiDeparturePoint(path, medial1, estPose1, estPose1, orderedPathIDs1, nodeID1)
			results1.append(result+(k,))


		results1 = sorted(results1, key=itemgetter(14))
		results1 = sorted(results1, key=itemgetter(12), reverse=True)				

		results2 = []		
		for k in range(len(splicedPaths2)):
			path = splicedPaths2[k]			
			result = getMultiDeparturePoint(path, medial2, estPose2, estPose2, orderedPathIDs2, nodeID2)
			results2.append(result+(k,))

		results2 = sorted(results2, key=itemgetter(14))
		results2 = sorted(results2, key=itemgetter(12), reverse=True)				


		kIndex = results1[0][15]

		return splicedPaths1[kIndex]

	@logFunction
	def selectLocalCommonOrigin(self, globalPath, medial1, estPose1):
		poseOrigin = Pose(estPose1)
		
		" FIXME:  change function so that it always returns a common pair, biasing towards it's current location "
		" alternately, make a new function that gives the option to enforce locality to common pair "
		
		
		#globalMedial = []
		#for p in medial1:
		#	globalMedial.append(poseOrigin.convertLocalToGlobal(p))
		
		#medialSpline1 = SplineFit(globalMedial, smooth=0.1)

		globalSpline = SplineFit(globalPath, smooth=0.1)
		medialSpline1 = SplineFit(medial1, smooth=0.1)


		globalSamples = globalSpline.getUniformSamples(spacing = 0.04)
		medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
		
		globalMedialSamples = []
		for p in medialSamples:
			result = poseOrigin.convertLocalOffsetToGlobal(p)	
			globalMedialSamples.append(result)
		
		" compute the local variance of the angle "
		globalVar = computePathAngleVariance(globalSamples)
		medialVar = computePathAngleVariance(globalMedialSamples)

	  
		"""
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
		"""


		" now lets find closest points and save their local variances "			
		closestPairs = []
		allPairs = []
		TERM_DIST = 20

		pathRail = int(len(globalSamples) / 20.0)
		medialRail = int(len(globalMedialSamples) / 20.0)

		print "rails:", len(globalSamples), len(globalMedialSamples), pathRail, medialRail

		
		for i in range(pathRail, len(globalSamples)-pathRail):
			pG = globalSamples[i]
			minDist = 1e100
			minJ = -1
			for j in range(medialRail, len(globalMedialSamples)-medialRail):
				pM = globalMedialSamples[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				
				
				if dist < minDist:
					minDist = dist
					minJ = j
			
			#poseDist = math.sqrt((pG[0]-estPose1[0])**2 + (pG[1]-estPose1[1])**2)
			#if poseDist < 0.3:
			if True:
				allPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))
					
			#if minDist < 0.1:
			#	closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

		for j in range(medialRail, len(globalMedialSamples)-medialRail):
			pM = globalMedialSamples[j]
			minDist = 1e100
			minI = -1
			for i in range(pathRail, len(globalSamples)-pathRail):
				pG = globalSamples[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				
				
				if dist < minDist:
					minDist = dist
					minI = i

			#pG = globalSamples[minI]
			#poseDist = math.sqrt((pG[0]-estPose1[0])**2 + (pG[1]-estPose1[1])**2)
			#if poseDist < 0.3:
			if True:		
				allPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))
			
			#if minDist < 0.1:
			#	closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))



		" remove duplicates "
		allPairs = list(set(allPairs))



		
		allPairs = sorted(allPairs, key=itemgetter(2))
		
		maxDistThresh = allPairs[-1][2]
		minDistThresh = 0.1
		
		print "minDistThresh,maxDistThresh =", allPairs[0][2], allPairs[-1][2]


		if len(allPairs) == 0:
			raise

		
		originU2 = 0.5
		originU1 = 0.5
		
		while minDistThresh <= maxDistThresh:
			closestPairs = []
			for pair in allPairs:			
				if pair[2] < minDistThresh:				
					closestPairs.append(pair)
			
			" sort by lowest angular variance"
			closestPairs = sorted(closestPairs, key=itemgetter(5,6))

			print len(closestPairs), "closest pairs for dist", minDistThresh

			if len(closestPairs) > 0:
				originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
				originU1 = globalSpline.findU(globalSamples[closestPairs[0][0]])
				
				break
			
			minDistThresh += 0.1
		
		u2 = originU2
		u1 = originU1
		angGuess = 0.0
		
		return u1, u2
			

	@logFunction
	def mergeSiblings(self, resultOffset1, resultOffset2, pathID1, pathID2, lastCost, matchCount):

		#tempID = pathID1
		#pathID1 = pathID2
		#pathID2 = tempID
		
		#tempOffset = resultOffset1
		#resultOffset1 = resultOffset2
		#resultOffset2 = tempOffset

		" verify that this path still exists before we try to merge it "
		allPathIDs = self.paths.getPathIDs()
		if pathID1 in allPathIDs and pathID2 in allPathIDs:
	
			parentID1 = self.paths.pathClasses[pathID1]["parentID"] 
			parentID2 = self.paths.pathClasses[pathID2]["parentID"] 
	
			if parentID1 != parentID2:
				print "ERROR mergeSiblings:", pathID1, pathID2, parentID1, parentID2, "not siblings"
				return
	
			
			" capture transform between pathID1 and pathID2 under correction "
	
			#cOffset = resultOffset1
			#rootID = pathID1
			#lessID = pathID2
			
			
			mergeNodes1 = copy(self.paths.getNodes(pathID1))
			mergeNodes2 = copy(self.paths.getNodes(pathID2))
			
			print "mergeNodes1:", mergeNodes1
			print "mergeNodes2:", mergeNodes2
			
			duplicateNodes = []
			" look to see if a node is in both paths "
			for nodeID1 in mergeNodes1:
				for nodeID2 in mergeNodes2:
					
					if nodeID1 == nodeID2:
						duplicateNodes.append(nodeID1)
						
			" if any duplicates, lets choose one path or the other "
			for nodeID in duplicateNodes:
				orderedPathIDs = self.paths.getOrderedOverlappingPaths(nodeID)
				print "duplicate node", nodeID, "is in paths:", orderedPathIDs
				if pathID1 in orderedPathIDs:
					mergeNodes2.remove(nodeID)
				elif pathID2 in orderedPathIDs:
					mergeNodes1.remove(nodeID)
				else:
					print "duplicate not in either sibling path"
					raise

			print "mergeNodes1:", mergeNodes1
			print "mergeNodes2:", mergeNodes2
				
			
			poseOrigin1 = Pose(resultOffset1)
			poseOrigin2 = Pose(resultOffset2)

			" move a node only once, keep track of whether we've moved a node or not "
			hasMoved = {}
			for nodeID in mergeNodes1:
				hasMoved[nodeID] = False
			for nodeID in mergeNodes2:
				hasMoved[nodeID] = False

			
			" capture transform between pathID2 and member nodes "
			" add compute new node locations under corrected transform "
			self.drawConstraints(self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull()
						
			" offset between paths "
			for nodeID in mergeNodes1:
	
				if not hasMoved[nodeID]:
					print "modifying node", nodeID, "in path", pathID1 
					" adjust node poses to relative path positions "
					nodePose = self.nodeHash[nodeID].getGlobalGPACPose()
					guessPose = poseOrigin1.convertLocalOffsetToGlobal(nodePose)
					self.nodeHash[nodeID].setGPACPose(guessPose)
					
					self.drawConstraints(self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull()
					
					hasMoved[nodeID] = True

			" offset between paths "
			for nodeID in mergeNodes2:
	
				if not hasMoved[nodeID]:
					print "modifying node", nodeID, "in path", pathID2 
					" adjust node poses to relative path positions "
					nodePose = self.nodeHash[nodeID].getGlobalGPACPose()
					guessPose = poseOrigin2.convertLocalOffsetToGlobal(nodePose)
					self.nodeHash[nodeID].setGPACPose(guessPose)
					
					self.drawConstraints(self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull()

					hasMoved[nodeID] = True

			self.paths.generatePaths()

			self.drawConstraints(self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull()
			self.drawTrimmedPaths(self.paths.trimmedPaths)


			for nodeID in mergeNodes2:
	
				hull, medial = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)
				estPose = self.nodeHash[nodeID].getGlobalGPACPose()
	
				orderedPathIDs = [parentID1, pathID1]
				splicedPaths1 = self.paths.splicePathIDs(orderedPathIDs)	
				print "received", len(splicedPaths1), "spliced paths from path IDs", orderedPathIDs
				" departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum, angDiff2 |"
		

				" exclusion not required since paths are siblings "
				try:
					self.consistentFit(nodeID, estPose)
				except IndexError:
					print "failed to consistentFit node", nodeID
					pass
				#self.consistentFit(nodeID, estPose, excludePathIDs = [pathID2])
		

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				
				"""
				results1 = []		
				for k in range(len(splicedPaths1)):
					path = splicedPaths1[k]			
					result = getMultiDeparturePoint(path, medial, estPose, estPose, orderedPathIDs, nodeID)
					results1.append(result+(k,))
		
		
				results1 = sorted(results1, key=itemgetter(12))
				results1 = sorted(results1, key=itemgetter(10), reverse=True)					
		
				kIndex = results1[0][13]
		
				fitSplice = splicedPaths1[kIndex]
				"""
	
	
	
	
			for nodeID in mergeNodes2:

				print "adding node", nodeID, "to path", pathID1
				self.paths.addNode(nodeID, pathID1)
				print "deleting node", nodeID, "from path", pathID2
				self.paths.delNode(nodeID, pathID2)

			print "A: self.paths.pathClasses[" + str(pathID1) + "] = " + repr(self.paths.pathClasses[pathID1])
			print "A: self.paths.pathClasses[" + str(pathID2) + "] = " + repr(self.paths.pathClasses[pathID2])


			" move nodes to pathID1, delete pathID2 "
			self.paths.delPath(pathID2, pathID1)

			print "B: self.paths.pathClasses[" + str(pathID1) + "] = " + repr(self.paths.pathClasses[pathID1])
			
			self.paths.generatePaths()

			print "C: self.paths.pathClasses[" + str(pathID1) + "] = " + repr(self.paths.pathClasses[pathID1])

			self.drawConstraints(self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull()
			self.drawTrimmedPaths(self.paths.trimmedPaths)


			""" Now we check to see if any branch events have occurred from this merge """
			" departure events for node 1 "
			for nodeID1 in mergeNodes2:
				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
				contig1 = []		

				orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
				#orderedPathIDs1.insert(0,0)
				#orderedPathIDs1.append(0)

				#orderedPathIDs1.append(0)
				print nodeID1, "orderedPathIDs1:", orderedPathIDs1

				" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.paths.getDeparturePoint(self.trimmedPaths[pathID], nodeID1, plotIter = True)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
					contig1.append((contigFrac, overlapSum))

				print "node pathIDs   ", nodeID1, ":", orderedPathIDs1
				print "node departures", nodeID1, ":", departures1
				print "node interiors ", nodeID1, ":", interiors1
				print "node depPoints ", nodeID1, ":", depPoints1
				print "node distances ", nodeID1, ":", distances1
				print "node depAngles ", nodeID1, ":", depAngles1
				print "node contiguity", nodeID1, ":", contig1
					
				" new junction finding logic "
				" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
				" if a terminal departure exists that is internal, than we have a new junction "
				DISC_THRESH = 0.2

				" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
				frontExist1 = departures1[0][0]
				backExist1 =  departures1[-1][1]
				frontInterior1 = interiors1[0][0]
				backInterior1 = interiors1[-1][1]
				
				depAngle1 = depAngles1[0][0]
				depPoint1 = depPoints1[0][0]
				parentPathID1 = orderedPathIDs1[0]

				foreTerm1 = frontInterior1 and frontExist1
				backTerm1 = backInterior1 and backExist1
				
				" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
				discForeTerm1 = distances1[0][0] > DISC_THRESH
				discBackTerm1 = distances1[-1][1] > DISC_THRESH
				
				foreAngle1 = depAngles1[0][0]
				backAngle1 = depAngles1[-1][1]
				
				isFront1 = self.nodeHash[nodeID1].faceDir
				if isFront1:
					dirFlag = 0
				elif not isFront1:
					dirFlag = 1
				else:
					print isFront1
					raise

				isBranch, pathBranchID, isNew = self.paths.determineBranchSingle(nodeID1, frontExist1, frontInterior1, depAngle1, depPoint1, parentPathID1, dirFlag)
				
				print "determineBranchSingle:"
				print isBranch
				print pathBranchID

				isNew1 = isNew

				if isBranch:
					orderedPathIDs1.insert(0,pathBranchID)
					departures1.insert(0, [False, False])
					interiors1.insert(0, [False, False])
					depPoints1.insert(0, [None, None])

				backExist1 = departures1[-1][1]
				backInterior1 = interiors1[-1][1]
				depAngle1 = depAngles1[-1][1]
				depPoint1 = depPoints1[-1][1]
				parentPathID1 = orderedPathIDs1[-1]

				isFront1 = self.nodeHash[nodeID1].faceDir
				if isFront1:
					dirFlag = 1
				elif not isFront1:
					dirFlag = 0
				else:
					print isFront1, isFront2
					raise
				
				isBranch, pathBranchID, isNew = self.paths.determineBranchSingle(nodeID1, backExist1, backInterior1, depAngle1, depPoint1, parentPathID1, dirFlag)
				print "determineBranchSingle:"
				print isBranch
				print pathBranchID


				isNew1 = isNew1 or isNew

				if isBranch:
					orderedPathIDs1.append(pathBranchID)
					departures1.append([False, False])
					interiors1.append([False, False])
					depPoints1.append([None, None])

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

				self.paths.generatePaths()
				self.trimmedPaths = self.paths.trimPaths(self.paths.paths)						

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				self.drawTrimmedPaths(self.paths.trimmedPaths)

			#for k, v in hasMoved.iteritems():
			#	 pass
			
	

	@logFunction
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
		toBeSibMerged = []
		
		delIndex = []
		for k in range(len(toBeMerged)):
			
			mergeThis = toBeMerged[k]

			if len(mergeThis) > 4:
				delIndex.append(k)
				toBeSibMerged.append(mergeThis)
		
		delIndex.reverse()
		for index in delIndex:
			toBeMerged.pop(index)

		
		sorted(toBeSibMerged, key=itemgetter(5))
		
		for mergeThis in toBeSibMerged:
			print "mergeThis:", mergeThis[0], mergeThis[1], mergeThis[2], mergeThis[4], mergeThis[5], mergeThis[6]				

			pathID1 = mergeThis[0]
			pathID2 = mergeThis[1]
			offset = mergeThis[2]
			pathSplices = mergeThis[3]
			offset2 = mergeThis[4]
			lastCost = mergeThis[5]
			matchCount = mergeThis[6]
			
			allPathIDs = self.paths.getPathIDs()
			if pathID1 in allPathIDs and pathID2 in allPathIDs:
				
				self.mergeSiblings(offset, offset2, pathID1, pathID2, mergeThis[5], mergeThis[6])			
			
		for mergeThis in toBeMerged:		

			print "mergeThis:", mergeThis[0], mergeThis[1], mergeThis[2]


			pathID1 = mergeThis[0]
			pathID2 = mergeThis[1]
			offset = mergeThis[2]
			pathSplices = mergeThis[3]
			
			" verify that this path still exists before we try to merge it "
			allPathIDs = self.paths.getPathIDs()
			if pathID1 in allPathIDs and pathID2 in allPathIDs:

				
				" capture transform between pathID1 and pathID2 under correction "

				cOffset = offset
				rootID = pathID1
				lessID = pathID2
				
				
				mergeNodes = copy(self.paths.getNodes(lessID))
				
				
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
					
					print "add path constraints:"
					print "pathNodes:", pathNodes
					print "targetNodeID:", targetNodeID
					
					nodeHasMerged[nodeID] = False
			

				

					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID].getGlobalGPACPose()[:2]					
					print "adjusting pose", nodeID, "from", self.nodeHash[nodeID].getGlobalGPACPose()
					
					" if siblings, include mutual parent in splice"
					#if 

					#splicedPaths1 = self.paths.splicePathIDs([rootID])

					try:
						self.consistentFit(nodeID, self.nodeHash[nodeID].getGlobalGPACPose(), excludePathIDs = [lessID])
					except IndexError:
						print "failed to consistentFit node", nodeID
						pass

					self.drawConstraints(self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull()

					
					
					self.paths.addNode(nodeID, rootID)
					self.paths.delNode(nodeID, lessID)
					
					if len(self.paths.getNodes(lessID)) == 0:
						self.paths.delPath(lessID, rootID)

					
					
					self.paths.generatePaths()
					

					self.drawConstraints(self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull()


				print "relaxed constraints",  self.statePlotCount
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				
				" constrain nodes to merge to pathID1 "
				
				" move nodes to pathID1, delete pathID2 "
				self.paths.delPath(lessID, rootID)
	
				self.paths.generatePaths()
				
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
			

	@logFunction
	def consistentFit(self, nodeID, estPose, numGuesses = 11, excludePathIDs = []):

		splicedPaths1, spliceTerms, splicePathIDs = self.paths.getSplicesByNearJunction(nodeID)


		print "consistentFit(", nodeID

		initGuesses = []
		orientedPaths = []
		
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

		resultsBySplice = []

		for spliceIndex in range(len(splicedPaths1)):
			
			path = splicedPaths1[spliceIndex]

	
			pathSpline = SplineFit(path, smooth=0.1)
			pathU1 = pathSpline.findU(estPose1[:2])
			pathForeU = pathSpline.getUOfDist(originU1, 1.0, distIter = 0.001)
			pathBackU = pathSpline.getUOfDist(originU1, -1.0, distIter = 0.001)
			
			orientedPath = orientPath(path, globalMedial)
			orientedPathSpline = SplineFit(orientedPath, smooth=0.1)
	
			#medialSpline1 = SplineFit(medial1, smooth=0.1)
	
	
			globalSamples = orientedPathSpline.getUniformSamples(spacing = 0.04)
			medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
			
			globalMedialSamples = []
			for p in medialSamples:
				result = poseOrigin.convertLocalOffsetToGlobal(p)	
				globalMedialSamples.append(result)
			
			
			globalVar = computePathAngleVariance(globalSamples)
			medialVar = computePathAngleVariance(globalMedialSamples)

			"""
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
			""" 
	
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
			

			time1 = time.time()
	
			" add guesses in the neighborhood of the closest pathU "
			" 5 forward, 5 backward "
			
	
			" while loop "
			halfNum = numGuesses/2
			dists = [(k-halfNum)*0.2 for k in range(numGuesses)]	
			for dist in dists:
				
				
				#pathUGuess = i*0.05
				pathUGuess = orientedPathSpline.getUOfDist(pathU1, dist)
				initGuesses.append([pathUGuess, originU2, 0.0, spliceIndex])
				
			orientedPaths.append(orientedPath)

		print len(orientedPaths), "orientedPaths"
		print len(initGuesses), "initGuesses"
		results = batchGlobalMultiFit(initGuesses, orientedPaths, medial1, estPose1, [], nodeID)

		time2 = time.time()
		print "time:", time2 - time1

		resultsBySplice += results
		
		print "consistentFit node", nodeID
		
		print len(resultsBySplice), "results"
		
		for k in range(len(resultsBySplice)):
			result = resultsBySplice[k]
			resultPose = result[0]
			dist = sqrt((estPose1[0]-resultPose[0])**2 + (estPose1[1] - resultPose[1])**2)
			resultsBySplice[k] = result + (dist, k,)


		print "unfilteredResults:"
		for result in resultsBySplice:
			resultPose = result[0]
			isInterior1 = result[5]
			isExist1 = result[6]
			isInterior2 = result[11]
			isExist2 = result[12]
			contigFrac = result[15]
			angDiff2 = result[17]
			spliceIndex = result[19]
			dist = sqrt((estPose1[0]-resultPose[0])**2 + (estPose1[1] - resultPose[1])**2)
			lastCost = result[1]
			matchCount = result[2]
			overlapSum = result[16]
			print spliceIndex, [int(isInterior1), int(isInterior2), int(isExist1), int(isExist2)], contigFrac, dist, angDiff2, lastCost, matchCount, overlapSum, spliceTerms[spliceIndex], splicePathIDs[spliceIndex], resultPose

		
		" resultPose, lastCost, matchCount, departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2, isDeparture, spliceIndex "


		"""
		1) departures:	(0,0) good, (1,0) okay, (1,1) bad 
		2) contigFrac:	1.0 good, 0.9 okay, 0.5 bad, 0.0 reject
		3) dist:   0.0 great, 0.5 good, 1.0 okay, 2.0 bad, 3.0 reject
		4) angDiff2:  0.0 great, 0.2 good, 0.5 bad/reject
		5) splicePathIDs:  exID in REJECT, currPathID in GOOD
		"""
		
		#(0.5-angDiff2)/0.5 + (contigFrac-0.5)/0.5 + (3.0-dist)/3.0 + 1-isDeparture*0.5
		
		filteredResults = []
		for result in resultsBySplice:
			resultPose = result[0]
			isInterior1 = result[5]
			isExist1 = result[6]
			isInterior2 = result[11]
			isExist2 = result[12]
			contigFrac = result[15]
			angDiff2 = result[17]
			isDeparture = result[18]
			spliceIndex = result[19]
			dist = result[20]
			
			
			" rejection filter "
			isReject = False
			
			if isInterior1 or isInterior2:
				isReject = True
			

			if isExist1 and isExist2:
				isReject = True
			
			for exID in excludePathIDs:
				if exID in splicePathIDs[spliceIndex]:
					isReject = True
		
			if contigFrac <= 0.5:
				isReject = True
				
			if dist >= 3.0:
				isReject = True
			
			if angDiff2 > 0.7:
				isReject = True
				
			if not isReject:
				filteredResults.append(result)

			"""			
			if isExist1 or isExist2:
				depVal = 50.0
			else:
				depVal = 100.0
			
			
			if not isInterior1 and not isInterior2 and not isExist1 and not isExist2:
				exclude = False
				for exID in excludePathIDs:
					if exID in splicePathIDs[spliceIndex]:
						exclude = True
				
				if dist < 3.0 and contigFrac > 0.7 and angDiff2 < 0.5 and not exclude:
					filteredResults.append(result)
			"""
			

		#for result in filteredResults:
		#	resultPose = result[0]
		#	spliceIndex = result[18]
			#getMultiDeparturePoint(orientedPaths[spliceIndex], medial1, estPose1, resultPose, [], nodeID, pathPlotCount = self.multiDepCount, plotIter = True)
			#self.multiDepCount += 1
		
		print "filteredResults:"
		for result in filteredResults:
			resultPose = result[0]
			isInterior1 = result[5]
			isExist1 = result[6]
			isInterior2 = result[11]
			isExist2 = result[12]
			contigFrac = result[15]
			angDiff2 = result[17]
			spliceIndex = result[19]
			dist = result[20]
			lastCost = result[1]
			matchCount = result[2]
			overlapSum = result[16]
			resultIndex = result[21]
			print resultIndex, spliceIndex, [int(isInterior1), int(isInterior2), int(isExist1), int(isExist2)], contigFrac, dist, angDiff2, lastCost, matchCount, overlapSum, spliceTerms[spliceIndex], splicePathIDs[spliceIndex], resultPose



		utilResults = []
		for k in range(len(filteredResults)):
			result = filteredResults[k]
			angDiff2 = result[17]
			contigFrac = result[15]
			dist = result[20]
			isDeparture = result[18]
			utilVal = (0.5-angDiff2)/0.5 + (contigFrac-0.5)/0.5 + (3.0-dist)/3.0 + 1-isDeparture*0.5
			
			filteredResults[k] = result + (utilVal,)
			

		sortedResults = sorted(filteredResults, key=itemgetter(15), reverse=True)
		contigSort = []
		for k in range(len(sortedResults)):
			result = sortedResults[k]
			contigSort.append((result[15],result[21]))
		
		sortedResults = sorted(filteredResults, key=itemgetter(20))
		distSort = []
		for k in range(len(sortedResults)):
			result = sortedResults[k]
			distSort.append((result[20],result[21]))
		
		sortedResults = sorted(filteredResults, key=itemgetter(17))
		angSort = []
		for k in range(len(sortedResults)):
			result = sortedResults[k]
			angSort.append((result[17],result[21]))
		
		sortedResults = sorted(filteredResults, key=itemgetter(18))
		depSort = []
		for k in range(len(sortedResults)):
			result = sortedResults[k]
			depSort.append((result[18],result[21]))


		sortedResults = sorted(filteredResults, key=itemgetter(22), reverse = True)
		utilSort = []
		for k in range(len(sortedResults)):
			result = sortedResults[k]
			utilSort.append((result[22],result[21]))
		
		
		print "contigSort:", contigSort
		print "distSort:", distSort
		print "angSort:", angSort
		print "depSort:", depSort
		print "utilSort:", utilSort


		#filteredResults = sorted(filteredResults, key=itemgetter(16))
		filteredResults = sorted(filteredResults, key=itemgetter(22), reverse=True)

		guessPose = filteredResults[0][0]
		self.nodeHash[nodeID].setGPACPose(guessPose)

	
	@logFunction
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


	@logFunction
	def makeGlobalMedialOverlapConstraint(self, nodeID, globalPath, globalJunctionPoint, globalDeparturePoint, fullOverlap = False ):

		" compute the medial axis for each pose "
		
		node1 = self.nodeHash[nodeID]
		posture1 = node1.getStableGPACPosture()
		#hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = True)
		hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = False)
		
		

		
		estPose1 = node1.getGlobalGPACPose()		
		poseOrigin = Pose(estPose1)
		
		globalMedial = []
		for p in medial1:
			globalMedial.append(poseOrigin.convertLocalToGlobal(p))
		
		medialSpline1 = SplineFit(medial1, smooth=0.1)

		orientedGlobalPath = orientPath(globalPath, globalMedial)
		globalSpline = SplineFit(orientedGlobalPath, smooth=0.1)


		globalSamples = globalSpline.getUniformSamples(spacing = 0.04)
		medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
		
		globalMedialSamples = []
		for p in medialSamples:
			result = poseOrigin.convertLocalOffsetToGlobal(p)	
			globalMedialSamples.append(result)
		

		globalVar = computePathAngleVariance(globalSamples)
		medialVar = computePathAngleVariance(globalMedialSamples)
		
		"""
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
		"""


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
		#s = sorted(student_objects, key=attrgetter('age'))		# sort on secondary key
		#sorted(s, key=attrgetter('grade'), reverse=True)  
		#closestPairs.sort()

		print len(closestPairs), "closest pairs"
		#for val in closestPairs:
		#	print val


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
		
		
		if fullOverlap:
			resultPose, lastCost, matchCount = gen_icp.globalOverlapICP_GPU2([u1,u2,angGuess], orientedGlobalPath, medial1, plotIter = True, n1 = nodeID, n2 = -1, minMatchDist = 2.0)
		else:
			resultPose, lastCost, matchCount = gen_icp.globalOverlapICP_GPU2([u1,u2,angGuess], orientedGlobalPath, medial1, plotIter = True, n1 = nodeID, n2 = -1)
		


		print "estimating branch pose", nodeID, "at",  resultPose[0], resultPose[1], resultPose[2]

		return resultPose, lastCost


	@logFunction
	def makeMultiJunctionMedialOverlapConstraint(self, nodeID1, nodeID2, isMove = True, isForward = True, inPlace = False, uRange = 1.5):

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
						medial2Reverse = deepcopy(medial2)
						medial2Reverse.reverse()
						orientedMedial2 = deepcopy(medial2Reverse)
					else:
						orientedMedial2 = deepcopy(medial2)
					
					#orientedMedial2 = orientPath(medial2, medial1)
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
					
					result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)
						
					transform = matrix([[result[0]], [result[1]], [result[2]]])
					covE =	self.E_overlap
					
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
						medial2Reverse = deepcopy(medial2)
						medial2Reverse.reverse()
						orientedMedial2 = deepcopy(medial2Reverse)
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
					covE =	self.E_overlap
					
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

		#print "Sorted:"
		#for guess in totalGuesses:
		#	print guess

		transform = results[selectedIndex][7]
		covE = results[selectedIndex][8]
		hist = results[selectedIndex][9]
		
		return transform, covE, hist


	@logFunction
	def makeMedialOverlapConstraint(self, i, j, isMove = True, isForward = True, inPlace = False, uRange = 0.1 ):

		#print "recomputing hulls and medial axis"
		" compute the medial axis for each pose "
		
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

		#distEst = 0.3
		distEst = 0.5
		#distEst = 0.8
		#distEst = 1.0

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
						u2 = medialSpline2.getUOfDist(originU2, distEst, distIter = 0.001)
				else:	
					u2 = medialSpline2.getUOfDist(originU2, distEst, distIter = 0.001)
					
			else:

				#u2 = medialSpline2.getUOfDist(originU2, -0.3, distIter = 0.001)

				if len(node2.backProbeError) > 0:
	
					backSum = 0.0
					backProbeError = node2.backProbeError
					for n in backProbeError:
						backSum += n
					backAvg = backSum / len(backProbeError)
									
					if backAvg >= 1.4:
						u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
					else:	
						u2 = medialSpline2.getUOfDist(originU2, -distEst, distIter = 0.001)
				else:	
					u2 = medialSpline2.getUOfDist(originU2, -distEst, distIter = 0.001)
				
				
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


		result, hist = gen_icp.overlapICP_GPU2(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)

		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE =	self.E_overlap
		
		print "making overlap constraint:", result[0], result[1], result[2]

		return transform, covE, hist
	
	@logFunction
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

	@logFunction
	def getNearestPathPoint(self, originPoint):

		newPoint = self.paths.getNearestPathPoint(originPoint)

		return newPoint

	@logFunction
	def computeNavigationPath(self, startPose, endPose):
		return self.paths.computeNavigationPath(startPose, endPose)


	@logFunction
	def plotEnv(self):
		
		#walls = self.probe.getWalls()
		walls = self.walls
	
		for wall in walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = copy(wall[i])
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g', zorder=0)
			
	
	@logFunction
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


		xP = []
		yP = []
		for pose in poses:
			xP.append(pose[0])
			yP.append(pose[1])		
		pylab.scatter(xP,yP, color='k', linewidth=1, zorder=10)

		self.plotEnv()

		#pylab.title("Corner Constraint %d ---> %d" % (id1, id2))
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		pylab.title("%d Poses" % self.numNodes)
		pylab.xlim(-10,10)
		pylab.ylim(-10,10)

			
		if id == []:
			print "plotEstimate:", self.numNodes
			printStack()
			pylab.savefig("plotEstimate%04u.png" % self.numNodes)
		else:
			print "plotEstimate:", id
			printStack()
			pylab.savefig("plotEstimate%04u.png" % id)

	@logFunction
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

			globJuncPose = self.paths.getGlobalJunctionPose(k)
			if globJuncPose != None:
				pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='k')
			
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		self.plotEnv()


		print "saving trimmedPath_%04u.png" % self.trimCount
		pylab.title("Trimmed Paths, numNodes = %d" % self.numNodes)
		pylab.savefig("trimmedPath_%04u.png" % self.trimCount)
		self.trimCount += 1


	@logFunction
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

		for k in pathIDs:

			globJuncPose = self.paths.getGlobalJunctionPose(k)
			if globJuncPose != None:
				pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='k', zorder=10)

		self.plotEnv()
		
		print "pathAndHull:", self.pathDrawCount
		printStack()

		pylab.axis("equal")
		#pylab.xlim(-10, 12)
		#pylab.ylim(-10, 10)
		pylab.title("paths: %s numNodes: %d %d" % (repr(self.paths.getPathIDs()), self.numNodes, highestNodeID))
		pylab.savefig("pathAndHull_%04u.png" % self.pathDrawCount)

		self.pathDrawCount += 1
			

