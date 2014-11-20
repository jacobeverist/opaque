

from Pose import Pose
import math
import ctypes, os, sys
import multiprocessing as processing
from copy import copy, deepcopy
from functions import *
import random
from subprocess import Popen, PIPE
from PIL import Image
from medialaxis import computeMedialAxis
import graph
from LocalNode import getLongestPath
import gen_icp
from SplineFit import SplineFit
import pylab
import numpy
from operator import itemgetter
import hashlib
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath, getTipAngles, orientPathLean
#from MapProcess import selectLocalCommonOrigin, selectCommonOrigin
from ParticleFilter import multiParticleFitSplice, batchLocalizeParticle, batchDisplaceParticles, Particle
import time
import traceback

import alphamod
from itertools import product

from shoots import computeShootSkeleton, spliceSkeletons, computeGlobalControlPoses, batchJointBranch, batchBranch, trimBranch, getBranchPoint, ensureEnoughPoints, getInitSkeletonBranchPoint, getSkeletonBranchPoint
from landmarks import *
#from shoots import *

pylab.ioff()


" 1) probability that a path is the same as another path "
" 2) new ID for each new path "
" 3) test fitness of being its own path and being the same as other path "
" 4) collapse path to old path if pass threshold "

"""
Path Data Structure

parentID: int, the parent path ID
branchNodeID: int, the node from which the branching point is computed
localJunctionPose: [float*3], the pose of the branching point in local coordinates of branchNodeID
sameProb: dict of floats for each pathID indicating probability this path is the same as older path 

"""

" Likeness Test "
" 1) % overlap of two paths via ICP "
" 2) successful overlap of member nodes translated to old path "
" 2a) maintain separate set of intra-path constraints for each path "
" 2b) measure X_theta error for each path hypothesis "

" New Path Condition "
" 1) any new departure automatically creates a new path "
" 2) a new departure that does not overlap correctly creates a new path "


" Overlap State "
" 1) perform best overlap of each splice "
" 2) perform ICP of each splice "
" 3) maintain path-state machine to keep track of overlaps and directions and likely configurations "
" 4) extract path-state from particle filter motion estimation and subsequent ICP fit on splice "


" FIXME:  Branch direction does not seem to be used properly yet.  Only used in getOverlapDeparture() "

"""

Things that are happening with the data set with branches 
-----------

Second branch detection fails because is deemed "similar" in departure angle to path 1.
The difference is about 1.01 radians but is just under the difference threshold.
Paradoxically, this node is not added to the similar path and is instead added to path 0.
This is because the similar path 1 is not part of the orderedOverlap set.

When the next instance of branch detection occurs, it is determined to be successful because of both dist and ang diff.
However, this new branch eventually is constrained to the junction where it actually belongs and on top of the node 
with the previous branch detection fail.  The mechanism of how this happens is?   NOTE:  This happened because of a bug in my code
that allowed cross-path constraints between member nodes (node path membership always returned True).  Remarkably, there were no map errors from this permissiveness.

In the event that the first node is successfully branch detected, subsequent branch nodes are added to existing path,
but not properly merged or constrained.

"""

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

			
" FIXME: branch from set of ordered overlap pathIDs (0,2) selects 2 as the parent, but should really be 0 "
" How did it overlap with 2 in the first place?  It departed from 0, but selected 2 as parent since it's the terminating pathID "
" 2 was marked positive by overlapCondition since it was parallel to 2 within the minMatchDist of 1.0 "

@logFunction
def computePathAngleVariance(pathSamples):
	pathVar = []
	
	" compute the local variance of the angle "
	VAR_WIDTH = 40
	for i in range(len(pathSamples)):
		
		lowK = i - VAR_WIDTH/2
		if lowK < 0:
			lowK = 0
			
		highK = i + VAR_WIDTH/2
		if highK >= len(pathSamples):
			highK = len(pathSamples)-1
		
		localSamp = []
		for k in range(lowK, highK+1):
			localSamp.append(pathSamples[k][2])
		
		sum = 0.0
		for val in localSamp:
			sum += val
			
		meanSamp = sum / float(len(localSamp))
		

		sum = 0
		for k in range(len(localSamp)):
			sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
	
		varSamp = sum / float(len(localSamp))
		
		pathVar.append((meanSamp, varSamp))		 

	return pathVar


@logFunction
def getOverlapCondition(medial2, estPose2, supportLine, nodeID, plotIter = False, overlapPlotCount = 0):

	if len(supportLine) == 0:
		return 1e100

	#medial2 = self.poseData.medialAxes[nodeID]

	#estPose2 = self.nodePoses[nodeID]
	
	minMatchDist2 = 0.2
		
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
	
	support_pairs = []
	for i in range(len(vecPoints2)):
		
		p_2 = vecPoints2[i]

		" for every transformed point of A, find it's closest neighbor in B "
		try:
			p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)

			if minDist <= minMatchDist2:
				C2 = points2[i][2]
				C1 = supportPoints[minI][2]

				" we store the untransformed point, but the transformed covariance of the A point "
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
		pylab.savefig("overlapCost_%04u.png" % overlapPlotCount)
		#overlapPlotCount += 1
	
	if len(support_pairs) == 0:
		return 1e100

	return cost




class MapState:
	
	def __init__(self, poseData, hypothesisID, args = None):
		
		""" raw sensor data products """
		self.poseData = deepcopy(poseData)
		
		""" pose particles for different combinations of poses and branch points """
		#self.numPoseParticles = 100
		#self.numPoseParticles = 40

		if args:
			self.numPoseParticles = args.numPoseParticles
		else:
			self.numPoseParticles = 40


		#self.numPoseParticles = 1
		self.poseParticles = {}
		self.poseParticles["numParticles"] = self.numPoseParticles
		self.poseParticles["updateCount"] = 0
		self.poseParticles["snapshots2"] = {0 : ()}

		self.parentParticle = None

		""" max index of particle """
		self.currMaxIndex = 0


		""" pose information for each of the spatial map nodes """
		self.isNodeBranching = {}
		self.nodeRawPoses = {}
		self.nodePoses = {}
		self.gndPoses = {}
		self.gndRawPoses = {}
		self.origPoses = {}

		""" evaluation result of map consistency """
		self.utility = 0.0

		""" if a branching decision was made, and this is the no-branch state """
		self.isNotBranched = False

		""" criteria for rejecting this map state """
		self.mapOverlapSum = 0.0
		self.isNoLocalize = False

		""" unique ID for this MapState instantation """
		self.hypothesisID = hypothesisID

		""" the alpha hull and maximum length path of shoot skeleton """
		self.paths = {0 : []}
		self.hulls = {0 : []}
		self.localPaths = {0 : []}
		self.localHulls = {0 : []}
		self.localLandmarks = {0 : {}}
		self.branchDiverges = {0 : True}

		""" intermediate data structure for computing the shoot skeleton """
		self.localLeaf2LeafPathJunctions = {}

		""" shoot spine only, separated from the skeleton at the branch point """
		self.trimmedPaths  = {}

		""" membership, parent-relationship, branch point, and control point for a shoot """
		self.pathClasses = {}
		self.pathClasses[0] = {"parentID" : None,
					"branchNodeID" : None,
					"localJunctionPose" : None, 
					"localDivergencePose" : None, 
					"sameProb" : {},
					"nodeSet" : [],
					"localNodePoses" : {},
					"globalJunctionPose" : None,
					"controlPose" : [0.0,0.0,0.0] }		

		""" number of shoots in the current map state """
		self.pathIDs = 1

		""" travel destinations on the shoot map """
		self.pathTermsVisited = {-1: False, 0: False}

		""" number of branch hypotheses and spacing between them for a branch point distribution """
		self.DIV_LEN = 0.2
		self.NUM_BRANCHES = 5
		#self.DIV_LEN = 0.1
		#self.NUM_BRANCHES = 1


		""" the entrance point of the environment """
		self.rootPoint = [-3.2, 0.0]
		
		""" data structures for computing splices on the shoot map """
		self.pathGraph = graph.graph()
		self.joins = []
		self.junctions = {}
		self.terminals = {}
		self.allSplices = {}
		self.isChanged = True

		self.shootIDs = 1

		""" plotting colors """
		self.colors = []

		self.colors.append([127, 127, 255])
		self.colors.append([127, 255, 127])
		self.colors.append([255, 127, 127])
		self.colors.append([255, 127, 255])
		self.colors.append([255, 255, 127])
		self.colors.append([127, 255, 255])

		for color in self.colors:
			color[0] = float(color[0])/256.0
			color[1] = float(color[1])/256.0
			color[2] = float(color[2])/256.0

		for i in range(1000):
			self.colors.append((random.random(),random.random(),random.random()))


		""" plot counts """
		self.alphaPlotCount = 0
		self.medialCount = 0
		self.topCount = 0
		self.spliceCount = 0
		self.multiDepCount = 0
		self.overlapPlotCount = 0
		self.tempCount = 0
		self.pathPlotCount = 0
		self.pathPlotCount2 = 0
		self.overlapPlotCount2 = 0
		self.medialSoupCount = 0
		self.posePlotCount = 0



	def getNodePose(self, nodeID):
		return copy(self.nodePoses[nodeID])

	def getNodeRawPose(self, nodeID):
		return copy(self.nodeRawPoses[nodeID])
	
	def setNodePose(self, nodeID, newPose):

		print "setNodePose(", nodeID, ",", newPose, ")"

		oldGPACPose = self.getNodePose(nodeID)

		self.nodePoses[nodeID] = newPose

		gpacProfile = Pose(oldGPACPose)
		
		localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[nodeID])
		
		" go back and convert this from GPAC pose to estPose "
		newProfile = Pose(newPose)
		newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
		
		self.nodeRawPoses[nodeID] = newEstPose


		allPathIDs = self.getPathIDs()

		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		particleDist2 = self.poseParticles["snapshots2"][0]

		for pathID in allPathIDs:

			if nodeID in self.pathClasses[pathID]["nodeSet"]:

				if pathID != 0:

					""" convert the controlPose to coordinates local to the parent frame """
					shootControlPose = globalControlPoses[pathID]
					currFrame = Pose(shootControlPose)
					localNodePose = currFrame.convertGlobalPoseToLocal(newPose)

					self.pathClasses[pathID]["localNodePoses"][nodeID] = localNodePose

					for part in particleDist2:
						part.junctionData[pathID]["localNodePoses"][nodeID] = localNodePose


				else:
					self.pathClasses[pathID]["localNodePoses"][nodeID] = newPose

					for part in particleDist2:
						part.junctionData[pathID]["localNodePoses"][nodeID] = newPose


	@logFunction
	def updateMaxParticle(self, maxIndex):

		print "updateMaxParticle(", maxIndex, ")"


		particleDist2 = self.poseParticles["snapshots2"][0]
		maxParticle = particleDist2[maxIndex]
		allPathIDs = self.getPathIDs()
		parentPathIDs = self.getParentHash()

		for pathID in allPathIDs:
			
			allNodes = self.pathClasses[pathID]["nodeSet"]
			for nodeID in allNodes:

				newPose = maxParticle.junctionData[pathID]["localNodePoses"][nodeID]
				self.pathClasses[pathID]["localNodePoses"][nodeID] = newPose


		""" change to the maximum likelihood branch position as well """

		""" maximum likelihood pose particle """
		maxParticle = particleDist2[maxIndex]
		partBranchTupleIndex = maxParticle.maxLikelihoodBranch
		branchArcDists = maxParticle.branchArcDists
		branchControls = maxParticle.branchControls

		branchPathIDs = deepcopy(allPathIDs)
		branchPathIDs.remove(0)
		branchPathIDs.sort()

		if partBranchTupleIndex != None:
			branchResult = self.jointBranchEvaluations[partBranchTupleIndex] 
		else:
			branchResult = None

		for k in range(len(branchPathIDs)):

			pathID = branchPathIDs[k]

			if partBranchTupleIndex != None:

				""" maximum likelihood branch point within maximum likelihood pose particle """
				""" update of canonical branch point from maximum likelihood branch point """
				arcDist = partBranchTupleIndex[k]
				controlPose_P = self.branchControlPoses[pathID][arcDist]

				origControlPose = copy(self.pathClasses[pathID]["controlPose"])
				print "changing path", pathID, "control position from", origControlPose, "to", controlPose_P

				self.pathClasses[pathID]["controlPose"] = controlPose_P

				branchPose_L = branchResult["branchPoses_L"][pathID]

				""" FIXME:  was recycling the angle, but now we are regenerating it """
				self.setLocalJunctionPose(pathID, branchPose_L)


		""" update the recent nodes if we are displaced """


		""" update the visible nodes, self.nodePoses """
		controlPoses = deepcopy(self.getControlPoses())
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		for pathID in allPathIDs:

			allNodes = self.pathClasses[pathID]["nodeSet"]

			print "local to global controlPose:", controlPoses[pathID], globalControlPoses[pathID]

			for nodeID in allNodes:

				shootControlPose_G = globalControlPoses[pathID]
				currFrame_G = Pose(shootControlPose_G)

				newPose_C = self.pathClasses[pathID]["localNodePoses"][nodeID]

				nodePose_G = currFrame_G.convertLocalOffsetToGlobal(newPose_C)

				""" get the relationship between current gpac and raw poses """
				oldGPACPose = self.getNodePose(nodeID)
				gpacProfile = Pose(oldGPACPose)
				localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[nodeID])
				
				" go back and convert this from GPAC pose to estPose "
				newProfile = Pose(nodePose_G)
				newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
				self.nodeRawPoses[nodeID] = newEstPose

				self.nodePoses[nodeID] = nodePose_G




	@logFunction
	def initializePoseParticles(self):

	
		pose0 = self.getNodePose(0)
		pose1 = self.getNodePose(1)

		path0 = deepcopy(self.paths[0])

		termPoint = self.terminals[0][1]

		termDist0 = sqrt((path0[0][0]-termPoint[0])**2 + (path0[0][1]-termPoint[1])**2)
		termDistN = sqrt((path0[-1][0]-termPoint[0])**2 + (path0[-1][1]-termPoint[1])**2)

		if termDist0 > termDistN: 
			path0.reverse()

		pathSpline = SplineFit(path0, smooth=0.1)

		minDist0, newU0, newP0 = pathSpline.findClosestPoint(pose0)
		minDist1, newU1, newP1 = pathSpline.findClosestPoint(pose1)

		dist0 = pathSpline.dist_u(newU0)
		dist1 = pathSpline.dist_u(newU1)

		#newAng0 = pathSpline.angle_u(newU0)
		#newAng1 = pathSpline.angle_u(newU1)
		newAng0 = pose0[2]
		newAng1 = pose1[2]

		totalLen = pathSpline.length()

		pathID = 0
		avgDist = (dist0 + dist1)/2.0
		spliceID = ('t0', 't1')

		part = (pathID, avgDist, ('t0','t1'))

		numParticles = self.poseParticles["numParticles"]
		initDist2 = []

		" create the initial particle distibution "
		while len(initDist2) < numParticles:
			hypDist = random.gauss(avgDist, 0.5)

			if hypDist > 0.0 and hypDist <= totalLen:
				hypPoint = pathSpline.getPointOfDist(hypDist)
				hypPose0 = copy(hypPoint)
				hypPose0[2] = newAng0
				hypPose1 = copy(hypPoint)
				hypPose1[2] = newAng1

				particleObj = Particle(hypPose0, hypPose1, pathID, hypDist, 0.0, self.hypothesisID)
				particleObj.spliceCurve = deepcopy(self.paths[0])
				particleObj.addNode(0,0, hypPose0)
				particleObj.addNode(1,0, hypPose1)
				initDist2.append(particleObj)

		self.poseParticles["snapshots2"][0] = initDist2

	@logFunction
	def batchDisplaceParticles(self, nodeID0, nodeID1):

		batchJobs = []
		pathIDs = self.getPathIDs()
		parentPathIDs = self.getParentHash()


		controlPoses = deepcopy(self.getControlPoses())
		controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

		candLandmarks_G = nodeToGlobalLandmarks(self.getControlPoses(), self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])

		targetNodeLandmarks_N = {nodeID0 : None, nodeID1 : None}
		targetNodeLandmarks_N[nodeID0] = getNodeLandmark(nodeID0, self.poseData)
		targetNodeLandmarks_N[nodeID1] = getNodeLandmark(nodeID1, self.poseData)

		""" collect landmarks that we can localize against """
		#candLandmarks_G = []
		#for pathID in pathIDs:

		#	nodeSet = self.getNodes(pathID)
		#	currFrame_L = Pose(controlPoses[pathID])

		#	for nodeID in nodeSet:

		#		landmarkPoint_N, pointThresh, pointName = self.nodeLandmarks[pathID][nodeID]

		#		currFrame_G = Pose(controlPoses_G[pathID])
		#		nodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
		#		poseFrame_L = Pose(nodePose_L)

		#		if landmarkPoint_N != None and nodeID != nodeID0 and nodeID != nodeID1:
		#			landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
		#			landmarkPoint_G = currFrame_G.convertLocalToGlobal(landmarkPoint_L)
		#			candLandmarks_G.append(landmarkPoint_G)

		particleDist2 = self.poseParticles["snapshots2"][0]
		for particleIndex in range(len(particleDist2)):

			part = particleDist2[particleIndex]

			hypPose0 = part.pose0
			hypPose1 = part.pose1
			prevPose0 = part.prevPose0
			prevPose1 = part.prevPose1

			staticSplicedPaths0, spliceTerms0, staticSplicePathIDs0 = self.getSplicesByNearJunction(hypPose0)
			staticSplicedPaths1, spliceTerms1, staticSplicePathIDs1 = self.getSplicesByNearJunction(hypPose1)



			batchJobs.append([self.poseData, part, particleIndex, nodeID1, prevPose0, prevPose1, hypPose0, hypPose1, self.paths[0], staticSplicedPaths0, staticSplicedPaths1, candLandmarks_G, targetNodeLandmarks_N])

		results = batchDisplaceParticles(batchJobs)

		for result in results:
			particleIndex = result[0]
			part = particleDist2[particleIndex]
			part.pose0 = result[1]
			part.pose1 = result[2]

			part.displacePose(result[1], result[2])

			part.displaceProb0 = result[3]
			part.displaceProb1 = result[4]


		#targetNodeLandmarks_N[nodeID0] = getNodeLandmark(nodeID0, self.poseData)
		#targetNodeLandmarks_N[nodeID1] = getNodeLandmark(nodeID1, self.poseData)

		""" evaluate the landmark consistency to select the max hypothesis """
		landmarkPoint0_N = targetNodeLandmarks_N[nodeID0]
		landmarkPoint1_N = targetNodeLandmarks_N[nodeID1]

		minSum = 1e100
		minPartIndex = self.currMaxIndex
		maxLandmarkSum = 0.0
		maxPartIndex = 0

		if landmarkPoint0_N != None or landmarkPoint1_N != None: 

			for partIndex in range(len(particleDist2)):

				part = particleDist2[partIndex]

				currLandmarks = deepcopy(candLandmarks_G)

				currPose0 = part.pose0
				currPose1 = part.pose1

				poseFrame0_G = Pose(currPose0)
				poseFrame1_G = Pose(currPose1)

				if landmarkPoint0_N != None:
					landmarkPoint0_G = (poseFrame0_G.convertLocalToGlobal(landmarkPoint0_N[0]), landmarkPoint0_N[1], landmarkPoint0_N[2])
					currLandmarks.append(landmarkPoint0_G)

					print "landmark0:", landmarkPoint0_G

				if landmarkPoint1_N != None:
					landmarkPoint1_G = (poseFrame1_G.convertLocalToGlobal(landmarkPoint1_N[0]), landmarkPoint1_N[1], landmarkPoint1_N[2])
					currLandmarks.append(landmarkPoint1_G)

					print "landmark1:", landmarkPoint1_G

				currSum = computeConsistency(currLandmarks)

				part.landmarkSum = currSum


				print "poseSum:", partIndex, currSum

				if currSum < minSum:
					minSum = currSum
					minPartIndex = partIndex

				if currSum > maxLandmarkSum:
					maxLandmarkSum = currSum
					maxPartIndex = partIndex

		maxProb = 0.0
		maxPart = 0

		for partIndex in range(len(particleDist2)):
			part = particleDist2[partIndex]
			currProb0 = part.displaceProb0
			currProb1 = part.displaceProb1

			poseLandmarkSum = part.landmarkSum

			""" check to avoid divide by zero """
			if maxLandmarkSum > 0.0:
				if maxLandmarkSum > poseLandmarkSum:
					#newProb0 = currProb0 * ((maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum)
					newProb0 = (maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum
				else:
					""" maximum landmark cost, means utility is zeroed """
					newProb0 = 0.0
			else:
				newProb0 = currProb0

			if newProb0 > maxProb:
				maxProb = newProb0
				maxPart = partIndex

		#print "setting max particle:", minPartIndex
		print "setting max particle:", maxPart, maxProb

		#self.currMaxIndex = minPartIndex
		self.currMaxIndex = maxPart

		""" change to the maximum likelihood branch position as well """
		self.updateMaxParticle(self.currMaxIndex)
		self.setNodePose(nodeID0, deepcopy(particleDist2[self.currMaxIndex].pose0))
		self.setNodePose(nodeID1, deepcopy(particleDist2[self.currMaxIndex].pose1))


	@logFunction
	def displacePoseParticles(self, nodeID0, nodeID1, travelDist0, travelDist1):

		#self.isChanged = True

		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]

		updateCount += 1
		self.poseParticles["updateCount"] = updateCount

		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]


		medialSpline = SplineFit(medial0, smooth=0.1)

		newPartDist2 = []
		for part in particleDist2:


			hypPose0 = part.pose0
			hypPose1 = part.pose1
			pathID = part.memberPaths[0]
			hypDist = part.hypDist
			spliceID = part.spliceName

			currSplice0 = part.spliceCurve

			poseOrigin0 = Pose(hypPose0)
			globalMedial0 = []
			for p in medial0:
				globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

			orientedSplicePath = orientPath(currSplice0, globalMedial0)				
			pathSpline = SplineFit(orientedSplicePath, smooth=0.1)


			minDist0, oldU0, oldP0 = pathSpline.findClosestPoint(hypPose0)
			minDist1, oldU1, oldP1 = pathSpline.findClosestPoint(hypPose1)

			thisDist0 = travelDist0
			thisDist1 = travelDist1

			moveChance = random.random()
			" 80% chance we move forward, 20% we move backward "    
			if moveChance >= 0.1:
				thisDist0 = travelDist0
				thisDist1 = travelDist1
			else:
				thisDist0 = -travelDist0
				thisDist1 = -travelDist1

			thisDist0 = random.gauss(thisDist0,0.6)
			thisDist1 = random.gauss(thisDist1,0.6)


			newU0 = pathSpline.getUOfDist(oldU0, thisDist0)
			newU1 = pathSpline.getUOfDist(oldU1, thisDist1)

			newDist0 = pathSpline.dist_u(newU0)
			newP0 = pathSpline.point_u(newU0)
			newAng0 = pathSpline.angle_u(newU0)
			oldAng0 = hypPose0[2]
			newPose0 = copy(newP0)

			newDist1 = pathSpline.dist_u(newU1)
			newP1 = pathSpline.point_u(newU1)
			newAng1 = pathSpline.angle_u(newU1)
			oldAng1 = hypPose1[2]
			newPose1 = copy(newP1)

			particleObj = part.copy()
			particleObj.displacePose(newPose0, newPose1)
			particleObj.newDist0 = newDist0

			newPartDist2.append(particleObj)


		self.poseParticles["snapshots2"][0] = newPartDist2

	#@logFunction
	def localizePoseParticles(self, nodeID0, nodeID1):

		#cProfile.run('runTest(probe)', 'test_prof')

		JOINT_DIST = 3.0

		poseData = self.poseData

		""" smoothed and vectored root path """
		rootSpline = SplineFit(self.paths[0], smooth=0.1)
		rootSplice = rootSpline.getUniformSamples()

		self.isNoLocalize = False
		self.resetBranches()

		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]

		nodePose0 = self.getNodePose(nodeID0)
		staticSplicedPaths, spliceTerms, staticSplicePathIDs = self.getSplicesByNearJunction(nodePose0)

		print "numStaticSplicedPaths =", len(staticSplicedPaths)

		pathIDs = self.getPathIDs()

		medial0 = poseData.medialAxes[nodeID0]
		medial1 = poseData.medialAxes[nodeID1]

		if nodeID0-2 >= 0:
			prevMedial0 = poseData.medialAxes[nodeID0-2]
			prevMedial1 = poseData.medialAxes[nodeID1-2]
			prevMedialSpline0 = SplineFit(prevMedial0, smooth=0.1)
			prevMedialSpline1 = SplineFit(prevMedial1, smooth=0.1)
			prevMedial0_vec = prevMedialSpline0.getUniformSamples()
			prevMedial1_vec = prevMedialSpline1.getUniformSamples()
		else:
			prevMedial0 = None
			prevMedial1 = None
			prevMedial0_vec = None
			prevMedial1_vec = None


		medialSpline0 = SplineFit(medial0, smooth=0.1)
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		minMedialDist0, oldMedialU0, oldMedialP0 = medialSpline0.findClosestPoint([0.0,0.0,0.0])
		minMedialDist1, oldMedialU1, oldMedialP1 = medialSpline1.findClosestPoint([0.0,0.0,0.0])

		medial0_vec = medialSpline0.getUniformSamples()
		medial1_vec = medialSpline1.getUniformSamples()




		""" FIXME:  convert particle to control matrix and switch to local controls """

		"""
		1.  First, for given pose, compute distance to each junction point for each branch point distribution
		2.  if branch point of shoot is beyond threshold distance, only include the most likely shoot point, all others same
		3.  For each branch point for each shoot, compute all possible combinations
		4.  Each combination is a problem set



		"""


		"""
		probSets = []
		for particleIndex in range(len(particleDist2)):
			part = particleDist2[particleIndex]
			hypPose0 = part.pose0
			pathIDs = self.getPathIDs()

			branchDists = {}

			for pathID in pathIDs:
				if pathID != 0:

					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]
					controlPose_L = part.junctionData[pathID]["controlPose"]

					branchPose_L = part.junctionData[pathID]["localJunctionPose"]

					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)

					branchDists[pathID] = dist1

					if dist1 < 3.0:
						probSets.append((pathID, branchPose_L, controlPose_L))
		"""


		if True:

			jointComboControlSets = []
			pathSplines = {}
			allPathIDs = self.getPathIDs()
			parentPathIDs = self.getParentHash()

			for pathID in allPathIDs:
				pathSpline = SplineFit(self.localPaths[pathID])
				pathSplines[pathID] = pathSpline

			branchPathIDs = deepcopy(allPathIDs)
			branchPathIDs.remove(0)
			branchPathIDs.sort()
			
			for particleIndex in range(len(particleDist2)):
				part = particleDist2[particleIndex]
				#part.maxLikelihoodBranch = None 
				hypPose0 = part.pose0


				origBranchDists = {}
				origBranchControls = {}
				origBranchArcDists = {}

				for pathID in branchPathIDs:

					parentID = parentPathIDs[pathID]
					pathSpline = pathSplines[parentID]
					controlPose = part.junctionData[pathID]["controlPose"]
					minDist, controlUVal, newControlPose = pathSpline.findClosestPoint(controlPose)
					arcDist = pathSpline.dist_u(controlUVal)

					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]
					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)

					origBranchDists[pathID] = dist1
					origBranchControls[pathID] = controlPose
					origBranchArcDists[pathID] = arcDist


				newBranchArcDists = {}
				newBranchControls = {}
				newBranchComboControls = {}


				if len(branchPathIDs) > 0:

					#self.branchControlPoses = {}
					#self.branchArcDists = {}


					""" bin the arc distances, and then compute discrete problem set """
					for pathID in branchPathIDs:

						minArcDiff = 1e100
						minIndex = 0

						origArcDist = origBranchArcDists[pathID]

						for k in range(len(self.branchArcDists[pathID])):

							currArcDist = self.branchArcDists[pathID][k]

							diff = abs(currArcDist-origArcDist)

							if diff < minArcDiff:
								minArcDiff = diff
								minIndex = k

						""" determine the range of values for this particular branch """
						newArcDist = self.branchArcDists[pathID][minIndex]
						if origBranchDists[pathID] < JOINT_DIST:

							if minIndex < 2:
								minIndex = 2

							if minIndex > len(self.branchArcDists[pathID]) - 3:
								minIndex = len(self.branchArcDists[pathID]) - 3


							arcList = []
							controlList = []
							comboList = []
							for k in range(minIndex-2,minIndex+3):
								newArcDist = self.branchArcDists[pathID][k]
								arcList.append(newArcDist)
								controlList.append(self.branchControlPoses[pathID][newArcDist])
								comboList.append((newArcDist, tuple(self.branchControlPoses[pathID][newArcDist])))

							newBranchArcDists[pathID] = arcList
							newBranchControls[pathID] = controlList
							newBranchComboControls[pathID] = comboList

						else:
							newBranchArcDists[pathID] = [newArcDist,]
							newBranchControls[pathID] = [self.branchControlPoses[pathID][newArcDist],]
							newBranchComboControls[pathID] = [(newArcDist, tuple(self.branchControlPoses[pathID][newArcDist])),]

					argSet2 = []
					for pathID in branchPathIDs:
						argSet2.append(newBranchComboControls[pathID])

					combIterator = product(*argSet2)
					for comb in combIterator:
						jointComboControlSets.append(comb)

				""" set the range of branches for this particular particle so we can reference them later """
				part.branchArcDists = newBranchArcDists
				part.branchControls = newBranchControls

		print "jointComboControl:", len(jointComboControlSets)
		jointComboControlSets = list(set(jointComboControlSets))
		print "jointComboControl:", len(jointComboControlSets)

		time1 = time.time()

		""" precompute the evaluation for branches and cache result """
		""" each branch point is the max likelihood of each pose particle """
		self.batchPrecomputeBranches(jointComboControlSets)

		time2 = time.time()
		print "TIME computeJointBranch", self.hypothesisID, "=", time2-time1

		parentPathIDs = self.getParentHash()
		localizeJobs = []
		allSplicedPaths = []
		spliceCount = 0


		""" collect landmarks that we can localize against for current nodes """
		landmark0_N = None
		landmark1_N = None
		#landmark0_L = None
		#landmark1_L = None

		for pathID in pathIDs:

			nodeSet = self.getNodes(pathID)
			for nodeID in nodeSet:

				landmarkPoint_N = self.nodeLandmarks[pathID][nodeID]

				nodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
				poseFrame_L = Pose(nodePose_L)

				if nodeID == nodeID0:
					landmark0_N = self.nodeLandmarks[pathID][nodeID]
					#if landmark0_N != None:
					#	landmark0_L = poseFrame_L.convertLocalToGlobal(landmark0_N)

				if nodeID == nodeID1:
					landmark1_N = self.nodeLandmarks[pathID][nodeID]
					#if landmark1_N != None:
					#	landmark1_L = poseFrame_L.convertLocalToGlobal(landmark1_N)


		for particleIndex in range(len(particleDist2)):

			#branchArcDists = part.branchArcDists
			#branchControls = part.branchControls

			#probDist = []
			#branchPoseDist = []
			#controlPoseDist = []
			#branchSplices = []
			#for branchSamp in samples:
			#	"""
			#	8 = initial probability value of this branch position, initialized to 0.0 in this function
			#	1 = pose of the branch point
			#	9 = set of splices including the branch at this position
			#	"""
			#	probDist.append(branchSamp["initProb"])
			#	branchPoseDist.append(branchSamp["modJuncPose"])
			#	controlPoseDist.append(branchSamp["modControlPose"])
			#	branchSplices.append(branchSamp["newSplices"])

			#	print "initProb, modJuncPose, modControlPose:", branchSamp["initProb"], branchSamp["modJuncPose"], branchSamp["modControlPose"]


			#probStr = ""
			#for probVal in probDist:
			#	probStr += "%1.2f " % probVal

			#print self.hypothesisID, particleIndex, pathID, "branch probDist:", probStr

			##self.poseParticles["snapshots2"][updateCount][particleIndex].junctionData[pathID] = {}
			#self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["probDist"] = probDist
			#self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["branchPoseDist"] = branchPoseDist
			#self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["controlPoseDist"] = controlPoseDist
			#self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["branchSplices"] = branchSplices


			time1 = time.time()

			part = particleDist2[particleIndex]

			hypPose0 = part.pose0
			hypPose1 = part.pose1
			hypDist = part.hypDist
			prevHypPose0 = part.prevPose0
			prevHypPose1 = part.prevPose1

			resultsBySplice = []

			thisSplicedPaths = []

			""" get control set of particle's branches indexed by arc distance """
			# part.branchArcDists = newBranchArcDists
			# part.branchControls = newBranchControls

			""" build indexing tuples for this particle """
			argSet = []
			for pathID in branchPathIDs:
				argSet.append(part.branchArcDists[pathID])

			arcIndexes = []
			if len(argSet) > 0:
				combIterator = product(*argSet)
				for comb in combIterator:
					arcIndexes.append(tuple(comb))


			# TODO: consider single path longestPaths , including root
			""" minimum 5 tuples """
			for arcTuple in arcIndexes:
				print "arcIndex:", arcTuple
				branchResult = self.jointBranchEvaluations[arcTuple]
				totalMatchCount = branchResult["totalMatchCount"]
				normMatchCount = branchResult["normMatchCount"]
				normLandmark = branchResult["normLandmark"]
				landmarks_G = branchResult["landmarks_G"]

				controlPoses = branchResult["controlSet"]
				controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)


				""" collect landmarks that we can localize against """
				#candLandmarks_G = []
				#for pathID in pathIDs:

				#	nodeSet = self.getNodes(pathID)
				#	currFrame_L = Pose(controlPoses[pathID])

				#	for nodeID in nodeSet:

				#		landmarkPoint_N = self.nodeLandmarks[pathID][nodeID]

				#		currFrame_G = Pose(controlPoses_G[pathID])
				#		nodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
				#		poseFrame_L = Pose(nodePose_L)

				#		if landmarkPoint_N != None and nodeID != nodeID0 and nodeID != nodeID1:
				#			landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
				#			landmarkPoint_G = currFrame_G.convertLocalToGlobal(landmarkPoint_L)
				#			candLandmarks_G.append(landmarkPoint_G)

				candLandmarks_G = nodeToGlobalLandmarks(controlPoses, self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])

				#LANDMARK_THRESH = 1e100
				#distSum = 0.0
				#for i in range(len(landmarks_G)):
				#	p1 = landmarks_G[i]
				#	for j in range(i+1, len(landmarks_G)):
				#		p2 = landmarks_G[j]
				#		dist = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

				#		if dist < LANDMARK_THRESH:
				#			distSum += dist


				splices_G = branchResult["splices_G"]
				for splice in splices_G:
					#thisSplicedPaths.append((arcTuple, normMatchCount*normLandmark, splice, None, []))
					thisSplicedPaths.append((arcTuple, normMatchCount*normLandmark, splice, None, candLandmarks_G))

				#splices_G = branchResult["splices_G"]
				#for pathID in branchPathIDs:

				#	""" 3 splices for this branch """
				#	numSplices = len(splices_G[pathID])

				#	for k in range(numSplices):
				#		splice = splices_G[pathID][k]
				#		thisSplicedPaths.append((arcTuple, normMatchCount, splice, pathID))

			#""" static splices have 1.0 probability """
			if len(thisSplicedPaths) == 0:

				controlPoses = deepcopy(self.getControlPoses())
				controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

				""" collect landmarks that we can localize against """
				#candLandmarks_G = []
				#for pathID in pathIDs:

				#	nodeSet = self.getNodes(pathID)
				#	currFrame_L = Pose(controlPoses[pathID])

				#	for nodeID in nodeSet:

				#		landmarkPoint_N = self.nodeLandmarks[pathID][nodeID]

				#		currFrame_G = Pose(controlPoses_G[pathID])
				#		nodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
				#		poseFrame_L = Pose(nodePose_L)

				#		if landmarkPoint_N != None and nodeID != nodeID0 and nodeID != nodeID1:
				#			landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
				#			landmarkPoint_G = currFrame_G.convertLocalToGlobal(landmarkPoint_L)
				#			candLandmarks_G.append(landmarkPoint_G)

				candLandmarks_G = nodeToGlobalLandmarks(controlPoses, self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])

				thisSplicedPaths.append((None, 1.0, rootSplice, None, candLandmarks_G))

			print "particle:", particleIndex, ",",  len(thisSplicedPaths), "localize jobs"
			poseFrame = Pose(hypPose0)
			globalMedial0 = []
			medialSpline0 = SplineFit(medial0, smooth=0.1)
			points0 = medialSpline0.getUniformSamples()
			#for p in points0:
			for p in medial0_vec:
				globalMedial0.append(poseFrame.convertLocalOffsetToGlobal(p))

			globalMedialP0 = poseFrame.convertLocalOffsetToGlobal(oldMedialP0)	
			globalMedialP1 = poseFrame.convertLocalOffsetToGlobal(oldMedialP1)	

			for spliceIndex in range(len(thisSplicedPaths)):
				
				branchSampleIndex = thisSplicedPaths[spliceIndex][0]
				probVal = thisSplicedPaths[spliceIndex][1]
				path = thisSplicedPaths[spliceIndex][2]
				landmarks_G = thisSplicedPaths[spliceIndex][4]

				orientedPath0 = orientPathLean(path, globalMedial0)
				
				localizeJobs.append([oldMedialP0, oldMedialU0, 0.0, oldMedialP1, oldMedialU1, 0.0, branchSampleIndex, spliceCount, orientedPath0, medial0_vec, medial1_vec, deepcopy(hypPose0), deepcopy(hypPose1), prevMedial0_vec, prevMedial1_vec, prevHypPose0, prevHypPose1, [], nodeID0, nodeID1, particleIndex, updateCount, self.hypothesisID, probVal, landmarks_G, landmark0_N, landmark1_N])

				self.pathPlotCount2 += 1
				spliceCount += 1

			allSplicedPaths += thisSplicedPaths


		print len(localizeJobs), "total localize jobs"

		time1 = time.time()

		results = batchLocalizeParticle(localizeJobs)

		time2 = time.time()

		print "TIME batchLocalizeParticle", self.hypothesisID, "=", time2-time1 
		print len(results[0]), "result arguments"

		maxLandmarkSum = 0.0
		for index in range(len(results)):
			part = results[index]
			poseLandmarkSum = part[51]
			if poseLandmarkSum > maxLandmarkSum:
				maxLandmarkSum = poseLandmarkSum

		print "maxLandmarkSum =", maxLandmarkSum

		isAllReject = True
		rejectDivergence = True

		""" if all rejected, try again without divergence filter """
		while isAllReject:

			""" set the rejection criteria """
			for index in range(len(results)):

				part = results[index]

				particleID = part[0]

				utilVal = part[45]
				spliceIndex = part[47]
				branchIndex = allSplicedPaths[spliceIndex][0]
				normMatchCount = allSplicedPaths[spliceIndex][1]
				spliceCurve = allSplicedPaths[spliceIndex][2]
				#pathID = allSplicedPaths[spliceIndex][3]

				initPose0 = part[48]
				initPose1 = part[49]

				overlapSum = part[50]
				poseLandmarkSum = part[51]


				newPose0 = part[2]

				isInterior1_0 = part[9]
				isExist1_0 = part[10]
				isInterior2_0 = part[15]
				isExist2_0 = part[16]
				contigFrac_0 = part[19]
				overlapSum_0 = part[20]

				newPose1 = part[24]
				isInterior1_1 = part[31]
				isExist1_1 = part[32]
				isInterior2_1 = part[37]
				isExist2_1 = part[38]
				contigFrac_1 = part[41]
				overlapSum_1 = part[42]


				isReject = False
				""" divergence is prohibited """
				if rejectDivergence and (isInterior1_0 or isInterior2_0 or isInterior1_1 or isInterior2_1):
					print "reject because divergence"
					isReject = True

				""" only a single extension is permitted, but not two """
				if isExist1_0 and isExist2_0 or isExist1_1 and isExist2_1:
					print "reject because double extension"
					isReject = True
				
				""" horrifically low contiguity is rejected out of hand """
				if contigFrac_0 <= 0.5 or contigFrac_1 <= 0.5:
					print "reject because low contiguity"
					isReject = True

				""" contiguity between current and previous pose """
				if overlapSum > 1e10:
					print "reject because no overlap"
					isReject = True

				#ANG_THRESH = 2.0*pi/3.0
				ANG_THRESH = 1.0*pi/3.0

				#if fabs(diffAngle(initPose0[2],newPose0[2])) > 2.0*pi/3.0 or fabs(diffAngle(initPose1[2],newPose1[2])) > 2.0*pi/3.0:

				if fabs(diffAngle(initPose0[2],newPose0[2])) > ANG_THRESH or fabs(diffAngle(initPose1[2],newPose1[2])) > ANG_THRESH:
					print "reject because change in pose angle is too different"
					isReject = True
				else:
					pass
					#print "localize angle diff:", particleID, initPose0[2], newPose0[2], initPose1[2], newPose1[2], fabs(diffAngle(initPose0[2],newPose0[2])), fabs(diffAngle(initPose1[2],newPose1[2]))

				normLandmarkSum = 0.0
				normAngDiff = pi-fabs(diffAngle(initPose0[2],newPose0[2])) 

				if isReject:
					newUtilVal = -1e100
					listCopy = list(part)
					listCopy[45] = newUtilVal
					tupleCopy = tuple(listCopy)

					results[index] = tupleCopy
				else:

					isAllReject = False

					""" check to avoid divide by zero """
					if maxLandmarkSum > 0.0:
						if maxLandmarkSum > poseLandmarkSum:
							normLandmarkSum = (maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum
							#newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0  * normLandmarkSum 
							#newProb = contigFrac_0 * normLandmarkSum 
							newProb = normLandmarkSum 
						else:
							""" maximum landmark cost, means utility is zeroed """
							newProb = 0.0
					else:
						#newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0  * 0.0
						newProb = contigFrac_0

					#newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0 * contigFrac_0 / overlapSum_0

					""" case if splice is root path """
					if normMatchCount != None:
						newProb *= normMatchCount

					listCopy = list(part)
					listCopy[45] = newProb
					tupleCopy = tuple(listCopy)

					results[index] = tupleCopy
				
				print "%d %d %d %d %d normMatchCount, newUtilVal, overlapSum:" % (nodeID0, self.hypothesisID, particleID, index, spliceIndex), normMatchCount, results[index][45], overlapSum, contigFrac_0, contigFrac_1, initPose0[2], newPose0[2], int(isExist1_0), int(isExist2_0), int(isExist1_1), int(isExist2_1), int(isInterior1_0), int(isInterior2_0), int(isInterior1_1), int(isInterior2_1), isReject, branchIndex, poseLandmarkSum, len(landmarks_G), normLandmarkSum, normAngDiff 

			""" break condition, do not accept divergence if we have already branched """
			if not rejectDivergence or self.isNodeBranching[nodeID0] or self.isNodeBranching[nodeID1]:
				#isAllReject = False
				break

			""" turn off divergence filter """
			rejectDivergence = False





		""" sort by pose particle index followed by utility value """
		sortedResults = sorted(results, key=itemgetter(0,45), reverse=True)

		print "particle sort"
		thisParticleID = -1
		distMin = 1e100
		filteredParticles = []
		probResults = {}
		for res in sortedResults:
			try:
				probResults[res[0]]
			except:
				probResults[res[0]] = []

			utilVal = res[45]
			probResults[res[0]].append(utilVal)
			if res[0] != thisParticleID:
				thisParticleID = res[0]
				distMin = res[1]
				filteredParticles.append(res)

		filteredParticles.reverse()

		print "localize results:"
		for key, values in probResults.iteritems():
			print key, values
			
		print "len(filteredParticles) =", len(filteredParticles)

		newParticleDist2 = []

		print "particle evaluation:", nodeID0, self.hypothesisID, updateCount
		for particleIndex in range(len(filteredParticles)):

			#if duplicateMatches[particleIndex] == particleIndex:

			part = filteredParticles[particleIndex]

			particleID = part[0]

			utilVal = part[45]
			spliceIndex = part[47]
			branchTupleIndex = allSplicedPaths[spliceIndex][0]
			normMatchCount = allSplicedPaths[spliceIndex][1]
			spliceCurve = allSplicedPaths[spliceIndex][2]
			#pathID = allSplicedPaths[spliceIndex][3]

			initPose0 = part[48]
			initPose1 = part[49]
			overlapSum = part[50]


			newPose0 = part[2]
			newDist0 = 0.5

			contigFrac_0 = part[19]
			overlapSum_0 = part[20]

			newPose1 = part[24]
				
			""" probability is the contigFrac squared """ 
			newProb = 0.0
			if utilVal > 0.0:
				newProb = utilVal
				#newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0 * contigFrac_0 / overlapSum_0
				#newProb *= float(normMatchCount)/float(self.maxMatchCount)
				#utilVal0 = (1.0-contigFrac_0) + (isExist1_0 or isExist2_0) + (1.0-contigFrac_1) + (isExist1_1 or isExist2_1)
			else:
				newProb = 0.0
			
			print "%d %d %d %d %d normMatchCount, utilVal, newProb" % (nodeID0, self.hypothesisID, particleID, particleIndex, spliceIndex), normMatchCount, utilVal, newProb, contigFrac_0, overlapSum_0

			print "localize angle diff:", particleID, initPose0[2], newPose0[2], initPose1[2], newPose1[2], fabs(diffAngle(initPose0[2],newPose0[2])), fabs(diffAngle(initPose1[2],newPose1[2]))

			particleObj = particleDist2[particleIndex].copy()

			#particleObj.updatePose(pose0, pose1)
			particleObj.pose0 = newPose0
			particleObj.pose1 = newPose1
			particleObj.hypDist = newDist0
			particleObj.weightVal = newProb
			particleObj.spliceCurve = spliceCurve



			""" Change branch if we are within a junction, otherwise, keep the most recent branch index """
			if branchTupleIndex != None:
				#particleObj.junctionData[pathID]["maxLikelihoodBranch"] = branchIndex
				particleObj.maxLikelihoodBranch = branchTupleIndex

				""" set branch points and control poses """
				branchResult = self.jointBranchEvaluations[branchTupleIndex] 
				controlSet = branchResult["controlSet"]
				branchPoses_L = branchResult["branchPoses_L"]

				controlPoses_G = computeGlobalControlPoses(controlSet, self.getParentHash())

				for pathID in branchPathIDs:
					branchPose_L = branchPoses_L[pathID]
					controlPose_G = controlPoses_G[pathID]
					currFrame = Pose(controlPose_G)
					branchPose_G = currFrame.convertLocalOffsetToGlobal(branchPose_L)
					particleObj.junctionData[pathID]["controlPose"] = controlSet[pathID]
					particleObj.junctionData[pathID]["localJunctionPose"] = branchPose_L
					particleObj.junctionData[pathID]["globalJunctionPose"] = branchPose_G

			#else:
			#	particleObj.maxLikelihoodBranch = None 



			newParticleDist2.append(particleObj)


		updateCount += 1
		self.poseParticles["updateCount"] = updateCount
		self.poseParticles["snapshots2"][0] = newParticleDist2
		numParticles = self.poseParticles["numParticles"]

		self.drawPoseParticles()

		""" now resample the particles """
		particleDist2 = self.poseParticles["snapshots2"][0]
		probSum = 0.0
		for part in particleDist2:
			probVal = part.weightVal
			probSum += probVal

		probParticles = []
		for k in range(len(particleDist2)):
			part = particleDist2[k]

			poseProbVal = part.weightVal

			""" if all 0.0 probabilities, make uniform distribution, set no localization flag """
			if probSum > 0.0:
				probParticles.append(poseProbVal/probSum)
			else:
				probParticles.append(float(1)/float(numParticles))
				self.isNoLocalize = True

		print "probParticles:", self.isNoLocalize, probParticles 

		""" now resample """
		resampledParticles2 = []
		numParticles = self.poseParticles["numParticles"]

		""" find the maximum likelihood particle.  Take average if tied for max """
		maxVal = 0.0
		maxIndex = 0
		numMax = 0
		for k in range(len(probParticles)):
			if probParticles[k] > maxVal:
				maxVal = probParticles[k]
				maxIndex = k 
				numMax = 1
			elif probParticles[k] == maxVal:
				numMax += 1

		print "maxIndex, numMax, maxVal:", maxIndex, numMax, maxVal

		""" if more than one max, than we don't change the node pose, we accept the localize on the canonical pose """
		if numMax >= 1:

			print "setting max likelihood poses", maxIndex, nodeID0, nodeID1, self.getNodePose(nodeID0), self.getNodePose(nodeID1)

			#self.nodePoses[nodeID0] = deepcopy(particleDist2[maxIndex].pose0)
			#self.nodePoses[nodeID1] = deepcopy(particleDist2[maxIndex].pose1)

			self.currMaxIndex = maxIndex

			""" change to the maximum likelihood branch position as well """
			self.updateMaxParticle(maxIndex)

			""" we update the canonical poses afterward because for now, we don't let node poses vary between particles """
			""" if we let them vary, then we would need to generate unique shoot maps for each particle """

			""" update the poses across all particles based on this max displacement result """
			self.setNodePose(nodeID0, deepcopy(particleDist2[maxIndex].pose0))
			self.setNodePose(nodeID1, deepcopy(particleDist2[maxIndex].pose1))

			""" using the results of localizePair() to position the nodes for generating map poses """
			self.generatePaths()


		else:
			print "no maximum particle!", probParticles
			raise

		""" FIXME:  do we really do nothing here when we tie for maximum? """


		print "particle evaluation:", nodeID0, self.hypothesisID, updateCount
		for i in range(numParticles):
			sampVal = random.random()

			k = 0
			probSum = probParticles[k]

			while sampVal > probSum:
				k += 1
				probSum += probParticles[k]

			oldPart2 = particleDist2[k]
			
			myKeys = oldPart2.junctionData.keys()
			printStr = ""
			for key in myKeys:
				printStr += repr(oldPart2.junctionData[key]["globalJunctionPose"]) 
				printStr += " " + repr(self.pathClasses[key]["globalJunctionPose"]) 

			#print "resampling particle", k, probParticles[k], printStr
			print "resampling particle", k, probParticles[k]
			print oldPart2.prevPose0, oldPart2.prevPose1
			print oldPart2.dispPose0, oldPart2.dispPose1
			print oldPart2.pose0, oldPart2.pose1

			resampledParticles2.append(oldPart2.copy())

		updateCount += 1
		self.poseParticles["updateCount"] = updateCount
		self.poseParticles["snapshots2"][0] = resampledParticles2

		self.drawPoseParticles()

		return


	@logFunction
	def drawPoseParticles(self):

		updateCount = self.poseParticles["updateCount"] 
		#particleDist = self.poseParticles["snapshots"][updateCount]
		particleDist2 = self.poseParticles["snapshots2"][0]


		localPathSets = {}

		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:

			if pathID != 0:

				junctionDetails = self.localLeaf2LeafPathJunctions[pathID]

				origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
				origJuncPose[2] = 0.0
				origJuncOrigin = Pose(origJuncPose)

				localPathSegs = junctionDetails["localSegments"]

				localPathSets[pathID] = localPathSegs


		#print "particles:", particleDist2

		pylab.ioff()
		print pylab.isinteractive()
		pylab.clf()
		#pylab.axis("equal")

		for k,path in  self.trimmedPaths.iteritems():
			print "path has", len(path), "points"
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color = self.colors[k], linewidth=4)


		nodeID0 = self.poseData.numNodes-2
		nodeID1 = self.poseData.numNodes-1
		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]

		# newPart = (newPose0, newPose1, 0, newDist0, ('t0', 't1'), newProb)
		hypPointsX_0 = []
		hypPointsY_0 = []
		hypPointsX_1 = []
		hypPointsY_1 = []
		for part in particleDist2:
			hypPose0 = part.pose0
			hypPose1 = part.pose1
			hypPointsX_0.append(hypPose0[0])
			hypPointsY_0.append(hypPose0[1])
			hypPointsX_1.append(hypPose1[0])
			hypPointsY_1.append(hypPose1[1])

			thisSplice = part.spliceCurve




			poseOrigin0 = Pose(hypPose0)

			xP = []
			yP = []
			for p in medial0:
				p1 = poseOrigin0.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'r', alpha = 0.2, zorder=9)


			allPathIDs = self.getPathIDs()
			branchPathIDs = copy(allPathIDs)
			branchPathIDs.remove(0)
			branchPathIDs.sort()

			controlPoses_L = {}
			for pathID in allPathIDs:
				controlPoses_L[pathID] = deepcopy(part.junctionData[pathID]["controlPose"])

			print "checkA:", controlPoses_L, self.getParentHash()
			controlPoses_G = computeGlobalControlPoses(controlPoses_L, self.getParentHash())

			junctionPoses = {}
			for pathID in allPathIDs:
				if pathID != 0:
					junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

			xP = []
			yP = []
			for pathID in allPathIDs:
				currFrame = Pose(controlPoses_G[pathID])
				for point_L, pointThresh, pointName in self.localLandmarks[pathID]:
					point_G = currFrame.convertLocalToGlobal(point_L)

					xP.append(point_G[0])
					yP.append(point_G[1])

			if len(xP) > 0:
				pylab.scatter(xP, yP, color='k', linewidth=1, zorder=9, alpha=0.4)

			xP = []
			yP = []
			for pathID, branchPoint in junctionPoses.iteritems():
				xP.append(branchPoint[0])
				yP.append(branchPoint[1])

			if len(xP) > 0:
				pylab.scatter(xP, yP, color='r', linewidth=1, zorder=10, alpha=0.4)

			xP = []
			yP = []
			for pathID in allPathIDs:


				if pathID != 0:


					#if pathID in part.junctionData.keys():
					if False:

						probDist = part.junctionData[pathID]["probDist"]
						branchPoseDist = part.junctionData[pathID]["branchPoseDist"]
						controlPoseDist = part.junctionData[pathID]["controlPoseDist"]

						""" FIXME:  convert to control matrix and local controls """

						maxProb = -1e100
						maxIndex = 0
						for k in range(len(probDist)):
							modPose0 = deepcopy(controlPoseDist[k])

							xP.append(modPose0[0])
							yP.append(modPose0[1])

							if probDist[k] > maxProb:
								maxProb = probDist[k]
								maxIndex = k

						modPose0 = deepcopy(controlPoseDist[maxIndex])
						modPose0 = deepcopy(part.junctionData[pathID]["controlPose"])

						print "probDist:", probDist
						print "branchPoseDist:", branchPoseDist
						print "controlPoseDist:", controlPoseDist

					else:
						#modPose0 = deepcopy(part.junctionData[pathID]["globalJunctionPose"])
						modPose0 = deepcopy(controlPoses_G[pathID])
						xP.append(modPose0[0])
						yP.append(modPose0[1])


					#print "modPose0:", modPose0
					#modPose0[2] = 0.0
					#modOrigin0 = Pose(modPose0)

			#if len(xP) > 0:
			#	pylab.scatter(xP, yP, color='k', linewidth=1, zorder=10, alpha=0.4)




			poseOrigin1 = Pose(hypPose1)

			xP = []
			yP = []
			for p in medial1:
				p1 = poseOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'b', alpha = 0.2, zorder=9)

		#pylab.scatter(hypPointsX_0, hypPointsY_0, color='r', linewidth=1, zorder=10, alpha=0.2)
		#pylab.scatter(hypPointsX_1, hypPointsY_1, color='b', linewidth=1, zorder=10, alpha=0.2)

		pose0 = self.getNodePose(nodeID0)
		pose1 = self.getNodePose(nodeID1)
		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]

		poseOrigin0 = Pose(pose0)
		poseOrigin1 = Pose(pose1)

		xP = []
		yP = []
		for p in medial0:
			p1 = poseOrigin0.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		#pylab.plot(xP,yP,color='k', zorder=9, alpha=0.2)
		pylab.plot(xP,yP,color='k', zorder=10, linewidth=1)
		#pylab.plot(xP,yP,color='k', zorder=8, linewidth=1)

		xP = []
		yP = []
		for p in medial1:
			p1 = poseOrigin1.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		#pylab.plot(xP,yP,color='k', zorder=9, alpha=0.2)
		pylab.plot(xP,yP,color='k', zorder=10, linewidth=1)
		#pylab.plot(xP,yP,color='k', zorder=8, linewidth=1)

		
		#pylab.xlim(-6,16.48)
		
		#pylab.ylim(-16,16)
		pylab.axis("equal")
		#pylab.xlim(-16,16)
		#pylab.title("%d particles" % len(particleDist))
		pylab.title("nodeIDs %d %d hypID %d" % (nodeID0, nodeID1, self.hypothesisID))


		#pylab.savefig("moveEstimate2_%04u_%04u.png" % (self.posePlotCount, self.hypothesisID))
		pylab.savefig("moveEstimate2_%04u_%04u.png" % (self.hypothesisID, self.posePlotCount))
		print "moveEstimate2_%04u_%04u.png" % (self.hypothesisID, self.posePlotCount)

		self.posePlotCount += 1

	@logFunction
	def updatePoseData(self, poseData):
		self.poseData = deepcopy(poseData)
		print "updated poseData:", poseData.numNodes, self.poseData.numNodes

	@logFunction
	def copy(self, hypothesisID):

		print "creating hypothesis", hypothesisID

		#newObj = deepcopy(self)
		#newObj.hypothesisID = hypothesisID

		#return newObj

		newObj = MapState(self.poseData, hypothesisID)

		newObj.numPoseParticles = deepcopy(self.numPoseParticles)

		newObj.poseParticles = deepcopy(self.poseParticles)

		updateCount = newObj.poseParticles["updateCount"] 
		particleDist2 = deepcopy(newObj.poseParticles["snapshots2"][0])
		for part in particleDist2:
			part.mapStateID = hypothesisID

		newObj.parentParticle = self.hypothesisID 

		newObj.pathClasses = deepcopy(self.pathClasses)
		newObj.pathTermsVisited = deepcopy(self.pathTermsVisited)
		newObj.topDict = deepcopy(self.topDict)
		newObj.pathIDs = deepcopy(self.pathIDs)
		newObj.paths = deepcopy(self.paths)
		newObj.hulls = deepcopy(self.hulls)
		newObj.trimmedPaths = deepcopy(self.trimmedPaths)
		newObj.nodeRawPoses = deepcopy(self.nodeRawPoses)
		newObj.nodePoses = deepcopy(self.nodePoses)
		newObj.gndPoses = deepcopy(self.gndPoses)
		newObj.gndRawPoses = deepcopy(self.gndRawPoses)
		newObj.origPoses = deepcopy(self.origPoses)
		newObj.utility = deepcopy(self.utility)
		newObj.isNodeBranching = deepcopy(self.isNodeBranching)

		newObj.alphaPlotCount = self.alphaPlotCount
		newObj.medialCount = self.medialCount
		newObj.topCount = self.topCount
		newObj.spliceCount = self.spliceCount
		newObj.multiDepCount = self.multiDepCount
		newObj.overlapPlotCount = self.overlapPlotCount
			
		newObj.pathPlotCount = self.pathPlotCount
		newObj.pathPlotCount2 = self.pathPlotCount2
		newObj.overlapPlotCount2 = self.overlapPlotCount2
		newObj.medialSouptCount = self.medialSoupCount

		
		newObj.rootPoint = deepcopy(self.rootPoint)
		
		newObj.pathGraph = deepcopy(self.pathGraph)
		newObj.joins = deepcopy(self.joins)
		newObj.junctions = deepcopy(self.junctions)
		newObj.terminals = deepcopy(self.terminals)
		newObj.topDict = deepcopy(self.topDict)
		newObj.allSplices = deepcopy(self.allSplices)

		newObj.orderedPathIDs1 = deepcopy(self.orderedPathIDs1)
		newObj.orderedPathIDs2 = deepcopy(self.orderedPathIDs2)

		" departure events for node 1 "
		newObj.departures1 = deepcopy(self.departures1) 
		newObj.interiors1 = deepcopy(self.interiors1)  
		newObj.depPoints1 = deepcopy(self.depPoints1) 
		newObj.distances1 = deepcopy(self.distances1)  
		newObj.depAngles1 = deepcopy(self.depAngles1)
		newObj.contig1 = deepcopy(self.contig1)

		" departure events for node 2 "
		newObj.departures2 = deepcopy(self.departures2) 
		newObj.interiors2 = deepcopy(self.interiors2) 
		newObj.depPoints2 = deepcopy(self.depPoints2)  
		newObj.distances2 = deepcopy(self.distances2) 
		newObj.depAngles2 = deepcopy(self.depAngles2)  
		newObj.contig2 = deepcopy(self.contig2)  

		newObj.departureResultSet1 = deepcopy(self.departureResultSet1)  
		newObj.departureResultSet2 = deepcopy(self.departureResultSet2)

		newObj.localPaths = deepcopy(self.localPaths)
		newObj.localHulls = deepcopy(self.localHulls)
		newObj.localLandmarks = deepcopy(self.localLandmarks)
		newObj.localLeaf2LeafPathJunctions = deepcopy(self.localLeaf2LeafPathJunctions) 

		newObj.branchDiverges = deepcopy(self.branchDiverges)

		return newObj

	@logFunction
	def computeEval(self):

		keys = self.trimmedPaths.keys()
		keys.sort()


		parentHash = self.getParentHash()

		totalSum = 0.0

		for j in keys:
			parentID = parentHash[j]
			if parentID != None:
				path1 = self.trimmedPaths[j]
				path2 = self.trimmedPaths[parentID]

				branchPose_G = self.branchPoses_G[j]

				#pathSpline1 = SplineFit(path1)
				#pathPoints1 = pathSpline1.getUniformSamples()
				pathSpline2 = SplineFit(path2)
				pathPoints2 = pathSpline2.getUniformSamples()

				minPoint, minI, minDist = gen_icp.findClosestPointInA(pathPoints2, branchPose_G)

				ang1 = branchPose_G[2]
				ang2 = pathPoints2[minI][2]

				ang1 = normalizeAngle(ang1)
				ang2 = normalizeAngle(ang2)
				branchDiffAngle1 = fabs(normalizeAngle(diffAngle(ang1, ang2)))
				branchDiffAngle2 = fabs(normalizeAngle(diffAngle(ang1, normalizeAngle(ang2 + math.pi))))

				print "computeEval:", self.hypothesisID, parentID, j, ang1, ang2, branchDiffAngle1, branchDiffAngle2, self.branchDiverges
				#if branchDiffAngle1 < 0.1 or branchDiffAngle2 < 0.1 or False in self.branchDiverges.values():
				if False in self.branchDiverges.values():
						totalSum = 1e100

		"""
		for j in keys:
			path1 = self.trimmedPaths[j]

			for k in keys:
				path2 = self.trimmedPaths[k]

				if j < k:


					resultSum = self.getPathOverlapSum(path1, path2, j, k, plotIter = True)
					print "computing sum of", j, "and", k, "=", resultSum
					totalSum += resultSum
		"""


		self.mapOverlapSum = totalSum
		#self.mapOverlapSum = 0.0

		return totalSum
	

	""" state save and restore methods """
	@logFunction
	def saveState(self, saveCount):
		
		saveFile = ""
		
		saveFile += "self.pathClasses = " + repr(self.pathClasses) + "\n"
		saveFile += "self.pathTermsVisited = " + repr(self.pathTermsVisited) + "\n"
		saveFile += "self.pathIDs = " + repr(self.pathIDs) + "\n"		 
		saveFile += "self.paths = " + repr(self.paths) + "\n"		 
		saveFile += "self.hulls = " + repr(self.hulls) + "\n"		 
		saveFile += "self.trimmedPaths = " + repr(self.trimmedPaths) + "\n"		   

		saveFile += "self.nodeRawPoses = " + repr(self.nodeRawPoses) + "\n"		 
		saveFile += "self.nodePoses = " + repr(self.nodePoses) + "\n"		 
		saveFile += "self.gndPoses = " + repr(self.gndPoses) + "\n"		 
		saveFile += "self.gndRawPoses = " + repr(self.gndRawPoses) + "\n"		 
		saveFile += "self.hypothesisID = " + repr(self.hypothesisID) + "\n"

		f = open("mapStateSave_%04u_%04u.txt" % (self.hypothesisID, saveCount), 'w')
		f.write(saveFile)
		f.close()
		
	@logFunction
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/mapStateSave_%04u_%04u.txt" % (self.hypothesisID, numNodes+1)
		f = open(dirName + "/mapStateSave_%04u_%04u.txt" % (self.hypothesisID, numNodes+1), 'r')		 
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		exec(saveStr)
		

	""" member access methods """

	@logFunction
	def getControlPose(self, pathID):

		controlPose = self.pathClasses[pathID]["controlPose"]
		return controlPose

	@logFunction
	def setLocalJunctionPose(self, pathID, localJuncPose):

		if pathID != 0:

			self.pathClasses[pathID]["localJunctionPose"] = localJuncPose

			globalControlPoses = computeGlobalControlPoses(self.getControlPoses(), self.getParentHash())


			""" convert the controlPose to coordinates local to the parent frame """
			shootControlPose = globalControlPoses[pathID]
			currFrame = Pose(shootControlPose)
			globalJuncPose = currFrame.convertLocalOffsetToGlobal(localJuncPose)

			self.pathClasses[pathID]["globalJunctionPose"] = globalJuncPose

	@logFunction
	def getGlobalJunctionPose(self, pathID):

		try:
			globJuncPose = self.pathClasses[pathID]["globalJunctionPose"]
		except:
			
			branchNodeID = self.pathClasses[pathID]["branchNodeID"]
			
			if branchNodeID == None:
				return None
			
			localDivergencePose = self.pathClasses[pathID]["localDivergencePose"]
			
			poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
			globJuncPose = poseOrigin.convertLocalOffsetToGlobal(localDivergencePose)

		return globJuncPose
		
	@logFunction
	def getNodes(self, pathID):
		return self.pathClasses[pathID]["nodeSet"]
	
	@logFunction
	def getPathIDs(self):
		allPathIDs = self.pathClasses.keys()
		allPathIDs.sort()
		return deepcopy(allPathIDs)

	@logFunction
	def getParentHash(self):

		parentPathIDs = {}
		pathIDs = self.getPathIDs()

		for pathID in pathIDs:
			cPath = self.getPath(pathID)
			parentPathID = cPath["parentID"]
			parentPathIDs[pathID] = parentPathID

		return parentPathIDs

	@logFunction
	def getControlPoses(self):
		controlPoses = {}
		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]

		return controlPoses
	

	@logFunction
	def getFinalControlPose(self, pathID):
		
		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)


		return globalControlPoses[pathID]

	@logFunction
	def getParentPathID(self, pathID):
		parentID = self.pathClasses[pathID]["parentID"]
		return parentID
	
	@logFunction
	def getPath(self, pathID):
		return self.pathClasses[pathID]
	
	def getAllSplices(self, plotIter = False):
		if not self.isChanged:
			print "returning without computing, not isChanged"
			return self.allSplices, self.terminals, self.junctions
		else:
			print "recomputing allSplices"
			self.generatePaths()
			return self.allSplices, self.terminals, self.junctions

	#@logFunction
	def computeAllSplices(self, plotIter = False):

		print "paths:", self.paths
		print "junctions:", self.junctions
		print "terminals:", self.terminals

		#if not self.isChanged:
		#	print "returning without computing, not isChanged"
		#	return self.allSplices, self.terminals, self.junctions
		#else:
		#	print "computing allSplices"



		" determine which paths are leaves "
		
		pathIDs = self.getPathIDs()
		isAParent = {}
		parents = {}
		pathIDGraph = graph.graph()

		for pathID in pathIDs:
			isAParent[pathID] = False
			pathIDGraph.add_node(pathID, [])
		
		for pathID in pathIDs:	  
			parentPathID = self.getPath(pathID)["parentID"]
			parents[pathID] = parentPathID
			if parentPathID != None:
				isAParent[parentPathID] = True
				pathIDGraph.add_edge(pathID, parentPathID)
		
		leafCount = 0
		for pathID in pathIDs:
			if not isAParent[pathID]:
				leafCount += 1


		" build topological graph "
		topGraph = graph.graph()
		
		topNodes = []
					   
		for k, term in self.terminals.iteritems():
			topGraph.add_node(term[0], term[1])
			print "add_node", term[0], term[1]
			topNodes.append(term[0])
		for k, junc in self.junctions.iteritems():
			topGraph.add_node(junc[2], junc[3])
			print "add_node", junc[2], junc[3]
			topNodes.append(junc[2])


		for pathID in pathIDs:
			pathIndex = []
			if pathID == 0:
				term = self.terminals[0]
				graphNodeKey = term[0]
				pathIndex.append(graphNodeKey[1])
				
				term = self.terminals[1]
				graphNodeKey = term[0]
				pathIndex.append(graphNodeKey[1])
			else:
				term = self.terminals[pathID+1]
				graphNodeKey = term[0]
				pathIndex.append(graphNodeKey[1])
				
			for k, junc in self.junctions.iteritems():
				graphNodeKey = junc[2]
				if graphNodeKey[0] == pathID:
					pathIndex.append(graphNodeKey[1])
						
			pathIndex.sort()
			
			print "pathIndex:", pathIndex
			
			for k in range(len(pathIndex)-1):
				print "add_edge", (pathID,pathIndex[k]),(pathID,pathIndex[k+1])
				topGraph.add_edge((pathID,pathIndex[k]),(pathID,pathIndex[k+1]))

			if pathID != 0:
				graphNodeKey = self.junctions[pathID][2]
				print "add_edge", graphNodeKey,(pathID,pathIndex[0])
				topGraph.add_edge(graphNodeKey,(pathID,pathIndex[0]))


		results = {}
		for i in range(len(pathIDs)):
			pathID1 = pathIDs[i]
			for j in range(i,len(pathIDs)):
				
				pathID2 = pathIDs[j]
				
				if pathID1 == pathID2:
					"term1 to term 2"

					termPaths = []
					if pathID1 == 0:
						termIDs = [self.topDict["t%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
					else:
						termIDs = [self.topDict["j%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
					termPaths.append(termIDs)  

					orderedPathIDs = [pathID1]
				else:
					
					startPathID = pathID1
					endPathID = pathID2
		
					shortestPathSpanTree, shortestDist = pathIDGraph.shortest_path(endPathID)
					nextPathID = shortestPathSpanTree[startPathID]
			
					orderedPathIDs = [startPathID, nextPathID]
					while nextPathID != endPathID:
						nextPathID = shortestPathSpanTree[nextPathID]
						orderedPathIDs.append(nextPathID)			 

					
					" start at first path ID "
					termPaths = []
					if parents[orderedPathIDs[0]] == orderedPathIDs[1]:
						" going to parent "
						termIDs = [self.topDict["t%u" % (orderedPathIDs[0]+1)]]
						termPaths.append(termIDs)
					
					elif parents[orderedPathIDs[1]] == orderedPathIDs[0]:
						" going to child "
						termIDs = [self.topDict["t%u" % (orderedPathIDs[0]+1)]]
						termPaths.append(termIDs)

						if orderedPathIDs[0] == 0:
							termIDs = [self.topDict["t%u" % 0]]
							termPaths.append(termIDs)
						else:
							termIDs = [self.topDict["j%u" % orderedPathIDs[0]]]
							termPaths.append(termIDs)
					else:
						print parents
						print orderedPathIDs
						print termPaths
						raise
						
					" transitions between middle path IDs "
					for k in range(len(orderedPathIDs)-1):
						if parents[orderedPathIDs[k]] == orderedPathIDs[k+1]:
							" child -> parent:	prev -> junction "
							for ids in termPaths:
								ids.append(self.topDict["j%u" % orderedPathIDs[k]])

						elif parents[orderedPathIDs[k+1]] == orderedPathIDs[k]:
							" parent -> child: prev -> junction "
							for ids in termPaths:
								ids.append(self.topDict["j%u" % orderedPathIDs[k+1]])

						else:
							print parents
							print orderedPathIDs
							print termPaths
							print k
							raise
					
					" last path ID "
					if parents[orderedPathIDs[-2]] == orderedPathIDs[-1]:
						" child -> parent:	prev -> junction "
						oldTermPaths = termPaths
						termPaths = []
						" split into 2 possibilities into parent path ID "
						for ids in oldTermPaths:
							ids1 = deepcopy(ids)
							ids2 = deepcopy(ids)
							if orderedPathIDs[-1] == 0:
								ids1.append(self.topDict["t%u" % orderedPathIDs[-1]])
								ids2.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])								  
							else:
								ids1.append(self.topDict["j%u" % orderedPathIDs[-1]])
								ids2.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])
							
							termPaths.append(ids1)
							termPaths.append(ids2)

					elif parents[orderedPathIDs[-1]] == orderedPathIDs[-2]:
						" parent -> child: prev -> junction "
						for ids in termPaths:
							ids.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])
										
					else:
						print parents
						print orderedPathIDs
						print termPaths
						raise
				
				#globJuncPose = self.pathClasses[pathID2]["globalJunctionPose"]

				#""" the terminals of the child shoot """
				#dist1 = sqrt((globJuncPose[0]-newPath3[0][0])**2 + (globJuncPose[1]-newPath3[0][1])**2)
				#dist2 = sqrt((globJuncPose[0]-newPath3[-1][0])**2 + (globJuncPose[1]-newPath3[-1][1])**2)

				#if dist1 < dist2:
				#	childTerm = newPath3[-1]
				#else:
				#	childTerm = newPath3[0]

				#dist3 = sqrt((childTerm[0]-particlePath[0][0])**2 + (childTerm[1]-particlePath[0][1])**2)
				#dist4 = sqrt((childTerm[0]-particlePath[-1][0])**2 + (childTerm[1]-particlePath[-1][1])**2)

				#if dist3 < dist4:
				#	childTerm2 = particlePath[-1]
				#else:
				#	childTerm2 = particlePath[0]

				#""" two different splices between child and parent terminals """
				#termPaths = [(parentTerm1, childTerm), (parentTerm2, childTerm), (childTerm2, childTerm)]


				finalResults = []


				print "termPaths:", termPaths
				
				for termPath in termPaths:

					joinPairs = []

					startKey = termPath[0]
					endKey = termPath[-1]

					print "startKey, endKey:", startKey, endKey

					""" use splice skeleton to find shortest path """
					startPose = self.pathGraph.get_node_attributes(startKey)
					endPose = self.pathGraph.get_node_attributes(endKey)

					print "startPose, endPose:", startPose, endPose

					minStartDist = 1e100
					minStartNode = None
					minEndDist = 1e100
					minEndNode = None
					
					for edge in self.spliceSkeleton.edges():
					
						globalNodePoint1 = edge[0]
						globalNodePoint2 = edge[1]

						dist1 = sqrt((globalNodePoint1[0]-startPose[0])**2 + (globalNodePoint1[1]-startPose[1])**2)
						dist2 = sqrt((globalNodePoint2[0]-startPose[0])**2 + (globalNodePoint2[1]-startPose[1])**2)

						if dist1 < minStartDist:
							minStartDist = dist1
							minStartNode = globalNodePoint1

						if dist2 < minStartDist:
							minStartDist = dist2
							minStartNode = globalNodePoint2

						dist1 = sqrt((globalNodePoint1[0]-endPose[0])**2 + (globalNodePoint1[1]-endPose[1])**2)
						dist2 = sqrt((globalNodePoint2[0]-endPose[0])**2 + (globalNodePoint2[1]-endPose[1])**2)

						if dist1 < minEndDist:
							minEndDist = dist1
							minEndNode = globalNodePoint1

						if dist2 < minEndDist:
							minEndDist = dist2
							minEndNode = globalNodePoint2

			
					startNode = minStartNode
					endNode = minEndNode


					print "nodes path from", startNode, "to", endNode
					shortestSpliceTree, shortestSpliceDist = self.spliceSkeleton.shortest_path(endNode)
					currNode = shortestSpliceTree[startNode]					 
					splicedSkel = [startNode]
					while currNode != endNode:
						print "currNode:", currNode
						splicedSkel.append(currNode)
						nextNode = shortestSpliceTree[currNode]
						currNode = nextNode
					splicedSkel.append(currNode)

					splicedSkel = ensureEnoughPoints(splicedSkel, max_spacing = 0.08, minPoints = 5)
					spliceSpline1 = SplineFit(splicedSkel, smooth=0.1)
					splicePoints1 = spliceSpline1.getUniformSamples()


					""" find splice using the old method """
					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					splicedPath = []
					while currNode != endKey:
						splicedPath.append(self.pathGraph.get_node_attributes(currNode))
						nextNode = shortestPathSpanTree[currNode]
						
						for join in self.joins:
							if join[0] == currNode and join[1] == nextNode:
								joinPairs.append((len(splicedPath)-1,len(splicedPath)))
							elif join[1] == currNode and join[0] == nextNode:
								joinPairs.append((len(splicedPath)-1,len(splicedPath)))
							
						currNode = nextNode

					splicedPath.append(self.pathGraph.get_node_attributes(currNode))


					" if there are any joins, we must interpolate a smooth transition "
					lastIndex = 0
					newPath = []
					for pair in joinPairs:
						index1 = pair[0]
						index2 = pair[1]
						cPoints = [splicedPath[index1-1], splicedPath[index1], splicedPath[index2], splicedPath[index2+1]]				 
						spline = SplineFit(cPoints, smooth=0.0, kp=3)
						points = spline.getUniformSamples(spacing = 0.1)
						
						newPath += splicedPath[lastIndex:index1]
						newPath += points
						
						lastIndex = index2+2

					newPath += splicedPath[lastIndex:]


					sPath = {}
					sPath['orderedPathIDs'] = orderedPathIDs
					sPath['path'] = splicePoints1
					sPath['oldPath'] = newPath
					sPath['termPath'] = termPath
					sPath['skelPath'] = splicePoints1
					
					finalResults.append(sPath)
					
				results[(pathID1,pathID2)] = finalResults
		
		print len(results), "results:"
		for k, result in results.iteritems():
			print k, ":"
			for sPath in result:
				print sPath['orderedPathIDs'], sPath['termPath'], len(sPath['path'])


		if plotIter:
			print "plotting splicedPath"
			pylab.clf()

			
			#for k,path in  self.trimmedPaths.iteritems():
			#	print "path has", len(path), "points"
			#	xP = []
			#	yP = []
			#	for p in path:
			#		xP.append(p[0])
			#		yP.append(p[1])
			#	
			#	pylab.plot(xP,yP, color = self.colors[k])
	
	
			for k, result in results.iteritems():
				for sPath in result:
					#path = sPath['path']
					path = sPath['skelPath']
	
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])
					pylab.plot(xP,yP, color='b', zorder=2)


			for k in pathIDs:
				xP = []
				yP = []
				
				for p in self.paths[k]:
					xP.append(p[0])
					yP.append(p[1])
				pylab.plot(xP,yP, color=self.colors[k], linewidth=1, zorder=1)

			xP = []
			yP = []
			for k,path in  self.pathClasses.iteritems():
				if k != 0:
					globalJunctionPose = path["globalJunctionPose"]
					xP.append(globalJunctionPose[0])
					yP.append(globalJunctionPose[1])

			pylab.scatter(xP,yP, color='k', zorder=3)

			#self.drawWalls()
	
			print "saving splicedPath_%04u_%04u_%04u.png" % (self.hypothesisID, self.poseData.numNodes, self.spliceCount)
			pylab.title("Spliced Paths, pathIDs = %s" %  self.getPathIDs())
			pylab.savefig("splicedPath_%04u_%04u_%04u.png" % (self.hypothesisID, self.poseData.numNodes, self.spliceCount))
			self.spliceCount += 1

		self.allSplices = results

		#self.isChanged = False
		
		return self.allSplices, self.terminals, self.junctions


	@logFunction
	def isNodeExist(self, nodeID, pathID):
		return nodeID in self.pathClasses[pathID]["nodeSet"]

	
	@logFunction
	def localizeBranchPoint(self, orientedPathSoup, junctionPose1, parentPathID, pathID1, plotIter = False):

		parentPath = self.paths[parentPathID]

		juncOrigin1 = Pose(junctionPose1)
		
		" curves of both paths "

		globalPath2 = parentPath

		globalSpline2 = SplineFit(globalPath2, smooth=0.1)
		globalSamples2 = globalSpline2.getUniformSamples(spacing = 0.04)
		originU2 = globalSpline2.findU(junctionPose1)
		minDist, originU2, uPoint = globalSpline2.findClosestPoint(junctionPose1)
		
		u2 = originU2
		u1 = 0.5  # value is meaningless since data is not a spline curve
		#angGuess = 0.0
		#angGuess = -uPoint[2] + junctionPose1[2]
		angGuess = junctionPose1[2]
		
		resultPose1, lastCost1, matchCount1 = gen_icp.branchEstimateICP([u1,u2,angGuess], junctionPose1, orientedPathSoup, globalPath2, plotIter = plotIter, n1 = pathID1, n2 = parentPathID)

		print "matchCount,cost,result:", matchCount1, lastCost1, resultPose1
		
		return resultPose1, lastCost1, matchCount1
	 

	@logFunction
	def resetBranches(self):


		self.branchEvaluations = {}
		self.jointBranchEvaluations = {}

		parentPathIDs = self.getParentHash()

		#self.DIV_LEN = 0.2

		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:
			if self.pathClasses[pathID]["parentID"] != None:
				parentPathID = parentPathIDs[pathID]
				pathSpline = SplineFit(self.paths[parentPathID])

				totalDist = pathSpline.dist_u(1.0)

				branchSpace = {}

				currDist = 0.0
				branchSpace[currDist] = None
				while currDist <= totalDist:
					currDist += self.DIV_LEN
					branchSpace[currDist] = None

				self.branchEvaluations[pathID] = branchSpace


		self.branchControlPoses = {}
		self.branchArcDists = {}
		for pathID in allPathIDs:
			if pathID != 0:

				self.branchControlPoses[pathID] = {}
				self.branchArcDists[pathID] = []

				parentPathID = parentPathIDs[pathID]
				pathSpline = SplineFit(self.localPaths[parentPathID])

				totalDist = pathSpline.dist_u(1.0)

				""" arc distances """
				currDist = 0.0
				self.branchArcDists[pathID].append(currDist)

				""" control poses """
				controlPose_P = pathSpline.getPointOfDist(currDist)
				controlPose_P[2] = 0.0
				self.branchControlPoses[pathID][currDist] = controlPose_P


				while currDist <= totalDist:

					""" arc distances """
					currDist += self.DIV_LEN
					self.branchArcDists[pathID].append(currDist)

					""" control poses """
					controlPose_P = pathSpline.getPointOfDist(currDist)
					controlPose_P[2] = 0.0
					self.branchControlPoses[pathID][currDist] = controlPose_P
			

		self.jointBranchEvaluations = {}


	@logFunction
	def getPrecomputedBranch(self, pathID, arcDist):

		"""
		result values are as follows:
		0 = parent ID
		1 = pose of the branch point
		2 = arc distance on the parent curve of the branch point
		3 = child ID
		4 = match count from ICP algorithm
		5 = last cost from ICP algorithm
		6 = discrepancy distance of new pose from old pose, zeroed for this function
		7 = angle discrepancy difference of new pose from old pose
		8 = initial probability value of this branch position, initialized to 0.0 in this function
		9 = set of splices including the branch at this position
		10 = angle discrepancy 
		11 = dist discrepancy 
		"""

		if self.pathClasses[pathID]["parentID"] != None: 
			currBranchSpace = self.branchEvaluations[pathID]

			currKeys = currBranchSpace.keys()
			currKeys.sort()


			" find the bin for this arc distance "
			thisKey = currKeys[0]
			for k in range(len(currKeys)):
				val = currKeys[k]
				if val >= arcDist:
					break
				thisKey = val

			#print "precompute value:", arcDist, thisKey

			if currBranchSpace[thisKey] != None:
				#print "returned precomputation of branch", thisKey
				return deepcopy(currBranchSpace[thisKey])
			else:

				print "returned precomputation of branch", arcDist, thisKey, None

			raise

			return None


	@logFunction
	def batchPrecomputeBranches(self, jointComboControlSets):

		""" precompute the evaluation for branches and cache result """
		""" each branch point is the max likelihood of each pose particle """
		""" computes discretized samples around each maximum likelihood branch point """

		""" probSets:  [ (pathID, globJuncPose), ] """
		""" possible locations of branch point shoot pathID """

		jointComboControlSets = sorted(jointComboControlSets, reverse=True)

		allPathIDs = self.getPathIDs()
		branchPathIDs = deepcopy(allPathIDs)
		branchPathIDs.remove(0)
		branchPathIDs.sort()


		""" data needed for branch computation processes """
		localSkeletons = {}
		controlPoses = {}
		junctionPoses = {}
		localPathSegsByID = {}
		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			localSkeletons[pathID] = self.localLeaf2LeafPathJunctions[pathID]["skeletonGraph"]
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			junctionPoses[pathID] = self.pathClasses[pathID]["localJunctionPose"]
			localPathSegsByID[pathID] = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]

		parentPathIDs = self.getParentHash()

		""" construct the joint control pose problem sets """
		jointBranchJobs = []
		for controlSet in jointComboControlSets:

			thisControlPoses = deepcopy(controlPoses)
			thisArcDists = {}

			for k in range(len(controlSet)):

				arcDist = controlSet[k][0]
				cPose = controlSet[k][1]

				pathID = branchPathIDs[k]
				thisControlPoses[pathID] = cPose
				thisArcDists[pathID] = arcDist

			#self.localLandmarks
			jointBranchJobs.append((localPathSegsByID, self.localPaths, localSkeletons, thisControlPoses, junctionPoses, self.localLandmarks, parentPathIDs, thisArcDists, len(self.nodePoses)-1, self.hypothesisID ))


		jointResults = batchJointBranch(jointBranchJobs)

		"""
		results indexed by shoot ID

		branchResult["branchPathIDs"]
		branchResult["arcDists"]
		branchResult["controlSet"]
		branchResult["matchCounts"]
		branchResult["costSum"]
		branchResult["branchPoses_L"]
		branchResult["splices_G"]
		branchResult["trimmedPaths"]
		branchResult["longestPaths"] 


		#self.branchControlPoses
		#self.branchArcDists 
		#self.jointBranchEvaluations
		#self.branchArcDists[pathID].append(currDist)
		#self.branchControlPoses[pathID][currDist] = controlPose_P
			
		"""

		""" get the maximum value for each of our features """
		self.maxCost = -1e100
		self.maxMatchCount = -1000
		self.maxLandmarkCost = -1e100
		for k in range(len(jointResults)):
			branchResult = jointResults[k]

			totalMatchCount = branchResult["totalMatchCount"]
			totalCost = branchResult["totalCost"]
			landmarkCost = branchResult["landmarkCost"]

			if landmarkCost > self.maxLandmarkCost:
				self.maxLandmarkCost = landmarkCost

			if totalCost > self.maxCost:
				self.maxCost = totalCost

			if totalMatchCount > self.maxMatchCount:
				self.maxMatchCount = totalMatchCount

		""" add slight value so that maximum normalized landmark cost is greater than zero """
		self.maxLandmarkCost += 0.1

		for k in range(len(jointResults)):
			branchResult = jointResults[k]

			totalMatchCount = branchResult["totalMatchCount"]
			totalCost = branchResult["totalCost"]
			landmarkCost = branchResult["landmarkCost"]

			""" compute normalized cost and match count """
			if self.maxMatchCount > 0:
				normMatchCount = float(totalMatchCount) / float(self.maxMatchCount)
			else:
				normMatchCount = 1.0

			if self.maxMatchCount > 0:
				normCost = totalCost / self.maxCost
			else:
				normCost = 1.0

			if self.maxLandmarkCost > 0:
				normLandmark = (self.maxLandmarkCost-landmarkCost)/ self.maxLandmarkCost
			else:
				normLandmark = 1.0

			branchResult["normMatchCount"] = normMatchCount
			branchResult["normCost"] = normCost
			branchResult["normLandmark"] = normLandmark

			""" build indexing tuple """
			arcDists = branchResult["arcDists"]
			arcList = []
			for pathID in branchPathIDs:
				arcList.append(arcDists[pathID])

			arcTuple = tuple(arcList)

			self.jointBranchEvaluations[arcTuple] = branchResult

		
		if True:

			parentPathIDs = self.getParentHash()

			pylab.clf() 

			for arcTuple, branchResult in self.jointBranchEvaluations.iteritems():

				#modJuncPose = result["modJuncPose"]
				#modControlPose = result["modControlPose"]

				controlPoses_L = branchResult["controlSet"]
				controlPoses_G = computeGlobalControlPoses(controlPoses_L, parentPathIDs)

				branchPoses_L = branchResult["branchPoses_L"]


				allSplices = []
				splices_G = branchResult["splices_G"]
				allSplices = splices_G
				#for pathID in branchPathIDs:
				#	allSplices += splices_G[pathID]
				


				for path in allSplices:
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])

					pylab.plot(xP,yP,color='k', zorder=9, alpha=0.10)

				for k,path in self.trimmedPaths.iteritems():
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])

					pylab.plot(xP,yP, color = self.colors[k], linewidth=4)


				xC = []
				yC = []
				xB = []
				yB = []
				for pathID in branchPathIDs:
					controlPose_G = controlPoses_G[pathID]
					currFrame = Pose(controlPose_G)

					xC.append(controlPose_G[0])
					yC.append(controlPose_G[1])

					branchPose_L = branchPoses_L[pathID]
					branchPose_G = currFrame.convertLocalOffsetToGlobal(branchPose_L)

					xB.append(branchPose_G[0])
					yB.append(branchPose_G[1])
			
				pylab.scatter(xB, yB, color='y', zorder=8)
				pylab.scatter(xC, yC, color='k', zorder=8)


				"""

				pathID = result["pathID"]
				matchCount = result["matchCount"]
				lastCost = result["lastCost"]
				probVal = result["initProb"]
				junctionDetails = self.localLeaf2LeafPathJunctions[pathID]
				localPathSegs = junctionDetails["localSegments"]

				xP = []
				yP = []
				termPaths = result["termPaths"]
				for p1, p2 in termPaths:
					xP.append(p1[0])
					yP.append(p1[1])
					xP.append(p2[0])
					yP.append(p2[1])
				pylab.scatter(xP, yP, color='r', zorder=8)


				offsetOrigin1 = Pose(modControlPose)

				for k in range(len(localPathSegs)):
					localSeg = localPathSegs[k]
					xP1 = []
					yP1 = []
					for p in localSeg:
						p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
						xP1.append(p1[0])
						yP1.append(p1[1])

					#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.5)
					#if maxCost > 0:
					if maxMatchCount > 0:
						#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=probVal/maxProb)
						#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=lastCost/maxCost)
						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=float(matchCount)/float(maxMatchCount))
					else:
						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.10)
				"""



			pylab.title("nodeID: %d hyp: %d" % (len(self.nodePoses), self.hypothesisID))
			pylab.savefig("all_branch_plot_%04u_%02u_%04u.png" % (len(self.nodePoses), self.hypothesisID, self.tempCount) )
			self.tempCount += 1


	@logFunction
	def evaluateBranches(self, pathID, modControlPose, nodeID0, particleIndex):

		# pathID
		# parentID
		# curr path
		# parent path
		# origJuncPose
		# smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]
		# hypID


		if self.pathClasses[pathID]["parentID"] != None: 

			""" path information """
			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			pathSpline = SplineFit(self.paths[parentID])

			minDist, uVal, splinePoint = pathSpline.findClosestPoint(modControlPose)
			arcDist = pathSpline.dist_u(uVal)

			arcHigh = arcDist + self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
			arcLow = arcDist - self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)

			totalDist = pathSpline.dist_u(1.0)
		
			"""
			1) point that has the junction
			2) direction of junction
			3) set of long paths that include junction path  (only 2)
			"""

			""" initialize if this is the first time """
			branchSamples = []
			for k in range(self.NUM_BRANCHES):
				newArcDist = arcLow + k * self.DIV_LEN
				branchSample = self.getPrecomputedBranch(pathID, newArcDist)
				branchSamples.append(branchSample)

			""" get the maximum value for each of our features """
			maxCost = -1e100
			maxMatchCount = -1000
			maxAngDiff = -1e100
			maxDist = 10.0


			""" find maximum values for each metric """
			for k in range(len(branchSamples)):
				part = branchSamples[k]
				matchCount = part["matchCount"]
				lastCost = part["lastCost"]
				dist = part["distDisc"]
				angDiff = part["angDisc"]

				if matchCount > maxMatchCount:
					maxMatchCount = matchCount

				if lastCost > maxCost:
					maxCost = lastCost

				if dist > maxDist:
					maxDist = dist

				if angDiff > maxAngDiff:
					maxAngDiff = angDiff


			""" invert the distance, angDiff and cost features """
			totalProbSum = 0.0
			for k in range(len(branchSamples)):

				part = branchSamples[k]

				lastCost = part["lastCost"]
				dist = part["distDisc"]
				angDiff = part["angDisc"]
				matchCount = part["matchCount"]

				""" makes the maximum cost non-zero and inverted """
				matchCost = 0.1 + maxCost-lastCost

				nomDist = 0.0
				if maxDist > 0.0:
					nomDist = (maxDist-dist)/maxDist
				else:
					nomDist = 0.0


				""" angDiff feature times matchCount times dist """
				probVal = matchCount
				totalProbSum += probVal

				part2 = deepcopy(part)
				part2["lastCost"] = matchCost
				part2["distDisc"] = nomDist
				part2["initProb"] = probVal

				#part2 = (part[0], part[1], part[2], part[3], part[4], matchCost, nomDist, part[7], probVal, part[9], part[10], part[11])
				branchSamples[k] = part2

			""" normalize the probability values """
			maxProb = 0.0
			for k in range(len(branchSamples)):
				part = branchSamples[k]

				probVal = 0.0

				if totalProbSum > 0.0:

					probVal = part["initProb"] / totalProbSum

					if probVal > maxProb:
						maxProb = probVal
				else:
					probVal = 0.0

				part2 = deepcopy(part)
				part2["initProb"] = probVal

				branchSamples[k] = part2

				#part2 = (parentID, modJuncPose, newArcDist, pathID, matchCount1, lastCost1, distDisc, angDisc, initProb, newSplices, juncDiscAngle, juncDist)
				print "branchSample %04u %02u %02u %1.4f %03u %1.2f %1.2f %1.2f %1.2f %d %1.2f %1.2f" % (nodeID0, particleIndex, k, part2["modJuncPose"][2], part2["matchCount"], part2["lastCost"], part2["distDisc"], part2["angDisc"], part2["initProb"], len(part2["newSplices"]), part2["juncDiscAngle"], part2["juncDiscDist"])
			print "branchSample"

			""" cartesian distance """
			DISC_THRESH = 0.5

			""" 60 degree threshold """
			#ANG_THRESH = 0.523 # pi/6
			ANG_THRESH = 1.5708 # pi/2

			
			""" get splices for new branch position """
			""" reject branch locations whose branch angle discrepancy is too high """ 
			for k in range(len(branchSamples)):
				part = branchSamples[k]

				newProbVal = branchSamples[k]["initProb"] 

				juncDiscAngle = part["juncDiscAngle"]

				if fabs(juncDiscAngle) > ANG_THRESH:
					newProbVal = 0.0

				part2 = deepcopy(part)
				part2["initProb"] = newProbVal

				#part2 = (part[0], part[1], part[2], part[3], part[4], part[5], part[6], part[7], newProbVal, part[9], part[10], part[11])
				branchSamples[k] = part2


			junctionDetails = self.localLeaf2LeafPathJunctions[pathID]
			localPathSegs = junctionDetails["localSegments"]


			if True:
				pylab.clf() 
				for k,path in self.trimmedPaths.iteritems():
					print "path has", len(path), "points"
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])

					pylab.plot(xP,yP, color = self.colors[k], linewidth=4)

			if True:
				xP = []
				yP = []

				xP2 = []
				yP2 = []
				for part in branchSamples:

					modControlPose = part["modControlPose"]
					modJuncPose = part["modJuncPose"]

					controlOrigin1 = Pose(modControlPose)
					globalJuncPose = controlOrigin1.convertLocalOffsetToGlobal(modJuncPose)

					xP.append(globalJuncPose[0])
					yP.append(globalJuncPose[1])
					xP2.append(modControlPose[0])
					yP2.append(modControlPose[1])

					newSplices = part["newSplices"]

					for path in newSplices:
						xP1 = []
						yP1 = []
						for p in path:
							xP1.append(p[0])
							yP1.append(p[1])
						if maxProb > 0:
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part["initProb"]/maxProb)
						else:
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.10)


					"""
					offsetOrigin1 = Pose(modControlPose)

					for k in range(len(localPathSegs)):
						localSeg = localPathSegs[k]
						xP1 = []
						yP1 = []
						for p in localSeg:
							p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
							xP1.append(p1[0])
							yP1.append(p1[1])

						if maxProb > 0:
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part["initProb"]/maxProb)
						else:
							#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part[8])
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.10)
					"""


				pylab.scatter(xP, yP, color='y', zorder=8)
				pylab.scatter(xP2, yP2, color='k', zorder=8)
				pylab.title("nodeID: %d hyp: %d, branchSampleID: %d pathID: %d, localPathSegs %d" % (nodeID0, self.hypothesisID, particleIndex, pathID, len(localPathSegs)))
				
				pylab.savefig("bayes_plot_%04u_%02u_%03u_%04u.png" % (nodeID0, self.hypothesisID, particleIndex, self.tempCount) )
				self.tempCount += 1

		return branchSamples

	@logFunction
	def getJuncSample(self, pathID, sampleIndex):
		if self.pathClasses[pathID]["parentID"] == None:
			raise

		path2 = self.paths[pathID]

		origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
		origJuncPose[2] = 0.0

		origJuncOrigin = Pose(origJuncPose)

		globJuncPose = branchPart[1]
		offsetOrigin1 = Pose(globJuncPose)
		
		path2 = self.paths[pathID]

		localPath2 = []
		for p in path2:
			p1 = origJuncOrigin.convertGlobalToLocal(p)
			localPath2.append(p1)

		samplePath2 = []
		for p in localPath2:
			p1 = offsetOrigin1.convertLocalToGlobal(p)
			samplePath2.append(p1)

		return samplePath2


	""" update landmark data structures """
	def updateLandmarks(self):


		BEND_THRESH = 0.9
		ARCH_THRESH = 0.05
		BLOOM_THRESH = 0.05


		pathIDs = self.getPathIDs()
		controlPoses = {}
		for pathID in pathIDs:
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]

		""" store landmark points into their proper shoots """
		self.localLandmarks = {}
		self.nodeLandmarks = {}
		for pathID in pathIDs:
			self.localLandmarks[pathID] = []
			self.nodeLandmarks[pathID] = {}

			nodeSet = self.getNodes(pathID)
			for nodeID in nodeSet:

				currFrame_L = Pose(controlPoses[pathID])
				nodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
				poseFrame_L = Pose(nodePose_L)

				spatialFeature = self.poseData.spatialFeatures[nodeID][0]

				landmarkPoint_N = None
				landmarkPoint_L = None



				self.nodeLandmarks[pathID][nodeID] = None

				if spatialFeature["bloomPoint"] != None:
					print self.hypothesisID, "adding bloomPoint for node", nodeID, "in path", pathID
					landmarkPoint_N = spatialFeature["bloomPoint"]
					landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
					self.localLandmarks[pathID].append((landmarkPoint_L, BLOOM_THRESH, "bloomPoint"))
					self.nodeLandmarks[pathID][nodeID] = (landmarkPoint_N, BLOOM_THRESH, "bloomPoint")

				elif spatialFeature["archPoint"] != None:
					print self.hypothesisID, "adding archPoint for node", nodeID, "in path", pathID
					landmarkPoint_N = spatialFeature["archPoint"]
					landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
					self.localLandmarks[pathID].append((landmarkPoint_L, ARCH_THRESH, "archPoint"))
					self.nodeLandmarks[pathID][nodeID] = (landmarkPoint_N, ARCH_THRESH, "archPoint")

				elif spatialFeature["inflectionPoint"] != None:
					print self.hypothesisID, "adding inflectionPoint for node", nodeID, "in path", pathID
					landmarkPoint_N = spatialFeature["inflectionPoint"]
					landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
					self.localLandmarks[pathID].append((landmarkPoint_L, BEND_THRESH, "bendPoint"))
					self.nodeLandmarks[pathID][nodeID] = (landmarkPoint_N, BEND_THRESH, "bendPoint")

				#if landmarkPoint_L != None:
				#	#self.localLandmarks[pathID].append(landmarkPoint_L)
				#	self.localLandmarks[pathID].append((landmarkPoint_L, 0.3, "bloomPoint"))

				#self.nodeLandmarks[pathID][nodeID] = landmarkPoint_N


		print self.hypothesisID, "landmarks:", self.localLandmarks

	
				
	""" generate the paths from the node and pose information """
	@logFunction
	def generatePaths(self):


		self.updateLandmarks()

	
		" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
		self.paths = {}
		self.hulls = {}
		self.localLeaf2LeafPathJunctions = {}
		self.localPaths = {}
		self.localHulls = {}
		
		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			print "computing path for node set", pathID, ":", self.getNodes(pathID)
			localResults = computeShootSkeleton(self.poseData,
										pathID,
										self.pathClasses[pathID]["localJunctionPose"],
										self.getNodes(pathID),
										self.pathClasses[pathID]["localNodePoses"],
										self.localLandmarks[pathID],
										self.hypothesisID,
										self.colors[pathID],
										self.topCount,
										plotIter = True)

			self.localLeaf2LeafPathJunctions[pathID] = localResults[0]
			self.localPaths[pathID] = localResults[1]
			self.localHulls[pathID] = localResults[2]




			self.topCount += 1
			
		localSkeletons = {}
		controlPoses = {}
		junctionPoses = {}
		for pathID in pathIDs:
			localSkeletons[pathID] = self.localLeaf2LeafPathJunctions[pathID]["skeletonGraph"]
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

		parentPathIDs = self.getParentHash()
		
		finalPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		print "controlPoses:", controlPoses
		print "finalPoses:", finalPoses
		print "junctionPoses:", junctionPoses


		""" convert local to max likelihood global """
		for pathID, localPath in self.localPaths.iteritems():
			shootFrame = Pose(finalPoses[pathID])

			globalPath = []
			for p in localPath:
				globalP = shootFrame.convertLocalToGlobal(p)
				globalPath.append(globalP)

			self.paths[pathID] = globalPath

		""" convert local to max likelihood global """
		for pathID, localHull in self.localHulls.iteritems():
			shootFrame = Pose(finalPoses[pathID])

			globalHull = []
			for p in localHull:
				globalP = shootFrame.convertLocalToGlobal(p)
				globalHull.append(globalP)

			self.hulls[pathID] = globalHull


		#self.paths[pathID] = results[1]
		#self.hulls[pathID] = results[2]

		self.trimmedPaths = self.trimPaths(self.paths)


		#self.spliceSkeleton = spliceSkeletons(localSkeletons, controlPoses, junctionPoses, parentPathIDs)
		self.spliceSkeleton = spliceSkeletons(localSkeletons, finalPoses, junctionPoses, parentPathIDs)

		" for each path, attempt to join with its parent path "
		self.joins = []
		self.junctions = {}
		for pathID in pathIDs:

			cPath = self.getPath(pathID)
			
			parentPathID = cPath["parentID"]
			
			" parent does not concern us "
			if parentPathID == None:
				continue

			junctionPose = self.getGlobalJunctionPose(pathID)
			junctionPoint = [junctionPose[0],junctionPose[1]]

			print "node globJuncPose:",  junctionPose

			path1 = self.trimmedPaths[pathID]
			path2 = self.trimmedPaths[parentPathID]

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
			
			self.joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])

			" get junctions " 
			branchNodeID = cPath["branchNodeID"]
			
			self.junctions[pathID] = [branchNodeID, junctionPose, (parentPathID,minI2), path2[minI2], minI1]


		" create a tree with the node IDs and then stitch them together with joins "
		self.pathGraph = graph.graph()
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for k in range(len(path)):
				self.pathGraph.add_node((pathID, k), path[k])

			for k in range(len(path)-1):
				p1 = path[k]
				p2 = path[k+1]
				weight = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

				self.pathGraph.add_edge((pathID, k), (pathID, k+1), wt=weight)
			
			" parent does not concern us "
			cPath = self.getPath(pathID)			
			parentPathID = cPath["parentID"]
			
		" join with the junction in between the join points "
		for k in range(len(self.joins)):
			join = self.joins[k]

			pID1, k1 = join[0]
			pID2, k2 = join[1]

			p1 = self.trimmedPaths[pID1][k1]
			p2 = self.trimmedPaths[pID2][k2]
			weight = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

			self.pathGraph.add_edge(join[0], join[1], wt=weight)



		self.topDict = {}
		
		
		" get terminals "

		self.terminals = {}
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			
			if len(path) > 0:
				if pathID == 0:
					
					dist1 = sqrt((self.trimmedPaths[pathID][0][0] - self.rootPoint[0])**2 + (self.trimmedPaths[pathID][0][1] - self.rootPoint[1])**2)
					dist2 = sqrt((self.trimmedPaths[pathID][len(path)-1][0] - self.rootPoint[0])**2 + (self.trimmedPaths[pathID][len(path)-1][1] - self.rootPoint[1])**2)

					if dist1 < dist2:
						self.terminals[0] = [(pathID,0), path[0]]
						self.terminals[1] = [(pathID,len(path)-1), path[len(path)-1]]
						self.topDict["t%u" % 0] = (pathID,0)
						self.topDict["t%u" % 1] = (pathID,len(path)-1)
						self.rootPoint = self.trimmedPaths[pathID][0]

					else:
						self.terminals[1] = [(pathID,0), path[0]]
						self.terminals[0] = [(pathID,len(path)-1), path[len(path)-1]]
						self.topDict["t%u" % 1] = (pathID,0)
						self.topDict["t%u" % 0] = (pathID,len(path)-1)
						self.rootPoint = self.trimmedPaths[pathID][len(path)-1]
					
				else:
					
					self.topDict["j%u" % pathID] = self.junctions[pathID][2]
					minI1 = self.junctions[pathID][4]
					
					" determine which side is junction and which side is terminal"
					if minI1 > len(path)-1 - minI1:    
						self.terminals[pathID+1] = [(pathID,0), path[0]]
						self.topDict["t%u" % (pathID+1)] = (pathID,0)
					else:
						self.terminals[pathID+1] = [(pathID,len(path)-1), path[len(path)-1]]
						self.topDict["t%u" % (pathID+1)] = (pathID,len(path)-1)


		#if len(self.trimmedPaths[pathID]) > 0:
		""" don't care about return values, are stored as object member variable """
		self.computeAllSplices(plotIter = True)
		self.isChanged = False

	@logFunction
	def getOrderedOverlappingPaths(self, nodeID):

		nodePose = self.getNodePose(nodeID)
		splicePaths, spliceTerms, splicePathIDs = self.getSplicesByNearJunction(nodePose)

		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		medialSpline2 = SplineFit(medial2, smooth=0.1)
		medial2_vec = medialSpline2.getUniformSamples()

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.getNodePose(nodeID)

		poseOrigin2 = Pose(estPose2)
		globalMedial2 = []
		for p in medial2:
			globalMedial2.append(poseOrigin2.convertLocalToGlobal(p))

		resultSet = []

		for k in range(len(splicePaths)):
			
			print "comparing node", nodeID, "to splice", splicePathIDs[k], spliceTerms[k],  "plot", self.multiDepCount
			
			currPath = splicePaths[k]
			pathIDs = splicePathIDs[k]
			currPath = orientPath(currPath, globalMedial2)
			results = getMultiDeparturePoint(currPath, medial2_vec, estPose2, estPose2, pathIDs, nodeID, pathPlotCount = self.multiDepCount, hypID = self.hypothesisID, plotIter = False) 

			self.multiDepCount += 1

			resultSet.append(results+(k,))

			"departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2"

		for result in resultSet:
			print "result:", result
		
		""" overlapSum secondary, contigFrac primary """
		resultSet = sorted(resultSet, key=itemgetter(13), reverse=False)
		resultSet = sorted(resultSet, key=itemgetter(12), reverse=True)
		
		print "getOrderedOverlappingPaths() sorted node", nodeID, ":"
		for result in resultSet:
			print result[15], result[12], result[13]

		# FIXME:  temporary forced result 
		if len(resultSet) > 1:
			contigFrac0 = resultSet[0][12]
			contigFrac1 = resultSet[1][12]
			spliceIndex0 = resultSet[0][15]
			spliceIndex1 = resultSet[1][15]

			diff = fabs(contigFrac0-contigFrac1)

			if diff < 0.05:
				if len(splicePathIDs[spliceIndex0]) < len(splicePathIDs[spliceIndex1]): 
					result = resultSet[0]	 
					spliceIndex = result[15]
					overlappedPathIDs = splicePathIDs[spliceIndex]
				else:
					result = resultSet[1]	 
					spliceIndex = result[15]
					overlappedPathIDs = splicePathIDs[spliceIndex]
			else:
				result = resultSet[0]	 
				spliceIndex = result[15]
				overlappedPathIDs = splicePathIDs[spliceIndex]

		else:
			result = resultSet[0]	 
			spliceIndex = result[15]
			overlappedPathIDs = splicePathIDs[spliceIndex]

		"""
		for k in range(len(resultSet)):
			result = resultSet[k]	 
			spliceIndex = result[15]
			overlappedPathIDs = splicePathIDs[spliceIndex]
			
			print "checking overlap of", overlappedPathIDs
			isNotOverlapped = False
			for pathID in overlappedPathIDs:

				#medial2 = self.poseData.medialAxes[nodeID]
				#estPose2 = self.nodePoses[nodeID]
				sum1 = getOverlapCondition(self.poseData.medialAxes[nodeID], self.getNodePose(nodeID), self.trimmedPaths[pathID], nodeID, plotIter = True, overlapPlotCount = self.overlapPlotCount)
				self.overlapPlotCount += 1

				#sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID, plotIter = False)

				print "pathID", pathID, "returned", sum1, "for splice", spliceIndex
				if sum1 > 1e10:
					isNotOverlapped = True
					
			if not isNotOverlapped:
				print "selected", k, "'th result"
				break
		"""
			
		
		orderedPathIDs = self.getPathOrdering(nodeID, overlappedPathIDs)

		print "returned ordered", orderedPathIDs, "from unordered", overlappedPathIDs
		
		return orderedPathIDs
		

	#@logFunction
	def getSplicesByNearJunction(self, initPose2):
		
		allSplices, terminals, junctions = self.getAllSplices(plotIter = False)

		print "junctions:", junctions
		print "initPose2:", initPose2
		
		junctions2 = []
		for pathID, params in junctions.iteritems():
			junctionPose = params[1]
			
			print "checking", pathID, params
			
			dist2 = sqrt((junctionPose[0]-initPose2[0])**2 + (junctionPose[1]-initPose2[1])**2)
			
			if dist2 < 3.0:
				junctions2.append((pathID,params[2][0]))
			
		"self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]"
		
		closePathID2 = self.getClosestPath(initPose2)

		print "closePathID2:", closePathID2
		
		pathSet2 = [closePathID2]

		junctionSet = junctions2
		junctionSet = list(set(junctionSet))

		print "junctions2:", junctions2

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
			parentID = self.getParentPathID(pathID)

			#if parentID != None and not parentID in pathSet2:
			if parentID != None:
				termSet2.append(junctions[pathID][2])
			

		

		splicePaths = []
		spliceTerms = []
		splicePathIDs = []

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
						spliceTerms.append(termPath)
						splicePathIDs.append(sPath['orderedPathIDs'])

		return splicePaths, spliceTerms, splicePathIDs


	 
	 
	@logFunction
	def delNode(self, nodeID, pathID):
		print "deleting node", nodeID, "from path", pathID
		self.pathClasses[pathID]["nodeSet"].remove(nodeID)
		del self.pathClasses[pathID]["localNodePoses"][nodeID]
		
		self.isChanged = True		

	@logFunction
	def moveNode(self, nodeID, splicePathIDs):
		" change the shoot membership based on fitted splice "

		particleDist = self.poseParticles["snapshots2"][0]

		" determine which paths are leaves "
		pathIDs = self.getPathIDs()

		" remove node from other shoots first"
		for pathID in pathIDs:
			if nodeID in self.pathClasses[pathID]["nodeSet"]:
				print "removing node", nodeID, "from path", pathID
				self.pathClasses[pathID]["nodeSet"].remove(nodeID)
				del self.pathClasses[pathID]["localNodePoses"][nodeID]

				for p in particleDist:
					p.delNode(nodeID, pathID)


		isAParent = {}
		for k in pathIDs:
			isAParent[k] = False
		for k in splicePathIDs:
			currPath = self.getPath(k)
			currParent = currPath["parentID"]
			if currParent != None:
				isAParent[currParent] = True

		"add nodes to paths that are the leaves "
		for pathID in splicePathIDs:
			if not isAParent[pathID]:
				print "adding node", nodeID, "to path", pathID
				self.pathClasses[pathID]["nodeSet"].append(nodeID)
				#self.addNode(nodeID,pathID)

				""" convert the controlPose to coordinates local to the parent frame """
				parentPathIDs = self.getParentHash()

				""" control state of maximum likelihood particle determine's location of shoot frame """
				controlPoses = self.getControlPoses()
				globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)
				shootControlPose_G = globalControlPoses_G[pathID]
				currFrame = Pose(shootControlPose_G)

				""" convert from global to local """
				nodePose_G = self.getNodePose(nodeID)
				nodePose_C = currFrame.convertGlobalPoseToLocal(nodePose_G)
				self.pathClasses[pathID]["localNodePoses"][nodeID] = nodePose_C

				for p in particleDist:
					p.addNode(nodeID, pathID, nodePose_C)

		print "shoot", pathID, "has nodes", self.pathClasses[pathID]["nodeSet"]

		self.isChanged = True		

	@logFunction
	def addNode(self, nodeID, pathID):


		print "adding node", nodeID, "to path", pathID
		if not nodeID in self.pathClasses[pathID]["nodeSet"]:
			self.pathClasses[pathID]["nodeSet"].append(nodeID)


			""" convert the controlPose to coordinates local to the parent frame """
			parentPathIDs = self.getParentHash()
			controlPoses = self.getControlPoses()
			globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)
			shootControlPose_G = globalControlPoses_G[pathID]
			currFrame = Pose(shootControlPose_G)

			nodePose_G = self.getNodePose(nodeID)

			nodePose_C = currFrame.convertGlobalPoseToLocal(nodePose_G)
			
			self.pathClasses[pathID]["localNodePoses"][nodeID] = nodePose_C

		else:
			print "node", nodeID, "already present in path", pathID

		updateCount = self.poseParticles["updateCount"] 
		particleDist = self.poseParticles["snapshots2"][0]

		for p in particleDist:
			nodePose_C = self.pathClasses[pathID]["localNodePoses"][nodeID]
			p.addNode(nodeID, pathID, nodePose_C)

		print pathID, "has nodes", self.pathClasses[pathID]["nodeSet"]


		self.isChanged = True

	@logFunction
	def delPath(self, pathID, mergeTargetID):
		
		print "deleting path", pathID, "merging to path", mergeTargetID
		try: 
			del self.pathClasses[pathID]
			del self.pathTermsVisited[pathID]
			
			
			for k, pathClass in self.pathClasses.iteritems():
				
				if pathClass["parentID"] == pathID:
					pathClass["parentID"] = mergeTargetID
					
		except:
			pass
	
		self.isChanged = True		

	@logFunction
	def addPath(self, branchNodeID, localDivergencePose_R):

		"""

		Create a new shoot branch given the following starting information:
		
		1) the parent shoot ID
		2) the spatial node ID that begins the new shoot with a detected divergence
		3) the pose of the divergence point in the local coordinates of the spatial node


		From this compute the following information:

		1) the junction point, representing the center of a hypothetical junction
		2) a control point on the parent shoot that allows a distribution of possible branch locations
		3) the computed shoot branch parameters for each currently active pose particle

		"""


		print "addPath(", self.shootIDs, branchNodeID, localDivergencePose_R
		
		""" the raw pose and posture pose """
		nodePose_G = self.getNodePose(branchNodeID)
		rawPose_G = self.getNodeRawPose(branchNodeID)

		""" the localDivergencePose computed from raw local coordinates to global coordinates """
		rawPoseFrame = Pose(rawPose_G)
		globalJunctionPose_G = rawPoseFrame.convertLocalOffsetToGlobal(localDivergencePose_R)
	
		
		""" pathIDs from self.localLeaf2LeafPathJunctions structure because recent shoots may not be generated yet """
		#allPathIDs = self.getPathIDs()
		allPathIDs = self.localLeaf2LeafPathJunctions.keys()
		

		""" new shoot ID """
		#newPathID = self.pathIDs
		newPathID = self.shootIDs

		""" posture curve of local spatial map """
		medial0_L = self.poseData.medialAxes[branchNodeID]

		""" medial axis converted to global coordinates from local posture coordinates """
		nodeFrame = Pose(nodePose_G)
		globalMedial0_G = []
		for p in medial0_L:
			globalMedial0_G.append(nodeFrame.convertLocalToGlobal(p))

		#trimmedParent_G = self.trimmedPaths[parentID]

		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

		localPathSegsByID = {}
		for pathID in allPathIDs:
			localPathSegs = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]
			localPathSegsByID[pathID] = localPathSegs


		#self.localPaths[pathID] = localResults[1]
		#localPathSegsByID[newPathID] 


		controlPose_G, controlParentID, branchPose_G = getInitSkeletonBranchPoint(globalJunctionPose_G, newPathID, globalMedial0_G, parentPathIDs, localPathSegsByID, self.localPaths, globalControlPoses_G, plotIter = False, hypothesisID = self.hypothesisID, nodeID = branchNodeID)

		#try:
		#	branchPose_G, isNoDiverge = getSkeletonBranchPoint(oldBranchPose_G, pathID, parentPathIDs, localPathSegsByID, localPaths, controlPoses_G)
		#except:
		#	branchPose_G, isNoDiverge = getSkeletonBranchPoint(oldBranchPose_G, pathID, parentPathIDs, localPathSegsByID, localPaths, controlPoses_G, angThresh = pi)

		#""" compute branch point of posture curve diverging from parent shoot """
		#newGlobJuncPose1_G, controlPoint1_G, angDeriv1 = getBranchPoint(globalJunctionPose_G, parentID, newPathID, trimmedParent_G, globalMedial0_G, plotIter = False, hypothesisID = self.hypothesisID, nodeID = branchNodeID)

		#""" compute branch point of parent shoot diverging from posture curve """
		#newGlobJuncPose2_G, controlPoint2_G, angDeriv2 = getBranchPoint(globalJunctionPose_G, newPathID, parentID, globalMedial0_G, trimmedParent_G, plotIter = False, hypothesisID = self.hypothesisID, nodeID = branchNodeID)

		#""" the control point that is the flattest determines which curve is the diverging one """
		#""" the curve that is sharpest at the control point is the diverging curve """
		#if angDeriv1 > angDeriv2:
		#	newGlobJuncPose_G = newGlobJuncPose1_G
		#	angDeriv = angDeriv1
		#else:
		#	newGlobJuncPose_G = newGlobJuncPose2_G
		#	angDeriv = angDeriv2


		#""" regardless, take the branch angle of the child from parent shoot branch point """
		#print "branchPointProb:", angDeriv1, angDeriv2, newGlobJuncPose1_G, newGlobJuncPose2_G
		#newGlobJuncPose_G[2] = newGlobJuncPose1_G[2]
		#print "setting new", newGlobJuncPose_G
		
		newGlobJuncPose_G = branchPose_G


		""" choose a control point that is common between parent and child shoots """
		#commonU1, commonU2, cPoint1_G, cPoint2_G = selectCommonOrigin(trimmedParent_G, globalMedial0_G)	
		#controlPose_G = [cPoint1_G[0],cPoint1_G[1], 0.0]

		""" convert the controlPose to coordinates local to the parent frame """
		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)
		
		""" compute the control pose relative to the parent """
		parentControlPose_G = globalControlPoses_G[controlParentID]
		parentFrame = Pose(parentControlPose_G)
		controlPose_P = parentFrame.convertGlobalPoseToLocal(controlPose_G)

		nodePose_G = self.getNodePose(branchNodeID)
		currFrame = Pose(controlPose_G)
		nodePose_C = currFrame.convertGlobalPoseToLocal(nodePose_G)
		junctionPose_C = currFrame.convertGlobalPoseToLocal(newGlobJuncPose_G)

		""" basic shoot data structure """
		self.pathClasses[newPathID] = {"parentID" : controlParentID,
						"branchNodeID" : branchNodeID,
						"localDivergencePose" : localDivergencePose_R, 
						"localJunctionPose" : junctionPose_C, 
						"sameProb" : {},
						"nodeSet" : [branchNodeID,],
						"localNodePoses" : {branchNodeID : nodePose_C},
						"globalJunctionPose" : newGlobJuncPose_G,
						"controlPose" : controlPose_P }		

		print "new pathClass:", self.pathClasses[newPathID]


		self.pathTermsVisited[newPathID] = False
		
		self.pathIDs += 1
		self.shootIDs += 1


		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]

		""" FIXME:  Full path spline instead of just relevant section.  Does this always work? """
		#parentShoot_G = self.paths[parentID]
		parentShoot_L = self.localPaths[controlParentID]
		parentShoot_G = []
		for p in parentShoot_L:
			parentShoot_G.append(parentFrame.convertLocalToGlobal(p))

		pathSpline_G = SplineFit(parentShoot_G)


		""" add new path to particles """
		for part in particleDist2:

			nodePose_G = self.getNodePose(branchNodeID)
			rawPose_G = self.getNodeRawPose(branchNodeID)

			""" get local transform of GPAC to raw pose """
			nodeFrame = Pose(nodePose_G)
			transform_N_to_R = nodeFrame.convertGlobalPoseToLocal(rawPose_G)
			
			""" go back and convert this from GPAC pose to raw pose """
			particlePose_G = deepcopy(part.pose0)
			particlePose_G[2] = nodePose_G[2]
			particlePoseFrame = Pose(particlePose_G)
			particleRawFrame_G = particlePoseFrame.convertLocalOffsetToGlobal(transform_N_to_R)

			""" location of junction point from raw to global in the particle pose """
			rawParticlePoseFrame = Pose(particleRawFrame_G)
			particleJunctionPose_G = rawParticlePoseFrame.convertLocalOffsetToGlobal(localDivergencePose_R)


			""" compute the particle control pose based on the global location of the particle pose """
			controlPose_L = nodeFrame.convertGlobalPoseToLocal(controlPose_G)
			particleControlPose_G = particlePoseFrame.convertLocalOffsetToGlobal(controlPose_L)

			#print "control poses:", nodePose_G, particlePose_G, controlPose_G, particleControlPose_G
		
			""" location on the parent shoot closest to the locally computed branch point """
			minDist, uVal, newJunctionPose_G = pathSpline_G.findClosestPoint(particleJunctionPose_G)

			""" angle of junction corrected to the original orientation """
			modJunctionPose_G = copy(newJunctionPose_G)
			modJunctionPose_G[2] = globalJunctionPose_G[2]

			""" get the arc distance of the control point """
			minDist, controlUVal, newControlPose_G = pathSpline_G.findClosestPoint(particleControlPose_G)
			arcDist = pathSpline_G.dist_u(controlUVal)

			""" compute the rails of the branch point distribution """
			arcHigh = arcDist + self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
			arcLow = arcDist - self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)

			""" compute the sample points for the branch point distribution """
			arcDists = []
			for k in range(self.NUM_BRANCHES):
				newArcDist = arcLow + k * self.DIV_LEN
				arcDists.append((newPathID, newArcDist))

			""" convert the controlPose to coordinates local to the parent frame """
			""" compute the control pose relative to the parent """
			parentPathIDs = self.getParentHash()
			controlPoses = self.getControlPoses()
			globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)
			parentControlPose_G = globalControlPoses_G[controlParentID]
			parentFrame = Pose(parentControlPose_G)

			currFrame = Pose(newControlPose_G)
			localJunctionPose_C = currFrame.convertGlobalPoseToLocal(modJunctionPose_G)

			controlPoses_P = []
			for k in range(self.NUM_BRANCHES):
				globalControlPose_G = pathSpline_G.getPointOfDist(arcDists[k][1])
				localControlPose_P = parentFrame.convertGlobalPoseToLocal(globalControlPose_G)
				controlPoses_P.append(localControlPose_P)

			localParticleControlPose_P = parentFrame.convertGlobalPoseToLocal(particleControlPose_G)

			""" add the details of this junction given particle pose is true """
			part.addPath(newPathID, controlParentID, branchNodeID, nodePose_C, localDivergencePose_R, modJunctionPose_G, localJunctionPose_C, localParticleControlPose_P, self.NUM_BRANCHES, arcDists, controlPoses_P)

		print "newPath", newPathID, "=", self.pathClasses[newPathID]

		self.updateLandmarks()

		self.isChanged = True

		return newPathID
 
	@logFunction
	def trimPaths(self, foo = None):

		paths = self.paths
		localPaths = self.localPaths
		trimmedPaths = {}
		branchPoses_G = {}
		localTrimmedPaths = {}
		branchDiverges = {}
		
		print "path lengths:"
		for k,v in paths.iteritems():
			print k, len(v)
		
		
		if len(paths) <= 1:
			for k, v in paths.iteritems():
				trimmedPaths[k] = v
				localTrimmedPaths[k] = v
				
			self.trimmedPaths = trimmedPaths
			self.localTrimmedPaths = localTrimmedPaths
			return trimmedPaths

		" path0 is root, and path1 is a child of path0 "

		" information on parentage is required to determine which path belongs to which "

		" global junction point for parent 0 and child 1 "


		pathIDs = self.getPathIDs()

		childParents = {}
		isParentComputed = {}
		for pathID in pathIDs:
			isParentComputed[pathID] = False
			parentPathID = self.pathClasses[pathID]["parentID"]
			childParents[pathID] = parentPathID

		isParentComputed[0] = True

		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		localPathSegsByID = {}
		for pathID in pathIDs:
			#localSkeletons[pathID] = self.localLeaf2LeafPathJunctions[pathID]["skeletonGraph"]
			#controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			#junctionPoses[pathID] = self.pathClasses[pathID]["localJunctionPose"]
			localPathSegsByID[pathID] = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]

		while False in isParentComputed.values():
			print "trim parents:", isParentComputed, childParents

			for pathID in pathIDs:
		
				if isParentComputed[pathID]:
				
					if childParents[pathID] == None:
						trimmedPaths[pathID] = deepcopy(paths[pathID])
						localTrimmedPaths[pathID] = deepcopy(paths[pathID])
						branchDiverges[0] = True

					else:
						parentPathID = childParents[pathID]
						childPathID = pathID

						localJunctionPose = self.pathClasses[childPathID]["localJunctionPose"]
						localControlPose = self.pathClasses[childPathID]["controlPose"]
						globalJunctionPose = self.getGlobalJunctionPose(childPathID)
						globalControlPose = globalControlPoses[pathID]

						localNewPath3, branchPose_L, localParticlePath, isNoDiverge = trimBranch(childPathID, parentPathID, localControlPose, localJunctionPose, localPathSegsByID, self.localPaths, parentPathIDs, globalControlPoses, plotIter=True, hypothesisID=self.hypothesisID, nodeID=(self.poseData.numNodes), arcDist = childPathID)


						currFrame = Pose(globalControlPose)
						branchPose_G = currFrame.convertLocalOffsetToGlobal(branchPose_L)

						offsetOrigin1 = Pose(globalControlPose)
						globalPath3 = []
						for p in localNewPath3:
							globalPath3.append(offsetOrigin1.convertLocalToGlobal(p))

						trimmedPaths[pathID] = deepcopy(globalPath3)
						localTrimmedPaths[pathID] = deepcopy(localNewPath3)
						branchPoses_G[pathID] = branchPose_G
						branchDiverges[pathID] = not isNoDiverge

					for childID in pathIDs:
						if self.pathClasses[childID]["parentID"] == pathID:
							isParentComputed[childID] = True

		self.branchPoses_G = branchPoses_G
		self.trimmedPaths = trimmedPaths
		self.localTrimmedPaths = localTrimmedPaths
		self.branchDiverges = branchDiverges
			
		return trimmedPaths


	@logFunction
	def pathTermVisited(self, termID):
		
		try:
			self.pathTermsVisited[termID]
		except:
			pass
		else:
			self.pathTermsVisited[termID] = True

	@logFunction
	def getPathTermsVisited(self):
		return self.pathTermsVisited

	@logFunction
	def resetTerms(self):
		
		for k, v in self.pathTermsVisited.iteritems():
			self.pathTermsVisited[k] = False

	@logFunction
	def getPathTerms(self):

		terms = {}

		pathIDs = self.getPathIDs()

		for pathID in pathIDs:

			pID, pathI = self.topDict["t%u" % (pathID+1)]
			path = self.trimmedPaths[pID]
			terms[pathID] = [path[pathI][0],path[pathI][1], 0.0]

		""" root terminal point """
		pID, pathI = self.topDict["t%u" % 0]
		path = self.trimmedPaths[pID]
		terms[-1] = [path[pathI][0],path[pathI][1], 0.0]
			
		print "returning path terms:", terms

		#self.topDict["t%u" % (pathID+1)] = (pathID,0)
		return terms

	@logFunction
	def getClosestPath(self, locPose):
	
		pathIDs = self.getPathIDs()
		
		minDist = 1e100
		minPathID = 0
		for pathID in pathIDs:
			
			pathSeg = self.trimmedPaths[pathID]
			
			for p in pathSeg:
				dist = sqrt((p[0]-locPose[0])**2+(p[1]-locPose[1])**2)
				
				if dist < minDist:
					minDist = dist
					minPathID = pathID
		
		return minPathID

	@logFunction
	def getPathPath(self, startPathID, endPathID):
		" given start pathID and end pathID, what is the series of paths I need to follow to get from start to end "

		if startPathID == endPathID:
			return [startPathID]
		
		pathGraph = graph.graph()

		pathIDs = self.getPathIDs()

		for pathID in pathIDs:
			pathGraph.add_node(pathID, [])
		
		for pathID in pathIDs:	  
			parentPathID = self.getPath(pathID)["parentID"]
			if parentPathID != None:
				pathGraph.add_edge(pathID, parentPathID)
		
		shortestPathSpanTree, shortestDist = pathGraph.shortest_path(endPathID)
		nextPathID = shortestPathSpanTree[startPathID]

		splicedPaths = [startPathID, nextPathID]
		while nextPathID != endPathID:
			nextPathID = shortestPathSpanTree[nextPathID]
			splicedPaths.append(nextPathID)

		return splicedPaths

	@logFunction
	def getPathOrdering(self, nodeID, pathIDs):

		plotIter = False

		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.getNodePose(nodeID)
			
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
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			 
		
		pathSelect = []
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

				
				except:
					pass

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

	@logFunction
	def getPathDeparture(self, path1, path2, pathID1, pathID2, termPathIDs1 = [], termPathIDs2 = [], plotIter = False):
		
		isExist1 = False
		isInterior1 = False
		departurePoint1 = 0
		angle1 = 0.0
		isExist2 = False
		isInterior2 = False
		departurePoint2 = 0
		angle2 = 0.0
		
		if len(path1) == 0:
			return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
				
				
				
				
		orientedPath1 = orientPath(path1, path2)		   


		"""Assumption:  one section of the medial axis is closely aligned with the path """
		
		pathSpline2 = SplineFit(path2, smooth=0.1)
		pathPoints2 = pathSpline2.getUniformSamples()
		print "path2:", path2[0], pathPoints2[0], path2[-1], pathPoints2[-1]

		pathSpline1 = SplineFit(orientedPath1, smooth=0.1)
		pathPoints1 = pathSpline1.getUniformSamples()
		print "path1:", orientedPath1[0], pathPoints1[0], orientedPath1[-1], pathPoints1[-1]


		angle1, angle2 = getTipAngles(pathPoints2)

		distSum = 0.0
		contigCount = 0
		maxContig = 0
		distances = []
		indices = []
		for i in range(0,len(pathPoints2)):
			p_2 = pathPoints2[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
			distances.append(minDist)
			indices.append(i_1)
			distSum += minDist
			
			if minDist < 0.2:
				contigCount += 1
				if contigCount > maxContig:
					maxContig = contigCount
			else:
				contigCount = 0
		
		overlapSum = distSum / float(len(pathPoints2))
		
		print "maxContig,overlapSum:", maxContig, overlapSum

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


		frontAngleRefI = frontDepI
		while distances[frontAngleRefI] < 0.1:
			frontAngleRefI -= 1
			
			if frontAngleRefI < 0:
				frontAngleRefI = 0
				break
						
		forePathIndex = indices[frontAngleRefI]
		forePathAngle = pathPoints1[forePathIndex][2]
		forePathAngle = normalizeAngle(forePathAngle + pi)

		foreDiffAngle = diffAngle(angle1, forePathAngle)
		
		newFrontDepI = 4
		while fabs(diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]+pi), forePathAngle)) > pi/6.0:
			
			if newFrontDepI < frontDepI:
				newFrontDepI += 1
			else:
				break
		
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

		backAngleRefI = backDepI
		while distances[backAngleRefI] < 0.1:
			backAngleRefI += 1
			
			if backAngleRefI >= len(distances):
				backAngleRefI = len(distances)-1
				break
		
		backPathIndex = indices[backAngleRefI]
		backPathAngle = pathPoints1[backPathIndex][2]
		
		backDiffAngle = diffAngle(angle2, backPathAngle)
		
		newBackDepI = len(distances)-5
		while fabs(diffAngle(pathPoints2[newBackDepI][2], backPathAngle)) > pi/6.0:
			
			if newBackDepI > backDepI:
				newBackDepI -= 1
			else:
				break		 


		foreAngDiffs = []
		for i in range(0,len(pathPoints2)):
			foreAngDiffs.append(fabs(diffAngle(normalizeAngle(pathPoints2[i][2] + pi), forePathAngle)))

		backAngDiffs = []
		for i in range(0,len(pathPoints2)):
			backAngDiffs.append(fabs(diffAngle(pathPoints2[i][2], backPathAngle)))
			
		"reset to the maximum distances "
		maxFront = distances[0]
		maxBack = distances[-1]
		
		
		
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

		tipMatch1 = pathPoints1[indices[0]]

		" discrepancy distance between tip closest point and average departure point "
		dist1 = sqrt((tipMatch1[0]-pathPoints1[max1][0])**2 + (tipMatch1[1] - pathPoints1[max1][1])**2)

		tipMatch2 = pathPoints1[indices[-1]]

		" discrepancy distance between tip closest point and average departure point "
		dist2 = sqrt((tipMatch2[0]-pathPoints1[max2][0])**2 + (tipMatch2[1] - pathPoints1[max2][1])**2)


		DEP_THRESH = 0.2

		if maxFront > DEP_THRESH:
			departurePoint1 = pathPoints2[newFrontDepI]
			isExist1 = True

			if max1 == 0 or max1 == len(pathPoints1)-1:
				isInterior1 = False
			else:
				isInterior1 = True

		if maxBack > DEP_THRESH:
			departurePoint2 = pathPoints2[newBackDepI]
			isExist2 = True

			if max2 == 0 or max2 == len(pathPoints1)-1:
				isInterior2 = False
			else:
				isInterior2 = True

		" sum of closest points on front and back "
		" select the one with minimal cost "
		
		if False:
			pylab.clf()
			xP = range(len(pathPoints2))
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
			pylab.title("%1.2f %1.2f %d %d %d" % ( maxFront, maxBack, len(pathPoints1), max1, max2))
			pylab.savefig("distances_%04u.png" % self.pathPlotCount)

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in pathPoints2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')   

			if True:	
				xP = [pathPoints1[max1][0]]
				yP = [pathPoints1[max1][1]]
				pylab.scatter(xP,yP, color='b')		   
	
				xP = [tipMatch1[0]]
				yP = [tipMatch1[1]]
				pylab.scatter(xP,yP, color='r')		   

				xP = [pathPoints2[frontDepI][0]]
				yP = [pathPoints2[frontDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints1[forePathIndex][0]]
				yP = [pathPoints1[forePathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [pathPoints2[newFrontDepI][0]]
				yP = [pathPoints2[newFrontDepI][1]]
				pylab.scatter(xP,yP, color='m')		   
	
	
			if True:
				xP = [pathPoints1[max2][0]]
				yP = [pathPoints1[max2][1]]
				pylab.scatter(xP,yP, color='b')		   
	
				xP = [tipMatch2[0]]
				yP = [tipMatch2[1]]
				pylab.scatter(xP,yP, color='r')		   

				xP = [pathPoints2[backDepI][0]]
				yP = [pathPoints2[backDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints1[backPathIndex][0]]
				yP = [pathPoints1[backPathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [pathPoints2[newBackDepI][0]]
				yP = [pathPoints2[newBackDepI][1]]
				pylab.scatter(xP,yP, color='m')		   
	
	
			xP = []
			yP = []
			for p in pathPoints1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], (%d,%d)" % (maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, pathID1, pathID2))
			pylab.savefig("pathDeparture_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1
		
		
		return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2

	@logFunction
	def getDeparturePoint(self, currPath, nodeID, plotIter = False):
		
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
		
		
		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.getNodePose(nodeID)
		
		"Assumption:  one section of the medial axis is closely aligned with the path "
		poseOrigin = Pose(estPose2)
		
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		#points2 = medialSpline2.getUniformSamples(interpAngle=True)
		points2 = medialSpline2.getUniformSamples()

		points2_offset = []
		for p in points2:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points2_offset.append(result)


		pathPoints = orientPath(currPath, points2_offset)

		angle1, angle2 = getTipAngles(points2_offset)


		print "ang1:", angle1
		print "ang2:", angle2
		print "diff:", diffAngle(angle1, angle2)
		
		distSum = 0.0
		contigCount = 0
		maxContig = 0
		distances = []
		indices = []
		for i in range(0,len(points2_offset)):
			p_2 = points2_offset[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints, p_2)
			distances.append(minDist)
			indices.append(i_1)
			distSum += minDist
			
			if minDist < 0.2:
				contigCount += 1
				if contigCount > maxContig:
					maxContig = contigCount
			else:
				contigCount = 0
		
		overlapSum = distSum / float(len(points2_offset))
		
		print "maxContig,overlapSum:", maxContig, overlapSum
		
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


		frontAngleRefI = frontDepI
		while distances[frontAngleRefI] < 0.1:
			frontAngleRefI -= 1
			
			if frontAngleRefI < 0:
				frontAngleRefI = 0
				break
		
		forePathIndex = indices[frontAngleRefI]
		forePathAngle = pathPoints[forePathIndex][2]
		forePathAngle = normalizeAngle(forePathAngle + pi)

		foreDiffAngle = diffAngle(angle1, forePathAngle)
		
		newFrontDepI = 4
		while fabs(diffAngle(normalizeAngle(points2_offset[newFrontDepI][2]+pi), forePathAngle)) > pi/6.0:
			
			if newFrontDepI < frontDepI:
				newFrontDepI += 1
			else:
				break		 
		
		print "front:", nodeID, frontDepI, distances[frontDepI], frontAngleRefI, forePathIndex, forePathAngle, angle1, foreDiffAngle, newFrontDepI, diffAngle(normalizeAngle(points2_offset[newFrontDepI][2]+pi), forePathAngle) 
		

		" FIXME:  index out of bounds case "
		currI = 2
		try:
			while distances[-currI-3] < maxBack:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass
			""" raw sensor data products """

		backDepI = len(distances) - currI
		backPoint = [backDepI, distances[backDepI]]

		backAngleRefI = backDepI
		while distances[backAngleRefI] < 0.1:
			backAngleRefI += 1
			
			if backAngleRefI >= len(distances):
				backAngleRefI = len(distances)-1
				break
		
		backPathIndex = indices[backAngleRefI]
		backPathAngle = pathPoints[backPathIndex][2]
		
		backDiffAngle = diffAngle(angle2, backPathAngle)
		
		newBackDepI = len(distances)-5
		while fabs(diffAngle(points2_offset[newBackDepI][2], backPathAngle)) > pi/6.0:
			
			if newBackDepI > backDepI:
				newBackDepI -= 1
			else:
				break		 

		print "back:", nodeID, backDepI, distances[backDepI], backAngleRefI, backPathIndex, backPathAngle, angle2, backDiffAngle, newBackDepI, diffAngle(points2_offset[newBackDepI][2], backPathAngle) 


		foreAngDiffs = []
		for i in range(0,len(points2_offset)):
			foreAngDiffs.append(fabs(diffAngle(normalizeAngle(points2_offset[i][2] + pi), forePathAngle)))

		backAngDiffs = []
		for i in range(0,len(points2_offset)):
			backAngDiffs.append(fabs(diffAngle(points2_offset[i][2], backPathAngle)))



		"reset to the maximum distances "
		maxFront = distances[0]
		maxBack = distances[-1]



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

		tipMatch1 = pathPoints[indices[0]]

		" discrepancy distance between tip closest point and average departure point "
		dist1 = sqrt((tipMatch1[0]-pathPoints[max1][0])**2 + (tipMatch1[1] - pathPoints[max1][1])**2)

		tipMatch2 = pathPoints[indices[-1]]

		" discrepancy distance between tip closest point and average departure point "
		dist2 = sqrt((tipMatch2[0]-pathPoints[max2][0])**2 + (tipMatch2[1] - pathPoints[max2][1])**2)


		DEP_THRESH = 0.3

		if maxFront > DEP_THRESH:
			departurePoint1 = points2_offset[newFrontDepI]
			isExist1 = True

			if max1 == 0 or max1 == len(pathPoints)-1:
				isInterior1 = False
			else:
				isInterior1 = True

		if maxBack > DEP_THRESH:
			departurePoint2 = points2_offset[newBackDepI]
			isExist2 = True
			
			if max2 == 0 or max2 == len(pathPoints)-1:
				isInterior2 = False
			else:
				isInterior2 = True

		" sum of closest points on front and back "
		" select the one with minimal cost "
		
		if False:
			pylab.clf()
			xP = range(len(points2_offset))
			yP = distances
			pylab.plot(xP,yP, color ='b')
			
			yP = foreAngDiffs
			pylab.plot(xP,yP, color ='r')

			yP = backAngDiffs
			pylab.plot(xP,yP, color ='g')
			
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
			pylab.title("nodeID %d: %1.2f %1.2f %d %d %d %d %1.2f" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2, maxContig, overlapSum))
			pylab.savefig("distances_%04u.png" % self.pathPlotCount)

		if plotIter:
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

				xP = [points2_offset[frontDepI][0]]
				yP = [points2_offset[frontDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints[forePathIndex][0]]
				yP = [pathPoints[forePathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [points2_offset[newFrontDepI][0]]
				yP = [points2_offset[newFrontDepI][1]]
				pylab.scatter(xP,yP, color='m')		   
	
	
			if True:
				xP = [pathPoints[max2][0]]
				yP = [pathPoints[max2][1]]
				pylab.scatter(xP,yP, color='b')		   
	
				xP = [tipMatch2[0]]
				yP = [tipMatch2[1]]
				pylab.scatter(xP,yP, color='r')		   

				xP = [points2_offset[backDepI][0]]
				yP = [points2_offset[backDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints[backPathIndex][0]]
				yP = [pathPoints[backPathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [points2_offset[newBackDepI][0]]
				yP = [points2_offset[newBackDepI][1]]
				pylab.scatter(xP,yP, color='m')		   

			xP = [points2_offset[-1][0]]
			yP = [points2_offset[-1][1]]
			pylab.scatter(xP,yP, color=(0.5,0.5,0.5))		   

	
			xP = []
			yP = []
			for p in currPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2))
			pylab.savefig("departure_%04u_%04u.png" % (self.hypothesisID, self.pathPlotCount))
			
			self.pathPlotCount += 1
		
		print "departure %d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, frontDepI, backDepI)
		
		
		" if the medial axis does not overlap the path contiguously enough, mark a high discrepancy "
		
		maxContig, overlapSum
		contigFrac = float(maxContig)/float(len(points2_offset))
		print "returning:", contigFrac, overlapSum
		return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum
		

	@logFunction
	def getPathOverlapSum(self, path1, path2, pathID1, pathID2, plotIter = False):


		if len(path1) == 0:
			return 0.0
	
		medial2 = path2
			   
		minMatchDist2 = 0.1
		angThresh = math.pi/8.0
					
	
		supportSpline = SplineFit(path1, smooth=0.1)		
		vecPoints1 = supportSpline.getUniformSamples()
		supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			 
		
		support_pairs = []
		for i in range(len(vecPoints2)):
			
			p_2 = vecPoints2[i]
	
			try:
				p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, angThresh)
	
				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = supportPoints[minI][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					support_pairs.append([points2[i],supportPoints[minI],C2,C1])
			except:
				pass

		cost = 0.0
		if len(support_pairs) == 0:
			cost = 0.0
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
				
			#cost = sum1 / len(support_pairs)
			cost = sum1
			

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
				pylab.plot(xP,yP, color='k', alpha=0.3)
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d %d %d cost = %f, count = %d" % (self.hypothesisID, pathID1, pathID2, cost, len(support_pairs)))
			pylab.savefig("pathOverlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 0.0

		return cost


	@logFunction
	def checkUniqueBranch(self, parentPathID, nodeID1, depAngle, depPoint):


		print "checkUniqueBranch(", parentPathID, ",", nodeID1, ",", depAngle, ",", depPoint, ")"

		#foreTerm1 = frontInterior1 and frontExist1
		#foreTerm2 = frontInterior2 and frontExist2

		if depPoint == 0:
			print "REJECT, no departure point!", nodeID1, depPoint
			return False, -1

		" check cartesian distance to similar junction points from parent path "

		" BEING MORE PERMISSIVE IN CREATING NEW BRANCHES BECAUSE WE CAN MERGE LATER "

		" cartesian distance "
		DISC_THRESH = 0.5

		JOINT_DIST = 3.0

		" 60 degree threshold "
		#ANG_THRESH = 0.523 # pi/6
		ANG_THRESH = pi/8
		
		" maximum distance to the parent terminal point before creating a junction point "
		#TERM_THRESH = 0.8
		TERM_THRESH = 0.5
		
		pathIDs = self.getPathIDs()


		""" enforce too many neighbors constraint """
		neighborCount = 0

		for pathID in pathIDs:
			if pathID != 0:
				junctionPoint = self.getGlobalJunctionPose(pathID)
				neighDist = sqrt((depPoint[0]-junctionPoint[0])**2 + (depPoint[1]-junctionPoint[1])**2)

				if neighDist < JOINT_DIST:
					neighborCount += 1

		if neighborCount >= 2:
			print "REJECT, too many neighbors", nodeID1, pathID, neighborCount, depPoint
			return False, parentPathID


		""" check over all spatial landmarks to determine if we are close to a junction """

		controlPoses_L = {}
		for pathID in pathIDs:
			controlPoses_L[pathID] = self.pathClasses[pathID]["controlPose"]

		parentPathIDs = self.getParentHash()
		controlPoses_G = computeGlobalControlPoses(controlPoses_L, parentPathIDs)

		minLandmark = None
		minLandDist = 1e100

		for pathID in pathIDs:
			currFrame = Pose(controlPoses_G[pathID])
			for point_L, pointThresh, pointName in self.localLandmarks[pathID]:
				point_G = currFrame.convertLocalToGlobal(point_L)
				landmarkDist = sqrt((depPoint[0]-point_G[0])**2 + (depPoint[1]-point_G[1])**2)

				if landmarkDist < minLandDist:
					minLandDist = landmarkDist
					minLandmark = point_G

		#if minLandDist >= 0.5:
		if minLandDist >= 1.0:
			print "REJECT, proposed branch point is not close enough to a spatial landmark", nodeID1, pathID, minLandDist, minLandmark, neighborCount, depPoint
			return False, parentPathID

		#sF0 = self.poseData.spatialFeatures[nodeID1][0]
		#if sF0["bloomPoint"] == None and sF0["archPoint"] == None and sF0["inflectionPoint"] == None:
		#	print "REJECT, diverging pose has no spatial landmark feature", nodeID1, pathID, neighborCount, depPoint
		#	return False, parentPathID

		
		for pathID in pathIDs:

			if pathID == parentPathID:
				termPoint = self.terminals[pathID+1]
				
				print "termPoint =", termPoint
				termDist = sqrt((depPoint[0]-termPoint[1][0])**2 + (depPoint[1]-termPoint[1][1])**2)
				
				print "termDist =", termDist
				
				if termDist < TERM_THRESH:
					print "REJECT, proposed junction point is too close to parent pathID terminal", pathID, termDist, termPoint, depPoint
					return False, pathID

			
			path = self.getPath(pathID)
			
			if path["parentID"] == parentPathID:

				
				junctionNodeID = path["branchNodeID"]
				
				junctionPoint = self.getGlobalJunctionPose(pathID)
				print "existing branch pose:", junctionPoint
					
				dist = sqrt((depPoint[0]-junctionPoint[0])**2 + (depPoint[1]-junctionPoint[1])**2 )

				" check difference of tangential angle "
				angDiff = diffAngle(depAngle, junctionPoint[2])

				print "result = ", junctionNodeID, dist, angDiff
		
				#isNodeFeatureless = self.nodeHash[nodeID1].getIsFeatureless()
				isNodeFeatureless = self.poseData.isNodeFeatureless[nodeID1]
				print "node", nodeID1, "is featureless =", isNodeFeatureless
		
				if dist < DISC_THRESH:
					if fabs(angDiff) < ANG_THRESH:
						
						print "DUPLICATE junction of", junctionPoint, "rejecting with differences of", dist, angDiff
						return False, pathID
				
				" Explicit limitation that we wont be able to detect a T-Junction coming from left or right"
				#if isNodeFeatureless:
				#	print "REJECT new branch because the branching node is featureless"
				#	return False, -1
		
		return True, -1



	@logFunction
	def determineBranchSingle(self, nodeID1, frontExist1, frontInterior1, depAngle1, depPoint1, parentPathID1, dirFlag):
	   
		"""
		3) one node is departing only
			- is only one departure ( other node has no departure, or has external departure but not same angle )
			- is a false departure	( falseness happens if bad map data, ignore )
		
		Now, we have the departing angles, but we need to compare them with each other. 
		
		dirFlag specifies which node we want to check branches on.	We do not want to branch with the anchored end of a probe sweep   
		"""

		foreTerm1 = frontInterior1 and frontExist1

		pathBranchID = -1
		isBranch = False
		isNew = False
		
		if foreTerm1:

			isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				
			if isUnique1:
				   
				" foreTerm1 has unique departure "	  
				pathID = parentPathID1
				branchNodeID = nodeID1
				globalJunctionPoint = depPoint1
				depAng = depAngle1
				junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

				poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
				newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
				pathBranchID = newPathID
				isBranch = True
				isNew = True

					
			if duplicatePathID1 != -1:
				pathBranchID = duplicatePathID1
				isBranch = True

		else:
			" no new departure, just add nodes to the leaf paths "
			pathID = parentPathID1
			branchNodeID = nodeID1
			globalJunctionPoint = depPoint1
			
		
		return isBranch, pathBranchID, isNew
		

	@logFunction
	def determineBranchPair(self, nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag, isUnique1, isUnique2, duplicatePathID1, duplicatePathID2, shootIDs):
		
		"""
		1)	both nodes departing
			- are same departure  ( equal departing angle and are on top of each other)  1 new path
			- are different departures ( different departing angle and are not necessarily on top of each other) 2 new paths
			
		2) neither node is departing
			- are on same path
			
		3) one node is departing only
			- is only one departure ( other node has no departure, or has external departure but not same angle )
			- is a false departure	( falseness happens if bad map data, ignore )
			- are both departing ( should have at least external departure on both, compare departure angle )
		
		Now, we have the departing angles, but we need to compare them with each other. 
		
		dirFlag specifies which node we want to check branches on.	We do not want to branch with the anchored end of a probe sweep   
		"""

		print self.hypothesisID, "determineBranchPair(", nodeID1, ",", nodeID2, ",", frontExist1, ",", frontExist2, ",", frontInterior1, ",", frontInterior2, ",", depAngle1, ",", depAngle2, ",", depPoint1, ",", depPoint2, ",", parentPathID1, ",", parentPathID2, ",", dirFlag, ",", isUnique1, ",", isUnique2, ",", duplicatePathID1, ",", duplicatePathID2, ",", shootIDs, ")"
		
		self.shootIDs = shootIDs


		foreTerm1 = frontInterior1 and frontExist1
		foreTerm2 = frontInterior2 and frontExist2

		frontAngDiff = diffAngle(depAngle1, depAngle2)

		pathBranchIDs = [-1, -1]
		isBranch = [False, False]
		isNew = [False, False]
		
		" 60 degree threshold "
		ANG_THRESH = 1.047


		if foreTerm1 and foreTerm2:
			" both departing case: check if same or different "
			if fabs(frontAngDiff) < ANG_THRESH:
				" are on top of each other and are branching to the same path "

				#isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				#isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

				if isUnique1 and isUnique2:
					
					pathID = parentPathID1
					branchNodeID = nodeID1
					globalJunctionPoint = depPoint1
					depAng = depAngle1
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					
					newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID

					isBranch[0] = True
					isBranch[1] = True
					isNew[0] = True
					isNew[1] = True
					
					print "foreA"

				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True
						isNew[0] = True

				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True
						isNew[1] = True
					
				
				if duplicatePathID1 != -1:
					pathBranchIDs[0] = duplicatePathID1
					isBranch[0] = True
					

				if duplicatePathID2 != -1:
					pathBranchIDs[1] = duplicatePathID2
					isBranch[1] = True
					
				
			else:
				" these are branching to separate new paths "

				#isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				#isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

				if isUnique1 and isUnique2:

					if dirFlag == 0:
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID1 = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID1
						isBranch[0] = True
						isNew[0] = True


					if dirFlag == 1:
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID2 = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID2
						isBranch[1] = True
						isNew[1] = True

					print "foreB"

				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True
						isNew[0] = True

				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True
						isNew[1] = True

				if duplicatePathID1 != -1:
					pathBranchIDs[0] = duplicatePathID1
					isBranch[0] = True
					
				if duplicatePathID2 != -1:
					pathBranchIDs[1] = duplicatePathID2
					isBranch[1] = True
					
		
		elif foreTerm1 and not foreTerm2:

			#isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				
			if isUnique1:
				   
				if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
					" both have same departure "
					pathID = parentPathID1
					branchNodeID = nodeID1
					globalJunctionPoint = depPoint1
					depAng = depAngle1
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID
					isBranch[0] = True
					isBranch[1] = True
					isNew[0] = True
					isNew[1] = True

					print "foreC"
					
				else:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True
						isNew[0] = True

					print "foreD"
					
			if duplicatePathID1 != -1:
				pathBranchIDs[0] = duplicatePathID1
				isBranch[0] = True

		elif foreTerm2 and not foreTerm1:

			#isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

			if isUnique2:
				if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

					" both have same departure "
					pathID = parentPathID2
					branchNodeID = nodeID2
					globalJunctionPoint = depPoint2
					depAng = depAngle2
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID
					isBranch[0] = True
					isBranch[1] = True
					isNew[0] = True
					isNew[1] = True

					print "foreE"

				else:
					
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True
						isNew[1] = True

					print "foreF"

			if duplicatePathID2 != -1:
				pathBranchIDs[1] = duplicatePathID2
				isBranch[1] = True

		"""
		" CHECKING FOR BAD PATH MEMBERSHIP HERE.  May reuse later "
		else:
			print "foreG"
			" no new departure, just add nodes to the leaf paths "
			pathID = parentPathID1
			branchNodeID = nodeID1
			globalJunctionPoint = depPoint1
			
			if parentPathID1 != parentPathID2:
				print "departing path IDs of two colocated nodes are not the same"
				raise
		"""	
		
		return isBranch, pathBranchIDs, isNew, self.shootIDs
		
		
	@logFunction
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
			parent = self.getPath(currPathID)
			thisParentID = parent["parentID"]
			
			if thisParentID == None:
				break
			
			if pathIDs.count(thisParentID) == 0:
				break
			
			currPathID = thisParentID

		" currPathID is the root path "
		rootPathID = currPathID
		
			
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


			for join in self.joins:
				if join[0][0] == childPathID:
					childJunctionPoint = join[2]
					childI = join[0][1]

			if childI < abs(childI-len(childPath)-1):
				childTermI = len(childPath)-1
			else:
				childTermI = 0

			"find path between: (childPathID, childTermI) to (rootPathID, 0) and (rootPathID, len(rootPath)-1)"
			shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path((childPathID, childTermI))
			startNode1 = shortestPathSpanTree[(rootPathID, 0)]
			startNode2 = shortestPathSpanTree[(rootPathID, len(rootPath)-1)]
			
			currNode = startNode1
			splicedPath1 = []
			while currNode != (childPathID, childTermI):
				splicedPath1.append(self.pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath1.append(self.pathGraph.get_node_attributes(currNode))

			currNode = startNode2
			splicedPath2 = []
			while currNode != (childPathID, childTermI):
				splicedPath2.append(self.pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath2.append(self.pathGraph.get_node_attributes(currNode))

			newPaths.append(splicedPath1)
			newPaths.append(splicedPath2)
			
		else:
			"find entire path from rightPathID to leftPathID "
			" start from point farthest from child's junction point "
			
			rightPath = self.trimmedPaths[rightPathID]
			leftPath = self.trimmedPaths[leftPathID]
			
			" find the junction of this child's to its parent "    
			for join in self.joins:
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
			shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path((rightPathID, rightTermI))
			currNode = shortestPathSpanTree[(leftPathID, leftTermI)]

			splicedPath = []
			while currNode != (rightPathID, rightTermI):
				splicedPath.append(self.pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath.append(self.pathGraph.get_node_attributes(currNode))

			newPaths.append(splicedPath)





		return newPaths

	@logFunction
	def getNearestPathPoint(self, originPoint):

		" find closest point on each path "
		minI1 = 0
		minJ1 = 0
		minDist1 = 1e100
		
		pathIDs = self.getPathIDs()    
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for j in range(len(path)):
				p0 = path[j]
				dist = sqrt((originPoint[0]-p0[0])**2 + (originPoint[1]-p0[1])**2)
				
				if dist < minDist1:
					minDist1 = dist
					minI1 = pathID
					minJ1 = j
		
		newPoint = self.trimmedPaths[minI1][minJ1]
		
		return newPoint
		
	@logFunction
	def computeNavigationPath(self, startPose, endPose):
		
		" get the trimmed paths "
		
		" find closest point on each path "
		minI1 = 0
		minJ1 = 0
		minDist1 = 1e100
		minI2 = 0
		minJ2 = 0
		minDist2 = 1e100
		
		pathIDs = self.getPathIDs()    
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for j in range(len(path)):
				p0 = path[j]
				dist = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
				
				if dist < minDist1:
					minDist1 = dist
					minI1 = pathID
					minJ1 = j

				dist = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

				if dist < minDist2:
					minDist2 = dist
					minI2 = pathID
					minJ2 = j
		
		" select shortest one "
		
		" get origin and target path "
		
		" get splice of origin and target path "
		" splice minJ1 and minJ2 "
		orderPathIDs = self.getPathPath(minI1, minI2)		 
		
		#splicedPaths = self.splicePathIDs(orderPathIDs)
		
		
		allSplices, terminals, junctions = self.getAllSplices()

		splicedPaths = []
		for k, result in allSplices.iteritems():
			print "allSplices key:", k
			for sPath in result:
				#path = sPath['path']
				path = sPath['skelPath']
				splicedPaths.append(path)
	
		print len(splicedPaths), "spliced paths for path finding "
		
		" again, find closest point of start and end pose "
		
		minTotalDist = 3e100
		minSplicePathID = 0
		spliceIndex = (0,1)
		for i in range(len(splicedPaths)):

			minDist1 = 1e100
			minDist2 = 1e100

			minJ1 = 0
			minJ2 = 0

			
			path = splicedPaths[i]

			for j in range(len(path)):
				p0 = path[j]
				dist1 = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
				
				if dist1 < minDist1:
					minDist1 = dist1
					minJ1 = j

				dist2 = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

				if dist2 < minDist2:
					minDist2 = dist2
					minJ2 = j

			
			totalDist = minDist1 + minDist2
			if totalDist < minTotalDist:
				minTotalDist = totalDist
				minSplicePathID = i
				spliceIndex = (minJ1,minJ2)

			print i, path[0], path[-1], minJ1, minJ2, minDist1, minDist2, totalDist, minTotalDist, minSplicePathID, spliceIndex
		
		
		" return splicePath[startIndex:endIndex+1]"
		if spliceIndex[0] < spliceIndex[1]:
			return splicedPaths[minSplicePathID]
			#pathSpline = SplineFit(splicedPaths[minSplicePathID])
			#newPath = pathSpline.getUniformSamples()
			#return newPath
		else:
			path = splicedPaths[minSplicePathID]
			path.reverse()
			return path
			#pathSpline = SplineFit(path)
			#newPath = pathSpline.getUniformSamples()
			#return newPath


			

