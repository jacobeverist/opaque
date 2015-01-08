

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
from ParticleFilter import multiParticleFitSplice, batchLocalizeParticle, batchDisplaceParticles2, Particle, batchLocalizeLandmark
import time
import traceback
from uuid import uuid4
import matplotlib.pyplot as plt

import alphamod
from itertools import product

from shoots import computeShootSkeleton, spliceSkeletons, computeGlobalControlPoses, batchJointBranch, batchBranch, trimBranch, getBranchPoint, ensureEnoughPoints, getInitSkeletonBranchPoint, getSkeletonBranchPoint, getSkeletonPath, evaluateJointBranch, computeJointBranch
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
		#self.poseParticles = {}
		#self.poseParticles["numParticles"] = self.numPoseParticles
		#self.poseParticles["updateCount"] = 0
		#self.poseParticles["snapshots2"] = {0 : ()}

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
		self.branchDivergeCount = {0 : 0}
		self.branchSubsumeIDs = {0: []}
		self.branchTermDivergenceDist = {0: []}
		self.collisionLandmarks = {}

		""" intermediate data structure for computing the shoot skeleton """
		self.localLeaf2LeafPathJunctions = {}


		""" membership, parent-relationship, branch point, and control point for a shoot """
		self.pathClasses = {}
		self.pathClasses[0] = {"parentID" : None,
					"branchNodeID" : None,
					"localJunctionPose" : None, 
					"tipPoint_L" : None,
					"localDivergencePose" : None, 
					"sameProb" : {},
					"nodeSet" : [],
					"localNodePoses" : {},
					"globalJunctionPose" : None,
					"controlPose" : [0.0,0.0,0.0] }		

		""" number of shoots in the current map state """
		self.pathIDs = 1

		""" travel destinations on the shoot map """
		#self.pathTermsVisited = {-1: False, 0: False}
		self.pathTermsVisited = {}
		self.pathTermsData = {}
		self.pathTerms = {}

		""" number of branch hypotheses and spacing between them for a branch point distribution """
		#self.DIV_LEN = 0.2
		self.DIV_LEN = 0.1
		self.NUM_BRANCHES = 5
		#self.DIV_LEN = 0.1
		#self.NUM_BRANCHES = 1

		self.unfeaturedStepCount = 0


		""" the entrance point of the environment """
		self.rootPoint = [-3.2, 0.0]
		
		""" data structures for computing splices on the shoot map """
		self.allSplices = {}
		self.allSplices2 = {}
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

		print "parentPathIDs, controlPoses:", parentPathIDs, controlPoses


		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		#particleDist2 = self.poseParticles["snapshots2"][0]

		for pathID in allPathIDs:

			if nodeID in self.pathClasses[pathID]["nodeSet"]:

				if pathID != 0:

					""" convert the controlPose to coordinates local to the parent frame """
					shootControlPose = globalControlPoses[pathID]
					currFrame = Pose(shootControlPose)
					localNodePose = currFrame.convertGlobalPoseToLocal(newPose)

					self.pathClasses[pathID]["localNodePoses"][nodeID] = localNodePose

					#for part in particleDist2:
					#	part.junctionData[pathID]["localNodePoses"][nodeID] = localNodePose


				else:
					self.pathClasses[pathID]["localNodePoses"][nodeID] = newPose

					#for part in particleDist2:
					#	part.junctionData[pathID]["localNodePoses"][nodeID] = newPose


	@logFunction
	def updateMaxParticle2(self, maxIndex):

		print "updateMaxParticle2(", maxIndex, ")"


		maxParticle = self.stepResults[maxIndex]

		allPathIDs = self.getPathIDs()
		parentPathIDs = self.getParentHash()

		""" allow for the local node poses to change """
		#for pathID in allPathIDs:
			
		#	allNodes = self.pathClasses[pathID]["nodeSet"]
		#	for nodeID in allNodes:

		#		newPose = maxParticle.junctionData[pathID]["localNodePoses"][nodeID]
		#		self.pathClasses[pathID]["localNodePoses"][nodeID] = newPose


		""" change to the maximum likelihood branch position as well """

		""" maximum likelihood pose particle """
		partBranchTupleIndex = maxParticle["maxLikelihoodBranch"]
		branchArcDists = maxParticle["branchArcDists"]
		branchControls = maxParticle["branchControls"]

		branchPathIDs = deepcopy(allPathIDs)
		branchPathIDs.remove(0)
		branchPathIDs.sort()

		branchResult = None
		try:
			if partBranchTupleIndex != None:
				branchResult = self.jointBranchEvaluations[partBranchTupleIndex] 
		except:
			""" must have merged path, branch states are not updated yet """
			pass

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

				if branchResult != None:
					branchPose_L = branchResult["branchPoses_L"][pathID]
					tipPoint_L = branchResult["tipPoints_L"][pathID]
				else:
					branchPose_L = self.pathClasses[pathID]["localJunctionPose"]
					tipPoint_L = self.pathClasses[pathID]["tipPoint_L"]

				print "updating local tip point", self.hypothesisID, pathID, tipPoint_L

				""" FIXME:  was recycling the angle, but now we are regenerating it """
				self.setLocalJunctionPose(pathID, branchPose_L)
				self.setLocalTipPoint(pathID, tipPoint_L)


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

		termPoint = self.rootPoint

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

		#numParticles = self.poseParticles["numParticles"]
		#initDist2 = []

		" create the initial particle distibution "
		#while len(initDist2) < numParticles:
		#	hypDist = random.gauss(avgDist, 0.5)

		#	if hypDist > 0.0 and hypDist <= totalLen:
		#		hypPoint = pathSpline.getPointOfDist(hypDist)
		#		hypPose0 = copy(hypPoint)
		#		hypPose0[2] = newAng0
		#		hypPose1 = copy(hypPoint)
		#		hypPose1[2] = newAng1

		#		particleObj = Particle(hypPose0, hypPose1, pathID, hypDist, 0.0, self.hypothesisID)
		#		particleObj.spliceCurve = deepcopy(self.paths[0])
		#		particleObj.addNode(0,0, hypPose0)
		#		particleObj.addNode(1,0, hypPose1)
		#		initDist2.append(particleObj)

		#self.poseParticles["snapshots2"][0] = initDist2

		self.stepResults = []
		pathIDs = self.getPathIDs()
		parentPathIDs = self.getParentHash()

		nodeID0 = 0
		nodeID1 = 1

		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]

		medialSpline0 = SplineFit(medial0, smooth=0.1)
		medial0_vec = medialSpline0.getUniformSamples()
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medial1_vec = medialSpline1.getUniformSamples()

		allSplices = []
		sPaths0 = self.getAllSplices2()
		for sPath in sPaths0:
			allSplices.append(sPath['skelPath'])


		""" take all splices to initialize """
		matchedSplices = []
		for k in range(len(allSplices)):
			splice = allSplices[k]
			matchedSplices.append(splice)

		""" build sample locations along each splice, starting from the pose's origin """
		poseOrigin0 = Pose(pose0)
		globalMedial0 = []
		for p in medial0:
			globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

		self.stepDists = []
		for j in range(len(matchedSplices)):

			matchSplice = matchedSplices[j]


			orientedSplicePath = orientPath(matchSplice, globalMedial0)				
			pathSpline = SplineFit(orientedSplicePath, smooth=0.1)

			minDist0, oldU0, oldP0 = pathSpline.findClosestPoint(pose0)
			minDist1, oldU1, oldP1 = pathSpline.findClosestPoint(pose1)
			oldDist0 = pathSpline.dist_u(oldU0)
			oldDist1 = pathSpline.dist_u(oldU1)

			STEP_DIST = 0.2
			NUM_SAMPLES = 10


			""" sample points in the neighborhood """
			stepHigh = STEP_DIST * floor(NUM_SAMPLES/2.0)
			stepLow = -STEP_DIST * floor(NUM_SAMPLES/2.0)

			""" compute the sample points for the pose point distribution """
			for k in range(NUM_SAMPLES):
				newStepDist = stepLow + k * STEP_DIST

				newU0 = pathSpline.getUOfDist(oldU0, newStepDist)
				newU1 = pathSpline.getUOfDist(oldU1, newStepDist)

				newDist0 = pathSpline.dist_u(newU0)
				newPose0 = pathSpline.point_u(newU0)
				newAng0 = newPose0[2]

				newDist1 = pathSpline.dist_u(newU1)
				newPose1 = pathSpline.point_u(newU1)
				newAng1 = newPose1[2]

				self.stepDists.append((j, newDist0, newDist1, newPose0, newPose1))


		self.stepResults = []
		for result in self.stepDists:
			particleIndex = result[0]
			newPose0 = result[3]
			newPose1 = result[4]


			resultDict = {}
			resultDict["newPose0"] = newPose0
			resultDict["newPose1"] = newPose1
			resultDict["controlPoses_P"] = {}
			resultDict["displaceProb0"] = 0.0
			resultDict["displaceProb1"] = 0.0
			resultDict["branchPoses_G"] = {}
			resultDict["branchPoses_L"] = {}
			resultDict["landmarkSum"] = 0.0
			resultDict["maxLikelihoodBranch"] = None
			resultDict["branchArcDists"] = {}
			resultDict["branchControls"] = {}
			resultDict["motionBias"] = 0.0

			self.stepResults.append(resultDict)





		#for result in results:
		#	particleIndex = result[0]
		#	displaceProb0 = result[3]
		#	displaceProb1 = result[4]

		#	newPose0 = result[1]
		#	newPose1 = result[2]


		#	resultDict = {}
		#	resultDict["newPose0"] = newPose0
		#	resultDict["newPose1"] = newPose1
		#	resultDict["controlPoses_P"] = controlPoses
		#	resultDict["displaceProb0"] = displaceProb0
		#	resultDict["displaceProb1"] = displaceProb1
		#	resultDict["branchPoses_G"] = branchPoses_G
		#	resultDict["branchPoses_L"] = branchPoses_L
		#	resultDict["landmarkSum"] = 0.0
		#	resultDict["maxLikelihoodBranch"] = None
		#	resultDict["branchArcDists"] = ()
		#	resultDict["branchControls"] = ()

		#	self.stepResults.append(resultDict)


	@logFunction
	#def localizeLandmarkPose(self, pathID, nodeID0):
	def localizeLandmarkPose(self, landmarkNodeSets):
		#newNodePoses = self.localizeLandmarkPose(landmarkNodeSets)

		batchJobs = []

		maxProbs = {}
		maxSamples = {}

		for pathID, nodeSet in landmarkNodeSets.iteritems():
			maxProbs[pathID] = {}
			maxSamples[pathID] = {}

			longPaths_L = self.localLongPaths[pathID]
			allSplices = longPaths_L

			for nodeID in nodeSet:
				maxProbs[pathID][nodeID] = -1e100
				maxSamples[pathID][nodeID] = None



				medial0 = self.poseData.medialAxes[nodeID]

				medialSpline0 = SplineFit(medial0, smooth=0.1)
				medial0_vec = medialSpline0.getUniformSamples()

				oldNodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
				prevPose0 = oldNodePose_L


				""" find the splices that our poses match to so we can perform displacements """
				matchedSplices = []

				for k in range(len(allSplices)):

					splice = allSplices[k]

					results0 = getMultiDeparturePoint(splice, medial0_vec, prevPose0, prevPose0, [pathID,], nodeID, spliceIndex=k, plotIter=False)

					contigFrac0 = results0[12]

					print nodeID, "contigFrac0:", k, contigFrac0

					if contigFrac0 > 0.1:
						matchedSplices.append(splice)


				""" build sample locations along each splice, starting from the pose's origin """
				poseOrigin0 = Pose(prevPose0)
				globalMedial0 = []
				for p in medial0:
					globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

				""" compute the distribution spread based on how long we've gone unfeatured or unlandmarked """
				STEP_DIST = 0.1
				NUM_SAMPLES = 10 


				#print nodeID, "NUM_SAMPLES =", NUM_SAMPLES

				isFeatureless0 = self.poseData.isNodeFeatureless[nodeID]

				spatialFeature0 = self.poseData.spatialFeatures[nodeID][0]
				isSpatialFeature0 = spatialFeature0["bloomPoint"] != None or spatialFeature0["archPoint"] != None or spatialFeature0["inflectionPoint"] != None

				print nodeID, "feature state:", isFeatureless0, isSpatialFeature0

				sampleDists = []
				for j in range(len(matchedSplices)):

					matchSplice = matchedSplices[j]

					distEst = 0.0

					orientedSplicePath = orientPath(matchSplice, globalMedial0)				
					pathSpline = SplineFit(orientedSplicePath, smooth=0.1)

					minDist0, oldU0, oldP0 = pathSpline.findClosestPoint(prevPose0)
					dispU0 = pathSpline.getUOfDist(oldU0, distEst)

					""" sample points in the neighborhood """
					stepHigh = STEP_DIST * floor(NUM_SAMPLES/2.0)
					stepLow = -STEP_DIST * floor(NUM_SAMPLES/2.0)

					""" compute the sample points for the pose point distribution """
					for k in range(NUM_SAMPLES):
						newStepDist = stepLow + k * STEP_DIST

						newU0 = pathSpline.getUOfDist(dispU0, newStepDist)

						newDist0 = pathSpline.dist_u(newU0)
						newPose0 = pathSpline.point_u(newU0)

						sampleDists.append((j, newDist0, newPose0))


				""" collect landmarks that we can localize against """
				targetNodeLandmarks_N = {nodeID : None}
				targetNodeLandmarks_N[nodeID] = getNodeLandmark(nodeID, self.poseData)

				candLandmarks_L = nodeToLocalLandmarks(pathID, self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID,])

				#candLandmarks_L = deepcopy(self.localLandmarks[pathID])
				#del candLandmarks_L[nodeID]

				for stepJob in sampleDists:

					spliceIndex = stepJob[0]
					sampDist = stepJob[1]
					hypPose0 = stepJob[2]

					batchJobs.append([self.poseData, spliceIndex, pathID, nodeID, sampDist, hypPose0, [matchedSplices[spliceIndex],], candLandmarks_L, targetNodeLandmarks_N])



		print len(batchJobs), "total landmark jobs"

		time1 = time.time()
		results = batchLocalizeLandmark(batchJobs)
		time2 = time.time()
		print "TIME landmark batch", time2-time1

		sampleResults = []
		maxProb = -1e100
		maxSample = None
		for result in results:
			particleIndex = result[0]
			nodeID = result[1]
			pathID = result[2]
			newPose0 = result[3]
			displaceProb0 = result[4]


			resultDict = {}
			resultDict["newPose0"] = newPose0
			resultDict["displaceProb0"] = displaceProb0
			resultDict["pathID"] = pathID
			resultDict["nodeID"] = nodeID

			if displaceProb0 > maxProbs[pathID][nodeID] and displaceProb0 > 0.0:
				maxProbs[pathID][nodeID] = displaceProb0
				maxSamples[pathID][nodeID] = newPose0

			sampleResults.append(resultDict)

		return maxSamples


	@logFunction
	def batchDisplaceParticles(self, nodeID0, nodeID1):

		batchJobs = []
		pathIDs = self.getPathIDs()
		parentPathIDs = self.getParentHash()

		medial0 = self.poseData.medialAxes[nodeID0-2]
		medial1 = self.poseData.medialAxes[nodeID1-2]


		""" state of map:
		
		1) pose  = P
		2) fitted curve
		3) branch state = B

		utilVal = P( P | B ) * P(B)
		maxUtilVal = P( P | B_max ) * P(B_max)
		maxLikelihoodMap = P( P_max | B_max ) * P(B_max)


		"""


		medialSpline0 = SplineFit(medial0, smooth=0.1)
		medial0_vec = medialSpline0.getUniformSamples()
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medial1_vec = medialSpline1.getUniformSamples()

		prevPose0 = self.getNodePose(nodeID0-2)
		prevPose1 = self.getNodePose(nodeID1-2)

		allSplices = []
		sPaths0 = self.getAllSplices2()
		for sPath in sPaths0:
			allSplices.append(sPath['skelPath'])


		""" find the splices that our poses match to so we can perform displacements """
		matchedSplices = []

		#for splice in allSplices:
		for k in range(len(allSplices)):

			splice = allSplices[k]

			results0 = getMultiDeparturePoint(splice, medial0_vec, prevPose0, prevPose0, pathIDs, nodeID0, spliceIndex=k, plotIter=False)
			results1 = getMultiDeparturePoint(splice, medial1_vec, prevPose1, prevPose1, pathIDs, nodeID1, spliceIndex=k, plotIter=False)

			contigFrac0 = results0[12]
			contigFrac1 = results1[12]

			print nodeID0, nodeID1, "contigFrac0,contigFrac1:", k, contigFrac0, contigFrac1

			if contigFrac0 > 0.1 and contigFrac1 > 0.1:
				matchedSplices.append(splice)

		if len(matchedSplices) == 0:
			""" unable to match with any splice.  We have a problem! """
			raise


		""" build sample locations along each splice, starting from the pose's origin """
		poseOrigin0 = Pose(prevPose0)
		globalMedial0 = []
		for p in medial0:
			globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

		""" compute the distribution spread based on how long we've gone unfeatured or unlandmarked """
		#STEP_DIST = 0.2
		STEP_DIST = 0.4
		NUM_SAMPLES = 10 + self.unfeaturedStepCount * 2


		#print nodeID0, "NUM_SAMPLES =", NUM_SAMPLES

		isFeatureless0 = self.poseData.isNodeFeatureless[nodeID0]
		isFeatureless1 = self.poseData.isNodeFeatureless[nodeID1]

		spatialFeature0 = self.poseData.spatialFeatures[nodeID0][0]
		isSpatialFeature0 = spatialFeature0["bloomPoint"] != None or spatialFeature0["archPoint"] != None or spatialFeature0["inflectionPoint"] != None
		spatialFeature1 = self.poseData.spatialFeatures[nodeID1][0]
		isSpatialFeature1 = spatialFeature1["bloomPoint"] != None or spatialFeature1["archPoint"] != None or spatialFeature1["inflectionPoint"] != None

		print nodeID0, "feature state:", isFeatureless0, isFeatureless1, isSpatialFeature0, isSpatialFeature1

		if isFeatureless0 and isFeatureless1 and not isSpatialFeature0 and not isSpatialFeature1:
			self.unfeaturedStepCount += 1
		else:
			self.unfeaturedStepCount = 0


		""" forward direction , direction=True """
		direction = self.poseData.travelDirs[nodeID0]

		collisionState = {}
		collisionState[nodeID0] = [None, None]
		collisionState[nodeID1] = [None, None]

		try:
			self.collisionLandmarks
		except:
			self.collisionLandmarks = {}
		
		if direction:
			
			frontSum = 0.0
			frontProbeError = self.poseData.frontProbeError[nodeID0]
			for n in frontProbeError:
				frontSum += n
			foreAvg = frontSum / len(frontProbeError)

			print nodeID0, nodeID1, "fore collision avg:", foreAvg

			if foreAvg >= 1.1:
				collisionState[nodeID0][0] = True
				print nodeID0, nodeID1, "forward collision motion"

				self.collisionLandmarks[nodeID0] = medial0[0]

		else:

			backSum = 0.0
			backProbeError = self.poseData.backProbeError[nodeID1]
			for n in backProbeError:
				backSum += n
			backAvg = backSum / len(backProbeError)

			print nodeID0, nodeID1, "back collision avg:", backAvg

			if backAvg >= 1.1:
				collisionState[nodeID1][1] = True
				print nodeID0, nodeID1, "backward collision motion"

				self.collisionLandmarks[nodeID1] = medial1[-1]


		#travelDist0 = 0., travelDist1):
		#for matchSplice in matchedSplices:
		self.stepDists = []
		for j in range(len(matchedSplices)):

			matchSplice = matchedSplices[j]




			""" forward direction , direction=True """
			direction = self.poseData.travelDirs[nodeID0]
			
			if direction:

				if collisionState[nodeID0][0]:
					distEst = 0.0
				else:
					distEst = -0.8
			else:

				if collisionState[nodeID1][1]:
					distEst = 0.0
				else:
					distEst = 0.8

			print "distEst:", distEst

			#u2 = medialSpline2.getUOfDist(originU2, distEst, distIter = 0.001)

			orientedSplicePath = orientPath(matchSplice, globalMedial0)				
			pathSpline = SplineFit(orientedSplicePath, smooth=0.1)

			minDist0, oldU0, oldP0 = pathSpline.findClosestPoint(prevPose0)
			minDist1, oldU1, oldP1 = pathSpline.findClosestPoint(prevPose1)
			dispU0 = pathSpline.getUOfDist(oldU0, distEst)
			dispU1 = pathSpline.getUOfDist(oldU1, distEst)


			meanDist = pathSpline.dist_u(dispU0)

			#STEP_DIST = 0.2
			#NUM_SAMPLES = 10

			maxArcDist = pathSpline.dist_u(1.0)

			numSamples = int(pathSpline.dist_u(1.0)/STEP_DIST)

			print nodeID0, "numSamples =", numSamples, ", maxArcDist =", maxArcDist

			""" sample points in the neighborhood """
			stepHigh = STEP_DIST * floor(numSamples/2.0)
			stepLow = -STEP_DIST * floor(numSamples/2.0)

			SPARSE_DIST = 0.4
			DENSE_DIST = 0.4
			DENSE_ANG_THRESH = pi/6.0 # 30 degrees

			""" compute the sample points for the pose point distribution """
			#for k in range(numSamples):
			newStepDist = 0.0

			while newStepDist < maxArcDist:
				#newStepDist = stepLow + k * STEP_DIST

				#newU0 = pathSpline.getUOfDist(dispU0, newStepDist)
				#newU1 = pathSpline.getUOfDist(dispU1, newStepDist)
				newU0 = pathSpline.getUOfDist(0.0, newStepDist)
				newU1 = pathSpline.getUOfDist(0.0, newStepDist)

				newDist0 = pathSpline.dist_u(newU0)
				newPose0 = pathSpline.point_u(newU0)

				newDist1 = pathSpline.dist_u(newU1)
				newPose1 = pathSpline.point_u(newU1)


				#motionBias = 0.1 * gaussian(newDist0, meanDist, 0.5)
				motionBias = 0.25 * gaussian(newDist0, meanDist, 0.5)

				if abs(meanDist-newDist0) < 3.0:
					print "sample dist:", newStepDist, newDist0, newDist1, newU0, newU1, motionBias
					self.stepDists.append((j, newDist0, newDist1, newPose0, newPose1, newU0, newU1, motionBias))

				frontCurveU = pathSpline.getUOfDist(0.0, newStepDist + 2.5)
				backCurveU = pathSpline.getUOfDist(0.0, newStepDist - 2.5)
				frontPose = pathSpline.point_u(frontCurveU)
				backPose = pathSpline.point_u(backCurveU)

				if abs(newPose0[2]-frontPose[2]) > DENSE_ANG_THRESH or abs(newPose0[2]-backPose[2]) > DENSE_ANG_THRESH:
					newStepDist += DENSE_DIST
				else:
					newStepDist += SPARSE_DIST

		branchPoses_G = self.getGlobalBranchPoses()
		branchPoses_L = self.getLocalBranchPoses()

		controlPoses = deepcopy(self.getControlPoses())
		controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

		""" collect landmarks that we can localize against """
		candLandmarks_G = nodeToGlobalLandmarks(self.getControlPoses(), self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])

		targetNodeLandmarks_N = {nodeID0 : None, nodeID1 : None}
		targetNodeLandmarks_N[nodeID0] = getNodeLandmark(nodeID0, self.poseData)
		targetNodeLandmarks_N[nodeID1] = getNodeLandmark(nodeID1, self.poseData)


		batchJobs = []
		#for stepJob in self.stepDists:
		for jobIndex in range(len(self.stepDists)):

			stepJob = self.stepDists[jobIndex]

			#self.stepDists.append((j, newDist0, newDist1, newPose0, newPose1))
			spliceIndex = stepJob[0]
			hypPose0 = stepJob[3]
			hypPose1 = stepJob[4]

			batchJobs.append([self.poseData, jobIndex, spliceIndex, nodeID1, prevPose0, prevPose1, hypPose0, hypPose1, matchedSplices[spliceIndex], [matchedSplices[spliceIndex],], [], candLandmarks_G, targetNodeLandmarks_N])

		print len(batchJobs), "total displacement2 jobs"

		results = batchDisplaceParticles2(batchJobs)

		self.stepResults = []
		for result in results:
			jobIndex = result[0]
			spliceIndex = result[1]
			displaceProb0 = result[4]
			displaceProb1 = result[5]

			#return currPose2, currPose3, currProb2, currProb3, currAngDiff2, currAngDiff3, currContigFrac2, currContigFrac3
			
			arcDist0 = self.stepDists[jobIndex][1]
			arcDist1 = self.stepDists[jobIndex][2]
			spliceU0 = self.stepDists[jobIndex][5]
			spliceU1 = self.stepDists[jobIndex][6]
			motionBias = self.stepDists[jobIndex][7]

			#print "jobIndex:", jobIndex, arcDist0, arcDist1, spliceU0, spliceU1

			newPose0 = result[2]
			newPose1 = result[3]

			newAngDiff0 = result[6]
			newAngDiff1 = result[7]
			newContigFrac0 = result[8]
			newContigFrac1 = result[9]
			newOverlap0 = result[10]
			newOverlap1 = result[11]


			resultDict = {}
			resultDict["newPose0"] = newPose0
			resultDict["newPose1"] = newPose1
			resultDict["controlPoses_P"] = controlPoses
			resultDict["displaceProb0"] = displaceProb0
			resultDict["displaceProb1"] = displaceProb1
			resultDict["branchPoses_G"] = branchPoses_G
			resultDict["branchPoses_L"] = branchPoses_L
			resultDict["landmarkSum"] = 0.0
			resultDict["maxLikelihoodBranch"] = None
			resultDict["branchArcDists"] = {}
			resultDict["branchControls"] = {}

			resultDict["angDiff0"] = newAngDiff0
			resultDict["angDiff1"] = newAngDiff1
			resultDict["contigFrac0"] = newContigFrac0
			resultDict["contigFrac1"] = newContigFrac1
			resultDict["overlap0"] = newOverlap0
			resultDict["overlap1"] = newOverlap1
			resultDict["arcDist0"] = arcDist0
			resultDict["arcDist1"] = arcDist1
			resultDict["spliceU0"] = spliceU0
			resultDict["spliceU1"] = spliceU1
			resultDict["motionBias"] = motionBias

			resultDict["spliceIndex"] = spliceIndex

			self.stepResults.append(resultDict)


			#self.stepResults.append((newPose0, newPose1, displaceProb0, displaceProb1))

		#sortedResults = sorted(self.stepResults, key=itemgetter(2)) 

		""" evaluate the landmark consistency to select the max hypothesis """
		landmarkPoint0_N = targetNodeLandmarks_N[nodeID0]
		landmarkPoint1_N = targetNodeLandmarks_N[nodeID1]

		maxLandmarkSum = 0.0
		maxPartIndex = 0

		if landmarkPoint0_N != None or landmarkPoint1_N != None: 

			#for partIndex in range(len(particleDist2)):
			for partIndex in range(len(self.stepResults)):

				part = self.stepResults[partIndex]

				currLandmarks = deepcopy(candLandmarks_G)

				currPose0 = part["newPose0"]
				currPose1 = part["newPose1"]

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

				part["landmarkSum"] = currSum


				print "poseSum:", partIndex, currSum

				if currSum > maxLandmarkSum:
					maxLandmarkSum = currSum
					maxPartIndex = partIndex

		maxProb = -1e100
		maxPart2 = 0

		#for partIndex in range(len(particleDist2)):
		for partIndex in range(len(self.stepResults)):
			resultDict = self.stepResults[partIndex]
			currProb0 = resultDict["displaceProb0"]
			currProb1 = resultDict["displaceProb1"]
			angDiff0 = resultDict["angDiff0"]
			angDiff1 = resultDict["angDiff1"]
			contigFrac0 = resultDict["contigFrac0"]
			contigFrac1 = resultDict["contigFrac1"]
			overlap0 = resultDict["overlap0"]
			overlap1 = resultDict["overlap1"]
			arcDist0 = resultDict["arcDist0"]
			arcDist1 = resultDict["arcDist1"]
			spliceU0 = resultDict["spliceU0"]
			spliceU1 = resultDict["spliceU1"]
			poseLandmarkSum = resultDict["landmarkSum"]
			motionBias = resultDict["motionBias"]


			""" check to avoid divide by zero """
			newProb0 = 0.0
			if maxLandmarkSum > 0.0:
				if maxLandmarkSum > poseLandmarkSum:

					newProb0 = -2*overlap0 - 2*overlap1 - angDiff0/pi - angDiff1/pi + 2.0*(maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum


					#newProb0 = currProb0 * currProb1 * ((maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum)
					#newProb0 = currProb0 * currProb1 * ((maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum)
					#newProb0 = (maxLandmarkSum-poseLandmarkSum)/maxLandmarkSum
				else:
					""" maximum landmark cost, means utility is zeroed """
					newProb0 = 0.0
			else:
				#newProb0 = currProb0 * currProb1
				newProb0 = -overlap0 - overlap1 - angDiff0/pi - angDiff1/pi

			""" if the sample has been rejected, zero out evaluation """
			if contigFrac0 <= 0.0 or contigFrac1 <= 0.0:
				newProb0 = 0.0

			resultDict["evalProb"] = newProb0

	
			""" if the sample hasn't been rejected, compute evaluation """
			if contigFrac0 > 0.0 and contigFrac1 > 0.0:
				evalWithBias = newProb0 + motionBias

				""" find the maximum evaluation """
				if evalWithBias > maxProb:
					maxProb = evalWithBias
					maxPart2 = partIndex
			else:
				evalWithBias = 0.0

			resultDict["evalWithBias"] = evalWithBias

			print nodeID0, nodeID1, "displace sample:", partIndex, arcDist0, arcDist1, spliceU0, spliceU1, poseLandmarkSum, angDiff0, angDiff1, contigFrac0, contigFrac1, overlap0, overlap1, currProb0, currProb1, newProb0, motionBias, evalWithBias

			#if newProb0 > maxProb:
			#	maxProb = newProb0

		self.currMaxIndex = maxPart2

		""" change to the maximum likelihood branch position as well """
		self.setNodePose(nodeID0, deepcopy(self.stepResults[maxPart2]["newPose0"]))
		self.setNodePose(nodeID1, deepcopy(self.stepResults[maxPart2]["newPose1"]))
		self.updateMaxParticle2(self.currMaxIndex)


		#self.snapToParent()

	#@logFunction
	def localizePoseParticles2(self, nodeID0, nodeID1):

		#cProfile.run('runTest(probe)', 'test_prof')

		JOINT_DIST = 3.0

		poseData = self.poseData

		""" smoothed and vectored root path """


		staticSplicedPaths = []
		sPaths0 = self.getAllSplices2()
		for sPath in sPaths0:
			staticSplicedPaths.append(sPath['skelPath'])

		self.isNoLocalize = False
		self.resetBranches()

		#updateCount = self.poseParticles["updateCount"] 
		#particleDist2 = self.poseParticles["snapshots2"][0]

		#resultDict["controlPoses_P"] = controlPoses

		nodePose0 = self.getNodePose(nodeID0)

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



		if True:

			jointComboControlSets = []
			pathSplines = {}
			allPathIDs = self.getPathIDs()
			parentPathIDs = self.getParentHash()

			for pathID in allPathIDs:
				#pathSpline = SplineFit(self.localPaths[pathID])
				pathSpline = SplineFit(self.controlCurves[pathID])
				pathSplines[pathID] = pathSpline

			branchPathIDs = deepcopy(allPathIDs)
			branchPathIDs.remove(0)
			branchPathIDs.sort()

			#self.pathClasses[0] = {"parentID" : None,
			#			"branchNodeID" : None,
			#			"localJunctionPose" : None, 
			#			"tipPoint_L" : None,
			#			"localDivergencePose" : None, 
			#			"sameProb" : {},
			#			"nodeSet" : [],
			#			"localNodePoses" : {},
			#			"globalJunctionPose" : None,
			#			"controlPose" : [0.0,0.0,0.0] }		
			
			for particleIndex in range(len(self.stepResults)):
				part = self.stepResults[particleIndex]
				hypPose0 = part["newPose0"]


				origBranchDists = {}
				origBranchControls = {}
				origBranchArcDists = {}

				for pathID in branchPathIDs:

					parentID = parentPathIDs[pathID]
					pathSpline = pathSplines[pathID]

					""" this is the point where the control pose gets erroneously moved back onto the longest path """
					controlPose = part["controlPoses_P"][pathID]
					#print "oldControlPose:", particleIndex, pathID, controlPose
					minDist, controlUVal, newControlPose = pathSpline.findClosestPoint(controlPose)
					arcDist = pathSpline.dist_u(controlUVal)

					globJuncPose = part["branchPoses_G"][pathID]
					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)

					origBranchDists[pathID] = dist1
					origBranchControls[pathID] = newControlPose
					origBranchArcDists[pathID] = arcDist


				newBranchArcDists = {}
				newBranchControls = {}
				newBranchComboControls = {}

				print "origBranch:", origBranchArcDists, origBranchControls

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
						#print "newArcDist:", newArcDist
						#if origBranchDists[pathID] < JOINT_DIST:
						#if False:
						#if True:
						""" localize the shoots every 20 poses """
						if self.poseData.numNodes % 20 == 0:

							#if minIndex < 8:
							#	minIndex = 8

							#if minIndex > len(self.branchArcDists[pathID]) - 9:
							#	minIndex = len(self.branchArcDists[pathID]) - 9

							#arcHigh = arcDist + self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
							#arcLow = arcDist - self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
							if minIndex < 4:
								minIndex = 4

							if minIndex > len(self.branchArcDists[pathID]) - 5:
								minIndex = len(self.branchArcDists[pathID]) - 5

							arcList = []
							controlList = []
							comboList = []
							#for k in range(minIndex-8,minIndex+9):
							for k in range(minIndex-4,minIndex+5):
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

						#print "newBranchArcDists:", newBranchArcDists
						#print "newBranchControls:", newBranchControls

					argSet2 = []
					for pathID in branchPathIDs:
						argSet2.append(newBranchComboControls[pathID])

					combIterator = product(*argSet2)
					for comb in combIterator:
						jointComboControlSets.append(comb)

				""" set the range of branches for this particular particle so we can reference them later """
				part["branchArcDists"] = newBranchArcDists
				part["branchControls"] = newBranchControls

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

		print "self.nodeLandmarks =", self.nodeLandmarks

		for pathID in pathIDs:

			nodeSet = self.getNodes(pathID)
			for nodeID in nodeSet:

				print "pathID,nodeID", pathID, nodeID

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

		maxTuple = None
		maxNormCost = -1e10
		for arcTuple, branchResult in self.jointBranchEvaluations.iteritems():

			normLandmark = branchResult["normLandmark"]
			normMatchCount = branchResult["normMatchCount"]
			normCost = branchResult["normCost"] 

			#print "branch eval:", arcTuple, normLandmark, normMatchCount, normCost

			if normLandmark > maxNormCost:
				maxNormCost = normLandmark
				maxTuple = arcTuple

		maxSplicedPaths = []
		if maxTuple != None:

			branchResult = self.jointBranchEvaluations[maxTuple]

			#branchResult = self.jointBranchEvaluations[arcTuple]
			totalMatchCount = branchResult["totalMatchCount"]
			normMatchCount = branchResult["normMatchCount"]
			normLandmark = branchResult["normLandmark"]


			print "maxTuple:", maxTuple, normLandmark, maxNormCost, totalMatchCount, branchResult["landmarkCost"]

			#if normLandmark > maxNormCost:
			#	maxNormCost = normLandmark
			#	maxTuple = arcTuple

			landmarks_G = branchResult["landmarks_G"]

			controlPoses = branchResult["controlSet"]
			controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

			""" collect landmarks that we can localize against """
			candLandmarks_G = nodeToGlobalLandmarks(controlPoses, self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])


			splices_G = branchResult["splices_G"]
			for splice in splices_G:
				maxSplicedPaths.append((arcTuple, normLandmark, splice, None, candLandmarks_G))

		thisSplicedPaths = []
		if maxTuple != None:
			print "maxTuple:", maxTuple

			arcTuple = maxTuple
			# consider single path longestPaths , including root

			branchResult = self.jointBranchEvaluations[arcTuple]
			normLandmark = branchResult["normLandmark"]
			landmarks_G = branchResult["landmarks_G"]

			controlPoses = branchResult["controlSet"]
			controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

			""" collect landmarks that we can localize against """
			candLandmarks_G = nodeToGlobalLandmarks(controlPoses, self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])

			splices_G = branchResult["splices_G"]
			for splice in splices_G:
				thisSplicedPaths.append((arcTuple, normLandmark, splice, None, candLandmarks_G))

		else:
			#if len(thisSplicedPaths) == 0:
			""" static splices have 1.0 probability """

			controlPoses = deepcopy(self.getControlPoses())
			controlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

			""" collect landmarks that we can localize against """
			candLandmarks_G = nodeToGlobalLandmarks(controlPoses, self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses, exemptNodes = [nodeID0,nodeID1])

			for splice in staticSplicedPaths:
				thisSplicedPaths.append((None, 1.0, splice, None, candLandmarks_G))

		allSplicedPaths = thisSplicedPaths

		for particleIndex in range(len(self.stepResults)):

			time1 = time.time()

			part = self.stepResults[particleIndex]
			hypPose0 = part["newPose0"]
			hypPose1 = part["newPose1"]

			if nodeID0 > 0:
				prevHypPose0 = self.getNodePose(nodeID0-2)
				prevHypPose1 = self.getNodePose(nodeID1-2)
			else:
				prevHypPose0 = self.getNodePose(nodeID0)
				prevHypPose1 = self.getNodePose(nodeID1)

			resultsBySplice = []

			#thisSplicedPaths = []

			""" get control set of particle's branches indexed by arc distance """
			# part.branchArcDists = newBranchArcDists
			# part.branchControls = newBranchControls

			if maxTuple != None:
				for kIndex in range(len(branchPathIDs)):
					pathID = branchPathIDs[kIndex]
					maxDist = maxTuple[kIndex]

					""" compute the sample points for the branch point distribution """
					arcDists = [maxDist,]
					part["branchArcDists"][pathID] = arcDists
					print "arcDists:", pathID, arcDists

				part["maxLikelihoodBranch"] = maxTuple

			""" build indexing tuples for this particle """
			argSet = []
			for pathID in branchPathIDs:
				argSet.append(part["branchArcDists"][pathID])

			arcIndexes = []
			if len(argSet) > 0:
				combIterator = product(*argSet)
				for comb in combIterator:
					arcIndexes.append(tuple(comb))


			#print "particle:", particleIndex, ",",  len(thisSplicedPaths), "localize jobs"
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
				
				#localizeJobs.append([oldMedialP0, oldMedialU0, 0.0, oldMedialP1, oldMedialU1, 0.0, branchSampleIndex, spliceCount, orientedPath0, medial0_vec, medial1_vec, deepcopy(hypPose0), deepcopy(hypPose1), prevMedial0_vec, prevMedial1_vec, prevHypPose0, prevHypPose1, [], nodeID0, nodeID1, particleIndex, 0, self.hypothesisID, probVal, landmarks_G, landmark0_N, landmark1_N])
				localizeJobs.append([oldMedialP0, oldMedialU0, 0.0, oldMedialP1, oldMedialU1, 0.0, branchSampleIndex, spliceIndex, orientedPath0, medial0_vec, medial1_vec, deepcopy(hypPose0), deepcopy(hypPose1), prevMedial0_vec, prevMedial1_vec, prevHypPose0, prevHypPose1, [], nodeID0, nodeID1, particleIndex, 0, self.hypothesisID, probVal, landmarks_G, landmark0_N, landmark1_N])

				self.pathPlotCount2 += 1
				spliceCount += 1

			#allSplicedPaths += thisSplicedPaths


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
				
				branchNormLandmarkCost = 1.0
				if branchIndex != None:
					branchResult = self.jointBranchEvaluations[branchIndex]
					branchNormLandmarkCost = branchResult["normLandmark"]


				initPose0 = part[48]
				initPose1 = part[49]

				overlapSum = part[50]
				poseLandmarkSum = part[51]
				tipAngDiff0 = part[52]
				tipAngDiff1 = part[53]


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
				#ANG_THRESH = 1.0*pi/3.0
				#ANG_THRESH = 0.5*pi/3.0
				ANG_THRESH = pi/3.0


				#if fabs(diffAngle(initPose0[2],newPose0[2])) > 2.0*pi/3.0 or fabs(diffAngle(initPose1[2],newPose1[2])) > 2.0*pi/3.0:
				#if fabs(diffAngle(initPose0[2],newPose0[2])) > ANG_THRESH or fabs(diffAngle(initPose1[2],newPose1[2])) > ANG_THRESH:

				if tipAngDiff0 > ANG_THRESH or tipAngDiff1 > ANG_THRESH:
					print "reject because change in tip angle is too different"
					isReject = True
				else:
					pass
					#print "localize angle diff:", particleID, initPose0[2], newPose0[2], initPose1[2], newPose1[2], fabs(diffAngle(initPose0[2],newPose0[2])), fabs(diffAngle(initPose1[2],newPose1[2]))

				normLandmarkSum = 0.0
				normAngDiff = pi-fabs(diffAngle(initPose0[2],newPose0[2])) 


				if isReject:
					#newUtilVal = -1e100
					newUtilVal = 0.0
					listCopy = list(part)
					listCopy[45] = newUtilVal
					tupleCopy = tuple(listCopy)

					preBranchProb = 0.0

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
						#newProb = contigFrac_0
						newProb = contigFrac_0 * contigFrac_1

					#newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0 * contigFrac_0 / overlapSum_0

					""" case if splice is root path """
					#if normMatchCount != None:
					#	newProb *= normMatchCount

					preBranchProb = newProb

					""" only using one branch, so we comment this out for now """
					#newProb *= branchNormLandmarkCost


					""" add motion bias """
					try:
						motionBias = self.stepResults[particleID]["motionBias"]
					except:
						print self.stepResults[particleID]
						raise
					


					""" disable motion bias """
					#if preBranchProb > 0.0:
					#	newProb = preBranchProb + motionBias/3.0

					listCopy = list(part)
					listCopy[45] = newProb
					tupleCopy = tuple(listCopy)

					results[index] = tupleCopy

				print "%d %d %d %d %d isReject branchProb poseProb_B poseProb %d %1.2f %1.2f %1.2f" %  (nodeID0, self.hypothesisID, particleID, index, spliceIndex, isReject, branchNormLandmarkCost, results[index][45], preBranchProb), normMatchCount, overlapSum, contigFrac_0, contigFrac_1, initPose0[2], newPose0[2], int(isExist1_0), int(isExist2_0), int(isExist1_1), int(isExist2_1), int(isInterior1_0), int(isInterior2_0), int(isInterior1_1), int(isInterior2_1), branchIndex, poseLandmarkSum, len(landmarks_G), tipAngDiff0, tipAngDiff1
				
				#print "%d %d %d %d %d normMatchCount, newUtilVal, overlapSum:" % (nodeID0, self.hypothesisID, particleID, index, spliceIndex), normMatchCount, results[index][45], overlapSum, contigFrac_0, contigFrac_1, initPose0[2], newPose0[2], int(isExist1_0), int(isExist2_0), int(isExist1_1), int(isExist2_1), int(isInterior1_0), int(isInterior2_0), int(isInterior1_1), int(isInterior2_1), isReject, branchIndex, poseLandmarkSum, len(landmarks_G), normLandmarkSum, normAngDiff, branchNormLandmarkCost

			""" break condition, do not accept divergence if we have already branched 
				in reality, we should still select best diverging pose when branching because otherwise it is random	
			"""

			#if not rejectDivergence or self.isNodeBranching[nodeID0] or self.isNodeBranching[nodeID1]:
			if not rejectDivergence:
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

		print "particle evaluation:", nodeID0, self.hypothesisID
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

			contigFrac_1 = part[41]
			overlapSum_1 = part[42]

			newPose1 = part[24]

			poseLandmarkSum = part[51]
			tipAngDiff0 = part[52]
			tipAngDiff1 = part[53]
				
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


			particleDict = self.stepResults[particleID]

			particleDict["newPose0"] = newPose0
			particleDict["newPose1"] = newPose1
			particleDict["evalProb"] = newProb
			particleDict["spliceIndex"] = spliceIndex

			particleDict["overlap0"] = overlapSum_0
			particleDict["overlap1"] = overlapSum_1

			particleDict["contigFrac0"] = contigFrac_0
			particleDict["contigFrac1"] = contigFrac_1

			particleDict["angDiff0"] = tipAngDiff0
			particleDict["angDiff1"] = tipAngDiff1

			particleDict["landmarkSum"] = poseLandmarkSum

			""" Change branch if we are within a junction, otherwise, keep the most recent branch index """
			if branchTupleIndex != None:
				#particleObj.junctionData[pathID]["maxLikelihoodBranch"] = branchIndex
				#particleObj.maxLikelihoodBranch = branchTupleIndex
				particleDict["maxLikelihoodBranch"] = branchTupleIndex

				""" set branch points and control poses """
				branchResult = self.jointBranchEvaluations[branchTupleIndex] 
				controlSet = branchResult["controlSet"]
				branchPoses_L = branchResult["branchPoses_L"]

				controlPoses_G = computeGlobalControlPoses(controlSet, self.getParentHash())

				particleDict["controlPoses_P"] = controlSet

				for pathID in branchPathIDs:
					branchPose_L = branchPoses_L[pathID]
					controlPose_G = controlPoses_G[pathID]
					currFrame = Pose(controlPose_G)
					branchPose_G = currFrame.convertLocalOffsetToGlobal(branchPose_L)
					#particleObj.junctionData[pathID]["controlPose"] = controlSet[pathID]
					#particleObj.junctionData[pathID]["localJunctionPose"] = branchPose_L
					#particleObj.junctionData[pathID]["globalJunctionPose"] = branchPose_G
					particleDict["branchPoses_G"][pathID] = branchPose_G
					particleDict["branchPoses_L"][pathID] = branchPose_L



		""" update the stepResults structure to include arcDists """

		poseOrigin0 = Pose(self.getNodePose(nodeID0))
		globalMedial0 = []
		for p in medial0:
			globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

		poseOrigin1 = Pose(self.getNodePose(nodeID1))
		globalMedial1 = []
		for p in medial1:
			globalMedial1.append(poseOrigin1.convertLocalToGlobal(p))

		orientedSplices = []
		orientedPathSplines = []

		for spliceIndex in range(len(allSplicedPaths)):
			spliceCurve = allSplicedPaths[spliceIndex][2]

			orientedSplicePath = orientPath(spliceCurve, globalMedial0)				
			pathSpline = SplineFit(orientedSplicePath, smooth=0.1)

			orientedPathSplines.append(pathSpline)

		for k in range(len(self.stepResults)):
			part = self.stepResults[k]
			spliceIndex = part["spliceIndex"]
			poseProbVal = part["evalProb"]
			newPose0 = part["newPose0"]

			pathSpline = orientedPathSplines[spliceIndex]

			minDist0, newU0, newP0 = pathSpline.findClosestPoint(newPose0)
			newDist0 = pathSpline.dist_u(newU0)

			minDist1, newU1, newP1 = pathSpline.findClosestPoint(newPose1)
			newDist1 = pathSpline.dist_u(newU1)

			part["arcDist0"] = newDist0
			part["arcDist1"] = newDist1

			part["spliceU0"] = newU0
			part["spliceU1"] = newU1

			meanMinDist0, meanU0, meanP0 = pathSpline.findClosestPoint(self.getNodePose(nodeID0))
			meanDist0 = pathSpline.dist_u(meanU0)
			motionBias = 0.2 * gaussian(newDist0, meanDist0, 0.5)
			part["motionBias"] = motionBias
			if poseProbVal > 0.0:
				part["evalWithBias"] = poseProbVal + motionBias
			else:
				part["evalWithBias"] = 0.0




		#self.drawPoseParticles()

		""" now find the maximum pose """
		probSum = 0.0
		for part in self.stepResults:
			#probVal = part["evalProb"]
			probVal = part["evalWithBias"]
			probSum += probVal

		probParticles = []
		for k in range(len(self.stepResults)):
			part = self.stepResults[k]

			#poseProbVal = part["evalProb"]
			poseProbVal = part["evalWithBias"]

			""" if all 0.0 probabilities, make uniform distribution, set no localization flag """
			if probSum > 0.0:
				probParticles.append(poseProbVal/probSum)
			else:
				probParticles.append(float(1)/float(len(self.stepResults)))
				self.isNoLocalize = True

		print "probParticles:", self.isNoLocalize, probParticles 
		print "stepResults:", self.stepResults

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
			self.updateMaxParticle2(maxIndex)

			""" we update the canonical poses afterward because for now, we don't let node poses vary between particles """
			""" if we let them vary, then we would need to generate unique shoot maps for each particle """

			""" update the poses across all particles based on this max displacement result """
			#self.setNodePose(nodeID0, deepcopy(particleDist2[maxIndex].pose0))
			#self.setNodePose(nodeID1, deepcopy(particleDist2[maxIndex].pose1))
			self.setNodePose(nodeID0, deepcopy(self.stepResults[maxIndex]["newPose0"]))
			self.setNodePose(nodeID1, deepcopy(self.stepResults[maxIndex]["newPose1"]))

			""" using the results of localizePair() to position the nodes for generating map poses """
			self.generatePaths()


		else:
			print "no maximum particle!", probParticles
			raise

		""" FIXME:  do we really do nothing here when we tie for maximum? """
		#hypPose0 = samp["newPose0"]
		#hypPose1 = samp["newPose1"]
		#currProb0 = samp["displaceProb0"]
		#currProb1 = samp["displaceProb1"]
		#arcDist0 = samp["arcDist0"]
		#arcDist1 = samp["arcDist1"]
		#evalProb = samp["evalProb"]
		#spliceIndex = samp["spliceIndex"]

		#angDiff0 = samp["angDiff0"]
		#angDiff1 = samp["angDiff1"]
		#contigFrac0 = samp["contigFrac0"]
		#contigFrac1 = samp["contigFrac1"]
		#spliceU0 = samp["spliceU0"]
		#spliceU1 = samp["spliceU1"]
		#poseLandmarkSum = samp["landmarkSum"]
		#motionBias = samp["motionBias"]
		#evalWithBias = samp["evalWithBias"]
		

		self.drawDist()
		self.drawPoseParticles()

		return


	@logFunction
	def drawDist(self):

		allPathIDs = self.getPathIDs()


		pylab.ioff()
		print pylab.isinteractive()
		pylab.clf()
		#pylab.axis("equal")

		#allSplices = self.getAllSplices2()
		#numSplices = len(allSplices)

		#fig, (ax1, ax2, ax3) = plt.subplots(3,  sharex=False, sharey=False)
		#fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5,  sharex=False, sharey=False)
		#fig, (ax2, ax3, ax4, ax5) = plt.subplots(4,  sharex=True, sharey=False)
		#fig, (ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(6,  sharex=True, sharey=False)
		fig, (ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(7,  sharex=True, sharey=False)
		#fig.set_size_inches(16,12)
		#fig.tight_layout(pad=4.0)
		#fig.tight_layout()

		#fig, axesTuple = plt.subplots(1+numSplices,  sharex=False, sharey=False)

		#ax1 = axesTuple[0]

		#ax1.set_xlim(-3, 3)
		#ax1.set_aspect("equal")

		#ax2.set_ylim(0.0, 1.0)
		#ax2.set_aspect("equal")

		#ax3.set_ylim(-pi, pi)

		""" find number of splices from highest splice index """
		numSplices = 0
		for samp in self.stepResults:
			spliceIndex = samp["spliceIndex"]

			if spliceIndex+1 > numSplices:
				numSplices = spliceIndex+1



		nodeID0 = self.poseData.numNodes-2
		nodeID1 = self.poseData.numNodes-1
		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]

		controlPoses_G = self.getGlobalControlPoses()

		probLists = [[] for k in range(numSplices)]
		landmarkLists = [[] for k in range(numSplices)]
		angDiffLists0 = [[] for k in range(numSplices)]
		overlapLists0 = [[] for k in range(numSplices)]
		contigFracLists0 = [[] for k in range(numSplices)]
		angDiffLists1 = [[] for k in range(numSplices)]
		overlapLists1 = [[] for k in range(numSplices)]
		contigFracLists1 = [[] for k in range(numSplices)]
		biasLists = [[] for k in range(numSplices)]
		finalSumLists = [[] for k in range(numSplices)]


		for samp in self.stepResults:
			hypPose0 = samp["newPose0"]
			hypPose1 = samp["newPose1"]
			currProb0 = samp["displaceProb0"]
			currProb1 = samp["displaceProb1"]
			arcDist0 = samp["arcDist0"]
			arcDist1 = samp["arcDist1"]
			evalProb = samp["evalProb"]
			spliceIndex = samp["spliceIndex"]

			angDiff0 = samp["angDiff0"]
			angDiff1 = samp["angDiff1"]
			overlap0 = samp["overlap0"]
			overlap1 = samp["overlap1"]
			contigFrac0 = samp["contigFrac0"]
			contigFrac1 = samp["contigFrac1"]
			spliceU0 = samp["spliceU0"]
			spliceU1 = samp["spliceU1"]
			poseLandmarkSum = samp["landmarkSum"]
			motionBias = samp["motionBias"]
			evalWithBias = samp["evalWithBias"]

			"""
			poseOrigin0 = Pose(hypPose0)
			xP = []
			yP = []
			for p in medial0:
				p1 = poseOrigin0.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			ax1.plot(xP,yP, color = 'r', alpha = 0.2, zorder=7)

			poseOrigin1 = Pose(hypPose1)

			xP = []
			yP = []
			for p in medial1:
				p1 = poseOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			ax1.plot(xP,yP, color = 'b', alpha = 0.2, zorder=7)
			"""

			#ax2 = axesTuple[1+spliceIndex]

			probLists[spliceIndex].append((arcDist0, evalProb))
			landmarkLists[spliceIndex].append((arcDist0, poseLandmarkSum))
			angDiffLists0[spliceIndex].append((arcDist0, angDiff0))
			overlapLists0[spliceIndex].append((arcDist0, overlap0))
			contigFracLists0[spliceIndex].append((arcDist0, contigFrac0))
			angDiffLists1[spliceIndex].append((arcDist0, angDiff1))
			overlapLists1[spliceIndex].append((arcDist0, overlap1))
			contigFracLists1[spliceIndex].append((arcDist0, contigFrac1))
			biasLists[spliceIndex].append((arcDist0, motionBias))
			finalSumLists[spliceIndex].append((arcDist0, evalWithBias))

			#ax2.scatter([arcDist0, arcDist1], [currProb0, currProb1], color='k')
			#ax2.scatter([arcDist0,], [evalProb,], color='r', zorder=10)


			#resultDict = {}
			#resultDict["newPose0"] = newPose0
			#resultDict["newPose1"] = newPose1
			#resultDict["controlPoses_P"] = controlPoses
			#resultDict["displaceProb0"] = displaceProb0
			#resultDict["displaceProb1"] = displaceProb1
			#resultDict["branchPoses_G"] = branchPoses_G
			#resultDict["branchPoses_L"] = branchPoses_L
			#resultDict["landmarkSum"] = 0.0
			#resultDict["maxLikelihoodBranch"] = None
			#resultDict["branchArcDists"] = {}
			#resultDict["branchControls"] = {}

			#resultDict["angDiff0"] = newAngDiff0
			#resultDict["angDiff1"] = newAngDiff1
			#resultDict["contigFrac0"] = newContigFrac0
			#resultDict["contigFrac1"] = newContigFrac1
			#resultDict["arcDist0"] = arcDist0
			#resultDict["arcDist1"] = arcDist1
			#resultDict["spliceU0"] = spliceU0
			#resultDict["spliceU1"] = spliceU1
		

		#probLists = sorted(probLists, reverse=False)
		for splice in probLists:
			splice.sort()
		for splice in landmarkLists:
			splice.sort()
		for splice in angDiffLists0:
			splice.sort()
		for splice in overlapLists0:
			splice.sort()
		for splice in contigFracLists0:
			splice.sort()
		for splice in angDiffLists1:
			splice.sort()
		for splice in overlapLists1:
			splice.sort()
		for splice in contigFracLists1:
			splice.sort()
		for splice in biasLists:
			splice.sort()
		for splice in finalSumLists:
			splice.sort()

		#print "probLists:", probLists

		#probLists[spliceIndex].append((arcDist0, evalProb))
		#landmarkLists[spliceIndex].append((arcDist0, poseLandmarkSum))
		#angDiffLists[spliceIndex].append((arcDist0, angDiff0))
		#contigFracLists[spliceIndex].append((arcDist0, contigFrac0))
		#biasLists[spliceIndex].append((arcDist0, motionBias))
		#finalSumLists[spliceIndex].append((arcDist0, evalWithBias))

		#self.currMaxIndex = maxPart2
		xP = [self.stepResults[self.currMaxIndex]["arcDist0"],]
		yP = [self.stepResults[self.currMaxIndex]["evalProb"],]
		maxSpliceIndex = self.stepResults[self.currMaxIndex]["spliceIndex"]
		ax2.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		yP = [self.stepResults[self.currMaxIndex]["landmarkSum"],]
		ax3.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		yP = [self.stepResults[self.currMaxIndex]["angDiff0"],]
		ax4.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		yP = [self.stepResults[self.currMaxIndex]["overlap0"],]
		ax5.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		yP = [self.stepResults[self.currMaxIndex]["contigFrac0"],]
		ax6.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		yP = [self.stepResults[self.currMaxIndex]["motionBias"],]
		ax7.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		yP = [self.stepResults[self.currMaxIndex]["evalWithBias"],]
		ax8.scatter(xP,yP, color=self.colors[maxSpliceIndex], zorder=10)

		for k in range(len(probLists)):
			splice = probLists[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			print "plot probLists:", xP#, yP
			ax2.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

		for k in range(len(landmarkLists)):
			splice = landmarkLists[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			#ax3.plot(xP,yP, color = self.colors[k])
			ax3.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

		for k in range(len(angDiffLists0)):
			splice = angDiffLists0[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			ax4.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

			splice = angDiffLists1[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			ax4.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

		for k in range(len(overlapLists0)):
			splice = overlapLists0[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			ax5.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

			splice = overlapLists1[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			ax5.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

		for k in range(len(contigFracLists0)):
			splice = contigFracLists0[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			ax6.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

			splice = contigFracLists1[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			ax6.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

		for k in range(len(biasLists)):
			splice = biasLists[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			#ax6.plot(xP,yP, color = self.colors[k])
			ax7.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')

		for k in range(len(finalSumLists)):
			splice = finalSumLists[k]
			xP = []
			yP = []

			for valTuple in splice:

				xP.append(valTuple[0])
				yP.append(valTuple[1])

			#ax7.plot(xP,yP, color = self.colors[k])
			ax8.plot(xP,yP, color = self.colors[k], marker='x', linestyle='--')


		controlPoses_G = computeGlobalControlPoses(self.getControlPoses(), self.getParentHash())

		#for k, segs in self.globalSegments.iteritems():
		#	for seg in segs:
		#		xP = []
		#		yP = []
		#		for p in seg:
		#			xP.append(p[0])
		#			yP.append(p[1])
		#		ax1.plot(xP,yP, color = self.colors[k], linewidth=4)



		"""
		for sPath in allSplices:
			path = sPath['skelPath']

			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])
			#ax1.plot(xP,yP, color='g', zorder=2)


		xP = []
		yP = []
		for pathID in allPathIDs:
			currFrame = Pose(controlPoses_G[pathID])
			for point_L, pointThresh, pointName in self.localLandmarks[pathID]:
				point_G = currFrame.convertLocalToGlobal(point_L)

				xP.append(point_G[0])
				yP.append(point_G[1])

		if len(xP) > 0:
			pass
			#ax1.scatter(xP, yP, color='k', linewidth=1, zorder=9, alpha=0.9)


		xP = []
		yP = []
		for pathID in allPathIDs:
			terms = self.globalTerms[pathID]
			for term in terms:
				xP.append(term[0])
				yP.append(term[1])

		if len(xP) > 0:
			pass
			#ax1.scatter(xP, yP, color='g', linewidth=1, zorder=11, alpha=0.9)


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
		#ax1.plot(xP,yP,color='k', zorder=10, linewidth=1)

		xP = []
		yP = []
		for p in medial1:
			p1 = poseOrigin1.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		#ax1.plot(xP,yP,color='k', zorder=10, linewidth=1)
		"""

		ax2.set_ylabel('eval')
		ax3.set_ylabel('landmarkCost')
		ax4.set_ylabel('angDiff')
		ax5.set_ylabel('overlap')
		ax6.set_ylabel('contigFrac')
		ax7.set_ylabel('motion')
		ax8.set_ylabel('final')
		
		#pylab.axis("equal")
		#ax1.set_aspect("equal")
		#ax1.set_title("nodeIDs %d %d hypID %d" % (nodeID0, nodeID1, self.hypothesisID))
		fig.suptitle("nodeIDs %d %d hypID %d" % (nodeID0, nodeID1, self.hypothesisID))
		plt.savefig("distEstimate_%04u_%04u.png" % (self.hypothesisID, self.posePlotCount))
		print "distEstimate2_%04u_%04u.png" % (self.hypothesisID, self.posePlotCount)

		plt.clf()
		plt.close()


		self.posePlotCount += 1

	@logFunction
	def drawPoseParticles(self):

		allPathIDs = self.getPathIDs()

		pylab.ioff()
		print pylab.isinteractive()
		pylab.clf()
		#pylab.axis("equal")


		nodeID0 = self.poseData.numNodes-2
		nodeID1 = self.poseData.numNodes-1
		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]


		controlPoses_G = self.getGlobalControlPoses()
		#self.stepResults.append((newPose0, newPose1))

		#for samp in self.stepDists:
		for samp in self.stepResults:
			#self.stepDists.append((j, newDist0, newDist1, newPose0, newPose1))
			#hypPose0 = samp[3]
			#hypPose1 = samp[4]
			hypPose0 = samp["newPose0"]
			hypPose1 = samp["newPose1"]

			poseOrigin0 = Pose(hypPose0)
			xP = []
			yP = []
			for p in medial0:
				p1 = poseOrigin0.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'r', alpha = 0.2, zorder=7)

			poseOrigin1 = Pose(hypPose1)

			xP = []
			yP = []
			for p in medial1:
				p1 = poseOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'b', alpha = 0.2, zorder=7)

		for pathID in allPathIDs:

			xP = []
			yP = []

			if pathID != 0:
				parentPathID = self.pathClasses[pathID]["parentID"]
				currFrame = Pose(controlPoses_G[parentPathID])

				controlCurve = self.controlCurves[pathID]

				for p in controlCurve:
					p1 = currFrame.convertLocalToGlobal(p)
					xP.append(p1[0])
					yP.append(p1[1])

				pylab.plot(xP,yP, color = 'y', alpha = 0.9, zorder=13)
				# = getSkeletonPath(parentSkeleton, parentTerm1, parentTerm2)


		#pylab.scatter(hypPointsX_0, hypPointsY_0, color='r', linewidth=1, zorder=10, alpha=0.2)
		#pylab.scatter(hypPointsX_1, hypPointsY_1, color='b', linewidth=1, zorder=10, alpha=0.2)

		controlPoses_G = computeGlobalControlPoses(self.getControlPoses(), self.getParentHash())
		#landmarks_G = nodeToGlobalLandmarks(self.getControlPoses(), self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses)

		#xP = []
		#yP = []
		#for point_G, pointThresh, pointName in landmarks_G:
		#	xP.append(point_G[0])
		#	yP.append(point_G[1])

		#if len(xP) > 0:
		#	pylab.scatter(xP, yP, color='k', linewidth=1, zorder=9, alpha=0.9)

		for k, segs in self.globalSegments.iteritems():
			for seg in segs:
				xP = []
				yP = []
				for p in seg:
					xP.append(p[0])
					yP.append(p[1])
				pylab.plot(xP,yP, color = self.colors[k], linewidth=4)

		xP = []
		yP = []
		for pathID in allPathIDs:
			currFrame = Pose(controlPoses_G[pathID])
			for point_L, pointThresh, pointName in self.localLandmarks[pathID]:
				point_G = currFrame.convertLocalToGlobal(point_L)

				xP.append(point_G[0])
				yP.append(point_G[1])

		if len(xP) > 0:
			pylab.scatter(xP, yP, color='k', linewidth=1, zorder=9, alpha=0.9)

		xP = []
		yP = []
		for pathID in allPathIDs:
			if pathID != 0:
				pose_G = controlPoses_G[pathID]
				xP.append(pose_G[0])
				yP.append(pose_G[1])

		if len(xP) > 0:
			pylab.scatter(xP, yP, color='r', linewidth=1, zorder=10, alpha=0.9)


		xP = []
		yP = []
		for pathID in allPathIDs:
			terms = self.globalTerms[pathID]
			for term in terms:
				xP.append(term[0])
				yP.append(term[1])

		if len(xP) > 0:
			pylab.scatter(xP, yP, color='g', linewidth=1, zorder=11, alpha=0.9)

		junctionPoses = {}
		for pathID in allPathIDs:
			if pathID != 0:
				junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

		#xP = []
		#yP = []
		#for pathID, branchPoint in junctionPoses.iteritems():
		#	xP.append(branchPoint[0])
		#	yP.append(branchPoint[1])

		#if len(xP) > 0:
		#	pylab.scatter(xP, yP, color='r', linewidth=1, zorder=10, alpha=0.9)


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

		direction = self.poseData.travelDirs[nodeID0]

		if direction:
			"forward"

			p0 = poseOrigin0.convertLocalToGlobal(medial0[0])
			p1 = poseOrigin1.convertLocalToGlobal(medial1[0])

			pylab.scatter([p0[0],p1[0]],[p0[1],p1[1]], zorder=15, color='m')

		else:
			"backward"
			p0 = poseOrigin0.convertLocalToGlobal(medial0[-1])
			p1 = poseOrigin1.convertLocalToGlobal(medial1[-1])

			pylab.scatter([p0[0],p1[0]],[p0[1],p1[1]], zorder=15, color='m')

		xP = []
		yP = []
		for thisNodeID, colPoint in self.collisionLandmarks.iteritems():
			thisPose = self.getNodePose(thisNodeID)
			thisPoseFrame = Pose(thisPose)
			globalColPoint = thisPoseFrame.convertLocalToGlobal(colPoint)
			xP.append(globalColPoint[0])
			yP.append(globalColPoint[1])

		pylab.scatter(xP,yP, zorder=14, color='b')

		
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

		#newObj.poseParticles = deepcopy(self.poseParticles)

		#updateCount = newObj.poseParticles["updateCount"] 
		#particleDist2 = deepcopy(newObj.poseParticles["snapshots2"][0])
		#for part in particleDist2:
		#	part.mapStateID = hypothesisID

		newObj.parentParticle = self.hypothesisID 

		newObj.pathClasses = deepcopy(self.pathClasses)
		newObj.pathTermsVisited = deepcopy(self.pathTermsVisited)
		newObj.pathTermsData = deepcopy(self.pathTermsData)
		newObj.pathTerms = deepcopy(self.pathTerms)
		newObj.pathIDs = deepcopy(self.pathIDs)
		newObj.paths = deepcopy(self.paths)
		newObj.hulls = deepcopy(self.hulls)
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
		
		newObj.allSplices = deepcopy(self.allSplices)
		newObj.allSplices2 = deepcopy(self.allSplices2)

		#newObj.orderedPathIDs1 = deepcopy(self.orderedPathIDs1)
		#newObj.orderedPathIDs2 = deepcopy(self.orderedPathIDs2)

		" departure events for node 1 "
		#newObj.departures1 = deepcopy(self.departures1) 
		#newObj.interiors1 = deepcopy(self.interiors1)  
		#newObj.depPoints1 = deepcopy(self.depPoints1) 
		#newObj.distances1 = deepcopy(self.distances1)  
		#newObj.depAngles1 = deepcopy(self.depAngles1)
		#newObj.contig1 = deepcopy(self.contig1)

		" departure events for node 2 "
		#newObj.departures2 = deepcopy(self.departures2) 
		#newObj.interiors2 = deepcopy(self.interiors2) 
		#newObj.depPoints2 = deepcopy(self.depPoints2)  
		#newObj.distances2 = deepcopy(self.distances2) 
		#newObj.depAngles2 = deepcopy(self.depAngles2)  
		#newObj.contig2 = deepcopy(self.contig2)  

		""" branching primitive data """
		newObj.departureResultSet1 = deepcopy(self.departureResultSet1)  
		newObj.departureResultSet2 = deepcopy(self.departureResultSet2)
		newObj.overlapSplice1 = deepcopy(self.overlapSplice1)
		newObj.overlapSplice2 = deepcopy(self.overlapSplice2)
		newObj.memberShootIDs1 = deepcopy(self.memberShootIDs1)
		newObj.memberShootIDs2 = deepcopy(self.memberShootIDs2)

		newObj.localPaths = deepcopy(self.localPaths)
		newObj.localHulls = deepcopy(self.localHulls)
		newObj.localLandmarks = deepcopy(self.localLandmarks)
		newObj.localLeaf2LeafPathJunctions = deepcopy(self.localLeaf2LeafPathJunctions) 

		newObj.localLongPaths = deepcopy(self.localLongPaths)
		newObj.globalLongPaths = deepcopy(self.globalLongPaths)
		newObj.localSegments = deepcopy(self.localSegments)
		newObj.globalSegments = deepcopy(self.globalSegments)
		newObj.localTerms = deepcopy(self.localTerms)
		newObj.globalTerms = deepcopy(self.globalTerms)
		newObj.subsumedTerms = deepcopy(self.subsumedTerms)

		newObj.branchDiverges = deepcopy(self.branchDiverges)

		return newObj

	@logFunction
	def computeEval(self):

		parentHash = self.getParentHash()

		totalSum = 0.0

		if False in self.branchDiverges.values():
				totalSum = 1e100

		""" look for outlier landmarks to invalidate a map state """
		LANDMARK_THRESH = 7.0
		#OUTLIER_THRESH = 1.5
		OUTLIER_THRESH = 0.7
		#OUTLIER_THRESH = 0.4

		print "controlPoses_L:", self.getControlPoses()

		controlPoses_G = computeGlobalControlPoses(self.getControlPoses(), self.getParentHash())

		print "controlPoses_G:", controlPoses_G

		""" collect the current landmarks """
		landmarks_G = nodeToGlobalLandmarks(self.getControlPoses(), self.getPathIDs(), self.getParentHash(), self.nodeLandmarks, self.pathClasses)

		print "check for outliers:", landmarks_G

		print "computeEval:", self.hypothesisID



		allNearestNeighbors = []
		bloomNearestNeighbors = []

		for j in range(len(landmarks_G)):

			p1 = landmarks_G[j][0]
			thresh1 = landmarks_G[j][1]
			landmarkType1 = landmarks_G[j][2]

			thisLandmarksBloomNeighbors = []
			thisLandmarksNeighbors = []

			for k in range(0, len(landmarks_G)):
				if j != k:
					p2 = landmarks_G[k][0]
					thresh2 = landmarks_G[k][1]
					landmarkType2 = landmarks_G[k][2]

					dist = sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
					diffVal = sqrt(dist*dist/(thresh1*thresh1 + thresh2*thresh2))


					if (landmarkType1 == "archPoint" or landmarkType1 == "bloomPoint") and (landmarkType2 == "archPoint" or landmarkType2 == "bloomPoint"):
						print landmarkType1, landmarkType2
						""" bloom only nearest neighbors """
						thisLandmarksBloomNeighbors.append((j, k, dist, diffVal, landmarks_G[j][1], landmarks_G[k][1], landmarks_G[j][0], landmarks_G[k][0]))

					""" all nearest neighbors """
					thisLandmarksNeighbors.append((j, k, dist, diffVal, landmarks_G[j][1], landmarks_G[k][1], landmarks_G[j][0], landmarks_G[k][0]))

			thisLandmarksBloomNeighbors = sorted(thisLandmarksBloomNeighbors, key=itemgetter(2), reverse=False)
			thisLandmarksNeighbors = sorted(thisLandmarksNeighbors, key=itemgetter(2), reverse=False)
		
			if len(thisLandmarksNeighbors) > 0:
				allNearestNeighbors.append(thisLandmarksNeighbors)

			if len(thisLandmarksBloomNeighbors) > 0:
				bloomNearestNeighbors.append(thisLandmarksBloomNeighbors)


		isReject = False
		print self.hypothesisID, "bloom neighbors:"
		for result in bloomNearestNeighbors:

			print result

			if len(result) > 1:
				print result[0][2], result[1][2]
				if result[0][2] > 0.4 or result[1][2] > 0.4:
					isReject = True
			else:
				print result[0][2]
				if result[0][2] > 0.4:
					isReject = True

		if isReject:
			totalSum = 5e100

		print self.hypothesisID, "bloom/bend neighbors:"
		for result in allNearestNeighbors:

			print result

			if len(result) > 1:
				print result[0][2], result[1][2]
				if result[0][2] > 1.0 or result[1][2] > 1.0:
					isReject = True
			else:
				print result[0][2]
				if result[0][2] > 1.0:
					isReject = True

		if isReject:
			totalSum = 6e100


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
	def setLocalTipPoint(self, pathID, tipPoint_L):

		if pathID != 0:

			self.pathClasses[pathID]["tipPoint_L"] = tipPoint_L

			globalControlPoses = computeGlobalControlPoses(self.getControlPoses(), self.getParentHash())

			#""" convert the controlPose to coordinates local to the parent frame """
			shootControlPose = globalControlPoses[pathID]
			currFrame = Pose(shootControlPose)
			tipPoint_G = currFrame.convertLocalToGlobal(tipPoint_L)

			print "tipPoint_G =", self.hypothesisID, tipPoint_G

			#self.pathClasses[pathID]["globalJunctionPose"] = globalJuncPose

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

	def getChildHash(self):

		" determine which paths are descendants "
		pathIDs = self.getPathIDs()
		childPathIDs = {}
		for k in pathIDs:
			childPathIDs[k] = []

		parentPathIDs = self.getParentHash()

		for childID in pathIDs:
			parentID = parentPathIDs[childID]

			if parentID != None:
				childPathIDs[parentID].append(childID)

		return childPathIDs

	@logFunction
	def getGlobalBranchPoses(self):
		branchPoses_G = {}
		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:
			branchPoses_G[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

		return branchPoses_G

	@logFunction
	def getLocalBranchPoses(self):
		branchPoses_L = {}
		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:
			branchPoses_L[pathID] = self.pathClasses[pathID]["localJunctionPose"]

		return branchPoses_L

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
	def getGlobalControlPoses(self):
		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)
		return globalControlPoses


	@logFunction
	def getParentPathID(self, pathID):
		parentID = self.pathClasses[pathID]["parentID"]
		return parentID
	
	@logFunction
	def getPath(self, pathID):
		return self.pathClasses[pathID]
	
	def getAllSplices2(self, plotIter = False):
		if not self.isChanged:
			print "returning without computing, not isChanged"
			return self.allSplices2
		else:
			print "recomputing allSplices"
			self.generatePaths()
			return self.allSplices2

	def computeAllSplices2(self, plotIter = False):

		"""
		1) get all terminals of each skeleton
		2) eliminate terminals that are "subsumed" by other shoots
		3) find splices between the remaining terminals


		FIXME: how to eliminate so many possibilities when junctions are not aligned?
		FIXME: what if two terminals are very close to one another, how to select only one? not both or zero

		"""
		#DIST_THRESH = 0.2
		DIST_THRESH = 0.4
		DIVERGE_THRESH = 0.8

		
		""" of all the terms, find the ones that are not subsumed by other shoots """

		self.branchDivergeCount = {}
		self.branchSubsumeIDs = {}
		self.branchTermDivergenceDist = {}

		self.subsumedTerms_L = {}
		allTerms_L = {}
		allTerms_G = {}
		currKeys = self.globalSegments.keys()

		self.subsumptionMatrix = {}
		self.terminalSimilarity = {}
		for pathID1 in currKeys:
			self.subsumptionMatrix[pathID1] = {}
			self.terminalSimilarity[pathID1] = {}

			for pathID2 in currKeys:
				if pathID1 != pathID2:
					self.subsumptionMatrix[pathID1][pathID2] = 0
					self.terminalSimilarity[pathID1][pathID2] = 0





		controlPoses_G = self.getGlobalControlPoses()

		#for pathID, terms in self.localTerms.iteritems():
		#	self.globalTerms[pathID] = []
		#	shootFrame = Pose(finalPoses[pathID])
		#	for term in terms:
		#		globalP = shootFrame.convertLocalToGlobal(term)
		#		self.globalTerms[pathID].append(globalP)

		for currK1 in range(len(currKeys)): 

			pathID1 = currKeys[currK1]
			self.branchDivergeCount[pathID1] = 0
			self.branchSubsumeIDs[pathID1] = []
			self.branchTermDivergenceDist[pathID1] = []
			segs1 = self.globalSegments[pathID1]

			terms_L = self.localTerms[pathID1]


			shootFrame = Pose(controlPoses_G[pathID1])

			allTerms_L[pathID1] = []
			allTerms_G[pathID1] = []


			for term1 in terms_L:

				term1_G = shootFrame.convertLocalToGlobal(term1)

				minDist2 = 1e100
				minP2 = None
				subSkelID = None

				for currK2 in range(len(currKeys)): 

					if currK2 != currK1:

						pathID2 = currKeys[currK2]
						segs2 = self.globalSegments[pathID2]

						thisMinDist = 1e100

						for seg2 in segs2:
							for p in seg2:

								dist2 = sqrt((p[0]-term1_G[0])**2 + (p[1]-term1_G[1])**2)

								if dist2 < thisMinDist:
									thisMinDist = dist2

								if dist2 < minDist2:
									minDist2 = dist2
									minP2 = p
									subSkelID = pathID2
						
						if thisMinDist <= DIST_THRESH:
							self.subsumptionMatrix[pathID1][pathID2] += 1

						otherTerms_G = self.globalTerms[pathID2]
						for p in otherTerms_G:

							dist = sqrt((p[0]-term1_G[0])**2 + (p[1]-term1_G[1])**2)

							if dist < DIST_THRESH:
								self.terminalSimilarity[pathID1][pathID2] += 1

				
				self.branchTermDivergenceDist[pathID1].append(minDist2)
				print pathID1, minDist2, "terms:", term1_G
				if minDist2 > DIST_THRESH:
					allTerms_L[pathID1].append(term1)
					allTerms_G[pathID1].append(term1_G)
					#print pathID1, "terms:", term1
				else:
					#if pathID1 != 0:
					#	self.branchDivergeCount[pathID1] += 1
					self.branchDivergeCount[pathID1] += 1
					self.branchSubsumeIDs[pathID1].append(subSkelID)


					""" check if this terminal is close to another terminal in different shoot """
					isCoincident = False
					for currK2 in range(currK1+1, len(currKeys)): 
						pathID2 = currKeys[currK2]
						otherTerms_G = self.globalTerms[pathID2]
						for p in otherTerms_G:

							dist = sqrt((p[0]-term1_G[0])**2 + (p[1]-term1_G[1])**2)

							if dist < DIST_THRESH:
								isCoincident = True
					#if isCoincident:
					#	allTerms_L[pathID1].append(term1)
					#	allTerms_G[pathID1].append(term1_G)

				
					isUnique = True
					for pathID, otherTerms_G in allTerms_G.iteritems():
						for p in otherTerms_G:
							dist = sqrt((p[0]-term1_G[0])**2 + (p[1]-term1_G[1])**2)

							if dist < DIST_THRESH:
								print term1_G, "is similar to", p
								isUnique = False

					if isCoincident and isUnique:
						print "adding", term1_G
						allTerms_L[pathID1].append(term1)
						allTerms_G[pathID1].append(term1_G)





		self.subsumedTerms_L = allTerms_L


		""" previous terminal visitation status """
		prevTermsData = self.pathTermsData
		#prevTermsVisited = self.pathTermsVisited
		#self.pathTermsVisited = {}
		self.pathTermsData = {}

		allPathIDs = self.getPathIDs()

		""" previous terminals from last computation """
		prevTermList = []
		for pathID, terms_L in self.subsumedTerms_L.iteritems():

			shootFrame = Pose(controlPoses_G[pathID])


			for term_L in terms_L:

				term1_G = shootFrame.convertLocalToGlobal(term_L)

				prevTermList.append(term1_G)

				minDist = 1e100
				minMatchData = None
				minUUID = None

				#for termTuple, visitStatus in prevTermsVisited.iteritems():
				for termUUID, termData in prevTermsData.iteritems():
					termShootID = termData["frameID"]
					prevTermPoint = termData["term_L"]
					visitStatus = termData["isVisited"]

					if termShootID in allPathIDs:
						termFrame = Pose(controlPoses_G[termShootID])

						prevTermPoint_G = termFrame.convertLocalToGlobal(prevTermPoint)

						prevDist = sqrt((term1_G[0]-prevTermPoint_G[0])**2 + (term1_G[1]-prevTermPoint_G[1])**2)

						if prevDist < minDist:
							minDist = prevDist
							minMatchData = termData
							minUUID = termUUID

				if minDist < 1.0:

					try:
						""" check for prior matches """
						bestDist = self.pathTermsData[minUUID][4]
						if minDist < bestDist:

							""" overwrite previous match """
							dataVal = {}
							dataVal["frameID"] = pathID
							dataVal["term_L"] = tuple(term_L)
							dataVal["isVisited"] = minMatchData["isVisited"]
							dataVal["term_G"] = tuple(term1_G)
							dataVal["minDist"] = minDist
							self.pathTermsData[minUUID] = dataVal

						else:
							dataVal = {}
							dataVal["frameID"] = pathID
							dataVal["term_L"] = tuple(term_L)
							dataVal["isVisited"] = False
							dataVal["term_G"] = tuple(term1_G)
							dataVal["minDist"] = None
							self.pathTermsData[uuid4().int] = dataVal

					except:
						""" no other matches, create new one """
						dataVal = {}
						dataVal["frameID"] = pathID
						dataVal["term_L"] = tuple(term_L)
						dataVal["isVisited"] = minMatchData["isVisited"]
						dataVal["term_G"] = tuple(term1_G)
						dataVal["minDist"] = minDist
						self.pathTermsData[minUUID] = dataVal
						#(pathID, tuple(term_L), minMatchTuple[2], term1_G, minDist)

				else:
					dataVal = {}
					dataVal["frameID"] = pathID
					dataVal["term_L"] = tuple(term_L)
					dataVal["isVisited"] = False
					dataVal["term_G"] = tuple(term1_G)
					dataVal["minDist"] = None
					self.pathTermsData[uuid4().int] = dataVal
					#(pathID, tuple(term_L), False, term1_G, None)
					#self.pathTermsVisited[(pathID, tuple(term_L))] = False

		
		""" data structures for navigation """
		self.pathTermsVisited = {}
		self.pathTerms = {}
		for uuidVal, termData in self.pathTermsData.iteritems():
			self.pathTermsVisited[uuidVal] = termData["isVisited"]
			self.pathTerms[uuidVal] = termData["term_G"]

	
		print "self.pathTermsVisited:", self.pathTermsVisited
		print "self.pathTerms:", self.pathTerms
		print "self.pathTermsData:", self.pathTermsData

		termsVisited = self.getPathTermsVisited()
		print "self.getPathTerms():", self.getPathTerms()
		print "self.getPathTermsVisited:", self.getPathTermsVisited()

		termList = []
		pathIDs = allTerms_G.keys()

		for pathID in pathIDs:
			for term_G in allTerms_G[pathID]:
				termList.append(term_G)
		


		termCombos = []
		for j in range(len(termList)):
			for k in range(j+1, len(termList)):
				termCombos.append((termList[j], termList[k]))


		finalResults = []


		print "termCombos:", termCombos
		
		for termPath in termCombos:

			memberShootIDs = {}
			joinPairs = []

			startPose = termPath[0]
			endPose = termPath[1]

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

			nodeAttrs = self.spliceSkeleton.get_node_attributes(startNode)
			memberPathID = None
			for attr in nodeAttrs:
				if attr[0] == "pathID":
					memberPathID = attr[1]

			memberShootIDs[memberPathID] = None

			splicedSkel = [startNode]
			while currNode != endNode:
				#print "currNode:", currNode
				splicedSkel.append(currNode)

				nodeAttrs = self.spliceSkeleton.get_node_attributes(currNode)
				memberPathID = None
				for attr in nodeAttrs:
					if attr[0] == "pathID":
						memberPathID = attr[1]

				memberShootIDs[memberPathID] = None

				nextNode = shortestSpliceTree[currNode]
				currNode = nextNode
			splicedSkel.append(currNode)

			splicedSkel = ensureEnoughPoints(splicedSkel, max_spacing = 0.08, minPoints = 5)
			spliceSpline1 = SplineFit(splicedSkel, smooth=0.1)
			splicePoints1 = spliceSpline1.getUniformSamples()


			sPath = {}
			sPath['termPath'] = termPath
			sPath['skelPath'] = splicePoints1
			sPath['memberShootIDs'] = memberShootIDs.keys()
			
			finalResults.append(sPath)
			
		
		if plotIter:
			print "plotting splicedPath2"
			pylab.clf()

			
			xP = []
			yP = []
			for term in termList:
				xP.append(term[0])
				yP.append(term[1])

			pylab.scatter(xP, yP, color='g', linewidth=1, zorder=11, alpha=0.9)

	
	
			for sPath in finalResults:
				path = sPath['skelPath']

				xP = []
				yP = []
				for p in path:
					xP.append(p[0])
					yP.append(p[1])
				pylab.plot(xP,yP, color='b', zorder=2)


			#self.drawWalls()
	
			print "saving splicedPath_%04u_%04u_%04u.png" % (self.hypothesisID, self.poseData.numNodes, self.spliceCount)
			pylab.title("Spliced Paths, pathIDs = %s" %  self.getPathIDs())
			pylab.savefig("newSplicedPath_%04u_%04u_%04u.png" % (self.hypothesisID, self.poseData.numNodes, self.spliceCount))
			self.spliceCount += 1

		self.allSplices2 = finalResults


	@logFunction
	def isNodeExist(self, nodeID, pathID):
		return nodeID in self.pathClasses[pathID]["nodeSet"]

	
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
	 

	def snapToParent(self):

		allPathIDs = self.getPathIDs()
		parentPathIDs = self.getParentHash()

		controlPoses_G = self.getGlobalControlPoses()

		originPoint = (0.0,0.0)

		for pathID in allPathIDs:
			if self.pathClasses[pathID]["parentID"] != None:
				#childSkeleton = self.localSkeletons[pathID]

				segs1 = self.localSegments[pathID]

				minDist = 1e100
				minP = None

				for seg in segs1:
					for p in seg:
						dist = sqrt((originPoint[0]-p[0])**2+(originPoint[1]-p[1])**2)

						if dist < minDist:
							minDist = dist
							minP = p

				""" new origin of the current frame """
				newOriginPose = (minP[0], minP[1], 0.0)
				newOrigin = Pose(newOriginPose)

				""" existing frame """

				currFrame = Pose(controlPoses_G[pathID])

				for nodeID in self.pathClasses[pathID]["nodeSet"]:

					""" we just move the existing poses within the frame, and keep the frame as is """
					oldNodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]
					newNodePose_L = newOrigin.convertGlobalPoseToLocal(oldNodePose_L)
					self.pathClasses[pathID]["localNodePoses"][nodeID] = newNodePose_L


					newPose_G = currFrame.convertLocalOffsetToGlobal(newNodePose_L)


					""" for recomputing raw pose """
					oldGPACPose = self.nodePoses[nodeID]
					gpacProfile = Pose(oldGPACPose)
					localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[nodeID])

					

					""" new global pose """
					self.nodePoses[nodeID] = newPose_G
					
					
					" go back and convert this from GPAC pose to estPose "
					newProfile = Pose(newPose_G)
					newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
					
					self.nodeRawPoses[nodeID] = newEstPose

 

	def snapPoseToSkeleton(self, targetNodeIDs = []):

		allPathIDs = self.getPathIDs()
		controlPoses_G = self.getGlobalControlPoses()

		originPoint = (0.0,0.0)

		newNodePoses = {}

		landmarkNodeSets = {}

		for pathID in allPathIDs:

			""" existing frame """
			currFrame = Pose(controlPoses_G[pathID])

			segs1 = self.localSegments[pathID]
			longPaths_L = self.localLongPaths[pathID]

			#newNodePoses[pathID] = {}

			landmarkNodeSets[pathID] = []


			for nodeID in self.pathClasses[pathID]["nodeSet"]:

				if self.nodeLandmarks[pathID][nodeID] != None:

					if len(targetNodeIDs) == 0 or nodeID in targetNodeIDs:
						landmarkNodeSets[pathID].append(nodeID)


		#time1 = time.time()
		newNodePoses = self.localizeLandmarkPose(landmarkNodeSets)

		#newNodePose_L = self.localizeLandmarkPose(pathID, nodeID)
		#for pathID, nodePoses in newNodePoses.iteritems():
		#	for nodeID, nodePose_L in nodePoses.iteritems():
		#		if newNodePose_L != None:
		#			newNodePoses[pathID][nodeID] = newNodePose_L

		for pathID, nodePoses in newNodePoses.iteritems():
				
			for nodeID, newNodePose_L in nodePoses.iteritems():

				if newNodePose_L != None:

					# = (landmarkPoint_N, BLOOM_THRESH, "bloomPoint")
					oldNodePose_L = self.pathClasses[pathID]["localNodePoses"][nodeID]

					""" we just move the existing poses within the frame, and keep the frame as is """
					#newNodePose_L[2] = oldNodePose_L[2]

					self.pathClasses[pathID]["localNodePoses"][nodeID] = newNodePose_L

					newPose_G = currFrame.convertLocalOffsetToGlobal(newNodePose_L)

					""" for recomputing raw pose """
					oldGPACPose = self.nodePoses[nodeID]
					gpacProfile = Pose(oldGPACPose)
					localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[nodeID])

					""" new global pose """
					self.nodePoses[nodeID] = newPose_G
					
					
					" go back and convert this from GPAC pose to estPose "
					newProfile = Pose(newPose_G)
					newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
					
					self.nodeRawPoses[nodeID] = newEstPose



	@logFunction
	def resetBranches(self):


		self.branchArcBins = {}
		self.jointBranchEvaluations = {}

		parentPathIDs = self.getParentHash()



		#self.DIV_LEN = 0.2

		allPathIDs = self.getPathIDs()


		self.controlCurves = {}
		for pathID in allPathIDs:
			if self.pathClasses[pathID]["parentID"] != None:
				parentPathID = parentPathIDs[pathID]
				parentSkeleton = self.localSkeletons[parentPathID]
				parentTerm1 = self.pathClasses[pathID]["controlTerm1_P"]
				parentTerm2 = self.pathClasses[pathID]["controlTerm2_P"]

				self.controlCurves[pathID] = getSkeletonPath(parentSkeleton, parentTerm1, parentTerm2)
			else:
				self.controlCurves[pathID] = self.localPaths[pathID]




		for pathID in allPathIDs:
			if self.pathClasses[pathID]["parentID"] != None:
				parentPathID = parentPathIDs[pathID]
				#pathSpline = SplineFit(self.localPaths[parentPathID])
				pathSpline = SplineFit(self.controlCurves[pathID])

				totalDist = pathSpline.dist_u(1.0)

				branchSpace = {}

				currDist = 0.0
				branchSpace[currDist] = None
				while currDist <= totalDist:
					currDist += self.DIV_LEN
					branchSpace[currDist] = None

				self.branchArcBins[pathID] = branchSpace


		self.branchControlPoses = {}
		self.branchArcDists = {}
		for pathID in allPathIDs:
			if pathID != 0:

				self.branchControlPoses[pathID] = {}
				self.branchArcDists[pathID] = []

				parentPathID = parentPathIDs[pathID]
				#pathSpline = SplineFit(self.localPaths[parentPathID])
				pathSpline = SplineFit(self.controlCurves[pathID])

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
			currBranchSpace = self.branchArcBins[pathID]

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

		""" sorted list of pathIDs for indexing arc distances in order to make indexable tuples """
		allPathIDs = self.getPathIDs()
		branchPathIDs = deepcopy(allPathIDs)
		branchPathIDs.remove(0)
		branchPathIDs.sort()


		""" data needed for branch computation processes """
		localSkeletons = {}
		controlPoses = {}
		junctionPoses = {}
		tipPoints = {}
		localPathSegsByID = {}
		localTermsByID = {}
		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			localSkeletons[pathID] = self.localLeaf2LeafPathJunctions[pathID]["skeletonGraph"]
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			junctionPoses[pathID] = self.pathClasses[pathID]["localJunctionPose"]
			tipPoints[pathID] = self.pathClasses[pathID]["tipPoint_L"]
			localPathSegsByID[pathID] = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]
			localTermsByID[pathID] = self.localLeaf2LeafPathJunctions[pathID]["leafTerms"]
			#self.localTerms[pathID] = self.localLeaf2LeafPathJunctions[pathID]["leafTerms"]
			#self.localSegments[pathID] = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]

		parentPathIDs = self.getParentHash()

		""" construct the joint control pose problem sets """
		jointBranchJobs = []
		for controlSet in jointComboControlSets:

			""" values are overwritten by controlSet """
			thisControlPoses = deepcopy(controlPoses)
			thisArcDists = {}

			print "controlSet:", controlSet

			for k in range(len(controlSet)):

				arcDist = controlSet[k][0]
				cPose = controlSet[k][1]

				pathID = branchPathIDs[k]
				thisControlPoses[pathID] = cPose
				thisArcDists[pathID] = arcDist

			print "thisArcDists:", thisArcDists

			#self.localLandmarks
			jointBranchJobs.append((localPathSegsByID, localTermsByID, self.localPaths, localSkeletons, thisControlPoses, tipPoints, junctionPoses, self.localLandmarks, parentPathIDs, thisArcDists, len(self.nodePoses)-1, self.hypothesisID ))
			#jointBranchJobs.append((localPathSegsByID, localTermsByID, self.controlCurves, localSkeletons, thisControlPoses, tipPoints, junctionPoses, self.localLandmarks, parentPathIDs, thisArcDists, len(self.nodePoses)-1, self.hypothesisID ))

		
		#evaluateResults = []
		#for job in jointBranchJobs:
			#evaluateJointBranch(localPathSegsByID, localTerms, localPaths, localSkeletons, controlPoses, tipPoints, junctionPoses, landmarks, parentPathIDs, arcDists, numNodes=0, hypothesisID=0)
			#time1 = time.time()
			#evaluateResults.append(evaluateJointBranch(job[0], job[1], job[2], job[3], job[4], job[5], job[6], job[7], job[8], job[9], numNodes=job[10], hypothesisID=job[11]))
			#time2 = time.time()

			#print "TIME evaluateJointBranch:", time2-time1

		time1 = time.time()
		evaluateResults = batchJointBranch(jointBranchJobs)
		time2 = time.time()
		print "TIME batch evaluateJointBranch:", time2-time1

		print len(jointBranchJobs), "total branch jobs"
		#time1 = time.time()
		#jointResults = batchJointBranch(jointBranchJobs)
		#time2 = time.time()
		#print "TIME batch batchJointBranch:", time2-time1

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
		for k in range(len(evaluateResults)):
			branchResult = evaluateResults[k]

			totalMatchCount = branchResult["totalMatchCount"]
			totalCost = branchResult["totalCost"]
			landmarkCost = branchResult["landmarkCost"]

			if landmarkCost > self.maxLandmarkCost:
				self.maxLandmarkCost = landmarkCost

			if totalCost > self.maxCost:
				self.maxCost = totalCost

			if totalMatchCount > self.maxMatchCount:
				self.maxMatchCount = totalMatchCount

		print "max values:", self.maxLandmarkCost, self.maxCost, self.maxMatchCount

		""" add slight value so that maximum normalized landmark cost is greater than zero """
		self.maxLandmarkCost += 0.1

		for k in range(len(evaluateResults)):
			branchResult = evaluateResults[k]

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

			""" multiply the overlap with the landmark utility """
			normLandmark *= normMatchCount

			""" if branch is not diverging, then zero the landmark utility since we use this for evaluation """
			if self.maxMatchCount <= 0:
				normLandmark = 0.0

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

		maxTuple = None
		maxNormCost = -1e10
		maxMatchCount = -1e10
		maxLandmark = -1e10
		printList = []
		for arcTuple, branchResult in self.jointBranchEvaluations.iteritems():

			normLandmark = branchResult["normLandmark"]
			normMatchCount = branchResult["normMatchCount"]
			normCost = branchResult["normCost"] 
			controlSet = branchResult["controlSet"]
			landmarkCost = branchResult["landmarkCost"]

			matchCount = branchResult["totalMatchCount"]

			printList.append((self.poseData.numNodes, arcTuple, normLandmark, landmarkCost, normMatchCount, matchCount, normCost, controlSet))
			#print "branch eval:", self.poseData.numNodes, arcTuple, normLandmark, normMatchCount, normCost, controlSet

			#if normLandmark > maxNormCost:
			if normLandmark > maxLandmark:
				maxLandmark = normLandmark
				maxMatchCount = normMatchCount
				maxNormCost = normCost
				maxTuple = arcTuple

		printList.sort()
		for val in printList:
			print "branch eval:", val

		if maxTuple != None:
			maxBranchResult = self.jointBranchEvaluations[maxTuple]
			print "branch max:", self.poseData.numNodes, maxLandmark, maxBranchResult["landmarkCost"], maxMatchCount, maxNormCost, maxTuple, maxBranchResult["controlSet"], self.maxCost, self.maxMatchCount, self.maxLandmarkCost


		if maxTuple != None:
			""" only splice the maximum """
			maxBranchResult = self.jointBranchEvaluations[maxTuple]
			thisArcDists = maxBranchResult["arcDists"]
			thisControlPoses = maxBranchResult["controlSet"]
			maxResult = computeJointBranch(localPathSegsByID, localTermsByID, self.localPaths, localSkeletons, thisControlPoses, tipPoints, junctionPoses, self.localLandmarks, parentPathIDs, thisArcDists, len(self.nodePoses)-1, self.hypothesisID )

			maxResult["normMatchCount"] = maxMatchCount
			maxResult["normCost"] = maxNormCost
			maxResult["normLandmark"] = maxLandmark

			self.jointBranchEvaluations[maxTuple] = maxResult

		
		if False:

			parentPathIDs = self.getParentHash()

			pylab.clf() 

			for arcTuple, branchResult in self.jointBranchEvaluations.iteritems():

				#modJuncPose = result["modJuncPose"]
				#modControlPose = result["modControlPose"]

				controlPoses_L = branchResult["controlSet"]
				controlPoses_G = computeGlobalControlPoses(controlPoses_L, parentPathIDs)

				#branchPoses_L = branchResult["branchPoses_L"]


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

				#for k,path in self.trimmedPaths.iteritems():
				#	xP = []
				#	yP = []
				#	for p in path:
				#		xP.append(p[0])
				#		yP.append(p[1])

				#	pylab.plot(xP,yP, color = self.colors[k], linewidth=4)


				xC = []
				yC = []
				xB = []
				yB = []
				for pathID in branchPathIDs:
					controlPose_G = controlPoses_G[pathID]
					currFrame = Pose(controlPose_G)

					xC.append(controlPose_G[0])
					yC.append(controlPose_G[1])

					#branchPose_L = branchPoses_L[pathID]
					#branchPose_G = currFrame.convertLocalOffsetToGlobal(branchPose_L)

					#xB.append(branchPose_G[0])
					#yB.append(branchPose_G[1])
			
				#pylab.scatter(xB, yB, color='y', zorder=8)
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

				if spatialFeature["inflectionPoint"] != None:
					print self.hypothesisID, "adding inflectionPoint for node", nodeID, "in path", pathID
					landmarkPoint_N = spatialFeature["inflectionPoint"]
					landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
					self.localLandmarks[pathID].append((landmarkPoint_L, BEND_THRESH, "bendPoint"))
					self.nodeLandmarks[pathID][nodeID] = (landmarkPoint_N, BEND_THRESH, "bendPoint")

				elif spatialFeature["bloomPoint"] != None:
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

		self.localLongPaths = {}
		self.globalLongPaths = {}
		self.localSegments = {}
		self.globalSegments = {}
		self.localTerms = {}
		self.globalTerms = {}
		
		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			print "computing path for node set", pathID, ":", self.getNodes(pathID)
			localResults = computeShootSkeleton(self.poseData,
										pathID,
										self.pathClasses[pathID]["localJunctionPose"],
										self.getNodes(pathID),
										self.pathClasses[pathID]["localNodePoses"],
										self.localLandmarks[pathID],
										self.pathClasses[pathID]["tipPoint_L"],
										self.hypothesisID,
										self.colors[pathID],
										self.topCount,
										plotIter = True)

			self.localLeaf2LeafPathJunctions[pathID] = localResults[0]
			self.localPaths[pathID] = localResults[1]
			self.localHulls[pathID] = localResults[2]




			self.topCount += 1
			
		self.localSkeletons = {}
		controlPoses = {}
		junctionPoses = {}
		for pathID in pathIDs:
			self.localSkeletons[pathID] = self.localLeaf2LeafPathJunctions[pathID]["skeletonGraph"]
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]
			self.localTerms[pathID] = self.localLeaf2LeafPathJunctions[pathID]["leafTerms"]
			self.localSegments[pathID] = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]
			self.localLongPaths[pathID] = self.localLeaf2LeafPathJunctions[pathID]["longPaths"]


		#leaf2LeafPathJunctions["leafTerms"] = smoothLeafTerms
		#leaf2LeafPathJunctions["localSegments"] = localPathSegs


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

		for pathID, segs in self.localSegments.iteritems():
			self.globalSegments[pathID] = []
			shootFrame = Pose(finalPoses[pathID])
			for seg in segs:
				globalSeg = []
				for p in seg:
					globalP = shootFrame.convertLocalOffsetToGlobal(p)
					globalSeg.append(globalP)

				self.globalSegments[pathID].append(globalSeg)

		for pathID, longPaths in self.localLongPaths.iteritems():
			self.globalLongPaths[pathID] = []
			shootFrame = Pose(finalPoses[pathID])
			for path in longPaths:

				globalPath = []
				for p in path:
					globalP = shootFrame.convertLocalOffsetToGlobal(p)
					globalPath.append(globalP)

				self.globalLongPaths[pathID].append(globalPath)


		for pathID, terms in self.localTerms.iteritems():
			self.globalTerms[pathID] = []
			shootFrame = Pose(finalPoses[pathID])
			for term in terms:
				globalP = shootFrame.convertLocalToGlobal(term)
				self.globalTerms[pathID].append(globalP)

		#self.trimmedPaths = self.trimPaths(self.paths)


		self.spliceSkeleton = spliceSkeletons(self.localSkeletons, finalPoses, junctionPoses, parentPathIDs)

		""" find root terminal """
		rootPath = self.paths[0]
		dist1 = sqrt((rootPath[0][0] - self.rootPoint[0])**2 + (rootPath[0][1] - self.rootPoint[1])**2)
		dist2 = sqrt((rootPath[len(rootPath)-1][0] - self.rootPoint[0])**2 + (rootPath[len(rootPath)-1][1] - self.rootPoint[1])**2)
		
		if dist1 < dist2:
			self.rootPoint = rootPath[0]
		else:
			self.rootPoint = rootPath[len(path)-1]

		
		""" don't care about return values, are stored as object member variable """
		self.computeAllSplices2(plotIter = True)

		self.controlCurves = {}
		for pathID in pathIDs:
			if self.pathClasses[pathID]["parentID"] != None:
				parentPathID = parentPathIDs[pathID]
				parentSkeleton = self.localSkeletons[parentPathID]
				parentTerm1 = self.pathClasses[pathID]["controlTerm1_P"]
				parentTerm2 = self.pathClasses[pathID]["controlTerm2_P"]

				parentFrame = Pose(finalPoses[parentPathID])
				parentTerm1_G = parentFrame.convertLocalToGlobal(parentTerm1)
				parentTerm2_G = parentFrame.convertLocalToGlobal(parentTerm2)

				print "controlCurve:", pathID, parentPathID, parentTerm1_G, parentTerm2_G

				self.controlCurves[pathID] = getSkeletonPath(parentSkeleton, parentTerm1, parentTerm2)
			else:
				self.controlCurves[pathID] = self.localPaths[pathID]


		self.isChanged = False


	@logFunction
	def delNode(self, nodeID, pathID):
		print "deleting node", nodeID, "from path", pathID
		self.pathClasses[pathID]["nodeSet"].remove(nodeID)
		del self.pathClasses[pathID]["localNodePoses"][nodeID]
		
		self.isChanged = True		

	@logFunction
	#def moveNode(self, nodeID, splicePathIDs):
	def moveNode(self, nodeID, targetPathID):
		" change the shoot membership based on fitted splice "

		#particleDist = self.poseParticles["snapshots2"][0]

		" determine which paths are leaves "
		pathIDs = self.getPathIDs()

		" remove node from other shoots first"
		for pathID in pathIDs:
			if nodeID in self.pathClasses[pathID]["nodeSet"]:
				print "removing node", nodeID, "from path", pathID
				self.pathClasses[pathID]["nodeSet"].remove(nodeID)
				del self.pathClasses[pathID]["localNodePoses"][nodeID]

				#for p in particleDist:
				#	p.delNode(nodeID, pathID)



		print "adding node", nodeID, "to path", targetPathID
		self.pathClasses[targetPathID]["nodeSet"].append(nodeID)
		#self.addNode(nodeID,pathID)

		""" convert the controlPose to coordinates local to the parent frame """
		parentPathIDs = self.getParentHash()

		""" control state of maximum likelihood particle determine's location of shoot frame """
		controlPoses = self.getControlPoses()
		globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)
		shootControlPose_G = globalControlPoses_G[targetPathID]
		currFrame = Pose(shootControlPose_G)

		""" convert from global to local """
		nodePose_G = self.getNodePose(nodeID)
		nodePose_C = currFrame.convertGlobalPoseToLocal(nodePose_G)
		self.pathClasses[targetPathID]["localNodePoses"][nodeID] = nodePose_C

		#for p in particleDist:
		#	p.addNode(nodeID, targetPathID, nodePose_C)

		print "shoot", targetPathID, "has nodes", self.pathClasses[targetPathID]["nodeSet"]

		self.updateLandmarks()
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

		#updateCount = self.poseParticles["updateCount"] 
		#particleDist = self.poseParticles["snapshots2"][0]

		#for p in particleDist:
		#	nodePose_C = self.pathClasses[pathID]["localNodePoses"][nodeID]
		#	p.addNode(nodeID, pathID, nodePose_C)

		print pathID, "has nodes", self.pathClasses[pathID]["nodeSet"]


		self.updateLandmarks()
		self.isChanged = True


	@logFunction
	def mergePath(self, pathID, targetPathID = None):



		
		parentPathIDs = self.getParentHash()

		""" target shoot we want to merge to """
		if targetPathID == None:
			mergeTargetPathID = parentPathIDs[pathID]
		else:
			mergeTargetPathID = targetPathID

		print "merging", pathID, "to", mergeTargetPathID

		""" if the root shoot, return after doing nothing """
		if mergeTargetPathID == None:
			return

		""" reparent the direct children of the merged shoot """
		childPathIDs = self.getChildHash()
		tobeReparented = childPathIDs[pathID]

		nodeSet = deepcopy(self.pathClasses[pathID]["nodeSet"])

		for nodeID in nodeSet:
			self.moveNode(nodeID, mergeTargetPathID)

		self.delPath(pathID, mergeTargetPathID)

		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:
			if len(self.pathClasses[pathID]["nodeSet"]) == 0:
				self.delPath(pathID)

		#self.pathTermsData

		#for termUUID, termData in self.pathTermsData.iteritems():
		#	termShootID = termData["frameID"]

		#controlPoses = self.getControlPoses()
		#globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)
	

	@logFunction
	def delPath(self, pathID, mergeTargetID = None):

		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

		#particleDist2 = self.poseParticles["snapshots2"][0]

		""" default to parent merge """
		if mergeTargetID == None:
			mergeTargetID = parentPathIDs[pathID]


		""" book-keeping for branch tuple indexing """
		branchPathIDs = deepcopy(self.getPathIDs())
		branchPathIDs.remove(0)
		branchPathIDs.sort()

		partSetControlPoses_G = []

		print "deleting path", pathID, "merging to path", mergeTargetID
		try: 
			del self.pathClasses[pathID]
			#del self.pathTermsVisited[pathID]
			del self.localLeaf2LeafPathJunctions[pathID]
			
			#for part in particleDist2:
			for part in self.stepResults:

				""" temporarily disable this.  we don't really use the sample data """
				if False:
					partControlPoses = part["controlPoses_P"]
					partControlPoses_G = computeGlobalControlPoses(partControlPoses, parentPathIDs)

					partSetControlPoses_G.append(partControlPoses_G)

					del part["controlPoses_P"][pathID]
					del part["branchPoses_L"][pathID]
					del part["branchPoses_G"][pathID]
					del part["branchArcDists"][pathID]
					del part["branchControls"][pathID] 

					branchTupleIndex = part["maxLikelihoodBranch"]

					newTupleIndex = []
					for k in range(len(branchPathIDs)):
						indexID = branchPathIDs[k]
						if indexID != pathID:
							newTupleIndex.append(branchTupleIndex[k])

					newTupleIndex = tuple(newTupleIndex)

					if len(newTupleIndex) > 0:
						part["maxLikelihoodBranch"] = newTupleIndex
					else:
						part["maxLikelihoodBranch"] = ()

					print "oldTuple, newTuple:", branchTupleIndex, newTupleIndex
				

			reparentedPaths = []
			for childPathID, pathClass in self.pathClasses.iteritems():
				
				if pathClass["parentID"] == pathID:

					pathClass["parentID"] = mergeTargetID

					oldControlPose_G = globalControlPoses_G[childPathID]

					oldParentControlPose_G = globalControlPoses_G[pathID]
					oldParentFrame = Pose(oldParentControlPose_G)

					parentControlPose_G = globalControlPoses_G[mergeTargetID]
					newParentFrame = Pose(parentControlPose_G)

					""" compute terminals to control curve in new parent frame """
					oldParentTerm1 = self.pathClasses[childPathID]["controlTerm1_P"]
					oldParentTerm2 = self.pathClasses[childPathID]["controlTerm2_P"]

					oldParentTerm1_G = oldParentFrame.convertLocalToGlobal(oldParentTerm1)
					oldParentTerm2_G = oldParentFrame.convertLocalToGlobal(oldParentTerm2)

					newParentTerm1_P = newParentFrame.convertGlobalToLocal(oldParentTerm1_G)
					newParentTerm2_P = newParentFrame.convertGlobalToLocal(oldParentTerm2_G)

					self.pathClasses[childPathID]["controlTerm1_P"] = newParentTerm1_P
					self.pathClasses[childPathID]["controlTerm2_P"] = newParentTerm2_P


					""" compute the control pose relative to the parent """
					controlPose_P = newParentFrame.convertGlobalPoseToLocal(oldControlPose_G)
					pathClass["controlPose"] = controlPose_P

					currFrame = Pose(oldControlPose_G)

					#junctionPose_C = currFrame.convertGlobalPoseToLocal(globalJunctionPose_G)
					#pathClass["localJunctionPose"] = junctionPose_C 

					#tipPoint_C = currFrame.convertGlobalToLocal(tipPoint_G)
					#pathClass["tipPoint_L"] = tipPoint_C

					print "childPathID parent:", childPathID, pathID, mergeTargetID, pathClass["parentID"], self.pathClasses[childPathID]["parentID"]

					for nodeID in pathClass["nodeSet"]:

						nodePose_G = self.getNodePose(nodeID)

						#self.setNodePose(nodeID, nodePose_G)


					print "reparented pathClass:", childPathID, pathClass

					reparentedPaths.append(childPathID)

		
			for childPathID in reparentedPaths:

				pathClass = self.pathClasses[childPathID]

				for nodeID in pathClass["nodeSet"]:

					nodePose_G = self.getNodePose(nodeID)

					self.setNodePose(nodeID, nodePose_G)



		except:

			print "FAIL to delete path!"
			raise
			pass
	
		#part.addPath(newPathID, controlParentID, branchNodeID, nodePose_C, localDivergencePose_R, modJunctionPose_G, localJunctionPose_C, localParticleControlPose_P, self.NUM_BRANCHES, arcDists, controlPoses_P)

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

		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

		localPathSegsByID = {}
		for pathID in allPathIDs:
			localPathSegs = self.localLeaf2LeafPathJunctions[pathID]["localSegments"]
			localPathSegsByID[pathID] = localPathSegs


		#self.localPaths[pathID] = localResults[1]
		#localPathSegsByID[newPathID] 


		#self.localSkeletons = {}
		controlPose_G, controlParentID, tipPoint_G, branchPose_G, controlTerm1_P, controlTerm2_P = getInitSkeletonBranchPoint(globalJunctionPose_G, newPathID, globalMedial0_G, parentPathIDs, localPathSegsByID, self.localPaths, self.localTerms, self.localSkeletons, globalControlPoses_G, plotIter = True, hypothesisID = self.hypothesisID, nodeID = branchNodeID)

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
		tipPoint_C = currFrame.convertGlobalToLocal(tipPoint_G)

		""" basic shoot data structure """
		self.pathClasses[newPathID] = {"parentID" : controlParentID,
						"branchNodeID" : branchNodeID,
						"localDivergencePose" : localDivergencePose_R, 
						"localJunctionPose" : junctionPose_C, 
						"tipPoint_L" : tipPoint_C, 
						"sameProb" : {},
						"nodeSet" : [branchNodeID,],
						"localNodePoses" : {branchNodeID : nodePose_C},
						"globalJunctionPose" : newGlobJuncPose_G,
						"controlPose" : controlPose_P,
						"controlTerm1_P" : controlTerm1_P, 		
						"controlTerm2_P" : controlTerm2_P }		

		print "new pathClass:", self.pathClasses[newPathID]


		#self.pathTermsVisited[newPathID] = False
		
		self.pathIDs += 1
		self.shootIDs += 1

		#self.controlCurves = {}
		parentSkeleton = self.localSkeletons[controlParentID]
		parentTerm1 = self.pathClasses[newPathID]["controlTerm1_P"]
		parentTerm2 = self.pathClasses[newPathID]["controlTerm2_P"]
		self.controlCurves[newPathID] = getSkeletonPath(parentSkeleton, parentTerm1, parentTerm2)



		#updateCount = self.poseParticles["updateCount"] 
		#particleDist2 = self.poseParticles["snapshots2"][0]

		""" FIXME:  Full path spline instead of just relevant section.  Does this always work? """
		#parentShoot_G = self.paths[parentID]
		#parentShoot_L = self.localPaths[controlParentID]
		parentShoot_L = self.controlCurves[newPathID]
		parentShoot_G = []
		for p in parentShoot_L:
			parentShoot_G.append(parentFrame.convertLocalToGlobal(p))

		pathSpline_G = SplineFit(parentShoot_G)


		print "max control poses:", controlPose_G
		""" add new path to particles """
		#for part in particleDist2:
		#for k in range(len(particleDist2)):
		for partIndex in range(len(self.stepResults)):

			#part = particleDist2[k]
			partDict = deepcopy(self.stepResults[partIndex])



			nodePose_G = self.getNodePose(branchNodeID)
			rawPose_G = self.getNodeRawPose(branchNodeID)

			""" get local transform of GPAC to raw pose """
			nodeFrame = Pose(nodePose_G)
			transform_N_to_R = nodeFrame.convertGlobalPoseToLocal(rawPose_G)
			
			""" go back and convert this from GPAC pose to raw pose """
			particlePose_G = deepcopy(partDict["newPose0"])
			particlePose_G[2] = nodePose_G[2]
			particlePoseFrame = Pose(particlePose_G)
			particleRawFrame_G = particlePoseFrame.convertLocalOffsetToGlobal(transform_N_to_R)

			""" location of junction point from raw to global in the particle pose """
			rawParticlePoseFrame = Pose(particleRawFrame_G)
			particleJunctionPose_G = rawParticlePoseFrame.convertLocalOffsetToGlobal(localDivergencePose_R)


			""" compute the particle control pose based on the global location of the particle pose """
			controlPose_L = nodeFrame.convertGlobalPoseToLocal(controlPose_G)
			particleControlPose_G = particlePoseFrame.convertLocalOffsetToGlobal(controlPose_L)

			#print "control poses:", k, nodePose_G, particlePose_G, controlPose_G, particleControlPose_G
			#print "control poses:", partIndex, particleControlPose_G
		
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
			#resultDict = {}
			#resultDict["newPose0"] = newPose0
			#resultDict["newPose1"] = newPose1
			#resultDict["controlPoses_P"] = controlPoses
			#resultDict["displaceProb0"] = displaceProb0
			#resultDict["displaceProb1"] = displaceProb1
			#resultDict["branchPoses_G"] = branchPoses_G
			#resultDict["branchPoses_L"] = branchPoses_L
			#resultDict["landmarkSum"] = 0.0
			#resultDict["maxLikelihoodBranch"] = None
			#resultDict["branchArcDists"] = ()
			#resultDict["branchControls"] = ()

			partDict["controlPoses_P"][newPathID] = localParticleControlPose_P
			#print "addParticle control:", partIndex, localParticleControlPose_P
			partDict["branchPoses_L"][newPathID] = localJunctionPose_C
			partDict["branchPoses_G"][newPathID] = modJunctionPose_G
			partDict["branchArcDists"][newPathID] = arcDists
			partDict["branchControls"][newPathID] = controlPoses_P
			self.stepResults[partIndex] = partDict


			#part.addPath(newPathID, controlParentID, branchNodeID, nodePose_C, localDivergencePose_R, modJunctionPose_G, localJunctionPose_C, localParticleControlPose_P, self.NUM_BRANCHES, arcDists, controlPoses_P)

		#print "newPath", newPathID, "=", self.pathClasses[newPathID]
		#for k in range(len(self.stepResults)):
		#	part = self.stepResults[k]
		#	print "freshControl:", k, part["controlPoses_P"]

		self.updateLandmarks()

		self.isChanged = True

		return newPathID
 

	@logFunction
	def pathTermVisited(self, termID):
		
		try:
			self.pathTermsVisited[termID]
		except:
			pass
		else:
			self.pathTermsVisited[termID] = True

			self.pathTermsData[termID]["isVisited"] = True

		print "pathTermVisited(", termID, self.pathTermsVisited

	@logFunction
	def getPathTermsVisited(self):
		return self.pathTermsVisited

	@logFunction
	def resetTerms(self):
		
		for k, v in self.pathTermsVisited.iteritems():
			self.pathTermsVisited[k] = False
			self.pathTermsData[k]["isVisited"] = False

	@logFunction
	def getPathTerms(self):
		return self.pathTerms

	@logFunction
	def checkUniqueBranch2(self, parentPathID, nodeID1, depAngle, depPoint):


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
	def determineBranchPair2(self, nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, dirFlag, isUnique1, isUnique2, shootIDs):
		
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

		print self.hypothesisID, "determineBranchPair(", nodeID1, ",", nodeID2, ",", frontExist1, ",", frontExist2, ",", frontInterior1, ",", frontInterior2, ",", depAngle1, ",", depAngle2, ",", depPoint1, ",", depPoint2, ",", dirFlag, ",", isUnique1, ",", isUnique2, ",", ",", shootIDs, ")"
		
		self.shootIDs = shootIDs

		# depAngle1, depAngle2, frontInterior1, frontExist1, frontInterior2, frontExist2, depPoint1, depPoint2 , isUnique1, isUnique2


		foreTerm1 = frontInterior1 and frontExist1
		foreTerm2 = frontInterior2 and frontExist2

		frontAngDiff = diffAngle(depAngle1, depAngle2)
		
		" 60 degree threshold "
		ANG_THRESH = 1.047

		pathBranchIDs = [None, None]


		if foreTerm1 and foreTerm2:

			" both departing case: check if same or different "
			if fabs(frontAngDiff) < ANG_THRESH:

				" are on top of each other and are branching to the same path "
				if isUnique1 and isUnique2:
					
					branchNodeID = nodeID1
					globalJunctionPoint = depPoint1
					depAng = depAngle1
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					
					newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID

					print "foreA"

				if isUnique1 and not isUnique2:

					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[0] = newPathID

				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
					
						pathBranchIDs[1] = newPathID
				
			else:
				" these are branching to separate new paths "

				if isUnique1 and isUnique2:

					if dirFlag == 0:
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID1 = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[0] = newPathID1

					if dirFlag == 1:
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID2 = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[1] = newPathID2

					print "foreB"

				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[0] = newPathID

				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[1] = newPathID
		
		elif foreTerm1 and not foreTerm2:

			if isUnique1:
				   
				if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
					" both have same departure "
					branchNodeID = nodeID1
					globalJunctionPoint = depPoint1
					depAng = depAngle1
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID

					print "foreC"
					
				else:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[0] = newPathID

					print "foreD"
					
		elif foreTerm2 and not foreTerm1:

			if isUnique2:
				if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

					" both have same departure "
					branchNodeID = nodeID2
					globalJunctionPoint = depPoint2
					depAng = depAngle2
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
					print "foreE"

					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID

				else:
					
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

						pathBranchIDs[1] = newPathID

					print "foreF"

		
		#return isBranch, pathBranchIDs, isNew, self.shootIDs

		return self.shootIDs, pathBranchIDs
		

	@logFunction
	def getNearestPathPoint(self, originPoint):

		" find closest point out of all the global segments "

		pathIDs = self.getPathIDs()
		minP = None
		minDist = 1e100
		for pathID in pathIDs:
			segs1 = self.globalSegments[pathID]

			for seg in segs1:
				for p in seg:
					dist = sqrt((originPoint[0]-p[0])**2+(originPoint[1]-p[1])**2)

					if dist < minDist:
						minDist = dist
						minP = p
	
		return minP


	@logFunction
	def computeNavigationPath(self, startPose, endPose):
		
		" get the trimmed paths "
		
		" select shortest one "
		
		" get origin and target path "
		
		" get splice of origin and target path "
		" splice minJ1 and minJ2 "
		#orderPathIDs = self.getPathPath(minI1, minI2)		 
		
		#splicedPaths = self.splicePathIDs(orderPathIDs)
		
		
		#allSplices, terminals, junctions = self.getAllSplices()
		allSplices = self.getAllSplices2()

		#splicedPaths = allSplices

		splicedPaths = []
		#for k, result in allSplices.iteritems():
		for sPath in allSplices:
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


			

