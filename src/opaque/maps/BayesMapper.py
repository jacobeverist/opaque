
from LocalNode import getLongestPath, computeHullAxis
from StableCurve import StableCurve
from SplineFit import SplineFit
from Pose import Pose
from PoseData import PoseData
from MapProcess import movePath, addToPaths, localizePair, consistentFit, movePath, batchMovePath, batchLocalizePair, batchEval
from shoots import computeGlobalControlPoses
import gen_icp
from functions import *
from operator import itemgetter
import cProfile
import time
import traceback
import math
import cProfile

from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath
import pylab
import matplotlib.pyplot as plt
import graph

from MapState import MapState

import random

from guppy import hpy


def printStack():

	flist = traceback.format_stack()
	flist = flist[:-1]
	
	printStr = ""
	for line in flist:
		printStr += line
		
	print printStr


class BayesMapper:

	def __init__(self, walls):
		
		#self.probe = probe
		
		#self.walls = self.probe.getWalls()
		self.walls = walls

		self.poseData = PoseData()
		
		self.poseData.numNodes = 0


		""" initialize to a single map hypothesis """
		self.particleIDs = 0
		self.mapHyps = {}
		self.mapHyps[self.particleIDs] = MapState(self.poseData, self.particleIDs)
		self.particleIDs += 1
		#self.mapHyps[self.particleIDs] = self.mapHyps[0].copy(self.particleIDs)
		#self.mapHyps[self.particleIDs] = MapState(self.poseData, self.particleIDs)
		#self.particleIDs += 1

		self.shootIDs = 1


		self.activeHypID = None


		self.poseData.medialAxes = {}
		self.poseData.aHulls = {}
		self.poseData.isBowties = {}
		self.poseData.faceDirs = {}
		self.poseData.travelDirs = {}
		self.poseData.medialLongPaths = {}
		self.poseData.numLeafs = {}
		self.poseData.correctedPostures = {}
		self.poseData.isNodeFeatureless = {}
		self.poseData.frontProbeError = {}
		self.poseData.backProbeError = {}
		self.poseData.spatialFeatures = {}

	
		self.pathDrawCount = 0
		self.pathPlotCount = 0
		self.statePlotCount = 0
		self.trimCount = 0
		self.multiDepCount = 0
		self.profCount = 0


		self.tempCount = 0

		#self.overlapPlotCount = 0
		#self.pathDrawCount = 0
		#self.candidateCount = 0
		#self.commonOriginCount = 0

		self.colors = []
		self.colors.append([0, 0, 255])
		self.colors.append([0, 255, 0])
		self.colors.append([255, 0, 0])
		self.colors.append([255, 0, 255])
		self.colors.append([255, 255, 0])
		self.colors.append([0, 255, 255])
		"""
		self.colors.append([215, 48, 39])
		self.colors.append([252, 141, 89])
		self.colors.append([254, 224, 144])
		self.colors.append([224, 243, 248])
		self.colors.append([145, 191, 219])
		self.colors.append([69, 117, 180])
		"""

		for color in self.colors:
			color[0] = float(color[0])/256.0
			color[1] = float(color[1])/256.0
			color[2] = float(color[2])/256.0
		
		for i in range(1000):
			self.colors.append((random.random(),random.random(),random.random()))


		self.E_overlap = matrix([[ 0.1,  0.0, 0.0],
							[ 0.0,	0.01, 0.0],
							[0.0, 0.0,	0.02]])
		
		self.E_inplace = matrix([[ 0.05, 0.0, 0.0 ],
							[ 0.0, 0.05, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

	@property
	def numNodes(self):
		return self.poseData.numNodes

	@property
	def topHyp(self):
		hypIDs = self.mapHyps.keys()
		topID = hypIDs[0]
		return self.mapHyps[topID]

	@logFunction
	def loadNewNode(self, newNode):

		#self.currNode = newNode
		nodeID = self.poseData.numNodes
		#self.nodeHash[nodeID] = self.currNode

		print "incrementing numNodes"
		print self.poseData.numNodes
		self.poseData.numNodes += 1
		print self.poseData.numNodes
		
		hull1, medial1 = computeHullAxis(nodeID, newNode, tailCutOff = False)
		#hull1, medial1 = computeHullAxis(nodeID, newNode, tailCutOff = True)
		self.poseData.aHulls[nodeID] = hull1
		self.poseData.medialAxes[nodeID] = medial1
		self.poseData.numLeafs[nodeID] = newNode.getNumLeafs()
		self.poseData.faceDirs[nodeID] = newNode.faceDir
		self.poseData.isBowties[nodeID] = newNode.isBowtie			
		self.poseData.medialLongPaths[nodeID] = newNode.medialLongPaths
		self.poseData.correctedPostures[nodeID] = newNode.getStableGPACPosture()
		self.poseData.isNodeFeatureless[nodeID] = newNode.getIsFeatureless()
		self.poseData.frontProbeError[nodeID] = newNode.frontProbeError
		self.poseData.backProbeError[nodeID] = newNode.backProbeError
		self.poseData.travelDirs[nodeID] = newNode.travelDir
		self.poseData.spatialFeatures[nodeID] = newNode.spatialFeatures


		for mid, mapHyp in self.mapHyps.iteritems():
			mapHyp.updatePoseData(self.poseData)

		""" for each current map hypothesis, integrate the new node """
		currHyps = self.mapHyps
		for mid, mapHyp in currHyps.iteritems():

			print "loading", nodeID, "hyp", mapHyp.hypothesisID

			mapHyp.gndPoses[nodeID] = newNode.getGndGlobalGPACPose()
			mapHyp.gndRawPoses[nodeID] = newNode.getGndPose()
			#mapHyp.setNodePose(nodeID, newNode.getGlobalGPACPose())
			mapHyp.nodePoses[nodeID] = newNode.getGlobalGPACPose()

			" FIXME:  raw pose does not get updated yet with GPAC pose "

			#print "rawPoses:", mapHyp.nodeRawPoses[nodeID], newNode.getEstPose()

			mapHyp.nodeRawPoses[nodeID] = newNode.getEstPose()

			pathIDs = mapHyp.getPathIDs()
			for k in range(mapHyp.numPoseParticles):
				for pathID in pathIDs:

					if pathID != 0:
						particle = mapHyp.poseParticles["snapshots2"][0][k]

						controlPoseDist = particle.junctionData[pathID]["controlPoseDist"]
						branchPoseDist = particle.junctionData[pathID]["branchPoseDist"]

						print "hypID,pathID,particleIndex:", mapHyp.hypothesisID, pathID, k, controlPoseDist
					

		self.mapHyps = self.integrateNode(currHyps, nodeID)


	@logFunction
	def restoreNode(self, dirName, numNodes):
		
		print "loading" + dirName + "/stateSave_%04u.txt" % (numNodes)
		f = open(dirName + "/stateSave_%04u.txt" % (numNodes), 'r')		
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')

		nodeID = 'foo'

		existMapHyp = self.mapHyps
		existParticleIDs = self.particleIDs
		existWalls = self.walls

		self.travelDirs

		exec(saveStr)
		self.poseData.travelDirs = self.poseData.travelDirs[nodeID]
		self.particleIDs = existParticleIDs
		self.walls = existWalls

		self.mapHyps = existMapHyp
		nodeID = numNodes 
		self.poseData.numNodes = numNodes+1

		" get node poses from some map hypothesis "
		hid = mapHypIDs.values()[0]
		tempMapState = MapState(self.poseData, hid)
		tempMapState.restoreState(dirName, numNodes)

		for mid, mapHyp in self.mapHyps.iteritems():
			mapHyp.updatePoseData(self.poseData)

		""" for each current map hypothesis, integrate the new node """
		currHyps = self.mapHyps
		for mid, mapHyp in currHyps.iteritems():

			mapHyp.gndPoses[nodeID] = tempMapState.gndPoses[nodeID]
			mapHyp.gndRawPoses[nodeID] = tempMapState.gndRawPoses[nodeID]
			mapHyp.nodePoses[nodeID] = tempMapState.nodePoses[nodeID]
			#mapHyp.setNodePose(nodeID, newNode.getGlobalGPACPose())

			#print "rawPoses:", mapHyp.nodeRawPoses[nodeID], newNode.getEstPose()
			" FIXME:  raw pose does not get updated yet with GPAC pose "
			mapHyp.nodeRawPoses[nodeID] = tempMapState.nodeRawPoses[nodeID]

		self.mapHyps = self.integrateNode(currHyps, nodeID)

	@logFunction
	def integrateNode(self, hypSet, nodeID):

		" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
		#direction = newNode.travelDir
		direction = self.poseData.travelDirs[nodeID]

		" ensure the medial axes are computed before this check "
		#computeHullAxis(nodeID, newNode, tailCutOff = False)

		#currHyps = {}



		print "integrating node", nodeID

		if nodeID > 0:
			
			" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
			if self.poseData.numNodes >= 4:

				if nodeID % 2 == 1:
					#hypSet = batchMovePath(hypSet, nodeID, direction)


					for pID, currHyp in hypSet.iteritems():
						time1 = time.time()
						#movePath(currHyp, nodeID, direction, distEst = 1.0)
						currHyp.batchDisplaceParticles(nodeID-1, nodeID)


						"""
						particleDist2 = currHyp.poseParticles["snapshots2"][0]
						maxParticle = particleDist2[currHyp.currMaxIndex]

						currHyp.setNodePose(nodeID-1, maxParticle.pose0)
						currHyp.setNodePose(nodeID, maxParticle.pose1)
						"""

						#currHyp.updateMaxParticle(currHyp.currMaxIndex)


						currHyp.drawPoseParticles()
						time2 = time.time()
						print "TIME displace", currHyp.hypothesisID, "=", time2-time1 

				"""
				for pID, mapHyp in hypSet.iteritems():
					" Move node along path "
					#self.movePath(mapHyp, nodeID, direction)
					movePath(mapHyp, nodeID, direction)
					#self.drawConstraints(mapHyp, self.statePlotCount)
					#self.statePlotCount += 1
					#self.drawPathAndHull(mapHyp)
				"""
		
		nodeID1 = self.poseData.numNodes-2
		nodeID2 = self.poseData.numNodes-1

		print "nodeID1, nodeID2 =", nodeID1, nodeID2
		
		" CHECK FOR A BRANCHING EVENT "
		
		if self.poseData.numNodes >= 2 and self.poseData.numNodes % 2 == 0:

			for pID, mapHyp in hypSet.iteritems():

				print "entering node", nodeID, "from hypothesis", mapHyp.hypothesisID

				" DETECT BRANCHING EVENTS FOR THE 2 NODES OF LAST STOP "
				" AFTER ODD-NUMBER OVERLAP OF CURRENT STOP HAS IRONED OUT ERRORS "
				
				pathIDs = mapHyp.getPathIDs()

				print "hypothesis", mapHyp.hypothesisID, "paths =", pathIDs

				print "generating node", nodeID, "from hypothesis", mapHyp.hypothesisID

				" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
				#mapHyp.generatePaths()

				
				print "drawing node", nodeID, "from hypothesis", mapHyp.hypothesisID
				#self.drawConstraints(mapHyp, self.statePlotCount)
				#self.statePlotCount += 1
				#self.drawPathAndHull(mapHyp)


				mapHyp.isNodeBranching[nodeID1] = False
				mapHyp.isNodeBranching[nodeID2] = False


			#hypSet = self.addToPaths(hypSet, nodeID1, nodeID2)
			time1 = time.time()
			self.shootIDs, self.particleIDs, hypSet = addToPaths(self.shootIDs, self.particleIDs, hypSet, nodeID1, nodeID2)
			for pID, currHyp in hypSet.iteritems():
				currHyp.drawPoseParticles()
			time2 = time.time()
			print "TIME addToPaths =", time2-time1 

			#scurrHyps = self.addToPaths(mapHyp, nodeID1, nodeID2)

			#for pID, currHyp in hypSet.iteritems():
			#	localizePair(currHyp, nodeID1, nodeID2)

			#hypSet = batchLocalizePair(hypSet, nodeID1, nodeID2)
			#batchEval(hypSet)
			
			for pID, mapHyp in hypSet.iteritems():
				time1 = time.time()

				#self.propagateBranchStates(mapHyp, nodeID1)

				#cProfile.runctx("mapHyp.localizePoseParticles(nodeID1,nodeID2)", globals(), locals(), "localizeParticles_%d.prof" % self.profCount)	
				#self.profCount += 1
				
				mapHyp.localizePoseParticles(nodeID1, nodeID2)
				#mapHyp.drawPoseParticles()

				time2 = time.time()

				print "TIME localize", pID, "=", time2-time1

			for pID, currHyp in hypSet.iteritems():

				currHyp.computeEval()

				""" LOCALIZE NODE PAIR """
				#localizePair(currHyp, nodeID1, nodeID2)

				#self.propagateBranchStates(currHyp, nodeID1)

				#currHyp.generatePaths()
				self.drawPathAndHull2(currHyp)


				"""
				if nodeID1 >= 2:
					newPath = deepcopy(currHyp.paths[0])
					p0 = newPath[0]
					pN = newPath[-1]
					
					rootPose = [-2.0,0.0]
					
					dist0 = sqrt((rootPose[0]-p0[0])*(rootPose[0]-p0[0]) + (rootPose[1]-p0[1])*(rootPose[1]-p0[1]))
					distN = sqrt((rootPose[0]-pN[0])*(rootPose[0]-pN[0]) + (rootPose[1]-pN[1])*(rootPose[1]-pN[1]))
			
					if dist0 > distN:
						newPath.reverse()
								
					" new path "
					newSpline = SplineFit(newPath, smooth = 0.1)
		
		
					posNew = currHyp.nodePoses[nodeID1]
					posOld = currHyp.nodePoses[nodeID1-2]
					minDist, minU, closestPoint = newSpline.findClosestPoint(posNew)
					arcDistNew = newSpline.dist(0.0, minU)
		
					minDist, minU, closestPoint = newSpline.findClosestPoint(posOld)
					arcDistOld = newSpline.dist(0.0, minU)
		
					print "arcDistNew, arcDistOld, diff =", arcDistNew, arcDistOld, arcDistNew-arcDistOld
				"""
		
		""" remove defective maps """
		toDelete = []
		for pID, currHyp in hypSet.iteritems():
			print pID, "mapOverlapSum =", currHyp.mapOverlapSum, "isNoLocalize =", currHyp.isNoLocalize
			#if currHyp.mapOverlapSum > 0.0 or currHyp.isNoLocalize:
			#if currHyp.mapOverlapSum > 6 or currHyp.isNoLocalize:
			if currHyp.mapOverlapSum > 4.0 or currHyp.isNoLocalize:
				#if pID != self.activeHypID:
				#	toDelete.append(pID)
				toDelete.append(pID)
				if pID == self.activeHypID:
					self.activeHypID = None

				
		for pID in toDelete:
			del hypSet[pID]

		hp = hpy()
		print hp.heap()

		return hypSet


	@logFunction
	def propagateBranchStates(self, mapHyp, nodeID, numGuesses = 11, excludePathIDs = []):

		poseData = mapHyp.poseData

		allSplices, terminals, junctions = mapHyp.getAllSplices(plotIter = False)

		for pathID, pathDesc in mapHyp.pathClasses.iteritems():
			if pathID != 0:
				"""
				this is a branching path,
				now we need hypotheses about where it's branching from.

				pathDesc = {"parentID" : parentID,
						"branchNodeID" : branchNodeID,
						"localJunctionPose" : localJunctionPose,
						"sameProb" : {},
						"nodeSet" : [],
						"globalJunctionPose" : globalJunctionPose }		
				"""

				mapHyp.updateParticles(pathID)

		
				" draw env "
				" draw trimmed path "
				" draw points "

		" initial distribution of particles representing branching point, angle, and parent path "
		" do the particles first. "



		" now do the fitting of the path super-axis to the different splices "

		return


	@logFunction
	def saveState(self):

					
		saveFile = ""

		saveFile += "self.walls = " + repr(self.walls) + "\n"


		saveFile += "self.poseData.numNodes = " + repr(self.poseData.numNodes) + "\n"
		saveFile += "self.particleIDs = " + repr(self.particleIDs) + "\n"
		saveFile += "self.poseData.medialAxes = " + repr(self.poseData.medialAxes) + "\n"
		saveFile += "self.poseData.aHulls = " + repr(self.poseData.aHulls) + "\n"
		saveFile += "self.poseData.isBowties = " + repr(self.poseData.isBowties) + "\n"
		saveFile += "self.poseData.faceDirs = " + repr(self.poseData.faceDirs) + "\n"
		saveFile += "self.poseData.medialLongPaths = " + repr(self.poseData.medialLongPaths) + "\n"
		saveFile += "self.poseData.numLeafs = " + repr(self.poseData.numLeafs) + "\n"
		saveFile += "self.poseData.correctedPostures = " + repr(self.poseData.correctedPostures) + "\n"
		saveFile += "self.poseData.isNodeFeatureless = " + repr(self.poseData.isNodeFeatureless) + "\n"
		saveFile += "self.poseData.frontProbeError = " + repr(self.poseData.frontProbeError) + "\n"
		saveFile += "self.poseData.backProbeError = " + repr(self.poseData.backProbeError) + "\n"
		saveFile += "self.poseData.travelDirs = " + repr(self.poseData.travelDirs) + "\n"
		saveFile += "self.poseData.spatialFeatures = " + repr(self.poseData.spatialFeatures) + "\n"

		saveFile += "mapHypIDs = " + repr(self.mapHyps.keys()) + "\n"

		f = open("stateSave_%04u.txt" % (self.poseData.numNodes-1), 'w')
		f.write(saveFile)
		f.close()		



		" SAVE STATE "
		for k in self.mapHyps.keys():
			self.mapHyps[k].saveState(self.poseData.numNodes-1)

		
	@logFunction
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/stateSave_%04u.txt" % (numNodes)
		f = open(dirName + "/stateSave_%04u.txt" % (numNodes), 'r')		
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		self.poseData = PoseData()
		exec(saveStr)
		
		for hid in mapHypIDs:
			self.mapHyps[hid] = MapState(self.poseData, hid)
			self.mapHyps[hid].restoreState(dirName, numNodes)

		

	@logFunction
	def addNode(self, newNode):
		
		#self.currNode = newNode
		nodeID = self.poseData.numNodes
		self.nodeHash[nodeID] = newNode
		self.poseData.numNodes += 1


	@logFunction
	def selectSplice(self, mapHyp, nodeID1, nodeID2, medial1, medial2, estPose1, estPose2, orderedPathIDs1, orderedPathIDs2):
		
		splicedPaths1 = mapHyp.splicePathIDs(orderedPathIDs1)
		splicedPaths2 = mapHyp.splicePathIDs(orderedPathIDs2)

		
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medial1_vec = medialSpline1.getUniformSamples()
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		medial2_vec = medialSpline2.getUniformSamples()


		print "received", len(splicedPaths1), "spliced paths from path IDs", orderedPathIDs1
		print "received", len(splicedPaths2), "spliced paths from path IDs", orderedPathIDs2

		" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 |"
		" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 "

		poseOrigin1 = Pose(estPose1)
		globalMedial1 = []
		for p in medial1:
			globalMedial1.append(poseOrigin1.convertLocalToGlobal(p))

		results1 = []		
		for k in range(len(splicedPaths1)):
			path = splicedPaths1[k]			
			path = orientPath(path, globalMedial1)
			result = getMultiDeparturePoint(path, medial1_vec, estPose1, estPose1, orderedPathIDs1, nodeID1)
			results1.append(result+(k,))


		results1 = sorted(results1, key=itemgetter(14))
		results1 = sorted(results1, key=itemgetter(12), reverse=True)				

		poseOrigin2 = Pose(estPose2)
		globalMedial2 = []
		for p in medial2:
			globalMedial2.append(poseOrigin2.convertLocalToGlobal(p))

		results2 = []		
		for k in range(len(splicedPaths2)):
			path = splicedPaths2[k]			
			path = orientPath(path, globalMedial2)
			result = getMultiDeparturePoint(path, medial2_vec, estPose2, estPose2, orderedPathIDs2, nodeID2)
			results2.append(result+(k,))

		results2 = sorted(results2, key=itemgetter(14))
		results2 = sorted(results2, key=itemgetter(12), reverse=True)				


		""" NOTE:  We only discriminated based on node k but not node k+1 """
		kIndex = results1[0][15]

		return splicedPaths1[kIndex]

	@logFunction
	def mergeSiblings(self, mapHyp, resultOffset1, resultOffset2, pathID1, pathID2, lastCost, matchCount):

		" verify that this path still exists before we try to merge it "
		allPathIDs = mapHyp.getPathIDs()
		if pathID1 in allPathIDs and pathID2 in allPathIDs:
	
			parentID1 = mapHyp.pathClasses[pathID1]["parentID"] 
			parentID2 = mapHyp.pathClasses[pathID2]["parentID"] 
	
			if parentID1 != parentID2:
				print "ERROR mergeSiblings:", pathID1, pathID2, parentID1, parentID2, "not siblings"
				return
	
			
			" capture transform between pathID1 and pathID2 under correction "
	
			mergeNodes1 = copy(mapHyp.getNodes(pathID1))
			mergeNodes2 = copy(mapHyp.getNodes(pathID2))
			
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
				orderedPathIDs = mapHyp.getOrderedOverlappingPaths(nodeID)
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
			#self.drawConstraints(mapHyp, self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull(mapHyp)
			
			" offset between paths "
			for nodeID in mergeNodes1:
	
				if not hasMoved[nodeID]:
					print "modifying node", nodeID, "in path", pathID1 
					" adjust node poses to relative path positions "
					nodePose = mapHyp.getNodePose(nodeID)
					guessPose = poseOrigin1.convertLocalOffsetToGlobal(nodePose)
					#self.nodeHash[nodeID].setGPACPose(guessPose)
					mapHyp.setNodePose(nodeID, guessPose)
					#mapHyp.nodePoses[nodeID] = guessPose
					
					#self.drawConstraints(mapHyp, self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull(mapHyp)
					
					hasMoved[nodeID] = True

			" offset between paths "
			for nodeID in mergeNodes2:
	
				if not hasMoved[nodeID]:
					print "modifying node", nodeID, "in path", pathID2 
					" adjust node poses to relative path positions "
					nodePose = mapHyp.getNodePose(nodeID)
					guessPose = poseOrigin2.convertLocalOffsetToGlobal(nodePose)
					#self.nodeHash[nodeID].setGPACPose(guessPose)
					mapHyp.setNodePose(nodeID, guessPose)
					#mapHyp.nodePoses[nodeID] = guessPose
					
					#self.drawConstraints(mapHyp, self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull(mapHyp)

					hasMoved[nodeID] = True

			mapHyp.generatePaths()

			#self.drawConstraints(mapHyp, self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull(mapHyp)
			#self.drawTrimmedPaths(mapHyp)


			for nodeID in mergeNodes2:
	
				#hull, medial = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)
				hull = self.poseData.aHulls[nodeID]
				medial = self.poseData.medialAxes[nodeID]
				estPose = mapHyp.getNodePose(nodeID)
	
				orderedPathIDs = [parentID1, pathID1]
				splicedPaths1 = mapHyp.splicePathIDs(orderedPathIDs)	
				print "received", len(splicedPaths1), "spliced paths from path IDs", orderedPathIDs
				" departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum, angDiff2 |"
		

				" exclusion not required since paths are siblings "
				try:
					consistentFit(mapHyp, nodeID, estPose)
				except IndexError:
					print "failed to consistentFit node", nodeID
					pass
				#self.consistentFit(nodeID, estPose, excludePathIDs = [pathID2])
		

				#self.drawConstraints(mapHyp, self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull(mapHyp)
				

			for nodeID in mergeNodes2:

				print "adding node", nodeID, "to path", pathID1
				mapHyp.addNode(nodeID, pathID1)
				print "deleting node", nodeID, "from path", pathID2
				mapHyp.delNode(nodeID, pathID2)

			print "A: mapHyp.pathClasses[" + str(pathID1) + "] = " + repr(mapHyp.pathClasses[pathID1])
			print "A: mapHyp.pathClasses[" + str(pathID2) + "] = " + repr(mapHyp.pathClasses[pathID2])


			" move nodes to pathID1, delete pathID2 "
			mapHyp.delPath(pathID2, pathID1)

			print "B: mapHyp.pathClasses[" + str(pathID1) + "] = " + repr(mapHyp.pathClasses[pathID1])
			
			mapHyp.generatePaths()

			print "C: mapHyp.pathClasses[" + str(pathID1) + "] = " + repr(mapHyp.pathClasses[pathID1])

			#self.drawConstraints(mapHyp, self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull(mapHyp)
			#self.drawTrimmedPaths(mapHyp)


			""" Now we check to see if any branch events have occurred from this merge """
			" departure events for node 1 "
			for nodeID1 in mergeNodes2:
				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
				contig1 = []		

				orderedPathIDs1 = mapHyp.getOrderedOverlappingPaths(nodeID1)
				#orderedPathIDs1.insert(0,0)
				#orderedPathIDs1.append(0)

				#orderedPathIDs1.append(0)
				print nodeID1, "orderedPathIDs1:", orderedPathIDs1

				" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = mapHyp.getDeparturePoint(mapHyp.trimmedPaths[pathID], nodeID1, plotIter = False)
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
				
				isFront1 = self.poseData.faceDirs[nodeID1]
				if isFront1:
					dirFlag = 0
				elif not isFront1:
					dirFlag = 1
				else:
					print isFront1
					raise

				isBranch, pathBranchID, isNew = mapHyp.determineBranchSingle(nodeID1, frontExist1, frontInterior1, depAngle1, depPoint1, parentPathID1, dirFlag)
				
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

				isFront1 = self.poseData.faceDirs[nodeID1]
				if isFront1:
					dirFlag = 1
				elif not isFront1:
					dirFlag = 0
				else:
					print isFront1, isFront2
					raise
				
				isBranch, pathBranchID, isNew = mapHyp.determineBranchSingle(nodeID1, backExist1, backInterior1, depAngle1, depPoint1, parentPathID1, dirFlag)
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
				pathIDs = mapHyp.getPathIDs()
				isAParent = {}
				for k in pathIDs:
					isAParent[k] = False
				for k in orderedPathIDs1:
					print "index:", k
					currPath = mapHyp.getPath(k)
					currParent = currPath["parentID"]
					if currParent != None:
						isAParent[currParent] = True
				
				"add nodes to paths that are the leaves "
				for pathID in orderedPathIDs1:
					if not isAParent[pathID]:				
						mapHyp.addNode(nodeID1,pathID)

				mapHyp.generatePaths()
				mapHyp.trimPaths(mapHyp.paths)						

				#self.drawConstraints(mapHyp, self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull(mapHyp)
				#self.drawTrimmedPaths(mapHyp)

			#for k, v in hasMoved.iteritems():
			#	 pass
			
	

	@logFunction
	def mergePaths(self, mapHyp):


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

		

		toBeMerged = mapHyp.comparePaths()
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
			
			allPathIDs = mapHyp.getPathIDs()
			if pathID1 in allPathIDs and pathID2 in allPathIDs:
				
				self.mergeSiblings(mapHyp, offset, offset2, pathID1, pathID2, mergeThis[5], mergeThis[6])			
			
		for mergeThis in toBeMerged:		

			print "mergeThis:", mergeThis[0], mergeThis[1], mergeThis[2]


			pathID1 = mergeThis[0]
			pathID2 = mergeThis[1]
			offset = mergeThis[2]
			pathSplices = mergeThis[3]
			
			" verify that this path still exists before we try to merge it "
			allPathIDs = mapHyp.getPathIDs()
			if pathID1 in allPathIDs and pathID2 in allPathIDs:

				
				" capture transform between pathID1 and pathID2 under correction "

				cOffset = offset
				rootID = pathID1
				lessID = pathID2
				
				
				mergeNodes = copy(mapHyp.getNodes(lessID))
				
				
				" capture transform between pathID2 and member nodes "
				" add compute new node locations under corrected transform "
				#self.drawConstraints(mapHyp, self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull(mapHyp)
				
				nodeHasMerged = {}
				
				" offset between paths "
				poseOrigin = Pose(cOffset)
				for nodeID in mergeNodes:

					" adjust node poses to relative path positions "
					nodePose = mapHyp.getNodePose(nodeID)
					guessPose = poseOrigin.convertLocalOffsetToGlobal(nodePose)
					#self.nodeHash[nodeID].setGPACPose(guessPose)
					mapHyp.setNodePose(nodeID, guessPose)
					#mapHyp.nodePoses[nodeID] = guessPose
					
					#self.drawConstraints(mapHyp, self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull(mapHyp)
					
					pathNodes = mapHyp.getNodes(rootID)
					targetNodeID = nodeID
					
					print "add path constraints:"
					print "pathNodes:", pathNodes
					print "targetNodeID:", targetNodeID
					
					nodeHasMerged[nodeID] = False
			

				

					" control point nearest the GPAC origin "
					globalPoint1 = mapHyp.getNodePose(nodeID)[:2]
					print "adjusting pose", nodeID, "from", mapHyp.getNodePose(nodeID)
					
					" if siblings, include mutual parent in splice"

					try:
						consistentFit(mapHyp, nodeID, mapHyp.getNodePose(nodeID), excludePathIDs = [lessID])
					except IndexError:
						print "failed to consistentFit node", nodeID
						pass

					#self.drawConstraints(mapHyp, self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull(mapHyp)

					
					
					mapHyp.addNode(nodeID, rootID)
					mapHyp.delNode(nodeID, lessID)
					
					if len(mapHyp.getNodes(lessID)) == 0:
						mapHyp.delPath(lessID, rootID)

					
					
					mapHyp.generatePaths()
					

					#self.drawConstraints(mapHyp, self.statePlotCount)
					self.statePlotCount += 1
					self.drawPathAndHull(mapHyp)


				print "relaxed constraints",  self.statePlotCount
				#self.drawConstraints(mapHyp, self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull(mapHyp)

				
				" constrain nodes to merge to pathID1 "
				
				" move nodes to pathID1, delete pathID2 "
				mapHyp.delPath(lessID, rootID)
	
				mapHyp.generatePaths()
				
				#self.drawConstraints(mapHyp, self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull(mapHyp)
			

	@logFunction
	def getNearestPathPoint(self, originPoint):

		candPoints = {}
		minDist = 1e100
		minHypID = -1 

		for hypID, mapHyp in self.mapHyps.iteritems():
			newPoint = mapHyp.getNearestPathPoint(originPoint)
			candPoints[hypID] = newPoint

			dist = sqrt((newPoint[0]-originPoint[0])**2 + (newPoint[1]-originPoint[1])**2)

			if dist < minDist:
				minDist = dist
				minHypID = hypID

		return candPoints[minHypID]

	@logFunction
	def computeNavigationPath(self, startPose, endPose):
		return self.mapHyps[self.activeHypID].computeNavigationPath(startPose, endPose)

	@logFunction
	def plotEnv(self, axes=0):
		
		#walls = self.probe.getWalls()
		walls = self.walls
	
		"""
		for wall in walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = copy(wall[i])
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g', zorder=0)
		"""
			

		if axes == 0:
		
			for wall in walls:
				xP = []
				yP = []
				for i in range(len(wall)):
					p = copy(wall[i])
					xP.append(p[0])
					yP.append(p[1])
		
				pylab.plot(xP,yP, linewidth=2, color = 'g', zorder=0)

		else:

			for wall in walls:
				xP = []
				yP = []
				for i in range(len(wall)):
					p = copy(wall[i])
					xP.append(p[0])
					yP.append(p[1])
		
				axes.plot(xP,yP, linewidth=2, color = 'g',zorder=0)
			

	@logFunction
	def drawPathAndHull2(self, mapHyp):
		
		" 1) plot the pose of local splines and postures "
		" 2) plot the alpha shape of the union of pose alpha shapes and medial axis tree "
		" 3) plot the alpha shape of the union of pose alpha shapes and the long medial axis "
		" 4) plot the trimmed path "

		#fig = plt.figure()
		fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
		#fig, (ax1, ax2) = plt.subplots(2,  sharex=True, sharey=True)
		fig.set_size_inches(16,12)
		fig.tight_layout(pad=4.0)

		#fig = plt.figure()
		#ax1 = fig.add_subplot(2,1,1)
		#ax1.plot(range(10), 'b-')

		#ax2 = fig.add_subplot(2,1,2)
		#ax2.plot(range(20), 'r^')

		#ax1.plot(x)
		#ax1.set_xlim(-4, 4)                    
		#ax1.set_ylim(-3, 3)
		ax1.set_xlim(-10, 10)
		ax1.set_aspect("equal")

		#ax3 = ax2
		#ax4 = ax3

		#fig.savefig("foo.png")

		poses = []
		for i in range(self.poseData.numNodes):
			poses.append(mapHyp.getNodePose(i))

		#pylab.clf()
		for i in range(self.poseData.numNodes):

			hull = self.poseData.aHulls[i]
			currPose = mapHyp.getNodePose(i)

			currProfile = Pose(currPose)
			posture1 = self.poseData.correctedPostures[i]

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
			#ax1.plot(xP,yP, color=(255/256.,182/256.,193/256.))	
			#ax1.plot(xP,yP, color=(112/256.,147/256.,219/256.))	
			ax1.plot(xP,yP, color='b')	

			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			#ax1.plot(xP,yP, color=(112/256.,147/256.,219/256.))	
			#ax2.plot(xP,yP, color=(112/256.,147/256.,219/256.))	
			ax2.plot(xP,yP, color='b')	


		xP = []
		yP = []
		for pose in poses:
			xP.append(pose[0])
			yP.append(pose[1])		
		ax1.scatter(xP,yP, color='k', linewidth=1, zorder=10)

		#self.plotEnv()

		ax1.set_title("%d Poses" % self.poseData.numNodes)
		ax2.set_title("%d Poses" % self.poseData.numNodes)
		#pylab.xlim(-10,10)
		#pylab.ylim(-10,10)


		" plot 3" 
			
		pathIDs = mapHyp.getPathIDs()

		allNodes = []
		for k in pathIDs:
			nodeSet = mapHyp.getNodes(k)
			allNodes += copy(nodeSet)
		
		allNodes.sort() 
		
		if len(allNodes) > 0:
			highestNodeID = allNodes[-1]
		else:
			highestNodeID = 1e100		
			


		"""
		for path in mapHyp.theoryMedialLongPaths:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color='k')

			globJuncPose = self.getGlobalJunctionPose(pathID)
			if globJuncPose != None:
				pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='k')

		"""

		"""
		
		xP = []
		yP = []
		for p in vertices:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='r')
		
		sizes = []
		for path in longPaths:
			sizes.append(len(path))
		
		bufStr1 = ""
		for dist in juncDists:
			if dist != None:
				bufStr1 += "%1.2f " % dist


		bufStr2 = ""
		for juncAngs in juncAngSet:
			for angs in juncAngs:
				if angs != None:
					bufStr2 += "%1.2f %1.2f " % (angs[0],angs[1])
		
		pylab.axis("equal")
		pylab.title("Path %d %s %s %s" % (pathID, sizes,bufStr1,bufStr2))
		pylab.savefig("medialOut2_%04u.png" % self.topCount)
		print "saving medialOut2_%04u.png" % self.topCount
		self.topCount += 1
		"""

		""" convert the controlPose to coordinates local to the parent frame """
		parentPathIDs = mapHyp.getParentHash()

		""" control state of maximum likelihood particle determine's location of shoot frame """
		controlPoses = mapHyp.getControlPoses()
		globalControlPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)


		for k in pathIDs:
			"""
			xP = []
			yP = []
			
			for p in mapHyp.paths[k]:
				xP.append(p[0])
				yP.append(p[1])
			ax3.plot(xP,yP, color=self.colors[k], linewidth=4)
			"""

			shootControlPose_G = globalControlPoses_G[k]
			currFrame = Pose(shootControlPose_G)

			#for path in mapHyp.medialLongPaths[k]:
			for path in mapHyp.localLeaf2LeafPathJunctions[k]["longPaths"]:
				xP = []
				yP = []
				for p in path:
					p1 = currFrame.convertLocalToGlobal(p)
					xP.append(p1[0])
					yP.append(p1[1])

				ax3.plot(xP,yP, color=self.colors[k], linewidth=4)


			nodeSet = mapHyp.getNodes(k)

			print "drawing pathID", k, "for nodes:", nodeSet
			"""
			for nodeID in nodeSet:
				xP = []
				yP = []

				estPose1 = mapHyp.nodePoses[nodeID]
		
				if self.poseData.isBowties[nodeID]:			
					hull1 = self.poseData.aHulls[nodeID]
					#medial1 = self.poseData.medialAxes[nodeID]
					#hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
				else:
					#hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
					hull1 = self.poseData.aHulls[nodeID]
		
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
					ax3.plot(xP,yP, color=(0,0,0))
				elif nodeID == highestNodeID-1:
					ax3.plot(xP,yP, color=(0.5,0.5,0.5))
				else:
					ax3.plot(xP,yP, color=self.colors[k])
			"""
					
			
			xP = []
			yP = []
			for p in mapHyp.hulls[k]:
				xP.append(p[0])
				yP.append(p[1])
				
			#ax3.plot(xP,yP, '--', color=self.colors[k], linewidth=4)
			#ax4.plot(xP,yP, '--', color=self.colors[k], linewidth=4)
			ax3.plot(xP,yP, color=self.colors[k], linewidth=1)
			#ax4.plot(xP,yP, color=self.colors[k], linewidth=4)

		for k in pathIDs:

			globJuncPose = mapHyp.getGlobalJunctionPose(k)
			if globJuncPose != None:
				ax3.scatter([globJuncPose[0],], [globJuncPose[1],], color='k', zorder=10)


		trimmedPaths = mapHyp.trimmedPaths
		
		#for k in range(len(trimmedPaths)):
		for k,path in trimmedPaths.iteritems():
			#path = trimmedPaths[k]
			print "path has", len(path), "points"
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			ax4.plot(xP,yP, color = self.colors[k], linewidth=4)

			globJuncPose = mapHyp.getGlobalJunctionPose(k)
			if globJuncPose != None:
				ax4.scatter([globJuncPose[0],], [globJuncPose[1],], color='k')
			
		#pylab.xlim(-4,4)

		#self.plotEnv()
		self.plotEnv(ax1)
		self.plotEnv(ax2)
		self.plotEnv(ax3)
		self.plotEnv(ax4)
		
		print "quadPath:", self.pathDrawCount
		printStack()

		#pylab.xlim(-10, 12)
		#pylab.ylim(-10, 10)
		#ax3.set_title("paths: %s numNodes: %d %d, hyp %d %3.2f" % (repr(mapHyp.getPathIDs()), self.poseData.numNodes, highestNodeID, mapHyp.hypothesisID, mapHyp.utility))
		#ax4.set_title("paths: %s numNodes: %d %d, hyp %d %3.2f" % (repr(mapHyp.getPathIDs()), self.poseData.numNodes, highestNodeID, mapHyp.hypothesisID, mapHyp.utility))

		ax1.set_title("Pose Robot Static Postures")
		ax2.set_title("Pose Alpha Shapes")
		ax3.set_title("Pose Union Medial Axes")
		ax4.set_title("Trimmed Paths")

		#foo = "paths: %s, nodeID: %d, hyp %d %3.2f" % (repr(mapHyp.getPathIDs()), highestNodeID, mapHyp.hypothesisID, mapHyp.utility)
		fig.suptitle("paths: %s, nodeID: %d, hyp %d %3.2f" % (repr(mapHyp.getPathIDs()), highestNodeID, mapHyp.hypothesisID, mapHyp.mapOverlapSum), fontsize=18, y=0.99)
		plt.savefig("quadPath_%04u_%04u.png" % (highestNodeID, mapHyp.hypothesisID))

		plt.clf()
		plt.close()

		self.pathDrawCount += 1
			

	@logFunction
	def drawPathAndHull(self, mapHyp):
		
		" 1) plot the pose of local splines and postures "
		" 2) plot the alpha shape of the union of pose alpha shapes and medial axis tree "
		" 3) plot the alpha shape of the union of pose alpha shapes and the long medial axis "
		" 4) plot the trimmed path "

		#fig = plt.figure()
		#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
		#ax1.plot(x)
		
		pylab.clf()
		pathIDs = mapHyp.getPathIDs()

		allNodes = []
		for k in pathIDs:
			nodeSet = mapHyp.getNodes(k)
			allNodes += copy(nodeSet)
		
		allNodes.sort() 
		
		if len(allNodes) > 0:
			highestNodeID = allNodes[-1]
		else:
			highestNodeID = 1e100		
			
			
		for k in pathIDs:
			xP = []
			yP = []
			
			for p in mapHyp.paths[k]:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=self.colors[k], linewidth=4)


			nodeSet = mapHyp.getNodes(k)

			print "drawing pathID", k, "for nodes:", nodeSet
			for nodeID in nodeSet:
				xP = []
				yP = []

				estPose1 = mapHyp.getNodePose(nodeID)
		
				if self.poseData.isBowties[nodeID]:			
					hull1 = self.poseData.aHulls[nodeID]
					#medial1 = self.poseData.medialAxes[nodeID]
					#hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
				else:
					#hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
					hull1 = self.poseData.aHulls[nodeID]
		
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
			for p in mapHyp.hulls[k]:
				xP.append(p[0])
				yP.append(p[1])
				
			pylab.plot(xP,yP, '--', color=self.colors[k], linewidth=4)

		allPathIDs = mapHyp.getPathIDs()

		controlPoses_L = {}
		for pathID in pathIDs:
			controlPoses_L[pathID] = self.pathClasses[pathID]["controlPose"]

		parentPathIDs = mapHyp.getParentHash()
		controlPoses_G = computeGlobalControlPoses(controlPoses_L, parentPathIDs)


		xP = []
		yP = []
		for pathID in allPathIDs:
			currFrame = Pose(controlPoses_G[pathID])
			for point_L in mapHyp.localLandmarks[pathID]:
				point_G = currFrame.convertLocalToGlobal(point_L)
				xP.append(point_G[0])
				yP.append(point_G[1])

		pylab.scatter(xP, yP, zorder=9, color='k')


		for k in pathIDs:

			globJuncPose = mapHyp.getGlobalJunctionPose(k)
			if globJuncPose != None:
				pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='m', zorder=10)

		self.plotEnv()
		
		print "pathAndHull:", self.pathDrawCount
		printStack()

		pylab.axis("equal")
		#pylab.xlim(-10, 12)
		#pylab.ylim(-10, 10)
		pylab.title("paths: %s numNodes: %d %d, hyp %d %3.2f" % (repr(mapHyp.getPathIDs()), self.poseData.numNodes, highestNodeID, mapHyp.hypothesisID, mapHyp.utility))
		pylab.savefig("pathAndHull_%04u_%04u.png" % (self.pathDrawCount, mapHyp.hypothesisID))

		self.pathDrawCount += 1
