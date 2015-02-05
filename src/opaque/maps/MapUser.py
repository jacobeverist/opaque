import BayesMapper
from LocalNode import LocalNode, getLongestPath, computeHullAxis
from MapProcess import getInPlaceGuess, getStepGuess
from Splices import getMultiDeparturePoint, orientPath
from Pose import Pose
import pylab
import numpy
from math import pi
import matplotlib.patches as mpatches
from SplineFit import SplineFit
from math import *
from functions import *

import cPickle as pickle
from StablePose import StablePose
#from math import *

# Map Space Parameters
CONST_PIXELSIZE = 0.05
CONST_MIDJOINT = 19

class MapUser:

	def __init__(self, probe, contacts, isStable = False, args = None):
		
		#
		self.isStable = isStable
		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)

		self.args = args
		
		#self.mapAlgorithm = PoseGraph(self.probe, self.contacts)
		self.mapAlgorithm = BayesMapper.BayesMapper(self.probe.robotParam, self.probe.getWalls(), args = self.args)	

		self.initPose = self.probe.getActualJointPose(CONST_MIDJOINT)
		
		self.currNode = 0
		self.foreNode = 0
		self.backNode = 0
		self.localNodes = []

		self.pixelSize = CONST_PIXELSIZE

		self.navPlotCount = 0

		" whether or not the data needs to be synchronized "
		self.isDirty = False

	def drawConstraints(self, id = []):

		pass


	def drawWalls(self):
		pass

	def drawNavigation(self, path, goalPoint):

		pylab.clf()

		segWidth = self.probe.segWidth
		segLength = self.probe.segLength

		print "segWidth, segLength:", segWidth, segLength

		walls = self.probe.getWalls()

		for wall in walls:
			xP = []
			yP = []
			for p in wall:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color='k', alpha=0.5)

		for i in range(0,40):
			#pose = self.probe.getActualJointPose(i)
			pose = self.probe.getActualSegPose(i)
			xTotal = pose[0]
			zTotal = pose[1]
			#totalAngle = pose[2] - pi
			totalAngle = pose[2]

			print "pose", i, pose


			p1 = [xTotal + 0.5*segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - 0.5*segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - 0.5*segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
			#p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			#p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			#p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			#p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
			
			xP = []
			yP = []
			xP.append(p4[0])
			xP.append(p3[0])
			xP.append(p2[0])
			xP.append(p1[0])
			xP.append(p4[0])
			yP.append(p4[1])
			yP.append(p3[1])
			yP.append(p2[1])
			yP.append(p1[1])
			yP.append(p4[1])

			xy = numpy.array([p4,p3,p2,p1])

			# add a rectangle
			polygon = mpatches.Polygon(xy, closed=True, color='k', alpha=0.2)
			pylab.gca().add_patch(polygon)
			#pylab.plot(xP,yP, color='k', alpha=0.5)

		#rootPose = self.contacts.getAverageSegPose(CONST_MIDJOINT)
		rootPose = self.contacts.getAveragePose(CONST_MIDJOINT)
		mapHyps = self.mapAlgorithm.mapHyps

		for hypID, mapHyp in mapHyps.iteritems():

			#if len(mapHyp.nodePoses) > 0:
			if True:

				#rootPose = mapHyp.getNodeRawPose(len(mapHyp.nodePoses)-1)

				#xTotal = 0.0
				#zTotal = 0.0
				#totalAngle = 0.0
				xTotal = rootPose[0]
				zTotal = rootPose[1]
				totalAngle = rootPose[2]

				joints = range(-1,CONST_MIDJOINT)
				joints.reverse()

				for i in joints:

					totalAngle = totalAngle + self.probe.getServo(i+1)
					totalAngle = normalizeAngle(totalAngle)

					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

					xTotal = xTotal - segLength*cos(totalAngle)
					zTotal = zTotal - segLength*sin(totalAngle)

					xP = []
					yP = []
					xP.append(p4[0])
					xP.append(p3[0])
					xP.append(p2[0])
					xP.append(p1[0])
					xP.append(p4[0])
					yP.append(p4[1])
					yP.append(p3[1])
					yP.append(p2[1])
					yP.append(p1[1])
					yP.append(p4[1])

					xy = numpy.array([p4,p3,p2,p1])

					# add a rectangle
					polygon = mpatches.Polygon(xy, closed=True, color='b', alpha=0.2)
					pylab.gca().add_patch(polygon)
					#pylab.plot(xP,yP, color='k', alpha=0.5)
					#actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])

				joints = range(CONST_MIDJOINT, self.probe.numSegs-1)
				
				#xTotal = 0.0
				#zTotal = 0.0
				#totalAngle = 0.0
				xTotal = rootPose[0]
				zTotal = rootPose[1]
				totalAngle = rootPose[2]

				for i in joints:
					xTotal = xTotal + segLength*cos(totalAngle)
					zTotal = zTotal + segLength*sin(totalAngle)
				
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

					xP = []
					yP = []
					xP.append(p4[0])
					xP.append(p3[0])
					xP.append(p2[0])
					xP.append(p1[0])
					xP.append(p4[0])
					yP.append(p4[1])
					yP.append(p3[1])
					yP.append(p2[1])
					yP.append(p1[1])
					yP.append(p4[1])

					xy = numpy.array([p4,p3,p2,p1])

					# add a rectangle
					polygon = mpatches.Polygon(xy, closed=True, color='b', alpha=0.2)
					pylab.gca().add_patch(polygon)
					#pylab.plot(xP,yP, color='k', alpha=0.5)
					#actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])
					
					if i < self.probe.numSegs-2:
						totalAngle = totalAngle - self.probe.getServo(i+1)
						totalAngle = normalizeAngle(totalAngle)


		"""
		for i in range(0,39):

			#pose = self.probe.getActualJointPose(i)
			#pose = self.probe.getActualSegPose(i)
			#pose = self.probe.getJointPose(rootPose, CONST_MIDJOINT, i)
			#pose = self.probe.getJointPose([0.0,0.0,0.0], CONST_MIDJOINT, i)
			pose = self.currNode.getJointPose(i)
			#pose = self.contacts.getAverageSegPose(i)
			xTotal = pose[0]
			zTotal = pose[1]
			#totalAngle = pose[2] - pi
			totalAngle = pose[2]

			print "pose", i, pose

			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]


	
			xP = []
			yP = []
			xP.append(p4[0])
			xP.append(p3[0])
			xP.append(p2[0])
			xP.append(p1[0])
			xP.append(p4[0])
			yP.append(p4[1])
			yP.append(p3[1])
			yP.append(p2[1])
			yP.append(p1[1])
			yP.append(p4[1])

			xy = numpy.array([p4,p3,p2,p1])

			# add a rectangle
			polygon = mpatches.Polygon(xy, closed=True, color='b', alpha=0.2)
			pylab.gca().add_patch(polygon)
			#pylab.plot(xP,yP, color='k', alpha=0.5)
		"""

		if path != None and len(path) > 0:

			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			#pylab.plot(xP,yP, color='r', alpha=0.5, zorder=2)
			pylab.plot(xP,yP, color='r', linewidth=2, zorder=2)

		if goalPoint != None and len(goalPoint) > 0:
			pylab.scatter([goalPoint[0],], [goalPoint[1],], color='k', zorder=3)

		activeHypID = self.mapAlgorithm.activeHypID
		mapHyps = self.mapAlgorithm.mapHyps

		for hypID, mapHyp in mapHyps.iteritems():

			if len(mapHyp.nodePoses) > 0:
				allSplices = mapHyp.getAllSplices2()
				for sPath in allSplices:
					path = sPath['skelPath']
	
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])
					pylab.plot(xP,yP, color='g', linewidth=2, zorder=1)


		if self.numNodes > 0:
			nodeID = self.numNodes-1

			hull1 = self.mapAlgorithm.poseData.aHulls[nodeID]
			medial1 = self.mapAlgorithm.poseData.medialAxes[nodeID]

			estPose1 = self.mapAlgorithm.mapHyps[0].getNodePose(nodeID)

			poseOrigin = Pose(estPose1)

			medialSpline1 = SplineFit(medial1, smooth=0.1)
			points1 = medialSpline1.getUniformSamples(interpAngle=True)
		
			points1_offset = []
			xP = []
			yP = []
			for p in points1:
				result = poseOrigin.convertLocalOffsetToGlobal(p)
				points1_offset.append(result)
				xP.append(result[0])
				yP.append(result[1])

			pylab.plot(xP, yP, color='k')

		if activeHypID != None:
			pass

			#print "activeHypID:", activeHypID, mapHyps[activeHypID].pathClasses.keys()

			#allSplices, terminals, junctions = mapHyps[activeHypID].getAllSplices()

			#for k, result in allSplices.iteritems():
			#	for sPath in result:
			#		#path = sPath['path']
			#		path = sPath['skelPath']
	
			#		xP = []
			#		yP = []
			#		for p in path:
			#			xP.append(p[0])
			#			yP.append(p[1])
			#		pylab.plot(xP,yP, color='g', zorder=1)



			#for k,trimPath in mapHyps[activeHypID].trimmedPaths.iteritems():
			#	print "path has", len(trimPath), "points"
			#	xP = []
			#	yP = []
			#	for p in trimPath:
			#		xP.append(p[0])
			#		yP.append(p[1])

			#	pylab.plot(xP,yP, color = 'g', alpha=0.7)

		pylab.gca().invert_yaxis()
		pylab.axis("equal")
		#pylab.xlim(-4.0,2.0)
		#pylab.ylim(-2.0,2.0)
		pylab.title("%d poses" % self.mapAlgorithm.poseData.numNodes)
		pylab.savefig("navEstimate_%04u.png" % self.navPlotCount)
		self.navPlotCount += 1

	def localizePose(self):
		
		self.foreNode.saveMap()
		self.backNode.saveMap()
	
		
		self.mapAlgorithm.loadNewNode(self.foreNode)
		self.mapAlgorithm.loadNewNode(self.backNode)

		self.localNodes.append(self.foreNode)
		self.localNodes.append(self.backNode)
#
		with open("map_%04u.obj" % self.backNode.nodeID, 'wb') as output:
			pickle.dump(self.mapAlgorithm, output, pickle.HIGHEST_PROTOCOL)

	@property
	def numNodes(self):
		return self.mapAlgorithm.numNodes

	@property
	def nodeHash(self):
		return self.mapAlgorithm.nodeHash
		#return self.mapAlgorithm.localNodes


	def restorePickle(self, dirName, nodeID):
		print self

		#with open('map_%04u.obj' % nodeID, 'rb') as inputVal:
		with open(dirName + '/' + 'map_%04u.obj' % nodeID, 'rb') as inputVal:
		
			print pickle
			print inputVal
			print self
			print self.mapAlgorithm
			self.mapAlgorithm = pickle.load(inputVal)
			print self.mapAlgorithm
			self.mapAlgorithm.saveState()

			for val in self.mapAlgorithm.mapHyps.values():
				self.mapAlgorithm.drawPathAndHull2(val)


	def restoreSeries(self, dirName, num_poses):

		self.mapAlgorithm.restoreState(dirName, num_poses)

		for i in range(0, num_poses):

			print "loading node", i			
			self.currNode = LocalNode(self.probe, self.contacts, i, CONST_MIDJOINT, self.pixelSize)
			self.currNode.readFromFile2(dirName, i)			
			self.currNode.setGPACPose(self.mapAlgorithm.nodePoses[i])
			self.mapAlgorithm.nodeHash[i] = self.currNode

		
		foreNode = LocalNode(self.probe, self.contacts, num_poses, CONST_MIDJOINT, self.pixelSize)
		foreNode.readFromFile2(dirName, num_poses)
		backNode = LocalNode(self.probe, self.contacts, num_poses+1, CONST_MIDJOINT, self.pixelSize)
		backNode.readFromFile2(dirName, num_poses+1)
		self.mapAlgorithm.insertPose(foreNode, backNode, initLocation = foreNode.getEstPose())
		
		self.isDirty = True
		
		self.synch()
		
	def loadSeries(self, dirName, num_poses):
	
		numNodes = 0
		self.currNode = 0
			
		for i in range(0, num_poses):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, CONST_MIDJOINT, self.pixelSize)
			currNode.readFromFile(dirName, i)
			self.mapAlgorithm.loadNewNode(currNode)
			self.drawConstraints(i)

		self.isDirty = True

		self.synch()

	def loadFile(self, dirName, num_poses):
		
		numNodes = 0
		self.currNode = 0		
		
		nodeHash = {}
		
		for i in range(0,num_poses):
			print "loading node", i			
			self.currNode = LocalNode(self.probe, self.contacts, i, CONST_MIDJOINT, self.pixelSize)
			self.currNode.readFromFile(dirName, i)
			nodeHash[i] = self.currNode

			numNodes += 1

		" add the nodes to the graph "
		self.mapAlgorithm.loadNodes(nodeHash)
		
		" add the constraint edges to the graph "
		self.mapAlgorithm.loadConstraints(dirName, num_poses-1)

		self.isDirty = True

		#self.drawConstraints()

		" save the maps "		
		#self.synch()

	def newInPlaceNode(self, faceDir, travelDir):
		
		#if self.numNodes > 0:
		#	
		#	if self.currNode.isDirty():
		#		self.currNode.synch()
		#
		#	self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, CONST_MIDJOINT, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir)
		
		self.mapAlgorithm.addInPlaceNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

	def loadWalls(self, walls):
		self.walls = walls

	def getWalls(self):
		return self.walls
	
	def newNode(self, faceDir, travelDir):
				
		#if self.numNodes > 0:
		#	
		#	if self.currNode.isDirty():
		#		self.currNode.synch()
		#
		#	self.currNode.saveToFile()

		
		if faceDir:
			self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, CONST_MIDJOINT, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir, useBloom = self.args.bloomFeature, useBend = self.args.bendFeature)		
			self.foreNode = self.currNode
		else:
			self.currNode = LocalNode(self.probe, self.contacts, self.numNodes+1, CONST_MIDJOINT, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir, useBloom = self.args.bloomFeature, useBend = self.args.bendFeature)		
			self.backNode = self.currNode
		
		#self.mapAlgorithm.addNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
		
		#if self.mapAlgorithm.numNodes > 100:
		#	print "max number of nodes reached, returning"
		#	raise

	def pairDone(self):

		" pair this last node, and the node before as a paired set for a local pose "
		#node1 = self.nodeHash[self.numNodes-1]
		#node2 = self.nodeHash[self.numNodes-2]

		self.foreNode.setPartnerNodeID(self.backNode.nodeID)
		self.backNode.setPartnerNodeID(self.foreNode.nodeID)

		#node1.setPartnerNodeID(self.numNodes-2)
		#node2.setPartnerNodeID(self.numNodes-1)
		
		#self.mapAlgorithm.pairLastTwo()

	def insertPose(self, estPose, direction):

		print "estPose:", estPose

		self.foreNode = LocalNode(self.probe, self.contacts, self.numNodes, CONST_MIDJOINT, self.pixelSize, stabilizePose = self.isStable, faceDir = True, travelDir = direction)		
		self.foreNode.setEstPose(estPose)
		self.foreNode.setGndPose(estPose)

		print "stable:", self.contacts.getAverageSegPose(CONST_MIDJOINT)

		print "insert:", self.foreNode.getGlobalGPACPose()
		self.currNode = self.foreNode
		#self.forceUpdate(True)		

		self.stablePose.setDirection(True)
		self.currNode.update()

		self.backNode = LocalNode(self.probe, self.contacts, self.numNodes+1, CONST_MIDJOINT, self.pixelSize, stabilizePose = self.isStable, faceDir = False, travelDir = direction)		
		self.backNode.setEstPose(estPose)
		self.backNode.setGndPose(estPose)
		print "insert:", self.backNode.getGlobalGPACPose()
		self.currNode = self.backNode
		#self.forceUpdate(False)		
		self.stablePose.setDirection(False)
		self.currNode.update()

		print "insert:", self.foreNode.getGlobalGPACPose()
		print "insert:", self.backNode.getGlobalGPACPose()
		self.mapAlgorithm.insertNewNode(self.foreNode)
		self.mapAlgorithm.insertNewNode(self.backNode)
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)	

	def synch(self):
		
		self.currNode.synch()
		
		print "LocalNode Synched"

		self.isDirty = False


	def correctPosture(self):
		self.currNode.updateCorrectPosture()

	def forceUpdate(self, isForward=True):


		self.stablePose.setDirection(isForward)
		self.currNode.update()

		nodeID = self.currNode.nodeID
		direction = self.currNode.travelDir
		distEst = 1.0

		mapHyp = self.getMaxHyp()

		poseData = deepcopy(mapHyp.poseData)

		if nodeID > 0:

			foreNodeID = self.foreNode.nodeID

			hull1, medial1 = computeHullAxis(foreNodeID, self.foreNode, tailCutOff = False)
			poseData.aHulls[foreNodeID] = hull1
			poseData.medialAxes[foreNodeID] = medial1
			poseData.numLeafs[foreNodeID] = self.foreNode.getNumLeafs()
			poseData.faceDirs[foreNodeID] = self.foreNode.faceDir
			poseData.isBowties[foreNodeID] = self.foreNode.isBowtie			
			poseData.medialLongPaths[foreNodeID] = self.foreNode.medialLongPaths
			poseData.correctedPostures[foreNodeID] = self.foreNode.getStableGPACPosture()
			poseData.isNodeFeatureless[foreNodeID] = self.foreNode.getIsFeatureless()
			poseData.frontProbeError[foreNodeID] = self.foreNode.frontProbeError
			poseData.backProbeError[foreNodeID] = self.foreNode.backProbeError
			poseData.travelDirs[foreNodeID] = self.foreNode.travelDir
			poseData.spatialFeatures[foreNodeID] = self.foreNode.spatialFeatures
			poseData.rootPoses[foreNodeID] = self.foreNode.rootPose
			poseData.localPostures[foreNodeID] = self.foreNode.correctedPosture

			if nodeID % 2 == 0:

				estPose0 = mapHyp.getNodePose(nodeID-2)
				estPose2 = self.foreNode.getGlobalGPACPose()
			
				#estPose0 = mapHyp.getNodePose(nodeID-2)
				#estPose2 = mapHyp.getNodePose(nodeID)
				newPose2 = getStepGuess(poseData, nodeID-2, nodeID, estPose0, estPose2, direction)

				self.foreNode.setGPACPose(newPose2)

			if nodeID % 2 == 1:
				
				foreNodeID = self.foreNode.nodeID
				backNodeID = self.backNode.nodeID

				hull1, medial1 = computeHullAxis(backNodeID, self.backNode, tailCutOff = False)
				poseData.aHulls[backNodeID] = hull1
				poseData.medialAxes[backNodeID] = medial1
				poseData.numLeafs[backNodeID] = self.backNode.getNumLeafs()
				poseData.faceDirs[backNodeID] = self.backNode.faceDir
				poseData.isBowties[backNodeID] = self.backNode.isBowtie			
				poseData.medialLongPaths[backNodeID] = self.backNode.medialLongPaths
				poseData.correctedPostures[backNodeID] = self.backNode.getStableGPACPosture()
				poseData.isNodeFeatureless[backNodeID] = self.backNode.getIsFeatureless()
				poseData.frontProbeError[backNodeID] = self.backNode.frontProbeError
				poseData.backProbeError[backNodeID] = self.backNode.backProbeError
				poseData.travelDirs[backNodeID] = self.backNode.travelDir
				poseData.spatialFeatures[backNodeID] = self.backNode.spatialFeatures
				poseData.rootPoses[backNodeID] = self.backNode.rootPose
				poseData.localPostures[backNodeID] = self.backNode.correctedPosture

				" the guess process gives a meaningful guess for these node's poses "
				" before this, it is meaningless "

				estPose2 = self.foreNode.getGlobalGPACPose()
				estPose3 = self.backNode.getGlobalGPACPose()

				#estPose3 = mapHyp.getNodePose(nodeID)
				supportLine = mapHyp.paths[0]
				newPose3 = getInPlaceGuess(poseData, foreNodeID, backNodeID, estPose2, estPose3, supportLine, direction)

				self.backNode.setGPACPose(newPose3)


	def getDistanceToPath(self, point, wayPath):
		
		pathSpline = SplineFit(wayPath)
		minDist, uVal, uPoint = pathSpline.findClosestPoint(point)
		
		return minDist

	def update(self, isForward=True):

		self.isDirty = True

		self.stablePose.setDirection(isForward)

		#self.currNode.update()

		" if the configuration is a stable pose, the maximum displacement is not exceeded "
		if self.stablePose.isStable():
			self.currNode.update()
		else:
			pass
	
	def getNearestPathPoint(self, originPoint):
		
		newPoint = self.mapAlgorithm.getNearestPathPoint(originPoint)
		
		return newPoint
	
	def getMaxHyp(self):

		activeHypID = self.mapAlgorithm.activeHypID
		mapHyps = self.mapAlgorithm.mapHyps

		if activeHypID != None:
			return mapHyps[activeHypID]

		else:
			""" return first """
			for hypID, mapHyp in mapHyps.iteritems():
				return mapHyp

		raise
	
	def getMaxPose(self):

		nodeID = self.numNodes-1

		activeHypID = self.mapAlgorithm.activeHypID
		mapHyps = self.mapAlgorithm.mapHyps

		if activeHypID != None:
			""" return active """
			print "return active:", nodeID
			estPose1 = mapHyps[activeHypID].getNodeRawPose(nodeID)

			return estPose1

		else:
			""" return first """
			print "return first:", nodeID
			for hypID, mapHyp in mapHyps.iteritems():
				if len(mapHyp.nodePoses) > self.numNodes:
					return mapHyp.getNodeRawPose(nodeID)

		print "return currNode.getEstPose():", nodeID
		return self.currNode.getEstPose()
	
	def getDirection(self, wayPath):

		nodeID = self.numNodes-1

		hull1 = self.mapAlgorithm.poseData.aHulls[nodeID]
		medial1 = self.mapAlgorithm.poseData.medialAxes[nodeID]

		estPose1 = self.mapAlgorithm.mapHyps[self.mapAlgorithm.activeHypID].getNodePose(nodeID)

		poseOrigin = Pose(estPose1)

		medialSpline1 = SplineFit(medial1, smooth=0.1)
		points1 = medialSpline1.getUniformSamples(interpAngle=True)
	
		points1_offset = []
		for p in points1:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points1_offset.append(result)
	
		globalMedialSpline = SplineFit(points1_offset, smooth=0.1)
		globalMedialPoints = globalMedialSpline.getUniformSamples()			

		orientedWayPath = orientPath(wayPath, points1_offset)

		dist1 = sqrt((orientedWayPath[0][0]-wayPath[0][0])**2 + (orientedWayPath[0][1]-wayPath[0][1])**2)
		dist2 = sqrt((orientedWayPath[-1][0]-wayPath[0][0])**2 + (orientedWayPath[-1][1]-wayPath[0][1])**2)

		if dist1 < dist2:
			return False
		else:
			return True

		#" - get position of both tips "
		#frontTip = globalMedialPoints[0]
		#backTip = globalMedialPoints[-1]

		#results0 = getMultiDeparturePoint(splice, medial0_vec, currPose, currPose, [0,], nodeID, spliceIndex=k, plotIter=False)

	def isFollowing(self, wayPath, foreAvg = 0.0, backAvg = 0.0):

		" - get splice of wayPath "
		
		" - get GPAC of pose "
		nodeID = self.numNodes-1

		hull1 = self.mapAlgorithm.poseData.aHulls[nodeID]
		medial1 = self.mapAlgorithm.poseData.medialAxes[nodeID]

		estPose1 = self.mapAlgorithm.mapHyps[self.mapAlgorithm.activeHypID].getNodePose(nodeID)

		poseOrigin = Pose(estPose1)

		medialSpline1 = SplineFit(medial1, smooth=0.1)
		points1 = medialSpline1.getUniformSamples(interpAngle=True)
	
		points1_offset = []
		for p in points1:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points1_offset.append(result)
	
		globalMedialSpline = SplineFit(points1_offset, smooth=0.1)
		globalMedialPoints = globalMedialSpline.getUniformSamples()			
		
		" - get position of both tips "
		frontTip = globalMedialPoints[0]
		backTip = globalMedialPoints[-1]
		#frontTipPathPoint = self.getNearestPathPoint(frontTip)
		#backTipPathPoint = self.getNearestPathPoint(backTip)
		
		pathSpline = SplineFit(wayPath)
		#minDist, uVal, uPoint = pathSpline.findClosestPoint(dest)
		
		frontMinDist, frontUVal, frontUPoint = pathSpline.findClosestPoint(frontTip)
		backMinDist, backUVal, backUPoint = pathSpline.findClosestPoint(backTip)

		results0 = getMultiDeparturePoint(wayPath, points1, estPose1, estPose1, [0,], nodeID, spliceIndex=0, plotIter=False)
		contigFrac = results0[12]

		#frontMinDist, frontUVal, frontUPoint = pathSpline.findClosestPoint(frontTipPathPoint)
		#backMinDist, backUVal, backUPoint = pathSpline.findClosestPoint(backTipPathPoint)

		return frontMinDist, backMinDist, contigFrac
		
	def isDestReached(self, dest, wayPath, foreAvg = 0.0, backAvg = 0.0):
		#dest = self.localWayPoints[0]			

		" - get splice of wayPath "
		
		" - get GPAC of pose "
		nodeID = self.numNodes-1
		#hull1, medial1 = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull1 = self.mapAlgorithm.poseData.aHulls[nodeID]
		medial1 = self.mapAlgorithm.poseData.medialAxes[nodeID]

		#estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()
		estPose1 = self.mapAlgorithm.mapHyps[self.mapAlgorithm.activeHypID].getNodePose(nodeID)

		poseOrigin = Pose(estPose1)

		medialSpline1 = SplineFit(medial1, smooth=0.1)
		points1 = medialSpline1.getUniformSamples(interpAngle=True)
	
		points1_offset = []
		for p in points1:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points1_offset.append(result)
	
		globalMedialSpline = SplineFit(points1_offset, smooth=0.1)
		globalMedialPoints = globalMedialSpline.getUniformSamples()			
		
		" - get position of both tips "
		frontTip = globalMedialPoints[0]
		backTip = globalMedialPoints[-1]
		frontTipPathPoint = self.getNearestPathPoint(frontTip)
		backTipPathPoint = self.getNearestPathPoint(backTip)
		
		pathSpline = SplineFit(wayPath)
		minDist, uVal, uPoint = pathSpline.findClosestPoint(dest)
		
		frontMinDist, frontUVal, frontUPoint = pathSpline.findClosestPoint(frontTipPathPoint)
		backMinDist, backUVal, backUPoint = pathSpline.findClosestPoint(backTipPathPoint)
		
		destDist = pathSpline.dist_u(uVal)
		frontDist = pathSpline.dist_u(frontUVal)
		backDist = pathSpline.dist_u(backUVal)


		frontMag = fabs(destDist-frontDist)
		backMag = fabs(destDist-backDist)
		
		if frontMag < backMag:
			goalDist = frontMag
		else:
			goalDist = backMag

		# FIXME:  this collision checks should go in the localization section
		# which would automatically localize us to the terminal point
		if frontMag < 0.8 and foreAvg >= 1.1:
			destReached = True

		if backMag < 0.8 and backAvg >= 1.1:
			destReached = True

		destReached = False
		if frontDist >= destDist and destDist >= backDist:
			" objective achieved!"
			destReached = True
			
		elif frontDist <= destDist and destDist <= backDist:
			" objective achieved!"
			destReached = True

		elif goalDist < 0.5:
			" objective achieved!"
			destReached = True

		print "isDestReached:", destReached, frontMag, backMag, foreAvg, backAvg, goalDist
			
		#if destReached:
		#	goalDist = 0.0

		return destReached, goalDist
				
	def computeGraphPath(self, startPose, endPose):
		
		path = self.mapAlgorithm.computeNavigationPath(startPose, endPose)

		return path

	def getPathLength(self, startPose, endPose, path):
		
		
		minI1 = 0
		minDist1 = 1e100
		minI2 = 0
		minDist2 = 1e100
		for i in range(len(path)):

			p0 = path[i]
			dist = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
			
			if dist < minDist1:
				minDist1 = dist
				minI1 = i

			dist = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

			if dist < minDist2:
				minDist2 = dist
				minI2 = i
		
		" find distance between minI1 and minI2 "
		if minI1 < minI2:
			lowIndex = minI1
			highIndex = minI2+1
		else:
			lowIndex = minI2
			highIndex = minI1+1
		
		totalDist = 0.0		
		for i in range(lowIndex, highIndex-1):
			p1 = path[i]
			p2 = path[i+1]
			
			dist = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
			
			totalDist += dist

		
		return totalDist + minDist1 + minDist2

	def pathTermVisited(self, termID):
		for hypID, mapHyp in self.mapAlgorithm.mapHyps.iteritems():
			mapHyp.pathTermVisited(termID)
		
	def selectNextDestination(self):

		if self.isDirty:
			self.synch()

		allTerms = {}
		allTermsVisited = {}
		termToHypID = {}

		#targetPoint = [2.255,5.93]
		targetPoint = [3.296, 0.48]

		" get the termination point and orientation of each path "
		for hypID, mapHyp in self.mapAlgorithm.mapHyps.iteritems():
			terms = mapHyp.getPathTerms()
			termsVisited = mapHyp.getPathTermsVisited()

			print "selectNextDestination for hyp", hypID
			print "terms:", terms
			print "termsVisited:", termsVisited

			for key, val in terms.iteritems():
				allTerms[key] = val
				allTermsVisited[key] = termsVisited[key]
				termToHypID[key] = hypID

		print "allTerms:", allTerms
		print "allTermsVisited:", allTermsVisited
		
		#for k, term in allTerms.iteritems():

		#	dist = sqrt((term[0]-targetPoint[0])**2+(term[1]-targetPoint[1])**2)
		#	print "targetDist=", dist, term, targetPoint

		#	if dist < 0.5:
		#		print "selecting term", k, "for hyp ID", termToHypID[k]

		#		self.mapAlgorithm.activeHypID = termToHypID[k]
		#		return term,k
		
		for k, term in allTerms.iteritems():

			if not allTermsVisited[k]:
				print "selecting term", k, "for hyp ID", termToHypID[k]

				self.mapAlgorithm.activeHypID = termToHypID[k]
				return term,k
		
		" if all terms visited, reset and go to root "
		for hypID, mapHyp in self.mapAlgorithm.mapHyps.iteritems():
			mapHyp.resetTerms()
		
		minDist = 1e100
		minPnt = None
		minTermID = None
		for hypID, mapHyp in self.mapAlgorithm.mapHyps.iteritems():
			terms = mapHyp.getPathTerms()

			for key,val in terms.iteritems():

				rootPoint = mapHyp.rootPoint
				dist = sqrt((rootPoint[0]-val[0])**2 + (rootPoint[1]-val[1])**2)
				if dist < minDist:
					minDist = dist
					minPnt = val
					minTermID = key

		return minPnt, minTermID
	
	def getDestination(self, termID, termPoint):

		if self.isDirty:
			self.synch()

		allTerms = {}
		allTermsVisited = {}
		termToHypID = {}

		" get the termination point and orientation of each path "
		for hypID, mapHyp in self.mapAlgorithm.mapHyps.iteritems():
			terms = mapHyp.getPathTerms()
			termsVisited = mapHyp.getPathTermsVisited()

			print "selectNextDestination for hyp", hypID
			print "terms:", terms
			print "termsVisited:", termsVisited

			for key, val in terms.iteritems():
				allTerms[key] = val
				allTermsVisited[key] = termsVisited[key]
				termToHypID[key] = hypID

		newTermPoint = allTerms[termID] 
		return newTermPoint

	def computePathSplice(self, termID, termPoint):

		frontierPoint = self.getDestination(termID, termPoint)


		allSplices = []
		""" find all splices that terminate at the terminal destination """
		allSplices = self.mapAlgorithm.mapHyps[0].getAllSplices2()
		#for hypID, mapHyp in self.mapAlgorithm.mapHyps.iteritems():
		#	allSplices = mapHyp.getAllSplices2()
		#	break

		termSplices = []
		for spliceIndex in range(len(allSplices)):
			sPath = allSplices[spliceIndex]
			path = sPath['skelPath']
			termPoints = sPath['termPath']

			term1 = termPoints[0]
			term2 = termPoints[1]

			dist1 = sqrt((term1[0]-frontierPoint[0])**2 + (term1[1]-frontierPoint[1])**2)
			dist2 = sqrt((term2[0]-frontierPoint[0])**2 + (term2[1]-frontierPoint[1])**2)

			if dist1 < 0.1 or dist2 < 0.1:
				termSplices.append(path)

		
		""" find the splice that is coaligned with current pose's spatial curve """

		nodeID = self.mapAlgorithm.poseData.numNodes-1

		medial0 = self.mapAlgorithm.poseData.medialAxes[nodeID]

		medialSpline0 = SplineFit(medial0, smooth=0.1)
		medial0_vec = medialSpline0.getUniformSamples()

		currPose = self.mapAlgorithm.mapHyps[0].getNodePose(nodeID)

		maxContigFrac = 0.0
		maxSpliceIndex = 0
		for k in range(len(termSplices)):

			splice = termSplices[k]

			results0 = getMultiDeparturePoint(splice, medial0_vec, currPose, currPose, [0,], nodeID, spliceIndex=k, plotIter=False)

			contigFrac0 = results0[12]

			if contigFrac0 >= maxContigFrac:
				maxContigFrac = contigFrac0
				maxSpliceIndex = k

		pathSplice = deepcopy(termSplices[maxSpliceIndex])


		""" destination should be at the n'th index """
		dist1 = sqrt((pathSplice[0][0]-frontierPoint[0])**2 + (pathSplice[0][1]-frontierPoint[1])**2)
		dist2 = sqrt((pathSplice[-1][0]-frontierPoint[0])**2 + (pathSplice[-1][1]-frontierPoint[1])**2)

		if dist1 < dist2:
			pathSplice.reverse()

		frontierPathPoint = self.getNearestPathPoint(frontierPoint)

		return frontierPathPoint, pathSplice

		""" over smooth sharp corners """
		newCurve = []
		rangeWidth = 50

		pose1 = pathSplice[0]
		pose2 = pathSplice[rangeWidth]

		newX = pose1[0] + (pose2[0]-pose1[0]) / 2.0
		newY = pose1[1] + (pose2[1]-pose1[1]) / 2.0

		joinP = [newX,newY]

		dist = sqrt((pose1[0]-joinP[0])**2+(pose1[1]-joinP[1])**2)
		vec = [joinP[0]-pose1[0], joinP[1]-pose1[1]]
		vec[0] /= dist
		vec[1] /= dist
		max_spacing = 0.02
		numPoints = int(dist/max_spacing)

		for j in range(numPoints-1):
			newP = [j*max_spacing*vec[0] + pose1[0], j*max_spacing*vec[1] + pose1[1]]
			newCurve.append(newP)


		for k in range(0,len(pathSplice)-rangeWidth):
			pose1 = pathSplice[k]
			pose2 = pathSplice[k+rangeWidth]

			newX = pose1[0] + (pose2[0]-pose1[0]) / 2.0
			newY = pose1[1] + (pose2[1]-pose1[1]) / 2.0


			newCurve.append([newX,newY])


		pose1 = pathSplice[len(pathSplice)-rangeWidth-1]
		pose2 = pathSplice[len(pathSplice)-1]

		newX = pose1[0] + (pose2[0]-pose1[0]) / 2.0
		newY = pose1[1] + (pose2[1]-pose1[1]) / 2.0

		joinP = [newX,newY]

		dist = sqrt((pose2[0]-joinP[0])**2+(pose2[1]-joinP[1])**2)
		vec = [-joinP[0]+pose2[0], -joinP[1]+pose2[1]]
		vec[0] /= dist
		vec[1] /= dist
		max_spacing = 0.02
		numPoints = int(dist/max_spacing)

		for j in range(1,numPoints):
			newP = [j*max_spacing*vec[0] + pose2[0], j*max_spacing*vec[1] + pose2[1]]
			newCurve.append(newP)

		overSmoothedPath = SplineFit(newCurve, smooth=0.1).getUniformSamples()

		return frontierPathPoint, overSmoothedPath


	def recomputePath(self, termID, termPoint):

		frontierPoint = self.getDestination(termID, termPoint)

		frontPose = self.contacts.getAverageSegPose(0)
		backPose = self.contacts.getAverageSegPose(39)

		frontPathPoint = self.getNearestPathPoint(frontPose)
		backPathPoint = self.getNearestPathPoint(backPose)
		frontierPathPoint = self.getNearestPathPoint(frontierPoint)
		
		goalPath1 = self.computeGraphPath(frontPathPoint, frontierPathPoint)
		goalPath2 = self.computeGraphPath(backPathPoint, frontierPathPoint)

		len1 = self.getPathLength(frontPathPoint, frontierPathPoint, goalPath1)
		len2 = self.getPathLength(backPathPoint, frontierPathPoint, goalPath2)

		if len1 > len2:
			wayPoints = [frontierPathPoint]
			wayPaths = [goalPath1]
		else:
			wayPoints = [frontierPathPoint]
			wayPaths = [goalPath2]
				
		return wayPoints, wayPaths
								
	def getNodeOccMap(self, nodeNum):
		localNode = self.nodeHash[nodeNum]
		return localNode.getOccMap()

	def getNodePose(self, nodeNum):
		localNode = self.nodeHash[nodeNum]
		return localNode.getEstPose()
	
	def setNodePose(self, nodeNum, estPose):
		localNode = self.nodeHash[nodeNum]
		localNode.setEstPose(estPose)
	
	def obstCallBack(self, direction):
		self.currNode.obstCallBack(direction)
	
	def saveLocalMap(self):
		" save the local map "
		if self.foreNode != 0:
			self.foreNode.saveMap()

		if self.backNode != 0:
			self.backNode.saveMap()
			
	def saveState(self):

		self.mapAlgorithm.saveState()
		

	def getCurrentNode(self):
		return self.currNode
	
	" FIXME:  Perhaps this should be moved to the PokeWalls behavior since this the only thing that requires it"
	def keepStablePose(self):
		self.stablePose.update()
	
