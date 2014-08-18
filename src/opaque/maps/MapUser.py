import BayesMapper
from LocalNode import LocalNode, getLongestPath, computeHullAxis
from Pose import Pose
import pylab
from SplineFit import SplineFit

import cPickle as pickle
from StablePose import StablePose
#from math import *

# Map Space Parameters
PIXELSIZE = 0.05
MAPSIZE = 40.0
#MAPSIZE = 20.0

class MapUser:

	def __init__(self, probe, contacts, isStable = False):
		
		#
		self.isStable = isStable
		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)
		
		#self.mapAlgorithm = PoseGraph(self.probe, self.contacts)
		self.mapAlgorithm = BayesMapper.BayesMapper(self.probe.getWalls())	

		self.initPose = self.probe.getActualJointPose(19)
		
		self.currNode = 0
		self.foreNode = 0
		self.backNode = 0
		self.localNodes = []

		self.pixelSize = PIXELSIZE

		" whether or not the data needs to be synchronized "
		self.isDirty = False

	def drawConstraints(self, id = []):
		pass

	def drawWalls(self):
		pass

	def localizePose(self):
		
		self.foreNode.saveMap()
		self.backNode.saveMap()
	
		
		self.mapAlgorithm.loadNewNode(self.foreNode)
		self.mapAlgorithm.loadNewNode(self.backNode)

		self.localNodes.append(self.foreNode)
		self.localNodes.append(self.backNode)

	@property
	def numNodes(self):
		return self.mapAlgorithm.numNodes

	@property
	def nodeHash(self):
		return self.mapAlgorithm.localNodes

	def restoreSeries(self, dirName, num_poses):

		self.mapAlgorithm.restoreState(dirName, num_poses)

		for i in range(0, num_poses):

			print "loading node", i			
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, self.pixelSize)
			self.currNode.readFromFile2(dirName, i)			
			self.currNode.setGPACPose(self.mapAlgorithm.nodePoses[i])
			self.mapAlgorithm.nodeHash[i] = self.currNode

		
		foreNode = LocalNode(self.probe, self.contacts, num_poses, 19, self.pixelSize)
		foreNode.readFromFile2(dirName, num_poses)
		backNode = LocalNode(self.probe, self.contacts, num_poses+1, 19, self.pixelSize)
		backNode.readFromFile2(dirName, num_poses+1)
		self.mapAlgorithm.insertPose(foreNode, backNode, initLocation = foreNode.getEstPose())
		
		self.isDirty = True
		
		self.synch()
		
	def loadSeries(self, dirName, num_poses):
	
		numNodes = 0
		self.currNode = 0
			
		for i in range(0, num_poses):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, 19, self.pixelSize)
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
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, self.pixelSize)
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

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir)
		
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
			self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir)		
			self.foreNode = self.currNode
		else:
			self.currNode = LocalNode(self.probe, self.contacts, self.numNodes+1, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir)		
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


		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = True, travelDir = direction)		
		self.foreNode = self.currNode
		self.forceUpdate(True)		
		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes+1, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = False, travelDir = direction)		
		self.backNode = self.currNode
		self.forceUpdate(False)		

		self.mapAlgorithm.insertPose(self.foreNode, self.backNode, initLocation = estPose)
		#self.mapAlgorithm.addInitNode(self.foreNode, estPose)
		#self.mapAlgorithm.loadNewNode(self.foreNode)
		#self.mapAlgorithm.loadNewNode(self.backNode)

		
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
		
	def update(self, isForward=True):

		self.isDirty = True

		self.stablePose.setDirection(isForward)

		" if the configuration is a stable pose, the maximum displacement is not exceeded "
		if self.stablePose.isStable():
			self.currNode.update()
		else:
			pass
	
	def getNearestPathPoint(self, originPoint):
		
		newPoint = self.mapAlgorithm.getNearestPathPoint(originPoint)
		
		return newPoint
	
	
	def isDestReached(self, dest, wayPath):
		#dest = self.localWayPoints[0]			

		" - get splice of wayPath "
		
		" - get GPAC of pose "
		nodeID = self.numNodes-1
		hull1, medial1 = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)
		estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()

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

		destReached = False
		if frontDist >= destDist and destDist >= backDist:
			" objective achieved!"
			destReached = True
			
		elif frontDist <= destDist and destDist <= backDist:
			" objective achieved!"
			destReached = True
			
		if destReached:
			goalDist = 0.0

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

		
	def selectNextDestination(self):
		if self.isDirty:
			self.synch()
		
		" get the termination point and orientation of each path "
		terms = self.mapAlgorithm.topHyp.getPathTerms()
		termsVisited = self.mapAlgorithm.topHyp.getPathTermsVisited()
		
		print "Terms:", terms
		print "TermsVisited:", termsVisited
		
		for k, term in terms.iteritems():
			#term = terms[k]
			
			if not termsVisited[k]:
				#self.mapAlgorithm.topHyp.pathTermVisited(k)
				print "selecting term", k
				return term,k
		
		" if all terms visited, reset and go to root "
		self.mapAlgorithm.topHyp.resetTerms()
		
		return self.mapAlgorithm.paths.rootPoint, -1
		
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
	
