from PoseGraph import PoseGraph, computeBareHull
from Pose import Pose

from LocalNode import LocalNode
from VoronoiMap import VoronoiMap
from FrontierMap import FrontierMap
from OccupancyMap import OccupancyMap
from FreeSpaceBoundaryMap import FreeSpaceBoundaryMap
from NavRoadMap import NavRoadMap
from StablePose import StablePose
from SplineFit import SplineFit
from PoseGraph import computeHullAxis

import gen_icp


import pylab
from math import *
from random import random

import Image
import ImageDraw

# Map Space Parameters
PIXELSIZE = 0.05
MAPSIZE = 40.0
#MAPSIZE = 20.0

class MapGraph:

	def __init__(self, probe, contacts, isStable = False):
		
		#
		self.isStable = isStable
		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)
		
		self.poseGraph = PoseGraph(self.probe, self.contacts)
	
		self.initPose = self.probe.getActualJointPose(19)
		
		self.currNode = 0
		self.foreNode = 0
		self.backNode = 0

		self.pixelSize = PIXELSIZE
		self.mapSize = MAPSIZE
		self.numPixel = int(2.0 * self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize / 2.0
		self.divPix = floor((2.0 * self.mapSize / self.pixelSize) / self.mapSize)
		self.saveCount = 0
		self.fileName = "mapGraph%04u.png"
		self.boundIteration = 0


		" ground truth walls of the environment "
		self.gndMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.gndImage = self.gndMapImage.load()
		gndDraw = ImageDraw.Draw(self.gndMapImage)
		

		walls = self.probe.getWalls()
		self.loadWalls(walls)
		for wall in walls:
			wallPoints = []
			for p in wall:
				pGrid = self.realToGrid(p)
				wallPoints.append((pGrid[0],pGrid[1]))
			
			gndDraw.line(wallPoints, fill=255)

		self.gndMapImage.save("mapGndGraph.png")

		self.occMap = OccupancyMap(self.probe, self, self.mapSize)
		self.boundMap = FreeSpaceBoundaryMap(self.mapSize)
		self.frontierMap = FrontierMap(self, self.mapSize)
		self.voronoiMap = VoronoiMap(self, self.mapSize)
		
		" whether or not the data needs to be synchronized "
		self.isDirty = False

		self.colors = []
		for i in range(100):
			self.colors.append((random(),random(),random()))

		if False:
			for k in range(1):
				pylab.clf()
				
				#[-4.0, 0.2] ,[-14.0, 0.2]
				#[-4.0, -0.2],[-14.0, -0.2]
				
				#xP = [-4.,-14.0]
				#yP = [0.2, 0.2]
				#pylab.plot(xP, yP, color='k')
				#yP = [-0.2, -0.2]
				#pylab.plot(xP, yP, color='k')
	
				self.plotEnv()
	
				theta = k*0.05
				 
				for i in range(self.probe.numSegs):
					
		
					pose = self.probe.getActualSegPose(i, phi = theta)
					#pose = self.probe.getActualJointPose(i, phi = theta)
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					
					p1 = [xTotal + 0.5*segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal - 0.5*segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal - 0.5*segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1], p2[1], p3[1], p4[1], p1[1]]
	
					pylab.plot(xP,yP, color='b')
					
				#for i in range(-1,20):
				for i in range(-1,self.probe.numSegs-1):
					
					pose = self.probe.getActualJointPose(i, phi = theta)
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1], p2[1], p3[1], p4[1], p1[1]]
		
					if i == 19:
						pylab.plot(xP,yP, color='r', linewidth='2')
						pylab.scatter([pose[0]], [pose[1]], color='r', linewidth=8)
					else:
						pylab.plot(xP,yP, color='r')				
						pylab.scatter([pose[0]], [pose[1]], color='r', linewidth=1)
	
				rootPose = self.probe.getActualJointPose(19)
	
				for i in range(-1,self.probe.numSegs-1):
				#for i in range(-1,20):
					
					#pose = self.probe.getActualJointPose(i, phi = theta)
					pose = self.probe.getJointPose(rootPose, 19, i)
					
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1], p2[1], p3[1], p4[1], p1[1]]
	
					if i == 19:
						pylab.plot(xP,yP, color='k', linewidth='2')
						pylab.scatter([pose[0]], [pose[1]], color='k', linewidth=4)
					else:
						pylab.plot(xP,yP, color='k')
						pylab.scatter([pose[0]], [pose[1]], color='k', linewidth=1)
	
				for i in range(-1,self.probe.numSegs-1):
					
					#pose = self.probe.getJointWRTJointPose([0.0,0.0,0.0], 19, i)
					pose = self.probe.getJointPose([0.0,0.0,0.0], 19, i)
					
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1]+3.0, p2[1]+3.0, p3[1]+3.0, p4[1]+3.0, p1[1]+3.0]
	
					if i == 19:
						pylab.plot(xP,yP, color='g', linewidth='2')
						pylab.scatter([pose[0]], [pose[1]+3.0], color='g', linewidth=4)
					else:
						pylab.plot(xP,yP, color='k')
						pylab.scatter([pose[0]], [pose[1]+3.0], color='g', linewidth=1)
									
				pylab.xlim(-4,4)
				pylab.ylim(-4,4)
				pylab.savefig("testplot%04u.png" % k)
				pylab.show()


	def drawConstraints(self, id = []):
		
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())
	
		cornerTrans = []

		pylab.clf()
		for i in range(self.numNodes):

			if self.nodeHash[i].isBowtie:			
				hull = computeBareHull(self.nodeHash[i], sweep = False, static = True)
				hull.append(hull[0])
			else:
				hull = computeBareHull(self.nodeHash[i], sweep = False)
				hull.append(hull[0])

			#hull = computeBareHull(self.nodeHash[i], sweep = False)
			#hull.append(hull[0])

			node1 = self.nodeHash[i]
			currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	

			hull_trans = []
			for p in hull:
				hull_trans.append(gen_icp.dispOffset(p, currPose))	
			
			" extract corner points for each node "
			cornerCandidates = node1.cornerCandidates
			cornerPoints = []
			for cand in cornerCandidates:
				cornerPoints.append(cand[0])
			
			for p in cornerPoints:
				cornerTrans.append(gen_icp.dispOffset(p, currPose))		
			
			if i == (self.numNodes - 1):
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

		for i in range(len(self.poseGraph.cornerBins)):
			bin = self.poseGraph.cornerBins[i]
			clusterTrans = []
			for item in bin:
				nodeID, cornerID = item
				currPose = self.nodeHash[nodeID].getGlobalGPACPose()
				cornerP = self.nodeHash[nodeID].cornerCandidates[cornerID][0]
				clusterTrans.append(gen_icp.dispOffset(cornerP, currPose))

			if len(clusterTrans) > 0:
				xP = []
				yP = []
				for p in clusterTrans:
					xP.append(p[0])
					yP.append(p[1])
				pylab.scatter(xP,yP, color=self.colors[i], linewidth=1, zorder=10)	

		self.drawWalls()


			
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)

	def drawWalls(self):
		for wall in self.walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = wall[i]
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g')

	def localizePose(self):
		
		self.foreNode.saveMap()
		self.backNode.saveMap()
	
		
		self.poseGraph.loadNewNode(self.foreNode)
		self.poseGraph.loadNewNode(self.backNode)
		

	def __getattr__(self, name):

		if name == "numNodes":
			return self.poseGraph.numNodes

		if name == "nodeHash":
			return self.poseGraph.nodeHash

	def restoreSeries(self, dirName, num_poses):

		self.poseGraph.restoreState(dirName, num_poses)

		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", i			
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			self.currNode.readFromFile2(dirName, i)			
			self.currNode.setGPACPose(self.poseGraph.nodePoses[i])
			self.poseGraph.nodeHash[i] = self.currNode


								
		
		foreNode = LocalNode(self.probe, self.contacts, num_poses, 19, PIXELSIZE)
		foreNode.readFromFile2(dirName, num_poses)
		backNode = LocalNode(self.probe, self.contacts, num_poses+1, 19, PIXELSIZE)
		backNode.readFromFile2(dirName, num_poses+1)
		self.poseGraph.insertPose(foreNode, backNode, initLocation = foreNode.getEstPose())
		
		self.isDirty = True
		
		self.synch()
		
		self.saveMap()

	def loadSeries(self, dirName, num_poses):
	
		numNodes = 0
		self.currNode = 0
			
		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile(dirName, i)
			self.poseGraph.loadNewNode(currNode)
			self.drawConstraints(i)

		self.isDirty = True

		self.synch()
		self.saveMap()

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
		self.poseGraph.loadNodes(nodeHash)
		
		" add the constraint edges to the graph "
		self.poseGraph.loadConstraints(dirName, num_poses-1)

		self.isDirty = True

		#self.drawConstraints()

		" save the maps "		
		#self.synch()
		#self.saveMap()

	def newInPlaceNode(self, faceDir, travelDir):
		
		#if self.numNodes > 0:
		#	
		#	if self.currNode.isDirty():
		#		self.currNode.synch()
		#
		#	self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = faceDir, travelDir = travelDir)
		
		self.poseGraph.addInPlaceNode(self.currNode)

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
		
		#self.poseGraph.addNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
		
		#if self.poseGraph.numNodes > 100:
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
		
		#self.poseGraph.pairLastTwo()

	def insertPose(self, estPose, direction):


		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = True, travelDir = direction)		
		self.foreNode = self.currNode
		self.forceUpdate(True)		
		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes+1, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = False, travelDir = direction)		
		self.backNode = self.currNode
		self.forceUpdate(False)		

		self.poseGraph.insertPose(self.foreNode, self.backNode, initLocation = estPose)
		#self.poseGraph.addInitNode(self.foreNode, estPose)
		#self.poseGraph.loadNewNode(self.foreNode)
		#self.poseGraph.loadNewNode(self.backNode)

		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)	

	def synch(self):
		
		self.currNode.synch()
		
		print "LocalNode Synched"

		self.occMap.update()
		
		print "Occupancy Map Updated"
		
		self.boundMap.update(self.occMap)

		print "Boundary Map Updated"
		
		" 1. for each node, copy local image "
		" 2. for each point in image, transform appropriately and plot to main image "
		self.obstMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.obstImage = self.obstMapImage.load()

		for i in range(self.numNodes):
			localNode = self.nodeHash[i]
			localObstMap = localNode.getObstacleMap()
			localMap = localObstMap.getMap()
			localImage = localMap.load()
			
			mapSize = localMap.size
		
			for j in range(mapSize[0]):
				for k in range(mapSize[1]):
					if localImage[j, k] == 255:

						pnt = localObstMap.gridToReal([j, k])

						pnt = localNode.convertLocalToGlobal(pnt)

						indexX, indexY = self.realToGrid(pnt)

						self.obstImage[indexX, indexY] = 255
							
		print "Obstacle Map Updated"
		
		self.frontierMap.update()
		
		print "Frontier Map Updated"

		self.voronoiMap.update()

		print "Voronoi Map Updated"
		
		self.navRoadMap = NavRoadMap(self.mapSize, self.probe, self.voronoiMap.getGraph(), localNode = self.currNode)

		print "Navigation Road Map Updated"

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
		
		newPoint = self.poseGraph.getNearestPathPoint(originPoint)
		
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
		
		path = self.poseGraph.computeNavigationPath(startPose, endPose)

		return path

	def computePath(self, currPose, frontierPoint):
		if self.isDirty:
			self.synch()

		path = self.navRoadMap.computePath(currPose, frontierPoint)
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

		
	def computeHeadPath(self, currPose, frontierPoint, exploreRoot):
		if self.isDirty:
			self.synch()

		#self.poseGraph.computeNavigationPath(currPose, frontierPoint)

		vals = self.navRoadMap.computeHeadPath(currPose, frontierPoint, exploreRoot)
		return vals

	def selectNextFrontier(self):
		if self.isDirty:
			self.synch()
		return self.frontierMap.selectNextFrontier()

	def selectNextDestination(self):
		if self.isDirty:
			self.synch()
		
		" get the termination point and orientation of each path "
		terms = self.poseGraph.paths.getPathTerms()
		termsVisited = self.poseGraph.paths.getPathTermsVisited()
		
		print "Terms:", terms
		print "TermsVisited:", termsVisited
		
		for k, term in terms.iteritems():
			#term = terms[k]
			
			if not termsVisited[k]:
				#self.poseGraph.paths.pathTermVisited(k)
				print "selecting term", k
				return term,k
		
		" if all terms visited, reset and go to root "
		self.poseGraph.paths.resetTerms()
		
		return self.poseGraph.paths.rootPoint, -1
		
		#for term in terms:
		#	
		#	if self.frontierMap.isFrontier(loc = term):
		#		return term
		
		print "all terms reached, selecting frontier point"
		return self.frontierMap.selectNextFrontier(), -1

	def isFrontier(self):
		if self.isDirty:
			self.synch()
		return self.frontierMap.isFrontier()

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

		self.poseGraph.saveState()
		
	def saveMap(self):
		
		if self.isDirty:
			self.synch()

		" save the global maps to file "

		print "saving global occupancy map"
		
		self.occMap.saveMap()

		print "saving global boundary map"
		
		self.boundMap.saveMap()
		
		print "saving global obstacle map"

		self.obstMapImage.save("mapObstGraph%04u.png" % self.saveCount)	

		print "saving global frontier map"

		self.frontierMap.saveMap()			
		
		print "saving global voronoi map"

		self.voronoiMap.saveMap()			

		print "building global navigation road map"

		self.navRoadMap = NavRoadMap(self.mapSize, self.probe, self.voronoiMap.getGraph(), localNode = self.currNode)

		print "done building global maps"
			
		self.saveCount += 1	

	def getCurrentNode(self):
		return self.currNode
	
	" FIXME:  Perhaps this should be moved to the PokeWalls behavior since this the only thing that requires it"
	def keepStablePose(self):
		self.stablePose.update()
			
	def realToGrid(self, point):
		indexX = int(floor(point[0] * self.divPix)) + self.numPixel / 2 + 1
		indexY = int(floor(point[1] * self.divPix)) + self.numPixel / 2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0, (j - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0]
		return point		
	
