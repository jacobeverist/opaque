from PoseGraph import PoseGraph
from LocalNode import LocalNode
from VoronoiMap import VoronoiMap
from FrontierMap import FrontierMap
from OccupancyMap import OccupancyMap
from FreeSpaceBoundaryMap import FreeSpaceBoundaryMap
from NavRoadMap import NavRoadMap
from StablePose import StablePose

import pylab
from math import *

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
		for wall in walls:
			wallPoints = []
			for p in wall:
				pGrid = self.realToGrid(p)
				wallPoints.append((pGrid[0],pGrid[1]))
			
			gndDraw.line(wallPoints, fill=255)

		self.gndMapImage.save("mapGndGraph.png")

		self.occMap = OccupancyMap(self.probe, self, self.mapSize)
		self.boundMap = FreeSpaceBoundaryMap(self.occMap, self.mapSize)
		self.frontierMap = FrontierMap(self, self.mapSize)
		self.voronoiMap = VoronoiMap(self, self.mapSize)

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

	def localizeCurrentNode(self):
		self.poseGraph.localizeCurrentNode()

	def __getattr__(self, name):

		if name == "numNodes":
			return self.poseGraph.numNodes

		if name == "nodeHash":
			return self.poseGraph.nodeHash

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

		#self.drawConstraints()

		" save the maps "		
		#self.synch()
		#self.saveMap()

	def newInPlaceNode(self, direction):
		
		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, direction = direction)
		
		self.poseGraph.addInPlaceNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
	
	def newNode(self, stepDist, direction):

		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, direction = True)		
		
		self.poseGraph.addNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

	def initNode(self, estPose, direction):

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, direction = True)

		self.poseGraph.addInitNode(self.currNode, estPose)
		
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

	def correctPosture(self):
		self.currNode.updateCorrectPosture()

	def forceUpdate(self, isForward=True):
		self.stablePose.setDirection(isForward)
		self.currNode.update(isForward)
		
	def update(self, isForward=True):

		self.stablePose.setDirection(isForward)

		" if the configuration is a stable pose, the maximum displacement is not exceeded "
		if self.stablePose.isStable():
			self.currNode.update(isForward)
		else:
			pass
		
	def computeHeadPath(self, currPose, frontierPoint, exploreRoot):
		vals = self.navRoadMap.computeHeadPath(currPose, frontierPoint, exploreRoot)
		return vals

	def selectNextFrontier(self):
		return self.frontierMap.selectNextFrontier()

	def isFrontier(self):
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
		if self.currNode != 0:
			self.currNode.saveMap()
			
	def saveConstraints(self):

		self.poseGraph.saveConstraints()
		
	def saveMap(self):
		
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
