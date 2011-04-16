
import Image
import pylab
import scipy.optimize
from numpy import *
from subprocess import *
from random import *
from math import floor, cos, sin, sqrt, asin, acos, pi
from numpy import array, dot
from copy import deepcopy, copy
from functions import *

class LocalOccMap:

	#def __init__(self, probe, contacts, nodeID):
	def __init__(self, localNode):

		self.localNode = localNode

		self.prevRects = []
		
		self.pixelSize = self.localNode.pixelSize
		self.mapSize = self.localNode.mapSize
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		self.prevJoints = [0.0,0.0,0.0,0.0]
		self.rootNeighbors = [18,19,20,21]

		self.nodeID = self.localNode.getNodeID()
		
		self.probe = self.localNode.getProbe()
		self.contacts = self.localNode.getContacts()
		self.rootNode = self.localNode.getRootNode()
		#self.rootPose = self.localNode.getRootPose()
		self.numJoints = self.probe.numSegs-1
		self.changed = False

		self.gnd_rectangles = []
		#self.gndInitRoot = self.probe.getActualJointPose(self.rootNode)
		self.gndMapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.gndImage = self.gndMapImage.load()
		#da = self.gndInitRoot[2]	
		#self.gndT = array([[cos(da), -sin(da), 0.0],[sin(da), cos(da), 0.0],[0.0,0.0,1.0]])

		self.stable_rectangles = []
		self.stableMapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.stableImage = self.stableMapImage.load()


		self.rectangles = []

		self.gndFileName = "gndLocalOccMap%03u" % self.nodeID + "_%04u.png"
		self.fileName = "localOccMap%03u" % self.nodeID + "_%04u.png"
		self.stableFileName = "stableLocalOccMap%03u" % self.nodeID + "_%04u.png"
		self.saveCount = 0
		self.mapImage = 0

		#self.costFile = open("costFile%03u.txt" % self.nodeID, 'w')
		self.hullFile = "convexHull%03u_%04u.png"
		self.hullCount = 0
		
		self.T = array([[1,0,0],
			[0,1,0],
			[0,0,1]], float)

	def saveToFile(self, filename):
		
		#1. first line has the offset
		#2. next is the printout of the rectangles
		val = repr(self.rectangles)
		f = open(filename, 'w')
		f.write(val)
		f.close()

	def readFromFile(self, dirName):
		
		#self.fileName = dirName + "/stableLocalOccMap%03u" % self.nodeID + "_%04u.png"
		self.fileName = dirName + "/localOccMap%03u" % self.nodeID + "_%04u.png"
		
		self.mapImage = Image.open(self.fileName % 0)
		self.image = self.mapImage.load()
		
		self.changed = True
		
	def getRectangles(self):
		return deepcopy(self.rectangles)
	
	def setRectangles(self, rectangles):
		self.rectangles = rectangles
			
	def findClosestPoint(self, point, cHull):
		"""
		TWO CANDIDATES
		# 1. check each of the path points
		# 2. check orthogonal distance to each edge if less than both terminating points
		"""
		
		"  point-to-point distances "
		pointDist = []
		for p in cHull:
			dist = sqrt((point[0]-p[0])**2 + (point[1]-p[1])**2)
			pointDist.append(dist)

		" orthogonal distance to each segment "
		orthoDist = []
		closePoints = []

		for i in range(len(cHull)-1):
			p0 = copy(cHull[i])
			p1 = copy(cHull[i+1])

			# rotate each edge to 0 degrees, and compare point distance to the horizontal
			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			mag = sqrt(vec[0]**2 + vec[1]**2)
			if mag == 0:
				orthoDist.append( 1e100 )
				closePoints.append(copy(cHull[i]))
				print "degenerate"
				continue

			eAngle = acos(vec[0]/mag)
			if asin(vec[1]/mag) < 0:
				eAngle = -eAngle

			" orthogonal point "
			oPoint = [point[0]-p0[0], point[1]-p0[1]]
			oPoint = [oPoint[0] * cos(eAngle) + oPoint[1] * sin(eAngle), -oPoint[0] * sin(eAngle) + oPoint[1] * cos(eAngle)]

			dist = abs(oPoint[1])
			if oPoint[0] < 0.0:
				orthoDist.append( pointDist[i])
				closePoints.append(copy(cHull[i]))

			elif oPoint[0] > mag:
				orthoDist.append(pointDist[i+1])
				closePoints.append(copy(cHull[i+1]))

			else:
				orthoDist.append(dist)
				nearPoint = [oPoint[0],0.0]
				nearPoint = [nearPoint[0] * cos(eAngle) - nearPoint[1] * sin(eAngle), nearPoint[0] * sin(eAngle) + nearPoint[1] * cos(eAngle)]
				closePoints.append(copy(nearPoint))

		minDist = 1e100
		minIndex = 0
		for i in range(len(orthoDist)):
			if orthoDist[i] < minDist:
				minDist = orthoDist[i]
				minIndex = i
				
		return minDist, copy(closePoints[minIndex])
	
	def update(self, isForward = True):
		
		"""
		1. add a ground truth map to the LocalOccMap class
		2. requires setting of initial pose
		2. requires setting of ground pose at each iteration
		3. requires plotting of each ground rectangle to map, but skipping all corrections
		
		"""
		
		# compute ground free space first
		gndRects = self.gndComputeFreeSpace()
		gndRects.sort()
		self.gnd_rectangles.append(gndRects)
		
		#stableRects = self.computeStabilizedFreeSpace()
		#stableRects.sort()
		#self.stable_rectangles.append(stableRects)
		
		
		# compute free space rectangles
		#freeRects = self.computeFreeSpace(isForward)
		#freeRects.sort()

		freeRects = self.computeStabilizedFreeSpace()
		freeRects.sort()
		#self.stable_rectangles.append(freeRects)


		jointErr = 0.0
		currJoints = []
		for j in self.rootNeighbors:
			currJoints.append(self.probe.getServo(j))

		#for i in range(len(currJoints)):
		#	jointErr += abs(currJoints[i]-self.prevJoints[i])
		
		#newRects = []
		#for seg, rect, isActive in freeRects:	
		#	newRect = []	
		#	for p in rect:
		#		pVec = array([[p[0]],[p[1]],[1.0]])
		#		result = dot(self.T,pVec)
		#		newP = [result[0,0],result[1,0]]
		#		newRect.append(newP)
		#	newRects.append([seg,newRect, isActive])

		#self.prevRects = deepcopy(newRects)		
		#self.prevJoints = currJoints
		#self.rectangles.append(newRects)
		
		self.prevRects = deepcopy(freeRects)		
		self.prevJoints = currJoints

		self.rectangles.append(freeRects)

	def gndComputeFreeSpace(self):
							
		" Old motion shadow-based obstacle detection, not used "
		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		actualConfig = []
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		
		currGnd = self.probe.getActualJointPose(self.rootNode)

		globalVec = array([[currGnd[0]], [currGnd[1]]])
		tempVec = dot(self.localNode.gndBackR, globalVec)
		tempVec[0, 0] -= self.localNode.gndDist
		transVec = dot(self.localNode.gndForeR, tempVec)
		finalVec = dot(self.localNode.gndR, transVec)
		
		xTotal = finalVec[0, 0]
		zTotal = finalVec[1, 0]
		totalAngle = normalizeAngle(currGnd[2] - self.localNode.gndPose[2])

		#print xTotal, zTotal, totalAngle
		
		joints = range(-1,self.rootNode)
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

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])


		joints = range(self.rootNode, self.probe.numSegs-1)

		#xTotal = newOrigin[0]
		#zTotal = newOrigin[1]
		#totalAngle = newOrigin[2]
		xTotal = finalVec[0, 0]
		zTotal = finalVec[1, 0]
		totalAngle = normalizeAngle(currGnd[2] - self.localNode.gndPose[2])

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])
			
			if i < self.probe.numSegs-2:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)

		actualConfig.sort()
		return actualConfig	
	
	def computeFreeSpace(self, isForward = True):
							
		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		actualConfig = []
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth

		#xTotal = self.localNode.rootPose[0]
		#zTotal = self.localNode.rootPose[1]
		#totalAngle = self.localNode.rootPose[2]

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		joints = range(-1,self.rootNode)
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

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])

		joints = range(self.rootNode, self.probe.numSegs-1)

		#xTotal = self.localNode.rootPose[0]
		#zTotal = self.localNode.rootPose[1]
		#totalAngle = self.localNode.rootPose[2]
		
		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])
			
			if i < self.probe.numSegs-2:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)

		actualConfig.sort()
		return actualConfig		

	def computeStabilizedFreeSpace(self, isForward = True):
							
		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		actualConfig = []
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth

		xTotal = self.localNode.rootPose[0]
		zTotal = self.localNode.rootPose[1]
		totalAngle = self.localNode.rootPose[2]

		#xTotal = 0.0
		#zTotal = 0.0
		#totalAngle = 0.0

		joints = range(-1,self.rootNode)
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

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])

		joints = range(self.rootNode, self.probe.numSegs-1)

		xTotal = self.localNode.rootPose[0]
		zTotal = self.localNode.rootPose[1]
		totalAngle = self.localNode.rootPose[2]
		
		#xTotal = 0.0
		#zTotal = 0.0
		#totalAngle = 0.0

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])
			
			if i < self.probe.numSegs-2:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)

		actualConfig.sort()
		return actualConfig		
	
	def fillOpen(self, polygon, image):

		# vertices of the 4-sided convex polygon
		polyX = [polygon[0][0], polygon[1][0], polygon[2][0], polygon[3][0], polygon[0][0]]
		polyY = [polygon[0][1], polygon[1][1], polygon[2][1], polygon[3][1], polygon[0][1]]

		# bounding box of polygon
		maxX = -10000.0
		minX = 10000.0
		maxY = -10000.0
		minY = 10000.0
		for i in range(4):
			if polyX[i] > maxX:
				maxX = polyX[i]
			if polyX[i] < minX:
				minX = polyX[i]
			if polyY[i] > maxY:
				maxY = polyY[i]
			if polyY[i] < minY:
				minY = polyY[i]

		# set boundaries of map grid
		highIndexX, highIndexY = self.realToGrid([maxX, maxY])
		lowIndexX, lowIndexY = self.realToGrid([minX, minY])

		# fill map grid cell if occupied by arm
		for i in range(lowIndexX,highIndexX+1,1):
			for j in range(lowIndexY,highIndexY+1,1):
				point = self.gridToReal([i,j])
				if IsContained(polygon, point):
					# free space overwrites everything
					image[i,j] = 255

	def realToGrid(self, point):
		indexX = int(floor(point[0]*self.divPix)) + self.numPixel/2 + 1
		indexY = int(floor(point[1]*self.divPix)) + self.numPixel/2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0, (j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]
		return point

	def buildMap(self):
		
		# 1. create the map image
		if self.mapImage == 0:
			self.mapImage = Image.new('L', (self.numPixel,self.numPixel),127)
			self.image = self.mapImage.load()

		# 2. plot just the latest rectangles to the image
		for snapshot in self.rectangles:
			for rect in snapshot:
				self.fillOpen(rect[1], self.image)
		
		for snapshot in self.gnd_rectangles:
			for rect in snapshot:
				self.fillOpen(rect[1], self.gndMapImage.load())

		#for snapshot in self.stable_rectangles:
		#	for rect in snapshot:
		#		self.fillOpen(rect[1], self.stableMapImage.load())
	
		" so outside maps can trigger changes "
		self.changed = True

		" save the rectangles, and clear the list for new ones "
		self.rectangles = []
		self.gnd_rectangles = []
		self.stable_rectangles = []
	
	def saveMap(self, fileName = ""):
		
		self.buildMap()
		
		if fileName != "":
			self.mapImage.save(fileName)
		else:
			print "saving as", self.fileName % self.saveCount
			self.mapImage.save(self.fileName % self.saveCount)	
			self.gndMapImage.save(self.gndFileName % self.saveCount)
			#self.stableMapImage.save(self.stableFileName % self.saveCount)
			self.saveCount += 1
		
	def getMap(self):
		return self.mapImage
	
	def getGndMap(self):
		return self.gndMapImage
