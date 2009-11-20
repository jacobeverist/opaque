import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from Map import Map
import Image
import ImageDraw
import pylab

class LocalObstacleMap(Map):

	#def __init__(self, probe, boundMap, contacts):
	def __init__(self, localNode):
		Map.__init__(self)

		self.localNode = localNode
		self.contacts = self.localNode.getContacts()
		self.probe = self.localNode.getProbe()
		self.boundMap = self.localNode.getBoundMap()

		self.nodeID = self.localNode.getNodeID()
		self.fileName = "localObstacleMap%03u" % self.nodeID + "_%04u.png"

		self.pixelSize = self.boundMap.pixelSize
		self.mapSize = self.boundMap.mapSize
		self.numPixel = self.boundMap.numPixel
		self.halfPix = self.boundMap.halfPix
		self.divPix = self.boundMap.divPix
		self.resetMap()

		self.xMin = self.numPixel
		self.xMax = -1
		self.yMin = self.numPixel
		self.yMax = -1
			
		#self.poly = [[0.0,0.0],[0.0,0.0]]
		self.polygons = []
		self.oldPolygons = []
		self.update()

	def readFromFile(self):
		#self.nodeID = nodeID
		#fileName = "localObstacleMap%03u" % self.nodeID + "_0000.png"
		
		self.mapImage = Image.open("localObstacleMap%03u" % self.nodeID + ".png")
		self.loadImage(self.mapImage)
		
	def saveMap(self):

		self.plotPolygons()

		tempImage = self.mapImage.copy()
		
		"""
		draw = ImageDraw.Draw(tempImage)
		
		pixPoly = []
		for p in self.poly:
			x, y = self.realToGrid(p)
			pixPoly.append([x,y])
			
		
		num = len(self.poly)
		
		for i in range(num):
			draw.line((pixPoly[i][0], pixPoly[i][1], pixPoly[(i+1) % num][0], pixPoly[(i+1) % num][1]), fill = 255)
		"""
		tempImage.save(self.fileName % self.saveCount)	
		self.saveCount += 1
		
	def update(self):
		
		self.pixelSize = self.boundMap.pixelSize
		self.mapSize = self.boundMap.mapSize
		self.numPixel = self.boundMap.numPixel
		self.halfPix = self.boundMap.halfPix
		self.divPix = self.boundMap.divPix

		self.points = self.boundMap.getBoundaryPoints()
		
		self.plotPolygons()
		
		# clear out obstacle points that clearly reside within free space
		#self.clearFreeSpace()

		#self.mapImage = boundMap.copyImage()
		#self.image = self.mapImage.load()
	
	def clearFreeSpace(self):
		
		# clear obstacle points that are clearly within free space
		for i in range(self.xMin,self.xMax+1,1):
			for j in range(self.yMin,self.yMax+1,1):
				if self.image[i,j] == 255:
					if self.occMap.image[i,j] == 255:
						self.image[i,j] = 0

	def updateContact(self, isFront):

		
		" synchronize the maps "
		#self.localNode.synch()
		
		self.rootNode = self.localNode.getRootNode()
		#self.rootPose = self.localNode.getRootPose()
		self.rootPose = [0.0,0.0,0.0]		
		
		" find the obstacle points within the cone projected from tip "
		angle, x, y = self.getProbeProfile(isFront)

		" now mark the obstacles "
		self.markObstacle(x, y, angle)


		# recompute the boundary points just prior to marking obstacles
		#self.occMap.update()
		#self.boundMap.update(self.occMap)
		#self.update(self.boundMap)
			
		# if there is no collision, mark obstacles
		#if not self.checkBoundaryCollision(polygon):
		#	self.markObstacle(xTip, yTip, angle)
		
		# contact detection is good enough that I don't need to check for bad cases
		# just mark the obstacle pixels
		# if boundary points are within cone projected from tip, then they are contact points

			
		#self.occMap.saveMap()
		#self.boundMap.saveMap(self.poly)
		#self.saveMap()

	def getProbeProfile(self, isFront):

		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0
		
		if isFront:
			
			joints = range(-1,self.rootNode)
			joints.reverse()
	
			for i in joints:
	
				totalAngle = totalAngle + self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)

				xTotal = xTotal - segLength*cos(totalAngle)
				zTotal = zTotal - segLength*sin(totalAngle)
	
			# get actual configuration from rootNode reference
			#origin = self.contacts.getClosestPose(0)
			#xTotal = origin[0] 
			#zTotal = origin[1] 
			#totalAngle = origin[2]
			
			return normalizeAngle(totalAngle + pi), xTotal, zTotal
		
		else:
			
			joints = range(self.rootNode, self.probe.numSegs-1)
	
			xTotal = 0.0
			zTotal = 0.0
			totalAngle = 0.0
	
			for i in joints:
				xTotal = xTotal + segLength*cos(totalAngle)
				zTotal = zTotal + segLength*sin(totalAngle)
			
				if i < self.probe.numSegs-2:
					totalAngle = totalAngle - self.probe.getServo(i+1)
					totalAngle = normalizeAngle(totalAngle)

			# get actual configuration from rootNode reference
			#origin = self.contacts.getClosestPose(38)
			#xTotal = origin[0] 
			#zTotal = origin[1] 
			#totalAngle = origin[2]
			
			totalAngle = normalizeAngle(totalAngle)
			
			return totalAngle, xTotal, zTotal
	



	def markObstacle(self, x, y, angle):
		
		#radius = 0.05
		radius = 0.1
		polyX = [ x + radius, x - radius]
		polyY = [ y + radius, y - radius]
	
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth	
		#p1 = [x - 0.5*segWidth*sin(angle), y + 0.5*segWidth*cos(angle)]
		#p2 = [x + 0.5*segWidth*sin(angle), y - 0.5*segWidth*cos(angle)]
		#p3 = [x + segLength*cos(angle) + segWidth*sin(angle), y + segLength*sin(angle) - segWidth*cos(angle)]
		#p4 = [x + segLength*cos(angle) - segWidth*sin(angle), y + segLength*sin(angle) + segWidth*cos(angle)]
		p1 = [x - 2.0*segWidth*sin(angle), y + 2.0*segWidth*cos(angle)]
		p2 = [x + 2.0*segWidth*sin(angle), y - 2.0*segWidth*cos(angle)]
		p3 = [x + segLength*cos(angle) + 2.0*segWidth*sin(angle), y + segLength*sin(angle) - 2.0*segWidth*cos(angle)]
		p4 = [x + segLength*cos(angle) - 2.0*segWidth*sin(angle), y + segLength*sin(angle) + 2.0*segWidth*cos(angle)]
		#polygon = [p4,p3,p2,p1]
		polygon = [p1,p2,p3,p4]
		
		#print polygon
		
		self.polygons.append(copy(polygon))

	def plotPolygons(self):

		for polygon in self.polygons:
			
			polyX = [polygon[0][0], polygon[1][0], polygon[2][0], polygon[3][0], polygon[0][0]]
			polyY = [polygon[0][1], polygon[1][1], polygon[2][1], polygon[3][1], polygon[0][1]]
			
			# bounding box of circle
			maxX = -10000.0
			minX = 10000.0
			maxY = -10000.0
			minY = 10000.0
			for i in range(len(polyX)):
				if polyX[i] > maxX:
					maxX = polyX[i]
				if polyX[i] < minX:
					minX = polyX[i]
				if polyY[i] > maxY:
					maxY = polyY[i]
				if polyY[i] < minY:
					minY = polyY[i]

			# set boundaries of map grid
			highIndexX = int(floor(maxX*self.divPix)) + self.numPixel/2 + 1
			lowIndexX = int(floor(minX*self.divPix)) + self.numPixel/2 + 1
			highIndexY = int(floor(maxY*self.divPix)) + self.numPixel/2 + 1
			lowIndexY = int(floor(minY*self.divPix)) + self.numPixel/2 + 1
	
			# if needed, adjust the global bounding box
			if lowIndexX < self.xMin and lowIndexX >= 0:
				self.xMin = lowIndexX
	
			if highIndexX > self.xMax and highIndexX < self.numPixel:
				self.xMax = highIndexX
	
			if lowIndexY < self.yMin and lowIndexY >= 0:
				self.yMin = lowIndexY
	
			if highIndexY > self.yMax and highIndexY < self.numPixel:
				self.yMax = highIndexY
	
			for i in range(lowIndexX,highIndexX+1,1):
				for j in range(lowIndexY,highIndexY+1,1):
					xP, yP = self.gridToReal([i, j])
					
					result = IsContained(polygon, [xP, yP])
					
					if result:
						#print i, j
						if self.boundMap.image[i,j] == 255:
							# save this as an obstacle point
							self.image[i,j] = 255

		self.oldPolygons += self.polygons
		self.polygons = []

	def checkBoundaryCollision(self, polygon):
		
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
		highIndexX = int(floor(maxX*self.divPix)) + self.numPixel/2 + 1
		lowIndexX = int(floor(minX*self.divPix)) + self.numPixel/2 + 1
		highIndexY = int(floor(maxY*self.divPix)) + self.numPixel/2 + 1
		lowIndexY = int(floor(minY*self.divPix)) + self.numPixel/2 + 1

		# if needed, adjust the global bounding box
		if lowIndexX < self.xMin:
			self.xMin = lowIndexX

		if highIndexX > self.xMax:
			self.xMax = highIndexX

		if lowIndexY < self.yMin:
			self.yMin = lowIndexY

		if highIndexY > self.yMax:
			self.yMax = highIndexY

		for i in range(lowIndexX,highIndexX+1,1):
			for j in range(lowIndexY,highIndexY+1,1):
				
				if self.boundMap.image[i,j] == 255:
					xP, yP = self.gridToReal([i, j])
					result = IsContained(polygon, [xP, yP])
					
					if result:
						return True
				
		return False
						
	def getMap(self):
		return self.mapImage
