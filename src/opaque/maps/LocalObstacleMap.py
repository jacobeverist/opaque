
from Map import Map
import Image
import ImageDraw
import pylab
from math import cos, sin, floor, pi
from copy import copy


# determines signed area of 3 points (used for solving Point in Polygon problem)
def Area2(Ax,Ay,Bx,By,Cx,Cy):
	return (Bx - Ax) * (Cy - Ay) - (Cx - Ax)*(By - Ay)

# determines if point C is left of line segment AB
def LeftOn(Ax,Ay,Bx,By,Cx,Cy):
	return (Area2(Ax,Ay,Bx,By,Cx,Cy) >= 0)

# determine if point is located within a bounding box specified by vertices RectVert
def IsContained(RectVert, Point):
	for i in range(4):
		if not (LeftOn(RectVert[i%4][0],RectVert[i%4][1],
					   RectVert[(i+1)%4][0],RectVert[(i+1)%4][1], Point[0], Point[1])):
			return False
	return True

# this function converts the angle to its equivalent # in the range [-pi,pi]
def normalizeAngle(angle):
	
	while angle>pi:
		angle=angle-2*pi
	
	while angle<=-pi:
		angle=angle+2*pi
	
	return angle 

class LocalObstacleMap(Map):

	def __init__(self, localNode):
		#Map.__init__(self)

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
			
		self.polygons = []
		self.oldPolygons = []
		self.update()

	def readFromFile(self):

		self.fileName = "localObstacleMap%03u" % self.nodeID + "_%04u.png"
		
		self.mapImage = Image.open(self.fileName % 0)
					
		self.loadImage(self.mapImage)
		
	def saveMap(self):

		self.plotPolygons()

		tempImage = self.mapImage.copy()
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

	def getProbeProfile(self, isFront):

		segLength = self.probe.segLength

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
			
			totalAngle = normalizeAngle(totalAngle)
			
			return totalAngle, xTotal, zTotal


	def markObstacle(self, x, y, angle):
		
		" mark obstacle to be 2 widths and 1 length "	
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth	
		p1 = [x - 2.0*segWidth*sin(angle), y + 2.0*segWidth*cos(angle)]
		p2 = [x + 2.0*segWidth*sin(angle), y - 2.0*segWidth*cos(angle)]
		p3 = [x + segLength*cos(angle) + 2.0*segWidth*sin(angle), y + segLength*sin(angle) - 2.0*segWidth*cos(angle)]
		p4 = [x + segLength*cos(angle) - 2.0*segWidth*sin(angle), y + segLength*sin(angle) + 2.0*segWidth*cos(angle)]
		polygon = [p1,p2,p3,p4]
				
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

	def getMap(self):
		return self.mapImage
	
