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

class ObstacleMap(Map):

	def __init__(self, probe, boundMap, contacts):
		Map.__init__(self)

		self.contacts = contacts
		self.probe = probe
		self.fileName = "obstacleMap%04u.png"

		self.poly = [[0.0,0.0],[0.0,0.0]]
		self.update(boundMap)

	def saveMap(self):

		tempImage = self.mapImage.copy()
		draw = ImageDraw.Draw(tempImage)
		
		pixPoly = []
		for p in self.poly:
			x, y = self.realToGrid(p)
			pixPoly.append([x,y])
		
		num = len(self.poly)
		
		for i in range(num):
			draw.line((pixPoly[i][0], pixPoly[i][1], pixPoly[(i+1) % num][0], pixPoly[(i+1) % num][1]), fill = 255)
		
		tempImage.save(self.fileName % self.saveCount)	
		self.saveCount += 1
		
	def setFrontierMap(self, frontMap):
		self.frontMap = frontMap

	def update(self, boundMap):
		
		self.occMap = boundMap.occMap
		self.boundMap = boundMap
		
		self.xMin = boundMap.xMin
		self.xMax = boundMap.xMax
		self.yMin = boundMap.yMin
		self.yMax = boundMap.yMax

		self.points = boundMap.getBoundaryPoints()
		
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
		angle, x, y = self.getProbeProfile(isFront)
		self.setContacts(angle, x, y)

	def getProbeProfile(self, isFront):
		
		if isFront:
			# get actual configuration from rootNode reference
			origin = self.contacts.getClosestPose(0)
			#origin = self.probe.getActualJointPose(0)
			xTotal = origin[0] 
			zTotal = origin[1] 
			totalAngle = origin[2]
			
			totalAngle = totalAngle + self.probe.getServo(0)
			totalAngle = normalizeAngle(totalAngle)
			
			return normalizeAngle(totalAngle + pi), xTotal, zTotal
		else:
			# get actual configuration from rootNode reference
			#origin = self.probe.getActualJointPose(38)
			origin = self.contacts.getClosestPose(38)
			xTotal = origin[0] 
			zTotal = origin[1] 
			totalAngle = origin[2]
			
			totalAngle = normalizeAngle(totalAngle)
			
			return totalAngle, xTotal, zTotal
	
	def setContacts(self, angle, x, y):
		
		# find the obstacle points within the cone projected from tip
		segLength = self.probe.segLength
		xTip = x + segLength*cos(angle)
		yTip = y + segLength*sin(angle)

		# recompute the boundary points just prior to marking obstacles
		self.occMap.update()
		self.boundMap.update(self.occMap)
		self.update(self.boundMap)
			
		# if there is no collision, mark obstacles
		#if not self.checkBoundaryCollision(polygon):
		#	self.markObstacle(xTip, yTip, angle)
		
		# contact detection is good enough that I don't need to check for bad cases
		# just mark the obstacle pixels
		# if boundary points are within cone projected from tip, then they are contact points
		self.markObstacle(xTip, yTip, angle)
			
		self.update(self.boundMap)
		self.frontMap.update(self.boundMap, self)
		
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
		polyX = [polygon[0][0], polygon[1][0], polygon[2][0], polygon[3][0], polygon[0][0]]
		polyY = [polygon[0][1], polygon[1][1], polygon[2][1], polygon[3][1], polygon[0][1]]

		self.poly = copy(polygon)
		
		
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
				xP, yP = self.gridToReal([i, j])
				
				result = IsContained(polygon, [xP, yP])
				
				if result:
					if self.boundMap.image[i,j] == 255:
						# save this as an obstacle point
						self.image[i,j] = 255
					
				#dist = sqrt((x-xP)**2 + (y-yP)**2)
				#if dist <= radius:
				#	if self.boundMap.image[i,j] == 255:
					
						# save this as an obstacle point
						#self.image[i,j] = 255




		
		
		
		