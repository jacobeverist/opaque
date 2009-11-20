import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from Map import Map
import Image


class OccupancyMap(Map):

	def __init__(self, probe, rootNode, mapSize = MAPSIZE):
		Map.__init__(self, mapSize)

		self.probe = probe
		self.rootNode = rootNode
		self.numJoints = self.probe.numSegs-1

		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0

		self.fileName = "occupancyMap%04u.png"

	def resetMap(self):
		self.mapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.image = self.mapImage.load()
		
		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0
	
	# load a custom image into the class

	def setRootPose(self, pose):
		self.rootPose = pose
		
	def update(self, isForward = True):
	
		# compute free space rectangles
		freeRects = self.computeFreeSpace(isForward)
		for rectangle in freeRects:
			self.fillOpen(rectangle)

		# compute configuration shadow rectangles
		#shadowRects = self.computeConfigShadow()
		#for rectangle in shadowRects:
		#	self.fillShadow(rectangle)
		
	def computeFreeSpace(self, isForward = True):
		actualConfig = []

		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1

		# get actual configuration from rootNode reference
		#origin = self.probe.getActualJointPose(self.rootNode)
		#xTotal = origin[0] 
		#zTotal = origin[1] 
		#totalAngle = origin[2] 

		xTotal = self.rootPose[0]
		zTotal = self.rootPose[1]
		totalAngle = self.rootPose[2]

		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		if isForward:
			joints = range(0,self.rootNode+1)
			joints.reverse()
	
			for i in joints:
	
				totalAngle = totalAngle + self.probe.getServo(i)
				totalAngle = normalizeAngle(totalAngle)
	
				p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
				p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
				#p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
				#p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
				p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
				p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
				xTotal = xTotal - segLength*cos(totalAngle)
				zTotal = zTotal - segLength*sin(totalAngle)
	
				actualConfig.append([p4,p3,p2,p1])
				#actualConfig.append([p1,p2,p3,p4])
			
			return actualConfig
		
		else:
			joints = range(self.rootNode, self.probe.numSegs-1)
	
			for i in joints:
				xTotal = xTotal + segLength*cos(totalAngle)
				zTotal = zTotal + segLength*sin(totalAngle)
			
				p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
				p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
				#p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
				#p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
				p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
				p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
	
				actualConfig.append([p4,p3,p2,p1])
				#actualConfig.append([p1,p2,p3,p4])
				
				if i < self.probe.numSegs-2:
					totalAngle = totalAngle - self.probe.getServo(i+1)
					totalAngle = normalizeAngle(totalAngle)
			
			return actualConfig


	def computeConfigShadow(self):

		actAngles = []
		desAngles = []
		for i in range(0,self.numJoints):
			actAngles.append(self.probe.getServo(i))
			desAngles.append(self.probe.getServoCmd(i))

		segLength = self.probe.segLength
		segWidth = self.probe.segWidth

		joints = range(0,self.rootNode+1)
		joints.reverse()

		rectangles = []

		for i in joints:
			origin = self.probe.getActualJointPose(i)

			desAng = desAngles[i]
			actAng = actAngles[i]
			jointErr = desAng - actAng

			angleInc = []
			if jointErr > 0:
				angleInc = arange(actAng,desAng+INC_VAL, INC_VAL)
			else:
				angleInc = arange(desAng,actAng+INC_VAL, INC_VAL)

			distalJoints = range(max(0,i-DISTAL_RANGE),i)
			distalJoints.reverse()

			# compute the error increments, and the distal joint positions
			for angleVal in angleInc:

				xTotal = origin[0] 
				zTotal = origin[1] 
				totalAngle = origin[2] 

				totalAngle = totalAngle + angleVal

				p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
				p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
				#p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
				#p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
				p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
				p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
				xTotal = xTotal - segLength*cos(totalAngle)
				zTotal = zTotal - segLength*sin(totalAngle)

				totalAngle = normalizeAngle(totalAngle)

				#rect1 = [p1,p2,p3,p4]
				rect1 = [p4,p3,p2,p1]
				rectangles.append(rect1)

				# now compute distal joints with actual angles
				for j in distalJoints:

					totalAngle = totalAngle + actAngles[j]

					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					#p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					#p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

					xTotal = xTotal - segLength*cos(totalAngle)
					zTotal = zTotal - segLength*sin(totalAngle)

					totalAngle = normalizeAngle(totalAngle)

					#rect2 = [p1,p2,p3,p4]
					rect2 = [p4,p3,p2,p1]
					rectangles.append(rect2)

		return rectangles
	
	def fillOpen(self, polygon):

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

		# fill map grid cell if occupied by arm
		indices = []

		for i in range(lowIndexX,highIndexX+1,1):
			for j in range(lowIndexY,highIndexY+1,1):
				if IsContained(polygon, [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0,(j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]):
					# free space overwrites everything
					self.image[i,j] = 255


	def fillShadow(self, polygon):

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
				if IsContained(polygon, [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0,(j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]):

					# if an obstacle, don't overwrite free space
					val = self.image[i,j]
					if val != 255:
						val = 0
						self.image[i,j] = val
						
						