import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from Map import Map
import Image
import pylab
import ImageDraw


class FrontierMap(Map):

	def __init__(self, probe, boundMap, obstacleMap):
		Map.__init__(self)

		self.probe = probe
		self.fileName = "frontierMap%04u.png"
		
		self.update(boundMap, obstacleMap)
		
		self.checkSpline = 0
	
	def saveMap(self):

		tempImage = self.mapImage.copy()
		draw = ImageDraw.Draw(tempImage)
		
		pixPoly = []
		for p in self.obstacleMap.poly:
			x, y = self.realToGrid(p)
			pixPoly.append([x,y])
			
		num = len(self.obstacleMap.poly)
		
		for i in range(num):
			draw.line((pixPoly[i][0], pixPoly[i][1], pixPoly[(i+1) % num][0], pixPoly[(i+1) % num][1]), fill = 255)
		
		tempImage.save(self.fileName % self.saveCount)	
		self.saveCount += 1
		
	def isFrontier(self):
		
		centerNodes = [2, 6, 10, 14, 18, 22, 26, 30, 34]
		centerNodes.reverse()
		points = []
		for n in centerNodes:
			pnt = self.probe.getActualJointPose(n)
			points.append([pnt[0],pnt[1]])
			
		spline = SplineFit(points)
		self.checkSpline = spline
		
		vec = spline.getUVector(0.9)
		tip = spline.getU(1.0)
		
		#vecDist = 1.0
		#vecDist = 0.7
		vecDist = 0.75
		x = tip[0] + vec[0]*vecDist
		y = tip[1] + vec[1]*vecDist

		#x = point[0]
		#y = point[1]
		radius = 0.6
		
		polyX = [ x + radius, x - radius]
		polyY = [ y + radius, y - radius]

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
			
			
		upoints = arange(0,1.0,0.01)
		splPoints = self.checkSpline.getUSet(upoints)

		#Image.new('L', (self.numPixel,self.numPixel),0)
		im = self.mapImage.copy()
		imPix = im.load()
		draw = ImageDraw.Draw(im)
		draw.ellipse((lowIndexX, lowIndexY, highIndexX, highIndexY), outline=255)
		
		for p in splPoints:
			m,n = self.realToGrid(p)
			imPix[m,n] = 255
		
		filename = "checkMap%04u.png"
		im.save(filename % (self.saveCount-1))
		
		numPoints = 0
		
		for i in range(lowIndexX,highIndexX+1,1):
			for j in range(lowIndexY,highIndexY+1,1):
				xP, yP = self.gridToReal([i, j])
				
				dist = sqrt((x-xP)**2 + (y-yP)**2)
				if dist <= radius:
					# if there's any frontier points at all, return True
					if self.image[i,j] > 0:
						numPoints += 1
		
		
		print numPoints,"frontier points found: "
		if numPoints > 10:
			return True
		else:
			# no frontier points found in current region
			#print "frontier points NOT found: ", (self.saveCount-1)
			return False

	def update(self, boundMap, obstacleMap):
		
		self.occMap = boundMap.occMap
		self.boundMap = boundMap
		self.obstacleMap = obstacleMap
		
		self.xMin = boundMap.xMin
		self.xMax = boundMap.xMax
		self.yMin = boundMap.yMin
		self.yMax = boundMap.yMax

		self.points = boundMap.getBoundaryPoints()
		
		self.computeFrontiers()
		
		
	def computeFrontiers(self):
		
		
		# cycle through all boundary points in boundaryMap
		# for each:
		#       check if it is an obstacle on obstacleMap
		#       check if it is a shadow boundary or an unexplored boundary in occupancyMap
		#       non-obstacle + unexplored boundary has higher weight than only non-obstacle and shadow boundary
		
		# we don't consider frontier points next to the point of origin
		xBoundary, yBoundary = self.realToGrid((-8.0,0.0))
		
		self.resetMap()
		
		for i in range(self.xMin, self.xMax):
			if i > xBoundary:
				for j in range(self.yMin, self.yMax):
					
					if self.boundMap.image[i,j] == 255:
						if self.obstacleMap.image[i,j] == 0:
							
							# is it a shadow boundary or unexplored boundary?
							
							isShadow = False
							
							for m in range(-1,2):
								for n in range(-1,2):
									if self.occMap.image[i+m,j+n] == 0:
										isShadow = True
							
							if isShadow:
								self.image[i,j] = 127
							else:
								self.image[i,j] = 255
	
	
	def selectNextFrontier(self):
		# select for closeness and density of frontier points
		frontierDensity = Image.new('L', (self.numPixel,self.numPixel),0)
		fimg = frontierDensity.load()
		
		
		sumRange = 4
		num = (2*sumRange+1) * (2*sumRange+1)
		
		xMax = 0
		yMax = 0
		densityMax = 0
		
		for i in range(self.xMin, self.xMax):
			for j in range(self.yMin, self.yMax):
				
				fSum = 0
				for m in range(-sumRange, sumRange+1):
					for n in range(-sumRange, sumRange+1):
						fSum += self.image[i+m, j+n]
				
				density =  fSum/num
				fimg[i,j] = density
				
				if density > densityMax:
					xMax = i
					yMax = j
					densityMax = density

		draw = ImageDraw.Draw(frontierDensity)
		draw.ellipse((xMax-10, yMax-10, xMax+10, yMax+10), outline=255)
		
		filename = "densityMap%04u.png"
		frontierDensity.save(filename % (self.saveCount-1))
			
		x, y = self.gridToReal([xMax, yMax])
		
		return [x,y]
		
			