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


class LocalFrontierMap(Map):

	def __init__(self, localNode):

		mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0		
		Map.__init__(self, mapSize = mapSize)

		self.localNode = localNode
		self.probe = self.localNode.probe
		self.nodeID = self.localNode.getNodeID()
		self.fileName = "localFrontierMap%03u" % self.nodeID + "_%04u.png"
		
		self.checkSpline = 0
		
		self.inhibits = []
		self.inhibitRadius = 0.25
		self.densityThreshold = 30
	
	def saveMap(self):

		tempImage = self.mapImage.copy()
		draw = ImageDraw.Draw(tempImage)
		
		tempImage.save(self.fileName % self.saveCount)	
		#self.saveCount += 1
		
	def isFrontier(self):
		" checks if there are frontier points directly in front of snake "
		" if there are enough, it will return true "
		
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

	def update(self):
		
		self.occMap = self.localNode.getOccMap()
		self.boundMap = self.localNode.getBoundMap()
		self.obstacleMap = self.localNode.getObstacleMap()
		
		self.xMin = self.boundMap.xMin
		self.xMax = self.boundMap.xMax
		self.yMin = self.boundMap.yMin
		self.yMax = self.boundMap.yMax

		self.computeFrontiers()
		
		
	def computeFrontiers(self):
		
		"""
		cycle through all boundary points in boundaryMap
		  for each:
		       check if it is an obstacle on obstacleMap
		       otherwise, it is a frontier
		"""
		
		# we don't consider frontier points next to the point of origin
		# xBoundary, yBoundary = self.realToGrid((-8.0,0.0))
		
		self.resetMap()
		
		for i in range(self.xMin, self.xMax):
			for j in range(self.yMin, self.yMax):

				if self.boundMap.image[i,j] == 255:
					if self.obstacleMap.image[i,j] == 0:
										
						self.image[i,j] = 255
	
	def inhibitLocation(self, point):
		self.inhibits.append(copy(point))
	
	def selectNextFrontier(self, direction = True):
		
		# select for closeness and density of frontier points
		frontierDensity = Image.new('L', (self.numPixel,self.numPixel),0)
		fimg = frontierDensity.load()
		
		
		sumRange = 4
		num = (2*sumRange+1) * (2*sumRange+1)
		
		xMax = 0
		yMax = 0
		densityMax = 0
		
		
		cPoints = deepcopy(self.localNode.getCenterPoints())
		radii = [0.7 for i in range(len(cPoints))]

		print cPoints
		print radii
	
		if direction:
			cPoints = cPoints[1:]
			radii = radii[1:]
			radii[-1] = 1.0
		else:
			cPoints = cPoints[:-1]
			radii = radii[:-1]
			radii[0] = 1.0
		
		for i in range(self.xMin, self.xMax):
			for j in range(self.yMin, self.yMax):
				
				isTooClose = False

				xReal , yReal = self.gridToReal([i,j])
				
				" make sure we only consider points outside of the immediate proximity "
				for k in range(len(cPoints)):
					cP = cPoints[k]
					dist = sqrt((cP[0]-xReal)**2 + (cP[1]-yReal)**2)
					if dist < radii[k]:
						isTooClose = True
						break

				" make sure we only consider points that are outside inhibition circles "
				for k in range(len(self.inhibits)):
					iP = self.inhibits[k]
					dist = sqrt((iP[0]-xReal)**2 + (iP[1]-yReal)**2)
					if dist < self.inhibitRadius:
						isTooClose = True
						break
					
				if not isTooClose:
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
		
		print "maximum density =", densityMax

		draw = ImageDraw.Draw(frontierDensity)
		draw.ellipse((xMax-10, yMax-10, xMax+10, yMax+10), outline=255)
		
		for k in range(len(cPoints)):
			cP = cPoints[k]
			radius = radii[k]
			radius = int(radius/self.pixelSize)
			xGrid , yGrid = self.realToGrid(cP)
			draw.ellipse((xGrid-radius, yGrid-radius, xGrid+radius, yGrid+radius), outline=255)

		fileName = "localDensityMap%03u" % self.nodeID + "_%04u.png"		
		frontierDensity.save(fileName % (self.saveCount))
		self.saveCount += 1
			
		x, y = self.gridToReal([xMax, yMax])

		if densityMax < self.densityThreshold:
			raise
		
		return [x,y]
		
		