
from Map import Map
import Image
import pylab
import ImageDraw
from SplineFit import SplineFit
from math import floor, sqrt
from numpy import arange
from copy import copy

class FrontierMap(Map):

	#def __init__(self, boundMap, obstacleMap, saveCount):
	def __init__(self, mapGraph, mapSize):
		Map.__init__(self, mapSize)
				
		self.fileName = "mapFrontierGraph%04u.png"

		self.mapGraph = mapGraph

		self.checkSpline = 0


		#self.inhibitRadius = 1.0
		#self.inhibits = [[-3.0,0.0], [-2.0,0.0], [-1.0,0.0]]

		self.inhibitRadius = 0.5
		self.inhibits = []
						
		self.densityThreshold = 30

		#self.update()
	
	def saveMap(self):

		tempImage = self.mapImage.copy()
		
		tempImage.save(self.fileName % self.saveCount)	
		self.saveCount += 1
		
	def isFrontier(self):
		" checks if there are frontier points directly in front of snake "
		" if there are enough, it will return true "
		
		
		points = self.mapGraph.getCenterPoints()
		
		print "centerPoints =", points
			
		spline = SplineFit(points, kp=2)
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
		"""
		if lowIndexX < self.xMin:
			self.xMin = lowIndexX

		if highIndexX > self.xMax:
			self.xMax = highIndexX

		if lowIndexY < self.yMin:
			self.yMin = lowIndexY

		if highIndexY > self.yMax:
			self.yMax = highIndexY
		"""
			
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
		
		filename = "mapCheckMap%04u.png"
		im.save(filename % (self.saveCount))
		
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
			self.inhibitLocation([x,y])
			return True
		else:
			# no frontier points found in current region
			#print "frontier points NOT found: ", (self.saveCount-1)
			return False

	def update(self):
		
		#self.xMin = self.boundMap.xMin
		#self.xMax = self.boundMap.xMax
		#self.yMin = self.boundMap.yMin
		#self.yMax = self.boundMap.yMax

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
		
		self.obstacleMap = self.mapGraph.obstMapImage
		
		self.resetMap()
		
		sumRange = 3
		num = (2*sumRange+1) * (2*sumRange+1)
		
		densityMax = 0
		densityMin = 1e100

		boundImg = self.mapGraph.boundMap.image
		obstImg = self.obstacleMap.load()

		for i in range(self.numPixel/2+20,self.numPixel-sumRange):
			for j in range(sumRange,self.numPixel-sumRange):
				if boundImg[i,j] == 255:
	
					isTooClose = False
					xReal , yReal = self.gridToReal([i,j])
					
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
								fSum += obstImg[i+m, j+n]
						
						density =  fSum/num
						density = 255 - density
						self.image[i,j] = density
						
						if density > densityMax:
							densityMax = density
							
						if density < densityMin:
							densityMin = density

		maxDen = densityMax - densityMin

		#for i in range(sumRange,self.numPixel-sumRange):
		for i in range(self.numPixel/2+20,self.numPixel-sumRange):
			for j in range(sumRange,self.numPixel-sumRange):
				if boundImg[i,j] == 255:

					isTooClose = False
					xReal , yReal = self.gridToReal([i,j])

					" make sure we only consider points that are outside inhibition circles "
					for k in range(len(self.inhibits)):
						iP = self.inhibits[k]
						dist = sqrt((iP[0]-xReal)**2 + (iP[1]-yReal)**2)
						if dist < self.inhibitRadius:
							isTooClose = True
							break

					if not isTooClose:
						density =  self.image[i,j]
						density = density - densityMin
						if maxDen > 0:
							result = 255*(1.0-float(maxDen-density)/float(maxDen))
							self.image[i,j] = result
						#print maxDen, density, result, self.image[i,j]

		#self.mapImage.save("mapFrontierGraph%04u.png" % self.saveCount)		
	
	
	def inhibitLocation(self, point):
		self.inhibits.append(copy(point))
	
	def selectNextFrontier(self, direction = True):
		
		# select for closeness and density of frontier points
		frontierDensity = Image.new('L', (self.numPixel,self.numPixel),0)
		fimg = frontierDensity.load()

		densityMax = 0
		xMax = 0
		yMax = 0
		frontWidth = 4
		normFactor = (2*frontWidth+1) * (2*frontWidth+1)

		#boundImg = self.boundMap.load()

		#for i in range(frontWidth,self.numPixel-frontWidth):
		for i in range(self.numPixel/2+20,self.numPixel-frontWidth):
			for j in range(frontWidth,self.numPixel-frontWidth):
				#if self.image[i,j] > 0:
				if True:
					sum = 0
					for k in range(i-frontWidth,i+frontWidth):
						for l in range(j-frontWidth,j+frontWidth):						
							sum += self.image[k,l]

					density =  sum/normFactor
					fimg[i,j] = density

					if density > densityMax:
						densityMax = density
						xMax = k
						yMax = l
					
		draw = ImageDraw.Draw(frontierDensity)
		draw.ellipse((xMax-5, yMax-5, xMax+5, yMax+5), outline=255)

		for k in range(len(self.inhibits)):
			cP = self.inhibits[k]
			radius = self.inhibitRadius
			radius = int(radius*self.divPix)
			xGrid , yGrid = self.realToGrid(cP)
			draw.ellipse((xGrid-radius, yGrid-radius, xGrid+radius, yGrid+radius), outline=255)
		
		frontierDensity.save("mapFrontierDensity%04u.png" % self.saveCount)		

		#self.saveCount += 1
			
		x, y = self.gridToReal([xMax, yMax])

		print "maximum density point"
		print x, y, densityMax

		if densityMax < self.densityThreshold:
			raise
		
		self.inhibitLocation([x,y])

		return [x,y]
		
			