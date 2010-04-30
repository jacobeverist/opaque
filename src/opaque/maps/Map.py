
import Image
from math import *

# Map Space Parameters
PIXELSIZE = 0.05


class Map:

	def __init__(self, mapSize):
		self.pixelSize = PIXELSIZE
		self.mapSize = mapSize
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		self.fileName = "generalMap%04u.png"
		self.saveCount = 0
		
		self.resetMap()
		
	def setMapSize(self, mapSize):
		self.pixelSize = PIXELSIZE
		self.mapSize = mapSize
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
	def resetMap(self):
		self.mapImage = Image.new('L', (self.numPixel,self.numPixel),0)
		self.image = self.mapImage.load()
		
	def loadImage(self, img):
		self.mapImage = img
		self.image = self.mapImage.load()

		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0
		
		for i in range(self.numPixel):
			for j in range(self.numPixel):
				val = self.image[i,j]
				if val == 0 or val == 255:
					if i < self.xMin:
						self.xMin = i
					if i > self.xMax:
						self.xMax = i
					if j < self.yMin:
						self.yMin = j
					if j > self.yMax:
						self.yMax = j
						
	def __getattr__(self, name):
		if name == 'size':
			return self.mapImage.size
		else:
			raise AttributeError

	def getImage(self):
		return self.mapImage

	def copyImage(self):
		return self.mapImage.copy()

	def realToGrid(self, point):
		indexX = int(floor(point[0]*self.divPix)) + self.numPixel/2 + 1
		indexY = int(floor(point[1]*self.divPix)) + self.numPixel/2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0, (j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]
		return point

	def saveMap(self, point = []):

		self.mapImage.save(self.fileName % self.saveCount)	
		self.saveCount += 1
