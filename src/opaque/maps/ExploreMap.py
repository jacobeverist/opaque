import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Map import Map
from copy import *
from math import *

class ExploreMap(Map):

	def __init__(self, boundaryPoints):
		Map.__init__(self)
		
		self.boundaryPoints = boundaryPoints
		self.fileName = "exploreMap%04u.png"

	# select nearest frontier point to the robot's position
	def selectExplorationPoint(self, point):

		# convert point to grid coordinates
		indexX, indexY = self.realToGrid(point)

		if indexX >= self.numPixel or indexX < 0:
			print "point X out of range of map", point[0], indexX
			raise

		if indexY >= self.numPixel or indexY < 0:
			print "point Y out of range of map", point[1], indexY
			raise

		searchComplete = False
		targetFound = False
		targetIndex = (0,0)

		dist = 1e100 
		minIndex = [0,0]
		for i in range(self.numPixel):
			for j in range(self.numPixel):
				if self.image[i,j] == 255:
					currDist = sqrt((indexX-i)**2 + (indexY-j)**2)
					if currDist < dist:
						dist = currDist
						minIndex = [i,j]

		if dist == 1e100:
			print "Error: No frontier point found"
			raise

		# convert indices to real coordinates
		targetPoint = self.gridToReal(minIndex)
		return targetPoint

	def markExplored(self, point):

		indX, indY = self.realToGrid(point)
		
		for i in range(indX-1, indX+2):
			for j in range(indY-1, indY+2):
				self.image[i,j] = 127
