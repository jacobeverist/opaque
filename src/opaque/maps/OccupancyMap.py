import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from Map import Map
import Image


class OccupancyMap(Map):

	def __init__(self, probe, mapGraph, mapSize = MAPSIZE):
		Map.__init__(self, mapSize)

		self.probe = probe
		self.numJoints = self.probe.numSegs-1
		self.mapGraph = mapGraph
		
		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0

		self.fileName = "mapOccupancyMap%04u.png"

	def resetMap(self):
		self.mapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.image = self.mapImage.load()
		
		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0
	
	def update(self):
		
		self.resetMap()
		
		for i in range(self.mapGraph.numNodes):
			localNode = self.mapGraph.poseGraph.get_node_attributes(i)
			localOccMap = localNode.getOccMap()
			localMap = localOccMap.getMap()
			localImage = localMap.load()
			
			mapSize = localMap.size
			
			for j in range(mapSize[0]):
				for k in range(mapSize[1]):
					if localImage[j, k] == 255:

						pnt = localOccMap.gridToReal([j, k])
						
						pnt = localNode.convertLocalToGlobal(pnt)

						indexX, indexY = self.realToGrid(pnt)
						#if indexX >= 0 and indexX < self.numPixel and indexY >= 0 and indexY < self.numPixel:
						#	self.image[indexX,indexY] = 255
						
						" fill in the interstitial spaces to remove aliasing "
						for m in range(indexX - 1, indexX + 2):
							for n in range(indexY - 1, indexY + 2):
								if m >= 0 and m < self.numPixel and n >= 0 and n < self.numPixel:
									self.image[m, n] = 255

		pass
	

						