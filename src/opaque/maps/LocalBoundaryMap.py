import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from Map import Map
import Image
import ImageDraw
import graph
import traceback

class LocalBoundaryMap(Map):

	#def __init__(self, occMap):
	def __init__(self, localNode):
		Map.__init__(self)
		
		self.localNode = localNode
		self.occMap = self.localNode.getOccMap()
		
		self.nodeCount = 0
		self.nodeID = self.localNode.getNodeID()
		self.fileName = "localBoundaryMap%03u" % self.nodeID + "_%04u.png"
		self.boundGraph = graph.graph()

		self.pixelSize = self.occMap.pixelSize
		self.mapSize = self.occMap.mapSize
		self.numPixel = self.occMap.numPixel
		self.halfPix = self.occMap.halfPix
		self.divPix = self.occMap.divPix

		self.resetMap()

		self.xMin = self.numPixel
		self.xMax = -1
		self.yMin = self.numPixel
		self.yMax = -1
		
		#self.update(occMap)

	def saveMap(self, poly = []):

		if len(poly) > 0:
			tempImage = self.mapImage.copy()
			draw = ImageDraw.Draw(tempImage)
			
			pixPoly = []
			for p in poly:
				x, y = self.realToGrid(p)
				pixPoly.append([x,y])
				
			
			num = len(poly)
			
			for i in range(num):
				draw.line((pixPoly[i][0], pixPoly[i][1], pixPoly[(i+1) % num][0], pixPoly[(i+1) % num][1]), fill = 255)
			
			tempImage.save(self.fileName % self.saveCount)	
		else:
			self.mapImage.save(self.fileName % self.saveCount)	
		
		self.saveCount += 1
			
	def update(self):

		self.occMap = self.localNode.getOccMap()
		
		#print "computing boundary"
		
		#print self.occMap.changed
		
		#print traceback.print_stack()

		
		if self.occMap.changed:
		
			self.resetMap()

			self.computeBoundary(self.occMap)
			
			self.getBoundaryPoints()
			
			self.occMap.changed = False
	
	def getBoundaryPoints(self):
	# return the boundary as a set of (x,y) points
		
		points = []
		
		for i in range(self.xMin-1,self.xMax+2,1):
			for j in range(self.yMin-1,self.yMax+2,1):
				#print self.image[i,j]
				if self.image[i,j] == 255:
					#print i,j
					p = self.gridToReal([i,j])
					points.append(p)
					
		return points
	
	def getBoundaryGraph(self):
		
		self.rootNode = ""
				
		for i in range(self.xMin-1,self.xMax+2,1):
			for j in range(self.yMin-1,self.yMax+2,1):
				if self.image[i,j] == 255:
					if self.rootNode == "":
						self.rootNode = repr((i,j))
						
					self.boundGraph.add_node(repr((i,j)))

		for i in range(self.xMin-1,self.xMax+2,1):
			for j in range(self.yMin-1,self.yMax+2,1):
				if self.image[i,j] == 255:
					neighbors = self.boundGraph.neighbors(repr((i,j)))
					
					for k in range(-1,2):
						for l in range(-1,2):
							if not ( k == 0 and l == 0):
								if self.image[i+k,j+l] == 255:
									if neighbors.count(repr((i+k,j+l))) == 0:
										self.boundGraph.add_edge(repr((i,j)),repr((i+k,j+l)))
										#self.boundGraph.add_edge((i+k,j+l),(i,j))

		return self.boundGraph

	def computeBoundary(self, occMap):
		
		occImage = occMap.getMap()
		occAccess = occImage.load()
		
		for i in range(1,self.numPixel-1):
			for j in range(1,self.numPixel-1):

				isBoundary = False

				# in shadow or unexplored, so it could be a boundary
				if occAccess[i,j] <= 127:
					#print i,j

					# check if any neighbors are free space
					for m in range(-1,2):
						for n in range(-1,2):
							if occAccess[i+m,j+n] == 255:
								isBoundary = True
				elif occAccess[i,j] == 255:
					#print i,j
					self.image[i,j] = 127

				# if a neighbor is free space, then this is a boundary point
				if isBoundary:
					#print i,j
					self.image[i,j] = 255
					
					if i > self.xMax:
						self.xMax = i
					if i < self.xMin:
						self.xMin = i
					if j > self.yMax:
						self.yMax = j
					if j < self.yMin:
						self.yMin = j
					
					