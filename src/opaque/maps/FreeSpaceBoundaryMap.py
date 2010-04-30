
from Map import Map
import Image
import graph

class FreeSpaceBoundaryMap(Map):

	def __init__(self, occMap, mapSize):
		Map.__init__(self, mapSize)
		
		self.nodeCount = 0
		self.fileName = "mapBoundaryMap%04u.png"
		self.boundGraph = graph.graph()
		self.update(occMap)
	
	def update(self, occMap):
		
		self.occMap = occMap
		#self.xMin = occMap.xMin
		#self.xMax = occMap.xMax
		#self.yMin = occMap.yMin
		#self.yMax = occMap.yMax
		
		self.xMin = 2
		self.xMax = self.numPixel - 3
		self.yMin = 2
		self.yMax = self.numPixel - 3
		
		self.resetMap()

		self.computeBoundary(occMap)
	
	def getBoundaryPoints(self):
	# return the boundary as a set of (x,y) points
		
		points = []
		
		for i in range(self.xMin-1,self.xMax+2,1):
			for j in range(self.yMin-1,self.yMax+2,1):
				if self.image[i,j] == 255:
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
		
		occImage = occMap.getImage()
		occAccess = occImage.load()
		
		for i in range(self.xMin-1,self.xMax+2,1):
			for j in range(self.yMin-1,self.yMax+2,1):

				isBoundary = False

				# in shadow or unexplored, so it could be a boundary
				if occAccess[i,j] <= 127:
					#print i,j

					# check if any neighbors are free space
					for m in range(-1,2):
						for n in range(-1,2):
							#print self.numPixel, i, m, j, n
							if occAccess[i+m,j+n] == 255:
								isBoundary = True
				elif occAccess[i,j] == 255:
					self.image[i,j] = 127

				# if a neighbor is free space, then this is a boundary point
				if isBoundary:
					self.image[i,j] = 255


					