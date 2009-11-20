import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from Map import Map
from voronoi import *
import Image
import ImageDraw
import graph


class LocalVoronoiMap(Map):

	def __init__(self, localNode):
		
		mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0
		
		Map.__init__(self, mapSize = mapSize)
		
		self.localNode = localNode
		self.nodeID = localNode.getNodeID()
		
		self.fileName = "localVoronoiMap%03u" % self.nodeID + "_%04u.png"

		#self.update()
		
	def update(self):
		
		self.boundaryPoints = self.localNode.boundaryMap.getBoundaryPoints()
		self.occMap = self.localNode.occMap
		
		" standard python graph "
		self.roadGraph = graph.graph()
		
		" compute voronoi diagram "
		self.computeVoronoi()
		
		" prune the edges that are outside of free space and put into graph structure "
		self.pruneEdges()
		
	def getGraph(self):
		return self.roadGraph
	
	def computeVoronoi(self):
		
		'''
		Compute the voronoi diagram using a 3rd party library.
		
		Voronoi diagram data structures
		----
		
		self.vertices, 2-tuples
		self.lines, 3-tuples: a,b,c for a*x + b*y = c
		self.edges, 3-tuples, (l, v1, v2)  line, vertex 1, vertex 2
		vertex -1 means edge extends to infinity
		'''
		
		print "computing voronoi"
		print len(self.boundaryPoints), "boundary points"
		
		sites = []
		for p in self.boundaryPoints:
			sites.append(Site(p[0],p[1]))
		
		if sites != []:
			self.vertices, self.lines, self.edges = computeVoronoiDiagram(sites)
		else:
			self.vertices = []
			self.lines = []
			self.edges = []
		
		print len(self.vertices)
	def pruneEdges(self):
		
		'''remove edges and vertices that are outside the explored area'''

		occupancyImage = self.occMap.getMap()
		
		#NEIGHBOR_WIDTH = 3
		NEIGHBOR_WIDTH = 0

		xSize, ySize = self.size
		occPix = occupancyImage.load()

		newVertices = {}
		newEdges = []

		# checking all edges to check within the free space
		for edge in self.edges:
			v1 = edge[1]
			v2 = edge[2]

			if v1 != -1 and v2 != -1:

				# get grid value of real point
				i1, j1 = self.realToGrid(self.vertices[v1])
				i2, j2 = self.realToGrid(self.vertices[v2])
				
				if i1 < xSize and i1 >= 0 and i2 < xSize and i2 >= 0 and j1 < ySize and j1 >= 0 and j2 < ySize and j2 >= 0 :
					seg = [(i1,j1),(i2,j2)]

					indX1 = seg[0][0]
					indY1 = seg[0][1]
					indX2 = seg[1][0]
					indY2 = seg[1][1]

					isOK = True

					for i in range(-NEIGHBOR_WIDTH,NEIGHBOR_WIDTH+1,1):
						for j in range(-NEIGHBOR_WIDTH,NEIGHBOR_WIDTH+1,1):

							if indX1+i > self.numPixel or indX1+i < 0 or indY1+j > self.numPixel or indY1+j < 0:
								isOK = False
							else:
								val1 = occPix[indX1 + i, indY1 + j]
								if val1 <= 127: # 0 or 127
									isOK = False

							if indX2+i > self.numPixel or indX2+i < 0 or indY2+j > self.numPixel or indY2+j < 0:
								isOK = False
							else:
								val2 = occPix[indX2 + i, indY2 + j]
								if val2 <= 127: # 0 or 127
									isOK = False

					if isOK:
						# 1. add node if not already added
						# 2. add edge

						newVertices[v1] = self.vertices[v1]
						newVertices[v2] = self.vertices[v2]

						newEdges.append([v1,v2])

		for key, attr in newVertices.iteritems():
			self.roadGraph.add_node(key, attr)

		for edge in newEdges:
			self.roadGraph.add_edge(edge[0], edge[1])

			# assign weights to the edges according to their length
			v1 = self.roadGraph.get_node_attributes(edge[0])
			v2 = self.roadGraph.get_node_attributes(edge[1])
			length = sqrt((v2[0]-v1[0])**2+(v2[1]-v1[1])**2)
			self.roadGraph.set_edge_weight(edge[0],edge[1],length)
			
	def saveMap(self):
		
		draw = ImageDraw.Draw(self.mapImage)

		edges = self.roadGraph.edges()
		for edge in edges:
			v1 = self.roadGraph.get_node_attributes(edge[0])
			v2 = self.roadGraph.get_node_attributes(edge[1])
			
			# get grid value of real point
			i1, j1 = self.realToGrid(v1)
			i2, j2 = self.realToGrid(v2)

			draw.line([(i1,j1),(i2,j2)], fill = 255)

		Map.saveMap(self)
		
		
	def draw(self):
		edges = self.roadGraph.edges()
		for edge in edges:
			v1 = self.roadGraph.get_node_attributes(edge[0])
			v2 = self.roadGraph.get_node_attributes(edge[1])

			xP = [v1[0], v2[0]]
			zP = [v1[1], v2[1]]
			pylab.plot(xP,zP, color='0.5')

