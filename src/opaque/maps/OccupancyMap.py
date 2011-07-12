
from Map import Map
import Image
import ImageDraw
from functions import point_inside_polygon


class OccupancyMap(Map):

	def __init__(self, probe, mapGraph, mapSize):
		Map.__init__(self, mapSize)

		self.probe = probe
		self.numJoints = self.probe.numSegs-1
		self.mapGraph = mapGraph
		
		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0

		self.fileName = "mapOccupancyMap%04u.png"
		self.gndFileName = "mapGndOccupancyMap%04u.png"

	def resetMap(self):
		self.mapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.image = self.mapImage.load()

		self.gndMapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.gndImage = self.gndMapImage.load()
				
		self.xMin = self.numPixel
		self.xMax = 0
		self.yMin = self.numPixel
		self.yMax = 0

	def saveMap(self, point = []):

		localNode = self.mapGraph.nodeHash[self.mapGraph.numNodes-1]

		estPose = localNode.getGlobalGPACPose()
		indexX, indexY = self.realToGrid([estPose[0],estPose[1]])
			
		draw = ImageDraw.Draw(self.mapImage)
		#draw.ellipse((indexX-5, indexY-5, indexX+5, indexY+5), outline=0)	

		walls = self.mapGraph.getWalls()
		for wall in walls:
			wallPoints = []
			for p in wall:
				pGrid = self.realToGrid(p)
				wallPoints.append((pGrid[0],pGrid[1]))
			
			draw.line(wallPoints, fill=255)


		self.mapImage.save(self.fileName % self.saveCount)	

		gndPose = localNode.getGndGlobalGPACPose()
		indexX, indexY = self.realToGrid([gndPose[0],gndPose[1]])
		
		draw = ImageDraw.Draw(self.gndMapImage)
		#draw.ellipse((indexX-5, indexY-5, indexX+5, indexY+5), outline=0)	

		walls = self.mapGraph.getWalls()
		for wall in walls:
			wallPoints = []
			for p in wall:
				pGrid = self.realToGrid(p)
				wallPoints.append((pGrid[0],pGrid[1]))
			
			draw.line(wallPoints, fill=255)

		
		self.gndMapImage.save(self.gndFileName % self.saveCount)	
		self.saveCount += 1
			
	def loadMap(self, polygon):

		self.resetMap()

		for j in range(self.numPixel):
			for k in range(self.numPixel):
				pnt = self.gridToReal([j,k])

				if point_inside_polygon(pnt[0],pnt[1],polygon):
					self.image[j, k] = 255

	
	def update(self):
		
		self.resetMap()
		
		for i in range(self.mapGraph.numNodes):
			localNode = self.mapGraph.nodeHash[i]
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
												
						" fill in the interstitial spaces to remove aliasing "
						for m in range(indexX - 1, indexX + 2):
							for n in range(indexY - 1, indexY + 2):
								if m >= 0 and m < self.numPixel and n >= 0 and n < self.numPixel:
									self.image[m, n] = 255
	
		for i in range(self.mapGraph.numNodes):
			localNode = self.mapGraph.nodeHash[i]
			localOccMap = localNode.getOccMap()
			localMap = localOccMap.getGndMap()
			localImage = localMap.load()
			
			mapSize = localMap.size
			
			for j in range(mapSize[0]):
				for k in range(mapSize[1]):
					if localImage[j, k] == 255:

						pnt = localOccMap.gridToReal([j, k])
						
						pnt = localNode.convertLocalToGndGlobal(pnt)

						indexX, indexY = self.realToGrid(pnt)
												
						" fill in the interstitial spaces to remove aliasing "
						for m in range(indexX - 1, indexX + 2):
							for n in range(indexY - 1, indexY + 2):
								if m >= 0 and m < self.numPixel and n >= 0 and n < self.numPixel:
									self.gndImage[m, n] = 255
	


						