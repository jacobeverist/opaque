import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from numpy import *
from scipy.optimize import *
import graph
import csv
from LocalNode import *
from VoronoiMap import *
from FrontierMap import *
from NavRoadMap import *
from OccupancyMap import *
from FreeSpaceBoundaryMap import *
from Pose import Pose
import pylab
from matplotlib.patches import Circle
from pose import *

from motion import gen_icp


class MapGraph:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)
		
		self.initPose = self.probe.getActualJointPose(19)
		self.poseGraph = graph.graph()
		self.numNodes = 0
		self.currNode = 0

		self.pixelSize = PIXELSIZE
		self.mapSize = 20.0
		self.numPixel = int(2.0 * self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize / 2.0
		self.divPix = floor((2.0 * self.mapSize / self.pixelSize) / self.mapSize)
		self.fileName = "mapGraph%04u.png"
		self.saveCount = 0


		" ground truth walls of the environment "
		self.gndMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.gndImage = self.gndMapImage.load()
		gndDraw = ImageDraw.Draw(self.gndMapImage)
		
		walls = self.probe.getWalls()
		for wall in walls:
			wallPoints = []
			for p in wall:
				pGrid = self.realToGrid(p)
				wallPoints.append((pGrid[0],pGrid[1]))
			
			gndDraw.line(wallPoints, fill=255)

		self.gndMapImage.save("mapGndGraph.png")

		self.occMap = OccupancyMap(self.probe, self, self.mapSize)
		self.boundMap = FreeSpaceBoundaryMap(self.occMap, self.mapSize)
		self.frontierMap = FrontierMap(self, self.mapSize)
		self.voronoiMap = VoronoiMap(self, self.mapSize)

		self.childNodes = []
		self.childEntities = []
		self.boundParentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("globalBoundNode")

		#self.newNode()

	def loadFile(self, num_poses):
		
		#self.saveCount = 0
		self.poseGraph = graph.graph()
		self.numNodes = 0
		self.currNode = 0
		
		#for i in range(0,22):
		#for i in range(0, 14):
		for i in range(0,num_poses):
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, inSim = False)
			self.currNode.readFromFile(i)
			
			self.poseGraph.add_node(i, self.currNode)

			if self.numNodes > 0:
				self.poseGraph.add_edge(i-1, i)

			self.numNodes += 1
	
		#self.saveMap()
		#self.correctPoses()
		#self.saveMap()
		
	def setCenterPoints(self, centerPoints):
		self.currNode.setCenterPoints(centerPoints)
		
	def getCenterPoints(self):
		return self.currNode.getCenterPoints()
	
	def getCurrentNode(self):
		return self.currNode
	
	def newNode(self):
		
		#if self.numNodes >= 2:
		#	self.correctPoses2()
		
		if self.numNodes > 0:
			if self.currNode.isDirty():
				self.currNode.synch()
				#self.synch()
			#self.saveMap()
			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19)
		
		self.poseGraph.add_node(self.numNodes, self.currNode)
		
		if self.numNodes > 0:
			self.poseGraph.add_edge(self.numNodes - 1, self.numNodes)
		
		self.numNodes += 1
		
		#self.contacts.resetPose()
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
	
	def keepStablePose(self):
		self.stablePose.update()
	
	def correctPoses2(self):

		if self.currNode.isDirty():
			self.currNode.synch()
			#self.synch()

		" TUNE ME:  threshold cost difference between iterations to determine if converged "
		costThresh = 0.1
	
		" TUNE ME:   minimum match distance before point is discarded from consideration "
		minMatchDist = 2.0
	
		" plot the best fit at each iteration of the algorithm? "
		plotIteration = True
		#plotIteration = False
	
		" initial guess for x, y, theta parameters "
		offset = [0.0,0.0,0.0]
	
		" Extract the data from the files and put them into arrays "
		estPoses = []
		gndPoses = []
		poseNumbers = []
		for i in range(0,self.numNodes):
		
			node1 = self.poseGraph.get_node_attributes(i)
				
			estPose1 = node1.getEstPose()	
			gndPose1 = node1.getGndPose()
	
			estPoses.append(estPose1)
			gndPoses.append(gndPose1)
			poseNumbers.append(i)


		def decimatePoints(points):
			result = []
		
			for i in range(len(points)):
				if i%2 == 0:
					result.append(points[i])
		
			return result
				
		" Read in data of Alpha-Shapes and add their associated covariances "
		a_hulls = []
		for i in range(len(estPoses)):

			node1 = self.poseGraph.get_node_attributes(i)
			node1.computeAlphaBoundary()			
			a_data = node1.getAlphaBoundary()
			a_data = decimatePoints(a_data)
	
			" treat the points with the point-to-line constraint "
			gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
	
			" treat the points with the distance-from-origin increasing error constraint "
			gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
	
			a_hulls.append(a_data)
	
		occMaps = []
		for m in range(0,len(estPoses)):
			offset = estPoses[m]
			occMap = self.getNodeOccMap(m)
	
			mapImage = occMap.getMap()
			image = mapImage.load()
			
			" 1. pick out the points "
			points = []
			for j in range(occMap.numPixel):
				for k in range(occMap.numPixel):
					if image[j,k] == 255:
						pnt = occMap.gridToReal([j,k])
						points.append(gen_icp.dispOffset(pnt, offset))
			
			occMaps.append(points)
			
	
		a_hull_trans = []
		for m in range(0,len(estPoses)):
			offset = estPoses[m]
	
			" transform the past poses "
			past_data = a_hulls[m]
			a_trans = []
			for p in past_data:
				a_trans.append(gen_icp.dispPoint(p, offset))
			
			a_hull_trans.append(a_trans)
			
		offsets = []
		for i in range(len(estPoses)-1):
			estPose1 = estPoses[i]
			estPose2 = estPoses[i+1]
			
			pose1 = Pose(estPose1)
			offset = pose1.convertGlobalPoseToLocal(estPose2)
			offsets.append(offset)
				
		k = len(estPoses)-1
	
		" 1. target pose to correct "
		targetPose = estPoses[-1]
		targetHull = a_hulls[-1]
	
		a_hull_trans = []
		estPoseOrigin = estPoses[-2]
		for m in range(0,len(estPoses)-1):
			estPose2 = estPoses[m]
			poseOrigin = Pose(estPoseOrigin)
			offset = poseOrigin.convertGlobalPoseToLocal(estPose2)

			" transform the past poses "
			past_data = a_hulls[m]
			a_trans = []
			for p in past_data:
				a_trans.append(gen_icp.dispPoint(p, offset))
			
			a_hull_trans.append(a_trans)

		pastCircles = []
		for m in range(0,len(estPoses)-1):
			hull = a_hull_trans[m]
			radius, center = gen_icp.computeEnclosingCircle(hull)
			pastCircles.append([radius,center])
			
		pastPose = estPoseOrigin
		
		pastHull = gen_icp.computeUnions(a_hull_trans)
		gen_icp.addPointToLineCovariance(pastHull, high_var=1.0, low_var=0.001)
		#gen_icp.addDistanceFromOriginCovariance(pastHull, tan_var=0.1, perp_var=0.01)

		" run generalized ICP (a plot is made for each iteration of the algorithm) "
		offset = gen_icp.gen_ICP_global(pastPose, targetPose, pastHull, targetHull, pastCircles, costThresh, minMatchDist, plotIteration)
		
		offsets[k-1] = offset
		
		" recompute the estimated poses with the new offset "
		newEstPoses = []
		newEstPoses.append(estPoses[0])
		
		for m in range(len(offsets)):
			pose1 = Pose(newEstPoses[m])
			offset = offsets[m]
			newEstPose2 = pose1.convertLocalOffsetToGlobal(offset)
			newEstPoses.append(newEstPose2)
			
			
		estPoses = newEstPoses
					
		" update the estimated poses "
		for m in range(0,len(estPoses)):
			self.setNodePose(m, estPoses[m])
			
			
		" update the current estimated pose in AverageContacts "
		self.contacts.resetPose(estPose = estPoses[-1])

		#self.saveMap()
		
	def correctPoses(self):
			
		# 2. compute the global position of all the poses wrt to first pose
	
		# TUNE ME:  threshold cost difference between iterations to determine if converged
		#costThresh = 0.004
		costThresh = 0.1
	
		# TUNE ME:   minimum match distance before point is discarded from consideration
		minMatchDist = 2.0
	
		# plot the best fit at each iteration of the algorithm?
		plotIteration = True
		#plotIteration = False
	
		# initial guess for x, y, theta parameters
		offset = [0.0,0.0,0.0]
	
		# sample data
		a_data = []
		b_data = []
	
		estPoses = []
		finalOffsets = []

		for i in range(1, self.numNodes-1):
			
			node1 = self.poseGraph.get_node_attributes(i)
			node2 = self.poseGraph.get_node_attributes(i+1)
			
			print "correcting", i , "and", i+1

			node1.computeAlphaBoundary()
			node2.computeAlphaBoundary()
			
			b_data = node1.getAlphaBoundary()
			a_data = node2.getAlphaBoundary()
			
			def decimatePoints(points):
				result = []
			
				for i in range(len(points)):
					if i%2 == 0:
						result.append(points[i])
			
				return result
			
			a_data = decimatePoints(a_data)
			b_data = decimatePoints(b_data)

			estPose1 = node1.getEstPose()
			estPose2 = node2.getEstPose()
			
			# get the offset of node 2 with respect to node 1
			pose1 = Pose(estPose1)
			offset = pose1.convertGlobalPoseToLocal(estPose2)

			# treat the points with the point-to-line constraint
			gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
			gen_icp.addPointToLineCovariance(b_data, high_var=1.0, low_var=0.001)

			gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
			gen_icp.addDistanceFromOriginCovariance(b_data, tan_var=0.1, perp_var=0.01)

			# plot the data without A transformed, plot 997
			gen_icp.draw(a_data, b_data, "rawData.png")

			# transform the points in A by 'offset'
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispOffset(p, offset))

			# plot the data with A transformed, plot 998
			gen_icp.draw(a_trans, b_data, "initialGuess.png") 

			# run generalized ICP (a plot is made for each iteration of the algorithm)
			offset = gen_icp.gen_ICP(offset, a_data, b_data, costThresh, minMatchDist, plotIteration)

			# transform the points of A with parameters determined by algorithm
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispPoint(p, offset))

			print "drawing: ", "finalOutput_%04u.png" % i

			# plot the final result, plot 999
			gen_icp.draw(a_trans, b_data, "finalOutput_%04u.png" % i, fileWrite = True) 

			finalOffsets.append(offset)
			estPoses.append(estPose1)

		"""
		1. get the relative offsets (done)
		2. for each node, recompute the estimated pose of all the subsequent poses
		3. plot them
		"""
		
		for i in range(len(finalOffsets)):
			
			est1 = estPoses[i]
	
			offset = finalOffsets[i]
	
			pose1 = Pose(est1)
	
			newEst2 = pose1.convertLocalOffsetToGlobal(offset)
	
			if i+1 < len(estPoses):
				estPoses[i+1] = newEst2
			else:
				estPoses.append(newEst2)

		pylab.clf()

		for i in range(len(estPoses)):
			est1 = estPoses[i]
			node1 = self.poseGraph.get_node_attributes(i+1)
			node1.setEstPose(est1)

	
		for i in range(1, self.numNodes-1):

			#est1 = estPoses[i-1]
			node1 = self.poseGraph.get_node_attributes(i)
			est1 = node1.getEstPose()
			
			a_data = node1.estOccPoints
			b_data = []
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispOffset(p, est1))
				
			gen_icp.draw( a_trans, b_data,"foo.png", fileWrite = False) 
	
	
			#gnd1 = gndPoses[i]
			#a_data = eval(val)
			#b_data = []
			#a_trans = []
			#for p in a_data:
			#	a_trans.append(gen_icp.dispOffset(p, gnd1))
			#gen_icp.draw(b_data, a_trans, "foo.png", fileWrite = False) 
			
		def plotEnv():
	
			WLEN = 3.0
			WLEN2 = 5.0
			wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
			wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			w1 = wall1[2]
			w2 = wall2[0]
			
			wall3 = [[w1[0] + 0.4*math.cos(math.pi/6), w1[1] + 0.4*math.sin(math.pi/6)], [0.4*math.cos(math.pi/6) - 4, 0.0], [w2[0] + 0.4*math.cos(math.pi/6), w2[1] - 0.4*math.sin(math.pi/6)], w2]
			lp = wall3[0]
			rp = wall3[2]
			
			wall6 = [lp, [lp[0] + WLEN*math.cos(math.pi/6), lp[1] + WLEN*math.sin(math.pi/6)]]
			wall6.append([wall6[1][0] + 0.4*math.cos(math.pi/3), wall6[1][1] - 0.4*math.sin(math.pi/3)])
			wall6.append([wall6[2][0] - WLEN*math.cos(math.pi/6), wall6[2][1] - WLEN*math.sin(math.pi/6)])
			wall6.append([wall6[3][0] + WLEN*math.cos(math.pi/3), wall6[3][1] - WLEN*math.sin(math.pi/3)])
			wall6.append([wall6[4][0] - 0.4*math.cos(math.pi/6), wall6[4][1] - 0.4*math.sin(math.pi/6)])
			wall6.append(w1)
			wall6.reverse()
		
			walls = [wall1, wall2, wall3, wall6]
		
			for wall in walls:
				xP = []
				yP = []
				for i in range(len(wall)):
					p = copy(wall[i])
					p[0] += 6.0
					wall[i] = p
					xP.append(p[0])
					yP.append(p[1])
		
				pylab.plot(xP,yP, linewidth=2, color = 'g')
		
		plotEnv()
	
		pylab.xlim(-3,6)
		pylab.ylim(-4,4)
	
		gen_icp.save("finalMap.png")	
	
	def drawEstBoundary(self):
		
		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []

		# remove all children
		self.boundParentNode.removeAllChildren()

		points = self.boundMap.getBoundaryPoints()

		if self.currNode != 0:
			onSegPose = Pose(self.contacts.getClosestPose(self.currNode.rootNode))
			onSegActPose = Pose( self.probe.getActualJointPose(self.currNode.rootNode))

		for i in range(len(points)):
			pnt = points[i]
			
			if self.currNode != 0:
				newPnt = copy(pnt)				
				localPnt = onSegPose.convertGlobalToLocal(pnt)
				pnt = onSegActPose.convertLocalToGlobal(localPnt)

			childNode = self.boundParentNode.createChildSceneNode("globalBoundPoint" + "_" + str(i))
			self.childNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("globalBoundEnt" + "_" + str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Black")
			self.childEntities.append(currEntity)
	
			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			childNode.setPosition(position)
	
			size = ogre.Vector3(0.03,0.03,0.03)
			childNode.setScale(size)
	
			childNode.attachObject(currEntity)	
	
	def synch(self):
		self.currNode.synch()

		self.occMap.update()
		self.boundMap.update(self.occMap)
		
		" 1. for each node, copy local image "
		" 2. for each point in image, transform appropriately and plot to main image "
		self.obstMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.obstImage = self.obstMapImage.load()

		for i in range(self.numNodes):
			localNode = self.poseGraph.get_node_attributes(i)
			#localNode.saveMap()
			localObstMap = localNode.getObstacleMap()
			localMap = localObstMap.getMap()
			localImage = localMap.load()
			
			mapSize = localMap.size
		
			for j in range(mapSize[0]):
				for k in range(mapSize[1]):
					if localImage[j, k] == 255:

						pnt = localObstMap.gridToReal([j, k])

						pnt = localNode.convertLocalToGlobal(pnt)

						indexX, indexY = self.realToGrid(pnt)
						#if indexX >= 0 and indexX < self.numPixel and indexY >= 0 and indexY < self.numPixel:
						#	self.image[indexX,indexY] = 255
						
						" fill in the interstitial spaces to remove aliasing "
						#for m in range(indexX-1,indexX+2):
						#	for n in range(indexY-1,indexY+2):
						#		if m >= 0 and m < self.numPixel and n >= 0 and n < self.numPixel:
						#			self.obstImage[m,n] = 255

						self.obstImage[indexX, indexY] = 255
							

		
		self.frontierMap.update()
		
		self.voronoiMap.update()

	def forceUpdate(self, isForward=True):

		self.stablePose.setDirection(isForward)
		self.currNode.update(isForward)

	def update(self, isForward=True):

		self.stablePose.setDirection(isForward)

		#if True:
		if self.stablePose.isStable():
			#print "stable"
			self.currNode.update(isForward)
		else:
			#print "not stable"
			pass
		
		#self.currNode.update(isForward)
		
		#self.saveMap()

	def computeHeadPath(self, currPose, frontierPoint, exploreRoot):
		vals = self.navRoadMap.computeHeadPath(currPose, frontierPoint, exploreRoot)

		#self.navRoadMap.draw()
		
		return vals

	def selectNextFrontier(self):
		return self.frontierMap.selectNextFrontier()

	def isFrontier(self):
		return self.frontierMap.isFrontier()

	def getNodeOccMap(self, nodeNum):
		localNode = self.poseGraph.get_node_attributes(nodeNum)
		return localNode.getOccMap()

	def getNodePose(self, nodeNum):
		localNode = self.poseGraph.get_node_attributes(nodeNum)
		return localNode.getEstPose()
	
	def setNodePose(self, nodeNum, estPose):
		localNode = self.poseGraph.get_node_attributes(nodeNum)
		localNode.setEstPose(estPose)
	
	def obstCallBack(self, direction):
		self.currNode.obstCallBack(direction)
	
	def saveLocalMap(self):
		
		if self.currNode != 0:
			self.currNode.saveMap()
		
		#for i in range(1, self.numNodes):
		#	localNode = self.poseGraph.get_node_attributes(i)
		#	localNode.saveMap()
			
	def saveMap(self):

		" 1. for each node, copy local image "
		" 2. for each point in image, transform appropriately and plot to main image "
		"""
		self.mapImage = Image.new('L', (self.numPixel, self.numPixel), 127)
		self.image = self.mapImage.load()
		
		for i in range(self.numNodes):
			localNode = self.poseGraph.get_node_attributes(i)
			#localNode.saveMap()
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
		
		self.mapImage.save(self.fileName % self.saveCount)	
		"""
		
		self.occMap.saveMap()
		#self.saveCount += 1	
		
		self.boundMap.saveMap()
		
		"""
		" build global boundary map "
		self.boundMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.boundImage = self.boundMapImage.load()
		
		image = self.occMap.image
		
		for i in range(1,self.numPixel-1):
			for j in range(1,self.numPixel-1):

				isBoundary = False

				# in shadow or unexplored, so it could be a boundary
				if image[i,j] <= 127:

					# check if any neighbors are free space
					for m in range(-1,2):
						for n in range(-1,2):
							if image[i+m,j+n] == 255:
								isBoundary = True
				
				elif image[i,j] == 255:
					self.boundImage[i,j] = 127

				# if a neighbor is free space, then this is a boundary point
				if isBoundary:
					self.boundImage[i,j] = 255		

		self.boundMapImage.save("mapBoundGraph%04u.png" % self.saveCount)	
		"""


		self.obstMapImage.save("mapObstGraph%04u.png" % self.saveCount)	

		self.frontierMap.saveMap()			
		
		"""
		self.frontMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.frontImage = self.frontMapImage.load()


		sumRange = 3
		num = (2*sumRange+1) * (2*sumRange+1)
		
		xMax = 0
		yMax = 0
		densityMax = 0
		densityMin = 1e100

		inhibitRadius = 1.0
		inhibits = [[-3.0,0.0], [-2.0,0.0], [-1.0,0.0]]
		
		for i in range(sumRange,self.numPixel-sumRange):
			for j in range(sumRange,self.numPixel-sumRange):
				if self.boundImage[i,j] == 255:
	
					isTooClose = False
					xReal , yReal = self.gridToReal([i,j])

					" make sure we only consider points that are outside inhibition circles "
					for k in range(len(inhibits)):
						iP = inhibits[k]
						dist = sqrt((iP[0]-xReal)**2 + (iP[1]-yReal)**2)
						if dist < inhibitRadius:
							isTooClose = True
							break

					if not isTooClose:

						fSum = 0
						for m in range(-sumRange, sumRange+1):
							for n in range(-sumRange, sumRange+1):
								fSum += self.obstImage[i+m, j+n]
						
						density =  fSum/num
						density = 255 - density
						self.frontImage[i,j] = density
						
						if density > densityMax:
							xMax = i
							yMax = j
							densityMax = density
							
						if density < densityMin:
							densityMin = density

		maxDen = densityMax - densityMin

		for i in range(sumRange,self.numPixel-sumRange):
			for j in range(sumRange,self.numPixel-sumRange):
				if self.boundImage[i,j] == 255:

					isTooClose = False
					xReal , yReal = self.gridToReal([i,j])

					" make sure we only consider points that are outside inhibition circles "
					for k in range(len(inhibits)):
						iP = inhibits[k]
						dist = sqrt((iP[0]-xReal)**2 + (iP[1]-yReal)**2)
						if dist < inhibitRadius:
							isTooClose = True
							break

					if not isTooClose:
						density =  self.frontImage[i,j]
						density = density - densityMin
						result = 255*(1.0-float(maxDen-density)/float(maxDen))
						self.frontImage[i,j] = result
						#print maxDen, density, result, self.frontImage[i,j]

		maxSum = 0
		xMax = 0
		yMax = 0
		frontWidth = 4

		for i in range(frontWidth,self.numPixel-frontWidth):
			for j in range(frontWidth,self.numPixel-frontWidth):

				sum = 0
				for k in range(i-frontWidth,i+frontWidth):
					for l in range(j-frontWidth,j+frontWidth):						
						sum += self.frontImage[k,l]

				if sum > maxSum:
					maxSum = sum
					xMax = k
					yMax = l
				

		#print "maximum density =", densityMax
		#print "minimum density =", densityMin

		draw = ImageDraw.Draw(self.frontMapImage)
		draw.ellipse((xMax-10, yMax-10, xMax+10, yMax+10), outline=255)

		for k in range(len(inhibits)):
			cP = inhibits[k]
			radius = inhibitRadius
			radius = int(radius*self.divPix)
			xGrid , yGrid = self.realToGrid(cP)
			draw.ellipse((xGrid-radius, yGrid-radius, xGrid+radius, yGrid+radius), outline=255)
		
		self.frontMapImage.save("mapFrontierGraph%04u.png" % self.saveCount)	
		"""						
		
		#self.voronoiMap = VoronoiMap(self.mapImage, self.boundMapImage)
		self.voronoiMap.saveMap()			

		self.navRoadMap = NavRoadMap(self.probe, self.voronoiMap.getGraph(), localNode = self.currNode)
			
		self.saveCount += 1	
			
	def realToGrid(self, point):
		indexX = int(floor(point[0] * self.divPix)) + self.numPixel / 2 + 1
		indexY = int(floor(point[1] * self.divPix)) + self.numPixel / 2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0, (j - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0]
		return point			

	def computeAlphaBoundary(self):
		
		" 1. pick out the points "
		numPixel = self.numPixel
		mapImage = self.mapImage
		image = mapImage.load()
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					pnt = self.gridToReal([j,k])
					points.append(pnt)

		print len(points), "points"
		
		if len(points) > 0:

			a_vert = self.computeAlpha2(points)
	
			" cut out the repeat vertex "
			a_vert = a_vert[:-1]

			a_vert = self.convertAlphaUniform(a_vert)
		
		else:
			a_vert = []
		
		return a_vert

	def convertAlphaUniform(self, a_vert):
		
		" make the vertices uniformly distributed "
		
		new_vert = []
		
		max_spacing = 0.04
		#max_spacing = 0.1
		
		for i in range(len(a_vert)):
			p0 = a_vert[i]
			p1 = a_vert[(i+1) % len(a_vert)]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			vec[0] /= dist
			vec[1] /= dist
			
			new_vert.append(copy(p0))
			
			if dist > max_spacing:
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					new_vert.append(newP)
					
		return new_vert
		
	def computeAlpha2(self, points):
		
		try:
			numPoints = len(points)
			inputStr = str(numPoints) + " "
			
			" radius 0.2 "
			#inputStr += str(0.8) + " "
			inputStr += str(0.2) + " "
			
			for p in points:
				p2 = copy(p)
				p2[0] += gauss(0.0,0.0001)
				p2[1] += gauss(0.0,0.0001)
	
				inputStr += str(p2[0]) + " " + str(p2[1]) + " "
			
			inputStr += "\n"
			
			" start the subprocess "
			
			" send input and receive output "
			sout, serr = subProc.communicate(inputStr)
		
			" convert string output to typed data "
			sArr = sout.split(" ")
	
			numVert = int(sArr[0])
			
			vertices = []
			for i in range(numVert+1):
				vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
			
			return vertices
		
		except:
			print "computeAlpha2() failed"
			print numPoints, "input points"
			print "rawOutput ="
			print sout
			print "numVert =", numVert
			print numVert, len(sArr)
			print sArr
			raise
		
		
		
		
				