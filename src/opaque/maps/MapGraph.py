
from numpy import *
from scipy.optimize import *
import graph
import csv

from Map import Map
from LocalNode import LocalNode
from VoronoiMap import VoronoiMap
from FrontierMap import FrontierMap
from OccupancyMap import OccupancyMap
from FreeSpaceBoundaryMap import FreeSpaceBoundaryMap
from NavRoadMap import NavRoadMap
from StablePose import StablePose
from Pose import Pose
import gen_icp
from functions import *


import pylab
from matplotlib.patches import Circle

from math import *
from copy import copy

import Image
import ImageDraw

# Map Space Parameters
PIXELSIZE = 0.05
MAPSIZE = 20.0

class MapGraph:

	def __init__(self, probe, contacts):
		
		# 
		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)
		
		self.initPose = self.probe.getActualJointPose(19)
		self.poseGraph = graph.digraph()
		self.numNodes = 0
		self.currNode = 0

		self.pixelSize = PIXELSIZE
		self.mapSize = MAPSIZE
		self.numPixel = int(2.0 * self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize / 2.0
		self.divPix = floor((2.0 * self.mapSize / self.pixelSize) / self.mapSize)
		self.saveCount = 0
		self.fileName = "mapGraph%04u.png"
		self.boundIteration = 0


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
		
		self.boundParentNode = 0

		#self.newNode()

	def loadNext(self):
		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, self.probe.robotParam, inSim = False)
		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, inSim = False)
		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize)

		self.currNode.readFromFile(self.numNodes)
		
		self.poseGraph.add_node(self.numNodes, self.currNode)
		if self.numNodes > 0:
			self.addConstraint(self.numNodes-1, self.numNodes)

		self.numNodes += 1	

	def loadFile(self, dirName, num_poses):
		
		self.poseGraph = graph.digraph()
		self.numNodes = 0
		self.currNode = 0
		
		for i in range(0,num_poses):
			
			print "loading file", i
			
			# def __init__(self, probe, contacts, nodeID, rootNode, pixelSize, inSim = True):
			#self.currNode = LocalNode(self.probe, self.contacts, i, 19, self.pixelSize, inSim = False)
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, self.pixelSize)
			self.currNode.readFromFile(dirName, i)
			
			print "adding node", i
			self.poseGraph.add_node(i, self.currNode)
			if self.numNodes > 0:
				self.addConstraint(i-1, i)

			self.numNodes += 1
			
	def addConstraint(self, i, j):

		if self.numNodes > 0:		
			
			print "adding constraint", i, j
			
			pose1 = copy(self.poseGraph.get_node_attributes(i).getEstPose())
			pose2 = copy(self.poseGraph.get_node_attributes(j).getEstPose())
		
			pose2[0] -= pose1[0]
			pose2[1] -= pose1[1]
			pose2[2] -= pose1[2]
			pose2[2] = normalizeAngle(pose2[2])
		
			ang = pose1[2]
		
			xRot = pose2[0] * cos(ang) + pose2[1] * sin(ang)
			yRot = -pose2[0] * sin(ang) + pose2[1] * cos(ang)
		
			pose2[0] = xRot
			pose2[1] = yRot
			
			transform = matrix([[pose2[0]], [pose2[1]], [pose2[2]]])
			covE = matrix([[ 0.01866183, 0.00099555, 0.0004255 ],
				[ 0.00099555, 0.00145845, -0.00017733],
				[ 0.0004255,  -0.00017733,  0.0003789 ]])
	
			
			self.poseGraph.add_edge(i, j, attrs = [transform, covE])

	def setCenterPoints(self, centerPoints):
		self.currNode.setCenterPoints(centerPoints)
		
	def getCenterPoints(self):
		return self.currNode.getCenterPoints()
	
	def getCurrentNode(self):
		return self.currNode
	
	def newNode(self):
		
		print "checkN"
		
		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		print "checkO"

		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19)
		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, inSim = False)
		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize)

		print "checkP"
		
		self.poseGraph.add_node(self.numNodes, self.currNode)
		if self.numNodes > 0:
			self.addConstraint(self.numNodes-1, self.numNodes)
	
		print "checkQ"
		
		self.numNodes += 1
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

		print "checkR"
	
	" FIXME:  Perhaps this should be moved to the PokeWalls behavior since this the only thing that requires it"
	def keepStablePose(self):
		self.stablePose.update()
	
	def correctPoses2(self):

		if self.numNodes < 2:
			return

		print "check_A"
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

		print len(estPoses), "poses"
		print "old estPose:", estPoses[-1]

		print "check_B"

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

		print "check_C"
	
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
			

		print "check_D"
	
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

		print "check_E"
	
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

		print "check_F"

		pastCircles = []
		for m in range(0,len(estPoses)-1):
			hull = a_hull_trans[m]
			radius, center = gen_icp.computeEnclosingCircle(hull)
			pastCircles.append([radius,center])
			
		pastPose = estPoseOrigin
		
		pastHull = gen_icp.computeUnions(a_hull_trans)
		gen_icp.addPointToLineCovariance(pastHull, high_var=1.0, low_var=0.001)
		#gen_icp.addDistanceFromOriginCovariance(pastHull, tan_var=0.1, perp_var=0.01)

		print "check_G"

		" run generalized ICP (a plot is made for each iteration of the algorithm) "
		offset = gen_icp.gen_ICP(pastPose, targetPose, pastHull, targetHull, pastCircles, costThresh, minMatchDist, plotIteration)
		
		offsets[k-1] = offset
		
		" recompute the estimated poses with the new offset "
		newEstPoses = []
		newEstPoses.append(estPoses[0])
		
		for m in range(len(offsets)):
			pose1 = Pose(newEstPoses[m])
			offset = offsets[m]
			newEstPose2 = pose1.convertLocalOffsetToGlobal(offset)
			newEstPoses.append(newEstPose2)
			
		print "check_H"
			
		estPoses = newEstPoses
		
		" update the estimated poses "
		for m in range(0,len(estPoses)):
			self.setNodePose(m, estPoses[m])
			
		print len(estPoses), "poses"
		print "new estPose:", estPoses[-1]
		
		" update the current estimated pose in AverageContacts "
		self.contacts.resetPose(estPose = estPoses[-1])

	def dijkstra_proj(self):
		
		paths = {}
		
		visited = [False for i in range(self.numNodes)]
		optimals = {}
		
		" initial node 0 "
		visited[0] = True
		
		" for each edge out of 0, add to path "
		neighbors = self.poseGraph.neighbors(0)
		incidents = self.poseGraph.incidents(0)
		
		for neigh in neighbors:
			if not visited[neigh]:
				transform, covE = self.poseGraph.get_edge_attributes(0,neigh)
				try:
					paths[neigh]
				except:
					paths[neigh] = []
				
				paths[neigh].append([transform, covE])

		for incid in incidents:
			if not visited[incid]:
				transform, covE = self.poseGraph.get_edge_attributes(incid, 0)
				try:
					paths[incid]
				except:
					paths[incid] = []
				
				" TODO: invert the transform and covariance "
				#paths[incid].append([transform, covE])
		
		
		" while some nodes unvisited "
		while visited.count(False) > 0:
			print visited
			" find the minimum uncertainty path p "
			minUnc = 1e100
			dest = 0
			minTrans = 0
			minCov = 0
			#print paths
			for key, value in paths.iteritems():
				#print "key:", key
				#print "value:", value
				if not visited[key]:
					for p in value:
						#print "p:", p
						uncertainty = linalg.det(p[1]) 
						if uncertainty < minUnc:
							minUnc = uncertainty
							dest = key
							minTrans = p[0]
							minCov = p[1]
					
			" mark this as visited and record the optimal path "
			visited[dest] = True
			optimals[dest] = [minTrans, minCov]
		
			" for all edges leaving dest, add the composed path "

			" for each edge out of dest, add to path "
			neighbors = self.poseGraph.neighbors(dest)
			incidents = self.poseGraph.incidents(dest)
			
			print "neighbors of", dest, "=", neighbors
			
			for neigh in neighbors:
				if not visited[neigh]:
					print dest, neigh
					transform, covE = self.poseGraph.get_edge_attributes(dest,neigh)
					try:
						paths[neigh]
					except:
						paths[neigh] = []

					T_b_a = optimals[dest][0]
					T_c_b = transform
					
					E_a_b = optimals[dest][1]
					E_b_c = covE
					
					x1 = T_b_a[0,0]
					y1 = T_b_a[1,0]
					p1 = T_b_a[2,0]
					
					x2 = T_c_b[0,0]
					y2 = T_c_b[1,0]
					p2 = T_c_b[2,0]
					
					T_c_a = matrix([[x1 + x2*cos(p1) - y2*sin(p1)], [y1 + x2*sin(p1) + y2*cos(p1)], [p1+p2]])
					
					J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]])
					J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]])
					
					E_a_c = J1 * E_a_b * J1.T + J2 * E_b_c * J2.T

					paths[neigh].append([T_c_a, E_a_c])
	
			for incid in incidents:
				if not visited[incid]:
					transform, covE = self.poseGraph.get_edge_attributes(incid, dest)
					try:
						paths[incid]
					except:
						paths[incid] = []
					
					" TODO: invert the transform and covariance "
					#paths[incid].append([transform, covE])
					
		
			" compute the set of pairwise poses to consider using djikstra projection "
		
		print paths	

		#for i in range(self.numNodes):
		#	for j in range(i+1, self.numNodes):
		#		trans, cov = self.findMinPath(i,j)

	def correctPoses3(self):

		if self.numNodes < 2:
			return

		print "check_A"
		if self.currNode.isDirty():
			self.currNode.synch()
			#self.synch()
			
		self.dijkstra_proj()
		
		return


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

		print len(estPoses), "poses"
		print "old estPose:", estPoses[-1]

		print "check_B"

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

		print "check_C"
	
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
			

		print "check_D"
	
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

		print "check_E"
	
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

		print "check_F"

		pastCircles = []
		for m in range(0,len(estPoses)-1):
			hull = a_hull_trans[m]
			radius, center = gen_icp.computeEnclosingCircle(hull)
			pastCircles.append([radius,center])
			
		pastPose = estPoseOrigin
		
		pastHull = gen_icp.computeUnions(a_hull_trans)
		gen_icp.addPointToLineCovariance(pastHull, high_var=1.0, low_var=0.001)
		#gen_icp.addDistanceFromOriginCovariance(pastHull, tan_var=0.1, perp_var=0.01)

		print "check_G"

		" run generalized ICP (a plot is made for each iteration of the algorithm) "
		offset = gen_icp.gen_ICP(pastPose, targetPose, pastHull, targetHull, pastCircles, costThresh, minMatchDist, plotIteration)
		
		offsets[k-1] = offset
		
		" recompute the estimated poses with the new offset "
		newEstPoses = []
		newEstPoses.append(estPoses[0])
		
		for m in range(len(offsets)):
			pose1 = Pose(newEstPoses[m])
			offset = offsets[m]
			newEstPose2 = pose1.convertLocalOffsetToGlobal(offset)
			newEstPoses.append(newEstPose2)
			
		print "check_H"
			
		estPoses = newEstPoses
		
		" update the estimated poses "
		for m in range(0,len(estPoses)):
			self.setNodePose(m, estPoses[m])
			
		print len(estPoses), "poses"
		print "new estPose:", estPoses[-1]
		
		" update the current estimated pose in AverageContacts "
		self.contacts.resetPose(estPose = estPoses[-1])


	def drawEstBoundary(self):

		if self.boundParentNode == 0:
			self.boundParentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("globalBoundNode")

		
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

		pylab.clf()

		print "check_Z"
		print "self.contacts.activeRef =", self.contacts.activeRef

		if self.currNode != 0 and self.contacts.numRef > 0:
			#onSegPose = Pose(self.contacts.getClosestPose(self.currNode.rootNode))
			onSegPose = Pose(self.contacts.getAveragePose(self.currNode.rootNode))
			onSegActPose = Pose( self.probe.getActualJointPose(self.currNode.rootNode))

		plotPoints = [[], []]
		for i in range(len(points)):
			pnt = points[i]
			
			if self.currNode != 0 and self.contacts.numRef > 0:
				newPnt = copy(pnt)				
				localPnt = onSegPose.convertGlobalToLocal(pnt)
				pnt = onSegActPose.convertLocalToGlobal(localPnt)

			plotPoints[0].append(pnt[0])
			plotPoints[1].append(pnt[1])

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

		if len(plotPoints[0]) > 0:
			pylab.scatter(plotPoints[0],plotPoints[1], linewidth=1, color='k')

		print "check_Y"
		print "self.contacts.activeRef =", self.contacts.activeRef
		
		plotPoints = [[], []]
		for i in range(39):
			pnt = self.contacts.getAveragePose(i)
			
			plotPoints[0].append(pnt[0])
			plotPoints[1].append(pnt[1])
						
			childNode = self.boundParentNode.createChildSceneNode("globalPosePoint" + "_" + str(i))
			self.childNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("globalPoseEnt" + "_" + str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Red")
			self.childEntities.append(currEntity)
	
			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			childNode.setPosition(position)
	
			size = ogre.Vector3(0.06,0.06,0.06)
			childNode.setScale(size)
	
			childNode.attachObject(currEntity)	

		if len(plotPoints[0]) > 0:
			pylab.scatter(plotPoints[0],plotPoints[1], linewidth=0, color=(1.0,0.6,0.6))

		print "check_X"
		print "self.contacts.activeRef =", self.contacts.activeRef

		plotPoints = [[], []]
		for i in range(39):
			pnt = self.probe.getActualJointPose(i)
			
			plotPoints[0].append(pnt[0])
			plotPoints[1].append(pnt[1])

		if len(plotPoints[0]) > 0:
			pylab.scatter(plotPoints[0],plotPoints[1], linewidth=0, color=(0.6,0.6,1.0))


		print "check_V"
		print "self.contacts.activeRef =", self.contacts.activeRef

		plotPoints = [[], []]
		for i in range(39):
			if self.contacts.activeRef[i]:
				pnt = self.contacts.getClosestPose(i)
				
				plotPoints[0].append(pnt[0])
				plotPoints[1].append(pnt[1])
		
		if len(plotPoints[0]) > 0:
			pylab.scatter(plotPoints[0],plotPoints[1], linewidth=0, color=(0.6,1.0,0.6))

		print "check_U"
		print "self.contacts.activeRef =", self.contacts.activeRef
			
		self.plotEnv()	

		pylab.xlim(-4,9)
		pylab.ylim(-7,5)
		pylab.savefig("inSituBoundary_%04u.png" % self.boundIteration)
		pylab.clf()
		
		self.boundIteration += 1
						
	def plotEnv(self):
		
		walls = self.probe.getWalls()
	
		for wall in walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = copy(wall[i])
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g')
			
			
	def synch(self):
		self.currNode.synch()
		
		print "checkA"

		self.occMap.update()
		print "checkB"
		self.boundMap.update(self.occMap)

		print "checkC"
		
		" 1. for each node, copy local image "
		" 2. for each point in image, transform appropriately and plot to main image "
		self.obstMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.obstImage = self.obstMapImage.load()

		for i in range(self.numNodes):
			localNode = self.poseGraph.get_node_attributes(i)
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

						self.obstImage[indexX, indexY] = 255
							
		print "checkD"

		
		self.frontierMap.update()
		
		print "checkE"

		self.voronoiMap.update()

		print "checkF"

	def forceUpdate(self, isForward=True):

		self.stablePose.setDirection(isForward)
		self.currNode.update(isForward)

	def update(self, isForward=True):

		self.stablePose.setDirection(isForward)

		" if the configuration is a stable pose, the maximum displacement is not exceeded "
		if self.stablePose.isStable():
			self.currNode.update(isForward)
		else:
			pass
		
	def computeHeadPath(self, currPose, frontierPoint, exploreRoot):
		vals = self.navRoadMap.computeHeadPath(currPose, frontierPoint, exploreRoot)
	
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

		" save the local map "
		if self.currNode != 0:
			self.currNode.saveMap()

	def saveMap(self):
		
		" save the global maps to file "

		print "checkG"
		
		self.occMap.saveMap()

		print "checkH"
		
		self.boundMap.saveMap()
		
		print "checkI"

		self.obstMapImage.save("mapObstGraph%04u.png" % self.saveCount)	

		print "checkJ"

		self.frontierMap.saveMap()			
		
		print "checkK"

		self.voronoiMap.saveMap()			

		print "checkL"

		self.navRoadMap = NavRoadMap(self.mapSize, self.probe, self.voronoiMap.getGraph(), localNode = self.currNode)

		print "checkM"
			
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

		
		
		
				