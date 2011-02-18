
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

		if False:
			for k in range(1):
				pylab.clf()
				
				#[-4.0, 0.2] ,[-14.0, 0.2]
				#[-4.0, -0.2],[-14.0, -0.2]
				
				#xP = [-4.,-14.0]
				#yP = [0.2, 0.2]
				#pylab.plot(xP, yP, color='k')
				#yP = [-0.2, -0.2]
				#pylab.plot(xP, yP, color='k')
	
				self.plotEnv()
	
				theta = k*0.05
				 
				for i in range(self.probe.numSegs):
					
		
					pose = self.probe.getActualSegPose(i, phi = theta)
					#pose = self.probe.getActualJointPose(i, phi = theta)
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					
					p1 = [xTotal + 0.5*segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal - 0.5*segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal - 0.5*segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1], p2[1], p3[1], p4[1], p1[1]]
	
					pylab.plot(xP,yP, color='b')
					
				#for i in range(-1,20):
				for i in range(-1,self.probe.numSegs-1):
					
					pose = self.probe.getActualJointPose(i, phi = theta)
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1], p2[1], p3[1], p4[1], p1[1]]
		
					if i == 19:
						pylab.plot(xP,yP, color='r', linewidth='2')
						pylab.scatter([pose[0]], [pose[1]], color='r', linewidth=8)
					else:
						pylab.plot(xP,yP, color='r')				
						pylab.scatter([pose[0]], [pose[1]], color='r', linewidth=1)
	
				rootPose = self.probe.getActualJointPose(19)
	
				for i in range(-1,self.probe.numSegs-1):
				#for i in range(-1,20):
					
					#pose = self.probe.getActualJointPose(i, phi = theta)
					pose = self.probe.getJointPose(rootPose, 19, i)
					
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1], p2[1], p3[1], p4[1], p1[1]]
	
					if i == 19:
						pylab.plot(xP,yP, color='k', linewidth='2')
						pylab.scatter([pose[0]], [pose[1]], color='k', linewidth=4)
					else:
						pylab.plot(xP,yP, color='k')
						pylab.scatter([pose[0]], [pose[1]], color='k', linewidth=1)
	
				for i in range(-1,self.probe.numSegs-1):
					
					#pose = self.probe.getJointWRTJointPose([0.0,0.0,0.0], 19, i)
					pose = self.probe.getJointPose([0.0,0.0,0.0], 19, i)
					
					xTotal = pose[0]
					zTotal = pose[1]
					totalAngle = pose[2]
					segWidth = self.probe.segWidth
					segLength = self.probe.segLength
					
					p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
					p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
					p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
					p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
	
					xP = [p1[0], p2[0], p3[0], p4[0], p1[0]]
					yP = [p1[1]+3.0, p2[1]+3.0, p3[1]+3.0, p4[1]+3.0, p1[1]+3.0]
	
					if i == 19:
						pylab.plot(xP,yP, color='g', linewidth='2')
						pylab.scatter([pose[0]], [pose[1]+3.0], color='g', linewidth=4)
					else:
						pylab.plot(xP,yP, color='k')
						pylab.scatter([pose[0]], [pose[1]+3.0], color='g', linewidth=1)
									
				pylab.xlim(-4,4)
				pylab.ylim(-4,4)
				pylab.savefig("testplot%04u.png" % k)
				pylab.show()

		#self.newNode()

	def loadNext(self):
		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, self.probe.robotParam, inSim = False)
		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, inSim = False)
		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize)

		self.currNode.readFromFile(self.numNodes)
		
		self.poseGraph.add_node(self.numNodes, self.currNode)
		if self.numNodes > 0:
			self.addMotionConstraint(self.numNodes-1, self.numNodes)

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
				#self.addMotionConstraint(i-1, i)
				self.addNaiveMotionConstraint(i-1, i, 0.14, True)

			self.numNodes += 1

		#pylab.clf()			
		#for i in range(num_poses):
			
			#currNode = self.poseGraph.get_node_attributes(i)
			
			#localPosture = currNode.localPosture
			#splPoints = currNode.splPoints

			#print localPosture
			#print splPoints
		
			#xP = []
			#yP = []
			#for p in localPosture:
			#	xP.append(p[0])
			#	yP.append(p[1])
			
			#pylab.plot(xP,yP)
			
			#xP = []
			#yP = []
			#for p in splPoints:
			#	xP.append(p[0])
			#	yP.append(p[1])
			#pylab.plot(xP,yP)

			#pylab.savefig("plotCenter%04u.png" % i)
			#pylab.clf()		
			
	def addMotionConstraint(self, i, j):

		if self.numNodes > 0:		
			
			print "adding constraint", i, j
			
			pose1 = copy(self.poseGraph.get_node_attributes(i).getEstPose())
			pose2 = copy(self.poseGraph.get_node_attributes(j).getEstPose())

			xA = pose1[0]
			yA = pose1[1]
			pA = pose1[2]
			
			xB = pose2[0]
			yB = pose2[1]
			pB = pose2[2]

			xT = cos(pA)*(xB-xA) + sin(pA)*(yB-yA)
			yT = -sin(pA)*(xB-xA) + cos(pA)*(yB-yA)
			pT = pB - pA
			
			#pose2[0] -= pose1[0]
			#pose2[1] -= pose1[1]
			#pose2[2] -= pose1[2]
			#pose2[2] = normalizeAngle(pose2[2])
		
			#ang = pose1[2]
		
			#xRot = pose2[0] * cos(ang) + pose2[1] * sin(ang)
			#yRot = -pose2[0] * sin(ang) + pose2[1] * cos(ang)
		
			#pose2[0] = xRot
			#pose2[1] = yRot
			
			#transform = matrix([[pose2[0]], [pose2[1]], [pose2[2]]])
			
			transform = matrix([[xT], [yT], [pT]])
			covE = matrix([[ 0.1, 0.0, 0.0 ],
							[ 0.0, 0.1, 0.0],
							[ 0.0, -0.0, pi/4.0 ]])

			#covE = matrix([[ 0.01866183, 0.00099555, 0.0004255 ],
			#	[ 0.00099555, 0.00145845, -0.00017733],
			#	[ 0.0004255,  -0.00017733,  0.0003789 ]])





	def addNaiveMotionConstraint(self, i, j, stepDist, direction):

		if self.numNodes > 0:		

			" do stepDist in the direction of the centerline "
			if direction:
				firstGuess = self.makeGuess(i, j, stepDist)
			else:
				firstGuess = self.makeGuess(i, j, -stepDist)

			print
			print "adding naive motion constraint", i, j
			print firstGuess
			print "dist =", sqrt(firstGuess[0]*firstGuess[0] + firstGuess[1]*firstGuess[1])

			pose1 = copy(self.poseGraph.get_node_attributes(i).getGndPose())
			pose2 = copy(self.poseGraph.get_node_attributes(j).getGndPose())

			xA = pose1[0]
			yA = pose1[1]
			pA = pose1[2]
			
			xB = pose2[0]
			yB = pose2[1]
			pB = pose2[2]

			#xOff = stepDist * cos(angle1)
			#yOff = stepDist * sin(angle1)

			xT = cos(pA)*(xB-xA) + sin(pA)*(yB-yA)
			yT = -sin(pA)*(xB-xA) + cos(pA)*(yB-yA)
			pT = pB - pA
			
			gndOffset = [xT, yT, pT]
			
			print gndOffset
			print "dist =", sqrt(xT*xT + yT*yT)
			
			" FIXME:  add covariance from vector in guess direction "
			transform = matrix([[firstGuess[0]], [firstGuess[1]], [firstGuess[2]]])
			covE = matrix([[ 0.1, 0.0, 0.0 ],
							[ 0.0, 0.1, 0.0],
							[ 0.0, -0.0, pi/4.0 ]])
			
			self.poseGraph.add_edge(i, j, attrs = [transform, covE])

			" ################ "

			pylab.clf()
			currNode1 = self.poseGraph.get_node_attributes(i)
			currNode2 = self.poseGraph.get_node_attributes(j)
			
			
			localPosture = currNode1.localPosture
			splPoints = currNode1.splPoints

			poseProfile1 = Pose(pose1)

			xP = []
			yP = []
			for p in localPosture:
				#nP = p
				nP = poseProfile1.convertLocalToGlobal(p)

				xP.append(nP[0])
				yP.append(nP[1])
			pylab.plot(xP,yP, color='b')
			
			xP = []
			yP = []
			for p in splPoints:
				#nP = p
				nP = poseProfile1.convertLocalToGlobal(p)

				xP.append(nP[0])
				yP.append(nP[1])
			pylab.plot(xP,yP, color='b')


			localPosture = currNode2.localPosture
			splPoints = currNode2.splPoints

			estPose2 = poseProfile1.convertLocalOffsetToGlobal(gndOffset) 
			poseProfile2 = Pose(estPose2)
			poseProfile2 = Pose(pose2)
			print pose2, estPose2


			xP = []
			yP = []
			for p in localPosture:
				#nP = gen_icp.dispOffset(p,firstGuess)
				#nP = gen_icp.dispOffset(p,gndOffset)
				nP = poseProfile2.convertLocalToGlobal(p)
				xP.append(nP[0])
				yP.append(nP[1])
			pylab.plot(xP,yP, color='r')
			
			xP = []
			yP = []
			for p in splPoints:
				#nP = gen_icp.dispOffset(p,firstGuess)
				#nP = gen_icp.dispOffset(p,gndOffset)
				nP = poseProfile2.convertLocalToGlobal(p)
				xP.append(nP[0])
				yP.append(nP[1])
			pylab.plot(xP,yP, color='r')

			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("plotCenter%04u.png" % self.numNodes)
			pylab.clf()		

	def setCenterPoints(self, centerPoints):
		self.currNode.setCenterPoints(centerPoints)
		
	def getCenterPoints(self):
		return self.currNode.getCenterPoints()
	
	def getCurrentNode(self):
		return self.currNode
	
	def newNode(self, stepDist, direction):
		
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
			self.addMotionConstraint(self.numNodes-1, self.numNodes)
			#self.addNaiveMotionConstraint(self.numNodes-1, self.numNodes, stepDist, direction)
	
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

	def dijkstra_proj(self, initNode = 0):
		
		paths = {}
		
		visited = [False for i in range(self.numNodes)]
		distances = [Inf for i in range(self.numNodes)]
		optimals = {}
		
		" initial node 0 "
		distances[initNode] = 0.0
		visited[initNode] = True
		
		" for each edge out of 0, add to path "
		neighbors = self.poseGraph.neighbors(initNode)
		incidents = self.poseGraph.incidents(initNode)
		
		#[[ 0.1         0.          0.        ]
		# [ 0.          0.1         0.        ]
		# [ 0.         -0.          0.78539819]]		
		
		
		for neigh in neighbors:
			if not visited[neigh]:
				transform, covE = self.poseGraph.get_edge_attributes(initNode,neigh)
				
				dist = linalg.det(covE)
				if dist < distances[neigh]:
					paths[neigh] = [transform, covE]
					distances[neigh] = dist


		for incid in incidents:
			if not visited[incid]:
				transform, covE = self.poseGraph.get_edge_attributes(incid, initNode)
				
				" TODO: invert the transform "
				#paths[incid].append([transform, covE])
				#dist = linalg.det(covE)
				#if dist < distances[incid]:
				#	distances[incid] = linalg.det(covE)
				
				#print "incid:", incid
				#print transform
				
				#x2 = transform[0,0]
				#y2 = transform[1,0]
				#p2 = transform[2,0]
				
				#x3 = -x2*cos(-p2) + y2*sin(-p2)
				#y3 = -x2*sin(-p2) - y2*cos(-p2)
				#p3 = -p2
										
				#newTransform = matrix([[x3], [y3], [p3]])

				xA = transform[0,0]
				yA = transform[1,0]
				pA = transform[2,0]
				
				x1 = -cos(pA)*xA - sin(pA)*yA
				y1 = sin(pA)*xA - cos(pA)*yA
				p1 = normalizeAngle(-pA)

				newTransform = matrix([[x1], [y1], [p1]])

				paths[incid] = [newTransform, covE]
				dist = linalg.det(covE)
				if dist < distances[incid]:
					paths[incid] = [newTransform, covE]
					distances[incid] = dist


		" while some nodes unvisited "
		while visited.count(False) > 0:

			" find the minimum uncertainty path p "

			minDist = Inf
			minDest = -1
			#print visited
			#print distances
			for i in range(self.numNodes):
				if not visited[i]:
					if distances[i] < minDist:
						minDist = distances[i]
						minDest = i

			
			" an unvisited node is unreachable "
			if minDest == -1:
				break
			else:
				dest = minDest
				minUnc = Inf
				minTrans = 0
				minCov = 0
				
				p = paths[dest]
				uncertainty = linalg.det(p[1]) 
				if uncertainty < minUnc:
					minUnc = uncertainty
					minTrans = p[0]
					minCov = p[1]
						
				" mark this as visited and record the optimal path "
				visited[dest] = True
				optimals[dest] = [minTrans, minCov]
			
				" for all edges leaving dest, add the composed path "
	
				" for each edge out of dest, add to path "
				neighbors = self.poseGraph.neighbors(dest)
				incidents = self.poseGraph.incidents(dest)
				
				#print "neighbors of", dest, "=", neighbors
				
				for neigh in neighbors:
					if not visited[neigh]:
						#print dest, neigh
						transform, covE = self.poseGraph.get_edge_attributes(dest,neigh)
	
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
						
						T_c_a = matrix([[x1 + x2*cos(p1) - y2*sin(p1)],
										[y1 + x2*sin(p1) + y2*cos(p1)],
										[p1+p2]])
						
						J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]])
						J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]])
						
						E_a_c = J1 * E_a_b * J1.T + J2 * E_b_c * J2.T
	
						dist = linalg.det(E_a_c)
						if dist < distances[neigh]:
							paths[neigh] = [T_c_a, E_a_c]
							distances[neigh] = dist
		
				for incid in incidents:
					if not visited[incid]:
						transform, covE = self.poseGraph.get_edge_attributes(incid, dest)

						T_b_a = optimals[dest][0]
						T_c_b = transform

						E_a_b = optimals[dest][1]
						E_b_c = covE
						
						x1 = T_b_a[0,0]
						y1 = T_b_a[1,0]
						p1 = T_b_a[2,0]

						xA = T_c_b[0,0]
						yA = T_c_b[1,0]
						pA = T_c_b[2,0]
						
						x2 = -cos(pA)*xA - sin(pA)*yA
						y2 = sin(pA)*xA - cos(pA)*yA
						p2 = -pA
					
						#J1 = matrix([[1,0, x2*sin(p1-p2) + y2*cos(p1-p2)],
						#			[0,1,-x2*cos(p1-p2) + y2*sin(p1-p2)],
						#			[0,0,1]])
						#J2 = matrix([[-cos(p1-p2), sin(p1-p2), -x2*sin(p1-p2) - y2*cos(p1-p2)],
						#			[-sin(p1-p2), -cos(p1-p2), x2*cos(p1-p2) - y2*sin(p1-p2)],
						#			[0, 0, -1]])

						T_c_a = matrix([[x1 + x2*cos(p1) - y2*sin(p1)],
										[y1 + x2*sin(p1) + y2*cos(p1)],
										[p1+p2]])
				
						J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]])
						J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]])

						#J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]])
						#J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]])
						
						E_a_c = J1 * E_a_b * J1.T + J2 * E_b_c * J2.T
	
						dist = linalg.det(E_a_c)
						if dist < distances[incid]:
							paths[incid] = [T_c_a, E_a_c]
							distances[incid] = dist
								
						" TODO: invert the transform and covariance "
						#paths[incid].append([transform, covE])
						#dist = linalg.det(E_a_c)
						#if dist < distances[incid]:
						#	distances[incid] = dist
						
			
				" compute the set of pairwise poses to consider using djikstra projection "
		
		for key, value in paths.iteritems():
			#print key, value
			offset = value[0]
			covE = value[1]
			#print key, ":", offset[0,0], offset[1,0], offset[2,0], linalg.det(covE)
			
		return paths

	def correctPoses3(self):

		DIST_THRESHOLD = 3.0

		def mahab_dist(c1, c2, r1, r2, E):

			" distance between centroids "
			d = c2 - c1
			
			
			s = max(0, sqrt(d[0,0]**2 + d[1,0]**2) - r1 - r2) * d / sqrt(d[0,0]**2 + d[1,0]**2)
			E_inv = linalg.inv(E)
			dist = s.T * E_inv * s
			return dist[0,0]

		if self.numNodes < 2:
			return

		print "check_A"
		if self.currNode.isDirty():
			self.currNode.synch()
			#self.synch()

		" perform dijkstra projection over all nodes"
		paths = []
		for i in range(self.numNodes):
			paths.append(self.dijkstra_proj(i))
		
		" select pairs to attempt a sensor constraint "
		pair_candidates = []
		print "CANDIDATES:"
		for i in range(self.numNodes):
			path_set = paths[i]
			
			ind = range(i+1,self.numNodes)
			#print i, self.numNodes, ind

			for j in ind:				
				loc2 = path_set[j][0]
				
				c1 = matrix([[0.0],[0.0]])
				c2 = matrix([[loc2[0,0]],[loc2[1,0]]])
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]])
				
				dist = mahab_dist(c1, c2, 0.25, 0.25, E)
				
				#print i, j, dist
				if dist <= DIST_THRESHOLD:
					pair_candidates.append([dist,i,j,loc2])
				
			ind = range(0,i)
			ind.reverse()
			
			for j in ind:
				
				loc2 = path_set[j][0]
				
				c1 = matrix([[0.0],[0.0]])
				c2 = matrix([[loc2[0,0]],[loc2[1,0]]])
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]])
				
				dist = mahab_dist(c1, c2, 0.25, 0.25, E)
				
				#print i, j, dist
				if dist <= DIST_THRESHOLD:
					pair_candidates.append([dist,i,j,loc2])

		" remove duplicates "
		pair_unique = []
		is_free = [True for i in range(len(pair_candidates))]
		for i in range(len(pair_candidates)):
			
			if is_free[i]:
				dest1 = pair_candidates[i][0]
				x1 = pair_candidates[i][1]
				x2 = pair_candidates[i][2]
				
				" find the opposing pair if it exists "
				for j in range(i+1, len(pair_candidates)):
					
					if pair_candidates[j][1] == x2 and pair_candidates[j][2] == x1:
						
						dest2 = pair_candidates[j][0]
						
						if dest1 < dest2:
							pair_unique.append(pair_candidates[i])
						else:
							pair_unique.append(pair_candidates[j])
						
						is_free[i] = False
						is_free[j] = False
						
						break

				if is_free[i]:
					pair_unique.append(pair_candidates[i])
					is_free[i] = False

		print "UNIQUE PAIRS"
		#for p in pair_unique:
		#	print p[0], p[1], p[2]

		def decimatePoints(points):
			result = []
		
			for i in range(len(points)):
				#if i%2 == 0:
				if i%4 == 0:
					result.append(points[i])
		
			return result
		
		def computeHull(i):
			
			" Read in data of Alpha-Shapes and add their associated covariances "
			node1 = self.poseGraph.get_node_attributes(i)
			node1.computeAlphaBoundary()			
			a_data = node1.getAlphaBoundary()
			a_data = decimatePoints(a_data)
			
			" treat the points with the point-to-line constraint "
			gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
	
			" treat the points with the distance-from-origin increasing error constraint "
			gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
			
			return a_data
		
		hull_computed = [False for i in range(self.numNodes)]
		a_hulls = [0 for i in range(self.numNodes)]
		
		for p in pair_unique:
			if not hull_computed[p[1]]:
				print "compute hull", p[1]
				a_hulls[p[1]] = computeHull(p[1])
				hull_computed[p[1]] = True

			if not hull_computed[p[2]]:
				print "compute hull", p[2]
				a_hulls[p[2]] = computeHull(p[2])
				hull_computed[p[2]] = True
										
			print p[1], p[2], p[0]

		
		#trans_hulls = [0 for i in range(self.numNodes)]
		#for i in range(self.numNodes):
		#	if hull_computed[i]:
		#		node1 = self.poseGraph.get_node_attributes(i)
		#		estPose1 = node1.getEstPose()
		#		hull1 = a_hulls[i]
		#		hull_trans = []
		#		for p in hull1:
		#			hull_trans.append(gen_icp.dispPoint(p, estPose1))
		#		
		#		trans_hulls[i] = hull_trans
					
		
		#hull_trans1 = []
		#for p in hull1:
		#	hull_trans1.append(gen_icp.dispPoint(p, estPose1))

		hypotheses = []
		
		for p in pair_unique:
			transform = p[3]
			n1 = p[1]
			n2 = p[2]
			offset = self.makeSensorConstraint(transform, n1, n2, a_hulls[n1], a_hulls[n2])
			#offset = self.makeSensorConstraint(transform, n1, n2, trans_hulls[n1], trans_hulls[n2])

			#old = [p[3][0,0],p[3][1,0],p[3][2,0]]
			#print old

			hypotheses.append([p[1],p[2],offset])

		" need more hypotheses "
		" add a more naive motion model constraint "

		for hyp in hypotheses:
			print hyp[0], hyp[1], hyp[2]
			
		
		
		
		#p = pair_unique[-1]
		#print self.makeSensorConstraint(p[0],p[1],p[2], a_hulls[p[1]], a_hulls[p[2]])
		#print p[3]

	def makeGuess(self, n1, n2, stepDist):
		
		node1 = self.poseGraph.get_node_attributes(n1)	
		node2 = self.poseGraph.get_node_attributes(n2)	
		
		" align the root vectors "
		vec1 = node1.centerCurve.getUVector(0.5)
		vec2 = node2.centerCurve.getUVector(0.5)
		
		angle1 = acos(vec1[0])
		if asin(vec1[1]) < 0:
			angle1 = -angle1

		angle2 = acos(vec2[0])
		if asin(vec2[1]) < 0:
			angle2 = -angle2
			
		angle = angle2 - angle1
		
		" center the root node positions "
		
		xOff = stepDist * cos(angle1)
		yOff = stepDist * sin(angle1)
		
		return [xOff, yOff, -angle]
	

	def makeSensorConstraint(self, transform, n1, n2, hull1, hull2):

		print "hypothesizing sensor constraint between", n1, "and", n2

		" TUNE: estimated step distance from empirical results "
		STEP_DIST = 0.145

		" initial guess for x, y, theta parameters "
		firstGuess = self.makeGuess(n1, n2, STEP_DIST)
		print "guess =", firstGuess

		" TUNE ME:  threshold cost difference between iterations to determine if converged "
		costThresh = 0.1
	
		" TUNE ME:   minimum match distance before point is discarded from consideration "
		minMatchDist = 2.0
	
		" plot the best fit at each iteration of the algorithm? "
		plotIteration = True
		#plotIteration = False
	
		#offset = [0.0,0.0,0.0]
	
		" Extract the data from the files and put them into arrays "
		
		node1 = self.poseGraph.get_node_attributes(n1)
		node2 = self.poseGraph.get_node_attributes(n2)
		
		estPose1 = node1.getEstPose()
		estPose2 = node2.getEstPose()
		
		hull_trans1 = []
		for p in hull1:
			hull_trans1.append(gen_icp.dispPoint(p, estPose1))

		hull_trans2 = []
		for p in hull2:
			hull_trans2.append(gen_icp.dispPoint(p, estPose2))

		#offset2 = [transform[0,0],transform[1,0],transform[2,0]]
				
		radius, center = gen_icp.computeEnclosingCircle(hull1)
		circle1 = [radius, center]
		
		offset = gen_icp.gen_ICP2(estPose1, firstGuess, hull1, hull2, [circle1], costThresh, minMatchDist, plotIteration, n1, n2)
		
		return offset
		
		"""
		pylab.clf()
		for i in range(7):
			path_set = paths[i]

			xP = [0.0]
			yP = [0.0]
			
			ind = range(i+1,7)

			for j in ind:
				loc = path_set[j][0]
				xP.append(loc[0,0])
				yP.append(loc[1,0])

			ind = range(0,i)
			ind.reverse()
			
			for j in ind:
				loc = path_set[j][0]
				xP.insert(0, loc[0,0])
				yP.insert(0, loc[1,0])
			
			pylab.plot(xP,yP)
		pylab.show()			
		"""
		
		
		


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
		return
		
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

		print "saving global occupancy map"
		
		self.occMap.saveMap()

		print "saving global boundary map"
		
		self.boundMap.saveMap()
		
		print "saving global obstacle map"

		self.obstMapImage.save("mapObstGraph%04u.png" % self.saveCount)	

		print "saving global frontier map"

		self.frontierMap.saveMap()			
		
		print "saving global voronoi map"

		self.voronoiMap.saveMap()			

		print "building global navigation road map"

		self.navRoadMap = NavRoadMap(self.mapSize, self.probe, self.voronoiMap.getGraph(), localNode = self.currNode)

		print "done building global maps"
			
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

		
		
		
				
