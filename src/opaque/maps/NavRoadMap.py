
from Map import Map
from copy import copy
from math import *
from Pose import Pose

nodeCount = 0

class NavRoadMap(Map):

	# take in vertices and edges of roadmap graph
	def __init__(self, mapSize, probe, roadGraph, localNode = 0):
		global nodeCount
		Map.__init__(self, mapSize)
		
		self.probe = probe
		self.fileName = "mapNavRoadMap%04u.png"
		
		self.update(roadGraph)
		
		self.cleanPath = []
		
		self.currNode = localNode
		
		# stuff for drawing purposes
		self.nodeCount = nodeCount
		self.childNodes = []
		self.childEntities = []

		self.parentNode = 0
		
		self.originPose = [0.0,0.0,0.0]
		nodeCount += 1
		
	def __del__(self):
		pass
		#self.clearDraw()

	def clearDraw(self):

		if self.parentNode == 0:
			self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("navRoadCurveRoot" + str(self.nodeCount))
		
		# remove all children
		self.parentNode.removeAllChildren()

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []
					
	def draw(self):
		
		#self.clearDraw()

		if self.currNode != 0:
			pose = Pose(self.currNode.getEstPose())
			pose2 = Pose( self.probe.getActualJointPose(self.currNode.rootNode))
		
		"""
		# draw the points in the simulation
		for i in range(len(self.cleanPath)):
			
			pnt = copy(self.cleanPath[i])
			
			" convert points from estimated points to actual points in view "
			if self.currNode != 0:
				localPnt = pose.convertGlobalToLocal(pnt)
				pnt = pose2.convertLocalToGlobal(localPnt)
			
			childNode = self.parentNode.createChildSceneNode("navRoadCurvePoint" + str(self.nodeCount) + "_" + str(i))
			self.childNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("navRoadCurveEnt" + str(self.nodeCount) + "_" +str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)
	
			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			childNode.setPosition(position)
	
			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)
	
			childNode.attachObject(currEntity)
		"""
		
	def update(self, roadGraph):
		self.roadGraph = roadGraph

	def computeHeadPath(self, start, goal, root):
		'''
		Compute path, backing from start, then going forward again towards goal
		Root is the base of the roadmap graph
		'''

		print "computeHeadPath(", start, goal, root, ")"

		self.originPose = copy(start)

		# the path from the starting point to the origin
		originPath = self.computePath(root, start)

		# the path from the origin to the goal
		goalPath = self.computePath(root, goal)

		# find the break point between the two paths
		breakPoint = copy(originPath[0])

		if len(originPath) < len(goalPath):
			minLen = len(originPath)
		else:
			minLen = len(goalPath)

		for i in range(minLen):
			if originPath[i] != goalPath[i]:
				break

			breakPoint = copy(originPath[i])

		originPath.reverse()

		return originPath, goalPath, breakPoint
	
	def computePath(self, start, goal):
		'''
		Create a path of nodes from a start position to a goal position.
		Adds an edge to the roadmap from the start and goal points.
		
		1. create an edge from start point to startEdge
		2. create an edge from goalEdge to goal point
		3. pick closest vertex on startEdge and closest vertex on goalEdge
		4. find optimal path from start vertex to goal vertex
		5. concatenate edges to plan of 1) and 2) plus an edge from points to vertices
		'''

		
		print "computing path from", start, "to", goal

		# find closest points on voronoi diagram,
		startDist, startPoint, startEdge = self.getClosestPoint(copy(start))
		goalDist, goalPoint, goalEdge = self.getClosestPoint(copy(goal))
		
		print "closest points are:", startPoint, startDist, startEdge
		print "and:", goalPoint, goalDist, goalEdge

		# pick a vertex on both the start and goal edges
		vStart = startEdge[0]
		vGoal = goalEdge[0]
		
		#print vStart, vGoal, startEdge, goalEdge

		# then compute path between two points on tree
		shortSpanTree, shortDist = self.roadGraph.shortest_path(vStart)

		# path to be built
		path = []

		vx = vGoal
		#print vx
		path.append(vx)

		while True:
			try:
				vx = shortSpanTree[vx]	
			except:
				print shortSpanTree
				print "Error: path not found! Disconnected graph?"
				raise
				# what to do if path not found?
				# FIXME: find a partial path.....

			if vx == None:
				break

			#print vx
			path.append(copy(vx))

		# put it into robot following order
		path.reverse()

		#print "path in NavRoadMap"
		#print path

		# check for degenerate paths first
		if len(path) >= 2:

			# if both vertices of an edge are on the path,
			# discard the first one if startEdge,

			if (path[0] == startEdge[0] and path[1] == startEdge[1]) or (path[0] == startEdge[1] and path[1] == startEdge[0]):
				path = path[1:]

			# discard the last one if goalEdge
			if (path[:-1] == goalEdge[0] and path[:-2] == goalEdge[1]) or (path[:-1] == goalEdge[1] and path[:-2] == goalEdge[0]):
				path = path[:-1]

		# convert the path to real coordinates
		realPath = []
		for v in path:
			vReal = self.roadGraph.get_node_attributes(v)
			realPath.append(copy(vReal))

		# add the closest points to the realPath
		realPath.insert(0,startPoint)
		realPath.append(goalPoint)

		# add edges from start and goal to the voronoi diagram
		realPath.insert(0,copy(start))
		realPath.append(copy(goal))

		# remove duplicates
		cleanPath = []
		for i in range(len(realPath)-1):
			dist = sqrt((realPath[i][0]-realPath[i+1][0])**2 + (realPath[i][0]-realPath[i+1][0])**2)
			if dist > 1e-4:
				cleanPath.append(realPath[i])
				
		self.cleanPath = cleanPath

		return cleanPath

	

	def getClosestEdgePoint(self, edge, point):
		'''
		Given point in the world, find the closest point on the edge
		
		- two possibilities for closest point:
		A. either of two end points
		B. somewhere in the interior of the edge
		'''

		#print "getClosestEdgePoint(", edge, point, ")"

		v1 = edge[0]
		v2 = edge[1]

		if v1 == v2:
			dist = sqrt((v1[0]-point[0])**2 + (v1[1]-point[1])**2)
			return dist, v1

		edgeVec1 = [v2[0]-v1[0],v2[1]-v1[1]]
		edgeMag1 = sqrt(edgeVec1[0]**2 + edgeVec1[1]**2)
		edgeVec1 = [edgeVec1[0]/edgeMag1, edgeVec1[1]/edgeMag1]

		edgeVec2 = [v1[0]-v2[0],v1[1]-v2[1]]
		edgeMag2 = sqrt(edgeVec2[0]**2 + edgeVec2[1]**2)
		edgeVec2 = [edgeVec2[0]/edgeMag2, edgeVec2[1]/edgeMag2]

		pointVec1 = [point[0]-v1[0], point[1]-v1[1]]
		pointMag1 = sqrt(pointVec1[0]**2 + pointVec1[1]**2)
		pointVec1 = [pointVec1[0]/pointMag1, pointVec1[1]/pointMag1]

		pointVec2 = [point[0]-v2[0], point[1]-v2[1]]
		pointMag2 = sqrt(pointVec2[0]**2 + pointVec2[1]**2)
		pointVec2 = [pointVec2[0]/pointMag2, pointVec2[1]/pointMag2]

		# compute the angle from dot product formula
		dotProduct1 = edgeVec1[0]*pointVec1[0] + edgeVec1[1]*pointVec1[1]
		dotProduct2 = edgeVec2[0]*pointVec2[0] + edgeVec2[1]*pointVec2[1]

		# correct for any numerical errors
		if dotProduct1 > 1.0:
			dotProduct1 = 1.0
		if dotProduct1 < -1.0:
			dotProduct1 = -1.0
		if dotProduct2 > 1.0:
			dotProduct2 = 1.0
		if dotProduct2 < -1.0:
			dotProduct2 = -1.0

		cornerAngle1 = acos(dotProduct1)
		cornerAngle2 = acos(dotProduct2)

		# if closer to somewhere in the middle,
		# compute the orthogonal distance and closest point
		if cornerAngle1 < pi/2 and cornerAngle2 < pi/2:
			dist = pointMag1*sin(cornerAngle1)
			tanDist = pointMag1*cos(cornerAngle1)
			contactPoint = [tanDist*edgeVec1[0] + v1[0], tanDist*edgeVec1[1] + v1[1]]
			return dist, contactPoint

		if pointMag1 < pointMag2:
			return pointMag1, v1
		else:
			return pointMag2, v2

	def getClosestPoint(self, point):
		''' 
		given point in world, find the closest point on roadmap
		returns edge, point on edge, and minimum distance to the roadmap
		'''

		minDist = 1e100
		minPoint = [0.0,0.0]
		minEdge = 0

		# find the closest edge to the point
		edges = self.roadGraph.edges()
		for edge in edges:
			v1 = self.roadGraph.get_node_attributes(edge[0])
			v2 = self.roadGraph.get_node_attributes(edge[1])

			dist, cPoint = self.getClosestEdgePoint([v1,v2], point)

			if dist < minDist:
				minDist = dist
				minPoint = cPoint
				minEdge = edge

		return minDist, minPoint, minEdge

		
		
		
		
