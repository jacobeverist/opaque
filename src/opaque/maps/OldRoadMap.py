import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *


import Image
import ImageDraw
import graph
import pylab
from copy import *
from math import *
from voronoi import *

class RoadMap:

	# take in vertices and edges of voronoi diagram

	def __init__(self, shadowMap):
		self.imageCount = 0
		self.size = shadowMap.size
		self.shadowTime = Image.new("L", self.size, 0)
		self.setImage(shadowMap)
	

	# select nearest frontier point to the robot's position
	def selectExplorationPoint(self, point):

		imPix = self.boundaryImage.load()
		#self.imageCount += 1
		#self.boundaryImage.save("boundaryImage%04u.png")
		#self.boundaryImage.save("boundaryImage%04u.png" % self.imageCount)

		# convert point to grid coordinates
		indexX, indexY = self.realToGrid(point)

		if indexX >= self.numPixel or indexX < 0:
			print "point X out of range of map", point[0], indexX
			raise

		if indexY >= self.numPixel or indexY < 0:
			print "point Y out of range of map", point[1], indexY
			raise

		#print "start: ", indexX, indexY

		searchComplete = False
		targetFound = False
		targetIndex = (0,0)

		dist = 1e100 
		minIndex = [0,0]
		for i in range(self.numPixel):
			for j in range(self.numPixel):
				if imPix[i,j] == 255:
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
		imPix = self.boundaryImage.load()
		#draw = ImageDraw.Draw(self.boundaryImage)

		indX, indY = self.realToGrid(point)
		
		for i in range(indX-1, indX+2):
			for j in range(indY-1, indY+2):
				imPix[i,j] = 127


	def getClosestEdgePoint(self, edge, point):

		# two possibilities for closest point:
		# A. either of two end points
		# B. somewhere in the middle of the edge

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
		# given point in world, find the closest point on voronoi diagram
		# edge, and point on edge

		minDist = 1e100
		minPoint = [0.0,0.0]
		minEdge = 0

		# find the closest edge to the point
		edges = self.gr.edges()
		for edge in edges:
			v1 = self.gr.get_node_attributes(edge[0])
			v2 = self.gr.get_node_attributes(edge[1])

			dist, cPoint = self.getClosestEdgePoint([v1,v2], point)

			if dist < minDist:
				minDist = dist
				minPoint = cPoint
				minEdge = edge

		return minDist, minPoint, minEdge

	def computeHeadPath(self, start, goal, root):

		# compute path, backing from start, then going forward again towards goal
		# root is the base of the roadmap tree

		# the path from the starting point to the origin
		#print "computing originPath"
		originPath = self.computePath(root, start)
		#print originPath

		# the path from the origin to the goal
		#print "computing goalPath"
		goalPath = self.computePath(root, goal)
		#print goalPath

		# find the break point between the two paths
		#print "origin goal"
		breakPoint = copy(originPath[0])

		if len(originPath) < len(goalPath):
			minLen = len(originPath)
		else:
			minLen = len(goalPath)


		for i in range(minLen):
			#print originPath[i], goalPath[i]

			if originPath[i] != goalPath[i]:
				break

			breakPoint = copy(originPath[i])

		originPath.reverse()

		#print "originPath:", originPath
		#print "goalPath:", goalPath

		return originPath, goalPath, breakPoint
	
	def computePath(self, start, goal):

		# find closest points on voronoi diagram,
		startDist, startPoint, startEdge = self.getClosestPoint(copy(start))
		goalDist, goalPoint, goalEdge = self.getClosestPoint(copy(goal))

		# pick a vertex on both the start and goal edges
		vStart = startEdge[0]
		vGoal = goalEdge[0]

		# then compute path between two points on tree
		shortSpanTree, shortDist = self.gr.shortest_path(vStart)

		# path to be built
		path = []

		vx = vGoal
		path.append(vx)

		#print "startPoint:", startPoint
		#print "goalPoint:", goalPoint
		#print "start:", start
		#print "goal:", goal

		while True:
			#print path

			try:
				vx = shortSpanTree[vx]	
			except:
				print "Error: path not found! Disconnected graph?"
				raise
				# what to do if path not found?
				# FIXME: find a partial path.....

			if vx == None:
				#path.append(vStart)
				break

			path.append(copy(vx))

		# put it into robot following order
		path.reverse()

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
			vReal = self.gr.get_node_attributes(v)
			realPath.append(copy(vReal))

		# add the closest points to the realPath
		realPath.insert(0,startPoint)
		realPath.append(goalPoint)

		# add edges from start and goal to the voronoi diagram
		realPath.insert(0,copy(start))
		realPath.append(copy(goal))


		#print "pre-cleaned realPath:", realPath

		# remove duplicates
		cleanPath = []
		for i in range(len(realPath)-1):
			dist = sqrt((realPath[i][0]-realPath[i+1][0])**2 + (realPath[i][0]-realPath[i+1][0])**2)
			if dist > 1e-4:
				cleanPath.append(realPath[i])

		return cleanPath

		# 1. create an edge from start point to startEdge
		# 2. create an edge from goalEdge to goal point
		# 3. pick closest vertex on startEdge and closest vertex on goalEdge
		# 4. find optimal path from start vertex to goal vertex
		# 5. concatenate edges to plan of 1) and 2) plus an edge from points to vertices

	

		
		
		
		