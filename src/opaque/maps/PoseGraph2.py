import sys

from numpy import *
from scipy.optimize import *
import scipy
import scipy.linalg
import numpy
import numpy.linalg
import random
import csv

from Map import Map
from LocalNode import LocalNode, getLongestPath
from StablePose import StablePose
from StableCurve import StableCurve
from GPACurve import GPACurve
from SplineFit import SplineFit
from Pose import Pose
import gen_icp
from functions import *
import toro
import scsgp
import bayes
import estimateMotion
from operator import itemgetter

import pylab
from matplotlib.patches import Circle

import socket

from math import *
from copy import copy

import Image
import ImageDraw

import cProfile

from subprocess import Popen, PIPE
from medialaxis import computeMedialAxis
import graph

renderCount = 0
topCount = 0

GND_PRIORITY = 10
PATH_PRIORITY = 6
INPLACE_PRIORITY = 5
CORNER_PRIORITY = 3
SHAPE_PRIORITY = 2
OVERLAP_PRIORITY = 1
MOTION_PRIORITY = 0
INPLACE_PRIORITY2 = -1
INPLACE_PRIORITY3 = -2


def decimatePoints(points):
	result = []

	for i in range(len(points)):
		#if i%2 == 0:
		if i%4 == 0:
			result.append(points[i])

	return result

def computeBareHull(node1, sweep = False, static = False):
	
	if static:
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		a_data = decimatePoints(a_data)

		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
	else:
				
		" Read in data of Alpha-Shapes without their associated covariances "
		node1.computeAlphaBoundary(sweep = sweep)
		a_data = node1.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
	
	return a_data_GPAC

def computeHull(node1, sweep = False, static = False):
	
	if static:
		
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		a_data = decimatePoints(a_data)

		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		" treat the points with the point-to-line constraint "
		gen_icp.addPointToLineCovariance(a_data_GPAC, high_var=1.0, low_var=0.001)
	
		" treat the points with the distance-from-origin increasing error constraint "
		gen_icp.addDistanceFromOriginCovariance(a_data_GPAC, tan_var=0.1, perp_var=0.01)

	else:			
		" Read in data of Alpha-Shapes and add their associated covariances "
		node1.computeAlphaBoundary(sweep = sweep)
		a_data = node1.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		" treat the points with the point-to-line constraint "
		gen_icp.addPointToLineCovariance(a_data_GPAC, high_var=1.0, low_var=0.001)
	
		" treat the points with the distance-from-origin increasing error constraint "
		gen_icp.addDistanceFromOriginCovariance(a_data_GPAC, tan_var=0.1, perp_var=0.01)
	
	return a_data_GPAC

def computeBoundary(node1, sweep = False):
	
	" Read in data of Alpha-Shapes and add their associated covariances "
	node1.computeAlphaBoundary(sweep)
	#node1.boundaryMap.getBoundaryPoints()
	#computeAlphaBoundary(sweep = sweep)
	a_data = node1.getAlphaBoundary(sweep = sweep)
	a_data = decimatePoints(a_data)
	
	" convert hull points to GPAC coordinates before adding covariances "
	localGPACPose = node1.getLocalGPACPose()
	localGPACProfile = Pose(localGPACPose)
	
	a_data_GPAC = []
	for pnt in a_data:
		a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
	
	" treat the points with the point-to-line constraint "
	gen_icp.addPointToLineCovariance(a_data_GPAC, high_var=1.0, low_var=0.001)

	" treat the points with the distance-from-origin increasing error constraint "
	gen_icp.addDistanceFromOriginCovariance(a_data_GPAC, tan_var=0.1, perp_var=0.01)
	
	return a_data_GPAC

def computeHullAxis(nodeID, node2, tailCutOff = False):
	
	if node2.isBowtie:			
		hull2 = computeBareHull(node2, sweep = False, static = True)
		hull2.append(hull2[0])
		medial2 = node2.getStaticMedialAxis()

	else:
		hull2 = computeBareHull(node2, sweep = False)
		hull2.append(hull2[0])
		medial2 = node2.getMedialAxis(sweep = False)

	" take the long length segments at tips of medial axis"
	edge1 = medial2[0:2]
	edge2 = medial2[-2:]
	
	frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
	backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
	
	frontVec[0] /= frontMag
	frontVec[1] /= frontMag
	backVec[0] /= backMag
	backVec[1] /= backMag
	
	" make a smaller version of these edges "
	newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
	newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

	edge1 = [newP1, edge1[1]]
	edge2 = [edge2[0], newP2]

	
	" find the intersection points with the hull "
	interPoints = []
	for k in range(len(hull2)-1):
		hullEdge = [hull2[k],hull2[k+1]]
		isIntersect1, point1 = Intersect(edge1, hullEdge)
		if isIntersect1:
			interPoints.append(point1)
			break

	for k in range(len(hull2)-1):
		hullEdge = [hull2[k],hull2[k+1]]
		isIntersect2, point2 = Intersect(edge2, hullEdge)
		if isIntersect2:
			interPoints.append(point2)
			break
	
	" replace the extended edges with a termination point at the hull edge "			
	medial2 = medial2[1:-2]
	if isIntersect1:
		medial2.insert(0, point1)
	if isIntersect2:
		medial2.append(point2)
	

	" cut off the tail of the non-sweeping side "
	TAILDIST = 0.5

	if tailCutOff:
		
		if nodeID % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]
	
		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
			
	return hull2, medial2


def doTransform(T1, T2, E1, E2):

	x1 = T1[0,0]	
	y1 = T1[1,0]
	p1 = T1[2,0]

	x2 = T2[0,0]	
	y2 = T2[1,0]
	p2 = T2[2,0]

	newT = matrix([[x1 + x2*cos(p1) - y2*sin(p1)],
		[y1 + x2*sin(p1) + y2*cos(p1)],
		[p1+p2]],dtype=float)

	J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]],dtype=float)
	J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]],dtype=float)
						
	newE = J1 * E1 * J1.T + J2 * E2 * J2.T
	
	return newT, newE


class PoseGraph:

	def __init__(self, probe, contacts):
		
		self.probe = probe
		self.contacts = contacts
		
		self.initPose = self.probe.getActualJointPose(19)
		self.nodeHash = {}
		self.edgeHash = {}
		
		self.numNodes = 0
		self.currNode = 0

		self.inplace1 = []
		self.inplace2 = []
		self.inplace3 = []

		self.corner_constraints = []
		self.backtrack_constraints = []
		self.overlap_constraints = []
		self.motion_constraints = []
		self.sensor_constraints = []
		self.inplace_constraints = []
		self.gnd_constraints = []		
		self.merged_constraints = []
		self.sensorHypotheses = []
		self.activeHypotheses = []
		self.cornerHypotheses = []
		self.cornerHash = {}
		self.badHypotheses = []
		self.badHypotheses1 = []
		self.badHypotheses2 = []

		self.goodHypotheses = []
		self.allCornerHypotheses = []
		
		self.sweepHulls = []
		self.nonSweepHulls = []
		
		self.edgePriorityHash = {}
		self.cornerBins = []
		self.c_hulls = []
		self.b_hulls = []
		self.a_hulls = []
		self.a_medial = []

		self.E_gnd = matrix([[ 0.00001, 0.0, 0.0 ],
							[ 0.0, 0.00001, 0.0],
							[ 0.0, 0.0, 0.0001 ]])

		self.E_junction = matrix([[ 0.0001, 0.0, 0.0 ],
							[ 0.0, 0.0001, 0.0],
							[ 0.0, 0.0, 0.02 ]])

		self.E_corner = matrix([[ 0.01, 0.0, 0.0 ],
							[ 0.0, 0.01, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

		self.E_inplace = matrix([[ 0.05, 0.0, 0.0 ],
							[ 0.0, 0.05, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

		self.E_overlap = matrix([[ 0.1,  0.0, 0.0],
							[ 0.0,  0.01, 0.0],
							[0.0, 0.0,  0.02]])

		self.E_featureless = matrix([[ 1.0,  0.0, 0.0],
							[ 0.0,  0.01, 0.0],
							[0.0, 0.0,  0.02]])
		#elf.E_featureless = matrix([[ 0.1,  0.0, 0.0],
		#					[ 0.0,  0.0001, 0.0],
		#					[0.0, 0.0,  0.02]])

		#self.E_overlap = matrix([[ 0.2,  0.0, 0.0],
		#					[ 0.0,  0.02, 0.0],
		#					[0.0, 0.0,  0.1]])

		self.E_motion = matrix([[ 0.06, 0.0, 0.0],
							[0.0,  5.0, 0.0],
							[0.0, 0.0,  0.1 ]])
		
		self.E_sensor = matrix([[0.1,0.0,0.0],
							[0.0,0.05,0.0],
							[0.0,0.0,0.02]])
		
		
		self.pathPlotCount = 0
		self.overlapPlotCount = 0
		self.pathDrawCount = 0
		#self.numPaths = 2
		#self.paths = {0 : [], 1 : []}
		#self.hulls = {0 : [], 1 : []}
		#self.nodeSets = {0 : [], 1 : []}
		self.numPaths = 1
		self.paths = {0 : []}
		self.hulls = {0 : []}
		self.nodeSets = {0 : []}
		self.currPath = 0
		self.junctionPoint = 0
		self.pathParents = [[None,None,None,None]]

		self.trimCount = 0
		self.spliceCount = 0
		self.orderCount = 0		
		self.stateCount = 0
		
		self.colors = []
		#for i in range(1000):
		#	self.colors.append((random.random(),random.random(),random.random()))

		self.colors = ['b','r','g','k','y']

	def saveState(self):
		
		saveFile = ""
		
		saveFile += "self.initPose = " + repr(self.initPose) + "\n"
		saveFile += "self.numNodes = " + repr(self.numNodes) + "\n"

		saveFile += "self.cornerBins = " + repr(self.cornerBins) + "\n"

		#self.cornerHash = {}
		"TODO: save localNode corner candidates "
		
		saveFile += "self.edgePriorityHash = " + repr(self.edgePriorityHash) + "\n"
		#saveFile += "self.edgeHash = " + repr(self.edgeHash) + "\n"
		
		#self.edgePriorityHash = {}
		
		saveFile += "self.numPaths = " + repr(self.numPaths) + "\n"
		saveFile += "self.nodeSets = " + repr(self.nodeSets) + "\n"
		saveFile += "self.junctionPoint = " + repr(self.junctionPoint) + "\n"

		saveFile += "self.currPath = " + repr(self.currPath) + "\n"

		saveFile += "self.pathParents = " + repr(self.pathParents) + "\n"
		

		f = open("stateSave_%04u.txt" % (self.numNodes-1), 'w')
		f.write(saveFile)
		f.close()		
		
		
		for nodeID in range(self.numNodes):
			
			f = open("cornerCandidates_%04u.txt" % (nodeID), 'w')
			f.write(repr(self.nodeHash[nodeID].cornerCandidates))
			f.close()
		
		
		
		
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/stateSave_%04u.txt" % (numNodes-1)
		f = open(dirName + "/stateSave_%04u.txt" % (numNodes-1), 'r')		
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		exec(saveStr)
		
		print self.numNodes
		print self.edgePriorityHash

		
	def splicePathIDs(self, pathIDs):
		
		if len(pathIDs) == 0:
			return []
		
		if len(pathIDs) == 1:
			return [self.trimmedPaths[pathIDs[0]]]

		" Assumption:  the pathIDs are connected paths with parent-child relations "


		" find the root "
		" pick any node, and go up the tree until we hit the root "
		currPathID = pathIDs[0]
		
		while True:
			parent = self.pathParents[currPathID]
			thisParentID = parent[0]
			
			if thisParentID == None:
				break
			
			if pathIDs.count(thisParentID) == 0:
				break
			
			currPathID = thisParentID

		" currPathID is the root path "
		rootPathID = currPathID
		
		" for each path, attempt to join with its parent path "
		joins = []
		for pathID in pathIDs:

			if self.pathParents[pathID] == None:
				continue
						
			parentPathID = self.pathParents[pathID][0]
			
			" parent does not concern us "
			if pathIDs.count(parentPathID) < 1:
				continue
			
			junctionNodeID = self.pathParents[pathID][1]
			localJunctionPoint = self.pathParents[pathID][2]
			poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
			junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)

			path1 = self.trimmedPaths[pathID]
			path2 = self.trimmedPaths[parentPathID]

			
			#print "junctionPoint =", junctionPoint
			#print "path1:", self.paths[0]
			#print "path2:", self.paths[1]
	
			minDist1 = 1e100
			minI1 = 0		
			for i in range(len(path1)):
				pnt = path1[i]
				dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
			
				if dist < minDist1:
					minDist1 = dist
					minI1 = i
	
			minDist2 = 1e100
			minI2 = 0		
			for i in range(len(path2)):
				pnt = path2[i]
				dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
			
				if dist < minDist2:
					minDist2 = dist
					minI2 = i
			
			joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])

		print "joins:", joins

		" create a tree with the node IDs and then stitch them together with joins "
		pathGraph = graph.graph()
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for k in range(len(path)):
				pathGraph.add_node((pathID, k), path[k])
			for k in range(len(path)-1):
				pathGraph.add_edge((pathID, k), (pathID, k+1))

		" join with the junction in between the join points "
		for k in range(len(joins)):
			join = joins[k]
			pathGraph.add_node(k, join[2])
			pathGraph.add_edge(join[0], k)
			pathGraph.add_edge(k, join[1])
			
		" if both terminal paths are child paths, 1 resultant spliced path "
		" if one terminal path is the root, then 2 resultant spliced paths "
		rightPathID = pathIDs[0]
		leftPathID = pathIDs[-1]

		newPaths = []		
		if rightPathID == rootPathID or leftPathID == rootPathID:
			
			print "one of paths is root make 2 resultant spliced paths"
			
			if rightPathID != rootPathID:
				childPathID = rightPathID
			else:
				childPathID = leftPathID
			
			" find path 1 from rootPathID to childPathID "
			" find path 2 from rootPathID to childPathID "
			rootPath = self.trimmedPaths[rootPathID]
			childPath = self.trimmedPaths[childPathID]


			for join in joins:
				if join[0][0] == childPathID:
					childJunctionPoint = join[2]
					childI = join[0][1]

			if childI < abs(childI-len(childPath)-1):
				childTermI = len(childPath)-1
			else:
				childTermI = 0

			"find path between: (childPathID, childTermI) to (rootPathID, 0) and (rootPathID, len(rootPath)-1)"
			shortestPathSpanTree, shortestDist = pathGraph.shortest_path((childPathID, childTermI))
			startNode1 = shortestPathSpanTree[(rootPathID, 0)]
			startNode2 = shortestPathSpanTree[(rootPathID, len(rootPath)-1)]
			
			print "len(shortestPathSpanTree) =", len(shortestPathSpanTree)
			#print "shorestPathSpanTree:", shortestPathSpanTree

			currNode = startNode1
			splicedPath1 = []
			while currNode != (childPathID, childTermI):
				splicedPath1.append(pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath1.append(pathGraph.get_node_attributes(currNode))

			currNode = startNode2
			splicedPath2 = []
			while currNode != (childPathID, childTermI):
				splicedPath2.append(pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath2.append(pathGraph.get_node_attributes(currNode))

			#print "splicedPath1:", splicedPath1
			#print "splicedPath2:", splicedPath2
			
			newPaths.append(splicedPath1)
			newPaths.append(splicedPath2)
			
		else:
			"find entire path from rightPathID to leftPathID "
			" start from point farthest from child's junction point "

			print "both paths are childs,  make 1 resultant spliced paths"

			rightPath = self.trimmedPaths[rightPathID]
			leftPath = self.trimmedPaths[leftPathID]
		
			" find the junction of this child's to its parent "	
			#joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])
			
			for join in joins:
				if join[0][0] == rightPathID:
					rightJunctionPoint = join[2]
					rightI = join[0][1]

				if join[0][0] == leftPathID:
					leftJunctionPoint = join[2]
					leftI = join[0][1]
					
			
			if leftI < abs(leftI-len(leftPath)-1):
				leftTermI = len(leftPath)-1
			else:
				leftTermI = 0

			if rightI < abs(rightI-len(rightPath)-1):
				rightTermI = len(rightPath)-1
			else:
				rightTermI = 0
				
			"find path between: (rightPathID, rightTermI) to (leftPathID, rightTermI) "
			shortestPathSpanTree, shortestDist = pathGraph.shortest_path((rightPathID, rightTermI))
			currNode = shortestPathSpanTree[(leftPathID, leftTermI)]

			splicedPath = []
			while currNode != (rightPathID, rightTermI):
				splicedPath.append(pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath.append(pathGraph.get_node_attributes(currNode))

			print "len(shortestPathSpanTree) =", len(shortestPathSpanTree)
			#print "splicedPath:", splicedPath
			#print "shorestPathSpanTree:", shortestPathSpanTree

			newPaths.append(splicedPath)


	
			pylab.clf()
	
			print "spliced path has", len(splicedPath), "points"
			xP = []
			yP = []
			for p in splicedPath:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP)
	
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
	
			pylab.title("Spliced Paths 1, numNodes = %d" % self.numNodes)
			pylab.savefig("splicedPath_%04u.png" % self.spliceCount)
			self.spliceCount += 1



		return newPaths
			


	def getTopology(self, nodes = 0):
		global topCount
		
		def convertAlphaUniform(a_vert, max_spacing = 0.04):
			
			" make the vertices uniformly distributed "
			
			new_vert = []
		
			for i in range(len(a_vert)):
				p0 = a_vert[i]
				p1 = a_vert[(i+1) % len(a_vert)]
				dist = math.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
	
				vec = [p1[0]-p0[0], p1[1]-p0[1]]
				vec[0] /= dist
				vec[1] /= dist
				
				new_vert.append(copy(p0))
				
				if dist > max_spacing:
					" cut into pieces max_spacing length or less "
					numCount = int(math.floor(dist / max_spacing))
					
					for j in range(1, numCount+1):
						newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
						new_vert.append(newP)
			
			return new_vert			

		medialPointSoup = []

		if nodes == 0:
			nodes = range(self.numNodes)
			
		if nodes == []:
			return [], []


		#pylab.clf()	

		for nodeID in nodes:
			estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()		
	
			if self.nodeHash[nodeID].isBowtie:			
				hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
			else:
				hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
	
			" set the origin of pose 1 "
			poseOrigin = Pose(estPose1)
	
			#xP = []
			#yP = []	
			for p in hull1:
				p1 = poseOrigin.convertLocalToGlobal(p)
				medialPointSoup.append(p1)
				#xP.append(p1[0])
				#yP.append(p1[1])

			#pylab.scatter(xP,yP)
		"""
		#pylab.xlim(-4.4)
		#pylab.ylim(-3,3)
		pylab.xlim(-4,4)
		pylab.ylim(-4,4)
		pylab.savefig("plotMedialSoup_%04u.png" % (topCount))
		topCount += 1
		"""
		


		nodeID = self.numNodes - 1
			
		radius = 0.2

		numPoints = len(medialPointSoup)
		inputStr = str(numPoints) + " "

		" alpha shape circle radius "
		inputStr += str(radius) + " "

		#for p in medialPointSoup:
		#	p2 = copy(p)
		#	" add a little bit of noise to avoid degenerate conditions in CGAL "
		#	inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		#
		#inputStr += "\n"


		isDone = False
		
		while not isDone:

			inputStr = str(numPoints) + " "
	
			" alpha shape circle radius "
			inputStr += str(radius) + " "
			
			for p in medialPointSoup:
				p2 = copy(p)
				" add a little bit of noise to avoid degenerate conditions in CGAL "
				p2[0] += random.gauss(0.0,0.000001)
				p2[1] += random.gauss(0.0,0.000001)
	
				inputStr += str(p2[0]) + " " + str(p2[1]) + " "
			
			inputStr += "\n"
			
			try:			
				" start the subprocess "
				if sys.platform == "win32":
					subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
				else:
					subProc = Popen(["./alpha2"], stdin=PIPE, stdout=PIPE)
					
				
				" send input and receive output "
				sout, serr = subProc.communicate(inputStr)
		
				" convert string output to typed data "
				sArr = sout.split(" ")
				
				numVert = int(sArr[0])
				
				sArr = sArr[1:]
				
				
				vertices = []
				for i in range(len(sArr)/2):
					vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
				isDone = True
			except:
				print "hull has holes!  retrying..."
				#print sArr	

		"""
		vertices = []
		try:			
			" start the subprocess "
			#subProc = Popen(["./alpha2.exe"], stdin=PIPE, stdout=PIPE)
			if sys.platform == "win32":
				subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
			else:
				subProc = Popen(["./alpha2"], stdin=PIPE, stdout=PIPE)
			
			#print subProc.stdin, subProc.stderr, subProc.stdout
			#print "input:", inputStr
			" send input and receive output "
			sout, serr = subProc.communicate(inputStr)

			" convert string output to typed data "
			sArr = sout.split(" ")
			#print sArr
			
			numVert = int(sArr[0])
			
			sArr = sArr[1:]
			
			vertices = []
			for i in range(len(sArr)/2):
				vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
		except:
			print "hull has holes!"
		"""
		
		" cut out the repeat vertex "
		vertices = vertices[:-1]
		
		vertices = convertAlphaUniform(vertices)

		vertices.append(vertices[0])
		
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in vertices:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]
	
	
		PIXELSIZE = 0.05
		mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
		pixelSize = PIXELSIZE
		numPixel = int(2.0*mapSize / pixelSize + 1.0)
		divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)
		
		def realToGrid(point):
			indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
			indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
			return indexX, indexY
	
		def gridToReal(indices):
			i = indices[0]
			j = indices[1]
			point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
			return point
	
		gridHull = []
		for i in range(len(vertices)):
			p = vertices[i]
			gridHull.append(realToGrid(p))
	
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in gridHull:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]

		resultImg = Image.new('L', (numPixel,numPixel))
		resultImg = computeMedialAxis(nodeID, numPixel,numPixel, resultImg, len(gridHull[:-2]), gridHull[:-2])

		imgA = resultImg.load()
	
		points = []
		for i in range(1, numPixel-1):
			for j in range(1, numPixel-1):
				if imgA[i,j] == 0:
					points.append((i,j))
	
		medialGraph = graph.graph()
		for p in points:
			medialGraph.add_node(p, [])
	
		builtGraph = {}
		for i in range(2, numPixel-2):
			for j in range(2, numPixel-2):
				if imgA[i,j] == 0:
					builtGraph[(i,j)] = []
					for k in range(i-1, i+2):
						for l in range(j-1, j+2):
							if imgA[k,l] == 0:
								builtGraph[(i,j)].append((k,l))
								medialGraph.add_edge((i,j), (k,l))
								
		mst = medialGraph.minimal_spanning_tree()
		
		uni_mst = {}
		isVisited = {}
		nodeSum = {}
		for k, v in mst.items():
			uni_mst[k] = []
			isVisited[k] = 0
			nodeSum[k] = 0
		
		for k, v in mst.items():
			if v != None:
				uni_mst[k].append(v)
				uni_mst[v].append(k)
		
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
	
	
		" find the longest path from each leaf"
		maxPairDist = 0
		maxPair = None
		maxPath = []
		for leaf in leaves:
	
			isVisited = {}
			nodeSum = {}
			nodePath = {}
			for k, v in uni_mst.items():
				isVisited[k] = 0
				nodeSum[k] = 0
				nodePath[k] = []
	
			getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
	
			maxDist = 0
			maxNode = None
			for k, v in nodeSum.items():
				if v > maxDist:
					maxNode = k
					maxDist = v
	
			if maxDist > maxPairDist:
				maxPairDist = maxDist
				maxPair = [leaf, maxNode]
				maxPath = nodePath[maxNode]

		" convert path points to real "
		realPath = []
		for p in maxPath:
			realPath.append(gridToReal(p))
		
		return realPath, vertices

	def checkSupport(self, nodeID1, nodeID2, offset, supportLine):

		node1 = self.nodeHash[nodeID1]
		node2 = self.nodeHash[nodeID2]

		#hull1, medial1 = computeHullAxis(nodeID1, node1, tailCutOff = True)
		hull2, medial2 = computeHullAxis(nodeID2, node2, tailCutOff = True)


		"""
		if node2.isBowtie:			
			hull2 = computeBareHull(node2, sweep = False, static = True)
			hull2.append(hull2[0])
			medial2 = node2.getStaticMedialAxis()

		else:
			hull2 = computeBareHull(node2, sweep = False)
			hull2.append(hull2[0])
			medial2 = node2.getMedialAxis(sweep = False)

		" take the long length segments at tips of medial axis"
		edge1 = medial2[0:2]
		edge2 = medial2[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]

		
		" find the intersection points with the hull "
		interPoints = []
		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect2, point2 = Intersect(edge2, hullEdge)
			if isIntersect2:
				interPoints.append(point2)
				break
		
		" replace the extended edges with a termination point at the hull edge "			
		medial2 = medial2[1:-2]
		if isIntersect1:
			medial2.insert(0, point1)
		if isIntersect2:
			medial2.append(point2)
		

		" cut off the tail of the non-sweeping side "
		TAILDIST = 0.5

		if nodeID2 % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]

		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
		"""
		
		estPose1 = node1.getGlobalGPACPose()		
			
		#minMatchDist2 = 0.5
		minMatchDist2 = 2.0
			
		" set the initial guess "
		poseOrigin = Pose(estPose1)
		
		localSupport = []
		for pnt in supportLine:
			localSupport.append(poseOrigin.convertGlobalToLocal(pnt))
	
		supportSpline = SplineFit(localSupport, smooth=0.1)		
		supportPoints = gen_icp.addGPACVectorCovariance(supportSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		points2 = gen_icp.addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)

	
		" transform pose 2 by initial offset guess "	
		" transform the new pose "
		points2_offset = []
		for p in points2:
			result = gen_icp.dispPoint(p, offset)		
			points2_offset.append(result)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])			
		
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
				
		support_pairs = []
		for i in range(len(poly2)):
			
			p_2 = poly2[i]
	
			" for every transformed point of A, find it's closest neighbor in B "
			p_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

			if gen_icp.isInCircle(p_1, radius2, center2):

				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = p_1[2]

					" we store the untransformed point, but the transformed covariance of the A point "
					support_pairs.append([points2[i],p_1,C2,C1])

		vals = []
		sum1 = 0.0
		for pair in support_pairs:
	
			a = pair[0]
			b = pair[1]
			Ca = pair[2]
			Cb = pair[3]
	
			ax = a[0]
			ay = a[1]		
			bx = b[0]
			by = b[1]
	
			c11 = Ca[0][0]
			c12 = Ca[0][1]
			c21 = Ca[1][0]
			c22 = Ca[1][1]
					
			b11 = Cb[0][0]
			b12 = Cb[0][1]
			b21 = Cb[1][0]
			b22 = Cb[1][1]	
		
			val = gen_icp.computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
			
			vals.append(val)
			sum1 += val
			
		return sum1


	def addPriorityEdge(self, edge, priority):

		try:
			self.edgePriorityHash[(edge[0],edge[1])]
		except:
			self.edgePriorityHash[(edge[0],edge[1])] = []
		
		" corner constraints with priority of 3 "
		self.edgePriorityHash[(edge[0],edge[1])].append([edge[2], edge[3], priority])

	def getEdges(self, nodeID1, nodeID2):

		resultEdges = []
		
		for k, v in self.edgePriorityHash.items():

			if k[0] == nodeID1 and k[1] == nodeID2 or k[0] == nodeID2 and k[1] == nodeID1:
				
				for edge in v:
					resultEdges.append(edge)
				
		return resultEdges

	def getCornerEdges(self):
		
		cornerEdges = {}
		
		for k, v in self.edgePriorityHash.items():

			hasCornerEdge = False			
			for edge in v:
				" if priority status is the corner status "
				if edge[2] == CORNER_PRIORITY:
					hasCornerEdge = True
		
			if hasCornerEdge:
				cornerEdges[k] = []
				
				for edge in v:
					" if priority status is the corner status "
					if edge[2] == CORNER_PRIORITY:
						cornerEdges[k].append(edge)
				
		return cornerEdges
	
	def getPriorityEdges(self, priorityLevel = -1):

		priorityEdges = {}
		
		if priorityLevel == -1:
			for k, v in self.edgePriorityHash.items():
				
				if len(v) > 0:
			
					id1 = k[0]
					id2 = k[1]
					
					" take the highest priority constraints only "
					maxPriority = -1
					for const in v:
						thisPriority = const[2]
						if thisPriority > maxPriority:
							maxPriority = thisPriority
			
					if maxPriority == -1:
						raise
					
					priorityEdges[k] = []
					for const in v:
						if const[2] == maxPriority:
							priorityEdges[k].append(const)

		else:
			
			for k, v in self.edgePriorityHash.items():
				
				if len(v) > 0:
			
					id1 = k[0]
					id2 = k[1]
					
					" take the highest priority constraints only "
					priorityEdges[k] = []
					for const in v:
						if priorityLevel == const[2]:
							
							priorityEdges[k].append(const)


		return priorityEdges

	def deleteAllPriority(self, priorityLevel):
		
		for k, v in self.edgePriorityHash.items():
			newV = []
			
			for const in v:
				if const[2] != priorityLevel:
					newV.append(const)
			self.edgePriorityHash[k] = newV


	def deleteAllEdges(self):

		self.edgePriorityHash = {}

	def computeCornerClusterError(self):
		

		totalError = 0.0
				
		for i in range(len(self.cornerBins)):
			bin = self.cornerBins[i]
			
			" check the internal error of this corner cluster "
			#tup1 = (n1,point1_i)
			#tup2 = (n2,point2_i)
			
			binError = 0.0
			for j in range(len(bin)):
				tup1 = bin[j]
				node1 = self.nodeHash[tup1[0]]
				point1_i = tup1[1]
				corner1 = node1.cornerCandidates[point1_i]
				point1 = corner1[0]
				estPose1 = node1.getGlobalGPACPose()
				poseProf1 = Pose(estPose1)
				
				globalPoint1 = poseProf1.convertLocalToGlobal(point1)
				
				for k in range(j+1, len(bin)):
					tup2 = bin[k]
					node2 = self.nodeHash[tup2[0]]
					point2_i = tup2[1]
					corner2 = node2.cornerCandidates[point2_i]
					point2 = corner2[0]
					estPose2 = node2.getGlobalGPACPose()
					poseProf2 = Pose(estPose2)

					globalPoint2 = poseProf2.convertLocalToGlobal(point2)


					dist = sqrt((globalPoint1[0]-globalPoint2[0])**2 + (globalPoint1[1]-globalPoint2[1])**2)
					
					binError += dist

			print i, binError
			totalError += binError
	
		return totalError
		
	def addCornerBin(self, tup1, tup2):
		
		" put this pair of corners into matched bins "			
		matchBins = []
		for i in range(len(self.cornerBins)):
			
			bin = self.cornerBins[i]
			" first check if either point exists already in a bin"
			if tup1 in bin:
				matchBins.append(i)
			elif tup2 in bin:
				matchBins.append(i)

		" if not shared, create new bin "
		if len(matchBins) == 0:
			self.cornerBins.append([tup1, tup2])
		else:						
			" we have to merge into the old bins "
			if len(matchBins) == 1:
				" only 1 bin, just add the tuples that don't already exist "
				bin_i = matchBins[0]
				
				if not tup1 in self.cornerBins[bin_i]:
					self.cornerBins[bin_i].append(tup1)
				if not tup2 in self.cornerBins[bin_i]:
					self.cornerBins[bin_i].append(tup2)
					
			else:
				" we have more than one matching bin, so now we need to merge "
				newSuperBin = []
				for i in matchBins:
					bin = self.cornerBins[i]
					
					for item in bin:
						if not item in newSuperBin:
							newSuperBin.append(item)

				newCornerBins = [newSuperBin]
				
				for i in range(len(self.cornerBins)):
					if not i in matchBins:
						newCornerBins.append(self.cornerBins[i])

				self.cornerBins = newCornerBins
	
		#print "cornerBins:"
		#for i in range(len(self.cornerBins)):
		#	print i, self.cornerBins[i]


	def addCornerConstraints(self, nodeID):

		global renderCount

		" try to perform corner constraints "
		" overwrite overlap/motion constraints"
		" corner constraints with past nodes "

		" compute djikstra projection "		
		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
	
		" find paired candidates by mahalanobis distance "
		cornerPairs = self.findCornerCandidates(paths, nodeID = nodeID)
	
		print "found", len(cornerPairs), "corner pairs"
	
		" compute hypotheses and reject unsuitable matches "
		hypotheses = []
		for k in range(len(cornerPairs)):
			
			p = cornerPairs[k]
			dist = p[0]
			n1 = p[1]
			n2 = p[2]
			transform = p[3]
			angDist = p[4]
			point1 = p[5]
			point2 = p[6]
			ang1 = p[7]
			ang2 = p[8]
			point1_i = p[9]
			point2_i = p[10]	
	
			hull1 = self.a_hulls[n1]
			hull2 = self.a_hulls[n2]
	
			offset, covar, cost, angDiff, oriDiff = self.makeCornerConstraint(n1, n2, point1, point2, ang1, ang2, hull1, hull2, self.b_hulls[n1], self.b_hulls[n2])
			transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
	
			" ratio of match error threshold"		
			if cost <= 0.6:
				hypotheses.append([n1, n2, transform, covar, cost, point1, point2, angDiff, oriDiff, point1_i, point2_i])
	
		print "kept", len(hypotheses), "corner hypotheses"
		
		" remove multiple constraints between two nodes.  Take best cost "
		cornerHash = {}
		for hyp in hypotheses:
			n1 = hyp[0]
			n2 = hyp[1]
			newCost = hyp[4]
			
			try:
				if n1 > n2:
					otherHyp = cornerHash[(n1,n2)]
				else:
					otherHyp = cornerHash[(n2,n1)]
			
				oldCost = otherHyp[4]
				if newCost < oldCost:
					if n1 > n2:
						cornerHash[(n1,n2)] = hyp
					else:
						cornerHash[(n2,n1)] = hyp
			except:				
				if n1 > n2:
					cornerHash[(n1,n2)] = hyp
				else:
					cornerHash[(n2,n1)] = hyp			
	
		newCorners = []	
		totalHypotheses = []
		for k, v in cornerHash.items():
			print "corner constraint:", v[0], "->", v[1]
			totalHypotheses.append([v[0], v[1], v[2], v[3]])

			self.renderCornerConstraint(v[0], v[1], v[2], v[5], v[6], renderCount)
			renderCount += 1

			" corner constraints with priority of 3 "
			self.addPriorityEdge([v[0],v[1],v[2],v[3]], CORNER_PRIORITY)
			
			" put this pair of corners into matched bins "			
			n1 = v[0]
			n2 = v[1]
			point1_i = v[9]
			point2_i = v[10]
			tup1 = (n1,point1_i)
			tup2 = (n2,point2_i)

			" add the corner pair to a merged bin "
			self.addCornerBin(tup1, tup2)

			newCorners.append(tup1)
			newCorners.append(tup2)
		
		print "cornerBins:"
		for i in range(len(self.cornerBins)):
			print i, self.cornerBins[i]

		" remove duplicates and non-new node corners"
		d = {}
		for x in newCorners:
			if x[0] == nodeID:
				d[x] = 1
		newCorners = list(d.keys())

		if len(newCorners) > 0:
			
			binIndices = [None for i in range(len(newCorners))]
			for i in range(len(self.cornerBins)):
				bin = self.cornerBins[i]
				
				for j in range(len(binIndices)):
					tup1 = newCorners[j]
					" first check if either point exists already in a bin"
					if tup1 in bin:
						binIndices[j] = i
	
		
			" attempt a corner constraint between a newCorner and any corner in the same bin "
			internalHyps = []
			for i in range(len(binIndices)):
				tup1 = newCorners[i]
				tup2 = None
				
				" do not add repeat corner constraints between nodes "
				for newTup in self.cornerBins[binIndices[i]]:
					isReject = False
					
					if (tup1[0],newTup[0]) in self.edgePriorityHash:
						edges = self.edgePriorityHash[(tup1[0],newTup[0])]
						for edge in edges:
							" if any constraints are corner constraints "
							if edge[2] == CORNER_PRIORITY:
								isReject = True
						
					
					if (newTup[0],tup1[0]) in self.edgePriorityHash:							
						edges = self.edgePriorityHash[(newTup[0],tup1[0])]
						for edge in edges:
							" if any constraints are corner constraints "
							if edge[2] == CORNER_PRIORITY:
								isReject = True
					
					if not isReject and tup1[0] != newTup[0]:
						tup2 = newTup
						break
				
				if tup2 != None:
					n1 = tup1[0]
					n2 = tup2[0]
					point1_i = tup1[1]
					point2_i = tup2[1]
					corner1 = self.nodeHash[n1].cornerCandidates[point1_i]
					corner2 = self.nodeHash[n2].cornerCandidates[point2_i]
					point1 = corner1[0]
					point2 = corner2[0]
					ang1 = corner1[1]
					ang2 = corner2[1]
	
					print "TRYING INTERNAL corner hypothesis:", tup1, tup2
					offset, covar, cost, angDiff, oriDiff = self.makeCornerConstraint(n1, n2, point1, point2, ang1, ang2, self.a_hulls[n1], self.a_hulls[n2], self.b_hulls[n1], self.b_hulls[n2])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost <= 0.6:
						internalHyps.append([n1, n2, transform, covar, cost, point1, point2, angDiff, oriDiff, point1_i, point2_i])

			"""
			externalHyps = []
			for i in range(len(binIndices)):
				tup1 = newCorners[i]
				tup2 = None
				
				for j in range(len(self.cornerBins)):
					if j != binIndices[i]:
											
						" do not add repeat corner constraints between nodes "
						for newTup in self.cornerBins[j]:
							isReject = False
							
							if (tup1[0],newTup[0]) in self.edgePriorityHash:
								edges = self.edgePriorityHash[(tup1[0],newTup[0])]
								for edge in edges:
									" if any constraints are corner constraints "
									if edge[2] == CORNER_PRIORITY:
										isReject = True
								
							
							if (newTup[0],tup1[0]) in self.edgePriorityHash:							
								edges = self.edgePriorityHash[(newTup[0],tup1[0])]
								for edge in edges:
									" if any constraints are corner constraints "
									if edge[2] == CORNER_PRIORITY:
										isReject = True
							
							if not isReject and tup1[0] != newTup[0]:
								tup2 = newTup
								break
						
						if tup2 != None:
							n1 = tup1[0]
							n2 = tup2[0]
							point1_i = tup1[1]
							point2_i = tup2[1]
							corner1 = self.nodeHash[n1].cornerCandidates[point1_i]
							corner2 = self.nodeHash[n2].cornerCandidates[point2_i]
							point1 = corner1[0]
							point2 = corner2[0]
							ang1 = corner1[1]
							ang2 = corner2[1]
			
							print "TRYING EXTERNAL corner hypothesis:", tup1, tup2
							offset, covar, cost, angDiff, oriDiff = self.makeCornerConstraint(n1, n2, point1, point2, ang1, ang2, self.a_hulls[n1], self.a_hulls[n2], self.b_hulls[n1], self.b_hulls[n2])
							transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
							if cost <= 0.6:
								externalHyps.append([n1, n2, transform, covar, cost, point1, point2, angDiff, oriDiff, point1_i, point2_i])
			"""


			newCorners = []	
			totalHypotheses = []
			for v in internalHyps:
				print "internal corner constraint:", v[0], "->", v[1]
				totalHypotheses.append([v[0], v[1], v[2], v[3]])

				self.renderCornerConstraint(v[0], v[1], v[2], v[5], v[6], renderCount)
				renderCount += 1
				
				" corner constraints with priority of 3 "
				self.addPriorityEdge([v[0],v[1],v[2],v[3]], CORNER_PRIORITY)
	
				" put this pair of corners into matched bins "			
				n1 = v[0]
				n2 = v[1]
				point1_i = v[9]
				point2_i = v[10]
				tup1 = (n1,point1_i)
				tup2 = (n2,point2_i)
				newCorners.append(tup1)
				newCorners.append(tup2)
				matchBins = []
				for i in range(len(self.cornerBins)):
					
					bin = self.cornerBins[i]
					" first check if either point exists already in a bin"
					if tup1 in bin:
						matchBins.append(i)
					elif tup2 in bin:
						matchBins.append(i)
		
				" if not shared, create new bin "
				if len(matchBins) == 0:
					self.cornerBins.append([tup1, tup2])
				else:						
					" we have to merge into the old bins "
					if len(matchBins) == 1:
						" only 1 bin, just add the tuples that don't already exist "
						bin_i = matchBins[0]
						
						if not tup1 in self.cornerBins[bin_i]:
							self.cornerBins[bin_i].append(tup1)
						if not tup2 in self.cornerBins[bin_i]:
							self.cornerBins[bin_i].append(tup2)
							
					else:
						" we have more than one matching bin, so now we need to merge "
						newSuperBin = []
						for i in matchBins:
							bin = self.cornerBins[i]
							
							for item in bin:
								if not item in newSuperBin:
									newSuperBin.append(item)
		
						newCornerBins = [newSuperBin]
						
						for i in range(len(self.cornerBins)):
							if not i in matchBins:
								newCornerBins.append(self.cornerBins[i])
		
						self.cornerBins = newCornerBins

			
			print "cornerBins2:"
			for i in range(len(self.cornerBins)):
				print i, self.cornerBins[i]

	def getPathTerms(self):

		terms = []

		numPaths = len(self.paths)
		for k in range(numPaths):

			path = self.paths[k]
			
			pathSpline = SplineFit(path, smooth=0.1)
			vecPoints = pathSpline.getUniformSamples()
			
			angleSum = 0.0
			for p in vecPoints[-11:]:
				angleSum += p[2]
				
			angleAvg = angleSum / 10.
			
			terms.append([vecPoints[-1][0],vecPoints[-1][1], angleAvg])

		return terms

	def trimPaths(self, paths):

		trimmedPaths = []
		
		if len(paths) <= 1:
			for k, v in paths.iteritems():
				trimmedPaths.append(v)
				
			return trimmedPaths

		" path0 is root, and path1 is a child of path0 "
		#parents = [None, [0, self.junctionNodeID, self.localJunctionPoint]]
		parents = self.pathParents

		" information on parentage is required to determine which path belongs to which "

		" global junction point for parent 0 and child 1 "

		print "pathParents:", parents
		
		for pathID in range(len(paths)):
			
			if parents[pathID][0] == None:
				trimmedPaths.append(deepcopy(paths[pathID]))

			else:

				path1 = paths[parents[pathID][0]]
				path2 = paths[pathID]
				
				parentPathID = parents[pathID][0]
				childPathID = pathID

				#secP1, secP2 = self.getOverlapDeparture(path1, path2)				
				secP1, secP2 = self.getOverlapDeparture(parentPathID, childPathID, paths)				

				minI_1 = 0		
				minI_2 = 0
				minDist_1 = 1e100
				minDist_2 = 1e100		
				for i in range(len(path2)):
					pnt = path2[i]
					dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
					dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
				
					if dist1 < minDist_1:
						minDist_1 = dist1
						minI_1 = i

					if dist2 < minDist_2:
						minDist_2 = dist2
						minI_2 = i

				term0 = path2[0]
				termN = path2[-1]

				" smallest distance is the terminal point "
				dist0_1 = sqrt((term0[0]-secP1[0])**2 + (term0[1]-secP1[1])**2)
				distN_1 = sqrt((termN[0]-secP1[0])**2 + (termN[1]-secP1[1])**2)
				dist0_2 = sqrt((term0[0]-secP2[0])**2 + (term0[1]-secP2[1])**2)
				distN_2 = sqrt((termN[0]-secP2[0])**2 + (termN[1]-secP2[1])**2)
				print "terminal distances:", dist0_1, distN_1, dist0_2, distN_2
				
				distList = [dist0_1, distN_1, dist0_2, distN_2]
				
				minI = -1
				minDist = 1e100
				for i in range(len(distList)):
					if distList[i] < minDist:
						minDist = distList[i]
						minI = i

				junctionNodeID = parents[pathID][1]
				localJunctionPoint = parents[pathID][2]

				poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
				junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)
				
				
				if minDist_1 < minDist_2:
					junctionPoint_K = secP1
					juncI_K = minI_1
				else:
					junctionPoint_K = secP2
					juncI_K = minI_2
						
				#junctionPoint = junctionPoint_K
				#juncI = juncI_K
				
				
				minDist2 = 1e100
				juncI = 0		
				for i in range(len(path2)):
					pnt = path2[i]
					dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
				
					if dist < minDist2:
						minDist2 = dist
						juncI = i
 				
 				
 				
				print "len(path1):", len(path1)
				print "len(path2):", len(path2)
				print "juncI:", juncI
				print "minDist:", minDist_1, minDist_2
				
				" now we have closest point to departure point. "
				" Which side is the departing side? "	
				#sum1 = self.getOverlapSum(path1, path2[0:minI2+1] + [junctionPoint])
				#sum2 = self.getOverlapSum(path1, [junctionPoint] + path2[minI2:])
				

						
				if minI == 0:
					"secP1 is terminal 0"
					index = juncI-10
					if index < 1:
						index = 1
					#newPath2 = path2[:index+1] + [junctionPoint]
					newPath2 = path2[:index+1]
					newPath2.reverse()
				
				elif minI == 1:
					"secP1 is terminal N"
					index = juncI+10
					if index >= len(path2):
						index = len(path2)-1
					#newPath2 = [junctionPoint] + path2[index:]
					newPath2 = path2[index:]
					
				elif minI == 2:
					"secP2 is terminal 0"
					index = juncI-10
					if index < 1:
						index = 1
					#newPath2 = path2[:index+1] + [junctionPoint]
					newPath2 = path2[:index+1]
					newPath2.reverse()
					
				elif minI == 3:
					"secP2 is terminal N"
					index = juncI+10
					if index >= len(path2):
						index = len(path2)-1
					#newPath2 = [junctionPoint] + path2[index:]
					newPath2 = path2[index:]
				
				else:
					print "no terminal found"
					raise
								
				#if minI_1 > minI_2:
				#	path2[minI_2:minI_1+1]
				#else:
				#	path2[minI_1:minI_2+1]
				
				
				#if startI == 0:
				#	newPath2 = path2[startI:endI+1] + [junctionPoint]
				#	newPath2.reverse()
				#else:
				#	newPath2 = [junctionPoint] + path2[startI:endI+1]
					

				#print "Sums:", sum1, sum2
				
				#if sum1 > sum2:
				#	index = minI2+1-10
				#	if index < 1:
				#		index = 1
				#	newPath2 = path2[0:index] + [junctionPoint]
				#	newPath2.reverse()
				#else:
				#	index = minI2+10
				#	if index >= len(path2):
				#		index = len(path2)-1
				#	newPath2 = [junctionPoint] + path2[index:]

				" convert path so that the points are uniformly distributed "
				
				
				max_spacing = 0.08
				#max_spacing = 0.04
				#max_spacing = 0.1
				newPath3 = []

				" make sure path is greater than 5 points "
				while len(newPath3) <= 5:
	
					max_spacing /= 2
					
					newPath3 = [copy(newPath2[0])]
										
					for i in range(len(newPath2)-1):
						p0 = newPath2[i]
						p1 = newPath2[(i+1)]
						dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
			
						vec = [p1[0]-p0[0], p1[1]-p0[1]]
						vec[0] /= dist
						vec[1] /= dist
						
						
						if dist > max_spacing:
							" cut into pieces max_spacing length or less "
							numCount = int(floor(dist / max_spacing))
							
							for j in range(1, numCount+1):
								newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
								newPath3.append(newP)
	
						newPath3.append(copy(p1))			
				
				trimmedPaths.append(deepcopy(newPath3))

		"""
		self.drawTrimmedPaths(self.trimmedPaths)
		"""
		
		
		return trimmedPaths


	def getOrderedOverlappingPaths(self, nodeID):

		overlappedPaths = []
		
		"FIXME:  add departure detection and junction continuity checking "
		"        to reject paths that we are only incident to "
	
		print "getOrderedOverlappingPaths(", nodeID, ")"
	
		overlapSums = {}
	
		for pathID in range(len(self.trimmedPaths)):
			sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
			
			print "overlap sum", pathID, "=", sum1
			
			if sum1 <= 1e10:
				overlappedPaths.append(pathID)
				overlapSums[pathID] = sum1

		print "overlapSums:", overlapSums

		orderedPathIDs = self.getPathOrdering(nodeID, overlappedPaths)

		print "orderedPathIDs:", orderedPathIDs

		departures = []
		interiors = []
		depPoints = []

		for pathID in orderedPathIDs:
			departurePoint1, depAngle2, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID)
			departures.append([isExist1,isExist2])
			interiors.append([isInterior1, isInterior2])
			depPoints.append([departurePoint1, departurePoint2])
		
		print "departures:", departures
		print "interiors:", interiors
			
		" determine which paths are leaves "
		isAParent = [False for k in range(self.numPaths)]
		for k in orderedPathIDs:
			currPath = self.pathParents[k]
			currParent = currPath[0]
			if currParent != None:
				isAParent[currParent] = True

		print "isAParent:", isAParent
		
		toBeRemoved = []
		seqSize = len(orderedPathIDs)
		
		" TODO: is this necessary still? "
		" check for proper parent-child relationship and "
		" that there is interior departure between paths "

		
		
		for k in range(seqSize-1):
			pathID1 = orderedPathIDs[k]
			pathID2 = orderedPathIDs[k+1]

			" there should be a interior departure between these paths on the parent path "
			if self.pathParents[pathID1][0] == pathID2:
				pass
				#parentPath = pathID2
				
			elif self.pathParents[pathID2][0] == pathID1:
				pass
				#parentPath = pathID1

					
			else:
				print "ordered path IDs could not determine relationship", pathID1, pathID2, orderedPathIDs
				raise
		
		"""		
		for k in range(seqSize-1):
			pathID1 = orderedPathIDs[k]
			pathID2 = orderedPathIDs[k+1]

			" there should be a interior departure between these paths on the parent path "
			if self.pathParents[pathID1][0] == pathID2:
				parentPath = pathID2

				" if the departure doesn't show, remove the least overlapped path section "
				if not interiors[k+1][0]:
					if overlapSums[pathID1] > overlapSums[pathID2]:
						print "remove pathID1", pathID1
						#toBeRemoved.append(k)
					else:
						print "remove pathID2", pathID2
						#toBeRemoved.append(k+1)
					
					if k != 0 and k != seqSize-1:
						print "paths in the middle violate a departure constraint!"
						print pathID1, pathID2, orderedPathIDs
						print departures, interiors
						print overlapSums
						print toBeRemoved
						raise
				
				
			elif self.pathParents[pathID2][0] == pathID1:
				parentPath = pathID1
				if not interiors[k][1]:
					if overlapSums[pathID1] > overlapSums[pathID2]:
						print "remove pathID1", pathID1
						#toBeRemoved.append(k)
					else:
						print "remove pathID2", pathID2
						#toBeRemoved.append(k+1)

					if k != 0 and k != seqSize-1:
						print "paths in the middle violate a departure constraint!"
						print pathID1, pathID2, orderedPathIDs
						print departures, interiors
						print overlapSums
						print toBeRemoved
						raise
					
			else:
				print "ordered path IDs could not determine relationship", pathID1, pathID2, orderedPathIDs
				raise
		"""
		
		print "toBeRemoved:", toBeRemoved

		finalPathIDs = []
		for k in range(seqSize):
			
			if toBeRemoved.count(k) == 0:
				finalPathIDs.append(orderedPathIDs[k])
			
		print "finalPathIDs:", finalPathIDs	
		return finalPathIDs

	def getPathOrdering(self, nodeID, pathIDs):

		plotIter = False

		"""
		for pathID in pathIDs:
			
			path = self.paths[pathID]
			hull = self.hulls[pathID]
			overlapCost1 = self.getOverlapCondition(path, nodeID1)
			overlapCost2 = self.getOverlapCondition(path, nodeID2)


			departurePoint1, isInterior1, isExist1 = self.getDeparturePoint(path, nodeID1, hull)
			departurePoint2, isInterior2, isExist2 = self.getDeparturePoint(path, nodeID2, hull)
		"""

		node2 = self.nodeHash[nodeID]

		hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		
		"""
		if node2.isBowtie:			
			hull2 = computeBareHull(node2, sweep = False, static = True)
			hull2.append(hull2[0])
			medial2 = node2.getStaticMedialAxis()

		else:
			hull2 = computeBareHull(node2, sweep = False)
			hull2.append(hull2[0])
			medial2 = node2.getMedialAxis(sweep = False)

		" take the long length segments at tips of medial axis"
		edge1 = medial2[0:2]
		edge2 = medial2[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]

		
		" find the intersection points with the hull "
		interPoints = []
		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect2, point2 = Intersect(edge2, hullEdge)
			if isIntersect2:
				interPoints.append(point2)
				break
		
		" replace the extended edges with a termination point at the hull edge "			
		medial2 = medial2[1:-2]
		if isIntersect1:
			medial2.insert(0, point1)
		if isIntersect2:
			medial2.append(point2)
		

		" cut off the tail of the non-sweeping side "
		TAILDIST = 0.5
		"""
		"""
		if nodeID % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]

		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
		"""
		
		estPose2 = node2.getGlobalGPACPose()		
			
		#minMatchDist2 = 0.2
		#minMatchDist2 = 0.5
		#minMatchDist2 = 2.0
			
		" set the initial guess "
		poseOrigin = Pose(estPose2)
						
		localizedPaths = []
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			
			localPath = []
			for pnt in path:
				localPath.append(poseOrigin.convertGlobalToLocal(pnt))
			
			localizedPaths.append(localPath)
			
		" find minimum distance and its associated path "					
		#supportSpline = SplineFit(localSupport, smooth=0.1)		
		#vecPoints1 = supportSpline.getUniformSamples()
		#supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			
		
		pathSelect = []
		#minMatchDist2 = 0.5
		minMatchDist2 = 1.0

		for i in range(len(vecPoints2)):

			minDist = 1e100
			minPathID = -1
			
			p_2 = vecPoints2[i]

			for j in range(len(localizedPaths)):


				path = localizedPaths[j]
				medialSpline1 = SplineFit(path, smooth=0.1)
				vecPoints1 = medialSpline1.getUniformSamples()
				supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
	
				" for every transformed point of A, find it's closest neighbor in B "
				try:
					p_1, minI, minDist_I = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
		
					if minDist_I <= minMatchDist2:
						C2 = points2[i][2]
						C1 = supportPoints[minI][2]

						ax = supportPoints[minI][0]
						ay = supportPoints[minI][1]		
						bx = points2[i][0]
						by = points2[i][1]
				
						c11 = C1[0][0]
						c12 = C1[0][1]
						c21 = C1[1][0]
						c22 = C1[1][1]
								
						b11 = C2[0][0]
						b12 = C2[0][1]
						b21 = C2[1][0]
						b22 = C2[1][1]	
						
						dist = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])

						if dist < minDist:
							minDist = dist
							minPathID = pathIDs[j]
						#" we store the untransformed point, but the transformed covariance of the A point "
						#support_pairs.append([points2[i],supportPoints[minI],C2,C1])

				
				except:
					#raise
					pass

				#print "minPathID:", minPathID

			if minPathID != -1:
				pathSelect.append(minPathID)
		
		print "pathSelect:", pathSelect
		
		" for each point on the medial axis, find it's closest pathID "
		"""
		pathSelect = []
		#for pnt2 in poly2:
		for pnt2 in points2:

			minDist = 1e100
			minPathID = -1
			
			for i in range(len(localizedPaths)):


				path = localizedPaths[i]

				medialSpline1 = SplineFit(path, smooth=0.1)
				
				
				#p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)

				
				points1 = gen_icp.addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
				
				#for pnt2 in path:
				for pnt1 in points1:
					
					#a = pnt1[0]
					#b = pnt2[0]
					Ca = pnt1[2]
					Cb = pnt2[2]
			
					ax = pnt1[0]
					ay = pnt1[1]		
					bx = pnt2[0]
					by = pnt2[1]
			
					c11 = Ca[0][0]
					c12 = Ca[0][1]
					c21 = Ca[1][0]
					c22 = Ca[1][1]
							
					b11 = Cb[0][0]
					b12 = Cb[0][1]
					b21 = Cb[1][0]
					b22 = Cb[1][1]	
					
					dist = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
										
					
					#dist = sqrt((pnt2[0]-pnt1[0])**2 + (pnt2[0]-pnt1[0])**2)
					if dist < minDist:
						minDist = dist
						minPathID = pathIDs[i]
			
			pathSelect.append(minPathID)
		"""
		
		" find the average position of each path ID"
		avgIDPosition = {}
		avgIDCount = {}
		for pathID in pathIDs:
			avgIDPosition[pathID] = 0
			avgIDCount[pathID] = 0
			
		for i in range(len(pathSelect)):
			avgIDPosition[pathSelect[i]] += i
			avgIDCount[pathSelect[i]] += 1


		for pathID in pathIDs:
			if avgIDCount[pathID] != 0:
				avgIDPosition[pathID] /= avgIDCount[pathID]
			else:
				del avgIDPosition[pathID]
			
		" now sort out the order "
		orderedPaths = sorted(avgIDPosition, key=avgIDPosition.__getitem__)		
		
		#print "orderedPaths:", orderedPaths
		
		
		if plotIter:
			pylab.clf()
			
			xP = []
			yP = []
			
			for p in poly2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP)
			
			
			for i in range(len(localizedPaths)):
				xP = []
				yP = []
				path = localizedPaths[i]
				for p in path:
					xP.append(p[0])
					yP.append(p[1])
				
				pylab.plot(xP,yP)
	
			pylab.title("%d: %s %s" % (nodeID, repr(pathIDs), repr(orderedPaths)))
			pylab.savefig("orderedPath_%04u.png" % self.orderCount)
			
			self.orderCount += 1
		
		return orderedPaths

	def getDeparturePoint(self, currPath, nodeID):
		
		isExist1 = False
		isInterior1 = False
		departurePoint1 = 0
		angle1 = 0.0
		isExist2 = False
		isInterior2 = False
		departurePoint2 = 0
		angle2 = 0.0
		
		if len(currPath) == 0:
			return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
		
		node2 = self.nodeHash[nodeID]

		hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		
		estPose2 = node2.getGlobalGPACPose()		
		
		"Assumption:  one section of the medial axis is closely aligned with the path "
		poseOrigin = Pose(estPose2)
		
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		points2 = medialSpline2.getUniformSamples()

		pathSpline = SplineFit(currPath, smooth=0.1)
		pathPoints = pathSpline.getUniformSamples()

		
		points2_offset = []
		for p in points2:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points2_offset.append(result)


		" tip angles "
		angSum1 = 0.0
		angSum2 = 0.0
		angs1 = []
		angs2 = []
		phi1 = normalizeAngle(points2_offset[0][2])
		phi2 = normalizeAngle(points2_offset[-1][2])
		for i in range(10):
			ang1 = normalizeAngle(points2_offset[i][2]-phi1)
			ang2 = normalizeAngle(points2_offset[-i-1][2]-phi2)
			angSum1 += ang1
			angSum2 += ang2
			
			angs1.append(ang1+phi1)
			angs2.append(ang2+phi2)

		angle1 = angSum1 / 10.0 + phi1
		angle2 = angSum2 / 10.0 + phi2

		" invert one angle so opposite tips have opposite angles "
		angle1 = normalizeAngle(angle1 + pi)

		print "ang1:", angle1, angs1
		print "ang2:", angle2, angs2
		print "diff:", diffAngle(angle1, angle2)

		distances = []
		indices = []
		for i in range(0,len(points2_offset)):
			p_2 = points2_offset[i]
			#p_1, minDist = gen_icp.findClosestPointInB(pathPoints, p_2, [0.0,0.0,0.0])
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints, p_2)
			distances.append(minDist)
			indices.append(i_1)
		
		
		" Compute the front and back departure points by finding the inflection point on the distance curve "
		" these indices become frontDepI and backDepI respectively "
		maxFront = distances[0]
		maxBack = distances[-1]

		currI = 1
		try:
			while distances[currI+3] < maxFront:
				maxFront = distances[currI]
				currI += 1
		except:
			pass
		
		frontDepI = currI
		frontPoint = [frontDepI, distances[frontDepI]]

		" FIXME:  index out of bounds case "
		currI = 2
		try:
			while distances[-currI-3] < maxBack:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass

		backDepI = len(distances) - currI
		backPoint = [backDepI, distances[backDepI]]

		"reset to the maximum distances "
		maxFront = distances[0]
		maxBack = distances[-1]
		
		" perform a line fit on the distance curves from the departure point to the tip"	
		#(ar1,br1)= scipy.polyfit(range(0,frontDepI+1),distances[0:frontDepI+1],1)
		#(ar2,br2)= scipy.polyfit(range(backDepI,len(distances)),distances[backDepI:],1)

		#t1 = range(0,frontDepI+1)
		#xr1 = scipy.polyval([ar1,br1],t1)			
		
		#t2 = range(backDepI,len(distances))
		#xr2 = scipy.polyval([ar2,br2],t2)			
		
		" compute the average of distances between departure and termination "
		#frontAvg = 0.0
		#for i in range(0,frontDepI+1):
		#	frontAvg += distances[i]
		#frontAvg /= (frontDepI+1)

		#backAvg = 0.0
		#for i in range(backDepI,len(distances)):
		#	backAvg += distances[i]
		#backAvg /= len(distances)-backDepI

		" for all the points matched from the local curve to the path curve "
		" count the number of times they are matched "
		foo = indices[0:frontDepI+1]
		d1 = {}
		for i in set(foo):
			d1[i] = foo.count(i)

		foo = indices[backDepI:]
		d2 = {}
		for i in set(foo):
			d2[i] = foo.count(i)

		" find the point that has the most matches "
		max1 = max(d1, key=d1.get)
		max2 = max(d2, key=d2.get)

		" compute the selected point index by average instead of by maximum"
		#frontI = 0.0
		#for i in range(0,frontDepI+1):
		#	frontI += indices[i]
		#frontI /= (frontDepI+1)

		#backI = 0.0
		#for i in range(backDepI,len(distances)):
		#	backI += indices[i]
		#backI /= len(distances)-backDepI

		tipMatch1 = pathPoints[indices[0]]

		" discrepancy distance between tip closest point and average departure point "
		dist1 = sqrt((tipMatch1[0]-pathPoints[max1][0])**2 + (tipMatch1[1] - pathPoints[max1][1])**2)

		tipMatch2 = pathPoints[indices[-1]]

		" discrepancy distance between tip closest point and average departure point "
		dist2 = sqrt((tipMatch2[0]-pathPoints[max2][0])**2 + (tipMatch2[1] - pathPoints[max2][1])**2)

		if maxFront > 0.5:
			departurePoint1 = pathPoints[max1]
			isExist1 = True



			if max1 == 0 or max1 == len(pathPoints)-1:
				isInterior1 = False
			else:
				isInterior1 = True

			# Test purposes only
			#if nodeID == 4 or nodeID == 5:
			#	isInterior1 = True


		if maxBack > 0.5:
			departurePoint2 = pathPoints[max2]
			isExist2 = True


			if max2 == 0 or max2 == len(pathPoints)-1:
				isInterior2 = False
			else:
				isInterior2 = True



			
		"""
		if maxFront > 0.5 and maxFront > maxBack:
			departurePoint = pathPoints[max1]
			isExist = True

			if max1 == 0 or max1 == len(pathPoints)-1:
				isInterior = False
			else:
				isInterior = True

		if maxBack > 0.5 and maxBack > maxFront:
			departurePoint = pathPoints[max2]
			isExist = True

			if max2 == 0 or max2 == len(pathPoints)-1:
				isInterior = False
			else:
				isInterior = True
		"""	
			
		" sum of closest points on front and back "
		" select the one with minimal cost "
		
		if False:
			pylab.clf()
			xP = range(len(points2_offset))
			yP = distances
			pylab.plot(xP,yP, color ='b')
			
			if maxFront > 0.5:
				xP = [frontPoint[0]]
				yP = [frontPoint[1]]
				pylab.scatter(xP,yP, color='b')
	
			if maxBack > 0.5:
				xP = [backPoint[0]]
				yP = [backPoint[1]]
				pylab.scatter(xP,yP, color='b')
			
			pylab.xlim(0,200)
			pylab.ylim(0,2)
			#pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
			pylab.title("nodeID %d: %1.2f %1.2f %d %d %d" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2))
			pylab.savefig("distances_%04u.png" % self.pathPlotCount)

		if True:
			pylab.clf()
			xP = []
			yP = []
			for p in points2_offset:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')
	
			if True:	
				xP = [pathPoints[max1][0]]
				yP = [pathPoints[max1][1]]
				pylab.scatter(xP,yP, color='b')		
	
				xP = [tipMatch1[0]]
				yP = [tipMatch1[1]]
				pylab.scatter(xP,yP, color='r')		
	
	
			if True:
				xP = [pathPoints[max2][0]]
				yP = [pathPoints[max2][1]]
				pylab.scatter(xP,yP, color='b')		
	
				xP = [tipMatch2[0]]
				yP = [tipMatch2[1]]
				pylab.scatter(xP,yP, color='r')		
	
	
			xP = []
			yP = []
			for p in currPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			#pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
			pylab.title("nodeID %d: %1.2f, %1.2f, %1.2f, %1.2f, %s %s" % (nodeID, dist1, dist2, angle1, angle2, repr([isExist1, isExist2]), repr([isInterior1, isInterior2])))
			pylab.savefig("departure_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1
		
		
		#return departurePoint, isInterior, isExist
		return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2
		

	" returns the endpoints of a subpath of path 2 that does not overlap path 1 "

	" get the trimmed version of child and parent paths that are overlapping in some fashion "
	def getOverlapDeparture(self, parentPathID, childPathID, paths):

		"Assumption:  one section of the medial axis is closely aligned with the path "		
			
		plotIter = False
		print "getOverlapDeparture():"
		
		path1 = paths[parentPathID]
		path2 = paths[childPathID]
		
		
		isExist1 = False
		isInterior1 = False
		departurePoint1 = 0
		isExist2 = False
		isInterior2 = False
		departurePoint2 = 0


		branchNodeID = self.pathParents[childPathID][1]
		poseOrigin = Pose(self.nodeHash[branchNodeID].getGlobalGPACPose())
		localJunctionPoint = self.pathParents[childPathID][2]
		
		
		poseOrigin2 = Pose(self.nodeHash[branchNodeID].getEstPose())
		globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)
		
		
		" orienting the medial axis of the branch node correctly "
		branchDir = self.pathParents[childPathID][3]
		
		branchHull, branchMedial = computeHullAxis(branchNodeID, self.nodeHash[branchNodeID], tailCutOff = False)

		branchMedialSpline = SplineFit(branchMedial, smooth=0.1)
		branchPoints = branchMedialSpline.getUniformSamples()
		
		globalBranchPoints = []
		for p in branchPoints:
			result = poseOrigin.convertLocalToGlobal(p)
			globalBranchPoints.append(result)
		
		if not branchDir:
			" backward "
			globalBranchPoints.reverse()
			
		globalBranchSpline = SplineFit(globalBranchPoints, smooth = 0.1)

			

		" return exception if we receive an invalid path "		
		if len(path1) == 0:
			print "path1 has zero length"
			raise
			#return departurePoint1, isInterior1, isExist1, departurePoint2, isInterior2, isExist2


		" make sure the overlap of both paths are oriented the same way "
		path1Spline = SplineFit(path1, smooth=0.1)
		path2Spline = SplineFit(path2, smooth=0.1)

		path2Reverse = deepcopy(path2)
		path2Reverse.reverse()
		path2SplineReverse = SplineFit(path2Reverse, smooth=0.1)
		
		overlapMatch = []
		angleSum1 = 0.0
		angleSum2 = 0.0
		for i in range(0,len(path1)):
			p_1 = path1[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(path2, p_1)
			if minDist < 0.1:
				overlapMatch.append((i,j,minDist))

				pathU1 = path1Spline.findU(p_1)	
				pathU2 = path2Spline.findU(p_2)	
				pathU2_R = path2SplineReverse.findU(p_2)	

				pathVec1 = path1Spline.getUVector(pathU1)
				pathVec2 = path2Spline.getUVector(pathU2)
				pathVec2_R = path2SplineReverse.getUVector(pathU2_R)

				val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
				if val1 > 1.0:
					val1 = 1.0
				elif val1 < -1.0:
					val1 = -1.0
				ang1 = acos(val1)
				
				val2 = pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1]
				if val2 > 1.0:
					val2 = 1.0
				elif val2 < -1.0:
					val2 = -1.0
				ang2 = acos(val2)
				
				angleSum1 += ang1
				angleSum2 += ang2
	
		" select global path orientation based on which has the smallest angle between tangent vectors "
		print "angle sums for orienting child path with parent path:", angleSum1, angleSum2

		if angleSum1 > angleSum2:
			orientedPath2 = path2Reverse
		else:
			orientedPath2 = path2
			
		path2Spline = SplineFit(orientedPath2, smooth=0.1)



		" for each point on the child path, find its closest pair on the parent path "
		
		pathPoints1 = path1Spline.getUniformSamples()
		pathPoints2 = path2Spline.getUniformSamples()

		distances = []
		indices = []
		for i in range(0,len(pathPoints2)):
			p_2 = pathPoints2[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
			
			" keep the distance information "
			distances.append(minDist)
			" and the associated index of the point on the parent path "
			indices.append(i_1)

		" match distances of the tip points of child path "		
		maxFront = distances[0]
		maxBack = distances[-1]
		print "match distances of tip points:", maxFront, maxBack

		" walk back from tip point until we have a non-monotic increase in match distance "
		" this becomes our departure point "
		
		"TODO:  why is there a 3-offset in this comparison? "
		currI = 1
		try:
			while distances[currI+3] < maxFront:
				maxFront = distances[currI]
				currI += 1
		except:
			pass
		
		" departure index on child path "
		frontDepI = currI
		
		" departure index of child path and match distance "
		frontPoint = [frontDepI, distances[frontDepI]]

		" FIXME:  index out of bounds case "
		currI = 2
		try:
			while distances[-currI-3] < maxBack:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass

		" departure index on child path "
		backDepI = len(distances) - currI

		" departure index of child path and match distance "
		backPoint = [backDepI, distances[backDepI]]

		print "lengths of parent and child paths:", len(pathPoints1), len(pathPoints2)
		print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint


		"reset to the tip match distance "
		maxFront = distances[0]
		maxBack = distances[-1]
		
		
		" count the number of times a matched point on the parent path is used "
		foo = indices[0:frontDepI+1]
		d1 = {}
		for i in set(foo):
			d1[i] = foo.count(i)

		foo = indices[backDepI:]
		d2 = {}
		for i in set(foo):
			d2[i] = foo.count(i)

		" select the match point that is used the most on the parent path "
		max1 = max(d1, key=d1.get)
		max2 = max(d2, key=d2.get)


		" NOTE:  This guarantees detection of a departure.  We've assumed that there exists a"
		" departure between these two paths "
		
		" compute two candidate departure points "
		
		" NOTE:  departure point on parent path is not used here "
		if True:

			#departurePoint1 = pathPoints1[max1]
			isExist1 = True

			" if the highest match point is either of the tips of the parent path, then this is an external departure "
			if max1 == 0 or max1 == len(pathPoints1)-1:
				isInterior1 = False
			else:
				isInterior1 = True

		if True:
			
			#departurePoint2 = pathPoints1[max2]
			isExist2 = True

			" if the highest match point is either of the tips of the parent path, then this is an external departure "
			if max2 == 0 or max2 == len(pathPoints1)-1:
				isInterior2 = False
			else:
				isInterior2 = True

		print "isExist1 =", isExist1, "isInterior1 =", isInterior1
		print "isExist2 =", isExist2, "isInterior2 =", isInterior2

		" sum of closest points on front and back "
		" select the one with minimal cost "		
		
		" now we compare our candidates to the known direction and location of the branching point "
		angleSum1 = 0.0
		overlapSum1 = 0.0
		matchCount1 = 0
		
		" path section for our front departure hypothesis "
		pathSec1 = pathPoints2[:frontDepI+1]
		pathSec1.reverse()
		
		if len(pathSec1) > 5:
			pathSec1Spline = SplineFit(pathSec1, smooth = 0.1)		
			
			
			for i in range(0,len(pathSec1)):
				p_1 = pathSec1[i]
				p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
				if minDist < 0.1:
					overlapMatch.append((i,j,minDist))
					matchCount1 += 1
					
					overlapSum1 += minDist
	
					pathU1 = pathSec1Spline.findU(p_1)	
					pathUB = globalBranchSpline.findU(p_2)	
	
					pathVec1 = pathSec1Spline.getUVector(pathU1)
					pathVecB = globalBranchSpline.getUVector(pathUB)
	
					val1 = pathVec1[0]*pathVecB[0] + pathVec1[1]*pathVecB[1]
					if val1 > 1.0:
						val1 = 1.0
					elif val1 < -1.0:
						val1 = -1.0
					ang1 = acos(val1)
					
					angleSum1 += ang1

			if matchCount1 > 0:
				angleSum1 /= matchCount1
				overlapSum1 /= matchCount1
		

		angleSum2 = 0.0
		overlapSum2 = 0.0
		matchCount2 = 0
		
		" path section for our back departure hypothesis "
		pathSec2 = pathPoints2[backDepI:]

		if len(pathSec2) > 5:
			pathSec2Spline = SplineFit(pathSec2, smooth = 0.1)		
	
	
			for i in range(0,len(pathSec2)):
				p_1 = pathSec2[i]
				p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
				if minDist < 0.1:
					overlapMatch.append((i,j,minDist))
					matchCount2 += 1
					
					overlapSum2 += minDist
	
					pathU2 = pathSec2Spline.findU(p_1)	
					pathUB = globalBranchSpline.findU(p_2)	
	
					pathVec2 = pathSec2Spline.getUVector(pathU2)
					pathVecB = globalBranchSpline.getUVector(pathUB)
	
					val2 = pathVec2[0]*pathVecB[0] + pathVec2[1]*pathVecB[1]
					if val2 > 1.0:
						val2 = 1.0
					elif val2 < -1.0:
						val2 = -1.0
					ang2 = acos(val2)
					
					angleSum2 += ang2

			if matchCount2 > 0:
				angleSum2 /= matchCount2
				overlapSum2 /= matchCount2

		print "pathSec1 hypothesis angle and overlap sum and match count:", angleSum1, overlapSum1, matchCount1
		print "pathSec2 hypothesis angle and overlap sum and match count:", angleSum2, overlapSum2, matchCount2

		" distance of departure point from known junction point "
		juncDist1 = -1
		juncDist2 = -1
		if len(pathSec1) > 0:
			p0 = pathSec1[0]
			juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

		if len(pathSec2) > 0:
			p0 = pathSec2[0]
			juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

		print "pathSec1 hypothesis discrepancy distance:", juncDist1
		print "pathSec2 hypothesis discrepancy distance:", juncDist2

		secP1 = []
		secP2 = []

		if matchCount1 > 0 and matchCount2 > 0:
			if juncDist1 < juncDist2:
				secP1 = pathPoints2[0]
				secP2 = pathPoints2[frontDepI]
			else:
				secP1 = pathPoints2[backDepI]
				secP2 = pathPoints2[len(distances)-1]

		elif matchCount1 > 0:
			secP1 = pathPoints2[0]
			secP2 = pathPoints2[frontDepI]
			
		elif matchCount2 > 0:
			secP1 = pathPoints2[backDepI]
			secP2 = pathPoints2[len(distances)-1]
		
		else:
			print "no overlap detected!"
			raise
		
		"""
		if isInterior1 and isInterior2:
			if juncDist1 < juncDist2:
				secP1 = pathPoints2[0]
				secP2 = pathPoints2[frontDepI]
			else:
				secP1 = pathPoints2[backDepI]
				secP2 = pathPoints2[len(distances)-1]
		
		elif isExist1 and isInterior1:
			secP1 = pathPoints2[0]
			secP2 = pathPoints2[frontDepI]
			
		elif isExist2 and isInterior2:
			secP1 = pathPoints2[backDepI]
			secP2 = pathPoints2[len(distances)-1]
		
		else:
			print "Warning: external departure found"
			if isExist1:
				secP1 = pathPoints2[0]
				secP2 = pathPoints2[frontDepI]
				
			elif isExist2:
				secP1 = pathPoints2[backDepI]
				secP2 = pathPoints2[len(distances)-1]
		"""

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in path2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(0.5,0.5,1.0))
	
			if True:	
				P1 = pathPoints2[0]
				P2 = pathPoints2[frontDepI]
				
				P3 = pathPoints2[backDepI]
				P4 = pathPoints2[-1]

				" 0, frontDepI "
				xP = [P1[0], P2[0]]
				yP = [P1[1], P2[1]]
				pylab.scatter(xP,yP, color='b')		

				" backDepI, -1 "
				xP = [P3[0], P4[0]]
				yP = [P3[1], P4[1]]
				pylab.scatter(xP,yP, color='g')		
	
				
				#xP = [P2[0], P3[0]]
				#yP = [P2[1], P3[1]]
				#xP = [P1[0], P2[0], P3[0], P4[0]]
				#yP = [P1[1], P2[1], P3[1], P4[1]]
	
			pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')		
	
			xP = []
			yP = []
			for p in path1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(1.0,0.5,0.5))
	
			#xP = []
			#yP = []
			#for p in globalBranchPoints:
			#	xP.append(p[0])
			#	yP.append(p[1])
			#pylab.plot(xP,yP, color='g')
	
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d: %d %d %d %d %d %3.2f %3.2f %3.2f %d %3.2f %3.2f %3.2f" % (self.numNodes, isExist1, isExist2, isInterior1, isInterior2, matchCount1, overlapSum1, angleSum1, juncDist1, matchCount2, overlapSum2, angleSum2, juncDist2))
			pylab.savefig("trimDeparture_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1
		
		

		#if isInterior1 and isInterior2:
		#	print "2 departures found"
		#	raise
		
		if len(secP1) == 0:
			print "no departures found"
			raise
		
			
		
		return secP1, secP2


	def getOverlapCondition(self, supportLine, nodeID):

		plotIter = False

		if len(supportLine) == 0:
			return 1e100

		node2 = self.nodeHash[nodeID]

		hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
	
		estPose2 = node2.getGlobalGPACPose()		
			
		#minMatchDist2 = 0.2
		#minMatchDist2 = 0.5
		minMatchDist2 = 1.0
		#minMatchDist2 = 2.0
			
		" set the initial guess "
		poseOrigin = Pose(estPose2)
		
		localSupport = []
		for pnt in supportLine:
			localSupport.append(poseOrigin.convertGlobalToLocal(pnt))
	
		supportSpline = SplineFit(localSupport, smooth=0.1)		
		vecPoints1 = supportSpline.getUniformSamples()
		supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			
		
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2)
				
		support_pairs = []
		#for i in range(len(poly2)):
		for i in range(len(vecPoints2)):
			
			p_2 = vecPoints2[i]
	
			" for every transformed point of A, find it's closest neighbor in B "
			#_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

			try:
				p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
	
				#if gen_icp.isInCircle(p_1, radius2, center2):
	
				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = supportPoints[minI][2]
					#C1 = p_1[2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					#support_pairs.append([points2[i],p_1,C2,C1])
					support_pairs.append([points2[i],supportPoints[minI],C2,C1])
			except:
				pass

		cost = 0.0
		if len(support_pairs) == 0:
			cost = 1e100
		else:
			vals = []
			sum1 = 0.0
			for pair in support_pairs:
		
				a = pair[0]
				b = pair[1]
				Ca = pair[2]
				Cb = pair[3]
		
				ax = a[0]
				ay = a[1]		
				bx = b[0]
				by = b[1]
		
				c11 = Ca[0][0]
				c12 = Ca[0][1]
				c21 = Ca[1][0]
				c22 = Ca[1][1]
						
				b11 = Cb[0][0]
				b12 = Cb[0][1]
				b21 = Cb[1][0]
				b22 = Cb[1][1]	
			
				val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
				
				vals.append(val)
				sum1 += val
				
			cost = sum1 / len(support_pairs)
			

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in supportPoints:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')
	
			xP = []
			yP = []
			for p in poly2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			for pair in support_pairs:
				p1 = pair[0]
				p2 = pair[1]
				xP = [p1[0],p2[0]]
				yP = [p1[1],p2[1]]
				pylab.plot(xP,yP)
					
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("nodeID %d, cost = %f, count = %d" % (nodeID, cost, len(support_pairs)))
			pylab.savefig("overlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 1e100

		return cost


	def restoreNode(self, dirName, newNode):

		self.currNode = newNode
		nodeID = newNode.nodeID
		self.nodeHash[nodeID] = self.currNode
		self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		self.c_hulls.append(computeHull(self.nodeHash[nodeID], static = True))

		f = open(dirName + "/cornerCandidates_%04u.txt" % (nodeID), 'r')		
		newNode.cornerCandidates = eval(f.read().rstrip())
		f.close()


	def correctNode2(self, nodeID):

		"FIXME:  add latest changes from loadNewNode() "

		self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		self.c_hulls.append(computeHull(self.nodeHash[nodeID], static = True))

		newNode = self.nodeHash[nodeID]		

		return self.integrateNode(newNode, nodeID)

						
	def loadNewNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		#self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		#self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		#self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		#self.c_hulls.append(computeHull(self.nodeHash[nodeID], static = True))
		#alphaHull = computeHull(self.nodeHash[nodeID], static = True)
		##longPaths, medialLongPaths, longMedialWidths = self.nodeHash[nodeID].computeFullSkeleton(alphaHull)
		
		
		return self.integrateNode(newNode, nodeID)


	def drawTrimmedPaths(self, trimmedPaths):
		
		pylab.clf()
		
		for k in range(len(trimmedPaths)):
			path = trimmedPaths[k]
			print "path has", len(path), "points"
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color = self.colors[k])

		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		pylab.title("Trimmed Paths, numNodes = %d" % self.numNodes)
		pylab.savefig("trimmedPath_%04u.png" % self.trimCount)
		self.trimCount += 1



	def drawPathAndHull(self):
		
		
		pylab.clf()
		for k in range(len(self.nodeSets)):
			xP = []
			yP = []
			
			for p in self.paths[k]:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=self.colors[k], linewidth=4)


			for nodeID in self.nodeSets[k]:
				xP = []
				yP = []

				estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()		
		
				if self.nodeHash[nodeID].isBowtie:			
					hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
				else:
					hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
		
				" set the origin of pose 1 "
				poseOrigin = Pose(estPose1)
		
				points = []
				for p in hull1:
					p1 = poseOrigin.convertLocalToGlobal(p)
					points.append(p1)
				
				for p in points:
					xP.append(p[0])
					yP.append(p[1])
					
				pylab.plot(xP,yP, color=self.colors[k])

			xP = []
			yP = []
			for p in self.hulls[k]:
				xP.append(p[0])
				yP.append(p[1])
				
			pylab.plot(xP,yP, '--', color=self.colors[k], linewidth=4)
			
		print "pathAndHull:", self.pathDrawCount, self.nodeSets.values()	
		pylab.title("%s" % repr(self.nodeSets.values()))
		pylab.savefig("pathAndHull_%04u.png" % self.pathDrawCount)

		self.pathDrawCount += 1

	def integrateNode(self, newNode, nodeID):

		" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
		direction = newNode.travelDir


		
		if nodeID > 0:
			
			" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
			if nodeID % 2 == 0:

				transform1, covE1, hist1 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )
				if hist1[2] > 0 or hist1[1] > 10 and hist1[2] == 0:
				
					transform2, covE2, hist2 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = not direction )

	
					if hist1[2] < hist2[2]:
						transform = transform1
						covE = covE1
					else:
						if hist1[1] <= hist2[1]:
							transform = transform1
							covE = covE1
						else:
							transform = transform2
							covE = covE2
				else:
					transform = transform1
					covE = covE1
					
				self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)
	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" ODD NUMBER POSES RECEIVE INPLACE CONSTRAINT WITH EVEN PAIR "
				" PERFORM MEDIAL OVERLAP WITH PREVIOUS ODD NUMBER POSE "
			else:
			
				" FIXME:  apply this support check to spliced paths "
				supportLine = self.paths[self.currPath]
				
				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				#self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY2)
				self.inplace1.append([nodeID-1,nodeID,transform,covE])
				transform1 = transform
				offset1 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE1 = covE
				
				if len(supportLine) == 0:
					resultSum1 = 1e100
				else:
					resultSum1 = self.checkSupport(nodeID-1, nodeID, offset1, supportLine)

				node1 = self.nodeHash[nodeID-1]
				node2 = self.nodeHash[nodeID]
				
				" create the ground constraints "
				gndGPAC1Pose = node1.getGndGlobalGPACPose()
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = node2.getGndGlobalGPACPose()
				offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
				transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
				offset2 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE2 = covE
				#resultSum2 = self.checkSupport(nodeID-1, nodeID, offset2, supportLine)

	
				#self.addPriorityEdge([nodeID-1,nodeID,transform,self.E_inplace], INPLACE_PRIORITY3)				
				self.inplace2.append([nodeID-1,nodeID,transform,covE])

				transform, covE, overHist = self.makeMedialOverlapConstraint(nodeID-1, nodeID, isMove = False, inPlace = True, isForward = direction)
				transform3 = transform
				offset3 = [transform[0,0], transform[1,0], transform[2,0]]			
				covE3 = covE
				if len(supportLine) == 0:
					resultSum3 = 1e100
				else:
					resultSum3 = self.checkSupport(nodeID-1, nodeID, offset3, supportLine)
				#self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY)
				self.inplace3.append([nodeID-1,nodeID,transform,covE])

				#print "INPLACE sums:", resultSum1, resultSum3

				if len(supportLine) > 0 and resultSum1 < resultSum3 and resultSum3 - resultSum1 > 100.0:
					self.addPriorityEdge([nodeID-1,nodeID,transform1,covE1], INPLACE_PRIORITY)
				else:
					self.addPriorityEdge([nodeID-1,nodeID,transform3,covE3], INPLACE_PRIORITY)
					


				if nodeID > 2:
					transform1, covE1, hist1 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction)

					if hist1[2] > 0 or hist1[1] > 10 and hist1[2] == 0:
						transform2, covE2, hist2 = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = not direction )

						if hist1[2] < hist2[2]:
							transform = transform1
							covE = covE1
						else:
							if hist1[1] <= hist2[1]:
								transform = transform1
								covE = covE1
							else:
								transform = transform2
								covE = covE2
					else:
						transform = transform1
						covE = covE1
						
					
					self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)

	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()


		" CHECK FOR A BRANCHING EVENT "
		
		if self.numNodes >= 4 and self.numNodes % 2 == 0:


			" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
			for k in range(len(self.nodeSets)):
				print "computing path for node set", k, ":", self.nodeSets[k]
				self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])
					
			self.drawPathAndHull()
			
			" DETECT BRANCHING EVENTS FOR THE NODES OF LAST STOP "
			" AFTER ODD-NUMBER OVERLAP OF CURRENT STOP HAS IRONED OUT ERRORS "
			nodeID1 = self.numNodes-4
			nodeID2 = self.numNodes-3
		
			" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
			if len(self.paths[0]) == 0:
				" first nodes in path 0 "					
				self.nodeSets[0].append(nodeID1)
				self.nodeSets[0].append(nodeID2)

				" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
			else:
				
				" departure events for node 1 "
				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
	
				" departure events for node 2 "
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []

				
				" raw medial axis of union of path nodes "
				paths = {}
				for k in range(len(self.paths)):
					path = self.paths[k]
					if len(path) > 0:
						paths[k] = path
	
				print "paths:", len(paths)
				for k in range(len(paths)):
					print k, len(paths[k])
	
	
				" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
				self.trimmedPaths = self.trimPaths(paths)
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)
	
				
				" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
					
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
				
				print "orderedPathIDs1:", orderedPathIDs1
				print "orderedPathIDs2:", orderedPathIDs2
				
				
				" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
				
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					
				print "node departures", nodeID1, ":", departures1
				print "node  interiors", nodeID1, ":", interiors1
				print "node departures", nodeID2, ":", departures2
				print "node  interiors", nodeID2, ":", interiors2
	
	
				" new junction finding logic "
				" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
				" if a terminal departure exists that is internal, than we have a new junction "
				DISC_THRESH = 0.2

				" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
				frontExist1 = departures1[0][0]
				backExist1 =  departures1[-1][1]
				frontInterior1 = interiors1[0][0]
				backInterior1 = interiors1[-1][1]

				foreTerm1 = frontInterior1 and frontExist1
				backTerm1 = backInterior1 and backExist1
				
				" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
				discForeTerm1 = distances1[0][0] > DISC_THRESH
				discBackTerm1 = distances1[-1][1] > DISC_THRESH
				
				foreAngle1 = depAngles1[0][0]
				backAngle1 = depAngles1[-1][1]
				
				" NODE2: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
				frontExist2 = departures2[0][0]
				backExist2 =  departures2[-1][1]
				frontInterior2 = interiors2[0][0]
				backInterior2 = interiors2[-1][1]

				foreTerm2 = frontInterior2 and frontExist2
				backTerm2 = backInterior2 and backExist2


				discForeTerm2 = distances2[0][0] > DISC_THRESH
				discBackTerm2 = distances2[-1][1] > DISC_THRESH

				foreAngle2 = depAngles2[0][0]
				backAngle2 = depAngles2[-1][1]

				frontAngDiff = diffAngle(foreAngle1, foreAngle2)
				backAngDiff = diffAngle(backAngle1, backAngle2)


				print "pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
				print "pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2

				" check for the cases that internal departure may have a departure point discrepancy, "
				" then we should perform a pose adjustment and recompute departures "
				print "checking discrepancy in departure:", nodeID1, nodeID2, foreTerm1, discForeTerm1, backTerm1, discBackTerm1, foreTerm2, discForeTerm2, backTerm2, discBackTerm2

				#if discForeTerm1 or discBackTerm1 or discForeTerm2 or discBackTerm2:
				if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

					print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2

					#if discForeTerm1 or discBackTerm1:
					if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1:


						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]

						splicedPaths1 = self.splicePathIDs(orderedPathIDs1)

						totalGuesses = []
						for path in splicedPaths1:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
							angDiff = abs(estPose[2]-guessPose1[2])
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID1].setGPACPose(guessPose)


						#if foreTerm1 and discForeTerm1:
						#	pathID = orderedPathIDs1[0]
						#elif backTerm1 and discBackTerm1:
						#	pathID = orderedPathIDs1[-1]
						
							
						#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
						#self.nodeHash[nodeID1].setGPACPose(guessPose1)
							
					#elif discForeTerm2 or discBackTerm2:
					elif foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

						" control point nearest the GPAC origin "
						globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]


						splicedPaths2 = self.splicePathIDs(orderedPathIDs2)

						totalGuesses = []
						for path in splicedPaths2:
		
							" make the path constraints "								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
							
							estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
							angDiff = abs(estPose[2]-guessPose1[2])
							totalGuesses.append((angDiff, cost1, guessPose1))
						
						totalGuesses.sort()
						guessPose = totalGuesses[0][2]
						self.nodeHash[nodeID2].setGPACPose(guessPose)

						#if foreTerm2 and discForeTerm2:
						#	pathID = orderedPathIDs2[0]
						#elif backTerm2 and discBackTerm2:
						#	pathID = orderedPathIDs2[-1]
							
						#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
						#self.nodeHash[nodeID2].setGPACPose(guessPose1)

					departures1 = []
					interiors1 = []
					depPoints1 = []
					distances1 = []
					depAngles1 = []
		
					departures2 = []
					interiors2 = []
					depPoints2 = []
					distances2 = []
					depAngles2 = []
					
					" the overlapping paths are computed from the initial guess of position "
					orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
					orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
					
					print "orderedPathIDs1:", orderedPathIDs1
					print "orderedPathIDs2:", orderedPathIDs2
					
					" now compute whether there are departure points after we have guessed a better position in synch with the paths "
					for pathID in orderedPathIDs1:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
						departures1.append([isExist1,isExist2])
						interiors1.append([isInterior1, isInterior2])
						depPoints1.append([departurePoint1, departurePoint2])
						distances1.append([discDist1, discDist2])
						depAngles1.append([depAngle1, depAngle2])
					
					for pathID in orderedPathIDs2:
						departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
						departures2.append([isExist1,isExist2])
						interiors2.append([isInterior1, isInterior2])
						depPoints2.append([departurePoint1, departurePoint2])
						distances2.append([discDist1, discDist2])
						depAngles2.append([depAngle1, depAngle2])
						
					print "node", nodeID1, ":", departures1, interiors1
					print "node", nodeID2, ":", departures2, interiors2
		
					" new junction finding logic "
					" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
					" if a terminal departure exists that is internal, than we have a new junction "
					
					frontExist1 = departures1[0][0]
					backExist1 =  departures1[-1][1]
					frontInterior1 = interiors1[0][0]
					backInterior1 = interiors1[-1][1]

					foreTerm1 = frontInterior1 and frontExist1
					backTerm1 = backInterior1 and backExist1
					
					
					discForeTerm1 = distances1[0][0] > DISC_THRESH
					discBackTerm1 = distances1[-1][1] > DISC_THRESH
					
					foreAngle1 = depAngles1[0][0]
					backAngle1 = depAngles1[-1][1]
					
					#foreTerm2 = departures2[0][0] and interiors2[0][0]
					#backTerm2 = departures2[-1][1] and interiors2[-1][1]
					frontExist2 = departures2[0][0]
					backExist2 =  departures2[-1][1]
					frontInterior2 = interiors2[0][0]
					backInterior2 = interiors2[-1][1]

					foreTerm2 = frontInterior2 and frontExist2
					backTerm2 = backInterior2 and backExist2

	
					discForeTerm2 = distances2[0][0] > DISC_THRESH
					discBackTerm2 = distances2[-1][1] > DISC_THRESH

					foreAngle2 = depAngles2[0][0]
					backAngle2 = depAngles2[-1][1]
					

					frontAngDiff = diffAngle(foreAngle1, foreAngle2)
					backAngDiff = diffAngle(backAngle1, backAngle2)

					print "pass2:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
					print "pass2:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2
					
				print "frontAngDiff:", frontAngDiff
				print "backAngDiff:", backAngDiff



				newPaths = []
	
				frontExist1 = departures1[0][0]
				frontInterior1 = interiors1[0][0]
		
				frontExist2 = departures2[0][0]
				frontInterior1 = interiors2[0][0]
				
				depAngle1 = depAngles1[0][0]
				depAngle2 = depAngles2[0][0]
				
				depPoint1 = depPoints1[0][0]
				depPoint2 = depPoints2[0][0]
				
				parentPathID1 = orderedPathIDs1[0]
				parentPathID2 = orderedPathIDs2[0]
	
	
				isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)

				print "determineBranch:"
				print isBranch
				print pathBranchIDs
	
				if isBranch[0]:
					orderedPathIDs1.insert(0,pathBranchIDs[0])
					departures1.insert(0, [False, False])
					interiors1.insert(0, [False, False])
					depPoints1.insert(0, [None, None])
	
				if isBranch[1]:
					orderedPathIDs2.insert(0,pathBranchIDs[1])
					departures2.insert(0, [False, False])
					interiors2.insert(0, [False, False])
					depPoints2.insert(0, [None, None])
					
				for pathID in pathBranchIDs:
					if pathID != -1 and newPaths.count(pathID) == 0:
						newPaths.append(pathID)




				backExist1 = departures1[-1][1]
				backInterior1 = interiors1[-1][1]
		
				backExist2 = departures2[-1][1]
				backInterior1 = interiors2[-1][1]
				
				depAngle1 = depAngles1[-1][1]
				depAngle2 = depAngles2[-1][1]
				
				depPoint1 = depPoints1[-1][1]
				depPoint2 = depPoints2[-1][1]
	
				parentPathID1 = orderedPathIDs1[-1]
				parentPathID2 = orderedPathIDs2[-1]
				print "self.pathParents:", self.pathParents
				print "self.nodeSets:", self.nodeSets
				print "newPaths:", newPaths
				print "self.numPaths:", self.numPaths
	
	
				isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)
				print "determineBranch:"
				print isBranch
				print pathBranchIDs
	
				if isBranch[0]:
					orderedPathIDs1.append(0,pathBranchIDs[0])
					departures1.append([False, False])
					interiors1.append([False, False])
					depPoints1.append([None, None])
	
				if isBranch[1]:
					orderedPathIDs2.append(pathBranchIDs[1])
					departures2.append([False, False])
					interiors2.append([False, False])
					depPoints2.append([None, None])
	
				for pathID in pathBranchIDs:
					if pathID != -1 and newPaths.count(pathID) == 0:
						newPaths.append(pathID)
				print "self.pathParents:", self.pathParents
				print "self.nodeSets:", self.nodeSets
				print "newPaths:", newPaths
				print "self.numPaths:", self.numPaths
				




				"""
				newPaths = []
				
				" 60 degree threshold "
				ANG_THRESH = 1.047
				#ANG_THRESH = 0.3
				
				if foreTerm1 or foreTerm2:

					if foreTerm1 and foreTerm2:
						" both departing case: check if same or different "
						if fabs(frontAngDiff) < ANG_THRESH:
							" are on top of each other and are branching to the same path "

							if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]) and self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):
	
								
								pathID = orderedPathIDs1[0]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[0][0]
								depAng = depAngles1[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.insert(0,self.numPaths)
								departures1.insert(0, [False, False])
								interiors1.insert(0, [False, False])
								depPoints1.insert(0, [None, None])
	
								orderedPathIDs2.insert(0,self.numPaths)
								departures2.insert(0, [False, False])
								interiors2.insert(0, [False, False])
								depPoints2.insert(0, [None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
	
								print "foreA"
						else:
							" these are branching to separate new paths "

							if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]) and self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):
	
	
								pathID = orderedPathIDs1[0]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[0][0]
								depAng = depAngles1[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.insert(0,self.numPaths)
								departures1.insert(0, [False, False])
								interiors1.insert(0, [False, False])
								depPoints1.insert(0, [None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
	
	
	
								pathID = orderedPathIDs2[0]
								branchNodeID = nodeID2
								globalJunctionPoint = depPoints2[0][0]
								depAng = depAngles2[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs2.insert(0,self.numPaths)
								departures2.insert(0, [False, False])
								interiors2.insert(0, [False, False])
								depPoints2.insert(0, [None, None])
	
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "foreB"
					
					elif foreTerm1 and not foreTerm2:
						
						if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]):
							
							if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
								" both have same departure "
								pathID = orderedPathIDs1[0]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[0][0]
								depAng = depAngles1[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
								#depAng1 = depAngles1[0][0]
								#depAng2 = depAngles2[0][0]
	
								orderedPathIDs1.insert(0,self.numPaths)
								departures1.insert(0, [False, False])
								interiors1.insert(0, [False, False])
								depPoints1.insert(0, [None, None])
	
								orderedPathIDs2.insert(0,self.numPaths)
								departures2.insert(0, [False, False])
								interiors2.insert(0, [False, False])
								depPoints2.insert(0, [None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "foreC"
								
							else:
								" foreTerm1 has unique departure "	
								pathID = orderedPathIDs1[0]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[0][0]
								depAng = depAngles1[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.insert(0,self.numPaths)
								departures1.insert(0, [False, False])
								interiors1.insert(0, [False, False])
								depPoints1.insert(0, [None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "foreD"

					elif foreTerm2 and not foreTerm1:

						if self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):
	
							if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:
	
								" both have same departure "
								pathID = orderedPathIDs2[0]
								branchNodeID = nodeID2
								globalJunctionPoint = depPoints2[0][0]
								depAng = depAngles2[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.insert(0,self.numPaths)
								departures1.insert(0, [False, False])
								interiors1.insert(0, [False, False])
								depPoints1.insert(0, [None, None])
	
								orderedPathIDs2.insert(0,self.numPaths)
								departures2.insert(0, [False, False])
								interiors2.insert(0, [False, False])
								depPoints2.insert(0, [None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "foreE"
	
							else:
								" foreTerm2 has unique departure "	
								pathID = orderedPathIDs2[0]
								branchNodeID = nodeID2
								globalJunctionPoint = depPoints2[0][0]
								depAng = depAngles2[0][0]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs2.insert(0,self.numPaths)
								departures2.insert(0, [False, False])
								interiors2.insert(0, [False, False])
								depPoints2.insert(0, [None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
	
								print "foreF"
					else:
						print "foreG"
						" no new departure, just add nodes to the leaf paths "
						pathID = orderedPathIDs1[0]
						branchNodeID = nodeID1
						globalJunctionPoint = depPoints1[0][0]
						
						if orderedPathIDs1[0] != orderedPathIDs2[0]:
							print "departing path IDs of two colocated nodes are not the same"
							raise

				if backTerm1 or backTerm2:

					if backTerm1 and backTerm2:
						" both departing case: check if same or different "
						if fabs(backAngDiff) < ANG_THRESH:
							" are on top of each other and are branching to the same path "

							if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]) and self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):
	
								pathID = orderedPathIDs1[-1]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[-1][1]
								depAng = depAngles1[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.append(self.numPaths)
								departures1.append([False, False])
								interiors1.append([False, False])
								depPoints1.append([None, None])
	
								orderedPathIDs2.append(self.numPaths)
								departures2.append([False, False])
								interiors2.append([False, False])
								depPoints2.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "backA"

						else:
							" these are branching to separate new paths "
							if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]) and self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):
	
								pathID = orderedPathIDs1[-1]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[-1][1]
								depAng = depAngles1[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.append(self.numPaths)
								departures1.append([False, False])
								interiors1.append([False, False])
								depPoints1.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
	
								pathID = orderedPathIDs2[-1]
								branchNodeID = nodeID2
								globalJunctionPoint = depPoints2[-1][1]
								depAng = depAngles2[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs2.append(self.numPaths)
								departures2.append([False, False])
								interiors2.append([False, False])
								depPoints2.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "backB"					
					elif backTerm1 and not backTerm2:
						
						if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]):
							
							if backExist2 and fabs(backAngDiff) < ANG_THRESH:
								" both have same departure "
								pathID = orderedPathIDs1[-1]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[-1][1]
								depAng = depAngles1[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.append(self.numPaths)
								departures1.append([False, False])
								interiors1.append([False, False])
								depPoints1.append([None, None])
	
								orderedPathIDs2.append(self.numPaths)
								departures2.append([False, False])
								interiors2.append([False, False])
								depPoints2.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "backC"
							else:
								" backTerm1 has unique departure "	
								pathID = orderedPathIDs1[-1]
								branchNodeID = nodeID1
								globalJunctionPoint = depPoints1[-1][1]
								depAng = depAngles1[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.append(self.numPaths)
								departures1.append([False, False])
								interiors1.append([False, False])
								depPoints1.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "backD"

					elif backTerm2 and not backTerm1:

						if self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):
	
							if backExist1 and fabs(backAngDiff) < ANG_THRESH:
	
								" both have same departure "
								pathID = orderedPathIDs2[-1]
								branchNodeID = nodeID2
								globalJunctionPoint = depPoints2[-1][1]
								depAng = depAngles2[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs1.append(self.numPaths)
								departures1.append([False, False])
								interiors1.append([False, False])
								depPoints1.append([None, None])
	
								orderedPathIDs2.append(self.numPaths)
								departures2.append([False, False])
								interiors2.append([False, False])
								depPoints2.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "backE"
							else:
								" backTerm2 has unique departure "	
								self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1])
	
								pathID = orderedPathIDs2[-1]
								branchNodeID = nodeID2
								globalJunctionPoint = depPoints2[-1][1]
								depAng = depAngles2[-1][1]
								junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
								orderedPathIDs2.append(self.numPaths)
								departures2.append([False, False])
								interiors2.append([False, False])
								depPoints2.append([None, None])
	
								poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
								self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
								self.nodeSets[self.numPaths] = []
								newPaths.append(self.numPaths)
								self.numPaths += 1
								print "backF"
					else:
						print "backG"
						" no new departure, just add nodes to the leaf paths "
						pathID = orderedPathIDs1[-1]
						branchNodeID = nodeID1
						globalJunctionPoint = depPoints1[-1][1]
						if orderedPathIDs1[-1] != orderedPathIDs2[-1]:
							print "departing path IDs of two colocated nodes are not the same"
							raise


				"""



				" determine which paths are leaves "
				print "self.pathParents:", self.pathParents
				isAParent = [False for k in range(self.numPaths)]
				for k in orderedPathIDs1:
					print "index:", k
					currPath = self.pathParents[k]
					currParent = currPath[0]
					if currParent != None:
						isAParent[currParent] = True
				
				"add nodes to paths that are the leaves "
				for pathID in orderedPathIDs1:
					if not isAParent[pathID]:				
						self.nodeSets[pathID].append(nodeID1)

				isAParent = [False for k in range(self.numPaths)]
				for k in orderedPathIDs2:
					print "index:", k
					currPath = self.pathParents[k]
					currParent = currPath[0]
					if currParent != None:
						isAParent[currParent] = True

				for pathID in orderedPathIDs2:
					if not isAParent[pathID]:				
						self.nodeSets[pathID].append(nodeID2)

				paths = {}
				for k in range(len(self.paths)):
					path = self.paths[k]
					if len(path) > 0:
						paths[k] = path	
				self.trimmedPaths = self.trimPaths(paths)

				" perform path constraints only if we have not created a new path "
				isNew1 = False
				isNew2 = False
				for pathID in newPaths:
					if orderedPathIDs1.count(pathID) > 0:
						isNew1 = True
					if orderedPathIDs2.count(pathID) > 0:
						isNew2 = True
				
	
				
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)



				"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
				"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
				"3)  select choice with the lowest cost "
				"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
				
				" nodeID1:  is it in a junction or a single path? "
				
				
				if not isNew1:
					self.constrainToPaths(nodeID1, orderedPathIDs1, departures1, interiors1, depPoints1)
					
				if not isNew2:
					self.constrainToPaths(nodeID2, orderedPathIDs2, departures2, interiors2, depPoints2)

				self.mergePriorityConstraints()

			for k in range(len(self.nodeSets)):
				print "computing path for node set", k, ":", self.nodeSets[k]
				self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])

			self.drawPathAndHull()
			
			return


			if False:				
				if not isNew1:
					splicedPaths1 = self.splicePathIDs(orderedPathIDs1)


					"""
					pylab.clf()
	
					print len(splicedPaths1), "spliced paths 1"
					#print splicedPaths1
	
					for path in splicedPaths1:
						print "spliced path has", len(path), "points"
						xP = []
						yP = []
						for p in path:
							xP.append(p[0])
							yP.append(p[1])
			
						pylab.plot(xP,yP)
	
					pylab.xlim(-4,4)
					pylab.ylim(-4,4)
			
					pylab.title("Spliced Paths 1, numNodes = %d" % self.numNodes)
					pylab.savefig("splicedPath_%04u.png" % self.spliceCount)
					self.spliceCount += 1
					"""
					if len(orderedPathIDs1) == 1:
	
						" in a single path "
	
						pathID = orderedPathIDs1[0]
						nodeSet = self.nodeSets[pathID]
						
						" check that we're not the only pose in the path "
						" does it contain at least one node that is not nodeID1 or nodeID2 "
						doConstraint = False
						if len(nodeSet) == 1:						
							if nodeSet.count(nodeID1) == 0:
								doConstraint = True
						
						elif len(nodeSet) == 2:						
							if nodeSet.count(nodeID1) == 0 or nodeSet.count(nodeID2) == 0:
								doConstraint = True
						
						elif len(nodeSet) > 2:
							doConstraint = True
	
						if doConstraint:
							
							" control point nearest the GPAC origin "
							globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]
								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
							self.nodeHash[nodeID1].setGPACPose(guessPose1)
	
							" make path constraint "
							print "pathID: addPathConstraints(", pathID, nodeID1
							self.addPathConstraints(self.nodeSets[pathID], nodeID1)
						
					else:
						" in at least one junction "
						
						" we don't intersect the junction sufficiently "
						" add constraints to most overlapped path"
						#if interiors1.count(True) == 0:
						" no recognized internal departure point "
						if [item for inner_list in interiors1 for item in inner_list].count(True) == 0:

		
							" determine which paths are leaves "
							isAChild = [True for k in range(self.numPaths)]
							for k in orderedPathIDs1:
								currPath = self.pathParents[k]
								currParent = currPath[0]
								if currParent != None and orderedPathIDs1.count(currParent) > 0:
									isAChild[currParent] = False

							print "isAChild:", isAChild

							minPathID = 0
							minSum = 1e100
							for pathID in orderedPathIDs1:
								
								if isAChild[pathID]:
									sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID1)
									
									if sum1 < minSum:
										minPathID = pathID
										minSum = sum1
									print "overlap sum", pathID, "=", sum1						

							print "maximum overlap path is", minPathID

							nodeSet = self.nodeSets[minPathID]
							
							" check that we're not the only pose in the path "
							" does it contain at least one node that is not nodeID1 or nodeID2 "
							doConstraint = False
							if len(nodeSet) == 1:						
								if nodeSet.count(nodeID1) == 0:
									doConstraint = True
							
							elif len(nodeSet) == 2:						
								if nodeSet.count(nodeID1) == 0 or nodeSet.count(nodeID2) == 0:
									doConstraint = True
							
							elif len(nodeSet) > 2:
								doConstraint = True
		
							if doConstraint:
								
								" control point nearest the GPAC origin "
								globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]

								totalGuesses = []
								
								for path in splicedPaths1:
		
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
									
									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose1 = totalGuesses[0][2]
								
									
								#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, self.trimmedPaths[minPathID], localPoint1, localPoint1)
								self.nodeHash[nodeID1].setGPACPose(guessPose1)
		
								" make path constraint "
								print "minPathID: addPathConstraints(", minPathID, nodeID1
								self.addPathConstraints(self.nodeSets[minPathID], nodeID1)






						
						elif len(orderedPathIDs1) == 2:
							
							" one junction point "
							pathID1 = orderedPathIDs1[0]
							pathID2 = orderedPathIDs1[1]
	
							#[pathID, branchNodeID, poseOrigin.convertGlobalToLocal(globalJunctionPoint)]
							
							" find parent-child relationship "
							val1 = self.pathParents[pathID1]
							val2 = self.pathParents[pathID2]
	
	
							if val1[0] == None:
								parentPath = pathID1
								childPath = pathID2
								
							elif val2[0] == None:					
								parentPath = pathID2
								childPath = pathID1
								
							else:
								parentID1 = val1[0]
								parentID2 = val2[0]
								
								if parentID1 == pathID2:
									parentPath = pathID2
									childPath = pathID1
								elif parentID2 == pathID1:
									parentPath = pathID1
									childPath = pathID2
								else:
									raise
	
							" global path point "
							localJunctionNodeID = self.pathParents[childPath][1]
							localPathJunctionPoint = self.pathParents[childPath][2]						
							
							localJunctionNode = self.nodeHash[localJunctionNodeID]
							poseOrigin = Pose(localJunctionNode.getEstPose())
							globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
						
							print departures1, interiors1
							print departures2, interiors2
							print pathID1, childPath, parentPath
						
							" new node departure point "
							" try childPath first "
							if childPath == pathID1:
								index1 = orderedPathIDs1.index(pathID1)			
								isDepExists1 = departures1[index1][1] and interiors1[index1][1]
								depPoint1 = depPoints1[index1][1]
								print "caseA:", index1, isDepExists1
							
							else:
								index2 = orderedPathIDs1.index(pathID2)
								isDepExists1 = departures1[index2][0] and interiors1[index2][0]
								depPoint1 = depPoints1[index2][0]
								print "caseB:", index2, isDepExists1
	
							" try parentPath first "
							if parentPath == pathID1:
								index1 = orderedPathIDs1.index(pathID1)
								isDepExists2 = departures1[index1][1] and interiors1[index1][1]
								depPoint2 = depPoints1[index1][1]
								print "caseC:", index1, isDepExists2
							else:
								index2 = orderedPathIDs1.index(pathID2)
								isDepExists2 = departures1[index2][0] and interiors1[index2][0]
								depPoint2 = depPoints1[index2][0]
								print "caseD:", index2, isDepExists2
							
							
							if isDepExists1:
								
								totalGuesses = []
								
								for path in splicedPaths1:
		
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPathJunctionPoint, globalPathJunctionPoint)
									
									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose = totalGuesses[0][2]
								
								self.nodeHash[nodeID1].setGPACPose(guessPose)
		
								" make path constraint "
								print "childPath: addPathConstraints(", childPath, nodeID1
								self.addPathConstraints(self.nodeSets[childPath], nodeID1)
	
							
							elif isDepExists2:
								totalGuesses = []
								
								for path in splicedPaths1:
	
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPathJunctionPoint, depPoint2)

									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose = totalGuesses[0][2]
	
								self.nodeHash[nodeID1].setGPACPose(guessPose)
		
								" make path constraint "
								print "childPath: addPathConstraints(", childPath, nodeID1
								self.addPathConstraints(self.nodeSets[childPath], nodeID1)
	
							
							else:
								print "no local junction point "
								raise
							
							#nodeID1
							#departures1.append([isExist1,isExist2])
							#interiors1.append([isInterior1, isInterior2])
							#depPoints1.append([departurePoint1, departurePoint2])
							
						
						else:
							" two or more junction points "
							pathID1_1 = orderedPathIDs1[0]
							pathID1_2 = orderedPathIDs1[1]
							pathID2_1 = orderedPathIDs1[-1]
							pathID2_2 = orderedPathIDs1[-2]
	
							" pick the newest path from where to sample the local junction "
							if pathID1_1 > pathID2_2:
								pathID1 = pathID1_1
								pathID2 = pathID1_2
							else:
								pathID1 = pathID2_1
								pathID2 = pathID2_2
								
							" find parent-child relationship "
							val1 = self.pathParents[pathID1]
							val2 = self.pathParents[pathID2]
	
	
							if val1[0] == None:
								parentPath = pathID1
								childPath = pathID2
	
							elif val2[0] == None:					
								parentPath = pathID2
								childPath = pathID1
	
							else:
								parentID1 = val1[0]
								parentID2 = val2[0]
								
								if parentID1 == pathID2:
									parentPath = pathID2
									childPath = pathID1
								elif parentID2 == pathID1:
									parentPath = pathID1
									childPath = pathID2
								else:
									raise
	
							if self.pathParents[childPath][0] != parentPath:
								print "WARNING: selected childPath", childPath, "does not have parent", parentPath, "but has instead", self.pathParents[childPath][0]
	
							" global path point "
							" we don't check, but we assume the parent here is the parent in our orderedPaths"
							localJunctionNodeID = self.pathParents[childPath][1]
							localPathJunctionPoint = self.pathParents[childPath][2]						
							
							localJunctionNode = self.nodeHash[localJunctionNodeID]
							poseOrigin = Pose(localJunctionNode.getEstPose())
							globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
						
							totalGuesses = []
							
							for path in splicedPaths1:
	
								" make the path constraints "								
								guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPathJunctionPoint, globalPathJunctionPoint)
								estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
								angDiff = abs(estPose[2]-guessPose1[2])
								totalGuesses.append((angDiff, cost1, guessPose1))
							
							totalGuesses.sort()
							
							guessPose = totalGuesses[0][2]
							
							self.nodeHash[nodeID1].setGPACPose(guessPose)
	
							" make path constraint "
							print "childPath: addPathConstraints(", childPath, nodeID1
							self.addPathConstraints(self.nodeSets[childPath], nodeID1)
	
							
					
				if not isNew2:

					splicedPaths2 = self.splicePathIDs(orderedPathIDs2)

					"""
					pylab.clf()
	
					print len(splicedPaths2), "spliced paths 2"
					#print splicedPaths2
	
					for path in splicedPaths2:
						print "spliced path has", len(path), "points"
						xP = []
						yP = []
						for p in path:
							xP.append(p[0])
							yP.append(p[1])
			
						pylab.plot(xP,yP)
	
					pylab.xlim(-4,4)
					pylab.ylim(-4,4)
			
					pylab.title("Spliced Paths 2, numNodes = %d" % self.numNodes)
					pylab.savefig("splicedPath_%04u.png" % self.spliceCount)
					self.spliceCount += 1
					"""
	
					" nodeID2:  is it in a junction or a single path? "
					if len(orderedPathIDs2) == 1:
	
						" in a single path "
	
						pathID = orderedPathIDs2[0]
						nodeSet = self.nodeSets[pathID]
						
						" check that we're not the only pose in the path "
						" does it contain at least one node that is not nodeID1 or nodeID2 "
						doConstraint = False
						if len(nodeSet) == 1:						
							if nodeSet.count(nodeID2) == 0:
								doConstraint = True
						
						elif len(nodeSet) == 2:						
							if nodeSet.count(nodeID1) == 0 or nodeSet.count(nodeID2) == 0:
								doConstraint = True
						
						elif len(nodeSet) > 2:
							doConstraint = True
	
						if doConstraint:
							
							" control point nearest the GPAC origin "
							localPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]
								
							guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[pathID], localPoint1, localPoint1)
							self.nodeHash[nodeID2].setGPACPose(guessPose1)
	
							" make path constraint "
							print "pathID: addPathConstraints(", pathID, nodeID2
							self.addPathConstraints(self.nodeSets[pathID], nodeID2)
						
					else:
						" in at least one junction "

						" we don't intersect the junction sufficiently, only constraint to most overlapped path"
						
						if [item for inner_list in interiors2 for item in inner_list].count(True) == 0:

							" determine which paths are leaves "
							isAChild = [True for k in range(self.numPaths)]
							for k in orderedPathIDs2:
								currPath = self.pathParents[k]
								currParent = currPath[0]
								if currParent != None and orderedPathIDs2.count(currParent) > 0:
									isAChild[currParent] = False

							print "isAChild:", isAChild

							minPathID = 0
							minSum = 1e100
							for pathID in orderedPathIDs2:
								
								if isAChild[pathID]:
									sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID2)
									
									if sum1 < minSum:
										minPathID = pathID
										minSum = sum1
									print "overlap sum", pathID, "=", sum1						

							print "maximum overlap path is", minPathID

							nodeSet = self.nodeSets[minPathID]
							
							" check that we're not the only pose in the path "
							" does it contain at least one node that is not nodeID1 or nodeID2 "
							doConstraint = False
							if len(nodeSet) == 1:						
								if nodeSet.count(nodeID2) == 0:
									doConstraint = True
							
							elif len(nodeSet) == 2:						
								if nodeSet.count(nodeID1) == 0 or nodeSet.count(nodeID2) == 0:
									doConstraint = True
							
							elif len(nodeSet) > 2:
								doConstraint = True
		
							if doConstraint:
								
								" control point nearest the GPAC origin "
								globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]

								totalGuesses = []
								
								for path in splicedPaths2:
		
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
									
									estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose1 = totalGuesses[0][2]
								
									
								#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[minPathID], localPoint1, localPoint1)
								self.nodeHash[nodeID2].setGPACPose(guessPose1)
		
								" make path constraint "
								print "minPathID: addPathConstraints(", minPathID, nodeID2
								self.addPathConstraints(self.nodeSets[minPathID], nodeID2)



						
						elif len(orderedPathIDs2) == 2:
							
							" one junction point "
							pathID1 = orderedPathIDs2[0]
							pathID2 = orderedPathIDs2[1]
	
							#[pathID, branchNodeID, poseOrigin.convertGlobalToLocal(globalJunctionPoint)]
							
							" find parent-child relationship "
							val1 = self.pathParents[pathID1]
							val2 = self.pathParents[pathID2]
	
	
							if val1[0] == None:
								parentPath = pathID1
								childPath = pathID2
								
							elif val2[0] == None:					
								parentPath = pathID2
								childPath = pathID1
								
							else:
								parentID1 = val1[0]
								parentID2 = val2[0]
								
								if parentID1 == pathID2:
									parentPath = pathID2
									childPath = pathID1
								elif parentID2 == pathID1:
									parentPath = pathID1
									childPath = pathID2
								else:
									raise
	
							" global path point "
							localJunctionNodeID = self.pathParents[childPath][1]
							localPathJunctionPoint = self.pathParents[childPath][2]						
							
							localJunctionNode = self.nodeHash[localJunctionNodeID]
							poseOrigin = Pose(localJunctionNode.getEstPose())
							globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
						
							print departures1, interiors1
							print departures2, interiors2
							print pathID1, childPath, parentPath
						
							" new node departure point "
							" try childPath first "
							if childPath == pathID1:
								index1 = orderedPathIDs2.index(pathID1)
								isDepExists1 = departures2[index1][1] and interiors2[index1][1]
								depPoint1 = depPoints2[index1][1]
								print "caseE:", index1, isDepExists1
							
							else:
								index2 = orderedPathIDs2.index(pathID2)
								isDepExists1 = departures2[index2][0] and interiors2[index2][0]
								depPoint1 = depPoints2[index2][0]
								print "caseF:", index2, isDepExists1
	
							" try parentPath first "
							if parentPath == pathID1:
								index1 = orderedPathIDs2.index(pathID1)
								isDepExists2 = departures2[index1][1] and interiors2[index1][1]
								depPoint2 = depPoints2[index1][1]
								print "caseG:", index1, isDepExists2
							else:
								index2 = orderedPathIDs2.index(pathID2)
								isDepExists2 = departures2[index2][0] and interiors2[index2][0]
								depPoint2 = depPoints2[index2][0]
								print "caseH:", index2, isDepExists2
							
							
							if isDepExists1:
								
								totalGuesses = []
								
								for path in splicedPaths2:
		
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPathJunctionPoint, depPoint1)
									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
									#totalGuesses.append((cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose = totalGuesses[0][2]
								
								self.nodeHash[nodeID2].setGPACPose(guessPose)
		
								" make path constraint "
								print "childPath: addPathConstraints(", childPath, nodeID2
								self.addPathConstraints(self.nodeSets[childPath], nodeID2)
	
							
							elif isDepExists2:
								totalGuesses = []
								
								for path in splicedPaths2:
	
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPathJunctionPoint, depPoint2)
									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
									#totalGuesses.append((cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose = totalGuesses[0][2]
	
								self.nodeHash[nodeID2].setGPACPose(guessPose)
		
								" make path constraint "
								print "childPath: addPathConstraints(", childPath, nodeID2
								self.addPathConstraints(self.nodeSets[childPath], nodeID2)
	
							
							else:
								print "no local junction point "
								raise
							
							#nodeID1
							#departures1.append([isExist1,isExist2])
							#interiors1.append([isInterior1, isInterior2])
							#depPoints1.append([departurePoint1, departurePoint2])

						else:
							" two or more junction points "
							pathID1_1 = orderedPathIDs2[0]
							pathID1_2 = orderedPathIDs2[1]
							pathID2_1 = orderedPathIDs2[-1]
							pathID2_2 = orderedPathIDs2[-2]
	
							" pick the newest path from where to sample the local junction "
							if pathID1_1 > pathID2_2:
								pathID1 = pathID1_1
								pathID2 = pathID1_2
							else:
								pathID1 = pathID2_1
								pathID2 = pathID2_2
								
							" find parent-child relationship "
							val1 = self.pathParents[pathID1]
							val2 = self.pathParents[pathID2]
	
	
							if val1[0] == None:
								parentPath = pathID1
								childPath = pathID2
	
							elif val2[0] == None:					
								parentPath = pathID2
								childPath = pathID1
	
							else:
								parentID1 = val1[0]
								parentID2 = val2[0]
								
								if parentID1 == pathID2:
									parentPath = pathID2
									childPath = pathID1
								elif parentID2 == pathID1:
									parentPath = pathID1
									childPath = pathID2
								else:
									raise
	
							if self.pathParents[childPath][0] != parentPath:
								print "WARNING: selected childPath", childPath, "does not have parent", parentPath, "but has instead", self.pathParents[childPath][0]
	
							" global path point "
							" we don't check, but we assume the parent here is the parent in our orderedPaths"
							localJunctionNodeID = self.pathParents[childPath][1]
							localPathJunctionPoint = self.pathParents[childPath][2]						
							
							localJunctionNode = self.nodeHash[localJunctionNodeID]
							poseOrigin = Pose(localJunctionNode.getEstPose())
							globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
						
							totalGuesses = []
							
							for path in splicedPaths2:
	
								" make the path constraints "								
								guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPathJunctionPoint, globalPathJunctionPoint)
								estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
								angDiff = abs(estPose[2]-guessPose1[2])
								totalGuesses.append((angDiff, cost1, guessPose1))
							
							totalGuesses.sort()
							
							guessPose = totalGuesses[0][2]
							
							self.nodeHash[nodeID2].setGPACPose(guessPose)
	
							" make path constraint "
							print "childPath: addPathConstraints(", childPath, nodeID2
							self.addPathConstraints(self.nodeSets[childPath], nodeID2)
								
						"""
						else:
							" two junction points "
							"FIXME:  take the first only, select it more intelligently"
							pathID1_1 = orderedPathIDs2[0]
							pathID1_2 = orderedPathIDs2[1]
							pathID2_1 = orderedPathIDs2[-1]
							pathID2_2 = orderedPathIDs2[-2]
	
	
							" find parent-child relationship "
							val1 = self.pathParents[pathID1_1]
							val2 = self.pathParents[pathID1_2]
	
	
							if val1[0] == None:
								parentPath = pathID1_1
								childPath = pathID1_2
	
							elif val2[0] == None:					
								parentPath = pathID1_2
								childPath = pathID1_1
	
							else:
								parentID1 = val1[0]
								parentID2 = val2[0]
								
								if parentID1 == pathID1_2:
									parentPath = pathID1_2
									childPath = pathID1_1
								elif parentID2 == pathID1_1:
									parentPath = pathID1_1
									childPath = pathID1_2
								else:
									raise
	
							" global path point "
							localJunctionNodeID = self.pathParents[childPath][1]
							localPathJunctionPoint = self.pathParents[childPath][2]						
							
							localJunctionNode = self.nodeHash[localJunctionNodeID]
							poseOrigin = Pose(localJunctionNode.getEstPose())
							globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
						
							" new node departure point "
							" try childPath first "
							if childPath == pathID1_1:
								index1_1 = orderedPathIDs2.index(pathID1_1)
								isDepExists1 = departures2[index1_1][1] and interiors2[index1_1][1]
								depPoint1 = depPoints2[index1_1][1]
							else:
								index1_2 = orderedPathIDs2.index(pathID1_2)
								isDepExists1 = departures2[index1_2][0] and interiors2[index1_2][0]
								depPoint1 = depPoints2[index1_2][0]
	
							" try parentPath first "
							if parentPath == pathID1_1:
								index1_1 = orderedPathIDs2.index(pathID1_1)
								isDepExists2 = departures2[index1_1][1] and interiors2[index1_1][1]
								depPoint2 = depPoints2[index1_1][1]
							else:
								index1_2 = orderedPathIDs2.index(pathID1_2)
								isDepExists2 = departures2[index1_2][0] and interiors2[index1_2][0]
								depPoint2 = depPoints2[index1_2][0]
								
							
							if isDepExists1:
								totalGuesses = []
								
								for path in splicedPaths2:
		
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPathJunctionPoint, depPoint1)
									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
									#totalGuesses.append((cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose = totalGuesses[0][2]
								
								self.nodeHash[nodeID2].setGPACPose(guessPose)
		
								" make path constraint "
								print "childPath: addPathConstraints(", childPath, nodeID2
								self.addPathConstraints(self.nodeSets[childPath], nodeID2)
	
							
							elif isDepExists2:
								totalGuesses = []
								
								for path in splicedPaths2:
	
									" make the path constraints "								
									guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPathJunctionPoint, depPoint2)
									estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
									angDiff = abs(estPose[2]-guessPose1[2])
									totalGuesses.append((angDiff, cost1, guessPose1))
									#totalGuesses.append((cost1, guessPose1))
								
								totalGuesses.sort()
								
								guessPose = totalGuesses[0][2]
	
								self.nodeHash[nodeID2].setGPACPose(guessPose)
		
								" make path constraint "
								print "childPath: addPathConstraints(", childPath, nodeID2
								self.addPathConstraints(self.nodeSets[childPath], nodeID2)
	
							
							else:
								print "no local junction point "
								raise
						"""	
	

				self.mergePriorityConstraints()


				#newPaths1 = self.splicePathIDs(orderedPathIDs1)
				#newPaths2 = self.splicePathIDs(orderedPathIDs2)
				#for k in range(len(newPaths1)):
				#	guessPose, cost = self.makeGlobalMedialOverlapConstraint(nodeID1, newPaths1[k], self.junctionPoint, departurePoint0)
			
				
				#self.paths[self.numPaths], self.hulls[self.numPaths] = self.getTopology(self.nodeSets[self.numPaths])
			
				#pass




			for k in range(len(self.nodeSets)):
				print "computing path for node set", k, ":", self.nodeSets[k]
				self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])

			self.drawPathAndHull()
			
			return

	def addNewPath(self, parentID, branchNodeID, localJunctionPose):
				
		" FIXME:  Branch direction does not seem to be used properly yet.  Only used in getOverlapDeparture() "
		self.pathParents.append([parentID, branchNodeID, localJunctionPose, True])

		newPathID = self.numPaths
		self.nodeSets[newPathID] = []
		self.numPaths += 1

		return newPathID
			
	def constrainToPaths(self, nodeID, orderedPathIDs, departures, interiors, depPoints):
			
		"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
		"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
		"3)  select choice with the lowest cost "
		"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
		
		" nodeID1:  is it in a junction or a single path? "
		
		splicedPaths1 = self.splicePathIDs(orderedPathIDs)

		if len(orderedPathIDs) == 1:

			" in a single path "

			pathID = orderedPathIDs[0]
			nodeSet = self.nodeSets[pathID]
			
			" check that we're not the only pose in the path "
			" does it contain at least one node that is not nodeID1 or nodeID2 "
			doConstraint = False
			
			" companion node iD "
			if nodeID % 2 == 0:
				pairNodeID = nodeID + 1
			else:
				pairNodeID = nodeID - 1
				 
			if len(nodeSet) == 1:						
				if nodeSet.count(nodeID) == 0:
					doConstraint = True
			
			elif len(nodeSet) == 2:						
				if nodeSet.count(nodeID) == 0 or nodeSet.count(pairNodeID) == 0:
					doConstraint = True
			
			elif len(nodeSet) > 2:
				doConstraint = True

			if doConstraint:
				
				" control point nearest the GPAC origin "
				globalPoint1 = self.nodeHash[nodeID].getGlobalGPACPose()[:2]
					
				guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
				self.nodeHash[nodeID].setGPACPose(guessPose1)

				" make path constraint "
				print "pathID: addPathConstraints(", pathID, nodeID
				self.addPathConstraints(self.nodeSets[pathID], nodeID)
			
		else:
			" in at least one junction "
			
			" we don't intersect the junction sufficiently "
			" add constraints to most overlapped path"
			#if interiors.count(True) == 0:
			" no recognized internal departure point "
			if [item for inner_list in interiors for item in inner_list].count(True) == 0:


				" determine which paths are leaves "
				isAChild = [True for k in range(self.numPaths)]
				for k in orderedPathIDs:
					currPath = self.pathParents[k]
					currParent = currPath[0]
					if currParent != None and orderedPathIDs.count(currParent) > 0:
						isAChild[currParent] = False

				print "isAChild:", isAChild

				minPathID = 0
				minSum = 1e100
				for pathID in orderedPathIDs:
					
					if isAChild[pathID]:
						sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
						
						if sum1 < minSum:
							minPathID = pathID
							minSum = sum1
						print "overlap sum", pathID, "=", sum1						

				print "maximum overlap path is", minPathID

				nodeSet = self.nodeSets[minPathID]
				
				" check that we're not the only pose in the path "
				" does it contain at least one node that is not nodeID or nodeID2 "
				doConstraint = False

				" companion node iD "
				if nodeID % 2 == 0:
					pairNodeID = nodeID + 1
				else:
					pairNodeID = nodeID - 1

				if len(nodeSet) == 1:						
					if nodeSet.count(nodeID) == 0:
						doConstraint = True
				
				elif len(nodeSet) == 2:						
					if nodeSet.count(nodeID) == 0 or nodeSet.count(pairNodeID) == 0:
						doConstraint = True
				
				elif len(nodeSet) > 2:
					doConstraint = True

				if doConstraint:
					
					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID].getGlobalGPACPose()[:2]

					totalGuesses = []
					
					for path in splicedPaths1:

						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					
					guessPose1 = totalGuesses[0][2]
					
						
					#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, self.trimmedPaths[minPathID], localPoint1, localPoint1)
					self.nodeHash[nodeID].setGPACPose(guessPose1)

					" make path constraint "
					print "minPathID: addPathConstraints(", minPathID, nodeID
					self.addPathConstraints(self.nodeSets[minPathID], nodeID)






			
			elif len(orderedPathIDs) == 2:
				
				" one junction point "
				pathID1 = orderedPathIDs[0]
				pathID2 = orderedPathIDs[1]

				#[pathID, branchNodeID, poseOrigin.convertGlobalToLocal(globalJunctionPoint)]
				
				" find parent-child relationship "
				val1 = self.pathParents[pathID1]
				val2 = self.pathParents[pathID2]


				if val1[0] == None:
					parentPath = pathID1
					childPath = pathID2
					
				elif val2[0] == None:					
					parentPath = pathID2
					childPath = pathID1
					
				else:
					parentID1 = val1[0]
					parentID2 = val2[0]
					
					if parentID1 == pathID2:
						parentPath = pathID2
						childPath = pathID1
					elif parentID2 == pathID1:
						parentPath = pathID1
						childPath = pathID2
					else:
						raise

				" global path point "
				localJunctionNodeID = self.pathParents[childPath][1]
				localPathJunctionPoint = self.pathParents[childPath][2]						
				
				localJunctionNode = self.nodeHash[localJunctionNodeID]
				poseOrigin = Pose(localJunctionNode.getEstPose())
				globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
			
				print departures, interiors
				print pathID1, childPath, parentPath
			
				" new node departure point "
				" try childPath first "
				if childPath == pathID1:
					index1 = orderedPathIDs.index(pathID1)			
					isDepExists1 = departures[index1][1] and interiors[index1][1]
					depPoint1 = depPoints[index1][1]
					print "caseA:", index1, isDepExists1
				
				else:
					index2 = orderedPathIDs.index(pathID2)
					isDepExists1 = departures[index2][0] and interiors[index2][0]
					depPoint1 = depPoints[index2][0]
					print "caseB:", index2, isDepExists1

				" try parentPath first "
				if parentPath == pathID1:
					index1 = orderedPathIDs.index(pathID1)
					isDepExists2 = departures[index1][1] and interiors[index1][1]
					depPoint2 = depPoints[index1][1]
					print "caseC:", index1, isDepExists2
				else:
					index2 = orderedPathIDs.index(pathID2)
					isDepExists2 = departures[index2][0] and interiors[index2][0]
					depPoint2 = depPoints[index2][0]
					print "caseD:", index2, isDepExists2
				
				
				if isDepExists1:
					
					totalGuesses = []
					
					for path in splicedPaths1:

						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPathJunctionPoint, globalPathJunctionPoint)
						
						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					
					guessPose = totalGuesses[0][2]
					
					self.nodeHash[nodeID].setGPACPose(guessPose)

					" make path constraint "
					print "childPath: addPathConstraints(", childPath, nodeID
					self.addPathConstraints(self.nodeSets[childPath], nodeID)

				
				elif isDepExists2:
					totalGuesses = []
					
					for path in splicedPaths1:

						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPathJunctionPoint, depPoint2)

						estPose = self.nodeHash[nodeID].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					
					guessPose = totalGuesses[0][2]

					self.nodeHash[nodeID].setGPACPose(guessPose)

					" make path constraint "
					print "childPath: addPathConstraints(", childPath, nodeID
					self.addPathConstraints(self.nodeSets[childPath], nodeID)

				
				else:
					print "no local junction point "
					raise
		
			else:
				" two or more junction points "
				pathID1_1 = orderedPathIDs[0]
				pathID1_2 = orderedPathIDs[1]
				pathID2_1 = orderedPathIDs[-1]
				pathID2_2 = orderedPathIDs[-2]

				" pick the newest path from where to sample the local junction "
				if pathID1_1 > pathID2_2:
					pathID1 = pathID1_1
					pathID2 = pathID1_2
				else:
					pathID1 = pathID2_1
					pathID2 = pathID2_2
					
				" find parent-child relationship "
				val1 = self.pathParents[pathID1]
				val2 = self.pathParents[pathID2]


				if val1[0] == None:
					parentPath = pathID1
					childPath = pathID2

				elif val2[0] == None:					
					parentPath = pathID2
					childPath = pathID1

				else:
					parentID1 = val1[0]
					parentID2 = val2[0]
					
					if parentID1 == pathID2:
						parentPath = pathID2
						childPath = pathID1
					elif parentID2 == pathID1:
						parentPath = pathID1
						childPath = pathID2
					else:
						raise

				if self.pathParents[childPath][0] != parentPath:
					print "WARNING: selected childPath", childPath, "does not have parent", parentPath, "but has instead", self.pathParents[childPath][0]

				" global path point "
				" we don't check, but we assume the parent here is the parent in our orderedPaths"
				localJunctionNodeID = self.pathParents[childPath][1]
				localPathJunctionPoint = self.pathParents[childPath][2]						
				
				localJunctionNode = self.nodeHash[localJunctionNodeID]
				poseOrigin = Pose(localJunctionNode.getEstPose())
				globalPathJunctionPoint = poseOrigin.convertLocalToGlobal(localPathJunctionPoint)
			
				totalGuesses = []
				
				for path in splicedPaths1:

					" make the path constraints "								
					guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID, path, globalPathJunctionPoint, globalPathJunctionPoint)
					estPose = self.nodeHash[nodeID].getGlobalGPACPose()
					angDiff = abs(estPose[2]-guessPose1[2])
					totalGuesses.append((angDiff, cost1, guessPose1))
				
				totalGuesses.sort()
				
				guessPose = totalGuesses[0][2]
				
				self.nodeHash[nodeID].setGPACPose(guessPose)

				" make path constraint "
				print "childPath: addPathConstraints(", childPath, nodeID
				self.addPathConstraints(self.nodeSets[childPath], nodeID)

	def determineBranch(self, nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2):
		
		"""
		1)  both nodes departing
			- are same departure  ( equal departing angle and are on top of each other)  1 new path
			- are different departures ( different departing angle and are not necessarily on top of each other) 2 new paths
			
		2) neither node is departing
			- are on same path
			
		3) one node is departing only
			- is only one departure ( other node has no departure, or has external departure but not same angle )
			- is a false departure  ( falseness happens if bad map data, ignore )
			- are both departing ( should have at least external departure on both, compare departure angle )
		
		Now, we have the departing angles, but we need to compare them with each other.	
		"""

		
		
		foreTerm1 = frontInterior1 and frontExist1
		foreTerm2 = frontInterior2 and frontExist2

		frontAngDiff = diffAngle(depAngle1, depAngle2)

		pathBranchIDs = [-1, -1]
		isBranch = [False, False]
		
		" 60 degree threshold "
		ANG_THRESH = 1.047
		#ANG_THRESH = 0.3
		
		if foreTerm1 or foreTerm2:

			if foreTerm1 and foreTerm2:
				" both departing case: check if same or different "
				if fabs(frontAngDiff) < ANG_THRESH:
					" are on top of each other and are branching to the same path "

					if self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1) and self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2):

						
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						pathBranchIDs[1] = newPathID

						isBranch[0] = True
						isBranch[1] = True

						print "foreA"
				else:
					" these are branching to separate new paths "

					if self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1) and self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2):


						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID1 = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID1
						isBranch[0] = True



						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID2 = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID2
						isBranch[1] = True
						print "foreB"
			
			elif foreTerm1 and not foreTerm2:
				
				if self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1):
					
					if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
						" both have same departure "
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						pathBranchIDs[1] = newPathID
						isBranch[0] = True
						isBranch[1] = True

						print "foreC"
						
					else:
						" foreTerm1 has unique departure "	
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True

						print "foreD"

			elif foreTerm2 and not foreTerm1:

				if self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2):

					if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

						" both have same departure "
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						pathBranchIDs[1] = newPathID
						isBranch[0] = True
						isBranch[1] = True

						print "foreE"

					else:
						" foreTerm2 has unique departure "	
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
						newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True

						print "foreF"
			else:
				print "foreG"
				" no new departure, just add nodes to the leaf paths "
				pathID = parentPathID1
				branchNodeID = nodeID1
				globalJunctionPoint = depPoint1
				
				if parentPathID1 != parentPathID2:
					print "departing path IDs of two colocated nodes are not the same"
					raise
		
		
		return isBranch, pathBranchIDs
		
		"""
		pathID = parentPathID1
		branchNodeID = nodeID1
		globalJunctionPoint = depPoint1
		depAng = depAngle1
		junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
	
		orderedPathIDs1.insert(0,self.numPaths)
		departures1.insert(0, [False, False])
		interiors1.insert(0, [False, False])
		depPoints1.insert(0, [None, None])
	
		orderedPathIDs2.insert(0,self.numPaths)
		departures2.insert(0, [False, False])
		interiors2.insert(0, [False, False])
		depPoints2.insert(0, [None, None])
	
		poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
		self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
		self.nodeSets[self.numPaths] = []
		newPaths.append(self.numPaths)
		self.numPaths += 1
		"""



	def checkForBranches(self, nodeID1, nodeID2):

		" CHECK FOR A BRANCHING EVENT FROM LEAF "

		" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
		if len(self.paths[0]) == 0:
			" first nodes in path 0 "					
			self.nodeSets[0].append(nodeID1)
			self.nodeSets[0].append(nodeID2)

			" NOT THE FIRST, NOW CHECK FOR BRANCHING FROM PATHS "			
		else:
			
			" departure events for node 1 "
			departures1 = []
			interiors1 = []
			depPoints1 = []
			distances1 = []
			depAngles1 = []

			" departure events for node 2 "
			departures2 = []
			interiors2 = []
			depPoints2 = []
			distances2 = []
			depAngles2 = []

			
			" raw medial axis of union of path nodes "
			paths = {}
			for k in range(len(self.paths)):
				path = self.paths[k]
				if len(path) > 0:
					paths[k] = path

			print "paths:", len(paths)
			for k in range(len(paths)):
				print k, len(paths[k])


			" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
			self.trimmedPaths = self.trimPaths(paths)
			self.drawTrimmedPaths(self.trimmedPaths)
			print "trimmed paths:", len(self.trimmedPaths)

			
			" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
				
			" the overlapping paths are computed from the initial guess of position "
			orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
			orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
			
			print "orderedPathIDs1:", orderedPathIDs1
			print "orderedPathIDs2:", orderedPathIDs2
			
			
			" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
			for pathID in orderedPathIDs1:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
				departures1.append([isExist1,isExist2])
				interiors1.append([isInterior1, isInterior2])
				depPoints1.append([departurePoint1, departurePoint2])
				distances1.append([discDist1, discDist2])
				depAngles1.append([depAngle1, depAngle2])
			
			for pathID in orderedPathIDs2:
				departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
				departures2.append([isExist1,isExist2])
				interiors2.append([isInterior1, isInterior2])
				depPoints2.append([departurePoint1, departurePoint2])
				distances2.append([discDist1, discDist2])
				depAngles2.append([depAngle1, depAngle2])
				
			print "node departures", nodeID1, ":", departures1
			print "node  interiors", nodeID1, ":", interiors1
			print "node departures", nodeID2, ":", departures2
			print "node  interiors", nodeID2, ":", interiors2


			" new junction finding logic "
			" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
			" if a terminal departure exists that is internal, than we have a new junction "
			DISC_THRESH = 0.2

			" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
			frontExist1 = departures1[0][0]
			backExist1 =  departures1[-1][1]
			frontInterior1 = interiors1[0][0]
			backInterior1 = interiors1[-1][1]

			foreTerm1 = frontInterior1 and frontExist1
			backTerm1 = backInterior1 and backExist1
			
			" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
			discForeTerm1 = distances1[0][0] > DISC_THRESH
			discBackTerm1 = distances1[-1][1] > DISC_THRESH
			
			foreAngle1 = depAngles1[0][0]
			backAngle1 = depAngles1[-1][1]
			
			" NODE2: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
			frontExist2 = departures2[0][0]
			backExist2 =  departures2[-1][1]
			frontInterior2 = interiors2[0][0]
			backInterior2 = interiors2[-1][1]

			foreTerm2 = frontInterior2 and frontExist2
			backTerm2 = backInterior2 and backExist2


			discForeTerm2 = distances2[0][0] > DISC_THRESH
			discBackTerm2 = distances2[-1][1] > DISC_THRESH

			foreAngle2 = depAngles2[0][0]
			backAngle2 = depAngles2[-1][1]

			frontAngDiff = diffAngle(foreAngle1, foreAngle2)
			backAngDiff = diffAngle(backAngle1, backAngle2)


			print "pass1:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
			print "pass1:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2

			" check for the cases that internal departure may have a departure point discrepancy, "
			" then we should perform a pose adjustment and recompute departures "
			print "checking discrepancy in departure:", nodeID1, nodeID2, foreTerm1, discForeTerm1, backTerm1, discBackTerm1, foreTerm2, discForeTerm2, backTerm2, discBackTerm2

			#if discForeTerm1 or discBackTerm1 or discForeTerm2 or discBackTerm2:
			if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1 or foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

				print "adjusting pose guess of node because of discrepancy:", nodeID1, nodeID2

				#if discForeTerm1 or discBackTerm1:
				if foreTerm1 and discForeTerm1 or backTerm1 and discBackTerm1:


					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID1].getGlobalGPACPose()[:2]

					splicedPaths1 = self.splicePathIDs(orderedPathIDs1)

					totalGuesses = []
					for path in splicedPaths1:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID1].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					guessPose = totalGuesses[0][2]
					self.nodeHash[nodeID1].setGPACPose(guessPose)


					#if foreTerm1 and discForeTerm1:
					#	pathID = orderedPathIDs1[0]
					#elif backTerm1 and discBackTerm1:
					#	pathID = orderedPathIDs1[-1]
					
						
					#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID1, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
					#self.nodeHash[nodeID1].setGPACPose(guessPose1)
						
				#elif discForeTerm2 or discBackTerm2:
				elif foreTerm2 and discForeTerm2 or backTerm2 and discBackTerm2:

					" control point nearest the GPAC origin "
					globalPoint1 = self.nodeHash[nodeID2].getGlobalGPACPose()[:2]


					splicedPaths2 = self.splicePathIDs(orderedPathIDs2)

					totalGuesses = []
					for path in splicedPaths2:
	
						" make the path constraints "								
						guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, path, globalPoint1, globalPoint1)
						
						estPose = self.nodeHash[nodeID2].getGlobalGPACPose()
						angDiff = abs(estPose[2]-guessPose1[2])
						totalGuesses.append((angDiff, cost1, guessPose1))
					
					totalGuesses.sort()
					guessPose = totalGuesses[0][2]
					self.nodeHash[nodeID2].setGPACPose(guessPose)

					#if foreTerm2 and discForeTerm2:
					#	pathID = orderedPathIDs2[0]
					#elif backTerm2 and discBackTerm2:
					#	pathID = orderedPathIDs2[-1]
						
					#guessPose1, cost1 = self.makeGlobalMedialOverlapConstraint(nodeID2, self.trimmedPaths[pathID], globalPoint1, globalPoint1)
					#self.nodeHash[nodeID2].setGPACPose(guessPose1)

				departures1 = []
				interiors1 = []
				depPoints1 = []
				distances1 = []
				depAngles1 = []
	
				departures2 = []
				interiors2 = []
				depPoints2 = []
				distances2 = []
				depAngles2 = []
				
				" the overlapping paths are computed from the initial guess of position "
				orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
				orderedPathIDs2 = self.getOrderedOverlappingPaths(nodeID2)
				
				print "orderedPathIDs1:", orderedPathIDs1
				print "orderedPathIDs2:", orderedPathIDs2
				
				" now compute whether there are departure points after we have guessed a better position in synch with the paths "
				for pathID in orderedPathIDs1:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID1)
					departures1.append([isExist1,isExist2])
					interiors1.append([isInterior1, isInterior2])
					depPoints1.append([departurePoint1, departurePoint2])
					distances1.append([discDist1, discDist2])
					depAngles1.append([depAngle1, depAngle2])
				
				for pathID in orderedPathIDs2:
					departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2 = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID2)
					departures2.append([isExist1,isExist2])
					interiors2.append([isInterior1, isInterior2])
					depPoints2.append([departurePoint1, departurePoint2])
					distances2.append([discDist1, discDist2])
					depAngles2.append([depAngle1, depAngle2])
					
				print "node", nodeID1, ":", departures1, interiors1
				print "node", nodeID2, ":", departures2, interiors2
	
				" new junction finding logic "
				" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
				" if a terminal departure exists that is internal, than we have a new junction "
				
				frontExist1 = departures1[0][0]
				backExist1 =  departures1[-1][1]
				frontInterior1 = interiors1[0][0]
				backInterior1 = interiors1[-1][1]

				foreTerm1 = frontInterior1 and frontExist1
				backTerm1 = backInterior1 and backExist1
				
				
				discForeTerm1 = distances1[0][0] > DISC_THRESH
				discBackTerm1 = distances1[-1][1] > DISC_THRESH
				
				foreAngle1 = depAngles1[0][0]
				backAngle1 = depAngles1[-1][1]
				
				#foreTerm2 = departures2[0][0] and interiors2[0][0]
				#backTerm2 = departures2[-1][1] and interiors2[-1][1]
				frontExist2 = departures2[0][0]
				backExist2 =  departures2[-1][1]
				frontInterior2 = interiors2[0][0]
				backInterior2 = interiors2[-1][1]

				foreTerm2 = frontInterior2 and frontExist2
				backTerm2 = backInterior2 and backExist2


				discForeTerm2 = distances2[0][0] > DISC_THRESH
				discBackTerm2 = distances2[-1][1] > DISC_THRESH

				foreAngle2 = depAngles2[0][0]
				backAngle2 = depAngles2[-1][1]
				

				frontAngDiff = diffAngle(foreAngle1, foreAngle2)
				backAngDiff = diffAngle(backAngle1, backAngle2)

				print "pass2:", frontExist1, backExist1, frontInterior1, backInterior1, foreTerm1, backTerm1, discForeTerm1, discBackTerm1, foreAngle1, backAngle1
				print "pass2:", frontExist2, backExist2, frontInterior2, backInterior2, foreTerm2, backTerm2, discForeTerm2, discBackTerm2, foreAngle2, backAngle2
				
			print "frontAngDiff:", frontAngDiff
			print "backAngDiff:", backAngDiff


			newPaths = []

			frontExist1 = departures1[0][0]
			frontInterior1 = interiors1[0][0]
	
			frontExist2 = departures2[0][0]
			frontInterior1 = interiors2[0][0]
			
			depAngle1 = depAngles1[0][0]
			depAngle2 = depAngles2[0][0]
			
			depPoint1 = depPoints1[0][0]
			depPoint2 = depPoints2[0][0]
			
			parentPathID1 = orderedPathIDs1[0]
			parentPathID2 = orderedPathIDs2[0]


			isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)

			if isBranch[0]:
				orderedPathIDs1.insert(0,pathBranchIDs[0])
				departures1.insert(0, [False, False])
				interiors1.insert(0, [False, False])
				depPoints1.insert(0, [None, None])

			if isBranch[1]:
				orderedPathIDs2.insert(0,pathBranchIDs[1])
				departures2.insert(0, [False, False])
				interiors2.insert(0, [False, False])
				depPoints2.insert(0, [None, None])
				
			for pathID in pathBranchIDs:
				if pathID != -1 and newPaths.count(pathID) == 0:
					newPaths.append(pathID)
			
			#for newPath in addThesePaths:
			#	self.pathParents.append(newPath)
			#	self.nodeSets[newPath[0]] = []
			#	newPaths.append(newPath[0])
			#	self.numPaths += 1


			backExist1 = departures1[-1][1]
			backInterior1 = interiors1[-1][1]
	
			backExist2 = departures2[-1][1]
			backInterior1 = interiors2[-1][1]
			
			depAngle1 = depAngles1[-1][1]
			depAngle2 = depAngles2[-1][1]
			
			depPoint1 = depPoints1[-1][1]
			depPoint2 = depPoints2[-1][1]

			parentPathID1 = orderedPathIDs1[-1]
			parentPathID2 = orderedPathIDs2[-1]


			isBranch, pathBranchIDs = self.determineBranch(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2)

			if isBranch[0]:
				orderedPathIDs1.append(0,pathBranchIDs[0])
				departures1.append([False, False])
				interiors1.append([False, False])
				depPoints1.append([None, None])

			if isBranch[1]:
				orderedPathIDs2.append(pathBranchIDs[1])
				departures2.append([False, False])
				interiors2.append([False, False])
				depPoints2.append([None, None])

			for pathID in pathBranchIDs:
				if pathID != -1 and newPaths.count(pathID) == 0:
					newPaths.append(pathID)
			
			#for newPath in addThesePaths:
			#	self.pathParents.append(newPath)
			#	self.nodeSets[newPath[0]] = []
			#	newPaths.append(newPath[0])
			#	self.numPaths += 1

	



			"""
			newPaths = []
			
			" 60 degree threshold "
			ANG_THRESH = 1.047
			#ANG_THRESH = 0.3
			
			if foreTerm1 or foreTerm2:

				if foreTerm1 and foreTerm2:
					" both departing case: check if same or different "
					if fabs(frontAngDiff) < ANG_THRESH:
						" are on top of each other and are branching to the same path "

						if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]) and self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):

							
							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1

							print "foreA"
					else:
						" these are branching to separate new paths "

						if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]) and self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):


							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1



							pathID = orderedPathIDs2[0]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[0][0]
							depAng = depAngles2[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])


							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreB"
				
				elif foreTerm1 and not foreTerm2:
					
					if self.checkUniqueBranch(orderedPathIDs1[0], nodeID1, depAngles1[0][0], depPoints1[0][0]):
						
						if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
							" both have same departure "
							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
							#depAng1 = depAngles1[0][0]
							#depAng2 = depAngles2[0][0]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreC"
							
						else:
							" foreTerm1 has unique departure "	
							pathID = orderedPathIDs1[0]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[0][0]
							depAng = depAngles1[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreD"

				elif foreTerm2 and not foreTerm1:

					if self.checkUniqueBranch(orderedPathIDs2[0], nodeID2, depAngles2[0][0], depPoints2[0][0]):

						if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

							" both have same departure "
							pathID = orderedPathIDs2[0]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[0][0]
							depAng = depAngles2[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.insert(0,self.numPaths)
							departures1.insert(0, [False, False])
							interiors1.insert(0, [False, False])
							depPoints1.insert(0, [None, None])

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "foreE"

						else:
							" foreTerm2 has unique departure "	
							pathID = orderedPathIDs2[0]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[0][0]
							depAng = depAngles2[0][0]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.insert(0,self.numPaths)
							departures2.insert(0, [False, False])
							interiors2.insert(0, [False, False])
							depPoints2.insert(0, [None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1

							print "foreF"
				else:
					print "foreG"
					" no new departure, just add nodes to the leaf paths "
					pathID = orderedPathIDs1[0]
					branchNodeID = nodeID1
					globalJunctionPoint = depPoints1[0][0]
					
					if orderedPathIDs1[0] != orderedPathIDs2[0]:
						print "departing path IDs of two colocated nodes are not the same"
						raise

			if backTerm1 or backTerm2:

				if backTerm1 and backTerm2:
					" both departing case: check if same or different "
					if fabs(backAngDiff) < ANG_THRESH:
						" are on top of each other and are branching to the same path "

						if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]) and self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):

							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backA"

					else:
						" these are branching to separate new paths "
						if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]) and self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):

							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1

							pathID = orderedPathIDs2[-1]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[-1][1]
							depAng = depAngles2[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backB"					
				elif backTerm1 and not backTerm2:
					
					if self.checkUniqueBranch(orderedPathIDs1[-1], nodeID1, depAngles1[-1][1], depPoints1[-1][1]):
						
						if backExist2 and fabs(backAngDiff) < ANG_THRESH:
							" both have same departure "
							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backC"
						else:
							" backTerm1 has unique departure "	
							pathID = orderedPathIDs1[-1]
							branchNodeID = nodeID1
							globalJunctionPoint = depPoints1[-1][1]
							depAng = depAngles1[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backD"

				elif backTerm2 and not backTerm1:

					if self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1]):

						if backExist1 and fabs(backAngDiff) < ANG_THRESH:

							" both have same departure "
							pathID = orderedPathIDs2[-1]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[-1][1]
							depAng = depAngles2[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs1.append(self.numPaths)
							departures1.append([False, False])
							interiors1.append([False, False])
							depPoints1.append([None, None])

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backE"
						else:
							" backTerm2 has unique departure "	
							self.checkUniqueBranch(orderedPathIDs2[-1], nodeID2, depAngles2[-1][1], depPoints2[-1][1])

							pathID = orderedPathIDs2[-1]
							branchNodeID = nodeID2
							globalJunctionPoint = depPoints2[-1][1]
							depAng = depAngles2[-1][1]
							junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

							orderedPathIDs2.append(self.numPaths)
							departures2.append([False, False])
							interiors2.append([False, False])
							depPoints2.append([None, None])

							poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
							self.pathParents.append([pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug), True])
							self.nodeSets[self.numPaths] = []
							newPaths.append(self.numPaths)
							self.numPaths += 1
							print "backF"
				else:
					print "backG"
					" no new departure, just add nodes to the leaf paths "
					pathID = orderedPathIDs1[-1]
					branchNodeID = nodeID1
					globalJunctionPoint = depPoints1[-1][1]
					if orderedPathIDs1[-1] != orderedPathIDs2[-1]:
						print "departing path IDs of two colocated nodes are not the same"
						raise
			"""





			" determine which paths are leaves "
			print "self.pathParents:", self.pathParents
			isAParent = [False for k in range(self.numPaths)]
			for k in orderedPathIDs1:
				print "index:", k
				currPath = self.pathParents[k]
				currParent = currPath[0]
				if currParent != None:
					isAParent[currParent] = True
			
			"add nodes to paths that are the leaves "
			for pathID in orderedPathIDs1:
				if not isAParent[pathID]:				
					self.nodeSets[pathID].append(nodeID1)

			isAParent = [False for k in range(self.numPaths)]
			for k in orderedPathIDs2:
				print "index:", k
				currPath = self.pathParents[k]
				currParent = currPath[0]
				if currParent != None:
					isAParent[currParent] = True

			for pathID in orderedPathIDs2:
				if not isAParent[pathID]:				
					self.nodeSets[pathID].append(nodeID2)

			paths = {}
			for k in range(len(self.paths)):
				path = self.paths[k]
				if len(path) > 0:
					paths[k] = path	
			self.trimmedPaths = self.trimPaths(paths)

			" perform path constraints only if we have not created a new path "
			isNew1 = False
			isNew2 = False
			for pathID in newPaths:
				if orderedPathIDs1.count(pathID) > 0:
					isNew1 = True
				if orderedPathIDs2.count(pathID) > 0:
					isNew2 = True
			


			self.drawTrimmedPaths(self.trimmedPaths)
			print "trimmed paths:", len(self.trimmedPaths)



			"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
			"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
			"3)  select choice with the lowest cost "
			"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
			
			" nodeID1:  is it in a junction or a single path? "

			if not isNew1:
				self.constrainToPaths(nodeID1, orderedPathIDs1, departures1, interiors1, depPoints1)
				
			if not isNew2:
				self.constrainToPaths(nodeID2, orderedPathIDs2, departures2, interiors2, depPoints2)

			self.mergePriorityConstraints()


		for k in range(len(self.nodeSets)):
			print "computing path for node set", k, ":", self.nodeSets[k]
			self.paths[k], self.hulls[k] = self.getTopology(self.nodeSets[k])

		self.drawPathAndHull()
		
		return
	
		


	def checkUniqueBranch(self, parentPathID, nodeID1, depAngle, depPoint):

		" check cartesian distance to similar junction points from parent path "

		" cartesian distance "
		DISC_THRESH = 1.0

		" 60 degree threshold "
		ANG_THRESH = 1.047
		
		for path in self.pathParents:
			if path[0] == parentPathID:
				
				junctionNodeID = path[1]
				localJunctionPoint = path[2]
				" check dist "

				poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
				junctionPoint = poseOrigin.convertLocalOffsetToGlobal(localJunctionPoint)
				
				dist = sqrt((depPoint[0]-junctionPoint[0])**2 + (depPoint[1]-junctionPoint[1])**2 )

				" check difference of tangential angle "
				angDiff = diffAngle(depAngle, junctionPoint[2])

				print "checkUniqueBranch(", parentPathID, ",", nodeID1, ",", depAngle, ",", depPoint, ") = ", junctionNodeID, dist, angDiff
		
				if dist < DISC_THRESH:
					if fabs(angDiff) < ANG_THRESH:
						
						print "DUPLICATE junction of", junctionPoint, "rejecting with differences of", dist, angDiff
						return False
		
		return True

	def correctNode(self, nodeID):

		global topCount

		#self.currNode = newNode
		#nodeID = self.numNodes
		#self.nodeHash[nodeID] = self.currNode
		#self.numNodes += 1
		self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		
		
		if nodeID > 0:
			if nodeID % 2 == 0:

				direction = self.currNode.travelDir

				transform, covE, hist = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )

				self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)
	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" try to perform corner constraints "
				" overwrite overlap/motion constraints"
				" corner constraints with past nodes "
				self.addCornerConstraints(nodeID)
			
				
				" merge the constraints "
				self.mergePriorityConstraints()


				" check that each bin is well-connected, otherwise, try to add constraints "
				for bin in self.cornerBins:
					pass
			else:
	
				direction = self.currNode.travelDir

				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY)

				if nodeID > 2:
					#transform, covE = self.makeOverlapConstraint2(nodeID-2, nodeID)
					transform, covE, hist = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction)
					self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)

	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" try to perform corner constraints "
				" overwrite overlap/motion constraints"
				" corner constraints with past nodes "
				self.addCornerConstraints(nodeID)

				self.mergePriorityConstraints()
				
				" check that each bin is well-connected, otherwise, try to add constraints "
				for bin in self.cornerBins:
					pass
		
		if nodeID > 0 and nodeID % 10 == 0:
			self.performOverlapConstraints()

				

	def resetGraphToGround(self):
		
		self.deleteAllEdges()

		for i in range(self.numNodes):
			" reset the estimated pose to ground pose "
			pastNode = self.nodeHash[i]
			gndPose = pastNode.getGndPose()
			#gndProf = Pose(gndPose)
			#gndRootPose = pastNode.gndRootPose
			#pastNode.setEstPose(gndProf.convertLocalOffsetToGlobal(gndRootPose))
			pastNode.setEstPose(gndPose)

		" now create ground constraints "
		" ground constraints are generated from ground truth data "
		self.gnd_constraints = []

		for i in range(1,self.numNodes):

			if self.numNodes > 0:
				
				node1 = self.nodeHash[i-1]
				node2 = self.nodeHash[i]
				
				" create the ground constraints "
				gndGPAC1Pose = node1.getGndGlobalGPACPose()
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = node2.getGndGlobalGPACPose()
				offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
				self.gnd_constraints.append([i-1, i, matrix([[offset[0]], [offset[1]], [offset[2]]]), self.E_gnd])

		for const in self.gnd_constraints:
			#self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)

			node1 = const[0]
			node2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([node1,node2,transform,covE], GND_PRIORITY)

	def loadNodeGroundPast(self, newNode):

		#if self.numNodes > 0:
		#	self.resetGraphToGround()						

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))
		self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))

		if nodeID > 0:

			if nodeID % 2 == 0:
				" 2nd node from a pose "
				" 1) make overlap+motion constraint with i-1"
				" 2) attempt overlap with past nodes ? "
				" 3) attempt corner constraints with past nodes, replacing overlaps and sensors "
				" 4) discard sensor constraints and reapply sensor constraints to non-corner edges"

				#transform, covE = self.makeMotionConstraint(nodeID-1, nodeID)
				#self.edgePriorityHash[(nodeID-1, nodeID)].append([transform, covE, MOTION_PRIORITY])
				
				transform, covE = self.makeOverlapConstraint2(nodeID-2, nodeID)
				#transform, covE = self.makeOverlapConstraint(nodeID-2, nodeID)
				self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)
	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" try to perform corner constraints "
				" overwrite overlap/motion constraints"
				" corner constraints with past nodes "
				#self.addCornerConstraints(nodeID)
				
				
				" merge the constraints "
				#self.mergePriorityConstraints()


				" check that each bin is well-connected, otherwise, try to add constraints "
				for bin in self.cornerBins:
					pass
									
			else:
						
				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY)
	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" try to perform corner constraints "
				" overwrite overlap/motion constraints"
				" corner constraints with past nodes "
				#self.addCornerConstraints(nodeID)

				#self.mergePriorityConstraints()
				


		if nodeID > 1:

			paths = []
			for i in range(self.numNodes):
				paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
			
			newConstraints = []
			newConstraints, newConstraints2 = self.addShapeConstraints(paths, targetNode = nodeID)
			
			" remove previous shape constraints "
			self.deleteAllPriority(SHAPE_PRIORITY)
				
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

			totalError1 = self.computeCornerClusterError()
			print "totalError1:", totalError1


			self.deleteAllPriority(SHAPE_PRIORITY)
			for const in newConstraints2:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

			totalError2 = self.computeCornerClusterError()
			print "totalError2:", totalError2

			" revert to first hypothesis set if other is better "
			if totalError2 > totalError1:
				self.deleteAllPriority(SHAPE_PRIORITY)
					
				for const in newConstraints:
					n1 = const[0]
					n2 = const[1]
					transform = const[2]
					covE = const[3]
		
					self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
		
				" merge the constraints "
				self.mergePriorityConstraints()

		if False:

			paths = []
			for i in range(self.numNodes):
				paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
			
			newConstraints = []
			self.newConstraints2 = []
			newConstraints, self.newConstraints2 = self.addShapeConstraints(paths)
			#newConstraints = self.addShapeConstraints(paths, targetNode = nodeID)
			
			" remove previous shape constraints "
			#deleteAll(self.edgePriorityHash, SHAPE_PRIORITY)
			self.deleteAllPriority(SHAPE_PRIORITY)
				
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

			totalError1 = self.computeCornerClusterError()

			self.deleteAllPriority(SHAPE_PRIORITY)
			for const in self.newConstraints2:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

			totalError2 = self.computeCornerClusterError()

			" revert to first hypothesis set if other is better "
			if totalError2 > totalError1:
				self.deleteAllPriority(SHAPE_PRIORITY)
					
				for const in newConstraints:
					n1 = const[0]
					n2 = const[1]
					transform = const[2]
					covE = const[3]
		
					self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
		
				" merge the constraints "
				self.mergePriorityConstraints()
		


	def performOverlapConstraints(self):
		
		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
		
		self.deleteAllPriority(SHAPE_PRIORITY)

		newConstraints = self.addOverlapConstraints(paths)

		print "adding", len(newConstraints), "overlap constraints"

		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()

	def performShapeConstraints(self):
		
		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		totalInitError = self.computeCornerClusterError()
		
		newConstraints = []
		self.newConstraints2 = []
		newConstraints, self.newConstraints2 = self.addShapeConstraints(paths)
		#newConstraints = self.addShapeConstraints(paths, targetNode = nodeID)
		
		" remove previous shape constraints "
		#deleteAll(self.edgePriorityHash, SHAPE_PRIORITY)
		self.deleteAllPriority(SHAPE_PRIORITY)
			
		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()

		totalError1 = self.computeCornerClusterError()

		self.deleteAllPriority(SHAPE_PRIORITY)
		for const in self.newConstraints2:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()

		totalError2 = self.computeCornerClusterError()

		" revert to first hypothesis set if other is better "
		if totalError2 > totalError1:
			self.deleteAllPriority(SHAPE_PRIORITY)
				
			for const in newConstraints:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

		
	def mergePriorityConstraints(self):
		
		totalConstraints = []
		
		#print "merging:"
		
		" merge the constraint "
		for k, v in self.edgePriorityHash.items():
			
			if len(v) > 0:
		
				id1 = k[0]
				id2 = k[1]
				
				" take the highest priority constraints only "
				maxPriority = -1
				for const in v:
					thisPriority = const[2]
					if thisPriority > maxPriority:
						maxPriority = thisPriority
		
				if maxPriority == -1:
					raise
				
				for const in v:
					if const[2] == maxPriority:
						transform = const[0]
						covE = const[1]
						#print id1, id2, maxPriority
						totalConstraints.append([id1,id2,transform,covE])	
		
		#print "merging totalConstraints:", totalConstraints
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")
		
		self.edgeHash = {}
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1, node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]
		
		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])
					
	def loadNodes(self, nodeHash):

		self.nodeHash = nodeHash
		
		self.numNodes = len(self.nodeHash)
		self.currNode = self.nodeHash[self.numNodes-1]

		" ground constraints are generated from ground truth data "
		self.gnd_constraints = []

		for i in range(1,self.numNodes):

			if self.numNodes > 0:
				
				node1 = self.nodeHash[i-1]
				node2 = self.nodeHash[i]
				
				" create the ground constraints "
				gndGPAC1Pose = node1.getGndGlobalGPACPose()
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = node2.getGndGlobalGPACPose()
				offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
				self.gnd_constraints.append([i-1, i, matrix([[offset[0]], [offset[1]], [offset[2]]]), matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])])

	def loadConstraints(self, dirName, index):
		
		self.inplace_constraints = []
		self.overlap_constraints = []
		self.motion_constraints = []
		self.sensor_constraints = []	
		self.merged_constraints = []
		self.backtrack_constraints = []

		" load the constraints from file "
		f = open(dirName + "/motion_constraints_%04u.txt" % index, 'r')
		self.motion_constraints = eval(f.read().rstrip())
		f.close()
		f = open(dirName + "/overlap_constraints_%04u.txt" % index, 'r')
		self.overlap_constraints = eval(f.read().rstrip())
		f.close()
		f = open(dirName + "/sensor_constraints_%04u.txt" % index, 'r')
		self.sensor_constraints = eval(f.read().rstrip())
		f.close()
		f = open(dirName + "/inplace_constraints_%04u.txt" % index, 'r')
		self.inplace_constraints = eval(f.read().rstrip())
		f.close()
		
		f = open(dirName + "/backtrack_constraints_%04u.txt" % index, 'r')
		self.backtrack_constraints = eval(f.read().rstrip())
		f.close()

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.motion_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.overlap_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])


		#print "totalConstraints:"
		#print totalConstraints
		
		" merge the constraints with TORO"
		#self.merged_constraints = self.toroMerge(totalConstraints)
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")
		
		#print "merged_constraints:", self.merged_constraints
		
		self.edgeHash = {}
		for i in range(len(self.merged_constraints)):
			#print i, self.merged_constraints[i]
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1, node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]
		
		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])		

	def makeHypotheses(self):
		
		global renderCount

		if self.numNodes < 2:
			return

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.motion_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.overlap_constraints[i][3]			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		#for i in range(len(self.activeHypotheses)):
		#	nodeID1 = self.activeHypotheses[i][0]
		#	nodeID2 = self.activeHypotheses[i][1]
		#	transform = self.activeHypotheses[i][2]
		#	covE = self.activeHypotheses[i][3]
		#	totalConstraints.append([nodeID1,nodeID2,transform,covE])


		" merge the constraints with TORO"
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")

		newEdgeHash = {}		
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			newEdgeHash[(node1,node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]

		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, newEdgeHash))

		#cornerPairs = self.findCornerCandidates(paths)
		

		" find candidates to perform sensor constraints "
		#sensor_pairs = self.findCornerCandidates(paths)
		#sensor_pairs = self.findSweepCandidates(paths)
		sensor_pairs = self.findAlphaHullCandidates(paths, targetNode = self.numNodes-1)
		
		" compute the alpha hulls of the selected candidates "
		hull_computed = [False for i in range(self.numNodes)]
		a_hulls = [0 for i in range(self.numNodes)]
		
		
		
		for i in range(self.numNodes):
			a_hulls[i] = computeHull(self.nodeHash[i], sweep = False)
			#a_hulls[i] = computeHull(self.nodeHash[i], sweep = True)
			hull_computed[i] = True
		
		" make hypothesized constraints out of the pairs "
		hypotheses = []		
		for i in range(len(sensor_pairs)):
			p = sensor_pairs[i]
			transform = p[3]
			n1 = p[1]
			n2 = p[2]
			offset, covar, cost = self.makeSensorConstraint(transform, n1, n2, a_hulls[n1], a_hulls[n2])
			transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
		
			if cost < 1.5:
				hypotheses.append([p[1],p[2],transform,covar])
				print "hypothesis", i, ":",  "sensor constraint between", n1, "and", n2

		" compute the consistency between each pair of hypotheses "

		self.sensorHypotheses += hypotheses
		
		totalHypotheses = self.sensor_constraints + self.sensorHypotheses
		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) < 2:
			print "too little hypotheses, rejecting"
			return [] 
		
		" set all proscribed hypothesis sets to maximum error "	
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]

		maxError = maxError*2

		for result in results:
			if result[0] == None:
				result[0] = maxError
	
		" create the consistency matrix "
		A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
		A.resize(len(totalHypotheses), len(totalHypotheses))
			
		" populate consistency matrix "					
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]

		" disable the truncation "
		set_printoptions(threshold=nan)		
		fn = open("A_sensor.txt",'w')
		fn.write(repr(A))
		fn.close()
		
		
		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			#w = scsgp.getIndicatorVector(e[0])
			w = scsgp.getIndicatorVector(e[1])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break
		
		#fp = open("w_sensor.txt", 'w')
		#fp.write(repr(w))
		#fp.close()

		
		matrix([[ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 1.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.],
		        [ 0.]])


		
		#print "INDICATOR:", repr(w)
		selected = []	
		for i in range(len(totalHypotheses)):
			if w[i,0] >= 1.0:
				selected.append(totalHypotheses[i])
		
		
		
		#print len(selected), "added hypotheses"
		#print len(self.merged_constraints), "merged constraints"
		
		self.activeHypotheses = selected

		" merge the constraints with TORO"
		#v_list, self.merged_constraints = self.doToro(self.merged_constraints + self.activeHypotheses, fileName = "probe_init")
		
		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.motion_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.overlap_constraints[i][3]			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.activeHypotheses)):
			nodeID1 = self.activeHypotheses[i][0]
			nodeID2 = self.activeHypotheses[i][1]
			transform = self.activeHypotheses[i][2]
			covE = self.activeHypotheses[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])


		" merge the constraints with TORO"
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")



		self.edgeHash = {}		
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1,node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]

		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])		


		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		cornerPairs = self.findCornerCandidates(paths)

		b_hulls = [0 for i in range(self.numNodes)]		
		for i in range(self.numNodes):
			b_hulls[i] = computeHull(self.nodeHash[i], sweep = True)


		self.goodHypotheses = []
		self.badHypotheses = []
		hypotheses = []
		for k in range(len(cornerPairs)):
			
			p = cornerPairs[k]
			dist = p[0]
			n1 = p[1]
			n2 = p[2]
			transform = p[3]
			angDist = p[4]
			point1 = p[5]
			point2 = p[6]
			ang1 = p[7]
			ang2 = p[8]
			point1_i = p[9]
			point2_i = p[10]	

			hull1 = a_hulls[n1]
			hull2 = a_hulls[n2]

			offset, covar, cost, angDiff, oriDiff = self.makeCornerConstraint(n1, n2, point1, point2, ang1, ang2, hull1, hull2, b_hulls[n1], b_hulls[n2])
 
			transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
	
			#print "hypothesis", k, ":",  "sensor constraint between", n1, "and", n2, " cost =", cost

			self.allCornerHypotheses.append([n1, n2, transform, covar, cost, point1, point2, angDiff, point1_i, point2_i])

			" ratio of match error threshold"		
			if cost <= 0.6:
				hypotheses.append([n1, n2, transform, covar, cost, point1, point2])
				self.goodHypotheses.append([n1, n2, transform, covar, cost, point1, point2, angDiff, point1_i, point2_i])

			else:
				self.badHypotheses.append([n1, n2, transform, covar, cost, point1, point2, angDiff, point1_i, point2_i])

		set_printoptions(threshold=nan)		
		fn = open("cornerHypotheses.txt",'w')
		fn.write(repr(self.allCornerHypotheses))
		fn.close()
				

		for hyp in hypotheses:
			n1 = hyp[0]
			n2 = hyp[1]
			newCost = hyp[4]
			
			try:
				if n1 > n2:
					otherHyp = self.cornerHash[(n1,n2)]
				else:
					otherHyp = self.cornerHash[(n2,n1)]
			
				oldCost = otherHyp[4]
				if newCost < oldCost:
					if n1 > n2:
						self.cornerHash[(n1,n2)] = hyp
					else:
						self.cornerHash[(n2,n1)] = hyp
			except:				
				if n1 > n2:
					self.cornerHash[(n1,n2)] = hyp
				else:
					self.cornerHash[(n2,n1)] = hyp			

		#self.cornerHypotheses += hypotheses

		#totalHypotheses = self.motion_constraints + self.overlap_constraints + self.inplace_constraints + self.backtrack_constraints + self.sensor_constraints + self.corner_constraints + self.cornerHypotheses
		#totalHypotheses =  self.sensor_constraints + self.corner_constraints + self.cornerHypotheses
		#totalHypotheses =  self.corner_constraints + self.cornerHypotheses
		
		totalHypotheses = []
		for k, v in self.cornerHash.items():
			print "corner constraint:", v[0], "->", v[1]
			totalHypotheses.append([v[0], v[1], v[2], v[3]])
			self.renderCornerConstraint(v[0], v[1], v[2], v[5], v[6], renderCount)
			renderCount += 1

		for i in range(len(self.goodHypotheses)):
			v = self.goodHypotheses[i]
			id1 = v[0]
			id2 = v[1]
			transform = v[2]
			covar = v[3]
			cost = v[4]
			point1 = v[5]
			point2 = v[6]
			angDiff = v[7]
			self.renderGoodAndBadConstraint(id1, id2, transform, covar, cost, point1, point2, angDiff, "goodConstraint", i)
			
		for i in range(len(self.badHypotheses)):
			v = self.badHypotheses[i]
			id1 = v[0]
			id2 = v[1]
			transform = v[2]
			covar = v[3]
			cost = v[4]
			point1 = v[5]
			point2 = v[6]
			angDiff = v[7]
			self.renderGoodAndBadConstraint(id1, id2, transform, covar, cost, point1, point2, angDiff, "badConstraint", i)

		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		#print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) < 2:
			print "too little hypotheses, rejecting"
			return []
		
		" set all proscribed hypothesis sets to maximum error "	
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]

		maxError = maxError*2

		for result in results:
			if result[0] == None:
				result[0] = maxError
	
		" create the consistency matrix "
		A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
		A.resize(len(totalHypotheses), len(totalHypotheses))
			
		" populate consistency matrix "					
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]
		
		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			#w = scsgp.getIndicatorVector(e[0])
			w = scsgp.getIndicatorVector(e[1])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break
			
		#print "INDICATOR:", w
		selected = []	
		for i in range(len(totalHypotheses)):
			selected.append(totalHypotheses[i])
			#if w[i,0] >= 1.0:
			#	selected.append(totalHypotheses[i])
		
		print len(selected), "added hypotheses"
		#print len(self.merged_constraints), "merged constraints"
		
		self.activeHypotheses = selected

		" merge the constraints with TORO"
		#v_list, self.merged_constraints = self.doToro(self.merged_constraints + self.activeHypotheses, fileName = "probe_init")
		
		totalConstraints = []
		"""
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.motion_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		"""
				
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.overlap_constraints[i][3]			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		"""
		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		"""

		for i in range(len(self.activeHypotheses)):
			nodeID1 = self.activeHypotheses[i][0]
			nodeID2 = self.activeHypotheses[i][1]
			transform = self.activeHypotheses[i][2]
			covE = self.activeHypotheses[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])


		" merge the constraints with TORO"
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")



		self.edgeHash = {}		
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1,node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]

		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])		



	def processCornerHypotheses(self):
		
		hypotheses = []

		for i in range(len(self.allCornerHypotheses)):
			v = self.allCornerHypotheses[i]
			n1 = v[0]
			n2 = v[1]
			transform = v[2]
			covar = v[3]
			cost = v[4]
			point1 = v[5]
			point2 = v[6]
			angDiff = v[7]
			corner1_i = v[8]
			corner2_i = v[9]
	
			isOkay, angDiff2, oriDiff = self.checkOrientationError([transform[0,0], transform[1,0], transform[2,0]], n1, n2)
			
			" ratio of match error threshold"		
			if cost <= 0.6 and isOkay:
				hypotheses.append([n1, n2, transform, covar, cost, point1, point2])
				self.goodHypotheses.append([n1, n2, transform, covar, cost, point1, point2, angDiff2, oriDiff])
			else:
				self.badHypotheses.append([n1, n2, transform, covar, cost, point1, point2, angDiff2, oriDiff])

		
		
		for hyp in hypotheses:
			n1 = hyp[0]
			n2 = hyp[1]
			newCost = hyp[4]
			
			try:
				if n1 > n2:
					otherHyp = self.cornerHash[(n1,n2)]
				else:
					otherHyp = self.cornerHash[(n2,n1)]
			
				oldCost = otherHyp[4]
				if newCost < oldCost:
					if n1 > n2:
						self.cornerHash[(n1,n2)] = hyp
					else:
						self.cornerHash[(n2,n1)] = hyp
			except:				
				if n1 > n2:
					self.cornerHash[(n1,n2)] = hyp
				else:
					self.cornerHash[(n2,n1)] = hyp			
		
		totalHypotheses = []
		for k, v in self.cornerHash.items():
			print "corner constraint:", v[0], "->", v[1]
			totalHypotheses.append([v[0], v[1], v[2], v[3]])
			#self.renderCornerConstraint(v[0], v[1], v[2], v[5], v[6], renderCount)
			#renderCount += 1

		for i in range(len(self.goodHypotheses)):
			v = self.goodHypotheses[i]
			id1 = v[0]
			id2 = v[1]
			transform = v[2]
			covar = v[3]
			cost = v[4]
			point1 = v[5]
			point2 = v[6]
			angDiff = v[7]
			oriDiff = v[8]
			self.renderGoodAndBadConstraint(id1, id2, transform, covar, cost, point1, point2, angDiff, oriDiff, "goodConstraint", i)
		
		"""
		for i in range(len(self.badHypotheses)):
			v = self.badHypotheses[i]
			id1 = v[0]
			id2 = v[1]
			transform = v[2]
			covar = v[3]
			cost = v[4]
			point1 = v[5]
			point2 = v[6]
			angDiff = v[7]
			self.renderGoodAndBadConstraint(id1, id2, transform, covar, cost, point1, point2, angDiff, "badConstraint", i)
		"""
		
		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		#print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) < 2:
			print "too little hypotheses, rejecting"
			return []
		
		" set all proscribed hypothesis sets to maximum error "	
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]

		maxError = maxError*2

		for result in results:
			if result[0] == None:
				result[0] = maxError
	
		" create the consistency matrix "
		A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
		A.resize(len(totalHypotheses), len(totalHypotheses))
			
		" populate consistency matrix "					
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]
		
		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			#w = scsgp.getIndicatorVector(e[0])
			w = scsgp.getIndicatorVector(e[1])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break
			
		#print "INDICATOR:", w
		selected = []	
		for i in range(len(totalHypotheses)):
			selected.append(totalHypotheses[i])
			#if w[i,0] >= 1.0:
			#	selected.append(totalHypotheses[i])
		
		print len(selected), "added hypotheses"
		#print len(self.merged_constraints), "merged constraints"
		
		self.activeHypotheses = selected

		" merge the constraints with TORO"
		#v_list, self.merged_constraints = self.doToro(self.merged_constraints + self.activeHypotheses, fileName = "probe_init")
		
		totalConstraints = []
		"""
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.motion_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		"""
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.overlap_constraints[i][3]			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		"""
		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		"""

		for i in range(len(self.activeHypotheses)):
			nodeID1 = self.activeHypotheses[i][0]
			nodeID2 = self.activeHypotheses[i][1]
			transform = self.activeHypotheses[i][2]
			covE = self.activeHypotheses[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])


		" merge the constraints with TORO"
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")



		self.edgeHash = {}		
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1,node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]

		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])		


					
	def computeCovar(self):

		" compute the average position, and the covariance of position "

		N = len(self.gnd_constraints)

		for i in range(N):
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]
			motionConst = self.motion_constraints[i][2]
			motionConst = [motionConst[0,0],motionConst[1,0],motionConst[2,0]]
			overlapConst = self.overlap_constraints[i][2]
			overlapConst = [overlapConst[0,0],overlapConst[1,0],overlapConst[2,0]]

		avgX = 0.0
		avgY = 0.0
		avgPx = 0.0
		avgPy = 0.0
		for i in range(N):
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]

			avgX += gndConst[0]
			avgY += gndConst[1]
			
			P = gndConst[2]
			avgPx += cos(P)
			avgPy += sin(P)
			
		avgX /= N
		avgY /= N
		avgPx /= N
		avgPy /= N
		avgP = atan2(avgPy, avgPx)
		
		#print "average x,y:", avgX, avgY
		#print "average dist:", avgDist
		#print "average angle:", avgP

		stdX = 0.0
		stdY = 0.0
		stdP = 0.0

		Exx = 0.0
		Exy = 0.0
		Exp = 0.0
		Eyp = 0.0
		Eyy = 0.0
		Epp = 0.0
		
		for i in range(N):
			
			motionConst = self.motion_constraints[i][2]
			motionConst = [motionConst[0,0],motionConst[1,0],motionConst[2,0]]
			
			x = motionConst[0]
			y = motionConst[1]
			p = motionConst[2]
			
			Exx += (x - avgX) * (x - avgX)
			Exy += (x - avgX) * (y - avgY)
			Exp += (x - avgX) * (p - avgP)
			Eyp += (y - avgY) * (p - avgP)
			Eyy += (y - avgY) * (y - avgY)
			Epp += (p - avgP) * (p - avgP)
			
		#print "pre:", Exx, Eyy, Epp, Exy, Exp, Eyp
		Exx /= N
		Exy /= N
		Exp /= N
		Eyp /= N
		Eyy /= N
		Epp /= N
		
		E_motion = numpy.matrix([[Exx, Exy, Exp],
					[Exy, Eyy, Eyp],
					[Exp, Eyp, Epp]])
		#print "covar:", Exx, Eyy, Epp, Exy, Exp, Eyp	

		stdX = 0.0
		stdY = 0.0
		stdP = 0.0

		Exx = 0.0
		Exy = 0.0
		Exp = 0.0
		Eyp = 0.0
		Eyy = 0.0
		Epp = 0.0
		
		for i in range(N):
			
			overlapConst = self.overlap_constraints[i][2]
			overlapConst = [overlapConst[0,0],overlapConst[1,0],overlapConst[2,0]]
						
			x = overlapConst[0]
			y = overlapConst[1]
			p = overlapConst[2]
			
			Exx += (x - avgX) * (x - avgX)
			Exy += (x - avgX) * (y - avgY)
			Exp += (x - avgX) * (p - avgP)
			Eyp += (y - avgY) * (p - avgP)
			Eyy += (y - avgY) * (y - avgY)
			Epp += (p - avgP) * (p - avgP)
			
		#print "pre:", Exx, Eyy, Epp, Exy, Exp, Eyp
		Exx /= N
		Exy /= N
		Exp /= N
		Eyp /= N
		Eyy /= N
		Epp /= N
		
		E_overlap = numpy.matrix([[Exx, Exy, Exp],
					[Exy, Eyy, Eyp],
					[Exp, Eyp, Epp]])

		return E_motion, E_overlap


	def computeCovar2(self):
		
		" compute the average error in position, and the covariance of error "
		N = len(self.gnd_constraints)

		for i in range(N):
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]
			motionConst = self.motion_constraints[i][2]
			motionConst = [motionConst[0,0],motionConst[1,0],motionConst[2,0]]
			overlapConst = self.overlap_constraints[i][2]
			overlapConst = [overlapConst[0,0],overlapConst[1,0],overlapConst[2,0]]

		avgX = 0.0
		avgY = 0.0
		avgPx = 0.0
		avgPy = 0.0
		for i in range(N):
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]
			motionConst = self.motion_constraints[i][2]
			motionConst = [motionConst[0,0],motionConst[1,0],motionConst[2,0]]

			avgX += motionConst[0] - gndConst[0]
			avgY += motionConst[1] - gndConst[1]
			
			P = motionConst[2] - gndConst[2]
			avgPx += cos(P)
			avgPy += sin(P)
			
		avgX /= N
		avgY /= N
		avgPx /= N
		avgPy /= N
		avgP = atan2(avgPy, avgPx)

		stdX = 0.0
		stdY = 0.0
		stdP = 0.0

		Exx = 0.0
		Exy = 0.0
		Exp = 0.0
		Eyp = 0.0
		Eyy = 0.0
		Epp = 0.0
		
		for i in range(N):
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]			
			motionConst = self.motion_constraints[i][2]
			motionConst = [motionConst[0,0],motionConst[1,0],motionConst[2,0]]
			
			x = motionConst[0] - gndConst[0]
			y = motionConst[1] - gndConst[1]
			p = motionConst[2] - gndConst[2]
			
			Exx += (x - avgX) * (x - avgX)
			Exy += (x - avgX) * (y - avgY)
			Exp += (x - avgX) * (p - avgP)
			Eyp += (y - avgY) * (p - avgP)
			Eyy += (y - avgY) * (y - avgY)
			Epp += (p - avgP) * (p - avgP)
			
		#print "pre:", Exx, Eyy, Epp, Exy, Exp, Eyp
		Exx /= N
		Exy /= N
		Exp /= N
		Eyp /= N
		Eyy /= N
		Epp /= N
		
		E_motion = numpy.matrix([[Exx, Exy, Exp],
					[Exy, Eyy, Eyp],
					[Exp, Eyp, Epp]])
		#print "covar:", Exx, Eyy, Epp, Exy, Exp, Eyp	


		avgX = 0.0
		avgY = 0.0
		avgPx = 0.0
		avgPy = 0.0
		for i in range(N):
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]
			overlapConst = self.overlap_constraints[i][2]
			overlapConst = [overlapConst[0,0],overlapConst[1,0],overlapConst[2,0]]

			avgX += overlapConst[0] - gndConst[0]
			avgY += overlapConst[1] - gndConst[1]
			
			P = overlapConst[2] - gndConst[2]
			avgPx += cos(P)
			avgPy += sin(P)
			
		avgX /= N
		avgY /= N
		avgPx /= N
		avgPy /= N
		avgP = atan2(avgPy, avgPx)





		stdX = 0.0
		stdY = 0.0
		stdP = 0.0

		Exx = 0.0
		Exy = 0.0
		Exp = 0.0
		Eyp = 0.0
		Eyy = 0.0
		Epp = 0.0
		
		for i in range(N):
			
			gndConst = self.gnd_constraints[i][2]
			gndConst = [gndConst[0,0],gndConst[1,0],gndConst[2,0]]			
			overlapConst = self.overlap_constraints[i][2]
			overlapConst = [overlapConst[0,0],overlapConst[1,0],overlapConst[2,0]]
						
			x = motionConst[0] - gndConst[0]
			y = motionConst[1] - gndConst[1]
			p = motionConst[2] - gndConst[2]
			
			Exx += (x - avgX) * (x - avgX)
			Exy += (x - avgX) * (y - avgY)
			Exp += (x - avgX) * (p - avgP)
			Eyp += (y - avgY) * (p - avgP)
			Eyy += (y - avgY) * (y - avgY)
			Epp += (p - avgP) * (p - avgP)
			
		#print "pre:", Exx, Eyy, Epp, Exy, Exp, Eyp
		Exx /= N
		Exy /= N
		Exp /= N
		Eyp /= N
		Eyy /= N
		Epp /= N
		
		E_overlap = numpy.matrix([[Exx, Exy, Exp],
					[Exy, Eyy, Eyp],
					[Exp, Eyp, Epp]])

		return E_motion, E_overlap

	def makeOverlapHypothesis(self, transform, i, j ):	
	
		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
			
		curve1 = node1.getGPACurve()
		curve2 = node2.getGPACurve()

		offset = [transform[0,0],transform[1,0],transform[2,0]]


		result = gen_icp.motionICP(curve1.getUniformSamples(), curve2.getUniformSamples(), offset, plotIter = False, n1 = i, n2 = j)

		" update the pose with the motion estimation "
		offset = [result[0], result[1], result[2]]
		
		" add covariance from vector in guess  "
		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE = self.E_overlap
		
		return transform, covE

	def makeOverlapConstraint2(self, i, j ):

		hull1 = computeBareHull(self.nodeHash[i], sweep = False)
		hull1.append(hull1[0])
		medial1 = self.nodeHash[i].getMedialAxis(sweep = False)

		hull2 = computeBareHull(self.nodeHash[j], sweep = False)
		hull2.append(hull2[0])
		medial2 = self.nodeHash[j].getMedialAxis(sweep = False)


		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()		
		#hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = True)
		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = True)


		edge1 = medial1[0:2]
		edge2 = medial1[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]

		
		interPoints = []
		for k in range(len(hull1)-1):
			hullEdge = [hull1[k],hull1[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		for k in range(len(hull1)-1):
			hullEdge = [hull1[k],hull1[k+1]]
			isIntersect2, point2 = Intersect(edge2, hullEdge)
			if isIntersect2:
				interPoints.append(point2)
				break
			
		medial1 = medial1[1:-2]
		if isIntersect1:
			medial1.insert(0, point1)
		if isIntersect2:
			medial1.append(point2)


		edge1 = medial2[0:2]
		edge2 = medial2[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]

		
		interPoints = []
		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect2, point2 = Intersect(edge2, hullEdge)
			if isIntersect2:
				interPoints.append(point2)
				break
			
		medial2 = medial2[1:-2]
		if isIntersect1:
			medial2.insert(0, point1)
		if isIntersect2:
			medial2.append(point2)
		
		curve1 = SplineFit(medial1, smooth=0.1)
		curve2 = SplineFit(medial2, smooth=0.1)
		
		motionT, motionCov = self.makeMotionConstraint(i,j)
		travelDelta = motionT[0,0]

		if travelDelta < 0.3:
			offset = estimateMotion.makeGuess2(curve1, curve2, 0.0, posture1, posture2, forward = True)
		else:
			offset = estimateMotion.makeGuess2(curve1, curve2, 0.0, posture1, posture2, forward = False)
			
		result = gen_icp.motionICP(curve1.getUniformSamples(), curve2.getUniformSamples(), offset, plotIter = False, n1 = i, n2 = j)
		#result = gen_icp.motionICP(curve1.getUniformSamples(), curve2.getUniformSamples(), offset, plotIter = True, n1 = i, n2 = j)

		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE =  self.E_overlap
		
		print "making overlap constraint:", result[0], result[1], result[2]
		
		return transform, covE				

	def computeMedialError(self, i, j, offset):


		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()		
		hull1, medial1 = computeHullAxis(i, node1, tailCutOff = True)
		hull2, medial2 = computeHullAxis(j, node2, tailCutOff = True)

		"""
		if self.nodeHash[i].isBowtie:			
			hull1 = computeBareHull(self.nodeHash[i], sweep = False, static = True)
			hull1.append(hull1[0])
			medial1 = self.nodeHash[i].getStaticMedialAxis()
		else:
			hull1 = computeBareHull(self.nodeHash[i], sweep = False)
			hull1.append(hull1[0])
			medial1 = self.nodeHash[i].getMedialAxis(sweep = False)

		if self.nodeHash[j].isBowtie:			
			hull2 = computeBareHull(self.nodeHash[j], sweep = False, static = True)
			hull2.append(hull2[0])
			medial2 = self.nodeHash[j].getStaticMedialAxis()

		else:
			hull2 = computeBareHull(self.nodeHash[j], sweep = False)
			hull2.append(hull2[0])
			medial2 = self.nodeHash[j].getMedialAxis(sweep = False)

		" take the long length segments at tips of medial axis"
		edge1 = medial1[0:2]
		edge2 = medial1[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]

		" find the intersection points with the hull "
		interPoints = []
		for k in range(len(hull1)-1):
			hullEdge = [hull1[k],hull1[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		for k in range(len(hull1)-1):
			hullEdge = [hull1[k],hull1[k+1]]
			isIntersect2, point2 = Intersect(edge2, hullEdge)
			if isIntersect2:
				interPoints.append(point2)
				break

		" replace the extended edges with a termination point at the hull edge "			
		medial1 = medial1[1:-2]
		if isIntersect1:
			medial1.insert(0, point1)
		if isIntersect2:
			medial1.append(point2)

		TAILDIST = 0.5

		" cut off the tail of the non-sweeping side "
		if i % 2 == 0:
			termPoint = medial1[-1]
			for k in range(len(medial1)):
				candPoint = medial1[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
				
			medial1 = medial1[:-k-1]
		else:
			termPoint = medial1[0]
			for k in range(len(medial1)):
				candPoint = medial1[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial1 = medial1[k:]


		" take the long length segments at tips of medial axis"
		edge1 = medial2[0:2]
		edge2 = medial2[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]

		
		" find the intersection points with the hull "
		interPoints = []
		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		for k in range(len(hull2)-1):
			hullEdge = [hull2[k],hull2[k+1]]
			isIntersect2, point2 = Intersect(edge2, hullEdge)
			if isIntersect2:
				interPoints.append(point2)
				break
		
		" replace the extended edges with a termination point at the hull edge "			
		medial2 = medial2[1:-2]
		if isIntersect1:
			medial2.insert(0, point1)
		if isIntersect2:
			medial2.append(point2)
		

		" cut off the tail of the non-sweeping side "

		if j % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]

		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
		"""
		
		
		estPose1 = node1.getGlobalGPACPose()		
		estPose2 = node2.getGlobalGPACPose()
		
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		
		
		#points1 = medialSpline1.getUniformSamples()
		#points2 = medialSpline2.getUniformSamples()

		#offset = [transform[0,0],transform[1,0],transform[2,0]]
		#		
		#points2_offset = []
		#for p in points2:
		#	result = gen_icp.dispOffset(p, offset)		
		#	points2_offset.append(result)


		points1 = gen_icp.addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
		points2 = gen_icp.addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)
	
		" transform pose 2 by initial offset guess "	
		" transform the new pose "
		points2_offset = []
		for p in points2:
			result = gen_icp.dispPoint(p, offset)		
			points2_offset.append(result)


		#costThresh = 0.004
		minMatchDist = 2.0
		#lastCost = 1e100
		
		poly1 = []		
		for p in points1:
			poly1.append([p[0],p[1]])	
		
		" find the matching pairs "
		match_pairs = []
		
		" transformed points without associated covariance "
		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])	
		
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
		radius1, center1 = gen_icp.computeEnclosingCircle(points1)
	
		for i in range(len(points2_offset)):
			p_2 = poly2[i]

			if gen_icp.isInCircle(p_2, radius1, center1):

				" for every transformed point of A, find it's closest neighbor in B "
				p_1, minDist = gen_icp.findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2_offset[i][2]
					C1 = p_1[2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i],p_1,C2,C1])

		for i in range(len(points1)):
			p_1 = poly1[i]
	
			if gen_icp.isInCircle(p_1, radius2, center2):
		
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2_offset[i_2][2]
					
					C1 = points1[i][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i_2],points1[i],C2,C1])

		vals = []
		sum1 = 0.0
		for pair in match_pairs:
	
			a = pair[0]
			b = pair[1]
			Ca = pair[2]
			Cb = pair[3]
	
			ax = a[0]
			ay = a[1]		
			bx = b[0]
			by = b[1]
	
			c11 = Ca[0][0]
			c12 = Ca[0][1]
			c21 = Ca[1][0]
			c22 = Ca[1][1]
					
			b11 = Cb[0][0]
			b12 = Cb[0][1]
			b21 = Cb[1][0]
			b22 = Cb[1][1]	
		
			val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
						
			vals.append(val)
			sum1 += val


		return sum1


	def makeGlobalMedialOverlapConstraint(self, nodeID, globalPath, globalJunctionPoint, globalDeparturePoint):

		" compute the medial axis for each pose "
		
		node1 = self.nodeHash[nodeID]
		posture1 = node1.getStableGPACPosture()
		hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = True)
		
		globalPathReverse = deepcopy(globalPath)
		globalPathReverse.reverse()
		

		
		estPose1 = node1.getGlobalGPACPose()		
		poseOrigin = Pose(estPose1)
		
		globalMedial = []
		for p in medial1:
			globalMedial.append(poseOrigin.convertLocalToGlobal(p))
		
		medialSpline1 = SplineFit(globalMedial, smooth=0.1)
		globalSpline = SplineFit(globalPath, smooth=0.1)
		globalSplineReverse = SplineFit(globalPathReverse, smooth=0.1)


		overlapMatch = []
		angleSum1 = 0.0
		angleSum2 = 0.0
		for i in range(0,len(globalMedial)):
			p_1 = globalMedial[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(globalPath, p_1)
			if minDist < 0.5:
				overlapMatch.append((i,j,minDist))

				pathU1 = medialSpline1.findU(p_1)	
				pathU2 = globalSpline.findU(p_2)	
				pathU2_R = globalSplineReverse.findU(p_2)	

				pathVec1 = medialSpline1.getUVector(pathU1)
				pathVec2 = globalSpline.getUVector(pathU2)
				pathVec2_R = globalSplineReverse.getUVector(pathU2_R)

				#path1Mag = sqrt(pathVec1[0]*pathVec1[0] + pathVec1[1]*pathVec1[1])
				#path2Mag = sqrt(pathVec2[0]*pathVec2[0] + pathVec2[1]*pathVec2[1])
				#path2MagR = sqrt(pathVec2_R[0]*pathVec2_R[0] + pathVec2_R[1]*pathVec2_R[1])

				val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
				if val1 > 1.0:
					val1 = 1.0
				elif val1 < -1.0:
					val1 = -1.0
				ang1 = acos(val1)
				
				val2 = pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1]
				if val2 > 1.0:
					val2 = 1.0
				elif val2 < -1.0:
					val2 = -1.0
				ang2 = acos(val2)
				

				#ang1 = acos(pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1])
				#ang2 = acos(pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1])
		
				angleSum1 += ang1
				angleSum2 += ang2
		
		" select global path orientation based on which has the smallest angle between tangent vectors "
		print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
		if angleSum1 > angleSum2:
			orientedGlobalPath = globalPathReverse
		else:
			orientedGlobalPath = globalPath
			
		globalSpline = SplineFit(orientedGlobalPath, smooth=0.1)
		medialSpline1 = SplineFit(medial1, smooth=0.1)


		globalSamples = globalSpline.getUniformSamples(spacing = 0.04)
		medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
		
		globalMedialSamples = []
		for p in medialSamples:
			result = poseOrigin.convertLocalOffsetToGlobal(p)	
			globalMedialSamples.append(result)
		
		
		globalVar = []
		medialVar = []
		
		" compute the local variance of the angle "
		VAR_WIDTH = 40
		for i in range(len(globalSamples)):
			
			lowK = i - VAR_WIDTH/2
			if lowK < 0:
				lowK = 0
				
			highK = i + VAR_WIDTH/2
			if highK >= len(globalSamples):
				highK = len(globalSamples)-1
			
			localSamp = []
			for k in range(lowK, highK+1):
				localSamp.append(globalSamples[k][2])
			
			sum = 0.0
			for val in localSamp:
				sum += val
				
			meanSamp = sum / float(len(localSamp))
			

			sum = 0
			for k in range(len(localSamp)):
				sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
		
			varSamp = sum / float(len(localSamp))
			
			globalVar.append((meanSamp, varSamp))		

		for i in range(len(globalMedialSamples)):
			
			lowK = i - VAR_WIDTH/2
			if lowK < 0:
				lowK = 0
				
			highK = i + VAR_WIDTH/2
			if highK >= len(globalMedialSamples):
				highK = len(globalMedialSamples)-1
			
			localSamp = []
			for k in range(lowK, highK+1):
				localSamp.append(globalMedialSamples[k][2])
			
			sum = 0.0
			for val in localSamp:
				sum += val
				
			meanSamp = sum / float(len(localSamp))
			

			sum = 0
			for k in range(len(localSamp)):
				sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
		
			varSamp = sum / float(len(localSamp))
			
			medialVar.append((meanSamp, varSamp))		


		" now lets find closest points and save their local variances "			
		closestPairs = []
		TERM_DIST = 20
		
		for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
			pG = globalSamples[i]
			minDist = 1e100
			minJ = -1
			for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
				pM = globalMedialSamples[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minJ = j
					
			if minDist < 0.1:
				closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

		for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
			pM = globalMedialSamples[j]
			minDist = 1e100
			minI = -1
			for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
				pG = globalSamples[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minI = i
					
			if minDist < 0.1:
				closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))

		" remove duplicates "
		closestPairs = list(set(closestPairs))
		#closestPairs = list(closestPairs)
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(5,6))
		#s = sorted(student_objects, key=attrgetter('age'))     # sort on secondary key
		#sorted(s, key=attrgetter('grade'), reverse=True)  
		#closestPairs.sort()

		print len(closestPairs), "closest pairs"
		#for val in closestPairs:
		#	print val


		"""
		overlapMatch = []
		for i in range(0,len(globalMedial)):
			p_1 = globalMedial[i]
			p_2, j, minDist = gen_icp.findClosestPointInA(orientedGlobalPath, p_1)
			if minDist < 0.5:
				overlapMatch.append((i,j,minDist))

				pathU1 = medialSpline1.findU(p_1)	
				pathU2 = globalSpline.findU(p_2)	
				
				pathVec1 = medialSpline1.getUVector(pathU1)
				pathVec2 = globalSpline.getUVector(pathU2)

				val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
				if val1 > 1.0:
					val1 = 1.0
				elif val1 < -1.0:
					val1 = -1.0
				ang1 = acos(val1)
		"""
		if len(closestPairs) > 0:
			originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
			originU1 = globalSpline.findU(globalSamples[closestPairs[0][0]])

		else:
			localDeparturePoint = poseOrigin.convertGlobalToLocal(globalDeparturePoint)		
			originU2 = medialSpline1.findU(localDeparturePoint)	
			originU1 = globalSpline.findU(globalJunctionPoint)
		
		u2 = originU2
		u1 = originU1
		angGuess = 0.0
		
		resultPose, lastCost = gen_icp.globalOverlapICP([u1, u2, angGuess], orientedGlobalPath, hull1, medial1,  plotIter = False, n1 = nodeID, n2 = 0)

		print "estimating branch pose", nodeID, "at",  resultPose[0], resultPose[1], resultPose[2]

		return resultPose, lastCost

	def makeMedialOverlapConstraint(self, i, j, isMove = True, isForward = True, inPlace = False, uRange = 1.5 ):

		print "recomputing hulls and medial axis"
		" compute the medial axis for each pose "
		
			

		#self.nodeHash[i].computeStaticAlphaBoundary()
		#print "static alpha computed"

		#hull2 = computeBareHull(self.nodeHash[j], sweep = False)
		#hull2 = computeBareHull(self.nodeHash[j], sweep = False, static = True)
		#hull2.append(hull2[0])
		#medial2 = self.nodeHash[j].getMedialAxis(sweep = False)
		#medial2 = self.nodeHash[j].getStaticMedialAxis()

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()
		hull1, medial1 = computeHullAxis(i, node1, tailCutOff = True)
		hull2, medial2 = computeHullAxis(j, node2, tailCutOff = True)

		estPose1 = node1.getGlobalGPACPose()		
		estPose2 = node2.getGlobalGPACPose()
		
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medialSpline2 = SplineFit(medial2, smooth=0.1)


		#samples = scipy.arange(0.0,1.0,0.01)

		#originU1 = medialSpline2.findU(node1.rootPose)	
		#originU2 = medialSpline2.findU(node2.rootPose)	
		originU1 = medialSpline1.findU([0.0,0.0])	
		originU2 = medialSpline2.findU([0.0,0.0])	

		#medialSpline1.getUOfDist(originU1, dist)

		if inPlace:
			originU1 = 0.6
			originU2 = 0.4
			u2 = originU2
			print "computed u2 =", u2, "from originU2 =", originU2

		elif isMove:
			
			if isForward:
				
				if len(node2.frontProbeError) > 0:
	
					frontSum = 0.0
					frontProbeError = node2.frontProbeError
					for n in frontProbeError:
						frontSum += n
					foreAvg = frontSum / len(frontProbeError)
									
					if foreAvg >= 1.4:
						u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
					else:	
						u2 = medialSpline2.getUOfDist(originU2, 0.3, distIter = 0.001)
				else:	
					u2 = medialSpline2.getUOfDist(originU2, 0.3, distIter = 0.001)
					
			else:
				u2 = medialSpline2.getUOfDist(originU2, -0.3, distIter = 0.001)
				#u2 = 0.4
			print "computed u2 =", u2, "from originU2 =", originU2
			 
		else:
			poseOrigin = Pose(estPose1)
			offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
			
			points2 = medialSpline2.getUniformSamples()
			p_1 = medialSpline1.getU(0.5)
			
			points2_offset = []
			for p in points2:
				result = gen_icp.dispOffset(p, offset)		
				points2_offset.append(result)
	
			p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
	
			u2 = medialSpline2.findU(points2[i_2])	
			
			if u2 > 0.9 or u2 < 0.1:
				raise

		#print "u2 =", u2


		#motionT, motionCov = self.makeMotionConstraint(i,j)
		#travelDelta = motionT[0,0]
		
		u1 = originU1
		#u1 = 0.5
		#u2 = 0.6
		angGuess = 0.0
	
		" create the ground constraints "
		gndGPAC1Pose = node1.getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)
		gndGPAC2Pose = node2.getGndGlobalGPACPose()
		gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
		
		#supportLine = self.getTopology(nodes = range(self.numNodes-2))

		#result = gen_icp.overlapICP(estPose1, [u1, u2, angGuess], hull1, hull2, medial1, medial2, node1.rootPose, node2.rootPose, plotIter = True, n1 = i, n2 = j)
		result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = i, n2 = j, uRange = uRange)
		#result, hist = gen_icp.overlapICP2(estPose1, [u1, u2, angGuess], hull1, hull2, medial1, medial2, supportLine, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = True, n1 = i, n2 = j)

		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE =  self.E_overlap
		
		print "making overlap constraint:", result[0], result[1], result[2]

		return transform, covE, hist


	def makeOverlapConstraint(self, i, j ):	

		motionT, motionCov = self.makeMotionConstraint(i,j)
		travelDelta = motionT[0,0]

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		
		#curve1 = node1.getGPACurve()
		#curve2 = node2.getGPACurve()
		#posture1 = node1.getGPACPosture()
		#posture2 = node2.getGPACPosture()
		
		curve1 = node1.getStableGPACurve()
		curve2 = node2.getStableGPACurve()
		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()


		#(points1, points2, offset, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):
		if travelDelta < 0.3:
			offset = estimateMotion.makeGuess2(curve1, curve2, 0.0, posture1, posture2, forward = True)
		else:
			offset = estimateMotion.makeGuess2(curve1, curve2, 0.0, posture1, posture2, forward = False)
			
		result = gen_icp.motionICP(curve1.getUniformSamples(), curve2.getUniformSamples(), offset, plotIter = False, n1 = i, n2 = j)

		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE =  self.E_overlap
		
		print "making overlap constraint:", result[0], result[1], result[2]
		
		return transform, covE
	
	def makeMotionConstraint(self, i, j):
		
		#print "making motion constraint", i, j

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		pose1 = node1.getGlobalGPACPose()
		pose2 = node2.getGlobalGPACPose()
		
		profile1 = Pose(pose1)
		
		localOffset = profile1.convertGlobalPoseToLocal(pose2)

		xT = localOffset[0]
		yT = localOffset[1]
		pT = localOffset[2]
		
		transform = matrix([[xT], [yT], [pT]])
		covE = self.E_motion
		
		return transform, covE
				
	
	def makeCornerConstraint(self, n1, n2, point1, point2, ang1, ang2, hull1, hull2, sweepHull1, sweepHull2):
		
		" TUNE: estimated step distance from empirical results "
		STEP_DIST = 0.145

		" initial guess for x, y, theta parameters "
		firstGuess = 0.0

		" TUNE ME:  threshold cost difference between iterations to determine if converged "
		costThresh = 0.1
	
		" TUNE ME:   minimum match distance before point is discarded from consideration "
		minMatchDist = 0.5
	
		" plot the best fit at each iteration of the algorithm? "
		#plotIteration = True
		plotIteration = False
		
		" Extract the data from the files and put them into arrays "		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]		
		
		" FIXME:  deal with relative poses from the motion constraints "
		estPose1 = node1.getGlobalGPACPose()
		estPose2 = node2.getGlobalGPACPose()
		
		radius, center = gen_icp.computeEnclosingCircle(hull1)
		circle1 = [radius, center]

		#def cornerICP(estPose1, angle, point1, point2, ang1, ang2, hull1, hull2, pastCircles, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):
		
		offset, cost = gen_icp.cornerICP(estPose1, firstGuess, point1, point2, ang1, ang2, hull1, hull2, sweepHull1, sweepHull2, [circle1], costThresh, minMatchDist, plotIteration, n1, n2)
		#print "made corner constraint: %d -> %d" %(n1,n2), "cost =", cost
		covar = deepcopy(self.E_corner)
		
		
		
		" 1) get closest points where GPA curves overlap "
		" 2) check relative orientation at those curves "
		" 3) if orientation is around 180 degrees, reject "
		#print "checking orientation error"
		
		isOkay, angDiff2, oriDiff = self.checkOrientationError(offset, n1, n2)
		
		#isOkay, angDiff = self.checkOrientationError(offset, n1, n2)
		if isOkay:
			return offset, covar, cost, angDiff2, oriDiff
		else:
			return offset, covar, 1e100, angDiff2, oriDiff
			
			 
		" 1) need to vary angle until we have wall fitting "
		" 2) compute alpha hulls and then perform maximum edge fitting, not overlap constraint"
		" 3) directionality is not important here because we could be walking the opposite direction through a junction "

	def computeConstraintOrientation(self, offset, n1, n2):

		def dispOffset(p, offset):
			xd = offset[0]
			yd = offset[1]
			theta = offset[2]
		
			px = p[0]
			py = p[1]
			pa = p[2]
			
			newAngle = normalizeAngle(pa + theta)
			
			p_off = [px*math.cos(theta) - py*math.sin(theta) + xd, px*math.sin(theta) + py*math.cos(theta) + yd, newAngle]
			
			return p_off
		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]
								
		curve1 = node1.getStableGPACurve()
		curve2 = node2.getStableGPACurve()

		points1 = curve1.getUniformSamples()
		points2 = curve2.getUniformSamples()

		" transform the new pose "
		points2_offset = []
		for p in points2:
			result = dispOffset(p, offset)		
			points2_offset.append(result)
		
		poly1 = []
		for p in points1:
			poly1.append([p[0],p[1]])	

		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])	
		
		" find the matching pairs "
		match_pairs = []
			
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
		radius1, center1 = gen_icp.computeEnclosingCircle(points1)

		minMatchDist = 0.05
			
		for i in range(len(points2_offset)):
			p_2 = poly2[i]

			if gen_icp.isValid(p_2, radius1, center1, poly1):

				" for every transformed point of A, find it's closest neighbor in B "
				p_1, minDist = gen_icp.findClosestPointInHull1(points1, p_2, [0.0,0.0,0.0])
	
				if minDist <= minMatchDist:

					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i],p_1])

		for i in range(len(points1)):
			p_1 = poly1[i]
	
			if gen_icp.isValid(p_1, radius2, center2, poly2):
		
				#print "selected", b_p, "in circle", radiusA, centerA, "with distance"
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = gen_icp.findClosestPointInHull2(points2_offset, p_1)
	
				if minDist <= minMatchDist:
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i_2],points1[i]])
					
		" if we made no matches, default to regular orientation "
		if len(match_pairs) == 0:
			return True, 0.0

		distances = []
		angSum = 0
		for i in range(len(match_pairs)):
			pair = match_pairs[i]
			p1 = pair[0]
			p2 = pair[1]
			
			angDist = normalizeAngle(p1[2]-p2[2])
			distances.append(angDist)
			angSum += angDist

		
		angDiff = angSum / len(distances)


		if fabs(angDiff) >= 0.7:
			return False, angDiff
		else:
			return True, angDiff
	

	def checkOrientationError(self, offset, n1, n2):

		def dispOffset(p, offset):
			xd = offset[0]
			yd = offset[1]
			theta = offset[2]
		
			px = p[0]
			py = p[1]
			pa = p[2]
			
			newAngle = normalizeAngle(pa + theta)
			#print pa, repr(theta), repr(newAngle)
			
			p_off = [px*math.cos(theta) - py*math.sin(theta) + xd, px*math.sin(theta) + py*math.cos(theta) + yd, newAngle]
			
			return p_off

		#normalizeAngle(pose[2] - self.estPose[2])
		
		#print "offset:", repr(offset)
		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]
			
		estPose1 = node1.getGlobalGPACPose()
		estPose2 = node2.getGlobalGPACPose()
		
		poseProfile = Pose(estPose1)
		localOffset2 = poseProfile.convertGlobalPoseToLocal(estPose2)
		
		#globalAngDiff = normalizeAngle(estPose2[2] - estPose1[2])
			
		curve1 = node1.getStableGPACurve()
		curve2 = node2.getStableGPACurve()
		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()

		points1 = curve1.getUniformSamples()
		points2 = curve2.getUniformSamples()

		" transform the new pose "
		points2_offset = []
		for p in points2:
			result = dispOffset(p, offset)		
			points2_offset.append(result)

		points2_offset_global = []
		for p in points2:
			result = dispOffset(p, localOffset2)		
			points2_offset_global.append(result)
		
		poly1 = []
		for p in points1:
			poly1.append([p[0],p[1]])	

		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])	
		
		" find the matching pairs "
		match_pairs = []
		match_pairs2 = []
			
		" get the circles and radii "
		radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
		radius1, center1 = gen_icp.computeEnclosingCircle(points1)

		minMatchDist = 0.05
			
		for i in range(len(points2_offset)):
			p_2 = poly2[i]

			if gen_icp.isValid(p_2, radius1, center1, poly1):

				" for every transformed point of A, find it's closest neighbor in B "
				p_1, minDist = gen_icp.findClosestPointInHull1(points1, p_2, [0.0,0.0,0.0])
	
				if minDist <= minMatchDist:

					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i],p_1])
					match_pairs2.append([points2_offset_global[i],p_1])

		for i in range(len(points1)):
			p_1 = poly1[i]
	
			if gen_icp.isValid(p_1, radius2, center2, poly2):
		
				#print "selected", b_p, "in circle", radiusA, centerA, "with distance"
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = gen_icp.findClosestPointInHull2(points2_offset, p_1)
	
				if minDist <= minMatchDist:
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2_offset[i_2],points1[i]])
					match_pairs2.append([points2_offset_global[i_2],points1[i]])
					
		if len(match_pairs) == 0:
			return False, 0.0, 0.0

		distances = []
		angSum = 0
		for i in range(len(match_pairs)):
			pair = match_pairs[i]
			p1 = pair[0]
			p2 = pair[1]
			
			dist = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
			#print p1, p2
			angDist = normalizeAngle(p1[2]-p2[2])
			distances.append(angDist)
			angSum += angDist

		
		angDiff = angSum / len(distances)
		oriDiff = normalizeAngle(offset[2] - localOffset2[2])

		if fabs(angDiff) >= 0.7 or fabs(oriDiff) >= 0.4:
			return False, angDiff, oriDiff
		else:
			return True, angDiff, oriDiff
		
	def makeSensorConstraint(self, transform, n1, n2, hull1, hull2):		

		" TUNE: estimated step distance from empirical results "
		STEP_DIST = 0.145

		" initial guess for x, y, theta parameters "
		#firstGuess = self.makeGuess(n1, n2, STEP_DIST)
		#print "guess =", firstGuess
		firstGuess = [0.0,0.0,0.0]

		" TUNE ME:  threshold cost difference between iterations to determine if converged "
		costThresh = 0.1
	
		" TUNE ME:   minimum match distance before point is discarded from consideration "
		minMatchDist = 2.0
		#minMatchDist = 0.5
	
		" plot the best fit at each iteration of the algorithm? "
		#plotIteration = True
		plotIteration = False
	
		#offset = [0.0,0.0,0.0]
	
		" Extract the data from the files and put them into arrays "
		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]		
		
		" FIXME:  deal with relative poses from the motion constraints "
		estPose1 = node1.getGlobalGPACPose()
		estPose2 = node2.getGlobalGPACPose()
		#estPose1 = node1.getEstPose()
		#estPose2 = node2.getEstPose()
		#curve1 = node1.getStableGPACurve()
		#curve2 = node2.getStableGPACurve()
			
		#posture1 = node1.getStableGPACPosture()
		#posture_trans = []
		#for p in posture1:
		#	posture_trans.append(gen_icp.dispOffset(p, currPose))	

		poly1 = []
		for p in hull1:
			poly1.append([p[0],p[1]])
		poly2 = []
		for p in hull2:
			poly2.append([p[0],p[1]])


		posture1_unstable = node1.getGPACPosture()
		posture1_stable = node1.getStableGPACPosture()
		posture2_unstable = node2.getGPACPosture()
		posture2_stable = node2.getStableGPACPosture()

		count1_unstable = 0
		count1_stable = 0
		for p in posture1_unstable:
			if point_inside_polygon(p[0],p[1],poly1):
				count1_unstable += 1
		for p in posture1_stable:
			if point_inside_polygon(p[0],p[1],poly1):
				count1_stable += 1
	
		posture1 = node1.getStableGPACPosture()
		if count1_stable < 38:
			if count1_unstable > count1_stable:
				posture1 = node1.getGPACPosture()

		count2_unstable = 0
		count2_stable = 0
		for p in posture2_unstable:
			if point_inside_polygon(p[0],p[1],poly2):
				count2_unstable += 1
		for p in posture2_stable:
			if point_inside_polygon(p[0],p[1],poly2):
				count2_stable += 1
	
		posture2 = node2.getStableGPACPosture()
		if count2_stable < 38:
			if count2_unstable > count2_stable:
				posture2 = node2.getGPACPosture()


		#uniform1 = posture1
		#uniform2 = posture2
		
		" convert hull points to GPAC coordinates before adding covariances "
		#localGPACPose1 = node1.getLocalGPACPose()
		#localGPACProfile1 = Pose(localGPACPose1)
		#localGPACPose2 = node2.getLocalGPACPose()
		#localGPACProfile2 = Pose(localGPACPose2)
		
		
		#curve1 = GPACurve(posture1)
		#curve2 = GPACurve(posture2)
		#uniform1 = curve1.getUniformSamples()
		#uniform2 = curve2.getUniformSamples()

		curve1 = StableCurve(posture1)
		curve2 = StableCurve(posture2)
		uniform1 = curve1.getPlot()
		uniform2 = curve2.getPlot()

	
		curve1 = StableCurve(posture1)
		curve2 = StableCurve(posture2)
		
		headPoint1, tailPoint1 = curve1.getTips()
		headPoint2, tailPoint2 = curve2.getTips()
		
		stablePoints1 = curve1.getPoints()
		stablePoints2 = curve2.getPoints()
		
		#headPoint1 = curve1.getU(0.0)
		#tailPoint1 = curve1.getU(1.0)
		#headPoint2 = curve2.getU(0.0)
		#tailPoint2 = curve2.getU(1.0)

		#uniform1 = curve1.getKnots()
		#uniform2 = curve2.getKnots()
		#uniform1 = curve1.getUniformSamples()
		#uniform2 = curve2.getUniformSamples()
		

		#hull_trans1 = []
		#for p in hull1:
		#	hull_trans1.append(gen_icp.dispPoint(p, estPose1))

		#hull_trans2 = []
		#for p in hull2:
			#hull_trans2.append(gen_icp.dispPoint(p, estPose2))

		#offset2 = [transform[0,0],transform[1,0],transform[2,0]]
		
		radius, center = gen_icp.computeEnclosingCircle(hull1)
		circle1 = [radius, center]
		
		isForward, angDiff = self.computeConstraintOrientation(firstGuess, n1, n2)

		
		offset, cost = gen_icp.shapeICP(estPose1, firstGuess, hull1, hull2, stablePoints1, stablePoints2, uniform1, uniform2, isForward, (headPoint1,tailPoint1,headPoint2,tailPoint2), [circle1], costThresh, minMatchDist, plotIteration, n1, n2)
		print "sensor constraint: %d -> %d" %(n1,n2), "cost =", cost
		covar = self.E_sensor
		
		return offset, covar, cost

	def makeShapeConstraintData2(self, transform, n1, n2, hull1, hull2, medial1, medial2):		

		" TUNE ME:  threshold cost difference between iterations to determine if converged "
		costThresh = 0.1
	
		" TUNE ME:   minimum match distance before point is discarded from consideration "
		minMatchDist = 2.0
	
		" plot the best fit at each iteration of the algorithm? "
		#plotIteration = True
		plotIteration = False
		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]		
		
		estPose1 = node1.getGlobalGPACPose()
		estPose2 = node2.getGlobalGPACPose()

		posture1_unstable = node1.getGPACPosture()
		posture1_stable = node1.getStableGPACPosture()
		posture2_unstable = node2.getGPACPosture()
		posture2_stable = node2.getStableGPACPosture()

		return (estPose1, estPose2, hull1, hull2, medial1, medial2, posture1_unstable, posture1_stable, posture2_unstable, posture1_stable, costThresh, minMatchDist, plotIteration, n1, n2)

	def makeShapeConstraintData(self, transform, n1, n2, hull1, hull2):		

		" TUNE: estimated step distance from empirical results "
		STEP_DIST = 0.145

		" initial guess for x, y, theta parameters "
		firstGuess = [0.0,0.0,0.0]

		" TUNE ME:  threshold cost difference between iterations to determine if converged "
		costThresh = 0.1
	
		" TUNE ME:   minimum match distance before point is discarded from consideration "
		minMatchDist = 2.0
	
		" plot the best fit at each iteration of the algorithm? "
		#plotIteration = True
		plotIteration = False
	
		" Extract the data from the files and put them into arrays "
		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]		
		
		" FIXME:  deal with relative poses from the motion constraints "
		estPose1 = node1.getGlobalGPACPose()
		estPose2 = node2.getGlobalGPACPose()

		#estPose1, estPose2
		#posture1_unstable, posture1_stable, posture2_unstable, posture2_stable


		poly1 = []
		for p in hull1:
			poly1.append([p[0],p[1]])
		poly2 = []
		for p in hull2:
			poly2.append([p[0],p[1]])


		posture1_unstable = node1.getGPACPosture()
		posture1_stable = node1.getStableGPACPosture()
		posture2_unstable = node2.getGPACPosture()
		posture2_stable = node2.getStableGPACPosture()

		count1_unstable = 0
		count1_stable = 0
		for p in posture1_unstable:
			if point_inside_polygon(p[0],p[1],poly1):
				count1_unstable += 1
		for p in posture1_stable:
			if point_inside_polygon(p[0],p[1],poly1):
				count1_stable += 1
	
		posture1 = posture1_stable
		if count1_stable < 38:
			if count1_unstable > count1_stable:
				posture1 = posture2_unstable

		count2_unstable = 0
		count2_stable = 0
		for p in posture2_unstable:
			if point_inside_polygon(p[0],p[1],poly2):
				count2_unstable += 1
		for p in posture2_stable:
			if point_inside_polygon(p[0],p[1],poly2):
				count2_stable += 1
	
		posture2 = posture2_stable
		if count2_stable < 38:
			if count2_unstable > count2_stable:
				posture2 = posture2_unstable


		curve1 = StableCurve(posture1)
		curve2 = StableCurve(posture2)
		uniform1 = curve1.getPlot()
		uniform2 = curve2.getPlot()

		headPoint1, tailPoint1 = curve1.getTips()
		headPoint2, tailPoint2 = curve2.getTips()
		
		stablePoints1 = curve1.getPoints()
		stablePoints2 = curve2.getPoints()
		
		radius, center = gen_icp.computeEnclosingCircle(hull1)
		circle1 = [radius, center]
		
		isForward, angDiff = self.computeConstraintOrientation(firstGuess, n1, n2)

		return (estPose1, firstGuess, hull1, hull2, stablePoints1, stablePoints2, uniform1, uniform2, isForward, (headPoint1,tailPoint1,headPoint2,tailPoint2), [circle1], costThresh, minMatchDist, plotIteration, n1, n2)

	def makeInPlaceConstraint(self, nodeID1, nodeID2):

		" compute the new posture "
		" 1. compute the GPAC pose "
		" 2. perform correction of GPAC pose "
		" 3. compute location of rootPose from corrected GPAC pose "
		
		" node1 is the front poke node "
		" nodes is the back poke node "

		node1 = self.nodeHash[nodeID1]
		node2 = self.nodeHash[nodeID2]
		
		#print "rootPose1", node1.rootPose
		#print "rootPose2", node2.rootPose
		#print "estPose1", node1.estPose
		#print "estPose2", node2.estPose

		#originPosture = node1.localPosture
		originPosture = node1.correctedPosture
		originCurve = StableCurve(originPosture)
		originForePose, originBackPose = originCurve.getPoses()
		
		originForeProfile = Pose(originForePose)
		originBackProfile = Pose(originBackPose)
		
		originForePosture = []
		originBackPosture = []
		for i in range(len(originPosture)):
			originForePosture.append(originForeProfile.convertGlobalPoseToLocal(originPosture[i]))
		for i in range(len(originPosture)):
			originBackPosture.append(originBackProfile.convertGlobalPoseToLocal(originPosture[i]))

		#newPosture = node2.localPosture
		newPosture = node2.correctedPosture
		#newPosture = []
		#for j in range(self.probe.numSegs-1):
		#	newPosture.append(self.probe.getJointPose([0.0,0.0,0.0], node1.rootNode, j))


		localPosture = newPosture
		localCurve = StableCurve(localPosture)

		localForePose, localBackPose = localCurve.getPoses()
		 
		localForeProfile = Pose(localForePose)
		localBackProfile = Pose(localBackPose)
		
		localForePosture = []
		localBackPosture = []
		for i in range(len(localPosture)):
			localForePosture.append(localForeProfile.convertGlobalPoseToLocal(localPosture[i]))
		for i in range(len(localPosture)):
			localBackPosture.append(localBackProfile.convertGlobalPoseToLocal(localPosture[i]))
	
		foreCost = 0.0
		for j in range(0,20):
			p1 = originForePosture[j]
			p2 = localForePosture[j]
			
			foreCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)		   
		
		backCost = 0.0
		for j in range(20,39):
			p1 = originBackPosture[j]
			p2 = localBackPosture[j]
			
			backCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)		   

		angle = 0.0

		correctedGPACPose1 = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		correctedGPACPose2 = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		correctedProfile1 = Pose(correctedGPACPose1)
		correctedProfile2 = Pose(correctedGPACPose2)

		localRootOffset3_1 = correctedProfile1.convertGlobalPoseToLocal([0.0,0.0,0.0])
		localRootOffset3_2 = correctedProfile2.convertGlobalPoseToLocal([0.0,0.0,0.0])

		offset1 = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3_1)
		offset2 = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3_2)
		
		" except the case where one node has not updated "
		if self.nodeHash[nodeID1].occMap.mapImage == 0 or self.nodeHash[nodeID2].occMap.mapImage == 0:
			cost1 = 1.0
			cost2 = 0.0
		else:
			cost1 = self.computeMedialError(nodeID1, nodeID2, offset1)
			cost2 = self.computeMedialError(nodeID1, nodeID2, offset2)
					
		#if foreCost > backCost:
		if cost1 < cost2:
			#plotPosture(originBackPosture, localBackPosture)
			correctedGPACPose = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		else:
			#plotPosture(originForePosture, localForePosture)					
			correctedGPACPose = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])

		correctedProfile = Pose(correctedGPACPose)

		" offset from 0,0 to newGPACPose "
		" rootPose from neg. offset from correctedGPACPose "
		localRootOffset3 = correctedProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		#localRootOffset3 = correctedProfile.convertGlobalPoseToLocal(node2.rootPose)
		
		
		#if foreCost > backCost:
		#if foreCost > backCost:
		if cost1 < cost2:
			offset = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3)
		else:
			offset = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3)
		
		transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
		covE = self.E_inplace
		
		return transform, covE

	def addInPlaceNode(self, newNode):
		
		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		#self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		#self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		#self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		#self.c_hulls.append(computeHull(self.nodeHash[nodeID], static = True))
		
		
		#self.currNode = newNode
		
		self.nodeHash[nodeID] = self.currNode
		if nodeID > 0:
	
			node1 = self.nodeHash[nodeID-1]
			node2 = self.nodeHash[nodeID]
				
			transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
			
			self.inplace_constraints.append([nodeID-1, nodeID, transform, covE])	


			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(self.nodeHash[nodeID-1].getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			self.nodeHash[nodeID].setGPACPose(currPose)
						
		self.saveConstraints()
		
		#self.numNodes += 1
					
	def addNode(self, newNode):
		
		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
				
	def addInitNode(self, newNode, estPose):
		
		self.currNode = newNode

		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		#self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		#self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		#self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		#self.c_hulls.append(computeHull(self.nodeHash[nodeID], static = True))


		if nodeID > 0:

			transform, covE = self.makeOverlapConstraint(0, nodeID)
			self.overlap_constraints.append([0, nodeID, transform, covE])

			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(self.nodeHash[0].getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			self.nodeHash[nodeID].setGPACPose(currPose)

			self.edgeHash[(0,nodeID)] = [transform, covE]
		
		self.saveConstraints()
		
		#self.numNodes += 1

	
	def correctFullAlpha(self):

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
		
			node1 = self.nodeHash[i]
				
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

			node1 = self.nodeHash[i]
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

	def findBestOverlap(self, paths, targetNode):

		SENSOR_RADIUS = 0.0
		DIST_THRESHOLD = 6.0
		
		" select pairs to attempt a sensor constraint "
		pair_candidates = []

		path_set = paths[targetNode]
		
		ind = range(targetNode+2,self.numNodes)

		for j in ind:	
			loc2 = path_set[j][0]	
			
			c1 = matrix([[0.0],[0.0]],dtype=float)
			c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
			E_param = path_set[j][1]
			E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
			
			dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
			
			print "distance:", targetNode, j, dist

			if dist <= DIST_THRESHOLD:
				pair_candidates.append([dist,targetNode,j,loc2])
			
		ind = range(0,targetNode-1)
		ind.reverse()

		
		for j in ind:
			
			loc2 = path_set[j][0]
			
			c1 = matrix([[0.0],[0.0]],dtype=float)
			c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
			E_param = path_set[j][1]
			E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
			
			dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
			
			print "distance:", targetNode, j, dist
			if dist <= DIST_THRESHOLD:
				pair_candidates.append([dist,targetNode,j,loc2])

		if len(pair_candidates) == 0:
			return []
		
		minDist = 1e100
		bestCandidate = 0
		for i in range(len(pair_candidates)):
			if dist <= minDist:
				minDist = dist
				bestCandidate = i

		" make hypothesized constraint out of the pair "
		p = pair_candidates[bestCandidate]
		transform = p[3]
		n1 = p[1]
		n2 = p[2]
		transform, covar = self.makeOverlapHypothesis(transform, n1, n2 )
		return [p[1],p[2],transform,covar]

	def findOverlapCandidates(self, paths, targetNode = -1):

		SENSOR_RADIUS = 0.0
		#DIST_THRESHOLD = 6.0
		DIST_THRESHOLD = 0.8
		
		" select pairs to attempt a sensor constraint "
		pair_candidates = []

		if targetNode >= 0:

			path_set = paths[targetNode]
			
			#print "path_set:", path_set
			
			ind = range(targetNode+2,self.numNodes)

			for j in ind:	
				loc2 = path_set[j][0]	
				
				c1 = matrix([[0.0],[0.0]],dtype=float)
				c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
				
				dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
				
				print "distance:", targetNode, j, dist

				if dist <= DIST_THRESHOLD:
					pair_candidates.append([dist,targetNode,j,loc2])
				
			ind = range(0,targetNode-1)
			ind.reverse()
			
			for j in ind:
				
				loc2 = path_set[j][0]
				
				c1 = matrix([[0.0],[0.0]],dtype=float)
				c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
				
				dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
				
				print "distance:", targetNode, j, dist
				if dist <= DIST_THRESHOLD:
					pair_candidates.append([dist,targetNode,j,loc2])

		
		else:
			
			for i in range(self.numNodes):
				path_set = paths[i]
				
				ind = range(i+2,self.numNodes)
				#ind = range(i+1,self.numNodes)
				#print i, self.numNodes, ind
	
				for j in ind:	
					loc2 = path_set[j][0]	
					
					c1 = matrix([[0.0],[0.0]],dtype=float)
					c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
					E_param = path_set[j][1]
					E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
					
					dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
					
					print "distance:", i, j, dist
	
					if dist <= DIST_THRESHOLD and abs(i-j) < 5:
						pair_candidates.append([dist,i,j,loc2])
					
				ind = range(0,i-1)
				#ind = range(0,i)
				ind.reverse()
				
				for j in ind:
					
					loc2 = path_set[j][0]
					
					c1 = matrix([[0.0],[0.0]],dtype=float)
					c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
					E_param = path_set[j][1]
					E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
					
					dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
					
					print "distance:", i, j, dist
					if dist <= DIST_THRESHOLD and abs(i-j) < 5:
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

		pair_unique.sort()
		
		#print "UNIQUE PAIRS"
		for p in pair_unique:
			i = p[1]
			j = p[2]
			
			
			loc2 = paths[i][j][0]
			c1 = matrix([[0.0],[0.0]],dtype=float)
			c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
		
			E_param = paths[i][j][1]
			E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
			
			dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
			
			#print "distance:", i, j, dist

		return pair_unique

	def findCornerCandidates(self, paths, nodeID = None):

		SENSOR_RADIUS = 0.0
		DIST_THRESHOLD = 2.0
		ANG_THRESHOLD = 2.0
		CORNER_DIFF = pi/4.
		
		" select pairs to attempt a sensor constraint "
		pair_candidates = []

		if nodeID != None:
			print "searching for corner pairs from node", nodeID
			print "node has", len(self.nodeHash[nodeID].cornerCandidates), "corners"
			path_set = paths[nodeID]
			
			ind = range(nodeID+2,self.numNodes)

			for j in ind:	
				loc2 = path_set[j][0]
				
				" c1 and c2 are centroids of sweep regions "
				node1 = self.nodeHash[nodeID]
				node2 = self.nodeHash[j]
				corners1 = deepcopy(node1.cornerCandidates)
				corners2 = deepcopy(node2.cornerCandidates)

				for m in range(len(corners1)):
					point1, cornerAngle1, ang1 = corners1[m]
						
					for n in range(len(corners2)):
						point2, cornerAngle2, ang2 = corners2[n]

						if fabs(cornerAngle2-cornerAngle1) < CORNER_DIFF:
	
							sweepOffset1 = point1
							sweepOffset2 = point2				
	
							node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
							sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)

							oriPose = node2Profile.convertLocalOffsetToGlobal([0.0,0.0,ang2])
							ang2 = oriPose[2]

							
							c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
							c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
							E_param = path_set[j][1]
							E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
							
							dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
							angDist = bayes.mahabAngDist(ang1, ang2, E_param[2,2])

							#print "ang1,ang2,E:", ang1, ang2, E_param[2,2]
							#print "c1,c2, dist =", [c1[0,0],c1[1,0]], [c2[0,0],c2[1,0]]
							#print "corner:", nodeID, j, dist, angDist
			
							#if dist <= DIST_THRESHOLD and abs(i-j) < 5:
							#if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD and abs(i-j) < 5:
							if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD:
								#print "distance:", i, j, dist								
								#print "angDist:", i, j, angDist

								pair_candidates.append([dist,nodeID,j,loc2, angDist, point1, point2, ang1, ang2, m, n])
				
			ind = range(0,nodeID-1)
			ind.reverse()
			
			for j in ind:
				
				loc2 = path_set[j][0]

				
				" c1 and c2 are centroids of sweep regions "
				node1 = self.nodeHash[nodeID]
				node2 = self.nodeHash[j]
				corners1 = deepcopy(node1.cornerCandidates)
				corners2 = deepcopy(node2.cornerCandidates)
				
				#print nodeID, j, "corners:", len(corners1), len(corners2)

				for m in range(len(corners1)):
					point1, cornerAngle1, ang1 = corners1[m]

					for n in range(len(corners2)):
						point2, cornerAngle2, ang2 = corners2[n]

						if fabs(cornerAngle2-cornerAngle1) < CORNER_DIFF:
	
							sweepOffset1 = point1						
							sweepOffset2 = point2				
	
							node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
							sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
							
							oriPose = node2Profile.convertLocalOffsetToGlobal([0.0,0.0,ang2])
							ang2 = oriPose[2]

							
							c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
							c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
							
							E_param = path_set[j][1]
							E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
							
							dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
							angDist = bayes.mahabAngDist(ang1, ang2, E_param[2,2])
							#print "ang1,ang2,E:", ang1, ang2, E_param[2,2]
							#print "c1,c2, dist =", [c1[0,0],c1[1,0]], [c2[0,0],c2[1,0]]
							#print "corner:", nodeID, j, dist, angDist
	
							#if dist <= DIST_THRESHOLD and abs(i-j) < 5:
							#if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD and abs(i-j) < 5:
							if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD:
								#print "distance:", i, j, dist
								#print "angDist:", i, j, angDist						
	
								pair_candidates.append([dist,nodeID,j,loc2, angDist, point1, point2, ang1, ang2, m, n])
		else:
			
			for i in range(self.numNodes):
				path_set = paths[i]
				
				ind = range(i+2,self.numNodes)
	
				for j in ind:	
					loc2 = path_set[j][0]
					
					" c1 and c2 are centroids of sweep regions "
					node1 = self.nodeHash[i]
					node2 = self.nodeHash[j]
					corners1 = deepcopy(node1.cornerCandidates)
					corners2 = deepcopy(node2.cornerCandidates)
	
					for m in range(len(corners1)):
						point1, cornerAngle1, ang1 = corners1[m]
							
						for n in range(len(corners2)):
							point2, cornerAngle2, ang2 = corners2[n]
	
							if fabs(cornerAngle2-cornerAngle1) < CORNER_DIFF:
		
								sweepOffset1 = point1
								sweepOffset2 = point2				
		
								node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
								sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
	
								oriPose = node2Profile.convertLocalOffsetToGlobal([0.0,0.0,ang2])
								ang2 = oriPose[2]
	
								
								c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
								c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
								E_param = path_set[j][1]
								E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
								
								dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
								angDist = bayes.mahabAngDist(ang1, ang2, E_param[2,2])
	
								
				
								#if dist <= DIST_THRESHOLD and abs(i-j) < 5:
								#if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD and abs(i-j) < 5:
								if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD:
									#print "distance:", i, j, dist								
									#print "angDist:", i, j, angDist
	
									pair_candidates.append([dist,i,j,loc2, angDist, point1, point2, ang1, ang2, m, n])
					
				ind = range(0,i-1)
				ind.reverse()
				
				for j in ind:
					
					loc2 = path_set[j][0]
	
					
					" c1 and c2 are centroids of sweep regions "
					node1 = self.nodeHash[i]
					node2 = self.nodeHash[j]
					corners1 = deepcopy(node1.cornerCandidates)
					corners2 = deepcopy(node2.cornerCandidates)
	
					for m in range(len(corners1)):
						point1, cornerAngle1, ang1 = corners1[m]
	
						for n in range(len(corners2)):
							point2, cornerAngle2, ang2 = corners2[n]
	
							if fabs(cornerAngle2-cornerAngle1) < CORNER_DIFF:
		
								sweepOffset1 = point1						
								sweepOffset2 = point2				
		
								node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
								sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
								
								oriPose = node2Profile.convertLocalOffsetToGlobal([0.0,0.0,ang2])
								ang2 = oriPose[2]
	
								
								c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
								c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
								
								E_param = path_set[j][1]
								E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
								
								dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
								angDist = bayes.mahabAngDist(ang1, ang2, E_param[2,2])
		
								#if dist <= DIST_THRESHOLD and abs(i-j) < 5:
								#if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD and abs(i-j) < 5:
								if dist <= DIST_THRESHOLD and angDist <= ANG_THRESHOLD:
									#print "distance:", i, j, dist
									#print "angDist:", i, j, angDist						
		
									pair_candidates.append([dist,i,j,loc2, angDist, point1, point2, ang1, ang2, m, n])
		
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

		pair_unique.sort()
		
		#print "UNIQUE PAIRS"
		for k in range(len(pair_unique)):
			p = pair_unique[k]
			dist = p[0]
			i = p[1]
			j = p[2]
			transform = p[3]
			angDist = p[4]
			point1 = p[5]
			point2 = p[6]

			#print "distance:", i, j, dist
			#print "angDist:", i, j, angDist
			
			#self.renderCornerConstraint(i, j, transform, point1, point2, k)
	
		return pair_unique



	def findSweepCandidates(self, paths):

		SENSOR_RADIUS = 0.0
		DIST_THRESHOLD = 1.0
		
		" select pairs to attempt a sensor constraint "
		pair_candidates = []
		
		for i in range(self.numNodes):
			path_set = paths[i]
			
			ind = range(i+2,self.numNodes)

			for j in ind:	
				loc2 = path_set[j][0]
				
				" c1 and c2 are centroids of sweep regions "
				node1 = self.nodeHash[i]
				sweepOffset1 = node1.getSweepCentroid()
				
				node2 = self.nodeHash[j]
				sweepOffset2 = node2.getSweepCentroid()				
				node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
				sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
				
				c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
				c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
				
				dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
				
				print "distance:", i, j, dist

				if dist <= DIST_THRESHOLD and abs(i-j) < 5:
					pair_candidates.append([dist,i,j,loc2])
				
			ind = range(0,i-1)
			ind.reverse()
			
			for j in ind:
				
				loc2 = path_set[j][0]
				
				" c1 and c2 are centroids of sweep regions "
				node1 = self.nodeHash[i]
				sweepOffset1 = node1.getSweepCentroid()
				
				node2 = self.nodeHash[j]
				sweepOffset2 = node2.getSweepCentroid()				
				node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
				sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
				
				c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
				c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
				
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
				
				dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
				
				print "distance:", i, j, dist
				if dist <= DIST_THRESHOLD and abs(i-j) < 5:
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

		pair_unique.sort()
		
		print "UNIQUE PAIRS"
		for p in pair_unique:
			i = p[1]
			j = p[2]
			
			
			loc2 = paths[i][j][0]

			" c1 and c2 are centroids of sweep regions "
			node1 = self.nodeHash[i]
			sweepOffset1 = node1.getSweepCentroid()
			
			node2 = self.nodeHash[j]
			sweepOffset2 = node2.getSweepCentroid()				
			node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
			sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
			
			c1 = matrix([[sweepOffset1[0]],[sweepOffset1[1]]],dtype=float)
			c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
			
			
			E_param = paths[i][j][1]
			E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
			
			dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
			
			print "distance:", i, j, dist
	
		return pair_unique

	def findAlphaHullCandidates(self, paths, targetNode = -1):

		SENSOR_RADIUS = 0.5
		#DIST_THRESHOLD = 1.0
		DIST_THRESHOLD = 3.0
		
		" select pairs to attempt a sensor constraint "
		pair_candidates = []
		

		if targetNode >= 0:

			path_set = paths[targetNode]
			
			ind = range(targetNode+2,self.numNodes)

			for j in ind:	
				loc2 = path_set[j][0]	
				
				c1 = matrix([[0.0],[0.0]],dtype=float)
				c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
				
				dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
				
				print targetNode, j, dist
				if dist <= DIST_THRESHOLD:
					pair_candidates.append([dist,targetNode,j,loc2])
				
			ind = range(0,targetNode-1)
			ind.reverse()
			
			for j in ind:
				
				loc2 = path_set[j][0]
				
				c1 = matrix([[0.0],[0.0]],dtype=float)
				c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
				E_param = path_set[j][1]
				E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
				
				dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
				
				print targetNode, j, dist
				if dist <= DIST_THRESHOLD:
					pair_candidates.append([dist,targetNode,j,loc2])

		else:
			
			
			for i in range(self.numNodes):
				path_set = paths[i]
				
				ind = range(i+2,self.numNodes)
	
				for j in ind:	
					loc2 = path_set[j][0]
					c1 = matrix([[0.],[0.]],dtype=float)
					c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
					E_param = path_set[j][1]
					E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
					
					dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
					
					#print "distance:", i, j, dist
	
					if dist <= DIST_THRESHOLD:
						pair_candidates.append([dist,i,j,loc2])
					
				ind = range(0,i-1)
				ind.reverse()
				
				for j in ind:
									
					loc2 = path_set[j][0]				
					c1 = matrix([[0.],[0.]],dtype=float)
					c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
					E_param = path_set[j][1]
					E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
					
					dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
					
					#print "distance:", i, j, dist
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

		pair_unique.sort()
		
		#print "UNIQUE PAIRS"
		for p in pair_unique:
			i = p[1]
			j = p[2]
				
			loc2 = paths[i][j][0]
			c1 = matrix([[0.],[0.]],dtype=float)
			c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
			
			E_param = paths[i][j][1]
			E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
			
			dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
			
			#print "distance:", i, j, dist
	
		return pair_unique

	def correctPoses3(self):

		""" 
		things to convert for map correction to the GPAC coordinate system 
		1. estPose to GPAC estPose
		2. convert hull points to GPAC coordinates before adding covariances
		3. convert edges from root end points to GPAC centroid endpoints (motion constraints)
		
		
		"""
		" FIXME:  check case where there are less than one hypotheses and don't perform hypothesis weeding "  

		" max difference between node numbers "
		NODE_SEP = 5

		if self.numNodes < 2:
			return

		print "check_A"
		if self.currNode.isDirty():
			self.currNode.synch()
		
		nodes = []
		for k, v in self.nodeHash.items():
			nodes.append(k)
		nodes.sort()
		
		edges = []
		for k, v in self.edgeHash.items():
			edges.append(k)
	
		
		edge_attrs = []
		for i in range(len(edges)):
			edge_attrs.append( self.edgeHash[edges[i]] )
	
		" save the graph to file minus the local map content "
		fp = open("graph_out.txt", 'w')
		fp.write(repr((nodes, edges, edge_attrs)))
		fp.close()

		" perform dijkstra projection over all nodes"
		paths1 = []
		for i in range(self.numNodes):
			paths1.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		edgeHash1 = self.edgeHash
		
		" create new consistent overlap hypotheses"
		#added_hypotheses = self.addOverlaps(paths1)

		added_hypotheses = []

		" merge and add them to the graph "
		v_list, new_constraints = self.doToro(self.merged_constraints + added_hypotheses, fileName = "probe")
		
		self.edgeHash = {}
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]

		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])				

		" perform dijkstra projection over all nodes"
		paths2 = []
		for i in range(self.numNodes):
			paths2.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
		
		edgeHash2 = self.edgeHash
		added_hypotheses = self.addShapeConstraints(paths2)

		v_list, new_constraints = self.doToro(self.merged_constraints + added_hypotheses, fileName = "probe")

		self.edgeHash = {}		
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])		
		
		"""
		" print out the edges and the dijkstra projection paths "
		for i in range(self.numNodes):
			for j in range(self.numNodes):
				
				if i != j:
					offset1 = paths1[i][j][0]
					offset2 = paths2[i][j][0]
	
					loc1 = paths1[i][j][0]
					loc2 = paths2[i][j][0]
					
					" c1 and c2 are centroids of sweep regions "
					#originNode = self.nodeHash[i]
					#sweepOffset = originNode.getSweepCentroid()
					
					#node1 = self.nodeHash[j]
					#sweepOffset1 = node1.getSweepCentroid()				
					#node1Profile = Pose([loc1[0,0], loc1[1,0], loc1[2,0]])
					#sweep1Centroid = node1Profile.convertLocalToGlobal(sweepOffset1)
	
					#node2 = self.nodeHash[j]
					#sweepOffset2 = node2.getSweepCentroid()				
					#node2Profile = Pose([loc2[0,0], loc2[1,0], loc2[2,0]])
					#sweep2Centroid = node2Profile.convertLocalToGlobal(sweepOffset2)
					
					sweepOffset = [0,0]
					sweep1Centroid = [loc1[0,0], loc1[1,0]]
					sweep2Centroid = [loc2[0,0], loc2[1,0]]
					
					cO = matrix([[sweepOffset[0]],[sweepOffset[1]]],dtype=float)
					c1 = matrix([[sweep1Centroid[0]],[sweep1Centroid[1]]],dtype=float)
					c2 = matrix([[sweep2Centroid[0]],[sweep2Centroid[1]]],dtype=float)
	
					E_param1 = paths1[i][j][1]
					E_param2 = paths2[i][j][1]
				
					E1 = matrix([[E_param1[0,0], E_param1[0,1]],[E_param1[1,0], E_param1[1,1]]],dtype=float)
					E2 = matrix([[E_param2[0,0], E_param2[0,1]],[E_param2[1,0], E_param2[1,1]]],dtype=float)
					
					dist1 = bayes.mahab_dist(cO, c1, 0, 0, E1)
					dist2 = bayes.mahab_dist(cO, c2, 0, 0, E2)
					
					#print "distance:", i, j, dist
	
					print i, j, ":", dist1, offset1[0,0], offset1[1,0], offset1[2,0], E1[0,0], E1[1,1]
					print i, j, ":", dist2, offset2[0,0], offset2[1,0], offset2[2,0], E2[0,0], E2[1,1]
					print 
					
		"""

		"""
		print "EDGES1:"
		for k, v in edgeHash1.items():
			offset = [v[0][0,0], v[0][1,0], v[0][2,0]]
			covE = [v[1][0,0], v[1][1,1], v[1][2,2]]
			print k, offset, covE
		print "EDGES2:"
		for k, v in edgeHash2.items():
			offset = [v[0][0,0], v[0][1,0], v[0][2,0]]
			covE = [v[1][0,0], v[1][1,1], v[1][2,2]]
			print k, offset, covE
		"""

	def relaxCorrect(self):

		" take the current graph, relax the constraints and start performing new hypotheses to try and create more consistent map "

		if self.numNodes < 2:
			return

		if self.currNode.isDirty():
			self.currNode.synch()
		
		" edges between i and i+1 where i is even should remain unrelaxed "

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.E_motion
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.E_overlap			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])


		" merge the constraints with TORO"
		#self.merged_constraints = self.toroMerge(totalConstraints)
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")

		for const in self.merged_constraints:
			covE = const[3]
			covE[0,0] = 0.5
			covE[1,1] = 0.5
			covE[2,2] = 0.25
			covE[0,1] = covE[1,0] = covE[1,2] = covE[2,1] = 0.0
			const[3] = covE

		newEdgeHash = {}		
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			newEdgeHash[(node1,node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]

			
		" perform dijkstra projection over all nodes"
		paths1 = []
		for i in range(self.numNodes):
			paths1.append(bayes.dijkstra_proj(i, self.numNodes, newEdgeHash))
			print "paths", i, paths1[-1]
			
		

		edgeHash1 = newEdgeHash
		
		" create new consistent overlap hypotheses"
		added_hypotheses = self.addOverlaps(paths1)

		" merge and add them to the graph "
		#new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])
		v_list, new_constraints = self.doToro(self.merged_constraints + added_hypotheses, fileName = "probe")
		
		self.edgeHash = {}
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]

		" perform dijkstra projection over all nodes"
		paths2 = []
		for i in range(self.numNodes):
			paths2.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
		
		edgeHash2 = self.edgeHash
		added_hypotheses = self.addShapeConstraints(paths2)

		v_list, new_constraints = self.doToro(self.merged_constraints + added_hypotheses, fileName = "probe")

		self.edgeHash = {}		
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])				
		
	def localizeCurrentNode(self):

		" take the current graph, relax the constraints and start performing new hypotheses to try and create more consistent map "

		if self.numNodes < 2:	
			return

		if self.currNode.isDirty():
			self.currNode.synch()

		" edges between i and i+1 where i is even should remain unrelaxed "

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = self.E_motion
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = self.E_overlap			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.backtrack_constraints)):
			nodeID1 = self.backtrack_constraints[i][0]
			nodeID2 = self.backtrack_constraints[i][1]
			transform = self.backtrack_constraints[i][2]
			covE = self.backtrack_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])


		" merge the constraints with TORO"
		#self.merged_constraints = self.toroMerge(totalConstraints)
		v_list, self.merged_constraints = self.doToro(totalConstraints, fileName = "probe_init")

		newEdgeHash = {}		
		for i in range(len(self.merged_constraints)):
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			newEdgeHash[(node1,node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]
			
		" perform dijkstra projection from latest node"		

		paths1 = []
		for i in range(self.numNodes):
			paths1.append(bayes.dijkstra_proj(i, self.numNodes, newEdgeHash))
			
		" create new consistent overlap hypotheses"
		#added_hypotheses = self.addOverlaps(paths1, targetNode = self.numNodes-1)

		#new_overlap = self.findBestOverlap(paths1, self.numNodes-1)
		eTup = (self.numNodes-2, self.numNodes-1)
		new_overlap = [eTup[0], eTup[1], newEdgeHash[eTup][0], newEdgeHash[eTup][1] ]
		
		

		
		" make hypothesized constraints out of the pairs "
		if new_overlap != []:
			n1 = new_overlap[0]
			n2 = new_overlap[1]
			hull1 = computeHull(self.nodeHash[n1], sweep = False)
			hull2 = computeHull(self.nodeHash[n2], sweep = False)
			transform = new_overlap[2]
			offset, covar, cost = self.makeSensorConstraint(transform, n1, n2, hull1, hull2)
			transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
			
			sensor_constraint = [n1,n2,transform,covar]

			new_overlap = sensor_constraint
			
		#print "CONSTRAINTS BEFORE doToro()"

		#print "added_hypotheses:", added_hypotheses
		#print "self.merged_constraints:", self.merged_constraints

		" merge and add them to the graph "
		#new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])
		#v_list, new_constraints = self.doToro(self.merged_constraints + added_hypotheses, fileName = "probe")

		if new_overlap != []:
			v_list, new_constraints = self.doToro(self.merged_constraints + [new_overlap], fileName = "probe")
		else:
			v_list = v_list
			new_constraints = self.merged_constraints
		
		self.edgeHash = {}
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		" FIXME:  save this new overlap constraint "
		if new_overlap != []:
			#self.backtrack_constraints.append(new_overlap)
			self.sensor_constraints.append(new_overlap)

		for v in v_list:
			node = self.nodeHash[v[0]]
			node.setGPACPose(v[1])				
		
	def addOverlaps(self, paths, targetNode = -1):

		print "numNodes:", self.numNodes
		
		print len(paths), "paths"
		for path in paths:
			print len(path)

		" find candidates to perform constraints "
		if targetNode >= 0:
			overlap_pairs = self.findOverlapCandidates(paths, targetNode)
		else:
			overlap_pairs = self.findOverlapCandidates(paths)
			
		
		" make hypothesized constraints out of the pairs "
		hypotheses = []
		for i in range(len(overlap_pairs)):
			p = overlap_pairs[i]
			transform = p[3]
			n1 = p[1]
			n2 = p[2]
			transform, covar = self.makeOverlapHypothesis(transform, n1, n2 )
			hypotheses.append([p[1],p[2],transform,covar])


		" compute the consistency between each pair of hypotheses "
		results = []
		for i in range(len(hypotheses)):
			for j in range(i+1, len(hypotheses)):
				
				hyp1 = hypotheses[i]
				hyp2 = hypotheses[j]

				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "
				
				#covE = matrix([[ 0.1, 0.0, 0.0 ], [ 0.0, 0.1, 0.0], [ 0.0, 0.0, pi/4.0 ]],dtype=float) 
				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				
				#precision = [invMat[0,0], invMat[0,1], invMat[1,1], invMat[2,2], invMat[0,2], invMat[1,2]]

				results.append([err, i, j])
				
				" max difference between node numbers "
				
				"FIXME:  this should be hop count, not difference in node ID "
				#if abs(m2-n1) < NODE_SEP and abs(n2-m1) < NODE_SEP:
				#	#results.append([err, i, j, [result4[0,0], result4[1,0], result4[2,0]], precision])
				#	results.append([err, i, j])
				#else:
				#	results.append([None,i,j])
					
		results.sort()
		
		print "len(hypotheses)", len(hypotheses)
		if len(hypotheses) < 2:
			"too little hypotheses, rejecting"
			f = open("overlap_constraints.txt", 'w')
			f.write(repr([]))
			f.close()
		
			return []
		

		" set all proscribed hypothesis sets to maximum error "
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]
				
		maxError = maxError*2
		for result in results:
			if result[0] == None:
				result[0] = maxError

		" create the consistency matrix "
		A = matrix(len(hypotheses)*len(hypotheses)*[0.0], dtype=float)
		A.resize(len(hypotheses), len(hypotheses))

		" populate consistency matrix "							
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]
		

		" disable the truncation "
		set_printoptions(threshold=nan)		
		fn = open("A_overlap.txt",'w')
		fn.write(repr(A))
		fn.close()
		
		
		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			#w = scsgp.getIndicatorVector(e[0])
			w = scsgp.getIndicatorVector(e[1])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break
			else:
				print
		#print "INDICATOR:", w
		
		selected = []	
		for i in range(len(hypotheses)):
			if w[i,0] >= 1.0:
				selected.append(hypotheses[i])

		#print len(selected), "added hypotheses"
		#print len(self.merged_constraints), "merged constraints"

		f = open("overlap_constraints.txt", 'w')
		#f.write(repr(added_hypotheses))
		f.write(repr(selected))
		f.close()
		
		#return added_hypotheses
		return selected


	def addPathConstraints(self, pathNodes, targetNodeID):

		print "addPathConstraints:"
		print "pathNodes:", pathNodes
		print "targetNodeID:", targetNodeID

		paths = []
		for i in range(self.numNodes):
			paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

		cart_distances = []
		targetNode = self.nodeHash[targetNodeID]
		targetPose = targetNode.getGlobalGPACPose()
		for i in range(len(self.nodeSets[0])):
			
			candNode = self.nodeHash[self.nodeSets[0][i]]
			candPose = candNode.getGlobalGPACPose()
			
			dist = sqrt((candPose[0]-targetPose[0])**2 + (candPose[1]-targetPose[1])**2)
			cart_distances.append(dist)

		print "non-candidate node distances:"
		for i in range(len(self.nodeSets[0])):
			print self.nodeSets[0][i], cart_distances[i]

		print "candidate node distances", targetNodeID, ":"

		" compute cartesian distances "
		cart_distances = []
		angle_diffs = []
		targetNode = self.nodeHash[targetNodeID]
		targetPose = targetNode.getGlobalGPACPose()
		for i in range(len(pathNodes)):
			
			candNode = self.nodeHash[pathNodes[i]]
			candPose = candNode.getGlobalGPACPose()
			
			dist = sqrt((candPose[0]-targetPose[0])**2 + (candPose[1]-targetPose[1])**2)
			angDelta = diffAngle(candPose[2],targetPose[2])
			print pathNodes[i], dist, angDelta, targetPose, candPose
			cart_distances.append(dist)

			angle_diffs.append(angDelta)
		#for i in range(len(pathNodes)):
		#	print pathNodes[i], cart_distances[i]

		constraintNodes = deepcopy(pathNodes)

		" remove the node that is the target from consideration "
		delIndex = pathNodes.index(targetNodeID)
		constraintNodes.pop(delIndex)
		cart_distances.pop(delIndex)
		angle_diffs.pop(delIndex)

		" remove the paired inplace node if it exists "
		if targetNodeID % 2 == 0:
			pairNodeID = targetNodeID + 1
		else:
			pairNodeID = targetNodeID - 1

		try:
			delIndex = pathNodes.index(pairNodeID)
			constraintNodes.pop(delIndex)
			cart_distances.pop(delIndex)
			angle_diffs.pop(delIndex)
		except:
			pass

		#cart_distances.remove(pairNodeID)



		#constraintNodes = pathNodes[:-2]
		#constraintNodes = pNodes
		constraintPairs = []
		
		PATH_MATCH_DIST = 1.0
		#PATH_MATCH_DIST = 0.5
		

		for i in range(len(constraintNodes)):
			nodeID = constraintNodes[i]
			dist = cart_distances[i]
			angDelta = angle_diffs[i]

			if dist < PATH_MATCH_DIST:
				state1 = self.nodeHash[targetNodeID].getDeparture()
				state2 = self.nodeHash[nodeID].getDeparture()
				if state1 == state2:
					print targetNodeID, "and", nodeID, "have same departure state", state1, state2
					
					" check if a constraint already exists between these two nodes "
					oldConstraints = self.getEdges(targetNodeID, nodeID)
					isExist = False
					for constraint in oldConstraints:
						isExist = True

					if not isExist:
						constraintPairs.append([targetNodeID,nodeID,angDelta])
					else:
						print "constraint already exists between", targetNodeID, "and", nodeID
						
				else:
					print targetNodeID, "and", nodeID, "have different departure state", state1, state2 
			else:	
				print targetNodeID, "and", nodeID, "are too far away", dist 

		print "received", len(constraintPairs), "possible pairs"

		" now, take only the 5 oldest nodes to constrain with "
		constraintPairs.sort()
		if len(constraintPairs) > 3:
			constraintPairs = constraintPairs[:3]
			print "reducing to only", len(constraintPairs)

		PAIR_ANG_DELTA = 0.3

		" make hypothesized constraints out of the pairs "
		constraintResults = []
		for i in range(len(constraintPairs)):
			p = constraintPairs[i]
			n1 = p[0]
			n2 = p[1]
			angDiff = p[2]

			try:

				if self.nodeHash[n1].getIsFeatureless():	
					
					if self.nodeHash[n2].getIsFeatureless():	
						transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)				
	
						if hist[1] > 0 or hist[2] > 0:
							" disregard if the fit is poor "
							pass
						else:
							
							" compare this new orientation to the old, is it too different?  If it is, throw it out."
							angDelta = diffAngle(transform[2,0], angDiff)
							
							if fabs(angDelta) < PAIR_ANG_DELTA:
								#constraintResults.append([n1, n2, transform, deepcopy(self.E_featureless), angDiff])
								constraintResults.append([n1, n2, transform, deepcopy(self.E_featureless), angDelta])
							else:
								print "rejecting", (n1,n2), "constraint because angDelta is", angDelta, "with transform angle", transform[2,0]

				else:
					if not self.nodeHash[n2].getIsFeatureless():

						transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False, uRange = 3.0)				
	
						if hist[1] > 0 or hist[2] > 0:
							" disregard if the fit is poor "
							pass
						else:
	
							" compare this new orientation to the old, is it too different?  If it is, throw it out."
							angDelta = diffAngle(transform[2,0], angDiff)
							
							if fabs(angDelta) < PAIR_ANG_DELTA:								
								constraintResults.append([n1, n2, transform, deepcopy(self.E_junction), angDelta])
							else:
								print "rejecting", (n1,n2), "constraint because angDelta is", angDelta, "with transform angle", transform[2,0]
							
						" TODO: "
						" 1) find the overlapping region curves "
						" 2) determine if the curves are straight and featureless "
						" 3) if so, add high variance along the tangent "
				
			except:
				pass
		
		print "received", len(constraintResults), "constraint results"
		for result in constraintResults:
			print result
		#if len(constraintResults) > 0:
		#	print constraintResults[0]	


		totalHypotheses = constraintResults
		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) == 0:
			print "no hypotheses!  returning "
			newConstraints = []
			#return
		
		elif len(totalHypotheses) == 1:
			newConstraints = totalHypotheses

		else: 
					
			" set all proscribed hypothesis sets to maximum error "	
			maxError = 0.0
			for result in results:
				if result[0] != None and result[0] > maxError:
					maxError = result[0]
	
			maxError = maxError*2
	
			for result in results:
				if result[0] == None:
					result[0] = maxError
		
			" create the consistency matrix "
			A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
			A.resize(len(totalHypotheses), len(totalHypotheses))
				
			" populate consistency matrix "					
			for result in results:
				i = result[1]
				j = result[2]
				A[i,j] = maxError - result[0]
				A[j,i] = maxError - result[0]
	
			" do graph clustering of consistency matrix "
			w = []
			e = []
			for i in range(100):
				e, lmbda = scsgp.dominantEigenvectors(A)
				w = scsgp.getIndicatorVector(e[0])
				w2 = scsgp.getIndicatorVector(e[1])
				if len(e) <= 1:
					break
	
				" threshold test l1 / l2"
				ratio = lmbda[0][0,0] / lmbda[1][0,0]
				if ratio >= 2.0:
					break
	
			selected = []	
			for i in range(len(totalHypotheses)):
				if w[i,0] >= 1.0:
					selected.append(totalHypotheses[i])
	
			selected2 = []	
			for i in range(len(totalHypotheses)):
				if w2[i,0] >= 1.0:
					selected2.append(totalHypotheses[i])
	
			newConstraints = selected2

			newConstraints = totalHypotheses
		print "adding", len(newConstraints), "medial overlap constraints on path"

		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], CORNER_PRIORITY)

		" merge the constraints "
		#self.mergePriorityConstraints()
		
		#return selected, selected2

		return

		return selected2

		
		#self.makeMedialOverlapConstraint(targetNodeID, nodeID, isMove = False)
		
		
		return
		
		"""		
		self.deleteAllPriority(SHAPE_PRIORITY)

		newConstraints = self.addOverlapConstraints(paths)

		
		print "adding", len(newConstraints), "overlap constraints"

		for const in newConstraints:
			n1 = const[0]
			n2 = const[1]
			transform = const[2]
			covE = const[3]

			self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)

		" merge the constraints "
		self.mergePriorityConstraints()
		"""

		#sensor_pairs = self.findOverlapCandidates(paths)
		
		#return	
		

		



	def addOverlapConstraints(self, paths, targetNode = None, matchNode = None):

		if targetNode != None:
			if matchNode != None:
				sensor_pairs = [[0.0,targetNode,matchNode,paths[targetNode][matchNode][0]]]
			else:
				sensor_pairs = self.findOverlapCandidates(paths, targetNode = targetNode)
		else:
			sensor_pairs = self.findOverlapCandidates(paths)
		
		print "received", len(sensor_pairs), "possible pairs"
		#return	
		
		" make hypothesized constraints out of the pairs "
		constraintResults = []
		for i in range(len(sensor_pairs)):
			p = sensor_pairs[i]
			transform = p[3]
			n1 = p[1]
			n2 = p[2]

			try:
				transform, covE, hist = self.makeMedialOverlapConstraint(n1, n2, isMove = False)				
				constraintResults.append([n1, n2, transform, covE])
			except:
				pass

		print "received", len(constraintResults), "constraint results"
		if len(constraintResults) > 0:
			print constraintResults[0]	
		
		totalHypotheses = constraintResults
		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) < 2:
			print "too little hypotheses, rejecting"
			return []
		
		" set all proscribed hypothesis sets to maximum error "	
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]

		maxError = maxError*2

		for result in results:
			if result[0] == None:
				result[0] = maxError
	
		" create the consistency matrix "
		A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
		A.resize(len(totalHypotheses), len(totalHypotheses))
			
		" populate consistency matrix "					
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]

		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			w = scsgp.getIndicatorVector(e[0])
			w2 = scsgp.getIndicatorVector(e[1])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break

		selected = []	
		for i in range(len(totalHypotheses)):
			if w[i,0] >= 1.0:
				selected.append(totalHypotheses[i])

		selected2 = []	
		for i in range(len(totalHypotheses)):
			if w2[i,0] >= 1.0:
				selected2.append(totalHypotheses[i])

		#return selected, selected2

		return selected2


	def addShapeConstraints(self, paths, targetNode = None, matchNode = None):

		if targetNode != None:
			if matchNode != None:
				sensor_pairs = [[0.0,targetNode,matchNode,paths[targetNode][matchNode][0]]]
			else:
				sensor_pairs = self.findAlphaHullCandidates(paths, targetNode = targetNode)
		else:
			sensor_pairs = self.findAlphaHullCandidates(paths)
		
		" compute the alpha hulls of the selected candidates "
		#hull_computed = [False for i in range(self.numNodes)]
		#a_hulls = [0 for i in range(self.numNodes)]
		
		a_hulls = self.a_hulls
		a_medial = self.a_medial
		
		#for i in range(self.numNodes):
		#	a_hulls[i] = computeHull(self.nodeHash[i], sweep = False)
		#	hull_computed[i] = True

		constraintData = []
		constraintData2 = []
		
		" make hypothesized constraints out of the pairs "

		constraintResults = []
		hypotheses = []		
		for i in range(len(sensor_pairs)):
			p = sensor_pairs[i]
			transform = p[3]
			n1 = p[1]
			n2 = p[2]
			print n1, n2
			print len(a_hulls), len(a_medial)
			#constraintData.append(self.makeShapeConstraintData(transform, n1, n2, a_hulls[n1], a_hulls[n2]))
			constraintData2.append(self.makeShapeConstraintData2(transform, n1, n2, a_hulls[n1], a_hulls[n2], a_medial[n1], a_medial[n2]))
			#return (estPose1, estPose2, hull1, hull2, posture1_unstable, posture1_stable, posture2_unstable, posture1_stable, costThresh, minMatchDist, plotIteration, n1, n2)
			
			if len(constraintData2) >= 16:
				constraintResults += self.computeAllShapeConstraints(constraintData2)
				constraintData2 = []

		if len(constraintData2) > 0:
			constraintResults += self.computeAllShapeConstraints(constraintData2)
			constraintData2 = []
			
		#constraintResults = self.computeAllShapeConstraints(constraintData)

		#print constraintData2[0]	
		print "received", len(constraintResults), "constraint results"
		if len(constraintResults) > 0:
			print constraintResults[0]	

		self.goodHypotheses = []
		self.badHypotheses1 = []
		self.badHypotheses2 = []
		
		for i in range(len(constraintResults)):
			p = sensor_pairs[i]
			result = constraintResults[i]
			offset = result[1][0]
			cost = result[1][1]
			status = result[1][2]
			covar = self.E_sensor		
			transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)

			if cost < 1.5:
				hypotheses.append([p[1],p[2],transform,covar])
				#self.goodHypotheses.append([p[1],p[2],transform,covar])
				print "accepted", i, ":", "shape constraint between", p[1], "and", p[2], "cost =", cost
			else:
				if status == 0:
					badStatus = 2
				else:
					badStatus = status
				self.badHypotheses1.append([p[1],p[2],transform,covar, badStatus])
				self.badHypotheses2.append([p[1],p[2],transform,covar, badStatus])
				print "rejected", i, ":", "shape constraint between", p[1], "and", p[2], "cost =", cost


		#hypotheses = []		
		#for i in range(len(sensor_pairs)):
		#	p = sensor_pairs[i]
		#	transform = p[3]
		#	n1 = p[1]
		#	n2 = p[2]
		#	constraintData.append(self.makeShapeConstraintData(transform, n1, n2, a_hulls[n1], a_hulls[n2]))
		#	offset, covar, cost = self.makeSensorConstraint(transform, n1, n2, a_hulls[n1], a_hulls[n2])
		#	transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
		
		#	if cost < 1.5:
		#		hypotheses.append([p[1],p[2],transform,covar])
		#		print "hypothesis", i, ":",  "sensor constraint between", n1, "and", n2

		#set_printoptions(threshold=nan)
		#f = open("shapeInitData.txt", 'w')
		#f.write(repr(constraintData))
		#f.close()

		" compute the consistency between each pair of hypotheses "

		self.sensorHypotheses += hypotheses
		
		#totalHypotheses = self.sensor_constraints + self.sensorHypotheses
		#totalHypotheses = self.sensorHypotheses
		
		totalHypotheses = hypotheses
		
		results = []
		for i in range(len(totalHypotheses)):
			for j in range(i+1, len(totalHypotheses)):
				
				hyp1 = totalHypotheses[i]
				hyp2 = totalHypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				results.append([err, i, j])

					
		results.sort()
		
		print "len(totalHypotheses) =", len(totalHypotheses)
				
		if len(totalHypotheses) < 2:
			print "too little hypotheses, rejecting"
			return [], []
		
		" set all proscribed hypothesis sets to maximum error "	
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]

		maxError = maxError*2

		for result in results:
			if result[0] == None:
				result[0] = maxError
	
		" create the consistency matrix "
		A = matrix(len(totalHypotheses)*len(totalHypotheses)*[0.0], dtype=float)
		A.resize(len(totalHypotheses), len(totalHypotheses))
			
		" populate consistency matrix "					
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]

		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			w = scsgp.getIndicatorVector(e[0])
			w2 = scsgp.getIndicatorVector(e[1])
			#w3 = scsgp.getIndicatorVector(e[2])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break

		selected = []	
		for i in range(len(totalHypotheses)):
			if w[i,0] >= 1.0:
				selected.append(totalHypotheses[i])
			else:
				self.badHypotheses1.append(totalHypotheses[i] + [3])

		selected2 = []	
		for i in range(len(totalHypotheses)):
			if w2[i,0] >= 1.0:
				selected2.append(totalHypotheses[i])
			else:
				self.badHypotheses2.append(totalHypotheses[i]+ [3])


		"""


		" find candidates to perform sensor constraints "
		#sensor_pairs = self.findSweepCandidates(paths)
		sensor_pairs = self.findAlphaHullCandidates(paths)
		
		" compute the alpha hulls of the selected candidates "
		hull_computed = [False for i in range(self.numNodes)]
		a_hulls = [0 for i in range(self.numNodes)]
		
		for i in range(self.numNodes):
			a_hulls[i] = computeHull(self.nodeHash[i], sweep = False)
			#a_hulls[i] = computeHull(self.nodeHash[i], sweep = True)
			hull_computed[i] = True
		
		" make hypothesized constraints out of the pairs "
		hypotheses = []		
		for i in range(len(sensor_pairs)):
			p = sensor_pairs[i]
			transform = p[3]
			n1 = p[1]
			n2 = p[2]
			offset, covar, cost = self.makeSensorConstraint(transform, n1, n2, a_hulls[n1], a_hulls[n2])
			transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
			
			if cost < 1.5:
				hypotheses.append([p[1],p[2],transform,covar])
				print "hypothesis", i, ":",  "sensor constraint between", n1, "and", n2

		" compute the consistency between each pair of hypotheses "
		results = []
				
		for i in range(len(hypotheses)):
			for j in range(i+1, len(hypotheses)):
				
				hyp1 = hypotheses[i]
				hyp2 = hypotheses[j]
				
				m1 = hyp1[0]
				m2 = hyp1[1]
				
				n1 = hyp2[0]
				n2 = hyp2[1]
				
				Tp1 = paths[m2][n1][0]
				Ep1 = paths[m2][n1][1]
				
				Tp2 = paths[n2][m1][0]
				Ep2 = paths[n2][m1][1]
				
				Th1 = hyp1[2]
				Th2 = hyp2[2]
				
				Ch1 = hyp1[3]
				Ch2 = hyp2[3]
				
				" m1->m2, m2->n1, n1->n2, n2->m1 "
				" Th1, Tp1, Th2, Tp2 "

				covE = Ch1
				result1 = Th1
				result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
				result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
				result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)
				
				invMat = scipy.linalg.inv(cov4)
				err = sqrt(numpy.transpose(result4) * invMat * result4)
				
				#precision = [invMat[0,0], invMat[0,1], invMat[1,1], invMat[2,2], invMat[0,2], invMat[1,2]]

				results.append([err, i, j])
				
				" max difference between node numbers "
				
				"FIXME:  this should be hop count, not difference in node ID "
				#if abs(m2-n1) < NODE_SEP and abs(n2-m1) < NODE_SEP:
				#	#results.append([err, i, j, [result4[0,0], result4[1,0], result4[2,0]], precision])
				#	results.append([err, i, j])
				#else:
				#	results.append([None,i,j])
					
		results.sort()
		
		print "len(hypotheses)", len(hypotheses)
		
		
		if len(hypotheses) < 2:
			"too little hypotheses, rejecting"
			f = open("sensor_constraints.txt", 'w')
			f.write(repr([]))
			f.close()
			return []
		
		" set all proscribed hypothesis sets to maximum error "	
		maxError = 0.0
		for result in results:
			if result[0] != None and result[0] > maxError:
				maxError = result[0]

		maxError = maxError*2

		for result in results:
			if result[0] == None:
				result[0] = maxError
				
		if len(hypotheses) == 0:
			set_printoptions(threshold=nan)		
			fn = open("A_sensor.txt",'w')
			fn.write(repr(matrix([[]])))
			fn.close()
			
			return []

		" create the consistency matrix "
		A = matrix(len(hypotheses)*len(hypotheses)*[0.0], dtype=float)
		A.resize(len(hypotheses), len(hypotheses))
			
		" populate consistency matrix "					
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = maxError - result[0]
			A[j,i] = maxError - result[0]

		" disable the truncation "
		set_printoptions(threshold=nan)		
		fn = open("A_sensor.txt",'w')
		fn.write(repr(A))
		fn.close()
		
		
		" do graph clustering of consistency matrix "
		w = []
		e = []
		for i in range(100):
			e, lmbda = scsgp.dominantEigenvectors(A)
			#w = scsgp.getIndicatorVector(e[0])
			w = scsgp.getIndicatorVector(e[1])
			if len(e) <= 1:
				break

			" threshold test l1 / l2"
			ratio = lmbda[0][0,0] / lmbda[1][0,0]
			if ratio >= 2.0:
				break
		#print "INDICATOR:", w
		selected = []	
		for i in range(len(hypotheses)):
			if w[i,0] >= 1.0:
				selected.append(hypotheses[i])
		"""

		"""
		added_hypotheses = []
		added = {}
		for i in range(len(results)):
			result = results[i]
			h1 = result[1]
			h2 = result[2]
			val = result[0]
			hyp1 = hypotheses[h1]
			hyp2 = hypotheses[h2]

			print "(%d -> %d)" % (hyp1[0], hyp1[1]), "(%d -> %d)" % (hyp2[0], hyp2[1]) , val, w[h1,0], w[h2,0], e[0][h1,0], e[0][h2,0], A[h1, h2]
	
	
			try:
				added[(hyp1[0],hyp1[1])]
			except:
				added[(hyp1[0],hyp1[1])] = True
				added_hypotheses.append(hyp1)
				
			try:
				added[(hyp2[0],hyp2[1])]
			except:
				added[(hyp2[0],hyp2[1])] = True
				added_hypotheses.append(hyp2)
			
		print len(added_hypotheses), "added hypotheses"
		"""
		
		#print len(selected), "added hypotheses"
		#print len(self.merged_constraints), "merged constraints"

		#f = open("sensor_constraints.txt", 'w')
		#f.write(repr(added_hypotheses))
		#f.write(repr(selected))
		#f.close()
		
		#return added_hypotheses
		return selected, selected2

	def computeAllShapeConstraints(self, constraintData):

		#constraintResults = self.computeAllShapeConstraints(constraintData)


		def send_problem(address, data):
			"""Download a piece of poetry from the given address."""
		
			sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
			sock.connect(address)
		
			sock.sendall(data + '\n')
		
			results = ''
		
			while True:
		
				# This is the 'blocking' call in this synchronous program.
				# The recv() method will block for an indeterminate period
				# of time waiting for bytes to be received from the server.
		
				bytes = sock.recv(1024)
		
				if not bytes:
					sock.close()
					break
		
				results += bytes
		
			return results
		
		results = eval(send_problem(('127.0.0.1', 10000), repr(constraintData)))
		
		return results
	
	def doToro(self, constraints, fileName = "probe"):
		
		
		
		" vertices "
		v_list = []

		for k, v in self.nodeHash.items():
			node1 = v
			estPose1 = node1.getGlobalGPACPose()
			v_list.append([k, [estPose1[0], estPose1[1], estPose1[2]]])

		if len(constraints) == 0:
			return v_list, []
		
		e_list = []
		for const in constraints:
				transform, covE = const[2], const[3]
				node1 = const[0]
				node2 = const[1]
				offset = [transform[0,0],transform[1,0],transform[2,0]]

				invMat = scipy.linalg.inv(covE)				
				prec = [invMat[0,0], invMat[0,1], invMat[1,1], invMat[2,2], invMat[0,2], invMat[1,2]]
				
				# inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr 
				e_list.append([node1,node2,offset,prec])			
		
		#print "Running TORO with", len(e_list), "constraints"
		" add hypotheses to the TORO file "
		
		#print "e_list:", e_list
		#print "v_list:", v_list

		toro.writeConstrainedToroGraph(".", fileName + ".graph", v_list, e_list)
		toro.executeToro("./" + fileName + ".graph")
	
		finalFileName = fileName + "-treeopt-final.graph"
		v_list2, e_list2 = toro.readToroGraph("./" + finalFileName)		
	
		final_constraints = []
		
		#print "TORO RESULTS:"
		for edge in e_list2:
			node1 = edge[0]
			node2 = edge[1]
			offset = edge[2]
			prec = edge[3]
			#print node1, node2, offset, prec

			invMat = matrix([[0.,0.,0.],
							[0.,0.,0.],
							[0.,0.,0.]])
			invMat[0,0] = prec[0]
			invMat[0,1] = invMat[1,0] = prec[1]
			invMat[1,1] = prec[2]
			invMat[2,2] = prec[3]
			invMat[0,2] = invMat[2,0] = prec[4]
			invMat[1,2] = invMat[2,1] = prec[5]
			
			transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
			covE = scipy.linalg.inv(invMat)				
			
			final_constraints.append([node1, node2, transform, covE])			

		return v_list2, final_constraints

	def saveConstraints(self):

		f = open("overlap_constraints_%04u.txt" % self.numNodes, 'w')
		f.write(repr(self.overlap_constraints))
		f.write("\n")
		f.close()		
			
		f = open("motion_constraints_%04u.txt" % (self.numNodes-1), 'w')
		f.write(repr(self.motion_constraints))
		f.write("\n")
		f.close()		

		f = open("sensor_constraints_%04u.txt" % (self.numNodes-1), 'w')
		f.write(repr(self.sensor_constraints))
		f.write("\n")
		f.close()		

		f = open("inplace_constraints_%04u.txt" % (self.numNodes-1), 'w')
		f.write(repr(self.inplace_constraints))
		f.write("\n")
		f.close()		

		f = open("backtrack_constraints_%04u.txt" % (self.numNodes-1), 'w')
		f.write(repr(self.backtrack_constraints))
		f.write("\n")
		f.close()		

	def renderNode(self, i, isSweep = False):
		
		node1 = self.nodeHash[i]
		
		" create the ground constraints "
		gndGPAC1Pose = node1.getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)

		posture1 = node1.getStableGPACPosture()

		hull1 = computeHull(node1, sweep = isSweep)

		segWidth = self.probe.robotParam['segWidth']
		segLength = self.probe.robotParam['segLength']
				
		for p in posture1:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "black")

		xP = []
		yP = []
		for p in hull1:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull1[0][0])
		yP.append(hull1[0][1])
		pylab.plot(xP,yP, color = "black")

	def renderGoodAndBadConstraint(self, id1, id2, transform, covar, cost, point1, point2, angDiff, oriDiff, fileName, renderCount):
		
		pylab.clf()
			
		stdev = 1

		covE = covar
		
		node1 = self.nodeHash[id1]
		node2 = self.nodeHash[id2]
		
		offset = [transform[0,0], transform[1,0], transform[2,0]]

		cholE = numpy.linalg.cholesky(covE[0:2,0:2])
		M = matrix([[transform[0,0]],[transform[1,0]]])
		circle = [matrix([[cos(phi)], [sin(phi)]]) for phi in arange(0.0, 2*pi, 0.1)]
		ellipse = []
		for p in circle:
			X = p
			Y = M + stdev * cholE * X
			ellipse.append([Y[0,0], Y[1,0]])

		" create the ground constraints "
		gndGPAC1Pose = node1.getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)
		gndGPAC2Pose = node2.getGndGlobalGPACPose()
		gnd_offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)		

		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()
		
		#getBoundaryPoints

		hull1 = computeHull(node1, sweep = True)
		hull2 = computeHull(node2, sweep = True)
		hull2_trans = []
		hull2_gnd = []

		posture_trans = []
		pose2Profile = Pose(offset)
		for p in posture2:
			posture_trans.append(pose2Profile.convertLocalOffsetToGlobal(p))

		for p in hull2:
			hull2_trans.append(pose2Profile.convertLocalToGlobal(p))
			
			
		point2_trans = pose2Profile.convertLocalToGlobal(point2)

		posture_gnd = []
		gndProfile = Pose(gnd_offset)
		for p in posture2:
			posture_gnd.append(gndProfile.convertLocalOffsetToGlobal(p))

		for p in hull2:
			hull2_gnd.append(gndProfile.convertLocalToGlobal(p))

		segWidth = self.probe.robotParam['segWidth']
		segLength = self.probe.robotParam['segLength']
				
		for p in posture_gnd:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "black")

		xP = []
		yP = []
		for p in hull2_gnd:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull2_gnd[0][0])
		yP.append(hull2_gnd[0][1])
		pylab.plot(xP,yP, color = "black")

		for p in posture1:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "green")

		xP = []
		yP = []
		for p in hull1:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull1[0][0])
		yP.append(hull1[0][1])
		pylab.plot(xP,yP, color = "green")

		for p in posture_trans:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "red")

		xP = []
		yP = []
		for p in hull2_trans:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull2_trans[0][0])
		yP.append(hull2_trans[0][1])
		pylab.plot(xP,yP, color = "red")


		xP = []
		yP = []
		for p in ellipse:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, linewidth=3, color='blue')
	
		xP = [point1[0], point2_trans[0]]
		yP = [point1[1], point2_trans[1]]
		pylab.scatter(xP, yP, color='k')
		
		pylab.title("%d ---> %d, cost = %f, angDiff = %f, oriDiff = %f" % (id1, id2, cost, angDiff, oriDiff))
		pylab.xlim(-4,4)
		pylab.ylim(-4,4)
		pylab.savefig(fileName + "_%04u.png" % renderCount)


	def renderCornerConstraint(self, id1, id2, transform, point1, point2, count):

		
		pylab.clf()
			
		stdev = 1

		covE = self.E_corner
		
		node1 = self.nodeHash[id1]
		node2 = self.nodeHash[id2]
		
		offset = [transform[0,0], transform[1,0], transform[2,0]]

		cholE = numpy.linalg.cholesky(covE[0:2,0:2])
		M = matrix([[transform[0,0]],[transform[1,0]]])
		circle = [matrix([[cos(phi)], [sin(phi)]]) for phi in arange(0.0, 2*pi, 0.1)]
		ellipse = []
		for p in circle:
			X = p
			Y = M + stdev * cholE * X
			ellipse.append([Y[0,0], Y[1,0]])

		" create the ground constraints "
		gndGPAC1Pose = node1.getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)
		gndGPAC2Pose = node2.getGndGlobalGPACPose()
		gnd_offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)		

		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()
		
		#getBoundaryPoints

		hull1 = computeHull(node1, sweep = True)
		hull2 = computeHull(node2, sweep = True)
		hull2_trans = []
		hull2_gnd = []

		posture_trans = []
		pose2Profile = Pose(offset)
		for p in posture2:
			posture_trans.append(pose2Profile.convertLocalOffsetToGlobal(p))

		for p in hull2:
			hull2_trans.append(pose2Profile.convertLocalToGlobal(p))
			
			
		point2_trans = pose2Profile.convertLocalToGlobal(point2)

		posture_gnd = []
		gndProfile = Pose(gnd_offset)
		for p in posture2:
			posture_gnd.append(gndProfile.convertLocalOffsetToGlobal(p))

		for p in hull2:
			hull2_gnd.append(gndProfile.convertLocalToGlobal(p))

		segWidth = self.probe.robotParam['segWidth']
		segLength = self.probe.robotParam['segLength']
				
		for p in posture_gnd:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "black")

		xP = []
		yP = []
		for p in hull2_gnd:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull2_gnd[0][0])
		yP.append(hull2_gnd[0][1])
		pylab.plot(xP,yP, color = "black")

		for p in posture1:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "green")

		xP = []
		yP = []
		for p in hull1:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull1[0][0])
		yP.append(hull1[0][1])
		pylab.plot(xP,yP, color = "green")

		for p in posture_trans:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "red")

		xP = []
		yP = []
		for p in hull2_trans:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull2_trans[0][0])
		yP.append(hull2_trans[0][1])
		pylab.plot(xP,yP, color = "red")


		xP = []
		yP = []
		for p in ellipse:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, linewidth=3, color='blue')
	
		xP = [point1[0], point2_trans[0]]
		yP = [point1[1], point2_trans[1]]
		pylab.scatter(xP, yP, color='k')
		
		pylab.title("Corner Constraint %d ---> %d" % (id1, id2))
		pylab.xlim(-4,4)
		pylab.ylim(-4,4)
		pylab.savefig("cornerConstraint_%04u.png" % count)

	
	def renderSingleConstraint(self, const, isSweep = False):

		stdev = 1

		id1 = const[0]
		id2 = const[1]			
		transform = const[2]
		covE = const[3]
		
		node1 = self.nodeHash[id1]
		node2 = self.nodeHash[id2]
		
		offset = [transform[0,0], transform[1,0], transform[2,0]]

		cholE = numpy.linalg.cholesky(covE[0:2,0:2])
		M = matrix([[transform[0,0]],[transform[1,0]]])
		circle = [matrix([[cos(phi)], [sin(phi)]]) for phi in arange(0.0, 2*pi, 0.1)]
		ellipse = []
		for p in circle:
			X = p
			Y = M + stdev * cholE * X
			ellipse.append([Y[0,0], Y[1,0]])

		" create the ground constraints "
		gndGPAC1Pose = node1.getGndGlobalGPACPose()
		currProfile = Pose(gndGPAC1Pose)
		gndGPAC2Pose = node2.getGndGlobalGPACPose()
		gnd_offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)		

		posture1 = node1.getStableGPACPosture()
		posture2 = node2.getStableGPACPosture()

		hull1 = computeHull(node1, sweep = isSweep)
		hull2 = computeHull(node2, sweep = isSweep)
		hull2_trans = []
		hull2_gnd = []

		posture_trans = []
		pose2Profile = Pose(offset)
		for p in posture2:
			posture_trans.append(pose2Profile.convertLocalOffsetToGlobal(p))

		for p in hull2:
			hull2_trans.append(pose2Profile.convertLocalToGlobal(p))

		posture_gnd = []
		gndProfile = Pose(gnd_offset)
		for p in posture2:
			posture_gnd.append(gndProfile.convertLocalOffsetToGlobal(p))

		for p in hull2:
			hull2_gnd.append(gndProfile.convertLocalToGlobal(p))

		segWidth = self.probe.robotParam['segWidth']
		segLength = self.probe.robotParam['segLength']
				
		for p in posture_gnd:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "black")

		xP = []
		yP = []
		for p in hull2_gnd:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull2_gnd[0][0])
		yP.append(hull2_gnd[0][1])
		pylab.plot(xP,yP, color = "black")

		for p in posture1:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "green")

		xP = []
		yP = []
		for p in hull1:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull1[0][0])
		yP.append(hull1[0][1])
		pylab.plot(xP,yP, color = "green")

		for p in posture_trans:
			xP = []
			yP = []
			
			p1 = [p[0] - 0.5*segWidth*sin(p[2]), p[1] + 0.5*segWidth*cos(p[2])]
			p2 = [p[0] + 0.5*segWidth*sin(p[2]), p[1] - 0.5*segWidth*cos(p[2])]
			p3 = [p[0] + segLength*cos(p[2]) + 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) - 0.5*segWidth*cos(p[2])]
			p4 = [p[0] + segLength*cos(p[2]) - 0.5*segWidth*sin(p[2]), p[1] + segLength*sin(p[2]) + 0.5*segWidth*cos(p[2])]

			xP = [p4[0],p3[0],p2[0],p1[0],p4[0]]
			yP = [p4[1],p3[1],p2[1],p1[1],p4[1]]
			
			pylab.plot(xP,yP, linewidth=1, color = "red")

		xP = []
		yP = []
		for p in hull2_trans:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull2_trans[0][0])
		yP.append(hull2_trans[0][1])
		pylab.plot(xP,yP, color = "red")


		xP = []
		yP = []
		for p in ellipse:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, linewidth=3, color='blue')
	
	def renderConstraints(self):
		
		self.motion_constraints
		self.overlap_constraints
		self.inplace_constraints
		self.gnd_constraints
		
		stdev = 1

		#for k, v in self.nodeHash.items():
		#	self.nodeHash[k].computeAlphaBoundary()

		for k, v in self.nodeHash.items():
			pylab.clf()
			self.renderNode(k, isSweep = True)
			pylab.title("LocalNode %d" % k)
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("localNode_%04u.png" % k)	

		count = 0
		for const in self.motion_constraints:
			pylab.clf()
			self.renderSingleConstraint(const)
			
			pylab.title("Motion Constraint %d ---> %d" % (const[0], const[1]))
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("constraint_%04u.png" % count)
			count += 1
			
		for const in self.overlap_constraints:
			pylab.clf()
			self.renderSingleConstraint(const)
			
			pylab.title("Overlap Constraint %d ---> %d" % (const[0], const[1]))
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("constraint_%04u.png" % count)
			count += 1
			
		for const in self.inplace_constraints:
			pylab.clf()
			self.renderSingleConstraint(const)
			
			pylab.title("In-Place Constraint %d ---> %d" % (const[0], const[1]))
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("constraint_%04u.png" % count)
			count += 1

		for const in self.backtrack_constraints:
			pylab.clf()
			self.renderSingleConstraint(const)
			
			pylab.title("Backtrack Constraint %d ---> %d" % (const[0], const[1]))
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("constraint_%04u.png" % count)
			count += 1
						

		#for const in self.cornerHypotheses:

		#	pylab.clf()
		#	self.renderSingleConstraint(const)
			
		#	pylab.title("Corner Constraint %d ---> %d" % (const[0], const[1]))
		#	pylab.xlim(-4,4)
		#	pylab.ylim(-4,4)
		#	pylab.savefig("constraint_%04u.png" % count)
		#	count += 1

		#f = open("sensor_constraints.txt", 'r')
		#sensor_hypotheses = eval(f.read().rstrip())
		#f.close()

		#f = open("overlap_constraints.txt", 'r')
		#overlap_hypotheses = eval(f.read().rstrip())
		#f.close()

		#for const in overlap_hypotheses:
		#	pylab.clf()
		#	self.renderSingleConstraint(const)
			
		#	pylab.title("Overlap Hypotheses %d ---> %d" % (const[0], const[1]))
		#	pylab.xlim(-4,4)
		#	pylab.ylim(-4,4)
		#	pylab.savefig("constraint_%04u.png" % count)
		#	count += 1
					

		#for const in sensor_hypotheses:
		#	pylab.clf()
		#	self.renderSingleConstraint(const, isSweep = True)
		#	pylab.title("Sensor Hypotheses %d ---> %d" % (const[0], const[1]))
		#	pylab.xlim(-4,4)
		#	pylab.ylim(-4,4)
		#	pylab.savefig("constraint_%04u.png" % count)
		#	count += 1
		
		
		#self.merged_constraints 
		

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
			


	def drawConstraints(self, id = []):
		
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())

		pylab.clf()
		for i in range(self.numNodes):

			hull = computeBareHull(self.nodeHash[i], sweep = False)
			hull.append(hull[0])

			node1 = self.nodeHash[i]
			currPose = currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	

			hull_trans = []
			for p in hull:
				hull_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	

		for k, v in self.edgeHash.items():	
			node1 = self.nodeHash[k[0]]
			node2 = self.nodeHash[k[1]]
			pose1 = node1.getGlobalGPACPose()
			pose2 = node2.getGlobalGPACPose()

			xP = [pose1[0], pose2[0]]
			yP = [pose1[1], pose2[1]]

			#age = self.edgeAgeHash[k]
			age = 0

			fval = float(age) / 20.0
			if fval > 0.4:
				fval = 0.4
			color = (fval, fval, fval)
			pylab.plot(xP, yP, color=color, linewidth=1)


		xP = []
		yP = []
		for pose in poses:
			xP.append(pose[0])
			yP.append(pose[1])		
		pylab.scatter(xP,yP, color='k', linewidth=1)

		#pylab.title("Corner Constraint %d ---> %d" % (id1, id2))
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)

		pylab.title("%d Poses" % self.numNodes)
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)

			
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)

			
		
				