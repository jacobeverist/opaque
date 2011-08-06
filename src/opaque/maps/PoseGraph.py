
from numpy import *
from scipy.optimize import *
import scipy
import scipy.linalg
import numpy
import numpy.linalg
import random

import csv

from Map import Map
from LocalNode import LocalNode
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

import pylab
from matplotlib.patches import Circle

import socket

from math import *
from copy import copy

import Image
import ImageDraw

import cProfile

renderCount = 0

GND_PRIORITY = 10
INPLACE_PRIORITY = 5
CORNER_PRIORITY = 3
SHAPE_PRIORITY = 2
OVERLAP_PRIORITY = 0
MOTION_PRIORITY = 0


def decimatePoints(points):
	result = []

	for i in range(len(points)):
		#if i%2 == 0:
		if i%4 == 0:
			result.append(points[i])

	return result

def computeBareHull(node1, sweep = False):
	
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

def computeHull(node1, sweep = False):
	
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
	node1.boundaryMap.getBoundaryPoints()
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
		self.b_hulls = []
		self.a_hulls = []
		self.a_medial = []

		self.E_gnd = matrix([[ 0.00001, 0.0, 0.0 ],
							[ 0.0, 0.00001, 0.0],
							[ 0.0, 0.0, 0.0001 ]])

		self.E_corner = matrix([[ 0.01, 0.0, 0.0 ],
							[ 0.0, 0.01, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

		self.E_inplace = matrix([[ 0.05, 0.0, 0.0 ],
							[ 0.0, 0.05, 0.0],
							[ 0.0, 0.0, pi/16.0 ]])

		self.E_overlap = matrix([[ 0.1,  0.0, 0.0],
							[ 0.0,  0.01, 0.0],
							[0.0, 0.0,  0.02]])

		#self.E_overlap = matrix([[ 0.2,  0.0, 0.0],
		#					[ 0.0,  0.02, 0.0],
		#					[0.0, 0.0,  0.1]])

		self.E_motion = matrix([[ 0.06, 0.0, 0.0],
							[0.0,  5.0, 0.0],
							[0.0, 0.0,  0.1 ]])
		
		self.E_sensor = matrix([[0.1,0.0,0.0],
							[0.0,0.05,0.0],
							[0.0,0.0,0.02]])

	def addPriorityEdge(self, edge, priority):

		try:
			self.edgePriorityHash[(edge[0],edge[1])]
		except:
			self.edgePriorityHash[(edge[0],edge[1])] = []
		
		" corner constraints with priority of 3 "
		self.edgePriorityHash[(edge[0],edge[1])].append([edge[2], edge[3], priority])

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
		
	def matchExternalCorners(self):
		
		externalHyps = []
		
		" for every bin, attempt to match a random corner with every other bin "
		for i in range(len(self.cornerBins)):
			bin1 = self.cornerBins[i]

			tup1 = bin1[random.randint(0,len(bin1))]
			
			for j in range(i+1):
				bin2 = self.cornerBins[j]
				tup2 = bin2[random.randint(0,len(bin2))]
							
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

		return externalHyps

	def makeCornerBinConsistent(self):

		for k in range(len(self.cornerBins)):
			bin = self.cornerBins[k]
			attempts = []
			for i in range(len(bin)):
				tup1 = bin[i]
				for j in range(i+1, len(bin)):
					tup2 = bin[j]
					if tup1[0] != tup2[0]:
						attempts.append([tup1[0],tup2[0],tup1[1],tup2[1]])

			print "BIN", k
			internalHyps = []
			for attempt in attempts:
				tup1 = attempt[0]
				tup2 = attempt[1]
				n1 = attempt[0]
				n2 = attempt[1]
				point1_i = attempt[2]
				point2_i = attempt[3]
				corner1 = self.nodeHash[n1].cornerCandidates[point1_i]
				corner2 = self.nodeHash[n2].cornerCandidates[point2_i]
				point1 = corner1[0]
				point2 = corner2[0]
				ang1 = corner1[1]
				ang2 = corner2[1]
	
				offset, covar, cost, angDiff, oriDiff = self.makeCornerConstraint(n1, n2, point1, point2, ang1, ang2, self.a_hulls[n1], self.a_hulls[n2], self.b_hulls[n1], self.b_hulls[n2])
				print "CONSTRAINT:", tup1, tup2, cost, angDiff, oriDiff
				transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
				if cost <= 0.6:
					internalHyps.append([n1, n2, transform, covar, cost, point1, point2, angDiff, oriDiff, point1_i, point2_i])

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

	def loadNewNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))

		if nodeID < 26:
			direction = True
		elif nodeID < 34:
			direction = False
		elif nodeID < 38:
			direction = True
		elif nodeID < 42:
			direction = False
		elif nodeID < 46:
			direction = True
		elif nodeID < 60:
			direction = False
		elif nodeID < 64:
			direction = True
		else:
			direction = False

		

		print "direction =", direction
		#return

		if nodeID > 0:
			if nodeID % 2 == 0:
				" 2nd node from a pose "
				" 1) make overlap+motion constraint with i-1"
				" 2) attempt overlap with past nodes ? "
				" 3) attempt corner constraints with past nodes, replacing overlaps and sensors "
				" 4) discard sensor constraints and reapply sensor constraints to non-corner edges"

				#transform, covE = self.makeMotionConstraint(nodeID-1, nodeID)
				#self.edgePriorityHash[(nodeID-1, nodeID)].append([transform, covE, MOTION_PRIORITY])
				
				#transform, covE = self.makeOverlapConstraint(nodeID-2, nodeID)
				#transform, covE = self.makeOverlapConstraint2(nodeID-2, nodeID)
				transform, covE = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )

				

	
				#def wrapMedial(self, i,j,result):
				#		result.append(self.makeMedialOverlapConstraint(i,j))

				#result = []
				#cProfile.run('self.wrapMedial(nodeID-2, nodeID, result)', 'prof_result_%04u_%04u' % (nodeID-2, nodeID))
				#transform = result[0]
				#covE = result[1]
				

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
				
				" try to perform sensor constraints on past nodes "
				" overwrite overlap/motion constraints "
				" do not overwrite corner constraints "


				"""
				paths = []
				for i in range(self.numNodes):
					paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
				
				newConstraints = []
				if nodeID >= 2 and nodeID <= 24:
					newConstraints = self.addShapeConstraints(paths, targetNode = nodeID)

				for const in newConstraints:
					n1 = const[0]
					n2 = const[1]
					transform = const[2]
					covE = const[3]

					self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)

				" merge the constraints "
				self.mergePriorityConstraints()
				"""
				
				" remove all sensor constraints from hash "
				#deleteAll(self.edgePriorityHash, 3)
				
				"""
				
				if nodeID >= 2:
					if self.a_hulls[nodeID-2] == 0:
						self.a_hulls[nodeID-2] = computeHull(self.nodeHash[nodeID-2], sweep = False)
						self.b_hulls[nodeID-2] = computeHull(self.nodeHash[nodeID-2], sweep = True)
				if nodeID >= 4:
					if self.a_hulls[nodeID-4] == 0:
						self.a_hulls[nodeID-4] = computeHull(self.nodeHash[nodeID-4], sweep = False)
						self.b_hulls[nodeID-4] = computeHull(self.nodeHash[nodeID-4], sweep = True)
				if nodeID >= 6:
					if self.a_hulls[nodeID-6] == 0:
						self.a_hulls[nodeID-6] = computeHull(self.nodeHash[nodeID-6], sweep = False)
						self.b_hulls[nodeID-6] = computeHull(self.nodeHash[nodeID-6], sweep = True)

				newConstraints = []
				if nodeID >= 2:
					offset, covar, cost = self.makeSensorConstraint(paths[nodeID][nodeID-2][0], nodeID, nodeID-2, self.b_hulls[nodeID], self.a_hulls[nodeID-2])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost < 1.5:				
						newConstraints.append([nodeID,nodeID-2,transform,covar])			

				if nodeID >= 4:
					offset, covar, cost = self.makeSensorConstraint(paths[nodeID][nodeID-4][0], nodeID, nodeID-2, self.b_hulls[nodeID], self.a_hulls[nodeID-4])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost < 1.5:				
						newConstraints.append([nodeID,nodeID-4,transform,covar])			

				if nodeID >= 6:
					offset, covar, cost = self.makeSensorConstraint(paths[nodeID][nodeID-6][0], nodeID, nodeID-2, self.b_hulls[nodeID], self.a_hulls[nodeID-6])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost < 1.5:				
						newConstraints.append([nodeID,nodeID-6,transform,covar])			
				"""
					
	
			else:

						
				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY)

				if nodeID > 2:
					#transform, covE = self.makeOverlapConstraint2(nodeID-2, nodeID)
					transform, covE = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction)
					self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)

	
				" merge constraints by priority and updated estimated poses "
				self.mergePriorityConstraints()

				" try to perform corner constraints "
				" overwrite overlap/motion constraints"
				" corner constraints with past nodes "
				self.addCornerConstraints(nodeID)

				self.mergePriorityConstraints()
				
				" 1) add inplace constraint "
				" 2) attempt corner constraint with past nodes, replacing overlaps and sensors "
				" 3) apply sensor constraints to non-corner edges "				

				"""
				paths = []
				for i in range(self.numNodes):
					paths.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))

				newConstraints = []
				if nodeID >= 2 and nodeID <= 24:
					newConstraints = self.addShapeConstraints(paths, targetNode = nodeID)

				for const in newConstraints:
					n1 = const[0]
					n2 = const[1]
					transform = const[2]
					covE = const[3]

					self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)

				" merge the constraints "
				self.mergePriorityConstraints()
				"""

				"""
				
				if nodeID >= 2:
					if self.a_hulls[nodeID-2] == 0:
						self.a_hulls[nodeID-2] = computeHull(self.nodeHash[nodeID-2], sweep = False)
						self.b_hulls[nodeID-2] = computeHull(self.nodeHash[nodeID-2], sweep = True)
				if nodeID >= 4:
					if self.a_hulls[nodeID-4] == 0:
						self.a_hulls[nodeID-4] = computeHull(self.nodeHash[nodeID-4], sweep = False)
						self.b_hulls[nodeID-4] = computeHull(self.nodeHash[nodeID-4], sweep = True)
				if nodeID >= 6:
					if self.a_hulls[nodeID-6] == 0:
						self.a_hulls[nodeID-6] = computeHull(self.nodeHash[nodeID-6], sweep = False)
						self.b_hulls[nodeID-6] = computeHull(self.nodeHash[nodeID-6], sweep = True)

				newConstraints = []
				if nodeID >= 2:
					offset, covar, cost = self.makeSensorConstraint(paths[nodeID][nodeID-2][0], nodeID, nodeID-2, self.b_hulls[nodeID], self.a_hulls[nodeID-2])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost < 1.5:				
						newConstraints.append([nodeID,nodeID-2,transform,covar])			

				if nodeID >= 4:
					offset, covar, cost = self.makeSensorConstraint(paths[nodeID][nodeID-4][0], nodeID, nodeID-2, self.b_hulls[nodeID], self.a_hulls[nodeID-4])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost < 1.5:				
						newConstraints.append([nodeID,nodeID-4,transform,covar])			

				if nodeID >= 6:
					offset, covar, cost = self.makeSensorConstraint(paths[nodeID][nodeID-6][0], nodeID, nodeID-2, self.b_hulls[nodeID], self.a_hulls[nodeID-6])
					transform = matrix([ [offset[0]], [offset[1]], [offset[2]] ],dtype=float)
					if cost < 1.5:				
						newConstraints.append([nodeID,nodeID-6,transform,covar])		
						
				newConstraints = self.addShapeConstraints(paths, targetNode = nodeID)
				"""

		"""
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
		"""

		#if nodeID == 8:
		#if nodeID % 15 == 0 and nodeID > 0:
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
				

	

		
		"""
		#if nodeID == 26:
		if nodeID % 25 == 1 and nodeID > 1:

			self.deleteAllPriority(SHAPE_PRIORITY)
			for const in self.newConstraints2:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

		if nodeID % 50 == 2 and nodeID > 2:

			self.deleteAllPriority(SHAPE_PRIORITY)
			for const in self.newConstraints3:
				n1 = const[0]
				n2 = const[1]
				transform = const[2]
				covE = const[3]
	
				self.addPriorityEdge([n1,n2,transform,covE], SHAPE_PRIORITY)
	
			" merge the constraints "
			self.mergePriorityConstraints()

		print "Cluster Error:"
		totalError = self.computeCornerClusterError()
		"""

	def correctNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		self.a_hulls.append(computeHull(self.nodeHash[nodeID], sweep = False))
		self.b_hulls.append(computeHull(self.nodeHash[nodeID], sweep = True))

		self.a_medial.append(self.nodeHash[nodeID].getMedialAxis(sweep = False))
		
		
		if nodeID > 0:
			if nodeID % 2 == 0:

				direction = self.currNode.travelDir

				transform, covE = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction )

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
	
				transform, covE = self.makeInPlaceConstraint(nodeID-1, nodeID)
				self.addPriorityEdge([nodeID-1,nodeID,transform,covE], INPLACE_PRIORITY)

				if nodeID > 2:
					#transform, covE = self.makeOverlapConstraint2(nodeID-2, nodeID)
					transform, covE = self.makeMedialOverlapConstraint(nodeID-2, nodeID, isMove = True, isForward = direction)
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


	def makeMedialOverlapConstraint(self, i, j, isMove = True, isForward = True ):

		print "recomputing hulls and medial axis"
		" compute the medial axis for each pose "
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
		
		estPose1 = node1.getGlobalGPACPose()		
		estPose2 = node2.getGlobalGPACPose()
		
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		medialSpline2 = SplineFit(medial2, smooth=0.1)


		#samples = scipy.arange(0.0,1.0,0.01)

		if isMove:
			if isForward:
				u2 = 0.6
			else:
				u2 = 0.4
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

		print "u2 =", u2


		#motionT, motionCov = self.makeMotionConstraint(i,j)
		#travelDelta = motionT[0,0]
		
		u1 = 0.5
		#u2 = 0.6
		angGuess = 0.0

		result = gen_icp.overlapICP(estPose1, [u1, u2, angGuess], hull1, hull2, medial1, medial2, plotIter = False, n1 = i, n2 = j)

		transform = matrix([[result[0]], [result[1]], [result[2]]])
		covE =  self.E_overlap
		
		print "making overlap constraint:", result[0], result[1], result[2]
		
		return transform, covE				


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

	def makeInPlaceConstraint(self, i, j):

		" compute the new posture "
		" 1. compute the GPAC pose "
		" 2. perform correction of GPAC pose "
		" 3. compute location of rootPose from corrected GPAC pose "
		
		" node1 is the front poke node "
		" nodes is the back poke node "

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
		
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
		
		if foreCost > backCost:
			#plotPosture(originBackPosture, localBackPosture)
			correctedGPACPose = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		else:
			#plotPosture(originForePosture, localForePosture)					
			correctedGPACPose = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])

		correctedProfile = Pose(correctedGPACPose)
	
		" offset from 0,0 to newGPACPose "
		" rootPose from neg. offset from correctedGPACPose "
		localRootOffset3 = correctedProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		
		
		if foreCost > backCost:
			offset = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3)
		else:
			offset = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3)
		
		transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
		covE = self.E_inplace
		
		return transform, covE

	def addInPlaceNode(self, newNode):
		
		self.currNode = newNode
		
		self.nodeHash[self.numNodes] = self.currNode
		if self.numNodes > 0:
	
			node1 = self.nodeHash[self.numNodes-1]
			node2 = self.nodeHash[self.numNodes]
				
			transform, covE = self.makeInPlaceConstraint(self.numNodes-1, self.numNodes)
			
			self.inplace_constraints.append([self.numNodes-1, self.numNodes, transform, covE])	


			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(self.nodeHash[self.numNodes-1].getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			self.nodeHash[self.numNodes].setGPACPose(currPose)
						
		self.saveConstraints()
		
		self.numNodes += 1
			
	def addNode(self, newNode):
		
		self.currNode = newNode

		self.nodeHash[self.numNodes] = self.currNode
		if self.numNodes > 0:
						
			transform, covE = self.makeMotionConstraint(self.numNodes-1, self.numNodes)
			self.motion_constraints.append([self.numNodes-1, self.numNodes, transform, covE])
			
			transform, covE = self.makeOverlapConstraint(self.numNodes-1, self.numNodes)
			self.overlap_constraints.append([self.numNodes-1, self.numNodes, transform, covE])

			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(self.nodeHash[self.numNodes-1].getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			self.nodeHash[self.numNodes].setGPACPose(currPose)
		
			self.edgeHash[(0,self.numNodes)] = [transform, covE]				
		
		self.saveConstraints()
		
		self.numNodes += 1
		
	def addInitNode(self, newNode, estPose):
		
		self.currNode = newNode

		self.nodeHash[self.numNodes] = self.currNode
		if self.numNodes > 0:

			transform, covE = self.makeOverlapConstraint(0, self.numNodes)
			self.overlap_constraints.append([0, self.numNodes, transform, covE])

			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(self.nodeHash[0].getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			self.nodeHash[self.numNodes].setGPACPose(currPose)

			self.edgeHash[(0,self.numNodes)] = [transform, covE]
		
		self.saveConstraints()
		
		self.numNodes += 1

	
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
				transform, covE = self.makeMedialOverlapConstraint(n1, n2, isMove = False)
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
			
		f = open("motion_constraints_%04u.txt" % self.numNodes, 'w')
		f.write(repr(self.motion_constraints))
		f.write("\n")
		f.close()		

		f = open("sensor_constraints_%04u.txt" % self.numNodes, 'w')
		f.write(repr(self.sensor_constraints))
		f.write("\n")
		f.close()		

		f = open("inplace_constraints_%04u.txt" % self.numNodes, 'w')
		f.write(repr(self.inplace_constraints))
		f.write("\n")
		f.close()		

		f = open("backtrack_constraints_%04u.txt" % self.numNodes, 'w')
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
		pylab.savefig("constraint_%04u.png" % count)

	
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

			
		
				
