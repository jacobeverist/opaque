
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
from VoronoiMap import VoronoiMap
from FrontierMap import FrontierMap
from OccupancyMap import OccupancyMap
from FreeSpaceBoundaryMap import FreeSpaceBoundaryMap
from NavRoadMap import NavRoadMap
from StablePose import StablePose
from StableCurve import StableCurve
from Pose import Pose
import gen_icp
from functions import *
import toro
import scsgp
import bayes
import estimateMotion

import pylab
from matplotlib.patches import Circle

from math import *
from copy import copy

import Image
import ImageDraw

# Map Space Parameters
PIXELSIZE = 0.05
MAPSIZE = 40.0
#MAPSIZE = 20.0


def decimatePoints(points):
	result = []

	for i in range(len(points)):
		#if i%2 == 0:
		if i%4 == 0:
			result.append(points[i])

	return result

def computeHull(node1, sweep = False):
	
	" Read in data of Alpha-Shapes and add their associated covariances "
	node1.computeAlphaBoundary(sweep = sweep)
	a_data = node1.getAlphaBoundary()
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


class MapGraph:

	def __init__(self, probe, contacts, stabilizePose = False):
		
		#
		self.stabilizePose = stabilizePose
		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)
		
		self.initPose = self.probe.getActualJointPose(19)
		self.nodeHash = {}
		self.edgeHash = {}
		
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
		
		self.overlap_constraints = []
		self.motion_constraints = []
		self.sensor_constraints = []
		self.inplace_constraints = []
		self.gnd_constraints = []		
		self.merged_constraints = []
		
				
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

	def loadFile(self, dirName, num_poses):
		
		self.numNodes = 0
		self.currNode = 0
		
		self.inplace_constraints = []
		self.overlap_constraints = []
		self.motion_constraints = []
		self.sensor_constraints = []	
		self.gnd_constraints = []
		self.merged_constraints = []

		" error from average "
		E_motion = [[ 0.37854806, -0.0049439,  -0.01498673],
			[-0.0049439,   0.08228017, -0.01217169],
			[-0.01498673, -0.01217169,  0.03889908]]
		E_overlap = [[ 0.63050944,  0.00400762, -0.03559462],
			[ 0.00400762,  0.01275507, -0.00786715],
			[-0.03559462, -0.00786715,  0.02484153]]

		" positional error "
		E_motion = [[ 0.05582791, -0.00422392, -0.00551943],
			[-0.00422392,  0.06441896, -0.00346319],
			[-0.00551943, -0.00346319,  0.0218136 ]]
		E_overlap = [[ 0.52709836, -0.033402,   -0.03137499],
			[-0.033402,    0.01522771, -0.00650399],
			[-0.03137499, -0.00650399,  0.01903969]]	

		" theoretical "		
		E_motion = [[ 0.06, 0.0, 0.0],
			[0.0,  5.0, 0.0],
			[0.0, 0.0,  0.1 ]]
		#E_overlap = [[ 2.0, 0.0, 0.0],
		#	[0.0, 0.1, 0.0],
		#	[0.0, 0.0, 0.02]]	

		E_overlap = matrix([[ 0.2,  0.0, 0.0],
		        [ 0.0,  0.02, 0.0],
		        [0.0, 0.0,  0.1]])
		
		for i in range(0,num_poses):
			print "loading node", i			
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, self.pixelSize)
			self.currNode.readFromFile(dirName, i)
			self.nodeHash[i] = self.currNode

			if self.numNodes > 0:
				
				node1 = self.nodeHash[i-1]
				node2 = self.nodeHash[i]
		
				
				" create the ground constraints "
				gndGPAC1Pose = node1.getGndGlobalGPACPose()
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = node2.getGndGlobalGPACPose()
				offset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
				self.gnd_constraints.append([matrix([[offset[0]], [offset[1]], [offset[2]]]), matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])    ])

			self.numNodes += 1

		" load the constraints from file "
		index = num_poses - 1
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
		
		
		for const in self.overlap_constraints:
			const[3] =  E_overlap
		

		#print "covariance 1:"
		#E_motion, E_overlap = self.computeCovar()		
		#print E_motion
		#print E_overlap
		
		#print "covariance 2:"
		#E_motion, E_overlap = self.computeCovar2()		
		#print E_motion
		#print E_overlap

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = E_motion
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			#covE = E_overlap			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		print "totalConstraints:"
		print totalConstraints
		
		" merge the constraints with TORO"
		self.merged_constraints = self.toroMerge(totalConstraints)
		
		print "merged_constraints:", self.merged_constraints
		
		self.edgeHash = {}
		for i in range(len(self.merged_constraints)):
			print i, self.merged_constraints[i]
			node1 = self.merged_constraints[i][0]
			node2 = self.merged_constraints[i][1]
			self.edgeHash[(node1, node2)] = [self.merged_constraints[i][2], self.merged_constraints[i][3]]
			
		
		for i in range(self.numNodes-1):
			node1 = self.nodeHash[i]
			node2 = self.nodeHash[i+1]
			transform, covE = self.edgeHash[(i, i+1)]
			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(node1.getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)
			node2.setGPACPose(currPose)		

		#self.drawConstraints()

		" save the maps "		
		#self.synch()
		#self.saveMap()

		

	def computeCovar(self):

		" compute the average position, and the covariance of position "

		N = len(self.gnd_constraints)

		for i in range(N):
			gndConst = self.gnd_constraints[i][0]
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
			gndConst = self.gnd_constraints[i][0]
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
			gndConst = self.gnd_constraints[i][0]
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
			gndConst = self.gnd_constraints[i][0]
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
			gndConst = self.gnd_constraints[i][0]
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
			gndConst = self.gnd_constraints[i][0]
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
			
			gndConst = self.gnd_constraints[i][0]
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

		E_overlap = matrix([[ 0.2,  0.0, 0.0],
		        [ 0.0,  0.02, 0.0],
		        [0.0, 0.0,  0.1]])

		#def makeSensorConstraint(self, transform, n1, n2, hull1, hull2):
		
		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]
			
		curve1 = node1.getGPACurve()
		curve2 = node2.getGPACurve()

		offset = [transform[0,0],transform[1,0],transform[2,0]]


		result = gen_icp.motionICP(curve1.getUniformSamples(), curve2.getUniformSamples(), offset, plotIter = False, n1 = i, n2 = j)

		" update the pose with the motion estimation "
		offset = [result[0], result[1], result[2]]
		
		" add covariance from vector in guess direction "
		transform = matrix([[result[0]], [result[1]], [result[2]]])
		#covE = matrix([[ 0.4,  0.0, 0.0],
		#        [ 0.0,  0.1, 0.0],
		#        [0.0, 0.0,  0.1]])
		covE = E_overlap
		
		return transform, covE
				
	def makeOverlapConstraint(self, i, j ):	
		

		E_overlap = matrix([[ 0.2,  0.0, 0.0],
		        [ 0.0,  0.02, 0.0],
		        [0.0, 0.0,  0.1]])
				
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

		" update the pose with the motion estimation "
		#offset = [result[0], result[1], result[2]]
		#currProfile = Pose(node1.getGlobalGPACPose())
		#currPose = currProfile.convertLocalOffsetToGlobal(offset)	
		#node2.setGPACPose(currPose)
		
		" FIXME:  add covariance from vector in guess direction "
		
		transform = matrix([[result[0]], [result[1]], [result[2]]])
		#covE = matrix([[ 0.4,  0.0, 0.0],
		#        [ 0.0,  0.1, 0.0],
		#        [0.0, 0.0,  0.1]])
		covE =  E_overlap
		
		print "adding overlap constraint:", result[0], result[1], result[2]
		
		return transform, covE
	
	
	
	def makeMotionConstraint(self, i, j):
		E_motion = matrix([[ 0.06, 0.0, 0.0],
			[0.0,  5.0, 0.0],
			[0.0, 0.0,  0.1 ]])
		
		if self.numNodes > 0:		
			
			print "adding motion constraint", i, j
	
			node1 = self.nodeHash[i]
			node2 = self.nodeHash[j]
			pose1 = node1.getGlobalGPACPose()
			pose2 = node2.getGlobalGPACPose()
			
			profile1 = Pose(pose1)
			profile2 = Pose(pose2)

			localOffset = profile1.convertGlobalPoseToLocal(pose2)

			xT = localOffset[0]
			yT = localOffset[1]
			pT = localOffset[2]
			
			transform = matrix([[xT], [yT], [pT]])
			
			#covE = matrix([[ 0.1, 0.0, 0.0 ],
			#				[ 0.0, 0.1, 0.0],
			#				[ 0.0, -0.0, pi/4.0 ]])

			covE = E_motion
			
			return transform, covE
					
		
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
		#minMatchDist = 2.0
		minMatchDist = 0.5
	
		" plot the best fit at each iteration of the algorithm? "
		plotIteration = True
		#plotIteration = False
	
		#offset = [0.0,0.0,0.0]
	
		" Extract the data from the files and put them into arrays "
		
		node1 = self.nodeHash[n1]
		node2 = self.nodeHash[n2]		
		
		" FIXME:  deal with relative poses from the motion constraints "
		estPose1 = node1.getGlobalGPACPose()
		estPose2 = node2.getGlobalGPACPose()
		#estPose1 = node1.getEstPose()
		#estPose2 = node2.getEstPose()

		hull_trans1 = []
		for p in hull1:
			hull_trans1.append(gen_icp.dispPoint(p, estPose1))

		hull_trans2 = []
		for p in hull2:
			hull_trans2.append(gen_icp.dispPoint(p, estPose2))

		#offset2 = [transform[0,0],transform[1,0],transform[2,0]]
		
		radius, center = gen_icp.computeEnclosingCircle(hull1)
		circle1 = [radius, center]
		
		offset, cost = gen_icp.gen_ICP2(estPose1, firstGuess, hull1, hull2, [circle1], costThresh, minMatchDist, plotIteration, n1, n2)
		print "%d -> %d" %(n1,n2), "cost =", cost
		covar = matrix([[0.1,0.0,0.0],[0.0,0.1,0.0],[0.0,0.0,0.3]])
		
		return offset, covar, cost


	def makeInPlaceConstraint(self, i, j):

		" compute the new posture "
		" 1. compute the GPAC pose "
		" 2. perform correction of GPAC pose "
		" 3. compute location of rootPose from corrected GPAC pose "
		
		" node1 is the front poke node "
		" nodes is the back poke node "

		node1 = self.nodeHash[i]
		node2 = self.nodeHash[j]



		originPosture = node1.localPosture
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

		newPosture = []
		for j in range(self.probe.numSegs-1):
			newPosture.append(self.probe.getJointPose([0.0,0.0,0.0], node1.rootNode, j))


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
		covE = matrix([[ 0.05, 0.0, 0.0 ],
						[ 0.0, 0.05, 0.0],
						[ 0.0, 0.0, pi/16.0 ]])
		
		return transform, covE

	def getCurrentNode(self):
		return self.currNode

	def newInPlaceNode(self, direction):
		
		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.stabilizePose, direction = direction)
		
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
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

		print "checkR"
	
	
	def newNode(self, stepDist, direction):
		
		print "checkN"
		
		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		print "checkO"

		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19)
		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, inSim = False)
		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.stabilizePose, direction = True)

		print "checkP"
		
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
			#self.addNaiveMotionConstraint(self.numNodes-1, self.numNodes, stepDist, direction)
			"""
			node1 = self.poseGraph.get_node_attributes(self.numNodes-1)
			node2 = self.poseGraph.get_node_attributes(self.numNodes)
		
			curve1 = node1.getGPACurve()
			curve2 = node2.getGPACurve()
			posture1 = node1.getGPACPosture()
			posture2 = node2.getGPACPosture()
	
			#(points1, points2, offset, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):
			offset = estimateMotion.makeGuess2(curve1, curve2, 0.0, posture1, posture2)
			result = gen_icp.motionICP(curve1.getUniformSamples(), curve2.getUniformSamples(), offset, plotIter = False, n1 = (self.numNodes-1), n2 = self.numNodes)
	
			" update the pose with the motion estimation "
			offset = [result[0], result[1], result[2]]
			currProfile = Pose(node1.getGlobalGPACPose())
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			node2.setGPACPose(currPose)
	
				#result = constraints[i-1]
	
			self.gpacs.append([0.2, result])	
			#self.gpacs.append([0.2, offset])	
			constraints.append(result)
			
			" FIXME:  add covariance from vector in guess direction "
			
			transform = matrix([[result[0]], [result[1]], [result[2]]])
			#transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
			covE = matrix([[ 0.4,  0.0, 0.0],
			        [ 0.0,  0.1, 0.0],
			        [0.0, 0.0,  0.1]])
			print "adding motion constraint:", result[0], result[1], result[2]
			#print "adding motion constraint:", offset[0], offset[1], offset[2]
			self.poseGraph.add_edge(i-1, i, attrs = [transform, covE])
			"""
			
		print "checkQ"
		
		self.saveConstraints()
		
		self.numNodes += 1
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

		print "checkR"

	def initNode(self, estPose, direction):

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.stabilizePose, direction = True)

		print "checkP"
		
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
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
	
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

	def findOverlapCandidates(self, paths, targetNode = -1):

		SENSOR_RADIUS = 0.0
		DIST_THRESHOLD = 6.0
		
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
		
		print "UNIQUE PAIRS"
		for p in pair_unique:
			i = p[1]
			j = p[2]
			
			
			loc2 = paths[i][j][0]
			c1 = matrix([[0.0],[0.0]],dtype=float)
			c2 = matrix([[loc2[0,0]],[loc2[1,0]]],dtype=float)
		
			E_param = paths[i][j][1]
			E = matrix([[E_param[0,0], E_param[0,1]],[E_param[1,0], E_param[1,1]]],dtype=float)
			
			dist = bayes.mahab_dist(c1, c2, SENSOR_RADIUS, SENSOR_RADIUS, E)
			
			print "distance:", i, j, dist

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
			#self.synch()

		
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
		added_hypotheses = self.addOverlaps(paths1)

		" merge and add them to the graph "
		new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])
		
		self.edgeHash = {}
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		#self.drawConstraints()


		" save the maps "		
		#self.synch()
		#self.saveMap()

		" perform dijkstra projection over all nodes"
		paths2 = []
		for i in range(self.numNodes):
			paths2.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
		
		edgeHash2 = self.edgeHash
		added_hypotheses = self.addSensorConstraints(paths2)

		new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])

		self.edgeHash = {}		
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		#self.drawConstraints()


		" save the maps "		
		#self.synch()
		#self.saveMap()

	
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


	def relaxCorrect(self):

		" take the current graph, relax the constraints and start performing new hypotheses to try and create more consistent map "

		if self.numNodes < 2:
			return

		if self.currNode.isDirty():
			self.currNode.synch()

		" theoretical "		
		E_motion = [[ 0.06, 0.0, 0.0],
			[0.0,  5.0, 0.0],
			[0.0, 0.0,  0.1 ]]
		E_overlap = [[ 2.0, 0.0, 0.0],
			[0.0, 0.1, 0.0],
			[0.0, 0.0, 0.02]]	
		
		" edges between i and i+1 where i is even should remain unrelaxed "

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = E_motion
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = E_overlap			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		" merge the constraints with TORO"
		self.merged_constraints = self.toroMerge(totalConstraints)

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
		new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])
		
		self.edgeHash = {}
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		#self.drawConstraints()


		" save the maps "		
		#self.synch()
		#self.saveMap()

		" perform dijkstra projection over all nodes"
		paths2 = []
		for i in range(self.numNodes):
			paths2.append(bayes.dijkstra_proj(i, self.numNodes, self.edgeHash))
		
		edgeHash2 = self.edgeHash
		added_hypotheses = self.addSensorConstraints(paths2)

		new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])

		self.edgeHash = {}		
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
		#self.drawConstraints()


		" save the maps "		
		#self.synch()
		#self.saveMap()

	
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
		
	def localizeCurrentNode(self):

		" take the current graph, relax the constraints and start performing new hypotheses to try and create more consistent map "

		if self.numNodes < 2:	
			return

		if self.currNode.isDirty():
			self.currNode.synch()

		" theoretical "		
		E_motion = [[ 0.06, 0.0, 0.0],
			[0.0,  5.0, 0.0],
			[0.0, 0.0,  0.1 ]]
		E_overlap = [[ 2.0, 0.0, 0.0],
			[0.0, 0.1, 0.0],
			[0.0, 0.0, 0.02]]	
		
		" edges between i and i+1 where i is even should remain unrelaxed "

		totalConstraints = []
		for i in range(len(self.motion_constraints)):
			nodeID1 = self.motion_constraints[i][0]
			nodeID2 = self.motion_constraints[i][1]
			transform = self.motion_constraints[i][2]
			covE = E_motion
			totalConstraints.append([nodeID1,nodeID2,transform,covE])
		
		for i in range(len(self.overlap_constraints)):
			nodeID1 = self.overlap_constraints[i][0]
			nodeID2 = self.overlap_constraints[i][1]
			transform = self.overlap_constraints[i][2]
			covE = E_overlap			
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		for i in range(len(self.inplace_constraints)):
			nodeID1 = self.inplace_constraints[i][0]
			nodeID2 = self.inplace_constraints[i][1]
			transform = self.inplace_constraints[i][2]
			covE = self.inplace_constraints[i][3]
			totalConstraints.append([nodeID1,nodeID2,transform,covE])

		" merge the constraints with TORO"
		self.merged_constraints = self.toroMerge(totalConstraints)

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
		added_hypotheses = self.addOverlaps(paths1, targetNode = self.numNodes-1)

		print "CONSTRAINTS BEFORE doToro()"

		print "added_hypotheses:", added_hypotheses
		print "self.merged_constraints:", self.merged_constraints

		" merge and add them to the graph "
		new_constraints = self.doToro(added_hypotheses, self.merged_constraints, [])
		
		self.edgeHash = {}
		for i in range(len(new_constraints)):
			node1 = new_constraints[i][0]
			node2 = new_constraints[i][1]
			self.edgeHash[(node1,node2)] = [new_constraints[i][2], new_constraints[i][3]]
		
			" FIXME:  save this new overlap constraint "
			self.overlap_constraints.append([node1, node2, new_constraints[i][2], new_constraints[i][3]])
		
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
		print "INDICATOR:", w
		
		selected = []	
		for i in range(len(hypotheses)):
			if w[i,0] >= 1.0:
				selected.append(hypotheses[i])
		
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
		print len(selected), "added hypotheses"
		print len(self.merged_constraints), "merged constraints"

		f = open("overlap_constraints.txt", 'w')
		#f.write(repr(added_hypotheses))
		f.write(repr(selected))
		f.close()
		
		#return added_hypotheses
		return selected

	def addSensorConstraints(self, paths):

		" find candidates to perform sensor constraints "
		sensor_pairs = self.findSweepCandidates(paths)

		" compute the alpha hulls of the selected candidates "
		hull_computed = [False for i in range(self.numNodes)]
		a_hulls = [0 for i in range(self.numNodes)]
		
		for i in range(self.numNodes):
			a_hulls[i] = computeHull(self.nodeHash[i], sweep = True)
			hull_computed[i] = True
		
		for p in sensor_pairs:
			if not hull_computed[p[1]]:
				a_hulls[p[1]] = computeHull(self.nodeHash[p[1]], sweep = True)
				hull_computed[p[1]] = True

			if not hull_computed[p[2]]:
				a_hulls[p[2]] = computeHull(self.nodeHash[p[2]], sweep = True)
				hull_computed[p[2]] = True
		
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
		print "INDICATOR:", w
		selected = []	
		for i in range(len(hypotheses)):
			if w[i,0] >= 1.0:
				selected.append(hypotheses[i])

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
		
		print len(selected), "added hypotheses"
		print len(self.merged_constraints), "merged constraints"

		f = open("sensor_constraints.txt", 'w')
		#f.write(repr(added_hypotheses))
		f.write(repr(selected))
		f.close()
		
		#return added_hypotheses
		return selected

		
	def toroMerge(self, constraints):
		
		" vertices "
		v_list = []

		for i in range(self.numNodes):
			node1 = self.nodeHash[i]
			estPose1 = node1.getGlobalGPACPose()
			v_list.append([i, [estPose1[0], estPose1[1], estPose1[2]]])

		
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
				print "edge", node1, node2			
		
		print len(e_list), "constraints"
		" add hypotheses to the TORO file "

		toro.writeConstrainedToroGraph(".", "probe_init.graph", v_list, e_list)
		toro.executeToro("./" + "probe_init.graph")	
		finalFileName = "probe_init-treeopt-final.graph"
		v_list2, e_list2 = toro.readToroGraph("./" + finalFileName)		

		for item in v_list2:
			nodeID = item[0]
			pose1 = item[1]
			node1 = self.nodeHash[nodeID]
			node1.setGPACPose(pose1)
			
		merged_constraints = []
		
		for edge in e_list2:			
			node1 = edge[0]
			node2 = edge[1]
			offset = edge[2]
			prec = edge[3]
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
			covE = scipy.mat(scipy.linalg.inv(invMat))				

			merged_constraints.append([node1, node2, transform, covE])
		
		" return corrected graph "
		return merged_constraints

		
	def doToro(self, hypotheses, constraints, hulls = []):
		
		" vertices "
		v_list = []

		for i in range(self.numNodes):
			node1 = self.nodeHash[i]
			estPose1 = node1.getGlobalGPACPose()
			v_list.append([i, [estPose1[0], estPose1[1], estPose1[2]]])

		
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
		
		print len(e_list), "constraints"
		" add hypotheses to the TORO file "

		
		for k in range(len(hypotheses)):
			hyp = hypotheses[k]
			
			n1 = hyp[0]
			n2 = hyp[1]
			
			transform = hyp[2]
			covE = hyp[3]
			
			offset = [transform[0,0], transform[1,0], transform[2,0]]
			
			invMat = scipy.linalg.inv(covE)				
			precision = [invMat[0,0], invMat[0,1], invMat[1,1], invMat[2,2], invMat[0,2], invMat[1,2]]
			
			e_list.append([n1,n2,offset,precision])
			
		#print "e_list:"
		#print e_list
		
		pylab.clf()
		#toro.writeToroGraph(".", "probe.graph", v_list, e_list, [])
		toro.writeConstrainedToroGraph(".", "probe.graph", v_list, e_list)
		toro.plotToro(v_list, e_list, hulls, drawedge=True)
		pylab.savefig("uncorrected.png")
		pylab.clf()

		#exit()

		toro.executeToro("./" + "probe.graph")
	
		finalFileName = "probe-treeopt-final.graph"
		v_list2, e_list2 = toro.readToroGraph("./" + finalFileName)		
		
		toro.plotToro(v_list2, e_list2, hulls, drawedge=True)
		pylab.savefig("corrected.png")
		pylab.clf()
		
		" render uncorrect graph "
		
		" run correction "
		
		" render corrected graph "
		for item in v_list2:
			nodeID = item[0]
			pose1 = item[1]
			node1 = self.nodeHash[nodeID]
			node1.setGPACPose(pose1)
			
		final_constraints = []
		
		print "TORO RESULTS:"
		for edge in e_list2:
			node1 = edge[0]
			node2 = edge[1]
			offset = edge[2]
			prec = edge[3]
			print node1, node2, offset, prec

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


			#self.merged_constraints[node1] = [transform, covE]
			

		return final_constraints


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
			

		f = open("sensor_constraints.txt", 'r')
		sensor_hypotheses = eval(f.read().rstrip())
		f.close()

		f = open("overlap_constraints.txt", 'r')
		overlap_hypotheses = eval(f.read().rstrip())
		f.close()

		for const in overlap_hypotheses:
			pylab.clf()
			self.renderSingleConstraint(const)
			
			pylab.title("Overlap Hypotheses %d ---> %d" % (const[0], const[1]))
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("constraint_%04u.png" % count)
			count += 1
					

		for const in sensor_hypotheses:
			pylab.clf()
			self.renderSingleConstraint(const, isSweep = True)
			pylab.title("Sensor Hypotheses %d ---> %d" % (const[0], const[1]))
			pylab.xlim(-4,4)
			pylab.ylim(-4,4)
			pylab.savefig("constraint_%04u.png" % count)
			count += 1
		
		
		#self.merged_constraints 
		

	def drawConstraints(self):
		
		motion_poses = [self.nodeHash[0].getGlobalGPACPose()]

		for i in range(len(self.motion_constraints)):
			transform = self.motion_constraints[i][2]
			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(motion_poses[-1])
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			motion_poses.append(currPose)	
		
		overlap_poses = [self.nodeHash[0].getGlobalGPACPose()]

		for i in range(len(self.overlap_constraints)):
			transform = self.overlap_constraints[i][2]
			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(overlap_poses[-1])
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			overlap_poses.append(currPose)	

		merged_poses = [self.nodeHash[0].getGlobalGPACPose()]

		for i in range(len(self.merged_constraints)):
			transform = self.merged_constraints[i][2]
			offset = [transform[0,0], transform[1,0], transform[2,0]]
			currProfile = Pose(merged_poses[-1])
			currPose = currProfile.convertLocalOffsetToGlobal(offset)	
			merged_poses.append(currPose)	


		pylab.clf()
		for i in range(self.numNodes):
			node1 = self.nodeHash[i]
			currPose = currPose = node1.getGndGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(0.5,0.5,0.5))	
			
			
			"""
			node1 = self.nodeHash[i]
			currPose = motion_poses[i]
			currProfile = Pose(currPose)
			posture1 = node1.getGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')		

			node1 = self.nodeHash[i]
			currPose = overlap_poses[i]
			currProfile = Pose(currPose)
			posture1 = node1.getGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')		

			node1 = self.nodeHash[i]
			currPose = merged_poses[i]
			currProfile = Pose(currPose)
			posture1 = node1.getGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='k')	
			
			"""
			
			node1 = self.nodeHash[i]
			currPose = currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='purple')	

			
			pylab.xlim(-8,12)
			pylab.ylim(-10,10)
			pylab.savefig("plotEstimate%04u.png" % i)
		
				
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
		#return
		
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
			localNode = self.nodeHash[i]
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

	def correctPosture(self):
		self.currNode.updateCorrectPosture()

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
		localNode = self.nodeHash[nodeNum]
		return localNode.getOccMap()

	def getNodePose(self, nodeNum):
		localNode = self.nodeHash[nodeNum]
		return localNode.getEstPose()
	
	def setNodePose(self, nodeNum, estPose):
		localNode = self.nodeHash[nodeNum]
		localNode.setEstPose(estPose)
	
	def obstCallBack(self, direction):
		self.currNode.obstCallBack(direction)
	
	def saveLocalMap(self):

		" save the local map "
		if self.currNode != 0:
			self.currNode.saveMap()
			
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

		
		
		
				
