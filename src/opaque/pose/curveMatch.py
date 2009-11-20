
#!/usr/bin/python

import math
from math import *
from random import *
from scipy.optimize import *
import scipy.interpolate
import numpy
import pylab
import getopt, sys, os, shutil, csv
from copy import *
#import toro
import traceback
from PoseProfile import *


class SplineFit:

	def __init__ (self, points, smooth = 0.1):
		self.pointSet = points
		self.smoothNess = smooth

		# unzip the points
		intArray = [[],[]]
		for p in points:
			intArray[0].append(p[0])
			intArray[1].append(p[1])

		#print "smooth = ", self.smoothNess, "input = ", intArray

		dataOutput = open("dataFile.txt", 'w')
		dataOutput.write(repr(intArray))
		dataOutput.close()

		#print len(intArray[0]), len(intArray[1])
	
		#print "performing spline fit"
		self.tck, self.u = scipy.interpolate.splprep(intArray, s = self.smoothNess, k=5)
		#print "fit complete"

		#print "result =", self.u
		#print self.tck
		#unew = scipy.arange(0, 1.01, 0.01)
		#out = scipy.interpolate.splev(unew,tck)
	
	def findClosestFromPointSet(self):
		pass

	def getPointOfDist(self, dist):
		# return a point of dist distance along the spline from the origin

		# if spline not long enough!
		totalLen = self.length()
		if dist > totalLen:
			raise

		if dist == totalLen:
			point = self.getU(1.0)
			return point

		# spline distance should not be negative
		if dist < 0:
			raise

		if dist == 0.0:
			point = self.getU(0.0)
			return point

		estU = 0.5
		jumpVal = 0.5

		count = 0

		while True:

			count += 1
			if count > 10:
				return []

			estDist = self.dist(0.0, estU)
			jumpVal /= 2.0

			#print estDist, dist

			if abs(dist-estDist) < 0.01:
				point = self.getU(estU)
				return point

			# increase the estimate
			if estDist < dist:
				estU += jumpVal

			# decrease the estimate
			if estDist > dist:
				estU -= jumpVal

	# compute the distance along the curve between the points parameterized by u1 and u2
	def dist(self, u1, u2, iter = 0.0001):
		totalDist = 0.0

		if u1 >= u2:
			raise


		unew = scipy.arange(u1, u2 + iter, iter)
		out = scipy.interpolate.splev(unew,self.tck)

		points = []
		for i in range(len(out[0])):
			points.append([out[0][i],out[1][i]])

		#print "measuring", u1, "to", u2+iter, "with", len(points), "results"
		#print points

		if len(unew) <= 1:
			return 0.0

		origin = points.pop(0)
		while len(points) > 0:
			pnt = points[0]
			#print pnt, origin
			dist = sqrt( (pnt[0]-origin[0])**2 + (pnt[1]-origin[1])**2 )
			totalDist += dist
			origin = points.pop(0)

		return totalDist

	# compute the tangent vector to u1
	def getUVector(self, u1, iter = 0.01):

		if u1 < 0.0 or u1 > 1.0:
			print "ERROR: u =", u
			raise

		if u1 > 1.0 - iter:
			u1 = u1-iter

		unew1 = [u1]
		newPoint1 = scipy.interpolate.splev(unew1,self.tck)
		unew2 = [u1 + iter]
		newPoint2 = scipy.interpolate.splev(unew2,self.tck)

		# tangent vector
		vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]

		# normalize
		mag = sqrt(vec[0]**2 + vec[1]**2)
		vec = [vec[0]/mag, vec[1]/mag]

		return vec

	# return the length of the curve
	def length(self, iter = 0.0001):
		totalDist = 0

		unew = scipy.arange(0, 1.0 + iter, iter)
		out = scipy.interpolate.splev(unew,self.tck)

		points = []
		for i in range(len(out[0])):
			points.append([out[0][i],out[1][i]])

		origin = points.pop(0)
		while len(points) > 0:
			pnt = points[0]
			dist = sqrt( (pnt[0]-origin[0])**2 + (pnt[1]-origin[1])**2 )
			totalDist += dist
			origin = points.pop(0)

		return totalDist

	def matchSpline(self, rearSpline, plot = False, index = 0):
		# this spline is fore, the new one is rear
		# the assumption here is represented by the function call
		# whether this is forward, or the other is forward

		# parametric boundaries of the spline sections we are comparing
		u0_F = 0.0
		u0_R = 0.0
		u1_F = 1.0
		u1_R = 1.0

		# since this is forward, u = 0.0 for fore, needs to match with some u on rear
		# get the this point for u=0.0
		self.originF = self.getU(0.0)

		# now find the closest point on the rear spline to self.originF
		dist, u0_R, self.originR = rearSpline.findClosestPoint(self.originF)

		# these are now our origins: self.originF, self.originR

		# now find the terminator of fore from rear
		self.termR = rearSpline.getU(1.0)
		dist, u1_F, self.termF = self.findClosestPoint(self.termR)

		# now we sample ten points in the fore spline boundaries
		diff = u1_F - u0_F
		inc = diff/10.0
		#print "u0_F =", u0_F, "u1_F =", u1_F
		#print "u0_R =", u0_R, "u1_R =", u1_R
	
		# if the end points are the same, lets default to a full length comparison
		if u1_F == u0_F:
			samples = scipy.arange(0.0,1.0,0.1)
		else:
			samples = scipy.arange(u0_F,u1_F,inc)

		sample_points = self.getUSet(samples)

		# for each sample point, compute the min distance to the rear spline
		distances = []
		closePoints = []
		for p in sample_points:
			dist, u_value, closest_p = rearSpline.findClosestPoint(p)
			distances.append(dist)
			closePoints.append(closest_p)

		cost = sum(distances)

		# now we sample ten points in the rear spline boundaries
		diff = u1_R - u0_R
		inc = diff/10.0

		if u1_R == u0_R:
			samples2 = scipy.arange(0.0,1.0,0.1)
		else:
			samples2 = scipy.arange(u0_R,u1_R,inc)

		sample_points2 = rearSpline.getUSet(samples2)

		# for each sample point, compute the min distance to the rear spline
		distances = []
		closePoints2 = []
		for p in sample_points2:
			dist, u_value, closest_p = self.findClosestPoint(p)
			distances.append(dist)
			closePoints2.append(closest_p)

		cost += sum(distances)
		
		if plot:
			#pylab.clf()
			#self.drawSpline()
			#rearSpline.drawSpline()

			pylab.scatter([self.originF[0], self.originR[0]], [self.originF[1], self.originR[1]])
			pylab.scatter([self.termF[0], self.termR[0]], [self.termF[1], self.termR[1]])

			for i in range(0,len(sample_points)):
				xP = [sample_points[i][0], closePoints[i][0]]
				yP = [sample_points[i][1], closePoints[i][1]]
				pylab.plot(xP,yP,color='0.5')
				pylab.scatter(xP,yP,color='0.5',linewidth='1')

			for i in range(0,len(sample_points2)):
				xP = [sample_points2[i][0], closePoints2[i][0]]
				yP = [sample_points2[i][1], closePoints2[i][1]]
				pylab.plot(xP,yP,color='0.5')
				pylab.scatter(xP,yP,color='0.5',linewidth='1')

			#pylab.xlim(3,7)
			#pylab.ylim(1,4)

			#pylab.title("Cost = %f" % cost)
			#pylab.savefig("pics/fit%04u.png" % index)

		return cost

	def drawSpline(self, clr = '0.5'):
		unew = scipy.arange(0, 1.01, 0.01)
		out = scipy.interpolate.splev(unew,self.tck)


		#pylab.plot(out[0],out[1], color = clr)	

		# flipped to align with Ogre x-z coord system
		#pylab.plot(out[1],out[0], color = clr)	
		pylab.plot(out[0],out[1], color = clr)	

	def getU(self, u):
		if u < 0.0 or u > 1.0:
			print "ERROR: u =", u
			raise

		unew = [u]
		newPoint = scipy.interpolate.splev(unew,self.tck)
		return newPoint

	def getUSet(self, u_set):
		for u in u_set:
			if u < 0.0 or u > 1.0:
				print "ERROR: u =", u
				raise

		newPoints = scipy.interpolate.splev(u_set,self.tck)
		# zip the points together
		zipPoints = []
		for i in range(0,len(newPoints[0])):
			zipPoints.append([newPoints[0][i],newPoints[1][i]])

		return zipPoints

	def findClosestPoint(self, pos):
		# find closest location on spline to pos
		# we do this in three nested searches,
		# sampling 10 uniformly distributed points at each level	

		sample = scipy.arange(0.0,1.01,0.01)
		points = self.getUSet(sample)
		min, min_i = self.findMinDistancePoint(points, pos)

		#print min, min_i

		if min_i == 0:
			# sample above index zero
			high_u = sample[1]
			low_u = sample[0]

		elif min_i == len(sample)-1:
			# sample below index zero 
			high_u = sample[9]
			low_u = sample[8]

		else:
			# sample above and below min index
			high_u = sample[min_i+1]
			low_u = sample[min_i-1]

		diff = high_u - low_u
		#print "high_u =", high_u
		inc = diff/10.0
		sample2 = scipy.arange(low_u,high_u,inc)

		sample3 = []
		for samp in sample2:
			sample3.append(samp)
		if high_u == 1.0:
			sample3 += [1.0]

		points2 = self.getUSet(sample3)
		min, min_i = self.findMinDistancePoint(points2, pos)

		# return distance, u value, and point
		return min, sample2[min_i], points2[min_i]

	def findMinDistancePoint(self, points, pos):

		# compute distances of each point to pos
		min = 10e9
		min_i = 0
		for i in range(0,len(points)):
			p = points[i]
			dist = sqrt((p[0]-pos[0])**2 + (p[1]-pos[1])**2)
			if dist < min:
				min_i = i
				min = dist

		return min, min_i

class GroundedSpline(SplineFit):

	def __init__ (self, points, gndPoints, smooth = 0.1):

		SplineFit.__init__(self, points, smooth)

		self.gndPoints = gndPoints

	def computeGroundPose(self, pnt):

		# 1. find closest nodes to point
		minDist = 1e100
		minNode = 0
		minIndex = 0
		mX = 0
		mY = 0

		for i in range(len(self.pointSet)):
			point = self.pointSet[i]

			dist = sqrt((point[0]-pnt[0])**2 + (point[1]-pnt[1])**2)	
			if dist < minDist:
				minDist = dist
				minNode = point[2]
				minIndex = i
				mX = point[0]
				mY = point[1]

		# 2. compute the average transfrom from est to gnd nodes
		gndPoint = self.gndPoints[minIndex]
		xDiff = gndPoint[0] - mX
		yDiff = gndPoint[1] - mY

		if gndPoint[2] != minNode:
			raise

		# 3. use average to compute gnd of point of u
		newPoint = [pnt[0]+xDiff, pnt[1]+yDiff]
	
		return newPoint


def printPoints(pointSet, clr = '0.5'):

	xP = []
	yP = []

	for p in pointSet:
		xP.append(p[0])
		yP.append(p[1])

	pylab.scatter(xP,yP, color= clr)

def printPoseRef(pose, offset = [0.0,0.0,0.0], clr = '0.5'):

	# get raw points for pose
	L_pose = deepcopy(pose.getLPose())
	R_pose = deepcopy(pose.getRPose())
	I_pose = deepcopy(pose.getIPose())

	# translate all points by (x,y)
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			p[0] += offset[0]
			p[1] += offset[1]

	# rotate all points with respect to p0 of L_pose by offset[2]
	origin = L_pose[0]
	ang = offset[2]
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			xdiff = p[0] - origin[0]
			ydiff = p[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			p[0] = origin[0] + xnew
			p[1] = origin[1] + ynew

	printPoints(L_pose, clr= clr)
	printPoints(R_pose, clr= clr)
	printPoints(I_pose, clr= clr)


def printPose(pose, offset = [0.0,0.0,0.0], clr = '0.5'):

	# get raw points for pose
	L_pose = deepcopy(pose.getLPose())
	R_pose = deepcopy(pose.getRPose())
	I_pose = deepcopy(pose.getIPose())

	# translate all points by (x,y)
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			p[0] += offset[0]
			p[1] += offset[1]

	# rotate all points with respect to p0 of L_pose by offset[2]
	origin = L_pose[0]
	ang = offset[2]
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			xdiff = p[0] - origin[0]
			ydiff = p[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			p[0] = origin[0] + xnew
			p[1] = origin[1] + ynew

	# compute splines for the resultant offsets of pose
	L_spline = SplineFit(L_pose)
	R_spline = SplineFit(R_pose)
	I_spline = SplineFit(I_pose)

	L_spline.drawSpline(clr)
	R_spline.drawSpline(clr)
	I_spline.drawSpline(clr)

def printCenterLine(pose, offset = [0.0,0.0,0.0], clr = '0.5'):

	# get raw points for pose
	L_pose = deepcopy(pose.getLPose())
	R_pose = deepcopy(pose.getRPose())
	I_pose = deepcopy(pose.getIPose())

	# translate all points by (x,y)
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			p[0] += offset[0]
			p[1] += offset[1]

	# rotate all points with respect to p0 of L_pose by offset[2]
	origin = L_pose[0]
	ang = offset[2]
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			xdiff = p[0] - origin[0]
			ydiff = p[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			p[0] = origin[0] + xnew
			p[1] = origin[1] + ynew

	# compute splines for the resultant offsets of pose
	I_spline = SplineFit(I_pose)

	I_spline.drawSpline(clr)

def printPoseWall(pose, offset = [0.0,0.0,0.0], clr = '0.5'):

	# get raw points for pose
	L_pose = deepcopy(pose.getLPose())
	R_pose = deepcopy(pose.getRPose())
	I_pose = deepcopy(pose.getIPose())

	# translate all points by (x,y)
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			p[0] += offset[0]
			p[1] += offset[1]

	# rotate all points with respect to p0 of L_pose by offset[2]
	origin = L_pose[0]
	ang = offset[2]
	for pose in [L_pose, R_pose, I_pose]:
		for p in pose:
			xdiff = p[0] - origin[0]
			ydiff = p[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			p[0] = origin[0] + xnew
			p[1] = origin[1] + ynew

	# compute splines for the resultant offsets of pose
	L_spline = SplineFit(L_pose)
	R_spline = SplineFit(R_pose)
	I_spline = SplineFit(I_pose)

	L_spline.drawSpline(clr)
	R_spline.drawSpline(clr)
	#I_spline.drawSpline(clr)

class CostObject:

	def __init__(self):
		self.index = 0
		self.offsets = []

	def getOffsets(self):
		return self.offsets

	def clearHistory(self):
		self.index = 0
		self.offsets = []
		return

	def cost_func(self, offset, pose1, pose2):

		# 1. pose1 is our fixed reference, so compute the splines for each arc
		# 2. get raw points of pose2, and offset them
		# 3. compute splines with resultant points
		# 4. run the matchSpline cost function on the each of the 3 pairs
		# 5. return the sum of the 3

		"""
		if self.index > 1:
			self.offsets.append(offset)
			return 0.0
		"""

		self.index += 1

		self.offsets.append(offset)

		L_spline1 = SplineFit(pose1.getLPose())
		R_spline1 = SplineFit(pose1.getRPose())
		I_spline1 = SplineFit(pose1.getIPose())

		# get raw points for pose 2
		L_pose2 = deepcopy(pose2.getLPose())
		R_pose2 = deepcopy(pose2.getRPose())
		I_pose2 = deepcopy(pose2.getIPose())

		# translate all points by (x,y)
		for pose in [L_pose2, R_pose2, I_pose2]:
			for p in pose:
				p[0] += offset[0]
				p[1] += offset[1]

		# rotate all points with respect to p0 of L_pose2 by offset[2]
		#origin = L_pose2[0]
		origin = pose1.getRootPose()
		
		ang = offset[2]
		for pose in [L_pose2, R_pose2, I_pose2]:
			for p in pose:
				xdiff = p[0] - origin[0]
				ydiff = p[1] - origin[1]

				xnew = xdiff*cos(ang) - ydiff*sin(ang)
				ynew = xdiff*sin(ang) + ydiff*cos(ang)

				p[0] = origin[0] + xnew
				p[1] = origin[1] + ynew

		# compute splines for the resultant offsets of pose2
		L_spline2 = SplineFit(L_pose2)
		R_spline2 = SplineFit(R_pose2)
		I_spline2 = SplineFit(I_pose2)

		# now compute the cost of each	
		#cost1 = L_spline1.matchSpline(L_spline2, True, self.index)
		#pylab.clf()

		cost1 = L_spline1.matchSpline(L_spline2, plot = False, index = self.index)
		cost2 = R_spline1.matchSpline(R_spline2, plot = False, index = self.index)
		cost3 = I_spline1.matchSpline(I_spline2, plot = False, index = self.index)
		#cost1 = L_spline1.matchSpline(L_spline2)
		#cost2 = R_spline1.matchSpline(R_spline2)
		#cost3 = I_spline1.matchSpline(I_spline2)

		cost = cost1 + cost2 + cost3

		#pylab.title("Total Cost = %f" % cost)
		#pylab.savefig("test_0003/pics/fit%04u.png" % self.index)


		"""
		pylab.clf()
		pylab.title("Cost = %f" % cost)
		L_spline1.drawSpline('0.5')
		R_spline1.drawSpline('0.5')
		I_spline1.drawSpline('0.5')

		L_spline2.drawSpline('b')
		R_spline2.drawSpline('b')
		I_spline2.drawSpline('b')
		
		pylab.xlim(-5,5)
		pylab.ylim(-8,2)
		pylab.savefig("test_0003/pics/fit%04u.png" % self.index)
		self.index += 1
		"""

		return cost
	
def readInflectionPoints(filename):

	pointReader = csv.reader(file(filename), delimiter=' ')

	poses = []

	for row in pointReader:
		inflectionSet = []

		# drop first element since it is the index
		row = row[1:]

		#for i in range(0,len(row)/2):
			#point = [float(row[2*i]),float(row[2*i+1])]

		for i in range(0,len(row)/3):
			point = [float(row[3*i+1]),float(row[3*i+2]), int(row[3*i])]

			if point[2] != -1:
				inflectionSet.append(point)

		if inflectionSet != []:
			# create pose object
			# pose = PoseProfile(polySet)

			poses.append(inflectionSet)

	return poses

def readSaddlePoints(filename):

	pointReader = csv.reader(file(filename), delimiter=' ')

	poses = []
	rightPoses = []
	leftPoses = []

	for row in pointReader:
		saddleSet = []

		# drop first element since it is the index
		row = row[1:]
		#print row

		for i in range(0,len(row)/3):
			# the 3rd element is the node ID
			point = [float(row[3*i+1]),float(row[3*i+2]), int(row[3*i])]
			#point = [float(row[2*i]),float(row[2*i+1])]
			saddleSet.append(point)

		if saddleSet != []:
			# create pose object
			# pose = PoseProfile(polySet)

			leftSide = []
			rightSide = []

			for i in range(0,len(saddleSet)):

				point = saddleSet[i]

				if point[2] != -1:
					if ( i % 2 == 0 ):
						leftSide.append(point)
					else:
						rightSide.append(point)

			leftPoses.append(leftSide)
			rightPoses.append(rightSide)
			#poses.append(saddleSet)

	#return poses
	return leftPoses, rightPoses

def writeToroGraph(testDir, v_list, e_list, edgeConstraints):

	filename = testDir + "/snake_pipe_fix.graph"
	sf = open(filename, 'w')

	for vertex in v_list:
		sf.write("VERTEX " + str(vertex[3]) + " " + str(vertex[0]) + " " + str(vertex[1]) + " " + str(vertex[2]))
		sf.write("\n")

	for edge in e_list:
		sf.write("EDGE ")
		for item in edge:
			sf.write(str(item) + " ")
		sf.write("\n")

	#constraint = [nodeID1, nodeID2, xDiff, yDiff, rDiff, nomVar, nomVar, nomVar, nomVar, nomVar, nomVar]

	for constraint in edgeConstraints:
		sf.write("EDGE ")
		for item in constraint:
			sf.write(str(item) + " ")
		sf.write("\n")

	sf.close()


def readToroGraph(filename):

	pointReader = csv.reader(file(filename), delimiter=' ')

	# read data into lists
	vertexList = []
	edgeList = []

	for row in pointReader:

		# add to the vertex list
		if row[0] == "VERTEX":
			# x, y, theta
			vertexList.append([float(row[2]),float(row[3]),float(row[4]),int(row[1])])	

		# add to the edge list
		if row[0] == "EDGE":
			# node 1, node 2
			edge = [int(row[1]),int(row[2]),float(row[3]),float(row[4]),float(row[5]),
				float(row[6]),float(row[7]),float(row[8]),float(row[9]),float(row[10]),float(row[11])]	
			edgeList.append(edge)	

	for i in range(0,len(vertexList)):
		if i != vertexList[i][3]:
			raise

	return vertexList, edgeList


def printHistory(testDir, poses, v_list, v_list_gnd, offsetHistories):

	pylab.clf()

	# find the max number of function evaluations
	max = 0
	index = 0
	for values in offsetHistories:
		if len(values) > max:
			max = len(values)

	# pad each offset value up until the max
	for values in offsetHistories:
		final_offset = values[-1]
		while len(values) < max:
			values.append(copy(final_offset))

	# add zero offset pose at the very beginning of the histories
	#for values in offsetHistories:
	#	values.insert(0, [0.0, 0.0, 0.0])

	# now we print the pose of every data at every iteration
	# so we can see how the map in its entirety is corrected


	# nominal displacements between poses
	nomDisp = []
	nomDisp.append(poses[0].getLPose()[0])

	for i in range(1, len(poses)):
		o1 = poses[i-1].getLPose()[0]
		o2 = poses[i].getLPose()[0]
		nomDisp.append([o2[0]-o1[0], o2[1]-o1[1]])

	# create local versions of each pose
	# 1. add the translation offset first between the poses
	# 2. rotate displacement vector by the sum of the angles
	# 3. add in local rotation

	for j in range(0,len(offsetHistories[0])):

		origin = [0.0,0.0,0.0]
		refPoint = [0.0,0.0,0.0]
		angleSum = 0.0

		localPoses = []
		for i in range(0,len(poses)):
			newPose = poses[i].copy()
			nomOffset = newPose.getLPose()[0]
			newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
			localPoses.append(newPose)
		
		nomOffset = nomDisp[0]
		localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
		refPoint = localPoses[0].getLPose()[0]

		for i in range(1,len(localPoses)):

			# nominal + local + reference
			nomOffset = nomDisp[i]
			localOffset = offsetHistories[i-1][j]

			# nom + local should be rotated by angleSum
			rotateOffset = [0.0,0.0,0.0]
			xdiff = localOffset[0] + nomOffset[0]
			ydiff = localOffset[1] + nomOffset[1]

			rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
			rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

			totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

			# perform the translation
			localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

			# new reference point
			refPoint = localPoses[i].getLPose()[0]

			# perform the rotation
			angleSum += localOffset[2]
			localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

		pylab.clf()

		xP = []
		yP = []
		for v in v_list:
			pose = v[1]
			xP.append(pose[0])
			yP.append(pose[1])
		pylab.scatter(xP,yP, color='0.9', faceted=False)

		# print the ground references in background
		xP = []
		yP = []
		for v in v_list_gnd:
			pose = v[1]
			xP.append(pose[0])
			yP.append(pose[1])
		pylab.scatter(xP,yP, color='0.8', faceted=False)


		for lpose in localPoses:
			#print "printing from origin", lpose.getLPose()[0]
			printPose(lpose)

		#pylab.xlim(-5,8)
		#pylab.ylim(-5,5)
		pylab.ylim(-5,20)
		pylab.xlim(-15,15)
		#pylab.xlim(-4,0)
		#pylab.ylim(-4,0)
		pylab.savefig(testDir + "/pics/fitHist%04u.png" % index)
		index += 1

	#print "nom displacements are:", nomDisp

def printBeforeAfter(testDir, poses, v_list, v_list_gnd, offsetHistories):

	pylab.clf()

	# before and after offsets
	initOffsets = []
	finalOffsets = []

	# pad each offset value up until the max
	for values in offsetHistories:
		initOffsets.append([0.0, 0.0, 0.0])
		finalOffsets.append(copy(values[-1]))

	# nominal displacements between poses
	nomDisp = []
	nomDisp.append(poses[0].getLPose()[0])

	for i in range(1, len(poses)):
		o1 = poses[i-1].getLPose()[0]
		o2 = poses[i].getLPose()[0]
		nomDisp.append([o2[0]-o1[0], o2[1]-o1[1]])

	# create local versions of each pose
	# 1. add the translation offset first between the poses
	# 2. rotate displacement vector by the sum of the angles
	# 3. add in local rotation

	#for j in range(0,len(offsetHistories[0])):
	origin = [0.0,0.0,0.0]
	refPoint = [0.0,0.0,0.0]
	angleSum = 0.0

	# use these parameters to choose the best scale for the plot
	xHigh = -10e100
	xLow = 10e100
	yHigh = -10e100
	yLow = 10e100

	localPoses = []
	for i in range(0,len(poses)):
		newPose = poses[i].copy()
		nomOffset = newPose.getLPose()[0]
		newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
		localPoses.append(newPose)

	nomOffset = nomDisp[0]
	localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
	refPoint = localPoses[0].getLPose()[0]

	for i in range(1,len(localPoses)):

		# nominal + local + reference
		nomOffset = nomDisp[i]
		localOffset = initOffsets[i-1]

		# nom + local should be rotated by angleSum
		rotateOffset = [0.0,0.0,0.0]
		xdiff = localOffset[0] + nomOffset[0]
		ydiff = localOffset[1] + nomOffset[1]

		rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
		rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

		totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

		# perform the translation
		localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

		# new reference point
		refPoint = localPoses[i].getLPose()[0]

		# perform the rotation
		angleSum += localOffset[2]
		localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

		xHighLoc, xLowLoc, yHighLoc, yLowLoc = localPoses[i].getMaxMin()

		if xHighLoc > xHigh:
			xHigh = xHighLoc
		if xLowLoc < xLow:
			xLow = xLowLoc
		if yHighLoc > yHigh:
			yHigh = yHighLoc
		if yLowLoc < yLow:
			yLow = yLowLoc

	pylab.clf()

	# print the estimated and ground references in background
	toro.plotToro(v_list, clr='0.9')
	toro.plotToro(v_list_gnd, clr='0.8')

	# find the min and maximum
	for v in v_list:
		pose = v[1]
		if pose[0] > xHigh:
			xHigh = pose[0]
		if pose[0] < xLow:
			xLow = pose[0]
		if pose[1] > yHigh:
			yHigh = pose[1]
		if pose[1] < yLow:
			yLow = pose[1]

	for v in v_list_gnd:
		pose = v[1]
		if pose[0] > xHigh:
			xHigh = pose[0]
		if pose[0] < xLow:
			xLow = pose[0]
		if pose[1] > yHigh:
			yHigh = pose[1]
		if pose[1] < yLow:
			yLow = pose[1]

	for lpose in localPoses:
		#print "printing from origin", lpose.getLPose()[0]
		printPose(lpose)

	pylab.ylim(-5,20)
	pylab.xlim(-15,15)
	pylab.savefig(testDir + "/poseFitBefore.png")

	origin = [0.0,0.0,0.0]
	refPoint = [0.0,0.0,0.0]
	angleSum = 0.0

	localPoses = []
	for i in range(0,len(poses)):
		newPose = poses[i].copy()
		nomOffset = newPose.getLPose()[0]
		newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
		localPoses.append(newPose)
	
	nomOffset = nomDisp[0]
	localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
	refPoint = localPoses[0].getLPose()[0]

	for i in range(1,len(localPoses)):

		# nominal + local + reference
		nomOffset = nomDisp[i]
		localOffset = finalOffsets[i-1]

		# nom + local should be rotated by angleSum
		rotateOffset = [0.0,0.0,0.0]
		xdiff = localOffset[0] + nomOffset[0]
		ydiff = localOffset[1] + nomOffset[1]

		rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
		rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

		totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

		# perform the translation
		localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

		# new reference point
		refPoint = localPoses[i].getLPose()[0]

		# perform the rotation
		angleSum += localOffset[2]
		localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

	pylab.clf()

	# print the estimated and ground references in background
	toro.plotToro(v_list, clr='0.9')
	toro.plotToro(v_list_gnd, clr='0.8')

	for lpose in localPoses:
		#print "printing from origin", lpose.getLPose()[0]
		printPose(lpose)

	pylab.ylim(-5,20)
	pylab.xlim(-15,15)
	pylab.savefig(testDir + "/poseFitAfter.png")

	# return a point for new ERROR computation
	return localPoses[-1].getLPose()[0]

def getCorrectPoses(poses, offsets, color = 'b'):

	# nominal displacements between poses
	nomDisp = []
	nomDisp.append(poses[0].getLPose()[0])

	for i in range(1, len(poses)):
		o1 = poses[i-1].getLPose()[0]
		o2 = poses[i].getLPose()[0]
		nomDisp.append([o2[0]-o1[0], o2[1]-o1[1]])

	# create local versions of each pose
	# 1. add the translation offset first between the poses
	# 2. rotate displacement vector by the sum of the angles
	# 3. add in local rotation

	origin = [0.0,0.0,0.0]
	refPoint = [0.0,0.0,0.0]
	angleSum = 0.0

	localPoses = []
	for i in range(0,len(poses)):
		newPose = poses[i].copy()
		nomOffset = newPose.getLPose()[0]
		newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
		localPoses.append(newPose)
	
	nomOffset = nomDisp[0]
	localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
	refPoint = localPoses[0].getLPose()[0]

	for i in range(1,len(localPoses)):

		# nominal + local + reference
		nomOffset = nomDisp[i]
		localOffset = offsets[i-1]

		# nom + local should be rotated by angleSum
		rotateOffset = [0.0,0.0,0.0]
		xdiff = localOffset[0] + nomOffset[0]
		ydiff = localOffset[1] + nomOffset[1]

		rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
		rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

		totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

		# perform the translation
		localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

		# new reference point
		refPoint = localPoses[i].getLPose()[0]

		# perform the rotation
		angleSum += localOffset[2]
		localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

	#for lpose in localPoses:
	#	printPoseWall(lpose, clr=color )
	#	lpose.printPosePoints('#e3b1b1')

	return localPoses

def printPoseFit(testDir, poses, offsets, color = 'b'):

	# nominal displacements between poses
	nomDisp = []
	nomDisp.append(poses[0].getLPose()[0])

	for i in range(1, len(poses)):
		o1 = poses[i-1].getLPose()[0]
		o2 = poses[i].getLPose()[0]
		nomDisp.append([o2[0]-o1[0], o2[1]-o1[1]])

	# create local versions of each pose
	# 1. add the translation offset first between the poses
	# 2. rotate displacement vector by the sum of the angles
	# 3. add in local rotation

	origin = [0.0,0.0,0.0]
	refPoint = [0.0,0.0,0.0]
	angleSum = 0.0

	localPoses = []
	for i in range(0,len(poses)):
		newPose = poses[i].copy()
		nomOffset = newPose.getLPose()[0]
		newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
		localPoses.append(newPose)
	
	nomOffset = nomDisp[0]
	localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
	refPoint = localPoses[0].getLPose()[0]

	for i in range(1,len(localPoses)):

		# nominal + local + reference
		nomOffset = nomDisp[i]
		localOffset = offsets[i-1]

		# nom + local should be rotated by angleSum
		rotateOffset = [0.0,0.0,0.0]
		xdiff = localOffset[0] + nomOffset[0]
		ydiff = localOffset[1] + nomOffset[1]

		rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
		rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

		totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

		# perform the translation
		localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

		# new reference point
		refPoint = localPoses[i].getLPose()[0]

		# perform the rotation
		angleSum += localOffset[2]
		localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

	for lpose in localPoses:
		printPoseWall(lpose, clr=color )
		lpose.printPosePoints('#e3b1b1')

	return localPoses

def printPath(testDir, poses, offsets, color = 'b'):

	# nominal displacements between poses
	nomDisp = []
	nomDisp.append(poses[0].getLPose()[0])

	for i in range(1, len(poses)):
		o1 = poses[i-1].getLPose()[0]
		o2 = poses[i].getLPose()[0]
		nomDisp.append([o2[0]-o1[0], o2[1]-o1[1]])

	# create local versions of each pose
	# 1. add the translation offset first between the poses
	# 2. rotate displacement vector by the sum of the angles
	# 3. add in local rotation

	origin = [0.0,0.0,0.0]
	refPoint = [0.0,0.0,0.0]
	angleSum = 0.0

	localPoses = []
	for i in range(0,len(poses)):
		newPose = poses[i].copy()
		nomOffset = newPose.getLPose()[0]
		newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
		localPoses.append(newPose)
	
	nomOffset = nomDisp[0]
	localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
	refPoint = localPoses[0].getLPose()[0]

	for i in range(1,len(localPoses)):

		# nominal + local + reference
		nomOffset = nomDisp[i]
		localOffset = offsets[i-1]

		# nom + local should be rotated by angleSum
		rotateOffset = [0.0,0.0,0.0]
		xdiff = localOffset[0] + nomOffset[0]
		ydiff = localOffset[1] + nomOffset[1]

		rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
		rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

		totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

		# perform the translation
		localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

		# new reference point
		refPoint = localPoses[i].getLPose()[0]

		# perform the rotation
		angleSum += localOffset[2]
		localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

	for lpose in localPoses:
		printCenterLine(lpose, clr=color )
		#lpose.printPosePoints('#e3b1b1')

	return localPoses

def computeTermPoint(testDir, poses, v_list, v_list_gnd, offsetHistories):

	# before and after offsets
	initOffsets = []
	finalOffsets = []

	# pad each offset value up until the max
	for values in offsetHistories:
		initOffsets.append([0.0, 0.0, 0.0])
		finalOffsets.append(copy(values[-1]))

	# nominal displacements between poses
	nomDisp = []
	nomDisp.append(poses[0].getLPose()[0])

	for i in range(1, len(poses)):
		o1 = poses[i-1].getLPose()[0]
		o2 = poses[i].getLPose()[0]
		nomDisp.append([o2[0]-o1[0], o2[1]-o1[1]])

	# create local versions of each pose
	# 1. add the translation offset first between the poses
	# 2. rotate displacement vector by the sum of the angles
	# 3. add in local rotation

	# use these parameters to choose the best scale for the plot
	xHigh = -10e100
	xLow = 10e100
	yHigh = -10e100
	yLow = 10e100

	origin = [0.0,0.0,0.0]
	refPoint = [0.0,0.0,0.0]
	angleSum = 0.0

	localPoses = []
	for i in range(0,len(poses)):
		newPose = poses[i].copy()
		nomOffset = newPose.getLPose()[0]
		newPose.performOffset([-nomOffset[0],-nomOffset[1],0.0], origin)
		localPoses.append(newPose)
	
	nomOffset = nomDisp[0]
	localPoses[0].performOffset([nomOffset[0],nomOffset[1],0.0], origin)
	refPoint = localPoses[0].getLPose()[0]

	for i in range(1,len(localPoses)):

		# nominal + local + reference
		nomOffset = nomDisp[i]
		localOffset = finalOffsets[i-1]

		# nom + local should be rotated by angleSum
		rotateOffset = [0.0,0.0,0.0]
		xdiff = localOffset[0] + nomOffset[0]
		ydiff = localOffset[1] + nomOffset[1]

		rotateOffset[0] = xdiff*cos(angleSum) - ydiff*sin(angleSum)
		rotateOffset[1] = xdiff*sin(angleSum) + ydiff*cos(angleSum)

		totalOffset = [refPoint[0] + rotateOffset[0], refPoint[1] + rotateOffset[1]]

		# perform the translation
		localPoses[i].performOffset([totalOffset[0],totalOffset[1],0.0], origin)

		# new reference point
		refPoint = localPoses[i].getLPose()[0]

		# perform the rotation
		angleSum += localOffset[2]
		localPoses[i].performOffset([0.0, 0.0, angleSum], refPoint)

	# return a point for new ERROR computation
	return localPoses[-1].getLPose()[0]

def createPoses(testDir, L_saddlePoses, R_saddlePoses, inflectionPoses):
	
	# create a spline for each type of curve
	poses = []

	try:
		for pose in R_saddlePoses:
			newSpline = SplineFit(pose)
		for pose in L_saddlePoses:
			newSpline = SplineFit(pose)
		for pose in inflectionPoses:
			newSpline = SplineFit(pose)
	except:
		print "failed to fit spline"
		traceback.print_exc()
		raise

	for i in range(0, len(L_saddlePoses)):
		L_pose = L_saddlePoses[i]
		R_pose = R_saddlePoses[i]
		I_pose = inflectionPoses[i]
		pose = PoseProfile(L_pose, R_pose, I_pose)
		poses.append(pose)

	return poses

def readOffsets(testDir):

	fileName = testDir + "/offsets.txt"
	sf = open(fileName, 'r')

	offsetReader = csv.reader(file(fileName), delimiter=' ')
	offsetHistories = []

	for row in offsetReader:
		thisRelPose = []
		for i in range(0,len(row)/3):
			offset = [float(row[3*i]), float(row[3*i+1]), float(row[3*i+2])]
			thisRelPose.append(offset)
				
		offsetHistories.append(thisRelPose)

	return offsetHistories

def writeOffsets(testDir, offsetHistories):

	fileName = testDir + "/offsets.txt"
	sf = open(fileName, 'w')

	for relPose in offsetHistories:
		for offset in relPose:

			sf.write(str(offset[0]) + " " + str(offset[1]) + " " + str(offset[2]) + " ")

		sf.write("\n")

	sf.close()
	

def contourErrorCorrect(L_saddlePoses, R_saddlePoses, inflectionPoses, testDir = ".", v_list = [], e_list = []):

	# create a spline for each type of curve
	poses = []

	try:
		for pose in R_saddlePoses:
			newSpline = SplineFit(pose)
		for pose in L_saddlePoses:
			newSpline = SplineFit(pose)
		for pose in inflectionPoses:
			newSpline = SplineFit(pose)
	except:
		print "failed to fit spline"
		traceback.print_exc()
		raise

	for i in range(0, len(L_saddlePoses)):
		L_pose = L_saddlePoses[i]
		R_pose = R_saddlePoses[i]
		I_pose = inflectionPoses[i]
		pose = PoseProfile(L_pose, R_pose, I_pose)
		poses.append(pose)


	# make this initial guess
	# x, y, p
	initGuess = [0.0, 0.0, 0.0]
	offsets = []
	costs = []

	costObject = CostObject()

	offsetHistories = []

	rd = Random()
	
	for i in range(0,len(poses)-1):

		# now let's measure the cost of matching two splines
		# for each pair of splines, search for the lowest cost
		fmin_result = fmin(costObject.cost_func,initGuess,args=(poses[i], poses[i+1]))
		#fmin_result = fmin_powell(costObject.cost_func,initGuess,args=(poses[i], poses[i+1]))

		offsetHistories.append(costObject.getOffsets())
		costObject.clearHistory()

		cost = costObject.cost_func(fmin_result, poses[i],poses[i+1])

		costs.append(cost)
		offsets.append(fmin_result)
	
	# create the constraints between the contour matched poses
	#edgeConstraints = []
	#for i in range(0, len(poses)-1):
	#	constraints = poses[i].createConstraints(poses[i+1], v_list)
	#	edgeConstraints += constraints

	#writeToroGraph(testDir, v_list, e_list, edgeConstraints)
	#filename = testDir + "/snake_pipe_fix.graph"
	#toro.writeToroGraph(testDir, "snake_pipe_fix.graph", v_list, e_list, edgeConstraints)

	return poses, offsetHistories

	#printHistory(testDir, poses, v_list, v_list_gnd, offsetHistories)

	# print before and after
	#lastPoint = printBeforeAfter(testDir, poses, v_list, v_list_gnd, offsetHistories)





	## RUN OPTIMIZATION ALGORITHMS
	#fmin_result = fmin(cost_func,initGuess,args=(rect1, rect2))
	#powell_result = fmin_powell(cost_func,initGuess,args=(rect1, rect2))
	#cg_result = fmin_cg(cost_func,initGuess,args=(rect1, rect2))
	#bfgs_result = fmin_bfgs(cost_func,initGuess,args=(rect1, rect2))
	# we are offseting polygon 1

if __name__ == '__main__':

	contourErrorCorrect("test_0002")

