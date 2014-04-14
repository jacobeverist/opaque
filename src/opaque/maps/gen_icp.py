#!/usr/bin/python

"""
2D Generalized ICP

WHAT:

This is a python implementation of the generalized ICP algorithm described in:
"Generalized-ICP", Aleksandr V. Segal, Dirk Haehnel, and Sebastian Thrun, In Robotics: Science and Systems 2009

Paper available on Aleksander Segal's personal website:
http://www.stanford.edu/~avsegal/

TO INSTALL:

This code requires 3 modules to be installed on your system.  Numpy, Scipy, and the Python NIPALS PCA module by Henning Risvik.

http://www.scipy.org/kk

http://numpy.scipy.org/

http://folk.uio.no/henninri/pca_module/


TO RUN:

To run the example program, simply execute this file from the command-line with the following command:
python gen_icp.py

Copy the main section to your own code to have a working implementation and include this file as a module.

This code is provided under the GNU General Public License, Version 3 

Copyright 2009, by Jacob Everist
jacob.everist@gmail.com

http://jacobeverist.com/gen_icp

"""

import sys
import numpy
import scipy 
import scipy.linalg
import scipy.optimize
import math
import pylab
import pca_module
import functions
from functions import logFunction
from subprocess import *
import traceback
from math import cos, sin, pi
import cProfile

from matplotlib.patches import Circle
import matplotlib.pyplot as plt

from copy import copy
from copy import deepcopy

from Pose import Pose
from SplineFit import SplineFit

import Image

from icp import shapeCostC
from icp import computeMatchErrorP
from icp import matchPairs

import time

import multiprocessing as processing
import ctypes, os

numIterations = 0
globalPlotCount = 0

overlapPool = None
overlapIn = None
overlapOut = None
globalPool = None

import nelminICP


from scipy.spatial import cKDTree
from numpy import array


def __num_processors():
	
	return 4
	
	if os.name == 'nt': # Windows
		return int(os.getenv('NUMBER_OF_PROCESSORS'))
	else: # glibc (Linux, *BSD, Apple)
		get_nprocs = ctypes.cdll.libc.get_nprocs
		get_nprocs.restype = ctypes.c_int
		get_nprocs.argtypes = []
		return get_nprocs()

def __remote_ICP(rank, qin, qout):

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		#print rank, "received", nc, args
		# process data
		#knn = __do_nothing(data, nc, someArg2, someArg3)
		results = []
		for arg in args:
			u1 = arg[1]
			u2 = arg[2]
			angGuess = arg[3]
			junctionPose1 = arg[4]
			orientedPathSoup = arg[5]
			globalPath2 = arg[6]
			plotIter = arg[7]
			pathID1 = arg[8]
			parentID = arg[9]

			resultPose, lastCost, matchCount = branchEstimateICP([u1,u2,angGuess], junctionPose1, orientedPathSoup, globalPath2, plotIter = plotIter, n1 = pathID1, n2 = parentID)
			
			results.append([resultPose, lastCost, matchCount])
		# write to output queue
		qout.put((nc,results))


def batchGlobalICP(args):

 

	ndata = len(args)
	nproc = __num_processors()
	
	print "nproc =", nproc
	
	# compute chunk size
	chunk_size = ndata / nproc
	chunk_size = 2 if chunk_size < 2 else chunk_size
	
	print "chunk_size =", chunk_size
	print "max_size =", ndata/chunk_size
	
	# set up a pool of processes
	qin = processing.Queue(maxsize=ndata/chunk_size)
	qout = processing.Queue(maxsize=ndata/chunk_size)
	pool = [processing.Process(target=__remote_ICP,
				args=(rank, qin, qout))
					for rank in range(nproc)]
	for p in pool: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	while 1:
		_data = args[cur:cur+chunk_size]
		if len(_data) == 0: break
		print "put(", (nc,_data), ")"
		qin.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout.get()]
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp
	# terminate workers
	for p in pool: p.terminate()
	return knn


def computeOffset(point1, point2, ang1, ang2):

	" corner points and orientations "
	overlapPoint1Pose = [point1[0], point1[1], ang1]
	overlapPoint2Pose = [point2[0], point2[1], ang2]
	
	" convert the desired intersection point on curve 1 into global coordinates "
	poseProfile1 = Pose([0.0,0.0,0.0])
	
	" now convert this point into a pose, and perform the inverse transform using corner2Pose "
	desGlobalPose2 = Pose(overlapPoint1Pose)
	
	" perform inverse offset from the destination pose "
	negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
	
	" relative pose between pose 1 and pose 2 to make corners coincide and same angle "
	resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
	localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
	
	return [localOffset[0], localOffset[1], localOffset[2]]


def computeBBox(points):
	
	xMax = -1e100
	xMin = 1e100
	yMax = -1e100
	yMin = 1e100
	
	
	for p in points:
		
		if p[0] > xMax:
			xMax = p[0]
		
		if p[0] < xMin:
			xMin = p[0]
			
		if p[1] > yMax:
			yMax = p[1]
		
		if p[1] < yMin:
			yMin = p[1]


	" pad the edges so not points on are on the boundary "
	xMax = xMax + 0.4
	xMin = xMin - 0.4
	yMax = yMax + 0.4
	yMin = yMin - 0.4
	
	bbox = [[xMin,yMin], [xMin,yMax], [xMax,yMax], [xMax,yMin]]
			
	return bbox

" displace the point by the offset plus modify it's covariance "
def dispPointWithAngle(p, offset):
	
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]
	
	px = p[0]
	py = p[1]
	
	tx = px*math.cos(theta) - py*math.sin(theta) + xd
	ty = px*math.sin(theta) + py*math.cos(theta) + yd
	
	p_off = [tx, ty]

	Cv = p[2]
	
	r11 = math.cos(theta)
	r12 = -math.sin(theta)
	r21 = math.sin(theta)
	r22 = r11

	c11 = Cv[0][0]
	c12 = Cv[0][1]
	c21 = Cv[1][0]
	c22 = Cv[1][1]
	
	res11 = r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
	res12 = r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
	res21 = r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
	res22 = r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)

	Ca = [[res11, res12], [res21, res22]]
	p_off.append(Ca)

	pa = p[3]
	p_off.append(functions.normalizeAngle(pa + theta))
	
	return p_off

" displace the point by the offset plus modify it's covariance "
def dispPoint(p, offset):
	
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	px = p[0]
	py = p[1]
	
	tx = px*math.cos(theta) - py*math.sin(theta) + xd
	ty = px*math.sin(theta) + py*math.cos(theta) + yd
	
	p_off = [tx, ty]

	Cv = p[2]
	
	r11 = math.cos(theta)
	r12 = -math.sin(theta)
	r21 = math.sin(theta)
	r22 = r11

	c11 = Cv[0][0]
	c12 = Cv[0][1]
	c21 = Cv[1][0]
	c22 = Cv[1][1]
	
	res11 = r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
	res12 = r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
	res21 = r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
	res22 = r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)

	Ca = [[res11, res12], [res21, res22]]
	p_off.append(Ca)
	
	return p_off
	
" displace the point by the offset only.  No covariance "
def dispOffset(p, offset):
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	px = p[0]
	py = p[1]
	
	p_off = [px*math.cos(theta) - py*math.sin(theta) + xd, px*math.sin(theta) + py*math.cos(theta) + yd]
	
	return p_off

def disp(ai, bi, T):

	temp = T*ai
	result = bi-temp	

	result[2] = 1.0

	return result

def computeMatchError(offset, a, b, Ca, Cb):


	xd = offset[0]
	yd = offset[1]
	theta = offset[2]
	
	ax = a[0]
	ay = a[1]
	
	bx = b[0]
	by = b[1]

	tx = ax*math.cos(theta) - ay*math.sin(theta) + xd
	ty = ax*math.sin(theta) + ay*math.cos(theta) + yd
	dx = bx - tx
	dy = by - ty

	r11 = math.cos(theta)
	r12 = -math.sin(theta)
	r21 = math.sin(theta)
	r22 = r11

	c11 = Ca[0][0]
	c12 = Ca[0][1]
	c21 = Ca[1][0]
	c22 = Ca[1][1]
	
	b11 = Cb[0][0]
	b12 = Cb[0][1]
	b21 = Cb[1][0]
	b22 = Cb[1][1]
	
	res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
	res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
	res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
	res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)
	
	resDet = res22*res11 - res12*res21
	
	q11 = res22/resDet
	q12 = -res12/resDet
	q21 = -res21/resDet
	q22 = res11/resDet

	errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22)

	return errVal


def computeMatchErrorSimple(offset, a, b, Ca, Cb):

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	  [math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
			])

	R = numpy.matrix([	  [math.cos(theta), -math.sin(theta)],
		[math.sin(theta), math.cos(theta)] ])

	r00 = math.cos(theta)
	r01 = -math.sin(theta)
	r10 = math.sin(theta)
	r11 = math.cos(theta)

	ai = numpy.concatenate((a,numpy.matrix([1.0])))
	bi = numpy.concatenate((b,numpy.matrix([1.0])))
	temp = T*ai
	d_vec = bi-temp
	d_vec[2] = 1.0	  

	res = Cb + R * Cv * numpy.transpose(R)

	# remove the 3rd homogeneous dimension for inverse
	invMat = scipy.linalg.inv(res)

	# add homogeneous dimension back
	invMat = numpy.concatenate((invMat,numpy.matrix([[0.0],[0.0]])), 1)
	invMat = numpy.concatenate((invMat,numpy.matrix([0.0,0.0,0.0])))

	error = numpy.transpose(d_vec)*invMat*d_vec
	errVal = error[0,0]

	return errVal


def findLocalNormal(pnt,points):

	x_list = []
	y_list = []

	pnt_count = 0
	pnts_copy = deepcopy(points)

	while pnt_count < 10:
		# select 3 closest points
		minDist = 1e100
		minPoint = []
		for p in pnts_copy:
			dist = math.sqrt((p[0]-pnt[0])**2 + (p[1]-pnt[1])**2)
		
			if dist < minDist:
				minDist = dist
				minPoint = p

		x_list.append(minPoint[0])
		y_list.append(minPoint[1])

		pnts_copy.remove(minPoint)
		pnt_count += 1

	x_list.append(pnt[0])
	y_list.append(pnt[1])

	cov_a = scipy.cov(x_list,y_list)

	loadings = []

	# NOTE:  seems to create opposing colinear vectors if data is colinear, not orthogonal vectors

	try:
		scores, loadings, E = pca_module.nipals_mat(cov_a, 2, 0.000001, False)

	except:
		raise

	if len(loadings) < 2:
		raise

	# return the second vector returned from PCA because this has the least variance (orthogonal to plane)
	return loadings[1]

def computeVectorCovariance(vec,x_var,y_var):

	c11 = x_var
	c12 = c21 = 0.0
	c22 = y_var

	mag = math.sqrt(vec[0]**2 + vec[1]**2)
	normVec = [vec[0]/mag, vec[1]/mag]

	if normVec[1] == 0:
		r22 = r11 = 1.0
		r12 = r21 = 0.0

	else:
		B = -1 / (normVec[1] + normVec[0]**2/normVec[1])
		A = -normVec[0]*B/normVec[1]

		r22 = r11 = A
		r12 = -B
		r21 = B

	res11 = r11*(c11*r11 + c12*r21) + r21*(c21*r11 + c22*r21)
	res12 = r11*(c11*r12 + c12*r22) + r21*(c21*r12 + c22*r22)
	res21 = r12*(c11*r11 + c12*r21) + r22*(c21*r11 + c22*r21)
	res22 = r12*(c11*r12 + c12*r22) + r22*(c21*r12 + c22*r22)
	
	Ca = [[res11, res12], [res21, res22]]
	
	return Ca


def findClosestPointInA(a_trans, b):

	minDist = 1e100
	minPoint = None
	min_i = 0

	for i in range(len(a_trans)):
		p = a_trans[i]

		dist = math.sqrt((p[0]-b[0])**2 + (p[1]-b[1])**2)
		if dist < minDist:
			minPoint = copy(p)
			minDist = dist
			min_i = i

	if minPoint != None:
		return minPoint, min_i, minDist
	else:
		raise

def findClosestPointWithAngle(points1, pnt, angleThresh):

	#print "findClosestPointWithAngle:", len(points1), pnt, angleThresh

	ax = pnt[0]
	ay = pnt[1]    
	a_theta = pnt[2]
		
	minDist = 1e100
	minPoint = None
	minI = 0

	for i in range(len(points1)):

		p = points1[i]
		
		angDiff = abs(functions.normalizeAngle(a_theta-p[2]))

		if angleThresh > angDiff or (math.pi - angleThresh) < angDiff:
	
			dist = math.sqrt((p[0]-ax)**2 + (p[1]-ay)**2)
			if dist < minDist:
				minPoint = copy(p)
				minDist = dist
				minI = i

	if minPoint != None:
		return minPoint, minI, minDist
	else:
		raise


# for point T*a, find the closest point b in B
def findClosestPointInB(b_data, a, offset):

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	ax = a[0]
	ay = a[1]	 
	a_off = [ax*math.cos(theta) - ay*math.sin(theta) + xd, ax*math.sin(theta) + ay*math.cos(theta) + yd]
	
	minDist = 1e100
	minPoint = None

	for p in b_data:

		dist = math.sqrt((p[0]-a_off[0])**2 + (p[1]-a_off[1])**2)
		if dist < minDist:
			minPoint = copy(p)
			minDist = dist


	if minPoint != None:
		return minPoint, minDist
	else:
		raise


def medialOverlapCostFunc(params, match_pairs, poses_1, poses_2, uHigh, uLow, u1):
	global numIterations

	currU = params[0]
	currAng = params[1]

	if currU < 0.0:
		return 1e5 * (0.15-currU)
	if currU > 1.0:
		return 1e5 * (currU - 0.85)

	if currU < uLow:
		return 1e5 * (uLow+0.05-currU)
	if currU > uHigh:
		return 1e5 * (currU - uHigh+0.05)

	" compute the point and angles from the parameters "
	if u1 >= 1.0:
		pose1 = poses_1[-1]
	elif u1 < 0.0:
		pose1 = poses_1[0]
	else:
		pose1 = poses_1[int(u1*100)]

	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]
		
	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]

	offset = computeOffset(point1, point2, ang1, ang2 + currAng)

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
	
		val = computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
		
		vals.append(val)
		sum1 += val
		
	return sum1


def overlapHistogram(params, match_pairs, medialSpline1, medialSpline2, u1):

	global numIterations

	currU = params[0]
	currAng = params[1]
	
	" compute the point and angles from the paramters "
	if u1+0.02 > 1.0:
		pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
	elif u1 < 0.0:
		pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
	else:
		pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

	if currU+0.02 > 1.0:
		pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
	elif currU < 0.0:
		pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
	else:
		pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]

	offset = computeOffset(point1, point2, ang1, ang2 + currAng)

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

		xd = offset[0]
		yd = offset[1]
		theta = offset[2]
	
		tx = ax*math.cos(theta) - ay*math.sin(theta) + xd
		ty = ax*math.sin(theta) + ay*math.cos(theta) + yd
		dx = bx - tx
		dy = by - ty
		
		val = math.sqrt(dx*dx + dy*dy)
	
		vals.append(val)
	
	lowCount = 0
	midCount = 0
	highCount = 0
	for val in vals:
		if val < 0.2:
			lowCount += 1
		elif val >= 0.2 and val <= 0.7:
			midCount += 1
		else:
			highCount += 1
	
	return lowCount, midCount, highCount

def addDistanceFromOriginCovariance(points, tan_var=0.1, perp_var=0.01):

	for p in points:

		# the covariance matrix that represents increasing uncertainty as distance from origin increases
		pVec = p
		dist = math.sqrt(pVec[0]**2 + pVec[1]**2)
		
		# 1. get the orthogonal vector, rotation 90 degrees
		rotVec = [0.,0.]
		rotVec[0] = -pVec[1]
		rotVec[1] = pVec[0]

		# 2. normalize vector	 
		rotVec[0] /= dist
		rotVec[1] /= dist

		# 3. computeVectorCovariance() with tan_var, perp_var proportional to dist, and tan_var larger
		C = computeVectorCovariance(rotVec,dist*tan_var,dist*perp_var)

		if len(p) <= 2:
			p.append(C)
		else:
			p[2][0][0] += C[0][0]
			p[2][0][1] += C[0][1]
			p[2][1][0] += C[1][0]
			p[2][1][1] += C[1][1]
			

def addGPACVectorCovariance(points, high_var=1.0, low_var=0.001):

	newPoints = []

	for p in points:
		C = [ [high_var, 0.], [0., high_var]]

		# the covariance matrix that enforces the point-to-plane constraint
		normVec = [math.cos(p[2]+math.pi/2.0), math.sin(p[2]+math.pi/2.0)]		  
		C = computeVectorCovariance(normVec,low_var,high_var)


		newPoints.append([p[0], p[1], C])

	return newPoints

def addGPACVectorCovarianceWithAngle(points, high_var=1.0, low_var=0.001):

	newPoints = []

	for p in points:
		C = [ [high_var, 0.], [0., high_var]]

		# the covariance matrix that enforces the point-to-plane constraint
		normVec = [math.cos(p[2]+math.pi/2.0), math.sin(p[2]+math.pi/2.0)]		  
		C = computeVectorCovariance(normVec,low_var,high_var)


		newPoints.append([p[0], p[1], C, p[2]])

	return newPoints

def addPointToLineCovariance(points, high_var=1.0, low_var=0.001):

	for p in points:
		C = [ [high_var, 0.], [0., high_var]]

		try:
			# the covariance matrix that enforces the point-to-plane constraint
			normVec = findLocalNormal(p,points)
			C = computeVectorCovariance(normVec,low_var,high_var)

		except:
			pass

		if len(p) <= 2:
			p.append(C)
		else:
			p[2][0][0] += C[0][0]
			p[2][0][1] += C[0][1]
			p[2][1][0] += C[1][0]
			p[2][1][1] += C[1][1]
			
" check if a point is contained in a circle "
def isInCircle(p, radius, center):
	
	dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
	if dist < radius:
		return True

	return False

" check if a point is contained in the polygon, use radius and center to quicken check "
def isValid(p, radius, center, poly):
	
	dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
	if dist < radius:
		if not functions.point_inside_polygon(p[0],p[1],poly):
			return True

	return False
	
	
def isValidPast(p, pastCircles, poly):
	
	for circle in pastCircles:
		
		radius, center = circle
		dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
	
		if dist < radius:
			if not functions.point_inside_polygon(p[0],p[1],poly):
				return True

	return False

def filterVertices(offset, radius, a_trans, b_data):

	circleA = [[offset[0],offset[1]], radius]
	circleB = [[0.0,0.0], radius]

	filteredAPoints = []
	filteredBPoints = []
	refilteredAPoints = []
	refilteredBPoints = []
	
	for p in b_data:
		dist = math.sqrt((p[0] - offset[0]) ** 2 + (p[1] - offset[1]) ** 2)
		
		if dist < radius:
			filteredBPoints.append(copy(p))

	for p in a_trans:
		dist = math.sqrt((p[0] - 0.0) ** 2 + (p[1] - 0.0) ** 2)
		
		if dist < radius:
			filteredAPoints.append(copy(p))

	for p in filteredBPoints:
		if not functions.point_inside_polygon(p[0],p[1],a_trans):
			refilteredBPoints.append(copy(p))

	for p in filteredAPoints:
		if not functions.point_inside_polygon(p[0],p[1],b_data):
			refilteredAPoints.append(copy(p))

	return refilteredAPoints, refilteredBPoints

def computeEnclosingCircle(a_data):
		
	maxA = 0.0
	a_p1 = []
	a_p2 = []
	for i in range(len(a_data)):
		p1 = a_data[i]

		for j in range(i+1, len(a_data)):
			p2 = a_data[j]

			dist = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

			if dist > maxA:
				maxA = dist
				a_p1 = copy(p1)
				a_p2 = copy(p2)

	radiusA = maxA/2.0

	centerA = [(a_p1[0] + a_p2[0])/2.0, (a_p1[1] + a_p2[1])/2.0]
	
	return radiusA, centerA

	
@logFunction
def overlapICP(estPose1, initGuess, medialPoints1, medialPoints2, rootPose1, rootPose2, inPlace = False, plotIter = False, n1 = 0, n2 = 0, uRange = 1.5):

	global numIterations

	u1 = initGuess[0]
	u2 = initGuess[1]
	
	currU = u2
	currAng = initGuess[2]
	medialSpline1 = SplineFit(medialPoints1, smooth=0.1)
	medialSpline2 = SplineFit(medialPoints2, smooth=0.1)

	uSet = [i*0.01 for i in range(100)]
	poses_1 = medialSpline1.getUVecSet(uSet)
	poses_2 = medialSpline2.getUVecSet(uSet)



	if u1 >= 1.0:
		pose1 = poses_1[-1]
	elif u1 < 0.0:
		pose1 = poses_1[0]
	else:
		pose1 = poses_1[int(u1*100)]

	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]

	
	uHigh = medialSpline2.getUOfDist(u2, uRange, distIter = 0.001)
	uLow = medialSpline2.getUOfDist(u2, -uRange, distIter = 0.001)

	if inPlace:
		uHigh = u2 + 0.08
		uLow = u2 - 0.08


	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]

	offset = computeOffset(point1, point2, ang1, ang2 + currAng)


	costThresh = 0.004
	minMatchDist = 2.0
	lastCost = 1e100
	
	startIteration = numIterations

	" set the initial guess "
	poseOrigin = Pose(estPose1)
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
	points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)

	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	points2_offset = []
	for p in points2:
		result = dispPoint(p, offset)		 
		points2_offset.append(result)


	#kdInput1 = []
	#for i in range(0,len(points1)):
	#	p_1 = points1[i]
	#	kdInput1.append(p_1[0:2])
	#kdTree1 = cKDTree(array(kdInput1))

	
	poly1 = []		  
	for p in points1:
		poly1.append([p[0],p[1]])	 

	while True:
		
		" find the matching pairs "
		match_pairs = []

		" transform the target Hull with the latest offset "
		points2_offset = []
		for p in points2:
			result = dispPoint(p, offset)		 
			#result = dispPointWithAngle(p, offset)
			points2_offset.append(result)
		
		" transformed points without associated covariance "
		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])	 
		
		" get the circles and radii "
		radius2, center2 = computeEnclosingCircle(points2_offset)
		radius1, center1 = computeEnclosingCircle(points1)
	
		for i in range(len(points2_offset)):
			p_2 = poly2[i]

			if isInCircle(p_2, radius1, center1):

				" for every transformed point of A, find it's closest neighbor in B "
				#p_1, minDist = findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
				p_1, i_1, minDist = findClosestPointInA(points1, p_2)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2[i][2]
					C1 = p_1[2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2[i],p_1,C2,C1])

		for i in range(len(points1)):
			p_1 = poly1[i]
	
			if isInCircle(p_1, radius2, center2):
		
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = findClosestPointInA(points2_offset, p_1)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2[i_2][2]
					
					C1 = points1[i][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2[i_2],points1[i],C2,C1])

		if plotIter:
			
			" set the origin of pose 1 "
			poseOrigin = Pose(estPose1)
	
			pylab.clf()
			pylab.axes()
			match_global = []
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_o = dispOffset(p1, offset)
				
				p1_g = poseOrigin.convertLocalToGlobal(p1_o)
				p2_g = poseOrigin.convertLocalToGlobal(p2)
				match_global.append([p1_g,p2_g])

			draw_matches(match_global, [0.0,0.0,0.0])

			
			xP = []
			yP = []
			for b in poly1:
				p1 = poseOrigin.convertLocalToGlobal(b)
	
				xP.append(p1[0])	
				yP.append(p1[1])
	
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			xP = []
			yP = []
			for b in points2:
				p = [b[0],b[1]]
				p = dispOffset(p,offset)
				
				p1 = poseOrigin.convertLocalToGlobal(p)
				xP.append(p1[0])	
				yP.append(p1[1])

			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

			point2_trans = dispOffset(point2, offset)
			p1 = poseOrigin.convertLocalToGlobal(point1)
			p2 = poseOrigin.convertLocalToGlobal(point2_trans)

			xP = [p1[0],p2[0]]
			yP = [p1[1],p2[1]]
			pylab.scatter(xP,yP,color='b')


			gpac2_trans = dispOffset(rootPose1, offset)
			gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
			gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

			xP = [gpac1[0],gpac2[0]]
			yP = [gpac1[1],gpac2[1]]
			pylab.scatter(xP,yP,color='r')


			plotEnv()		 
			
			lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)			
			pylab.title("%s %s, ang = %1.3f, hist = %d %d %d" % (repr(n1), repr(n2), currAng, lowCount, midCount, highCount))
			pylab.xlim(estPose1[0]-4, estPose1[0]+4)					
			pylab.ylim(estPose1[1]-3, estPose1[1]+3)
			pylab.savefig("ICP_plot_%06u_0.png" % numIterations)
			pylab.clf()			   
		

		newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, poses_1, poses_2, uHigh, uLow, u1], disp = 0)
		newU = newParam[0]
		newAng = newParam[1]
		
		# get the current cost
		newCost = medialOverlapCostFunc([newU, newAng], match_pairs, poses_1, poses_2, uHigh, uLow, u1)
		
		" set the current parameters "
		currAng = functions.normalizeAngle(newAng)
		currU = newU
		
		if inPlace:
			print "currU =", currU, "currAng =", currAng, "newCost =", newCost

		" compute offset from newU and newAng"

		if u1+0.02 > 1.0:
			pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
		elif u1 < 0.0:
			pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
		else:
			pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
	
		if currU+0.02 > 1.0:
			pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
		elif currU < 0.0:
			pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
		else:
			pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

	
		point1 = [pose1[0],pose1[1]]
		point2 = [pose2[0],pose2[1]]
		ang1 = pose1[2]
		ang2 = pose2[2]
	
		newOffset = computeOffset(point1, point2, ang1, ang2 + currAng)

		# save the current offset and cost
		offset = newOffset

	
		# check for convergence condition, different between last and current cost is below threshold
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			break



		" update after check for termination condition "
		lastCost = newCost

		# optionally draw the position of the points in current transform
		if plotIter:
			
			" set the origin of pose 1 "
			poseOrigin = Pose(estPose1)
	
			pylab.clf()
			pylab.axes()
			match_global = []
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_o = dispOffset(p1, offset)
				
				p1_g = poseOrigin.convertLocalToGlobal(p1_o)
				p2_g = poseOrigin.convertLocalToGlobal(p2)
				match_global.append([p1_g,p2_g])
			
			draw_matches(match_global, [0.0,0.0,0.0])
			
			xP = []
			yP = []
			for b in poly1:
				p1 = poseOrigin.convertLocalToGlobal(b)

				xP.append(p1[0])	
				yP.append(p1[1])

			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			xP = []
			yP = []
			for b in points2:
				p = [b[0],b[1]]
				p = dispOffset(p,offset)
				
				p1 = poseOrigin.convertLocalToGlobal(p)
				xP.append(p1[0])	
				yP.append(p1[1])
			
			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))


			point2_trans = dispOffset(point2, offset)
			p1 = poseOrigin.convertLocalToGlobal(point1)
			p2 = poseOrigin.convertLocalToGlobal(point2_trans)

			xP = [p1[0],p2[0]]
			yP = [p1[1],p2[1]]

			pylab.scatter(xP,yP,color='b')

			gpac2_trans = dispOffset(rootPose1, offset)
			gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
			gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

			xP = [gpac1[0],gpac2[0]]
			yP = [gpac1[1],gpac2[1]]
			pylab.scatter(xP,yP,color='r')

			plotEnv()		 
			
			lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)			
			pylab.title("%s %s, ang = %1.3f, hist = %d %d %d" % (repr(n1), repr(n2), currAng, lowCount, midCount, highCount))
			pylab.xlim(estPose1[0]-4, estPose1[0]+4)					
			pylab.ylim(estPose1[1]-3, estPose1[1]+3)
			pylab.savefig("ICP_plot_%06u_1.png" % numIterations)
			pylab.clf()			   

		
		numIterations += 1

		" reduce the minMatch distance for each step down to a floor value "
	
	histogram = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)			


	offset[2] =  functions.normalizeAngle(offset[2])	
	return offset, histogram



@logFunction
def overlapICP_GPU2(estPose1, initGuess, medialPoints1, medialPoints2, rootPose1, rootPose2, inPlace = False, plotIter = False, n1 = 0, n2 = 0, uRange = 1.5):

	global numIterations
	global globalPlotCount

	u1 = initGuess[0]
	u2 = initGuess[1]
	
	currU = u2
	currAng = initGuess[2]
	medialSpline1 = SplineFit(medialPoints1, smooth=0.1)
	medialSpline2 = SplineFit(medialPoints2, smooth=0.1)


	uSet = [i*0.01 for i in range(100)]
	poses_1 = medialSpline1.getUVecSet(uSet)
	poses_2 = medialSpline2.getUVecSet(uSet)

	uHigh = medialSpline2.getUOfDist(u2, uRange, distIter = 0.001)
	uLow = medialSpline2.getUOfDist(u2, -uRange, distIter = 0.001)
	
	if inPlace:
		uHigh = u2 + 0.08
		uLow = u2 - 0.08

	if u1 >= 1.0:
		pose1 = poses_1[-1]
	elif u1 < 0.0:
		pose1 = poses_1[0]
	else:
		pose1 = poses_1[int(u1*100)]

	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]

	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]

	offset = computeOffset(point1, point2, ang1, ang2 + currAng)

	costThresh = 0.004
	minMatchDist = 2.0
	lastCost = 1e100
	
	startIteration = numIterations

	" set the initial guess "
	poseOrigin = Pose(estPose1)
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "

	points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
	points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)

	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	points2_offset = []
	for p in points2:
		result = dispPoint(p, offset)		 
		points2_offset.append(result)


	
	#kdInput1 = []
	#for i in range(0,len(points1)):
	#	p_1 = points1[i]
	#	kdInput1.append(p_1[0:2])
	#kdTree1 = cKDTree(array(kdInput1))

	poly1 = []		  
	for p in points1:
		poly1.append([p[0],p[1]])	 


	while True:
		
		" find the matching pairs "
		match_pairs = []

		" transform the target Hull with the latest offset "
		points2_offset = []
		for p in points2:
			result = dispPoint(p, offset)		 
			points2_offset.append(result)
		
		" transformed points without associated covariance "
		poly2 = []
		for p in points2_offset:
			poly2.append([p[0],p[1]])	 
		
		" get the circles and radii "
		radius2, center2 = computeEnclosingCircle(points2_offset)
		radius1, center1 = computeEnclosingCircle(points1)
	
		for i in range(len(points2_offset)):
			p_2 = poly2[i]

			if isInCircle(p_2, radius1, center1):

				" for every transformed point of A, find it's closest neighbor in B "
				#p_1, minDist = findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
				p_1, i_1, minDist = findClosestPointInA(points1, p_2)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2[i][2]
					C1 = p_1[2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2[i],p_1,C2,C1])

		for i in range(len(points1)):
			p_1 = poly1[i]
	
			if isInCircle(p_1, radius2, center2):
		
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = findClosestPointInA(points2_offset, p_1)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points2[i_2][2]
					
					C1 = points1[i][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points2[i_2],points1[i],C2,C1])




		flatMatchPairs = []
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			C1 = pair[2]
			C2 = pair[3]

			flatMatchPairs.append(p1[0])
			flatMatchPairs.append(p1[1])
			flatMatchPairs.append(C1[0][0])
			flatMatchPairs.append(C1[0][1])
			flatMatchPairs.append(C1[1][0])
			flatMatchPairs.append(C1[1][1])
			flatMatchPairs.append(p2[0])
			flatMatchPairs.append(p2[1])
			flatMatchPairs.append(C2[0][0])
			flatMatchPairs.append(C2[0][1])
			flatMatchPairs.append(C2[1][0])
			flatMatchPairs.append(C2[1][1])
		c_poses_1 = [item for sublist in poses_1 for item in sublist]
		c_poses_2 = [item for sublist in poses_2 for item in sublist]
		newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), initGuess, uHigh, uLow, c_poses_1, c_poses_2, len(poses_1))

		newU = newParam[0]
		newAng = newParam[1]
		
		" set the current parameters "
		currAng = functions.normalizeAngle(newAng)
		currU = newU
		
		if inPlace:
			print "currU =", currU, "currAng =", currAng, "newCost =", newCost

		" compute offset from newU and newAng"

		if u1+0.02 > 1.0:
			pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
		elif u1 < 0.0:
			pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
		else:
			pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
	
		if currU+0.02 > 1.0:
			pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
		elif currU < 0.0:
			pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
		else:
			pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

		point1 = [pose1[0],pose1[1]]
		point2 = [pose2[0],pose2[1]]
		ang1 = pose1[2]
		ang2 = pose2[2]
	
		newOffset = computeOffset(point1, point2, ang1, ang2 + currAng)

		# save the current offset and cost
		offset = newOffset

	
		# check for convergence condition, different between last and current cost is below threshold
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			break



		" update after check for termination condition "
		lastCost = newCost
		
		numIterations += 1

	if False:

		" set the origin of pose 1 "
		poseOrigin = Pose(estPose1)

		fig1 = pylab.figure(1)
		fig1.clf()

		axes2 = fig1.add_subplot(111)

		xP = []
		yP = []
		for b in poly1:
			p1 = poseOrigin.convertLocalToGlobal(b)

			xP.append(p1[0])	
			yP.append(p1[1])

		axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

		xP = []
		yP = []
		for b in points2:
			p = [b[0],b[1]]
			p = dispOffset(p,offset)
			
			p1 = poseOrigin.convertLocalToGlobal(p)
			xP.append(p1[0])	
			yP.append(p1[1])

		axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

		point2_trans = dispOffset(point2, offset)
		p1 = poseOrigin.convertLocalToGlobal(point1)
		p2 = poseOrigin.convertLocalToGlobal(point2_trans)

		xP = [p1[0],p2[0]]
		yP = [p1[1],p2[1]]
		axes2.scatter(xP,yP,color='b')


		gpac2_trans = dispOffset(rootPose1, offset)
		gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
		gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

		xP = [gpac1[0],gpac2[0]]
		yP = [gpac1[1],gpac2[1]]
		axes2.scatter(xP,yP,color='r')

		plotEnv(axes2)		  
		
		axes2.set_title("%s %s, ang = %1.3f" % (repr(n1), repr(n2), currAng))
		axes2.set_xlim(estPose1[0]-4, estPose1[0]+4)					
		axes2.set_ylim(estPose1[1]-3, estPose1[1]+3)
		fig1.savefig("ICP_plot_%06u_0.png" % globalPlotCount)
		globalPlotCount += 1
		pylab.clf()			   
		pylab.figure(1)

		numIterations += 1


		
	histogram = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)			


	offset[2] =  functions.normalizeAngle(offset[2])	
	return offset, histogram


@logFunction
def globalOverlapICP_GPU2(initGuess, globalPath, medialPoints,plotIter = False, n1 = 0, n2 = 0, minMatchDist = 1.0, arcLimit = 0.5):

	global numIterations
	global globalPlotCount
	
	u1 = initGuess[0]
	u2 = initGuess[1]

	uHigh = u2 + 0.2
	uLow = u2 - 0.2

	currU = u2
	currAng = initGuess[2]
	
	globalSpline = SplineFit(globalPath, smooth=0.1)
	medialSpline = SplineFit(medialPoints, smooth=0.1)
	
	
	uSet = [i*0.01 for i in range(100)]
	poses_1 = globalSpline.getUVecSet(uSet)
	poses_2 = medialSpline.getUVecSet(uSet)
	
	if u1 >= 1.0:
		pose1 = poses_1[-1]
	elif u1 < 0.0:
		pose1 = poses_1[0]
	else:
		pose1 = poses_1[int(u1*100)]

	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]
	
	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]

	currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

	costThresh = 0.004
	lastCost = 1e100
	
	startIteration = numIterations

	" set the initial guess "
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
	points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)

	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	points_offset = []

	globalPoly = []		   
	for p in globalPoints:
		globalPoly.append([p[0],p[1]])	  


	while True:
		
		" find the matching pairs "
		match_pairs = []

		" transform the target Hull with the latest offset "
		points_offset = []
		for p in points:
			result = dispPoint(p, currPose)		   
			points_offset.append(result)
		

		match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)		  
		
		
		flatMatchPairs = []
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			C1 = pair[2]
			C2 = pair[3]

			flatMatchPairs.append(p1[0])
			flatMatchPairs.append(p1[1])
			flatMatchPairs.append(C1[0][0])
			flatMatchPairs.append(C1[0][1])
			flatMatchPairs.append(C1[1][0])
			flatMatchPairs.append(C1[1][1])
			flatMatchPairs.append(p2[0])
			flatMatchPairs.append(p2[1])
			flatMatchPairs.append(C2[0][0])
			flatMatchPairs.append(C2[0][1])
			flatMatchPairs.append(C2[1][0])
			flatMatchPairs.append(C2[1][1])
		c_poses_1 = [item for sublist in poses_1 for item in sublist]
		c_poses_2 = [item for sublist in poses_2 for item in sublist]
		

		
		newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), [u1,currU,currAng], uHigh, uLow, c_poses_1, c_poses_2, len(poses_1))

		newU = newParam[0]
		newAng = newParam[1]
		
		" set the current parameters "
		currAng = functions.normalizeAngle(newAng)
		currU = newU


		if currU >= 1.0:
			pose2 = poses_2[-1]
		elif currU < 0.0:
			pose2 = poses_2[0]
		else:
			pose2 = poses_2[int(currU*100)]


		point1 = [pose1[0],pose1[1]]
		point2 = [pose2[0],pose2[1]]
		ang1 = pose1[2]
		ang2 = pose2[2]
	
		currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

		# check for convergence condition, different between last and current cost is below threshold
		
		isTerminate = False
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			isTerminate = True
			print "terminating globalOverlap_GPU2:", lastCost, newCost, costThresh, numIterations-startIteration
		

		" update after check for termination condition "
		lastCost = newCost

		numIterations += 1
		
		if isTerminate:
			break

		" draw final position "
		if plotIter:
			
			" set the origin of pose 1 "
			poseOrigin = Pose(currPose)
			
			pylab.clf()
			pylab.axes()
			match_global = []
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_o = dispOffset(p1, currPose)
				
				p1_g = p1_o
				p2_g = p2
				match_global.append([p1_g,p2_g])
	
			draw_matches(match_global, [0.0,0.0,0.0])
	
			
			xP = []
			yP = []
			for b in globalPoly:
				p1 = b
				xP.append(p1[0])	
				yP.append(p1[1])
	
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
	
			pylab.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))
	
			xP = []
			yP = []
			for b in points:
				p = [b[0],b[1]]
				p1 = dispOffset(p,currPose)
				
				xP.append(p1[0])	
				yP.append(p1[1])
	
			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
	
			plotEnv()		 
			pylab.title("(%u,%u) u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, n2, u1, currU, currAng, newCost))
	
			pylab.xlim(currPose[0]-4, currPose[0]+4)					
			pylab.ylim(currPose[1]-3, currPose[1]+3)
			pylab.savefig("ICP_plot_%06u.png" % globalPlotCount)
			pylab.clf()
			
			globalPlotCount += 1

	return currPose, newCost, len(match_pairs)

@logFunction
def globalPathToNodeOverlapICP2(initGuess, globalPath, medialPoints, plotIter = False, n1 = 0, n2 = 0, minMatchDist = 1.0, arcLimit = 0.5, origPose = [0.0,0.0,0.0]):

	global numIterations
	global globalPlotCount
	
	u1 = initGuess[0]
	u2 = initGuess[1]
	u3 = u1

	
	currU = u3
	currAng = initGuess[2]

	globalSpline = SplineFit(globalPath, smooth=0.1)
	medialSpline = SplineFit(medialPoints, smooth=0.1)

	uHigh = globalSpline.getUOfDist(u3, arcLimit, distIter = 0.001)
	uLow = globalSpline.getUOfDist(u3, -arcLimit, distIter = 0.001)

	" keep the rails within the range of u parameterization "
	if uHigh >= 1.0:
		uHigh = 1.0
	
	if uLow <= 0.0:
		uLow = 0.0

	print "u1,u2,u3,uHigh,uLow,currU,currAng:", u1, u2, u3, uHigh, uLow, currU, currAng
	
	
	uSet = [i*0.01 for i in range(100)]
	poses_1 = globalSpline.getUVecSet(uSet)
	poses_2 = medialSpline.getUVecSet(uSet)
	
	if u1 >= 1.0:
		pose1 = poses_1[-1]
	elif u1 < 0.0:
		pose1 = poses_1[0]
	else:
		pose1 = poses_1[int(u1*100)]

	if u2 >= 1.0:
		pose2 = poses_2[-1]
	elif u2 < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(u2*100)]


	pose2_orig = pose2
		
	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]


	toLocalOffset = computeOffset(point2, point1, ang2, ang1 + currAng)
	poseRoot = Pose([0.0,0.0,0.0])
	deOffset = poseRoot.doInverse(toLocalOffset)
	globalPath2 = []
	pathLocal = []
	for p in globalPath:
		result = dispOffset(p, toLocalOffset)
		pathLocal.append(result)
		result2 = dispOffset(result, deOffset)
		globalPath2.append(result2)


	poses_3 = []
	p3Origin = Pose(toLocalOffset)
	p3Restore = Pose(deOffset)
	for p1 in poses_1:
		p3 = p3Origin.convertLocalOffsetToGlobal(p1)
		poses_3.append(p3)


	globalPoints = globalSpline.getUniformSamples()
	localPathPoints = []
	for p1 in globalPoints:
		p3 = p3Origin.convertLocalOffsetToGlobal(p1)
		localPathPoints.append(p3)

	localPathPoints = addGPACVectorCovariance(localPathPoints,high_var=0.05, low_var = 0.001)

	if u3 >= 1.0:
		pose3 = poses_3[-1]
	elif u3 < 0.0:
		pose3 = poses_3[0]
	else:
		pose3 = poses_3[int(u3*100)]


	point3 = [pose3[0],pose3[1]]
	point2 = [pose2[0],pose2[1]]
	ang3 = pose3[2]
	ang2 = pose2[2]

	currPose = computeOffset(point2, point3, ang2, ang3 + currAng)





	" transform the new pose "
	
	costThresh = 0.004
	lastCost = 1e100
	
	startIteration = numIterations

	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)

	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	points_offset = []

	globalPoly = []		   
	for p in globalPoints:
		globalPoly.append([p[0],p[1]])	  


	medialPoly = []		   
	for p in points:
		medialPoly.append([p[0],p[1]])	  


	while True:
		
		" find the matching pairs "
		match_pairs = []


		" transform the target Hull with the latest offset "
		localPathPoints_offset = []
		for p in localPathPoints:
			result = dispPoint(p, currPose)		   
			localPathPoints_offset.append(result)
 
		match_pairs = matchPairs(localPathPoints, localPathPoints_offset, points, minMatchDist)		  


		flatMatchPairs = []
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			C1 = pair[2]
			C2 = pair[3]

			flatMatchPairs.append(p1[0])
			flatMatchPairs.append(p1[1])
			flatMatchPairs.append(C1[0][0])
			flatMatchPairs.append(C1[0][1])
			flatMatchPairs.append(C1[1][0])
			flatMatchPairs.append(C1[1][1])
			flatMatchPairs.append(p2[0])
			flatMatchPairs.append(p2[1])
			flatMatchPairs.append(C2[0][0])
			flatMatchPairs.append(C2[0][1])
			flatMatchPairs.append(C2[1][0])
			flatMatchPairs.append(C2[1][1])
		c_poses_2 = [item for sublist in poses_2 for item in sublist]
		c_poses_3 = [item for sublist in poses_3 for item in sublist]

		#if plotIter and pose2[0] > 1.8:
		#	print "input:", len(match_pairs), [u2, currU, currAng], uHigh, uLow, pose2, len(poses_2)
		newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), [u2,currU,currAng], uHigh, uLow, c_poses_2, c_poses_3, len(poses_2))

		#if plotIter and pose2[0] > 1.8:
		#	print poses_2

		print "newCost =", newCost

		newU = newParam[0]
		newAng = newParam[1]

		" set the current parameters "
		currAng = functions.normalizeAngle(newAng)
		currU = newU
			
		if currU >= 1.0:
			pose3 = poses_3[-1]
		elif currU < 0.0:
			pose3 = poses_3[0]
		else:
			pose3 = poses_3[int(currU*100)]


		point3 = [pose3[0],pose3[1]]
		point2 = [pose2[0],pose2[1]]
		ang3 = pose3[2]
		ang2 = pose2[2]
	
		currPose = computeOffset(point2, point3, ang2, ang3 + currAng)

		" check for convergence condition, different between last and current cost is below threshold "
		
		isTerminate = False
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			isTerminate = True
			print "terminating globalOverlap_GPU2:", lastCost, newCost, costThresh, numIterations-startIteration
		
		" update after check for termination condition "
		lastCost = newCost

		numIterations += 1
		
		
		doPoseOffset = Pose(currPose)
		unPoseOffset = doPoseOffset.doInverse(currPose)
		doUnPoseOffset = Pose(unPoseOffset)
			

		if plotIter:

			trueCost = shapeCostC(currPose, match_pairs)
			numPairs = len(match_pairs)
			numPoses = len(poses_3)
			trueCost2, resultParam, resultOffset = nelminICP.ICPcost(flatMatchPairs, numPairs, [u2,currU,currAng], uHigh, uLow, c_poses_2, c_poses_3, numPoses)
			trueCost3 = shapeCostC(resultOffset, match_pairs)

			pylab.clf()

			pLocal = doUnPoseOffset.convertLocalOffsetToGlobal([0.0,0.0,0.0])
			resultPose = p3Restore.convertLocalOffsetToGlobal(pLocal)

			xP = []
			yP = []
			for b in medialPoly:
				p1 = dispOffset(b, resultPose)
				xP.append(p1[0])	
				yP.append(p1[1])
			

			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

			xP = []
			yP = []
			for b in medialPoly:
				p1 = dispOffset(b, origPose)
				xP.append(p1[0])	
				yP.append(p1[1])

			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0), alpha=0.5)

			
			if currU >= 1.0:
				newPose1 = poses_1[-1]
			elif currU < 0.0:
				newPose1 = poses_1[0]
			else:
				newPose1 = poses_1[int(currU*100)]

			pylab.scatter([resultPose[0]],[resultPose[1]],color=(0.0,0.0,1.0))
			pylab.scatter([newPose1[0]],[newPose1[1]],color='k')

			zeroPose = dispOffset([0.0,0.0,0.0], resultPose)
			pylab.scatter([zeroPose[0]],[zeroPose[1]],color='y')
			
			xP = []
			yP = []
			for b in globalPath:
				xP.append(b[0])    
				yP.append(b[1])
			
			pylab.title("(%u,%u) u2 = %1.3f, u3 = %1.3f, ang = %1.3f, iter = %d" % (n1, n2, u2, currU, currAng, numIterations))
			pylab.xlim(-5, 10)					   
			pylab.ylim(-8, 8)					   
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))	
			pylab.savefig("ICP_plot_%06u_000.png" % globalPlotCount)

			globalPlotCount += 1

		if isTerminate:
			break


	if False:

		if False:

			pylab.clf()

			" transform the target Hull with the latest offset "
			localPathPoints_offset = []
			for p in localPathPoints:
				result = dispPoint(p, currPose)
				localPathPoints_offset.append(result)
			

			
			" set the origin of pose 1 "
			poseOrigin = Pose(currPose)
			
			match_global = []
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_g = dispOffset(p1, deOffset)
				
				p2Local = dispOffset(p2, unPoseOffset)
				p2_g = dispOffset(p2Local, deOffset)				
				 
				match_global.append([p1_g,p2_g])
			
			draw_matches(match_global, [0.0,0.0,0.0])
			
			xP = []
			yP = []
			for b in globalPoly:
				p1 = b
				xP.append(p1[0])	
				yP.append(p1[1])
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
			
			
			pylab.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))
			
			xP = []
			yP = []
			for b in localPathPoints_offset:
				pLocal = dispOffset(b, unPoseOffset)
				p1 = dispOffset(pLocal, deOffset)				 
				xP.append(p1[0])	
				yP.append(p1[1])
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,1.0,0.0))
			
			p3 = p3Restore.convertLocalOffsetToGlobal(pose3)
			
			pylab.scatter([p3[0]],[p3[1]],color=(0.0,1.0,0.0))
			
			
			xP = []
			yP = []
			for b in points:
				p1 = [b[0],b[1]]
				
				
				p1Local = dispOffset(p1, unPoseOffset)
				p1_g = dispOffset(p1Local, deOffset)					 
				xP.append(p1_g[0])	  
				yP.append(p1_g[1])
			  
			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
			
			p2Local = dispOffset(pose2, unPoseOffset)
			p2 = dispOffset(p2Local, deOffset)	 
			pylab.scatter([p2[0]],[p2[1]],color=(1.0,0.0,0.0))
			
			plotEnv()		
			pylab.title("(%u,%u) u1 = %1.3f, U = %1.3f, A = %1.3f, cost = %1.3f, %1.3f, %1.3f" % (n1, n2, u1, currU, currAng, newCost, trueCost, trueCost2))

			pylab.xlim(-5, 10)					   
			pylab.ylim(-8, 8)					   

			pylab.savefig("ICP_plot_%06u_001.png" % globalPlotCount)


			pylab.clf()

			" set the origin of pose 1 "
			poseOrigin = Pose(currPose)
			
			match_global = []
			
			for pair in match_pairs:
				p2 = pair[0]
				p3 = pair[1]
				
				p2_o = dispOffset(p2, currPose)
				
				p2_g = p2_o
				p3_g = p3
				match_global.append([p2_g,p3_g])
			
			draw_matches(match_global, [0.0,0.0,0.0])
			
			
			xP = []
			yP = []
			for b in medialPoly:
				p1 = b
				xP.append(p1[0])	
				yP.append(p1[1])
			
			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
			
			pylab.scatter([pose2[0]],[pose2[1]],color=(0.0,0.0,1.0))
			
			xP = []
			yP = []
			for b in localPathPoints_offset:
				xP.append(b[0])    
				yP.append(b[1])
			
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
			
			
			plotEnv()		
			pylab.title("(%u,%u) u2 = %1.3f, u3 = %1.3f, ang = %1.3f, cost = %f" % (n1, n2, u2, currU, currAng, trueCost))
			
			pylab.xlim(-5, 10)					   
			pylab.ylim(-8, 8)					   
			pylab.savefig("ICP_plot_%06u_002.png" % globalPlotCount)

		globalPlotCount += 1

	" draw final position "
	if False:
		
		
		trueCost = shapeCostC(currPose, match_pairs)
		
		numPairs = len(match_pairs)
		numPoses = len(poses_3)
		trueCost2, resultParam, resultOffset = nelminICP.ICPcost(flatMatchPairs, numPairs, [u2,currU,currAng], uHigh, uLow, c_poses_2, c_poses_3, numPoses)
		trueCost3 = shapeCostC(resultOffset, match_pairs)
		
		print "currPose:", currPose
		print "resultOffset:", resultOffset
		print "trueCost:", trueCost
		print "trueCost2:", trueCost2
		print "trueCost3:", trueCost3
		
		
		figPlot, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
		
		#doPoseOffset = Pose(currPose)
		#unPoseOffset = doPoseOffset.doInverse(currPose)
		#doUnPoseOffset = Pose(unPoseOffset)
		
		" transform the target Hull with the latest offset "
		localPathPoints_offset = []
		for p in localPathPoints:
			result = dispPoint(p, currPose)
			localPathPoints_offset.append(result)
		
		
		" set the origin of pose 1 "
		poseOrigin = Pose(currPose)
		
		match_global = []
		
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			
			p1_g = dispOffset(p1, deOffset)
			
			p2Local = dispOffset(p2, unPoseOffset)
			p2_g = dispOffset(p2Local, deOffset)				
			 
			match_global.append([p1_g,p2_g])
		
		draw_matches(match_global, [0.0,0.0,0.0], ax1)
		
		xP = []
		yP = []
		for b in globalPoly:
			p1 = b
			xP.append(p1[0])	
			yP.append(p1[1])
		
		ax1.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
		
		
		ax1.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))
		
		xP = []
		yP = []
		for b in localPathPoints_offset:
			pLocal = dispOffset(b, unPoseOffset)
			p1 = dispOffset(pLocal, deOffset)				 
			xP.append(p1[0])	
			yP.append(p1[1])
		
		ax1.plot(xP,yP,linewidth=1, color=(0.0,1.0,0.0))
		
		p3 = p3Restore.convertLocalOffsetToGlobal(pose3)
		
		ax1.scatter([p3[0]],[p3[1]],color=(0.0,1.0,0.0))
		
		
		xP = []
		yP = []
		for b in points:
			p1 = [b[0],b[1]]
			
			
			p1Local = dispOffset(p1, unPoseOffset)
			p1_g = dispOffset(p1Local, deOffset)					 
			xP.append(p1_g[0])	  
			yP.append(p1_g[1])
		  
		ax1.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
		
		p2Local = dispOffset(pose2, unPoseOffset)
		p2 = dispOffset(p2Local, deOffset)	 
		ax1.scatter([p2[0]],[p2[1]],color=(1.0,0.0,0.0))
		
		plotEnv(ax1)		
		ax1.set_title("(%u,%u) u1 = %1.3f, U = %1.3f, A = %1.3f, cost = %1.3f, %1.3f, %1.3f" % (n1, n2, u1, currU, currAng, newCost, trueCost, trueCost2))

		ax1.set_xlim(-4, 4)					   
		ax1.set_ylim(-3, 3)

		" set the origin of pose 1 "
		poseOrigin = Pose(currPose)
		
		match_global = []
		
		for pair in match_pairs:
			p2 = pair[0]
			p3 = pair[1]
			
			p2_o = dispOffset(p2, currPose)
			
			p2_g = p2_o
			p3_g = p3
			match_global.append([p2_g,p3_g])
		
		draw_matches(match_global, [0.0,0.0,0.0], ax2)
		
		
		xP = []
		yP = []
		for b in medialPoly:
			p1 = b
			xP.append(p1[0])	
			yP.append(p1[1])
		
		ax2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
		
		ax2.scatter([pose2[0]],[pose2[1]],color=(0.0,0.0,1.0))
		
		xP = []
		yP = []
		for b in localPathPoints_offset:
			xP.append(b[0])    
			yP.append(b[1])
		
		
		ax2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
		
		
		plotEnv(ax2)		
		ax2.set_title("(%u,%u) u2 = %1.3f, u3 = %1.3f, ang = %1.3f, cost = %f" % (n1, n2, u2, currU, currAng, trueCost))
		
		ax2.set_xlim(-4, 4)					   
		ax2.set_ylim(-3, 3)
		plt.savefig("ICP_plot_%06u.png" % globalPlotCount)
		plt.clf()
		plt.close()
		
		
		globalPlotCount += 1

	pLocal = doUnPoseOffset.convertLocalOffsetToGlobal([0.0,0.0,0.0])
	resultPose = p3Restore.convertLocalOffsetToGlobal(pLocal)

	return resultPose, newCost, len(match_pairs), currAng, currU


@logFunction
def branchEstimateCost(initGuess, junctionPose, pathSoup, globalPath, plotIter = False, n1 = 0, n2 = 0):

	global numIterations
	global globalPlotCount
	
	u1 = initGuess[0]
	u2 = initGuess[1]

	uHigh = u2 + 0.2
	uLow = u2 - 0.2

	currU = u2
	currAng = initGuess[2]

	angNom = currAng
	angLim = pi/4.0
	#angLim = pi/64.0
	
	startU2 = u2
	
	globalSpline = SplineFit(globalPath, smooth=0.1)
	

	uSet = [i*0.01 for i in range(100)]
	poses_2 = globalSpline.getUVecSet(uSet)
	
	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]
	
	point1 = [junctionPose[0],junctionPose[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = junctionPose[2]
	ang2 = pose2[2]

	currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

	costThresh = 0.004
	minMatchDist = 0.5
	minMatchDist2 = 0.5
	lastCost = 1e100
	matchAngTol = math.pi/4.0
	
	startIteration = numIterations

	" set the initial guess "
	poseOrigin = Pose(currPose)
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	globalVecPoints = globalSpline.getUniformSamples()
	#globalCovPoints = addGPACVectorCovariance(globalVecPoints,high_var=0.50, low_var = 0.10)
	globalCovPoints = addGPACVectorCovariance(globalVecPoints,high_var=1.00, low_var = 1.00)

	pointCovSoup = []
	pointVecSoup = []
	for path in pathSoup:
		#points = addGPACVectorCovariance(path, high_var=0.50, low_var = 0.10)
		points = addGPACVectorCovariance(path, high_var=1.00, low_var = 1.00)
		pointCovSoup += deepcopy(points)
		pointVecSoup += deepcopy(path)
	
	" transform pose 2 by initial offset guess "	
	" transform the new pose "

	
	" find the matching pairs "
	match_pairs = []

	" transform the target Hull with the latest offset "
	poseOrigin = Pose(currPose)
	globalVecPoints_offset = []
	for p in globalVecPoints:
		result = poseOrigin.convertLocalOffsetToGlobal(p)
		globalVecPoints_offset.append(result)


	" get the circles and radii "
	match_pairs = []

	for i in range(len(globalVecPoints)):
		p_1 = globalVecPoints_offset[i]

		try:
			p_2, i_2, minDist = findClosestPointWithAngle(pointVecSoup, p_1, matchAngTol)

			if minDist <= minMatchDist:
	
				" add to the list of match pairs less than 1.0 distance apart "
				" keep A points and covariances untransformed "
				C1 = globalCovPoints[i][2]
				C2 = pointCovSoup[i_2][2]

				" we store the untransformed point, but the transformed covariance of the A point "
				match_pairs.append([globalVecPoints[i],pointVecSoup[i_2],C1,C2])

		except:
			pass

	for i in range(len(pointVecSoup)):
		p_1 = pointVecSoup[i]

		try:
			p_2, i_2, minDist = findClosestPointWithAngle(globalVecPoints_offset, p_1, matchAngTol)

			if minDist <= minMatchDist2:
	
				" add to the list of match pairs less than 1.0 distance apart "
				" keep A points and covariances untransformed "
				C1 = globalCovPoints[i_2][2]	 
				C2 = pointCovSoup[i][2]						 

				" we store the untransformed point, but the transformed covariance of the A point "
				match_pairs.append([globalVecPoints[i_2],pointVecSoup[i],C1,C2])
		except:
			pass

	flatMatchPairs = []
	for pair in match_pairs:
		p1 = pair[0]
		p2 = pair[1]
		C1 = pair[2]
		C2 = pair[3]

		flatMatchPairs.append(p1[0])
		flatMatchPairs.append(p1[1])
		flatMatchPairs.append(C1[0][0])
		flatMatchPairs.append(C1[0][1])
		flatMatchPairs.append(C1[1][0])
		flatMatchPairs.append(C1[1][1])
		flatMatchPairs.append(p2[0])
		flatMatchPairs.append(p2[1])
		flatMatchPairs.append(C2[0][0])
		flatMatchPairs.append(C2[0][1])
		flatMatchPairs.append(C2[1][0])
		flatMatchPairs.append(C2[1][1])

	c_pose1 = junctionPose
	c_poses_2 = [item for sublist in poses_2 for item in sublist]
	numPoses = len(poses_2)


	newCost, newParam, newOffset  = nelminICP.pointCost(flatMatchPairs, len(match_pairs), [u1,currU,currAng], angNom, angLim, uHigh, uLow, c_pose1, c_poses_2, numPoses)

	#newParam, newCost = nelminICP.ICPbyPose(flatMatchPairs, len(match_pairs), [u1,currU,currAng], angNom, angLim, uHigh, uLow, c_pose1, c_poses_2, numPoses)



	" draw final position "
	if plotIter:
		
		" set the origin of pose 1 "
		poseOrigin = Pose(currPose)
		
		pylab.clf()
		pylab.axes()
		match_global = []
		
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			
			p1_o = dispOffset(p1, currPose)
			
			p1_g = p1_o
			p2_g = p2
			match_global.append([p1_g,p2_g])

		draw_matches(match_global, [0.0,0.0,0.0])

		
		for path in pathSoup:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])	
				yP.append(p[1])

			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

		pylab.scatter([junctionPose[0]],[junctionPose[1]],color=(0.0,0.0,1.0))

		xP = []
		yP = []
		for b in globalVecPoints:
			p = [b[0],b[1]]
			p1 = dispOffset(p,currPose)
			
			xP.append(p1[0])	
			yP.append(p1[1])

		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

		plotEnv()		 
		pylab.title("%u -- %s u1 = %1.3f, u2 = %1.3f, ang = %1.3f, %.1f, %d %d F" % (n1, repr(n2), u1, currU, currAng, newCost, len(match_pairs), numIterations))

		pylab.xlim(-6, 10)					  
		pylab.ylim(-8, 8)
		pylab.savefig("ICP_plot_%06u.png" % globalPlotCount)
		pylab.clf()
		
		" save inputs "
		globalPlotCount += 1

	return currPose, newCost, len(match_pairs)


@logFunction
def branchEstimateICP(initGuess, junctionPose, pathSoup, globalPath, plotIter = False, n1 = 0, n2 = 0):

	global numIterations
	global globalPlotCount
	
	u1 = initGuess[0]
	u2 = initGuess[1]

	uHigh = u2 + 0.2
	uLow = u2 - 0.2

	currU = u2
	currAng = initGuess[2]

	angNom = currAng
	angLim = pi/4.0
	#angLim = pi/64.0
	
	startU2 = u2
	
	globalSpline = SplineFit(globalPath, smooth=0.1)
	

	uSet = [i*0.01 for i in range(100)]
	poses_2 = globalSpline.getUVecSet(uSet)
	
	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]
	
	point1 = [junctionPose[0],junctionPose[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = junctionPose[2]
	ang2 = pose2[2]

	currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

	costThresh = 0.004
	minMatchDist = 0.5
	minMatchDist2 = 0.5
	lastCost = 1e100
	matchAngTol = math.pi/4.0
	
	startIteration = numIterations

	" set the initial guess "
	poseOrigin = Pose(currPose)
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	globalVecPoints = globalSpline.getUniformSamples()
	#globalCovPoints = addGPACVectorCovariance(globalVecPoints,high_var=0.50, low_var = 0.10)
	globalCovPoints = addGPACVectorCovariance(globalVecPoints,high_var=1.00, low_var = 1.00)

	pointCovSoup = []
	pointVecSoup = []
	for path in pathSoup:
		#points = addGPACVectorCovariance(path, high_var=0.50, low_var = 0.10)
		points = addGPACVectorCovariance(path, high_var=1.00, low_var = 1.00)
		pointCovSoup += deepcopy(points)
		pointVecSoup += deepcopy(path)
	
	" transform pose 2 by initial offset guess "	
	" transform the new pose "

	" draw final position "
	if plotIter:
		
		" set the origin of pose 1 "
		poseOrigin = Pose(currPose)
		
		pylab.clf()
		pylab.axes()
		
		for path in pathSoup:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])	
				yP.append(p[1])

			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

		pylab.scatter([junctionPose[0]],[junctionPose[1]],color=(0.0,0.0,1.0))

		xP = []
		yP = []
		for b in globalVecPoints:
			p = [b[0],b[1]]
			p1 = dispOffset(p,currPose)
			
			xP.append(p1[0])	
			yP.append(p1[1])

		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

		plotEnv()		 
		pylab.title("%u -- %s u1 = %1.3f, u2 = %1.3f, ang = %1.3f, %d" % (n1, repr(n2), u1, currU, currAng, numIterations))

		pylab.xlim(-6, 10)					  
		pylab.ylim(-8, 8)
		pylab.savefig("ICP_plot_%06u.png" % globalPlotCount)
		pylab.clf()
		
		" save inputs "
		globalPlotCount += 1

	
	while True:
		
		" find the matching pairs "
		match_pairs = []

		" transform the target Hull with the latest offset "
		poseOrigin = Pose(currPose)
		globalVecPoints_offset = []
		for p in globalVecPoints:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			globalVecPoints_offset.append(result)


		" get the circles and radii "
		match_pairs = []

		for i in range(len(globalVecPoints)):
			p_1 = globalVecPoints_offset[i]
	
			try:
				p_2, i_2, minDist = findClosestPointWithAngle(pointVecSoup, p_1, matchAngTol)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C1 = globalCovPoints[i][2]
					C2 = pointCovSoup[i_2][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([globalVecPoints[i],pointVecSoup[i_2],C1,C2])

			except:
				pass

		for i in range(len(pointVecSoup)):
			p_1 = pointVecSoup[i]

			try:
				p_2, i_2, minDist = findClosestPointWithAngle(globalVecPoints_offset, p_1, matchAngTol)
	
				if minDist <= minMatchDist2:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C1 = globalCovPoints[i_2][2]	 
					C2 = pointCovSoup[i][2]						 

					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([globalVecPoints[i_2],pointVecSoup[i],C1,C2])
			except:
				pass

		flatMatchPairs = []
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			C1 = pair[2]
			C2 = pair[3]

			flatMatchPairs.append(p1[0])
			flatMatchPairs.append(p1[1])
			flatMatchPairs.append(C1[0][0])
			flatMatchPairs.append(C1[0][1])
			flatMatchPairs.append(C1[1][0])
			flatMatchPairs.append(C1[1][1])
			flatMatchPairs.append(p2[0])
			flatMatchPairs.append(p2[1])
			flatMatchPairs.append(C2[0][0])
			flatMatchPairs.append(C2[0][1])
			flatMatchPairs.append(C2[1][0])
			flatMatchPairs.append(C2[1][1])

		c_pose1 = junctionPose
		c_poses_2 = [item for sublist in poses_2 for item in sublist]
		numPoses = len(poses_2)


		newParam, newCost = nelminICP.ICPbyPose(flatMatchPairs, len(match_pairs), [u1,currU,currAng], angNom, angLim, uHigh, uLow, c_pose1, c_poses_2, numPoses)



		newU = newParam[0]
		newAng = newParam[1]

		print "old to new:", u1, currU, currAng, " - ", newU, newAng
		
		" set the current parameters "
		currAng = functions.normalizeAngle(newAng)
		currU = newU


		if currU >= 1.0:
			pose2 = poses_2[-1]
		elif currU < 0.0:
			pose2 = poses_2[0]
		else:
			pose2 = poses_2[int(currU*100)]


		point1 = [junctionPose[0],junctionPose[1]]
		point2 = [pose2[0],pose2[1]]
		ang1 = junctionPose[2]
		ang2 = pose2[2]
	
		#currPose = computeOffset(point1, point2, ang1, ang2 + currAng)
		currPose = computeOffset(point1, point2, ang1, ang1)

		# check for convergence condition, different between last and current cost is below threshold
		
		isTerminate = False
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			isTerminate = True


		" update after check for termination condition "
		lastCost = newCost

		numIterations += 1
		
		if isTerminate:
			break
		
		" draw final position "
		if plotIter:
			
			" set the origin of pose 1 "
			poseOrigin = Pose(currPose)
			
			pylab.clf()
			pylab.axes()
			match_global = []
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_o = dispOffset(p1, currPose)
				
				p1_g = p1_o
				p2_g = p2
				match_global.append([p1_g,p2_g])

			draw_matches(match_global, [0.0,0.0,0.0])

			
			for path in pathSoup:
				xP = []
				yP = []
				for p in path:
					xP.append(p[0])	
					yP.append(p[1])

				pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			pylab.scatter([junctionPose[0]],[junctionPose[1]],color=(0.0,0.0,1.0))

			xP = []
			yP = []
			for b in globalVecPoints:
				p = [b[0],b[1]]
				p1 = dispOffset(p,currPose)
				
				xP.append(p1[0])	
				yP.append(p1[1])

			pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

			plotEnv()		 
			pylab.title("%u -- %s u1 = %1.3f, u2 = %1.3f, ang = %1.3f, %.1f, %d %d" % (n1, repr(n2), u1, currU, currAng, newCost, len(match_pairs), numIterations))

			pylab.xlim(-6, 10)					  
			pylab.ylim(-8, 8)
			pylab.savefig("ICP_plot_%06u.png" % globalPlotCount)
			pylab.clf()
			
			" save inputs "
			globalPlotCount += 1

	" draw final position "
	if plotIter:
		
		" set the origin of pose 1 "
		poseOrigin = Pose(currPose)
		
		pylab.clf()
		pylab.axes()
		match_global = []
		
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			
			p1_o = dispOffset(p1, currPose)
			
			p1_g = p1_o
			p2_g = p2
			match_global.append([p1_g,p2_g])

		draw_matches(match_global, [0.0,0.0,0.0])

		
		for path in pathSoup:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])	
				yP.append(p[1])

			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

		pylab.scatter([junctionPose[0]],[junctionPose[1]],color=(0.0,0.0,1.0))

		xP = []
		yP = []
		for b in globalVecPoints:
			p = [b[0],b[1]]
			p1 = dispOffset(p,currPose)
			
			xP.append(p1[0])	
			yP.append(p1[1])

		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

		plotEnv()		 
		pylab.title("%u -- %s u1 = %1.3f, u2 = %1.3f, ang = %1.3f, %.1f, %d %d F" % (n1, repr(n2), u1, currU, currAng, newCost, len(match_pairs), numIterations))

		pylab.xlim(-6, 10)					  
		pylab.ylim(-8, 8)
		pylab.savefig("ICP_plot_%06u.png" % globalPlotCount)
		pylab.clf()
		
		" save inputs "
		globalPlotCount += 1

	return currPose, newCost, len(match_pairs)



@logFunction
def pathOverlapICP(initGuess, globalPath, medialPoints,plotIter = False, n1 = 0, n2 = 0):

	global numIterations
	global globalPlotCount
	
	u1 = initGuess[0]
	u2 = initGuess[1]

	uHigh = u2 + 0.2
	uLow = u2 - 0.2

	currU = u2
	currAng = initGuess[2]
	
	startU2 = u2
	
	globalSpline = SplineFit(globalPath, smooth=0.1)
	medialSpline = SplineFit(medialPoints, smooth=0.1)
	

	uSet = [i*0.01 for i in range(100)]
	poses_1 = globalSpline.getUVecSet(uSet)
	poses_2 = medialSpline.getUVecSet(uSet)
	

	if u1 >= 1.0:
		pose1 = poses_1[-1]
	elif u1 < 0.0:
		pose1 = poses_1[0]
	else:
		pose1 = poses_1[int(u1*100)]

	if currU >= 1.0:
		pose2 = poses_2[-1]
	elif currU < 0.0:
		pose2 = poses_2[0]
	else:
		pose2 = poses_2[int(currU*100)]
	
	point1 = [pose1[0],pose1[1]]
	point2 = [pose2[0],pose2[1]]
	ang1 = pose1[2]
	ang2 = pose2[2]

	currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

	costThresh = 0.004
	minMatchDist = 1.0
	minMatchDist2 = 1.0
	lastCost = 1e100
	
	startIteration = numIterations

	" set the initial guess "
	poseOrigin = Pose(currPose)
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	globalVecPoints = globalSpline.getUniformSamples()
	globalPoints = addGPACVectorCovariance(globalVecPoints,high_var=0.05, low_var = 0.01)
	medialVecPoints = medialSpline.getUniformSamples()
	points = addGPACVectorCovariance(medialVecPoints,high_var=0.05, low_var = 0.01)
	
	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	
	globalPoly = []		   
	for p in globalPoints:
		globalPoly.append([p[0],p[1]])	  
	
	#kdInput1 = []
	#for i in range(0,len(points1)):
	#	p_1 = points1[i]
	#	kdInput1.append(p_1[0:2])
	#kdTree1 = cKDTree(array(kdInput1))
	
	while True:
		
		
		
		" find the matching pairs "
		match_pairs = []

		" transform the target Hull with the latest offset "
		poseOrigin = Pose(currPose)
		vecPoints_offset = []
		for p in medialVecPoints:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			vecPoints_offset.append(result)
		
		
		" get the circles and radii "
		match_pairs = []

		for i in range(len(globalPoints)):
			p_1 = globalVecPoints[i]
	
			if True:
		
				" for every point of B, find it's closest neighbor in transformed A "
				try:
					p_2, i_2, minDist = findClosestPointInA(vecPoints_offset, p_1)
		
					if minDist <= minMatchDist:
			
						" add to the list of match pairs less than 1.0 distance apart "
						" keep A points and covariances untransformed "
						C2 = points[i_2][2]
						
						C1 = globalPoints[i][2]
		
						" we store the untransformed point, but the transformed covariance of the A point "
						match_pairs.append([points[i_2],globalPoints[i],C2,C1])
				except:
					pass

			p_1 = globalVecPoints[i]

		for i in range(len(vecPoints_offset)):
			p_1 = vecPoints_offset[i]

			try:
				" for every transformed point of A, find it's closest neighbor in B "
				p_2, i_2, minDist = findClosestPointInA(globalVecPoints, p_1)
	
				if minDist <= minMatchDist2:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "

					C2 = points[i][2]						 
					C1 = globalPoints[i_2][2]	 
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points[i],globalPoints[i_2],C2,C1])
			except:
				pass

		flatMatchPairs = []
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			C1 = pair[2]
			C2 = pair[3]

			flatMatchPairs.append(p1[0])
			flatMatchPairs.append(p1[1])
			flatMatchPairs.append(C1[0][0])
			flatMatchPairs.append(C1[0][1])
			flatMatchPairs.append(C1[1][0])
			flatMatchPairs.append(C1[1][1])
			flatMatchPairs.append(p2[0])
			flatMatchPairs.append(p2[1])
			flatMatchPairs.append(C2[0][0])
			flatMatchPairs.append(C2[0][1])
			flatMatchPairs.append(C2[1][0])
			flatMatchPairs.append(C2[1][1])
		c_poses_1 = [item for sublist in poses_1 for item in sublist]
		c_poses_2 = [item for sublist in poses_2 for item in sublist]
		
		
		
		
		resultSum, resultParam, resultOffset = nelminICP.ICPcost(flatMatchPairs, len(match_pairs), [u1,currU,currAng], uHigh, uLow, c_poses_1, c_poses_2, len(poses_1))
		
		newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), [u1,currU,currAng], uHigh, uLow, c_poses_1, c_poses_2, len(poses_1))
		



		newU = newParam[0]
		newAng = newParam[1]
		
		" set the current parameters "
		currAng = functions.normalizeAngle(newAng)
		currU = newU


		if currU >= 1.0:
			pose2 = poses_2[-1]
		elif currU < 0.0:
			pose2 = poses_2[0]
		else:
			pose2 = poses_2[int(currU*100)]


		point1 = [pose1[0],pose1[1]]
		point2 = [pose2[0],pose2[1]]
		ang1 = pose1[2]
		ang2 = pose2[2]
	
		currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

		# check for convergence condition, different between last and current cost is below threshold
		
		isTerminate = False
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			isTerminate = True


		" update after check for termination condition "
		lastCost = newCost

		numIterations += 1
		
		if isTerminate:
			break
		

	" draw final position "
	if plotIter:
		
		" set the origin of pose 1 "
		poseOrigin = Pose(currPose)
		
		pylab.clf()
		pylab.axes()
		match_global = []
		
		for pair in match_pairs:
			p1 = pair[0]
			p2 = pair[1]
			
			p1_o = dispOffset(p1, currPose)
			
			p1_g = p1_o
			p2_g = p2
			match_global.append([p1_g,p2_g])

		draw_matches(match_global, [0.0,0.0,0.0])

		
		xP = []
		yP = []
		for b in globalPoly:
			p1 = b
			xP.append(p1[0])	
			yP.append(p1[1])

		pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

		pylab.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))

		xP = []
		yP = []
		for b in points:
			p = [b[0],b[1]]
			p1 = dispOffset(p,currPose)
			
			#p1 = poseOrigin.convertLocalToGlobal(p)
			xP.append(p1[0])	
			yP.append(p1[1])

		pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

		#pylab.scatter([pose2[0]],[pose2[1]],color=(1.0,0.0,0.0))


		plotEnv()		 
		pylab.title("%u -- %s u1 = %1.3f, u2 = %1.3f, ang = %1.3f, %.1f, %d F" % (n1, repr(n2), u1, currU, currAng, newCost, len(match_pairs)))

		pylab.xlim(-6, 10)					  
		pylab.ylim(-8, 8)
		pylab.savefig("ICP_plot_%06u.png" % globalPlotCount)
		pylab.clf()
		
		" save inputs "
		#saveFile = ""
		#saveFile += "initGuess = " + repr(initGuess) + "\n"
		#saveFile += "globalPath = " + repr(globalPath) + "\n"
		#saveFile += "medialPoints = " + repr(medialPoints) + "\n"

		#f = open("icpInputSave_%04u.txt" % globalPlotCount, 'w')
		globalPlotCount += 1
		#f.write(saveFile)
		#f.close()		  
			
		#numIterations += 1


	return currPose, newCost, len(match_pairs)


def plotEnv(axes = 0):

	# Y - junction, test 1 				
	"""
	WLEN2 = 7.0
	wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN2*cos(pi/3), -0.2 - WLEN2*sin(pi/3)]]
	wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
	wall5 = [wall2[2],wall1[0]]
	w1 = wall1[-1]
	w2 = wall2[0]

	wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
	wall4 = [w1, [w1[0] + 0.4*cos(pi/3-pi/2), w1[1] - 0.4*sin(pi/3-pi/2)]]

	walls = [wall1, wall2, wall3, wall4, wall5]

	for wall in walls:
		for i in range(len(wall)):
			p = copy(wall[i])
			p[0] += 6.0
			wall[i] = p
	"""

	# cross - junctions, test 6

	wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-5.2], [-3.6,-5.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
	wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,5.2], [-3.6,5.2], [-3.6,0.2], [2.0,0.2]]
	wall2.reverse()
	wall5 = [wall2[-1],wall1[0]]
	walls = [wall1, wall2, wall5]

	for wall in walls:
		for i in range(len(wall)):
			p = copy(wall[i])
			p[0] += 6.0
			wall[i] = p


	if axes == 0:
	
		for wall in walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = copy(wall[i])
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g')

	else:

		for wall in walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = copy(wall[i])
				xP.append(p[0])
				yP.append(p[1])
	
			axes.plot(xP,yP, linewidth=2, color = 'g')
		
	
def draw_matches(match_pairs, offset, axes = 0):


	if axes == 0:
		for pair in match_pairs:
			a_off = dispOffset(pair[0], offset)
			b = pair[1]
	
			pylab.plot([a_off[0],b[0]], [a_off[1],b[1]],color='k')

	else:
		for pair in match_pairs:
			a_off = dispOffset(pair[0], offset)
			b = pair[1]
	
			axes.plot([a_off[0],b[0]], [a_off[1],b[1]],color='k')

def doTest():



	f = open("icpInputSave_input.txt", 'r')		   
	saveStr = f.read()
	f.close()
		
	saveStr = saveStr.replace('\r\n','\n')		  
	exec(saveStr)

	uSet = [i*0.01 for i in range(100)]

	#globalSpline = SplineFit(globalPath, smooth=0.1)
	#medialSpline = SplineFit(medialPoints, smooth=0.1)

	#poses_1 = globalSpline.getUVecSet(uSet)
	#poses_2 = medialSpline.getUVecSet(uSet)

	
	"""
	saveStr = ''
	saveStr += "initGuess = " + repr(initGuess) + "\n"
	saveStr += "globalPath = " + repr(globalPath) + "\n"
	saveStr += "medialPoints = " + repr(medialPoints) + "\n"
	saveStr += "poses_1 = " + repr(poses_1) + "\n"
	saveStr += "poses_2 = " + repr(poses_2) + "\n"

	f = open("cudaData.txt", 'w')
	f.write(saveStr)
	f.close()
	"""


	#import nelminICP
	#param, cost = nelminICP.ICPmin(matchPairs, numPairs, initGuess, poses_1, poses_2, numPoses)

	#uSet = uSet[:20]
	
	args = []
	for i in range(len(uSet)):

		u1 = uSet[i]
		#print "doing u1 =", u1
		initGuess[0] = u1
		initGuess[1] = 0.5
		args.append(copy(initGuess))
		
	results = batchGlobalICP(globalPath, medialPoints, args)
	
	return
	
	print "initGuess:", initGuess
		
	#for u1 in uSet:
	for i in range(len(uSet)/2):

		u1 = uSet[2*i]
		#print "doing u1 =", u1
		initGuess[0] = u1
		initGuess[1] = 0.5
		globalOverlapICP_GPU2(initGuess, globalPath, medialPoints)
		#globalOverlapICP(initGuess, globalPath, medialPoints)
		#globalOverlapICP_GPU(initGuess, globalPath, medialPoints, poses_1, poses_2)

	#print "doing u1 =", u1
	#globalOverlapICP_GPU2(initGuess, globalPath, medialPoints, poses_1, poses_2)


def doTest2():


	f = open("costTest.txt", 'r')		 
	saveStr = f.read()
	f.close()
		
	saveStr = saveStr.replace('\r\n','\n')		  
	exec(saveStr)

	flatMatchPairs = []
	for pair in match_pairs:
		p1 = pair[0]
		p2 = pair[1]
		C1 = pair[2]
		C2 = pair[3]

		flatMatchPairs.append(p1[0])
		flatMatchPairs.append(p1[1])
		flatMatchPairs.append(C1[0][0])
		flatMatchPairs.append(C1[0][1])
		flatMatchPairs.append(C1[1][0])
		flatMatchPairs.append(C1[1][1])
		flatMatchPairs.append(p2[0])
		flatMatchPairs.append(p2[1])
		flatMatchPairs.append(C2[0][0])
		flatMatchPairs.append(C2[0][1])
		flatMatchPairs.append(C2[1][0])
		flatMatchPairs.append(C2[1][1])
	c_poses_1 = [item for sublist in poses_1 for item in sublist]
	c_poses_2 = [item for sublist in poses_2 for item in sublist]

	u1 = initGuess[0]
	currU = initGuess[1]
	currAng = initGuess[2]

	uHigh = currU + 0.2
	uLow = currU - 0.2

	#newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), initGuess, c_poses_1, c_poses_2, len(poses_1))

	for i in range(1):
		cost1, param, offset = nelminICP.ICPcost(flatMatchPairs, len(match_pairs), initGuess, currU+0.2, currU-0.2, c_poses_1, c_poses_2, len(poses_1))
		print "offset =", offset
		print "param =", param
		cost2 = medialOverlapCostFunc([currU, currAng], match_pairs, poses_1, poses_2, uHigh, uLow, u1)			   
		#print cost1, cost2
		print cost1, cost2
		

	#globalOverlapICP_GPU2(initGuess, globalPath, medialPoints)
	#globalOverlapICP(initGuess, globalPath, medialPoints)
	#globalOverlapICP_GPU(initGuess, globalPath, medialPoints, poses_1, poses_2)


class ControlError(Exception):
	def __init__(self, value):
		self.value = value
		
	def __str__(self):
		return repr(self.value)
		

if __name__ == '__main__':


	numpy.set_printoptions(threshold=numpy.nan)

	#print __num_processors()
	
	#exit()

	try:
		time1 = time.time()

		" while loop "		  
		#cProfile.run('doTest()', 'test_prof')
		#doTest2()
		doTest()
		time2 = time.time()
		print time2 - time1
	
	except ControlError as inst:
		print inst.value

	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]
		
	

	exit()

