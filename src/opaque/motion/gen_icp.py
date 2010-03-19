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

import numpy
import scipy 
import scipy.linalg
import scipy.optimize
import math
import pylab
import pca_module
import functions
from subprocess import *

from matplotlib.patches import Circle

from copy import copy
from copy import deepcopy

from coord import Pose

numIterations = 0

fig = pylab.figure()


" displace the point by the offset plus modify it's covariance "
def dispPoint(p, offset):
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
			])

	p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
	temp = T*p_hom
	p_off = [temp[0,0],temp[1,0]]

	R = numpy.matrix([	[math.cos(theta), -math.sin(theta)],
		[math.sin(theta), math.cos(theta)] ])

	Cv = p[2]
	Ca = R * Cv * numpy.transpose(R)
	p_off.append(Ca)

	return p_off

" displace the point by the offset only.  No covariance "
def dispOffset(p, offset):
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
			])

	p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
	temp = T*p_hom
	p_off = [temp[0,0],temp[1,0]]

	return p_off

def dispOffsetMany(points, offset):
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
			])

	results = []

	for p in points:
		p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
		temp = T*p_hom
		p_off = [temp[0,0],temp[1,0]]
		results.append(p_off)

	return results

def disp(ai, bi, T):

	temp = T*ai
	result = bi-temp	

	result[2] = 1.0

	return result


def computeMatchError(offset, a, b, Ca, Cb):

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
			])

	R = numpy.matrix([	[math.cos(theta), -math.sin(theta)],
		[math.sin(theta), math.cos(theta)] ])

	Cv = Ca

	d_vec = disp(numpy.concatenate((a,numpy.matrix([1.0]))), numpy.concatenate((b,numpy.matrix([1.0]))), T)

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

	Cv = numpy.matrix([	[x_var, 0.0],
			[0.0, y_var]
			])

	mag = math.sqrt(vec[0]**2 + vec[1]**2)
	normVec = [vec[0]/mag, vec[1]/mag]

	if normVec[1] == 0:
		R = numpy.matrix([	[1.0, 0.0],
				[0.0, 1.0]
				])

	else:
		B = -1 / (normVec[1] + normVec[0]**2/normVec[1])
		A = -normVec[0]*B/normVec[1]
		R = numpy.matrix([	[A, -B],
				[B, A]
				])

	Ca = numpy.transpose(R) * Cv * R

	return Ca

def findClosestPointInA(a_trans, b):

	#a_p, a_i, minDist = findClosestPointInA(a_trans, b_p, [0.0,0.0,0.0])

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


# for point T*a, find the closest point b in B
def findClosestPointInB(b_data, a, offset):

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
			])


	a_hom = numpy.matrix([[a[0]],[a[1]],[1.0]])
	temp = T*a_hom
	a_off = [temp[0,0],temp[1,0]]

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


def cost_func(offset, match_pairs, a_data_raw = [], polyB = [], circles = []):
	global numIterations
	global fig

	fig.clf()
	ax = fig.add_subplot(111)
	
	if False and len(a_data_raw) > 0 and len(polyB) > 0:

		numIterations += 1
		
		#pylab.clf()
		#pylab.axes()
		draw_matches(match_pairs, offset)
		
		xP = []
		yP = []
		for b in polyB:
			xP.append(b[0])	
			yP.append(b[1])
		
		xP.append(polyB[0][0])	
		yP.append(polyB[0][1])
		
		pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
	
		xP = []
		yP = []
		for b in a_data_raw:
			p = [b[0],b[1]]
			p = dispOffset(p,offset)
			xP.append(p[0])	
			yP.append(p[1])
		
		p = [a_data_raw[0][0],a_data_raw[0][1]]
		p = dispOffset(p,offset)
		xP.append(p[0])	
		yP.append(p[1])
		
		pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
	
		#for circle in circles:
		#	radius, center = circle
		#	c1 = Circle([center[0], center[1]], radius)
		#	ax.add_artist(c1)
		#	c1.set_clip_box(ax.bbox)
		#	c1.set_alpha(0.2)
		#	c1.set_facecolor('r')
		
			#cir = pylab.Circle((center[0],center[1]), radius=radius,  fc='y')
			#pylab.gca().add_patch(cir)
	
		pylab.xlim(-4.5,4.5)
		pylab.ylim(-4,4)
		#pylab.axis('equal')
		pylab.savefig("ICP_plot_%04u.png" % numIterations)
		pylab.clf()		
	
	
	sum = 0.0
	for pair in match_pairs:

		a = numpy.matrix([[pair[0][0]],[pair[0][1]]])
		b = numpy.matrix([[pair[1][0]],[pair[1][1]]])
		sum += computeMatchError(offset, a, b, pair[2], pair[3])

		#polyA_trans = dispOffsetMany(polyA, offset)
		originA = dispOffset([0.0,0.0], [offset[0], offset[1], 0.0])
		originB = dispOffset([0.0,0.0], [-offset[0], -offset[1], 0.0])

		#print polyA_trans

		"""
		if functions.point_inside_polygon(originB[0],originB[1],polyA):
			dist = math.sqrt(originA[0]**2 + originA[1]**2)
			cost = 10000.0 * dist
			sum += cost

		if functions.point_inside_polygon(originA[0],originA[1],polyB):
			dist = math.sqrt(originA[0]**2 + originA[1]**2)
			cost = 10000.0 * dist
			sum += cost
		"""

		# NOTE:  standard point-to-point ICP
		#a = numpy.matrix([[pair[0][0]],[pair[0][1]],[1.0]])
		#b = numpy.matrix([[pair[1][0]],[pair[1][1]],[1.0]])

		#distVec = disp(a, b, T)
		#mag = distVec[0,0]**2 + distVec[1,0]**2
		#sum += mag

		
	return sum


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
			p[2] += C

def addPointToLineCovariance(points, high_var=1.0, low_var=0.001):

	for p in points:
		C = numpy.matrix([	[high_var, 0.0],
				[0.0, high_var]
				])

		try:
			# the covariance matrix that enforces the point-to-plane constraint
			normVec = findLocalNormal(p,points)
			C = computeVectorCovariance(normVec,low_var,high_var)

		except:
			pass

		if len(p) <= 2:
			p.append(C)
		else:
			p[2] += C

" check if a point is contained in the polygon, use radius and center to quicken check "
def isValid(p, radius, center, poly):
	
	#isValid(b_p, radiusA, centerA, polyA):
	
	#print p, radius, center
	dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
	#print "dist =", dist, "radius =", radius, "p =", p, "center =", center
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


def isValidA(p, radiusB, centerB, polyB):

	dist = math.sqrt((p[0] - centerB[0]) ** 2 + (p[1] - centerB[1]) ** 2)
	if dist < radiusB:
		if not functions.point_inside_polygon(p[0],p[1],polyB):
			return True

	return False

def isValidB2(p, radiusA, centerA, polyA):

	dist = math.sqrt((p[0] - centerA[0]) ** 2 + (p[1] - centerA[1]) ** 2)

	if dist < radiusA:
		if not functions.point_inside_polygon(p[0],p[1],polyA):
			return True

	return False


def isValidB(p, offset, radiusA, centerA, polyA):

	newCenter = dispOffset(centerA, offset)

	dist = math.sqrt((p[0] - newCenter[0]) ** 2 + (p[1] - newCenter[1]) ** 2)

	if dist < radiusA:
		if not functions.point_inside_polygon(p[0],p[1],polyA):
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

def gen_ICP_past(estPoses, a_hulls, costThresh = 0.004, minMatchDist = 2.0, plotIter = False):
	
	global numIterations
	#numIterations = 0

	" Determine which local poses are close enough and should be fit together "
	circles = []
	for i in range(len(estPoses)):

		est1 = estPoses[i]

		a_data = a_hulls[i]
		b_data = []
		a_trans = []
		for p in a_data:
			a_trans.append(dispOffset(p, est1))

		radius, center = computeEnclosingCircle(a_trans)
		circles.append([radius, center])

	distances = []
	for j in range(len(circles)-1):
		radius1 = circles[-1][0]
		radius2 = circles[j][0]
		
		center1 = circles[-1][1]
		center2 = circles[j][1]
		
		dist = math.sqrt((center1[0]-center2[0])**2 + (center1[1]-center2[1])**2)
		
		maxDist = radius1 + radius2
		
		if dist < maxDist:
			distances.append(1)
		else:
			distances.append(0)

	print "distances =", distances

	lastCost = 1e100

	circles = []
	a_hull_trans = []
	poly_trans = []
	estPoseOrigin = estPoses[-2]
	for i in range(len(estPoses)):
		estPose2 = estPoses[i]
		poseOrigin = Pose(estPoseOrigin)
		offset = poseOrigin.convertGlobalPoseToLocal(estPose2)

		" transform the past poses "
		a_data = a_hulls[i]
		a_trans = []
		for p in a_data:
			a_trans.append(dispPoint(p, offset))
		
		a_hull_trans.append(a_trans)

		" transformed points without associated covariance "
		polyA_trans = []
		for p in a_trans:
			polyA_trans.append([p[0],p[1]])	
			
		poly_trans.append(polyA_trans)

		" get the circles and radii "
		radius, center = computeEnclosingCircle(a_trans)
		circles.append([radius, center])

	#print len(a_hull_trans)
	#print len(poly_trans)
	#print len(distances)
	#print
	
	if True:
		pylab.clf()
		for i in range(len(a_hull_trans)):
			a_pnts = a_hull_trans[i]
	
			xP = []
			yP = []
			for a in a_pnts:
				xP.append(a[0])	
				yP.append(a[1])
					
			xP.append(a_pnts[0][0])	
			yP.append(a_pnts[0][1])
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
	
		#pylab.show()
		#pylab.savefig("test3.png")
	
	#exit()

	" set the initial guess "
	estPoseOrigin = estPoses[-2]
	estPose2 = estPoses[-1]
	poseOrigin = Pose(estPoseOrigin)
	offset = poseOrigin.convertGlobalPoseToLocal(estPose2)

	while True:
		" find the matching pairs "
		match_pairs = []
		for j in range(0,len(distances)):
			if distances[j] == 1:
				
				if True:
					" transform the past poses "
					a_data_raw = a_hulls[-1]
					a_data = []
					for p in a_data_raw:
						a_data.append(dispPoint(p, offset))
					
					" transformed points without associated covariance "
					polyA = []
					for p in a_data:
						polyA.append([p[0],p[1]])	
											
					radiusA, centerA = computeEnclosingCircle(a_data)

				else:
					a_data = a_hull_trans[-1]
					polyA = poly_trans[-1]
					radiusA = circles[-1][0]
					centerA = circles[-1][1]					
									
				a_data_raw = a_hulls[-1]
				b_data = a_hull_trans[j]			
				polyB = poly_trans[j]
				radiusB = circles[j][0]				
				centerB = circles[j][1]
				
				#print centerA, centerB
				#print len(a_data), len(b_data), len(polyA), len(polyB)
				
				for i in range(len(a_data)):
					a_p = polyA[i]
	
					#if isValidA(a_p, radiusB, centerB, polyB):
					if isValid(a_p, radiusB, centerB, polyB):
	
						" for every transformed point of A, find it's closest neighbor in B "
						b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])
			
						if minDist <= minMatchDist:
				
							" add to the list of match pairs less than 1.0 distance apart "
							" keep A points and covariances untransformed "
							#Ca = a_data[i][2]
							Ca = a_data_raw[i][2]
							Cb = b_p[2]
			
							" we store the untransformed point, but the transformed covariance of the A point "
							#match_pairs.append([a_data[i],b_p,Ca,Cb])
							match_pairs.append([a_data_raw[i],b_p,Ca,Cb])
				
				for i in range(len(b_data)):
					b_p = polyB[i]
			
					#if isValidB2(b_p, radiusA, centerA, polyA):
					if isValid(b_p, radiusA, centerA, polyA):
	
						# for every point of B, find it's closest neighbor in transformed A
						a_p, a_i, minDist = findClosestPointInA(a_data, b_p, [0.0,0.0,0.0])
			
						if minDist <= minMatchDist:
				
							" add to the list of match pairs less than 1.0 distance apart "
							" keep A points and covariances untransformed "
							#Ca = a_data[a_i][2]
							Ca = a_data_raw[a_i][2]
							
							Cb = b_data[i][2]
			
							" we store the untransformed point, but the transformed covariance of the A point "
							#match_pairs.append([a_data[a_i],b_p,Ca,Cb])
							match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])

		#print match_pairs
		if False:
			for pair in match_pairs:
				a = pair[0]
				b = pair[1]
			
				pylab.plot([a[0],b[0]], [a[1],b[1]])
			
			pylab.show()

		# optimize the match error for the current list of match pairs
		newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs])
	
		# get the current cost
		newCost = cost_func(newOffset, match_pairs)
	
		# check for convergence condition, different between last and current cost is below threshold
		if abs(lastCost - newCost) < costThresh:
			offset = newOffset
			lastCost = newCost
			break
	
		numIterations += 1
	
		# optionally draw the position of the points in current transform
		if plotIter:
			pylab.clf()
			draw_matches(match_pairs, offset)
			
			for i in range(len(a_hull_trans)-1):
				a_pnts = a_hull_trans[i]
		
				xP = []
				yP = []
				for a in a_pnts:
					xP.append(a[0])	
					yP.append(a[1])
						
				xP.append(a_pnts[0][0])	
				yP.append(a_pnts[0][1])
				
				pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
			
			draw(a_data,[], "ICP_plot_%04u.png" % numIterations)
	
		# save the current offset and cost
		offset = newOffset
		lastCost = newCost

	return offset

#def gen_ICP_global(estPoses, a_hulls, costThresh = 0.004, minMatchDist = 2.0, plotIter = False):
def gen_ICP_global(pastPose, targetPose, pastHull, targetHull, pastCircles, costThresh = 0.004, minMatchDist = 2.0, plotIter = False):

	global numIterations
	
	stepOne = True
	lastCost = 1e100
	
	startIteration = numIterations

	estPoseOrigin = pastPose
	
	#print "estPoseOrigin =", estPoseOrigin
	
	" set the initial guess "
	estPose2 = targetPose
	#print "estPose2 =", estPose2
	
	poseOrigin = Pose(estPoseOrigin)
	
	offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
	#print "offset =", offset

	#exit()
	" transform the past poses "
	a_data_raw = targetHull
	a_data = []
	for p in a_data_raw:
		result = dispPoint(p, offset)
		#print offset, p, result
		#print
		
		a_data.append(result)
	
	#print "targetHull =", targetHull
	#print "a_trans =", a_trans
	

	
	polyB = []		
	for p in pastHull:
		polyB.append([p[0],p[1]])	

	while True:
		" find the matching pairs "
		match_pairs = []

		a_data_raw = targetHull
		b_data = pastHull			

		" transform the target Hull with the latest offset "
		a_data = []
		for p in a_data_raw:
			result = dispPoint(p, offset)
			a_data.append(result)

		" transformed points without associated covariance "
		polyA = []
		for p in a_data:
			polyA.append([p[0],p[1]])	
		
		" get the circles and radii "
		radiusA, centerA = computeEnclosingCircle(a_data)
		#radiusB, centerB = computeEnclosingCircle(pastHull)
		
		if True:
			for i in range(len(a_data)):
				a_p = polyA[i]
	
				#if isValid(a_p, radiusB, centerB, polyB):
				if isValidPast(a_p, pastCircles, polyB):
	
					" for every transformed point of A, find it's closest neighbor in B "
					b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])
		
					if minDist <= minMatchDist:
			
						" add to the list of match pairs less than 1.0 distance apart "
						" keep A points and covariances untransformed "
						Ca = a_data_raw[i][2]
						Cb = b_p[2]
		
						" we store the untransformed point, but the transformed covariance of the A point "
						match_pairs.append([a_data_raw[i],b_p,Ca,Cb])

		if True:
			for i in range(len(b_data)):
				b_p = polyB[i]
		
				if isValid(b_p, radiusA, centerA, polyA):
			
					#print "selected", b_p, "in circle", radiusA, centerA, "with distance"
					" for every point of B, find it's closest neighbor in transformed A "
					a_p, a_i, minDist = findClosestPointInA(a_data, b_p)
		
					if minDist <= minMatchDist:
			
						" add to the list of match pairs less than 1.0 distance apart "
						" keep A points and covariances untransformed "
						Ca = a_data_raw[a_i][2]
						
						Cb = b_data[i][2]
		
						" we store the untransformed point, but the transformed covariance of the A point "
						match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])


		" plot the initial configuration "
		if startIteration == numIterations:

			numIterations += 1

			pylab.clf()
			pylab.axes()
			match_global = []
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_o = dispOffset(p1, offset)
				#p2_o = dispOffset(p2, offset)

				
				p1_g = poseOrigin.convertLocalToGlobal(p1_o)
				p2_g = poseOrigin.convertLocalToGlobal(p2)
				match_global.append([p1_g,p2_g])
			
			draw_matches(match_global, [0.0,0.0,0.0])
			
			xP = []
			yP = []
			for b in polyB:
				p1 = poseOrigin.convertLocalToGlobal(b)

				xP.append(p1[0])	
				yP.append(p1[1])
			
			p1 = poseOrigin.convertLocalToGlobal(polyB[0])
			xP.append(p1[0])	
			yP.append(p1[1])
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			xP = []
			yP = []
			for b in a_data_raw:
				p = [b[0],b[1]]
				p = dispOffset(p,offset)
				
				p1 = poseOrigin.convertLocalToGlobal(p)
				xP.append(p1[0])	
				yP.append(p1[1])
			
			p = [a_data_raw[0][0],a_data_raw[0][1]]
			p = dispOffset(p,offset)

			p1 = poseOrigin.convertLocalToGlobal(p)
			xP.append(p1[0])	
			yP.append(p1[1])

			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			plotEnv()		
			
			pylab.xlim(-4,9)
			pylab.ylim(-7,5)
			pylab.savefig("ICP_plot_%04u.png" % numIterations)
			pylab.clf()					
			
		if False:
			for pair in match_pairs:
				a = pair[0]
				b = pair[1]
			
				pylab.plot([a[0],b[0]], [a[1],b[1]],color=(1.0,0.0,0.0))
			
			pylab.show()
			
		#allCircles = deepcopy(pastCircles)
		#allCircles.append([radiusA,centerA])

		allCircles = [[radiusA,centerA]]

		# optimize the match error for the current list of match pairs
		newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs, a_data_raw, polyB, allCircles])
	
		# get the current cost
		newCost = cost_func(newOffset, match_pairs)
	
		# check for convergence condition, different between last and current cost is below threshold
		if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
			offset = newOffset
			lastCost = newCost
			break
	
		numIterations += 1

		offset = newOffset

		" reduce the minMatch distance for each step down to a floor value "
		"""
		if stepOne:
			stepOne = False
		else:
			minMatchDist /= 2
			if minMatchDist < 0.25:
				minMatchDist = 0.25
		"""
		minMatchDist /= 2
		if minMatchDist < 0.25:
			minMatchDist = 0.25
	
		# optionally draw the position of the points in current transform
		if plotIter:
			pylab.clf()
			pylab.axes()
			match_global = []
			
			#match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])
			
			for pair in match_pairs:
				p1 = pair[0]
				p2 = pair[1]
				
				p1_o = dispOffset(p1, offset)
				#p2_o = dispOffset(p2, offset)

				
				p1_g = poseOrigin.convertLocalToGlobal(p1_o)
				p2_g = poseOrigin.convertLocalToGlobal(p2)
				match_global.append([p1_g,p2_g])
			
			#draw_matches(match_pairs, offset)
			draw_matches(match_global, [0.0,0.0,0.0])
			
			xP = []
			yP = []
			for b in polyB:
				p1 = poseOrigin.convertLocalToGlobal(b)

				xP.append(p1[0])	
				yP.append(p1[1])
				#xP.append(b[0])	
				#yP.append(b[1])
			
			p1 = poseOrigin.convertLocalToGlobal(polyB[0])
			xP.append(p1[0])	
			yP.append(p1[1])
			#xP.append(polyB[0][0])	
			#yP.append(polyB[0][1])
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			xP = []
			yP = []
			for b in a_data_raw:
				p = [b[0],b[1]]
				p = dispOffset(p,offset)
				
				p1 = poseOrigin.convertLocalToGlobal(p)
				xP.append(p1[0])	
				yP.append(p1[1])
				
				#xP.append(p[0])	
				#yP.append(p[1])
			
			p = [a_data_raw[0][0],a_data_raw[0][1]]
			p = dispOffset(p,offset)

			p1 = poseOrigin.convertLocalToGlobal(p)
			xP.append(p1[0])	
			yP.append(p1[1])

			#xP.append(p[0])	
			#yP.append(p[1])
			
			pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

			#cir = Circle( (0,0), radius=0.5)
			#a.add_patch(cir)

			plotEnv()		
			
			#for circle in pastCircles:
			#	radius, center = circle
			#	cir = pylab.Circle((center[0],center[1]), radius=radius,  fc='r')
			#	pylab.gca().add_patch(cir)
				
			#cir = pylab.Circle((centerA[0],centerA[1]), radius=radiusA,  fc='r')
			#pylab.gca().add_patch(cir)
			#cir = pylab.Circle((.5,.5), radius=0.25, alpha =.2, fc='b')
			#pylab.gca().add_patch(cir)

			#pylab.xlim(-4.5,4.5)
			#pylab.ylim(-4,4)
			#pylab.xlim(-3,6)
			#pylab.ylim(-4,4)
			#pylab.xlim(-4,7)
			#pylab.ylim(-7,3)
			pylab.xlim(-4,9)
			pylab.ylim(-7,5)
			#pylab.axis('equal')
			pylab.savefig("ICP_plot_%04u.png" % numIterations)
			pylab.clf()			
						
		# save the current offset and cost
		offset = newOffset
		lastCost = newCost
	
	


	return offset


def gen_ICP(offset, a_data, b_data, costThresh = 0.004, minMatchDist = 2.0, plotIter = False):

	# 1. compute the local line for each point
	# 2. compute the covariance for point-to-line constraint
	# 3. find closest points in A of B
	# 4. discard matches beyond threshold d_max < | ai - T *bi |
	# 5. optimize T for the sum of computeMatchError
	# 6. if converged, stop, else go to 3 with b offset by new T

	# DONE:  store the covariance for each point untransformed so we don't have to repeat
	# DONE:  transform each covariance by the appropriate rotation at each point

	# compute the local covariance matrix for a point-to-plane contrainst

	numIterations = 0
	lastCost = 1e100



	radiusA, centerA = computeEnclosingCircle(a_data)
	radiusB, centerB = computeEnclosingCircle(b_data)
	

	print "radiusA =", radiusA, "center =", centerA
	print "radiusB =", radiusB, "center =", centerB

	#radiusA = 2.0
	#radiusB = 2.0


	polyA = []
	polyB = []

	for p in a_data:
		polyA.append([p[0],p[1]])	

	for p in b_data:
		polyB.append([p[0],p[1]])	

	while True:

		a_trans = []

		# pre-transform the A points and their associated covariances
		for p in a_data:
			a_trans.append(dispPoint(p, offset))

		polyA_trans = []
		for p in a_trans:
			polyA_trans.append([p[0],p[1]])	

		match_pairs = []
	
		for i in range(len(a_trans)):
			a_p = a_trans[i]

			#radiusB = 2.5
			#if isValidA(a_p, offset, b_data):
			if isValidA(a_p, radiusB, centerB, polyB):

				# for every transformed point of A, find it's closest neighbor in B
				b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])

				if minDist <= minMatchDist:
		
					# add to the list of match pairs less than 1.0 distance apart 
					# keep A points and covariances untransformed
					Ca = a_data[i][2]
					
					Cb = b_p[2]

					# we store the untransformed point, but the transformed covariance of the A point
					match_pairs.append([a_data[i],b_p,Ca,Cb])


		for i in range(len(b_data)):
			b_p = b_data[i]

			#radiusA = 2.5
			if isValidB(b_p, offset, radiusA, centerA, polyA_trans):

				# for every point of B, find it's closest neighbor in transformed A
				a_p, a_i, minDist = findClosestPointInA(a_trans, b_p, [0.0,0.0,0.0])

				if minDist <= minMatchDist:
		
					# add to the list of match pairs less than 1.0 distance apart 
					# keep A points and covariances untransformed
					Ca = a_data[a_i][2]
					
					Cb = b_p[2]

					# we store the untransformed point, but the transformed covariance of the A point
					match_pairs.append([a_data[a_i],b_p,Ca,Cb])


		
		# optimize the match error for the current list of match pairs
		newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs])

		# get the current cost
		newCost = cost_func(newOffset, match_pairs)

		# check for convergence condition, different between last and current cost is below threshold
		if abs(lastCost - newCost) < costThresh:
			offset = newOffset
			lastCost = newCost
			break

		numIterations += 1

		# optionally draw the position of the points in current transform
		if plotIter:
			draw_matches(match_pairs, offset)
			draw(a_trans,b_data, "ICP_plot_%04u.png" % numIterations)

		# save the current offset and cost
		offset = newOffset
		lastCost = newCost

	return offset

def plotEnv():

	WLEN = 3.0
	WLEN2 = 5.0
	wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*math.cos(math.pi/3), -0.2 - WLEN*math.sin(math.pi/3)]]
	wall2 = [[-4.0 + WLEN2*math.cos(math.pi/3), 0.2 + WLEN2*math.sin(math.pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
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
	
def draw_matches(match_pairs, offset):

	for pair in match_pairs:
		a_off = dispOffset(pair[0], offset)
		b = pair[1]

		pylab.plot([a_off[0],b[0]], [a_off[1],b[1]])


def draw(a_pnts, b_pnts, filename, fileWrite = True):

	xP = []
	yP = []
	for a in a_pnts:
		xP.append(a[0])	
		yP.append(a[1])	
		#pylab.scatter([a[0]],[a[1]],linewidth=1, color=(0.6,0.6,1.0), edgecolors='none')
	if len(xP) > 0:
		pylab.plot(xP,yP, linewidth=1, color='b')
		pass
		
	xP = []
	yP = []
	for b in b_pnts:
		xP.append(b[0])	
		yP.append(b[1])	
		#pylab.scatter([b[0]],[b[1]],linewidth=1, color=(1.0,0.6,0.6), edgecolors='none')
	if len(xP) > 0:
		pylab.plot(xP,yP, linewidth=1, color='r')
		pass

	if fileWrite:
		pylab.xlim(-4.5,4.5)
		pylab.ylim(-4,4)
		pylab.savefig(filename)
		pylab.clf()

def save(filename):
	pylab.savefig(filename)
	pylab.clf()


def computeUnions(point_sets):

	currPoly = point_sets[0]

	for i in range(1,len(point_sets)):
		currPoly = computeUnion(currPoly,point_sets[i])

	return currPoly

def computeUnion(points1, points2):
	
	try:
		
		" remove any degenerate segments "
		nPoints1 = []
		nPoints2 = []
		
		for i in range(len(points1)-1):
			if not (points1[i][0] == points1[i+1][0] and points1[i][1] == points1[i+1][1]):
				nPoints1.append(points1[i])
		nPoints1.append(points1[-1])

		for i in range(len(points2)-1):
			if not (points2[i][0] == points2[i+1][0] and points2[i][1] == points2[i+1][1]):
				nPoints2.append(points2[i])
		nPoints2.append(points2[-1])

		numPoints1 = len(nPoints1)
		numPoints2 = len(nPoints2)
	
		#print numPoints1, numPoints2
	
		inputStr = str(numPoints1) + " " + str(numPoints2) + " "
		
		" Removed Gaussians because were causing self-intersections ."
		" not necessary anymore because polygons come from CGAL and are not degenerate. "
		for p in nPoints1:
			p2 = copy(p)
			#p2[0] += gauss(0.0,0.0001)
			#p2[1] += gauss(0.0,0.0001)
	
			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
	
		for p in nPoints2:
			p2 = copy(p)
			#p2[0] += gauss(0.0,0.0001)
			#p2[1] += gauss(0.0,0.0001)
	
			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		
		inputStr += "\n"
		
		#print inputStr 
		
		f = open("input2.txt", 'w')
		f.write(inputStr)
		f.close()
		
		" start the subprocess "
		subProc = Popen(["./poly_union.exe"], stdin=PIPE, stdout=PIPE)
		
		" send input and receive output "
		sout, serr = subProc.communicate(inputStr)
	
	
		
		" convert string output to typed data "
		sArr = sout.split(" ")
	
		numVert = int(sArr[0])
		
		#print "numVert =", numVert
		#print numVert, len(sArr)
		#print sArr
		
		vertices = []
		for i in range(numVert):
			#print sArr[2*i+1], sArr[2*i+2]
			vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
		
		return vertices
	except:
		print "polygon union failed with the following data:"
		print "sending input to poly_union:  ", inputStr
		print "rawOutput ="
		print sout
	
		print serr		
		raise


if __name__ == '__main__':



	inputStr = "323 175 0.782056 0.654747 0.746317 0.706635 0.702085 0.773295 0.674216 0.811786 0.636855 0.862563 0.577827 0.84605 0.543567 0.802175 0.49582 0.737985 0.457926 0.684164 0.429763 0.642696 0.384887 0.576467 0.356848 0.534852 0.328793 0.493153 0.283956 0.426898 0.256036 0.385406 0.22816 0.343921 0.183368 0.277636 0.153638 0.237144 0.105648 0.173136 0.0677733 0.119273 0.0397852 0.0777528 0.0263176 0.0578137 0.0135987 0.0611008 -0.0295399 0.0357618 -0.0598928 -0.0111072 -0.0905972 -0.0579506 -0.121016 -0.104847 -0.151364 -0.151644 -0.181771 -0.198579 -0.212376 -0.245375 -0.24257 -0.292198 -0.273265 -0.339256 -0.30371 -0.386283 -0.334194 -0.433067 -0.364533 -0.479823 -0.39506 -0.526632 -0.425435 -0.5737 -0.455891 -0.620513 -0.486375 -0.667436 -0.516476 -0.706214 -0.532903 -0.686855 -0.57036 -0.670711 -0.613448 -0.738116 -0.64888 -0.793543 -0.692269 -0.860755 -0.728033 -0.916154 -0.768024 -0.985442 -0.793122 -1.02881 -0.818296 -1.07228 -0.858302 -1.14156 -0.884855 -1.18384 -0.928024 -1.25119 -0.963537 -1.30657 -1.00673 -1.3739 -1.04062 -1.43014 -1.06582 -1.47371 -1.1058 -1.543 -1.13087 -1.58638 -1.15599 -1.62983 -1.1575 -1.63247 -1.17686 -1.65925 -1.20714 -1.70648 -1.24877 -1.7748 -1.28131 -1.83188 -1.29765 -1.86184 -1.29955 -1.86499 -1.31375 -1.88969 -1.34695 -1.94408 -1.38121 -2.0002 -1.39092 -2.01617 -1.41818 -2.05639 -1.44636 -2.09791 -1.49108 -2.16424 -1.51926 -2.20581 -1.54738 -2.24723 -1.59207 -2.31359 -1.62019 -2.35527 -1.64828 -2.39689 -1.69321 -2.46308 -1.72261 -2.50339 -1.77028 -2.56764 -1.80805 -2.62132 -1.83637 -2.66299 -1.88125 -2.72922 -1.90929 -2.77077 -1.9373 -2.8123 -1.99986 -2.86216 -2.06242 -2.91201 -2.12499 -2.96187 -2.18755 -3.01173 -2.25011 -3.06159 -2.31267 -3.11145 -2.37524 -3.16131 -2.42599 -3.20175 -2.42762 -3.21243 -2.43229 -3.21474 -2.504 -3.25021 -2.57571 -3.28568 -2.64742 -3.32115 -2.71913 -3.35661 -2.79083 -3.39208 -2.86254 -3.42755 -2.87831 -3.43535 -2.9253 -3.43372 -2.98099 -3.43349 -3.03695 -3.43319 -3.09278 -3.42989 -3.17252 -3.42347 -3.23807 -3.42122 -3.29407 -3.42095 -3.34996 -3.42083 -3.40583 -3.42057 -3.46168 -3.42026 -3.5175 -3.42003 -3.57351 -3.41985 -3.62938 -3.4196 -3.68545 -3.41938 -3.74128 -3.41917 -3.79723 -3.41892 -3.85313 -3.41867 -3.89148 -3.4297 -3.89873 -3.43928 -3.95947 -3.43567 -4.0153 -3.4369 -4.09516 -3.43229 -4.1608 -3.42851 -4.21679 -3.42965 -4.27264 -3.43075 -4.35251 -3.4262 -4.41826 -3.42246 -4.47385 -3.42332 -4.52982 -3.42437 -4.60969 -3.41979 -4.67559 -3.41602 -4.73121 -3.41707 -4.78723 -3.4178 -4.84316 -3.41922 -4.88737 -3.44257 -4.88737 -3.44257 -4.89015 -3.44238 -4.97015 -3.44187 -5.01499 -3.44168 -5.09499 -3.44136 -5.13999 -3.44114 -5.21999 -3.44075 -5.26516 -3.44279 -5.3152 -3.44537 -5.3952 -3.44492 -5.44008 -3.44467 -5.52008 -3.44426 -5.56516 -3.44403 -5.64516 -3.44366 -5.68997 -3.44344 -5.76997 -3.44302 -5.81526 -3.44278 -5.89516 -3.4468 -5.94495 -3.44927 -5.99548 -3.45177 -6.03537 -3.48146 -6.07342 -3.53214 -6.09594 -3.58336 -6.12648 -3.6573 -6.14148 -3.71626 -6.15172 -3.78604 -6.14205 -3.84119 -6.06205 -3.84151 -6.01715 -3.84169 -5.93715 -3.8421 -5.89206 -3.84232 -5.81207 -3.84273 -5.76699 -3.84295 -5.68785 -3.85468 -5.60871 -3.8664 -5.52958 -3.87812 -5.45234 -3.88956 -5.39736 -3.87987 -5.34216 -3.87004 -5.28712 -3.86027 -5.23209 -3.85062 -5.15219 -3.84654 -5.0723 -3.84247 -5.05231 -3.84145 -4.97241 -3.8376 -4.92197 -3.83724 -4.84197 -3.83753 -4.7969 -3.83777 -4.71691 -3.83823 -4.67183 -3.83843 -4.59183 -3.83879 -4.54711 -3.839 -4.46711 -3.83937 -4.43143 -3.83379 -4.39914 -3.83567 -4.34327 -3.83463 -4.28733 -3.83347 -4.20747 -3.8381 -4.14184 -3.84191 -4.08594 -3.84088 -4.03013 -3.83987 -3.97418 -3.83873 -3.89431 -3.84334 -3.82854 -3.84715 -3.77267 -3.84609 -3.71672 -3.84493 -3.66095 -3.84413 -3.58108 -3.84874 -3.51534 -3.85253 -3.45951 -3.85129 -3.4035 -3.85035 -3.34754 -3.84918 -3.29182 -3.84822 -3.23607 -3.8471 -3.15609 -3.84869 -3.07611 -3.85029 -2.99612 -3.85188 -2.91614 -3.85347 -2.83615 -3.85506 -2.77004 -3.85638 -2.7588 -3.85766 -2.67893 -3.86218 -2.6287 -3.86503 -2.57858 -3.86787 -2.49871 -3.87233 -2.44876 -3.8702 -2.39361 -3.86089 -2.31765 -3.83972 -2.31481 -3.83873 -2.31037 -3.83875 -2.25433 -3.83906 -2.19845 -3.83598 -2.13808 -3.83116 -2.10701 -3.83127 -2.06513 -3.81597 -2.06033 -3.81398 -1.99629 -3.83316 -1.94736 -3.84302 -1.86817 -3.83793 -1.82208 -3.81877 -1.76936 -3.7586 -1.74054 -3.69552 -1.72471 -3.6171 -1.70888 -3.53868 -1.69306 -3.46026 -1.67723 -3.38184 -1.66141 -3.30343 -1.64779 -3.23598 -1.63821 -3.21955 -1.59622 -3.15146 -1.57041 -3.10316 -1.54224 -3.05478 -1.50037 -2.98662 -1.49527 -2.97711 -1.48681 -2.9639 -1.44506 -2.89566 -1.41077 -2.83947 -1.36915 -2.77115 -1.33485 -2.7148 -1.29325 -2.64646 -1.2565 -2.5919 -1.22808 -2.54401 -1.20403 -2.50012 -1.16229 -2.43188 -1.12803 -2.37586 -1.08638 -2.30755 -1.05204 -2.25122 -1.01035 -2.18295 -0.976081 -2.12683 -0.943624 -2.08139 -0.905314 -2.01115 -0.881258 -1.96736 -0.864868 -1.93757 -0.852933 -1.91983 -0.808275 -1.85346 -0.763617 -1.78708 -0.718959 -1.72071 -0.674301 -1.65433 -0.645067 -1.60963 -0.606052 -1.54938 -0.575654 -1.50266 -0.545161 -1.45578 -0.514778 -1.40875 -0.484331 -1.36195 -0.435266 -1.29876 -0.401757 -1.2556 -0.371163 -1.20875 -0.345223 -1.16836 -0.321475 -1.15082 -0.276606 -1.08458 -0.248607 -1.04315 -0.220534 -1.00159 -0.175639 -0.935372 -0.1459 -0.894693 -0.0980938 -0.830549 -0.0603663 -0.776854 -0.032353 -0.735423 0.0123763 -0.669096 0.0404848 -0.627549 0.0686508 -0.585952 0.113469 -0.519685 0.141469 -0.478052 0.169376 -0.436498 0.214248 -0.370267 0.24237 -0.328761 0.270532 -0.287198 0.315167 -0.220807 0.343199 -0.179228 0.371258 -0.13764 0.419252 -0.0736363 0.458505 -0.0212898 0.503151 0.0450937 0.531337 0.086775 0.559558 0.128448 0.604252 0.194799 0.632253 0.236319 0.660345 0.277961 0.705247 0.344171 0.731318 0.382803 0.734514 0.385505 0.734732 0.387866 0.761258 0.427221 0.806076 0.493488 0.825424 0.534521 0.829394 0.590255 -3.84460529197 -3.82900135089 -3.76460540964 -3.82913856068 -3.71495033922 -3.83117007124 -3.63505189559 -3.83519979889 -3.56902055147 -3.83658476609 -3.51892676239 -3.83667500634 -3.43892697064 -3.83685754034 -3.38882622381 -3.83703303216 -3.33869715289 -3.83722408104 -3.2586973238 -3.83738944761 -3.20839117506 -3.8375155713 -3.1584067702 -3.8376465555 -3.0784070987 -3.8378758163 -3.02813688751 -3.8379083647 -2.97814373788 -3.83791225956 -2.89814407845 -3.83814569522 -2.84786524893 -3.83831920567 -2.79809505061 -3.83849777794 -2.71809539593 -3.83873283143 -2.66900492785 -3.84935510983 -2.59186031575 -3.87053785809 -2.51471570365 -3.89172060636 -2.43757109155 -3.91290335462 -2.36042647945 -3.93408610288 -2.28328186736 -3.95526885115 -2.21592011365 -3.95770077147 -2.16307593811 -3.91521348968 -2.12708025114 -3.84376903209 -2.09108456417 -3.77232457451 -2.07815332484 -3.73292888855 -2.05362342379 -3.65678242471 -2.03460052889 -3.59773097665 -2.00748227285 -3.5224674384 -1.99950994572 -3.46065816842 -1.99936675208 -3.41048352814 -2.02696047509 -3.36191822725 -2.05460686723 -3.31334170307 -2.0962044665 -3.28529373996 -2.13762567848 -3.25770106524 -2.19308664517 -3.25055352901 -2.26810166159 -3.2783518619 -2.34311667802 -3.30615019478 -2.41813169444 -3.33394852767 -2.49314671086 -3.36174686055 -2.56816172728 -3.38954519344 -2.6431767437 -3.41734352633 -2.71354302609 -3.44341917671 -2.79354244643 -3.44311463518 -2.84375602929 -3.4429899247 -2.89370887589 -3.44288273927 -2.97370853605 -3.44264955437 -3.02399973154 -3.44253246532 -3.07419109443 -3.44242318371 -3.15419098987 -3.4422938382 -3.20429594883 -3.44225834405 -3.25444544011 -3.44223432718 -3.33444506631 -3.44198976963 -3.38460735199 -3.44179940065 -3.43478855318 -3.44159954985 -3.51478852488 -3.44153225349 -3.56489678094 -3.44142140448 -3.61505675509 -3.44129306251 -3.69505639907 -3.44105439283 -3.74513095058 -3.43904044074 -3.82503250852 -3.43507294278 -3.88545835061 -3.44071462704 -3.95486234214 -3.45437801521 -4.03479214956 -3.45102750987 -4.11472195697 -3.44767700452 -4.19465176439 -3.44432649918 -4.27458157181 -3.44097599383 -4.35451137923 -3.43762548849 -4.43444118665 -3.43427498315 -4.4818553633 -3.43228747114 -4.56185471357 -3.43196504842 -4.60654200111 -3.43883254971 -4.68652435287 -3.43715224976 -4.76650670462 -3.4354719498 -4.84648905638 -3.43379164985 -4.92647140813 -3.43211134989 -5.00645375988 -3.43043104994 -5.08643611164 -3.42875074998 -5.16641846339 -3.42707045003 -5.24640081514 -3.42539015007 -5.3263831669 -3.42370985012 -5.38853899809 -3.42309471871 -5.43857377542 -3.42290750312 -5.51857374387 -3.4228364482 -5.56855078121 -3.42254646646 -5.61889810939 -3.42219261746 -5.69889787103 -3.4223879059 -5.74904069333 -3.42224832631 -5.79906185715 -3.42204281239 -5.87906120845 -3.42172064468 -5.92933357241 -3.42157766291 -5.97934683101 -3.42145061243 -6.05925094485 -3.41753492705 -6.12481820616 -3.41432181631 -6.20481719041 -3.41391867807 -6.2553648537 -3.41377550739 -6.30526332591 -3.41366321382 -6.38526331167 -3.41361548099 -6.43545380999 -3.41356517604 -6.48533692631 -3.41351002441 -6.56533689772 -3.41344238596 -6.6154585591 -3.41135902962 -6.6953508716 -3.40720952335 -6.7613936331 -3.40577722896 -6.81134372749 -3.40567775068 -6.89134351581 -3.40549371962 -6.94158754726 -3.40537322491 -6.99155977241 -3.40525213044 -7.07155968949 -3.40513694929 -7.12179525101 -3.40496005024 -7.17193373625 -3.40475678492 -7.24140803634 -3.41858973932 -7.26930020062 -3.45999080177 -7.29682111845 -3.50155641234 -7.32636971025 -3.55928081314 -7.35699748917 -3.63318574638 -7.38762526808 -3.70709067963 -7.40861204281 -3.75773183745 -7.36169827837 -3.79448404428 -7.28180082377 -3.79853333409 -7.21620573024 -3.79992253968 -7.1659977426 -3.80003803477 -7.08599802119 -3.80024915972 -7.03605066192 -3.80031331068 -6.98584523617 -3.80036087924 -6.90584606208 -3.80072439656 -6.85565374297 -3.80275372855 -6.77575209027 -3.80671931756 -6.70998286143 -3.80815431498 -6.65988319544 -3.8083497464 -6.57997688387 -3.81222032408 -6.51434672589 -3.8153993799 -6.4344547492 -3.81124341307 -6.38945352434 -3.80890245713 -6.30945370867 -3.80907419415 -6.2592619437 -3.81108788911 -6.17936117792 -3.81507130831 -6.11832188645 -3.83480137123 -6.04734670986 -3.87171375815 -5.97637153327 -3.90862614506 -5.90539635668 -3.94553853198 -5.83442118009 -3.9824509189 -5.7634460035 -4.01936330582 -5.6924708269 -4.05627569273 -5.62149565031 -4.09318807965 -5.55116256828 -4.1297665295 -5.47488661224 -4.15089592967 -5.42647739228 -4.12327464324 -5.35842361949 -4.0812189182 -5.2903698467 -4.03916319316 -5.2223160739 -3.99710746812 -5.15426230111 -3.95505174308 -5.08620852832 -3.91299601804 -5.01815475553 -3.87094029299 -4.95010098274 -3.82888456795 -4.91247492605 -3.80563249787 -4.84311373223 -3.79183448947 -4.76843992134 -3.78557932998 -4.68845043522 -3.78687629199 -4.6084609491 -3.788173254 -4.52847146298 -3.78947021601 -4.44848197686 -3.79076717802 -4.36849249075 -3.79206414003 -4.30242177961 -3.79313542085 -4.22713693455 -3.82019446933 -4.16544060725 -3.8282156323 -4.11519277984 -3.82849551347 -4.03519296498 -3.82866762414 -3.9850251929 -3.82869939508 -3.93478107715 -3.82871182562 -3.85478148949 -3.82896867842"
	points = [0.782056, 0.654747, 0.746317, 0.706635, 0.702085, 0.773295, 0.674216, 0.811786, 0.636855, 0.862563, 0.577827, 0.84605, 0.543567, 0.802175, 0.49582, 0.737985, 0.457926, 0.684164, 0.429763, 0.642696, 0.384887, 0.576467, 0.356848, 0.534852, 0.328793, 0.493153, 0.283956, 0.426898, 0.256036, 0.385406, 0.22816, 0.343921, 0.183368, 0.277636, 0.153638, 0.237144, 0.105648, 0.173136, 0.0677733, 0.119273, 0.0397852, 0.0777528, 0.0263176, 0.0578137, 0.0135987, 0.0611008, -0.0295399, 0.0357618, -0.0598928, -0.0111072, -0.0905972, -0.0579506, -0.121016, -0.104847, -0.151364, -0.151644, -0.181771, -0.198579, -0.212376, -0.245375, -0.24257, -0.292198, -0.273265, -0.339256, -0.30371, -0.386283, -0.334194, -0.433067, -0.364533, -0.479823, -0.39506, -0.526632, -0.425435, -0.5737, -0.455891, -0.620513, -0.486375, -0.667436, -0.516476, -0.706214, -0.532903, -0.686855, -0.57036, -0.670711, -0.613448, -0.738116, -0.64888, -0.793543, -0.692269, -0.860755, -0.728033, -0.916154, -0.768024, -0.985442, -0.793122, -1.02881, -0.818296, -1.07228, -0.858302, -1.14156, -0.884855, -1.18384, -0.928024, -1.25119, -0.963537, -1.30657, -1.00673, -1.3739, -1.04062, -1.43014, -1.06582, -1.47371, -1.1058, -1.543, -1.13087, -1.58638, -1.15599, -1.62983, -1.1575, -1.63247, -1.17686, -1.65925, -1.20714, -1.70648, -1.24877, -1.7748, -1.28131, -1.83188, -1.29765, -1.86184, -1.29955, -1.86499, -1.31375, -1.88969, -1.34695, -1.94408, -1.38121, -2.0002, -1.39092, -2.01617, -1.41818, -2.05639, -1.44636, -2.09791, -1.49108, -2.16424, -1.51926, -2.20581, -1.54738, -2.24723, -1.59207, -2.31359, -1.62019, -2.35527, -1.64828, -2.39689, -1.69321, -2.46308, -1.72261, -2.50339, -1.77028, -2.56764, -1.80805, -2.62132, -1.83637, -2.66299, -1.88125, -2.72922, -1.90929, -2.77077, -1.9373, -2.8123, -1.99986, -2.86216, -2.06242, -2.91201, -2.12499, -2.96187, -2.18755, -3.01173, -2.25011, -3.06159, -2.31267, -3.11145, -2.37524, -3.16131, -2.42599, -3.20175, -2.42762, -3.21243, -2.43229, -3.21474, -2.504, -3.25021, -2.57571, -3.28568, -2.64742, -3.32115, -2.71913, -3.35661, -2.79083, -3.39208, -2.86254, -3.42755, -2.87831, -3.43535, -2.9253, -3.43372, -2.98099, -3.43349, -3.03695, -3.43319, -3.09278, -3.42989, -3.17252, -3.42347, -3.23807, -3.42122, -3.29407, -3.42095, -3.34996, -3.42083, -3.40583, -3.42057, -3.46168, -3.42026, -3.5175, -3.42003, -3.57351, -3.41985, -3.62938, -3.4196, -3.68545, -3.41938, -3.74128, -3.41917, -3.79723, -3.41892, -3.85313, -3.41867, -3.89148, -3.4297, -3.89873, -3.43928, -3.95947, -3.43567, -4.0153, -3.4369, -4.09516, -3.43229, -4.1608, -3.42851, -4.21679, -3.42965, -4.27264, -3.43075, -4.35251, -3.4262, -4.41826, -3.42246, -4.47385, -3.42332, -4.52982, -3.42437, -4.60969, -3.41979, -4.67559, -3.41602, -4.73121, -3.41707, -4.78723, -3.4178, -4.84316, -3.41922, -4.88737, -3.44257, -4.88737, -3.44257, -4.89015, -3.44238, -4.97015, -3.44187, -5.01499, -3.44168, -5.09499, -3.44136, -5.13999, -3.44114, -5.21999, -3.44075, -5.26516, -3.44279, -5.3152, -3.44537, -5.3952, -3.44492, -5.44008, -3.44467, -5.52008, -3.44426, -5.56516, -3.44403, -5.64516, -3.44366, -5.68997, -3.44344, -5.76997, -3.44302, -5.81526, -3.44278, -5.89516, -3.4468, -5.94495, -3.44927, -5.99548, -3.45177, -6.03537, -3.48146, -6.07342, -3.53214, -6.09594, -3.58336, -6.12648, -3.6573, -6.14148, -3.71626, -6.15172, -3.78604, -6.14205, -3.84119, -6.06205, -3.84151, -6.01715, -3.84169, -5.93715, -3.8421, -5.89206, -3.84232, -5.81207, -3.84273, -5.76699, -3.84295, -5.68785, -3.85468, -5.60871, -3.8664, -5.52958, -3.87812, -5.45234, -3.88956, -5.39736, -3.87987, -5.34216, -3.87004, -5.28712, -3.86027, -5.23209, -3.85062, -5.15219, -3.84654, -5.0723, -3.84247, -5.05231, -3.84145, -4.97241, -3.8376, -4.92197, -3.83724, -4.84197, -3.83753, -4.7969, -3.83777, -4.71691, -3.83823, -4.67183, -3.83843, -4.59183, -3.83879, -4.54711, -3.839, -4.46711, -3.83937, -4.43143, -3.83379, -4.39914, -3.83567, -4.34327, -3.83463, -4.28733, -3.83347, -4.20747, -3.8381, -4.14184, -3.84191, -4.08594, -3.84088, -4.03013, -3.83987, -3.97418, -3.83873, -3.89431, -3.84334, -3.82854, -3.84715, -3.77267, -3.84609, -3.71672, -3.84493, -3.66095, -3.84413, -3.58108, -3.84874, -3.51534, -3.85253, -3.45951, -3.85129, -3.4035, -3.85035, -3.34754, -3.84918, -3.29182, -3.84822, -3.23607, -3.8471, -3.15609, -3.84869, -3.07611, -3.85029, -2.99612, -3.85188, -2.91614, -3.85347, -2.83615, -3.85506, -2.77004, -3.85638, -2.7588, -3.85766, -2.67893, -3.86218, -2.6287, -3.86503, -2.57858, -3.86787, -2.49871, -3.87233, -2.44876, -3.8702, -2.39361, -3.86089, -2.31765, -3.83972, -2.31481, -3.83873, -2.31037, -3.83875, -2.25433, -3.83906, -2.19845, -3.83598, -2.13808, -3.83116, -2.10701, -3.83127, -2.06513, -3.81597, -2.06033, -3.81398, -1.99629, -3.83316, -1.94736, -3.84302, -1.86817, -3.83793, -1.82208, -3.81877, -1.76936, -3.7586, -1.74054, -3.69552, -1.72471, -3.6171, -1.70888, -3.53868, -1.69306, -3.46026, -1.67723, -3.38184, -1.66141, -3.30343, -1.64779, -3.23598, -1.63821, -3.21955, -1.59622, -3.15146, -1.57041, -3.10316, -1.54224, -3.05478, -1.50037, -2.98662, -1.49527, -2.97711, -1.48681, -2.9639, -1.44506, -2.89566, -1.41077, -2.83947, -1.36915, -2.77115, -1.33485, -2.7148, -1.29325, -2.64646, -1.2565, -2.5919, -1.22808, -2.54401, -1.20403, -2.50012, -1.16229, -2.43188, -1.12803, -2.37586, -1.08638, -2.30755, -1.05204, -2.25122, -1.01035, -2.18295, -0.976081, -2.12683, -0.943624, -2.08139, -0.905314, -2.01115, -0.881258, -1.96736, -0.864868, -1.93757, -0.852933, -1.91983, -0.808275, -1.85346, -0.763617, -1.78708, -0.718959, -1.72071, -0.674301, -1.65433, -0.645067, -1.60963, -0.606052, -1.54938, -0.575654, -1.50266, -0.545161, -1.45578, -0.514778, -1.40875, -0.484331, -1.36195, -0.435266, -1.29876, -0.401757, -1.2556, -0.371163, -1.20875, -0.345223, -1.16836, -0.321475, -1.15082, -0.276606, -1.08458, -0.248607, -1.04315, -0.220534, -1.00159, -0.175639, -0.935372, -0.1459, -0.894693, -0.0980938, -0.830549, -0.0603663, -0.776854, -0.032353, -0.735423, 0.0123763, -0.669096, 0.0404848, -0.627549, 0.0686508, -0.585952, 0.113469, -0.519685, 0.141469, -0.478052, 0.169376, -0.436498, 0.214248, -0.370267, 0.24237, -0.328761, 0.270532, -0.287198, 0.315167, -0.220807, 0.343199, -0.179228, 0.371258, -0.13764, 0.419252, -0.0736363, 0.458505, -0.0212898, 0.503151, 0.0450937, 0.531337, 0.086775, 0.559558, 0.128448, 0.604252, 0.194799, 0.632253, 0.236319, 0.660345, 0.277961, 0.705247, 0.344171, 0.731318, 0.382803, 0.734514, 0.385505, 0.734732, 0.387866, 0.761258, 0.427221, 0.806076, 0.493488, 0.825424, 0.534521, 0.829394, 0.590255, -3.84460529197, -3.82900135089, -3.76460540964, -3.82913856068, -3.71495033922, -3.83117007124, -3.63505189559, -3.83519979889, -3.56902055147, -3.83658476609, -3.51892676239, -3.83667500634, -3.43892697064, -3.83685754034, -3.38882622381, -3.83703303216, -3.33869715289, -3.83722408104, -3.2586973238, -3.83738944761, -3.20839117506, -3.8375155713, -3.1584067702, -3.8376465555, -3.0784070987, -3.8378758163, -3.02813688751, -3.8379083647, -2.97814373788, -3.83791225956, -2.89814407845, -3.83814569522, -2.84786524893, -3.83831920567, -2.79809505061, -3.83849777794, -2.71809539593, -3.83873283143, -2.66900492785, -3.84935510983, -2.59186031575, -3.87053785809, -2.51471570365, -3.89172060636, -2.43757109155, -3.91290335462, -2.36042647945, -3.93408610288, -2.28328186736, -3.95526885115, -2.21592011365, -3.95770077147, -2.16307593811, -3.91521348968, -2.12708025114, -3.84376903209, -2.09108456417, -3.77232457451, -2.07815332484, -3.73292888855, -2.05362342379, -3.65678242471, -2.03460052889, -3.59773097665, -2.00748227285, -3.5224674384, -1.99950994572, -3.46065816842, -1.99936675208, -3.41048352814, -2.02696047509, -3.36191822725, -2.05460686723, -3.31334170307, -2.0962044665, -3.28529373996, -2.13762567848, -3.25770106524, -2.19308664517, -3.25055352901, -2.26810166159, -3.2783518619, -2.34311667802, -3.30615019478, -2.41813169444, -3.33394852767, -2.49314671086, -3.36174686055, -2.56816172728, -3.38954519344, -2.6431767437, -3.41734352633, -2.71354302609, -3.44341917671, -2.79354244643, -3.44311463518, -2.84375602929, -3.4429899247, -2.89370887589, -3.44288273927, -2.97370853605, -3.44264955437, -3.02399973154, -3.44253246532, -3.07419109443, -3.44242318371, -3.15419098987, -3.4422938382, -3.20429594883, -3.44225834405, -3.25444544011, -3.44223432718, -3.33444506631, -3.44198976963, -3.38460735199, -3.44179940065, -3.43478855318, -3.44159954985, -3.51478852488, -3.44153225349, -3.56489678094, -3.44142140448, -3.61505675509, -3.44129306251, -3.69505639907, -3.44105439283, -3.74513095058, -3.43904044074, -3.82503250852, -3.43507294278, -3.88545835061, -3.44071462704, -3.95486234214, -3.45437801521, -4.03479214956, -3.45102750987, -4.11472195697, -3.44767700452, -4.19465176439, -3.44432649918, -4.27458157181, -3.44097599383, -4.35451137923, -3.43762548849, -4.43444118665, -3.43427498315, -4.4818553633, -3.43228747114, -4.56185471357, -3.43196504842, -4.60654200111, -3.43883254971, -4.68652435287, -3.43715224976, -4.76650670462, -3.4354719498, -4.84648905638, -3.43379164985, -4.92647140813, -3.43211134989, -5.00645375988, -3.43043104994, -5.08643611164, -3.42875074998, -5.16641846339, -3.42707045003, -5.24640081514, -3.42539015007, -5.3263831669, -3.42370985012, -5.38853899809, -3.42309471871, -5.43857377542, -3.42290750312, -5.51857374387, -3.4228364482, -5.56855078121, -3.42254646646, -5.61889810939, -3.42219261746, -5.69889787103, -3.4223879059, -5.74904069333, -3.42224832631, -5.79906185715, -3.42204281239, -5.87906120845, -3.42172064468, -5.92933357241, -3.42157766291, -5.97934683101, -3.42145061243, -6.05925094485, -3.41753492705, -6.12481820616, -3.41432181631, -6.20481719041, -3.41391867807, -6.2553648537, -3.41377550739, -6.30526332591, -3.41366321382, -6.38526331167, -3.41361548099, -6.43545380999, -3.41356517604, -6.48533692631, -3.41351002441, -6.56533689772, -3.41344238596, -6.6154585591, -3.41135902962, -6.6953508716, -3.40720952335, -6.7613936331, -3.40577722896, -6.81134372749, -3.40567775068, -6.89134351581, -3.40549371962, -6.94158754726, -3.40537322491, -6.99155977241, -3.40525213044, -7.07155968949, -3.40513694929, -7.12179525101, -3.40496005024, -7.17193373625, -3.40475678492, -7.24140803634, -3.41858973932, -7.26930020062, -3.45999080177, -7.29682111845, -3.50155641234, -7.32636971025, -3.55928081314, -7.35699748917, -3.63318574638, -7.38762526808, -3.70709067963, -7.40861204281, -3.75773183745, -7.36169827837, -3.79448404428, -7.28180082377, -3.79853333409, -7.21620573024, -3.79992253968, -7.1659977426, -3.80003803477, -7.08599802119, -3.80024915972, -7.03605066192, -3.80031331068, -6.98584523617, -3.80036087924, -6.90584606208, -3.80072439656, -6.85565374297, -3.80275372855, -6.77575209027, -3.80671931756, -6.70998286143, -3.80815431498, -6.65988319544, -3.8083497464, -6.57997688387, -3.81222032408, -6.51434672589, -3.8153993799, -6.4344547492, -3.81124341307, -6.38945352434, -3.80890245713, -6.30945370867, -3.80907419415, -6.2592619437, -3.81108788911, -6.17936117792, -3.81507130831, -6.11832188645, -3.83480137123, -6.04734670986, -3.87171375815, -5.97637153327, -3.90862614506, -5.90539635668, -3.94553853198, -5.83442118009, -3.9824509189, -5.7634460035, -4.01936330582, -5.6924708269, -4.05627569273, -5.62149565031, -4.09318807965, -5.55116256828, -4.1297665295, -5.47488661224, -4.15089592967, -5.42647739228, -4.12327464324, -5.35842361949, -4.0812189182, -5.2903698467, -4.03916319316, -5.2223160739, -3.99710746812, -5.15426230111, -3.95505174308, -5.08620852832, -3.91299601804, -5.01815475553, -3.87094029299, -4.95010098274, -3.82888456795, -4.91247492605, -3.80563249787, -4.84311373223, -3.79183448947, -4.76843992134, -3.78557932998, -4.68845043522, -3.78687629199, -4.6084609491, -3.788173254, -4.52847146298, -3.78947021601, -4.44848197686, -3.79076717802, -4.36849249075, -3.79206414003, -4.30242177961, -3.79313542085, -4.22713693455, -3.82019446933, -4.16544060725, -3.8282156323, -4.11519277984, -3.82849551347, -4.03519296498, -3.82866762414, -3.9850251929, -3.82869939508, -3.93478107715, -3.82871182562, -3.85478148949, -3.82896867842]
	points1 = points[:323*2]
	points2 = points[323*2:]
	
	nPoints1 = []
	nPoints2 = []
	
	print len(points1)
	print len(points2)


	for i in range(len(points1)/2-1):
		if False and points1[2*i] == points1[2*i+2] and points1[2*i+1] == points1[2*i+3]:
			print "points1 ", i
			print points1[2*i], points1[2*i+2], points1[2*i+1], points1[2*i+3]
		else:
			nPoints1.append([points1[2*i], points1[2*i+1]])

	nPoints1.append([points1[-2], points1[-1]])
			

	for i in range(len(points2)/2-1):
		if False and points2[2*i] == points2[2*i+2] and points2[2*i+1] == points2[2*i+3]:
			print "points2 ", i
		else:
			nPoints2.append([points2[2*i],points2[2*i+1]])

	nPoints2.append([points2[-2], points2[-1]])
			
	print len(nPoints1)
	print len(nPoints2)
	
	xP = []
	yP = []
	for i in range(len(points1)/2):
		
		xP.append(points1[i*2])
		yP.append(points1[i*2+1])

	pylab.plot(xP,yP)

	xP = []
	yP = []
	for i in range(len(points2)/2):
		xP.append(points2[i*2])
		yP.append(points2[i*2+1])

	pylab.plot(xP,yP)
	#pylab.show()
	
	computeUnion(nPoints1, nPoints2)
	
	exit()
	
	# TUNE ME:  threshold cost difference between iterations to determine if converged
	costThresh = 0.004

	# TUNE ME:   minimum match distance before point is discarded from consideration
	minMatchDist = 2.0

	# plot the best fit at each iteration of the algorithm?
	plotIteration = True

	# initial guess for x, y, theta parameters
	offset = [0.0,0.0,-math.pi/4]

	# sample data
	a_data = [[1.0,1.0],[1.1,1.1],[1.2,1.2],[1.3,1.31],[1.4,1.4],[1.51,1.5],[1.6,1.6]]
	b_data = [[0.3,1.0],[0.3,1.1],[0.3,1.2],[0.31,1.3],[0.3,1.4],[0.3,1.5],[0.3,1.6]]

	# treat the points with the point-to-line constraint
	addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
	addPointToLineCovariance(b_data, high_var=1.0, low_var=0.001)

	#addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
	#addDistanceFromOriginCovariance(b_data, tan_var=0.1, perp_var=0.01)

	# plot the data without A transformed, plot 997
	draw(a_data, b_data, "rawData.png")

	# transform the points in A by 'offset'
	a_trans = []
	for p in a_data:
		a_trans.append(dispOffset(p, offset))

	# plot the data with A transformed, plot 998
	draw(a_trans, b_data, "initialGuess.png") 

	# run generalized ICP (a plot is made for each iteration of the algorithm)
	offset = gen_ICP(offset, a_data, b_data, costThresh, minMatchDist, plotIteration)

	# transform the points of A with parameters determined by algorithm
	a_trans = []
	for p in a_data:
		a_trans.append(dispPoint(p, offset))

	# plot the final result, plot 999
	draw(a_trans, b_data, "finalOutput.png") 
