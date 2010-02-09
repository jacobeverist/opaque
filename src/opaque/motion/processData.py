#!/usr/bin/python

import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import traceback
import Image
import ImageDraw
from copy import *
from math import *
from numpy import arange
from maps import *
from pose import *
from robot import *

import numpy
import math
import gen_icp
import functions

import scipy
import scipy.optimize

from copy import copy
from copy import deepcopy
from coord import Pose
from SplineFit import SplineFit

from numpy import array
import pylab


plotCount = 0

def normalizeAngle(angle):

	while angle>math.pi:
		angle=angle-2*math.pi

	while angle<=-math.pi:
		angle=angle+2*math.pi

	return angle 

def sanitizeMotion(offset, a_cnt, b_cnt):

	# optimize the match error for the current list of match pairs
	newOffset = scipy.optimize.fmin(sanitizeCost, offset, [a_cnt, b_cnt])

	return newOffset

def sanitizeCost(offset, a_cnt, b_cnt):

	#offset = [0.0,0.0,0.0]

	# transform the points in A by 'offset'
	a_trans = []
	for p in a_cnt:
		a_trans.append(gen_icp.dispOffset(p, offset))

	# fit splines
	splineA = SplineFit(a_trans, kp=3)
	splineB = SplineFit(b_cnt, kp=3)

	cost1 = splineB.matchSpline(splineA, plot=False, index=0)

	return cost1

def decimatePoints(points):
	result = []

	for i in range(len(points)):
		if i%2 == 0:
			result.append(points[i])

	return result

def modifyGuess(est1, newOffset, control_points):
	
	pose1 = Pose(est1)

	" first node stays in its current position "
	cntPoints = deepcopy(control_points)

	xP = []
	yP = []
	cntPoints1 = []
	for j in range(len(cntPoints)):
		pnt = cntPoints[j]
		
		# convert local points to global coordinates for node1
		pnt = pose1.convertLocalToGlobal(pnt)

		xP.append(pnt[0])
		yP.append(pnt[1])
		
		cntPoints1.append(copy(pnt))

	spline1 = SplineFit(cntPoints1, kp=3)
	
	samples = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
	totalVec = [0.0, 0.0]
	for i in range(len(samples)):
	
		uvec = spline1.getUVector(samples[i])
		totalVec[0] += uvec[0]
		totalVec[1] += uvec[1]
	
	" normalize the vector "
	mag = math.sqrt(totalVec[0]**2 + totalVec[1]**2)
	totalVec[0] /= mag
	totalVec[1] /= mag
	
	" convert the totalVec into local coordinates "
	globalVec = numpy.array([[totalVec[0]], [totalVec[1]]])
	finalVec = numpy.dot(pose1.R, globalVec)
	
	" modified guess in local coordinates "
	newGuess = copy(newOffset)
	newGuess[0] += finalVec[0,0] * 1.0
	newGuess[1] += finalVec[1,0] * 1.0
	
	return newGuess

def plotEnv():
	WLEN = 3.0
	wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*math.cos(math.pi/3), -0.2 - WLEN*math.sin(math.pi/3)]]
	wall2 = [[-4.0 + WLEN*math.cos(math.pi/3), 0.2 + WLEN*math.sin(math.pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
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

	pylab.xlim(-3,6)
	pylab.ylim(-4,4)


def computeUnions(point_sets):

	currPoly = point_sets[0]

	for i in range(1,len(point_sets)):
		currPoly = computeUnion(currPoly,point_sets[i])

	return currPoly

def computeUnion(points1, points2):
	
	numPoints1 = len(points1)
	numPoints2 = len(points2)

	inputStr = str(numPoints1) + " " + str(numPoints2) + " "
	
	" Removed Gaussians because were causing self-intersections ."
	" not necessary anymore because polygons come from CGAL and are not degenerate. "
	for p in points1:
		p2 = copy(p)
		#p2[0] += gauss(0.0,0.0001)
		#p2[1] += gauss(0.0,0.0001)

		inputStr += str(p2[0]) + " " + str(p2[1]) + " "

	for p in points2:
		p2 = copy(p)
		#p2[0] += gauss(0.0,0.0001)
		#p2[1] += gauss(0.0,0.0001)

		inputStr += str(p2[0]) + " " + str(p2[1]) + " "
	
	inputStr += "\n"
	
	#print "sending input to Union:  ", inputStr
	
	" start the subprocess "
	subProc = Popen(["./poly_union.exe"], stdin=PIPE, stdout=PIPE)
	
	" send input and receive output "
	sout, serr = subProc.communicate(inputStr)

	#print "rawOutput ="
	#print sout

	print serr
	
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
				
if __name__ == '__main__':

	probe = FakeProbe(40,0.15,0.05)
	contacts = ContactReferences(probe)
	mapGraph = MapGraph(probe, contacts)


	#numPoses = 8
	#numPoses = 11
	numPoses = 34

	mapGraph.loadFile(numPoses)
	mapGraph.saveMap()
	#exit()

	" TUNE ME:  threshold cost difference between iterations to determine if converged "
	#costThresh = 0.004
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
	for i in range(0,numPoses):
	
		f = open('estpose%04u.txt' % i,'r')
		val = f.read()
		val = val.rstrip('\r\n')
		estPose1 = eval(val)

		f = open('gndpose%04u.txt' % i,'r')
		val = f.read()
		val = val.rstrip('\r\n')
		gndPose1 = eval(val)

		estPoses.append(estPose1)
		gndPoses.append(gndPose1)
		poseNumbers.append(i)
	
	pylab.clf()
	
	" Read in data of Alpha-Shapes and add their associated covariances "
	a_hulls = []
	for i in range(len(estPoses)):
		f = open('alpha_bound_%04u.txt' % i,'r')
		val = f.read()
		a_data = eval(val)
		a_data = decimatePoints(a_data)

		" treat the points with the point-to-line constraint "
		gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)

		" treat the points with the distance-from-origin increasing error constraint "
		gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)

		a_hulls.append(a_data)

	occMaps = []
	for m in range(0,len(estPoses)):
		offset = estPoses[m]
		occMap = mapGraph.getNodeOccMap(m)

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
		

	a_hull_trans = []
	for m in range(0,len(estPoses)):
		offset = estPoses[m]

		" transform the past poses "
		past_data = a_hulls[m]
		a_trans = []
		for p in past_data:
			a_trans.append(gen_icp.dispPoint(p, offset))
		
		a_hull_trans.append(a_trans)
		
	pylab.clf()
	plotEnv()

	for points in occMaps:
		xP = []
		yP = []
		for p in points:
			xP.append(p[0])
			yP.append(p[1])
		pylab.scatter(xP,yP, linewidth=1, color=(1.0,0.6,0.6))
	
	for hull in a_hull_trans:
		
		xP = []
		yP = []
		for p in hull:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull[0][0])
		yP.append(hull[0][1])
		
		pylab.plot(xP,yP,color='b')
		
	pylab.savefig("uncorrectedPoses.png")	
	
	offsets = []
	for i in range(len(estPoses)-1):
		estPose1 = estPoses[i]
		estPose2 = estPoses[i+1]
		
		pose1 = Pose(estPose1)
		offset = pose1.convertGlobalPoseToLocal(estPose2)
		offsets.append(offset)
	
	for k in range(1,len(estPoses)):
		print k

		currEstPoses = estPoses[:k+1]
		curr_a_hulls = a_hulls[:k+1]
		
		print len(currEstPoses),len(curr_a_hulls)
		
		
		" 1. target pose to correct "
		targetPose = currEstPoses[k]
		targetHull = curr_a_hulls[k]
		
		" 2. past poses in single alpha shape "
		mapGraph.loadFile(k)
		
		" update the estimated poses "
		for m in range(0,k):
			mapGraph.setNodePose(m, estPoses[m])
					
		mapGraph.saveMap()
	
		a_hull_trans = []
		estPoseOrigin = estPoses[k-1]
		for m in range(0,k):
			estPose2 = estPoses[m]
			poseOrigin = Pose(estPoseOrigin)
			offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
	
			" transform the past poses "
			past_data = a_hulls[m]
			a_trans = []
			for p in past_data:
				a_trans.append(gen_icp.dispPoint(p, offset))
			
			a_hull_trans.append(a_trans)

		pastCircles = []
		for m in range(0,k):
			hull = a_hull_trans[m]
			radius, center = gen_icp.computeEnclosingCircle(hull)
			pastCircles.append([radius,center])
			
		pastPose = estPoses[k-1]
		
		pastHull = computeUnions(a_hull_trans)
		gen_icp.addPointToLineCovariance(pastHull, high_var=1.0, low_var=0.001)
		#gen_icp.addDistanceFromOriginCovariance(pastHull, tan_var=0.1, perp_var=0.01)
		
		"""	
		globalHull = mapGraph.computeAlphaBoundary()		
		globalHull = decimatePoints(globalHull)
		pastPose = estPoses[k-1]
		
		" convert points to local coordinates "
		refPose = Pose(pastPose)
		pastHull = []
		for p in globalHull:
			localPoint = refPose.convertGlobalToLocal(p)
			pastHull.append(localPoint)
		
		gen_icp.addPointToLineCovariance(pastHull, high_var=1.0, low_var=0.001)
		gen_icp.addDistanceFromOriginCovariance(pastHull, tan_var=0.1, perp_var=0.01)
		"""
		
		print "pastHull =", len(pastHull)
		print "targetHull =", len(targetHull)
		
		print "pastPose =", pastPose
		print "targetPose =", targetPose
				
		" run generalized ICP (a plot is made for each iteration of the algorithm) "
		offset = gen_icp.gen_ICP_global(pastPose, targetPose, pastHull, targetHull, pastCircles, costThresh, minMatchDist, plotIteration)
		#offset = gen_icp.gen_ICP_past(currEstPoses, curr_a_hulls, costThresh, minMatchDist, plotIteration)
		
		offsets[k-1] = offset
		
		" recompute the estimated poses with the new offset "
		newEstPoses = []
		newEstPoses.append(estPoses[0])
		
		for m in range(len(offsets)):
			pose1 = Pose(newEstPoses[m])
			offset = offsets[m]
			newEstPose2 = pose1.convertLocalOffsetToGlobal(offset)
			newEstPoses.append(newEstPose2)
			
			
		estPoses = newEstPoses
		
		#mapGraph.saveMap()
	
	" 2. all poses in single alpha shape "
	mapGraph.loadFile(len(estPoses))
	
	" update the estimated poses "
	for m in range(0,len(estPoses)):
		mapGraph.setNodePose(m, estPoses[m])
				
	mapGraph.saveMap()

	occMaps = []
	for m in range(0,len(estPoses)):
		offset = estPoses[m]
		occMap = mapGraph.getNodeOccMap(m)

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
		

	a_hull_trans = []
	for m in range(0,len(estPoses)):
		offset = estPoses[m]

		" transform the past poses "
		past_data = a_hulls[m]
		a_trans = []
		for p in past_data:
			a_trans.append(gen_icp.dispPoint(p, offset))
		
		a_hull_trans.append(a_trans)

	pastHull = computeUnions(a_hull_trans)
	
	pylab.clf()
	plotEnv()
	
	xP = []
	yP = []
	for p in pastHull:
		xP.append(p[0])
		yP.append(p[1])
	xP.append(pastHull[0][0])
	yP.append(pastHull[0][1])
	
	pylab.plot(xP,yP,color='b')
	pylab.savefig("finalMap.png")
	
	pylab.clf()
	plotEnv()
	
	for points in occMaps:
		xP = []
		yP = []
		for p in points:
			xP.append(p[0])
			yP.append(p[1])
		pylab.scatter(xP,yP, linewidth=1, color=(1.0,0.6,0.6))
	
	for hull in a_hull_trans:
		
		xP = []
		yP = []
		for p in hull:
			xP.append(p[0])
			yP.append(p[1])
		xP.append(hull[0][0])
		yP.append(hull[0][1])
		
		pylab.plot(xP,yP,color='b')
		
	pylab.savefig("finalPoses.png")
	
	" relative distance error "
	" separate error into two components, the direct error and the lateral error "
	" MOTION MODEL "
	
	"""
	for i in range(len(estPoses)-1):
		estPose1 = estPoses[i]
		estPose2 = estPoses[i+1]
		
		gndPose1 = gndPoses[i]
		gndPose2 = gndPoses[i+1]

		estDist = math.sqrt((estPose1[0]-estPose2[0])**2 + (estPose1[1]-estPose2[1])**2)
		gndDist = math.sqrt((gndPose1[0]-gndPose2[0])**2 + (gndPose1[1]-gndPose2[1])**2)

		estVec = [estPose2[0]-estPose1[0], estPose2[1]-estPose1[1]]
		gndVec = [gndPose2[0]-gndPose1[0], gndPose2[1]-gndPose1[1]]
		
		" rotate both vectors by ground angle to 0 "
		groundAng = gndPose1[2]
		rotGndVec1 = [gndVec[0]*math.cos(groundAng) + gndVec[1]*math.sin(groundAng), -gndVec[0]*math.sin(groundAng) + gndVec[1]*math.cos(groundAng)]

		estAng = estPose1[2]
		rotEstVec1 = [estVec[0]*math.cos(estAng) + estVec[1]*math.sin(estAng), -estVec[0]*math.sin(estAng) + estVec[1]*math.cos(estAng)]
		
		" rotate estVec to [1,0] "
		normVec = [rotGndVec1[0]/gndDist, rotGndVec1[1]/gndDist]
		rotAng = math.acos(normVec[0])
		if normVec[1] < 0.0:
			rotAng = -rotAng
		
		rotEstVec2 = [rotEstVec1[0]*math.cos(rotAng) + rotEstVec1[1]*math.sin(rotAng), -rotEstVec1[0]*math.sin(rotAng) + rotEstVec1[1]*math.cos(rotAng)]
		rotGndVec2 = [rotGndVec1[0]*math.cos(rotAng) + rotGndVec1[1]*math.sin(rotAng), -rotGndVec1[0]*math.sin(rotAng) + rotGndVec1[1]*math.cos(rotAng)]
				
		#print gndDist, (rotEstVec2[0]-rotGndVec2[0]), (rotEstVec2[1]-rotGndVec2[1])

		rotEstVec3 = [rotEstVec1[0]*math.cos(estAng) - rotEstVec1[1]*math.sin(estAng), rotEstVec1[0]*math.sin(estAng) + rotEstVec1[1]*math.cos(estAng)]
		rotGndVec3 = [rotGndVec1[0]*math.cos(estAng) - rotGndVec1[1]*math.sin(estAng), rotGndVec1[0]*math.sin(estAng) + rotGndVec1[1]*math.cos(estAng)]
		
		xP = [estPose1[0],estPose1[0] + rotEstVec3[0]]
		yP = [estPose1[1],estPose1[1] + rotEstVec3[1]]		
		pylab.plot(xP,yP, color='r')

		xP = [estPose1[0],estPose1[0] + rotGndVec3[0]]
		yP = [estPose1[1],estPose1[1] + rotGndVec3[1]]		
		pylab.plot(xP,yP, color='b')
	"""

	"""
	for i in range(len(estPoses)-1):
		estPose1 = estPoses[i]
		estPose2 = estPoses[i+1]
		
		xP = [estPose1[0],estPose2[0]]
		yP = [estPose1[1],estPose2[1]]
		
		#pylab.plot(xP,yP, color='r')

	for i in range(len(estPoses)):
		estPose1 = estPoses[i]
		gndPose1 = gndPoses[i]
		
		xP = [estPose1[0],gndPose1[0]]
		yP = [estPose1[1],gndPose1[1]]
		
		#pylab.plot(xP,yP, color='b')
	"""

	"""
	" Read in data of Alpha-Shapes and add their associated covariances "
	a_hulls = []
	for i in range(len(estPoses)):
		f = open('alpha_bound_%04u.txt' % i,'r')
		val = f.read()
		a_data = eval(val)
		a_data = decimatePoints(a_data)

		" treat the points with the point-to-line constraint "
		gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)

		" treat the points with the distance-from-origin increasing error constraint "
		gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)

		a_hulls.append(a_data)
	
	offsets = []
	for i in range(len(estPoses)-1):
		estPose1 = estPoses[i]
		estPose2 = estPoses[i+1]
		
		pose1 = Pose(estPose1)
		offset = pose1.convertGlobalPoseToLocal(estPose2)
		offsets.append(offset)
	
	for k in range(1,len(estPoses)):
		print k

		currEstPoses = estPoses[:k+1]
		curr_a_hulls = a_hulls[:k+1]
		
		print len(currEstPoses),len(curr_a_hulls)
		
		" run generalized ICP (a plot is made for each iteration of the algorithm) "
		offset = gen_icp.gen_ICP_past(currEstPoses, curr_a_hulls, costThresh, minMatchDist, plotIteration)
		
		offsets[k-1] = offset
		
		" recompute the estimated poses with the new offset "
		newEstPoses = []
		newEstPoses.append(estPoses[0])
		
		for m in range(len(offsets)):
			pose1 = Pose(newEstPoses[m])
			offset = offsets[m]
			newEstPose2 = pose1.convertLocalOffsetToGlobal(offset)
			newEstPoses.append(newEstPose2)
			
			
		estPoses = newEstPoses
		
		for m in range(len(estPoses)):
			mapGraph.setNodePose(m, estPoses[m])
	
		mapGraph.saveMap()
	"""
	exit()
	
	plotEnv()	
	#pylab.show()

	" 1. build data structure holding the relative positions between poses "
	" 2. query and return local poses that are close enough to a particular local pose "
	" 3. perform a simultaneous iterative closest fit between all the closest local poses "
	" 4. when the position of a local pose changes, alter its relative pose and propagate its changes globally"
	" 5. cycle through all of the local poses like 10 times and see how the system changes over time "
	" 6. add in the cost of the motion model in the iterative closest fit "

	#exit()
	

	" 2. query and return local poses that are close enough to a particular local pose "
	circles = []
	for i in range(len(estPoses)):

		est1 = estPoses[i]

		f = open('alpha_bound_%04u.txt' % poseNumbers[i],'r')
		val = f.read()
		a_data = eval(val)
		b_data = []
		a_trans = []
		for p in a_data:
			a_trans.append(gen_icp.dispOffset(p, est1))

		radius, center = gen_icp.computeEnclosingCircle(a_trans)

		circles.append([radius, center])

	print circles
	
	
	print
	for i in range(len(circles)):
		distances = []
		for j in range(0,i):
			radius1 = circles[i][0]
			radius2 = circles[j][0]
			
			center1 = circles[i][1]
			center2 = circles[j][1]
			
			dist = math.sqrt((center1[0]-center2[0])**2 + (center1[1]-center2[1])**2)
			
			maxDist = radius1 + radius2
			
			if dist < maxDist:
				distances.append(1)
			else:
				distances.append(0)
		print distances

	print a_data
	exit()

	count = 0
	pylab.savefig("uncorrectedMap_%04u.png" % count)

	#count += 1

	for i in range(len(estPoses)):

		est1 = estPoses[i]

		f = open('alpha_bound_%04u.txt' % poseNumbers[i],'r')
		val = f.read()
		a_data = eval(val)
		b_data = []
		a_trans = []
		for p in a_data:
			a_trans.append(gen_icp.dispOffset(p, est1))

		gen_icp.draw(a_trans, b_data, "foo.png", fileWrite = False) 

		gnd1 = gndPoses[i]
		a_data = eval(val)
		b_data = []
		a_trans = []
		for p in a_data:
			a_trans.append(gen_icp.dispOffset(p, gnd1))
		gen_icp.draw(b_data, a_trans, "foo.png", fileWrite = False) 

		pylab.xlim(-3,6)
		pylab.ylim(-4,4)
		pylab.savefig("uncorrectedMap_%04u.png" % count)
		count += 1

	#print a_trans


	pylab.xlim(-3,6)
	pylab.ylim(-4,4)
	gen_icp.save("uncorrectedMap.png")

	exit()

	finalOffsets = []
	estPoses = []
	#gndPoses = []
	poseNumbers = []
	#raise
	#for i in range(1,29):
	#for i in range(5,6):
	for i in range(1,8):
	#for i in range(2,3):
	
		#if i % 2 == 1:
		if True:

			print "correcting", i , "and", i+1

			f = open('alpha_bound_%04u.txt' % i,'r')
			#f = open('gnd_alpha_bound_%04u.txt' % i,'r')
			val = f.read()
			b_data = eval(val)

			f = open('alpha_bound_%04u.txt' % (i+1),'r')
			#f = open('gnd_alpha_bound_%04u.txt' % (i+1),'r')
			val = f.read()
			a_data = eval(val)

			a_data = decimatePoints(a_data)
			b_data = decimatePoints(b_data)

			f = open('cntpose%04u.txt' % i,'r')
			val = f.read()
			val = val.rstrip('\r\n')
			b_cnt = eval(val)
			f = open('cntpose%04u.txt' % (i+1) ,'r')
			val = f.read()
			val = val.rstrip('\r\n')
			a_cnt = eval(val)

			f = open('estpose%04u.txt' % i,'r')
			val = f.read()
			val = val.rstrip('\r\n')
			estPose1 = eval(val)

			f = open('estpose%04u.txt' % (i+1),'r')
			val = f.read()
			val = val.rstrip('\r\n')
			estPose2 = eval(val)

			# get the offset of node 2 with respect to node 1
			pose1 = Pose(estPose1)
			offset = pose1.convertGlobalPoseToLocal(estPose2)

			# make a more realistic motion estimate
			#newOffset = sanitizeMotion(offset, a_cnt, b_cnt)
			#offset[2] = newOffset[2]

			#offset = modifyGuess(estPose1, newOffset, b_cnt)

			# treat the points with the point-to-line constraint
			gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
			gen_icp.addPointToLineCovariance(b_data, high_var=1.0, low_var=0.001)

			gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
			gen_icp.addDistanceFromOriginCovariance(b_data, tan_var=0.1, perp_var=0.01)

			# plot the data without A transformed, plot 997
			gen_icp.draw(a_data, b_data, "rawData.png")

			# transform the points in A by 'offset'
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispOffset(p, offset))

			# plot the data with A transformed, plot 998
			gen_icp.draw(a_trans, b_data, "initialGuess.png") 

			# run generalized ICP (a plot is made for each iteration of the algorithm)
			offset = gen_icp.gen_ICP(offset, a_data, b_data, costThresh, minMatchDist, plotIteration)

			# transform the points of A with parameters determined by algorithm
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispPoint(p, offset))

			print "drawing: ", "finalOutput_%04u.png" % i

			# plot the final result, plot 999
			gen_icp.draw(a_trans, b_data, "finalOutput_%04u.png" % i, fileWrite = True) 

			finalOffsets.append(offset)
			estPoses.append(estPose1)
			poseNumbers.append(i)


	print "finalOffsets"
	print finalOffsets
	print "estPoses"
	print estPoses

	# 1. get the relative offsets (done)
	# 2. for each node, recompute the estimated pose of all the subsequent poses
	# 3. plot them

	
	for i in range(len(finalOffsets)):
		
		est1 = estPoses[i]

		offset = finalOffsets[i]

		pose1 = Pose(est1)

		newEst2 = pose1.convertLocalOffsetToGlobal(offset)

		if i+1 < len(estPoses):
			estPoses[i+1] = newEst2
		else:
			estPoses.append(newEst2)
			#poseNumbers.append(poseNumbers[i]+2)
			poseNumbers.append(poseNumbers[i]+1)

	print "estPoses"
	print estPoses
	
	pylab.clf()

	for i in range(len(estPoses)):

		est1 = estPoses[i]

		f = open('alpha_bound_%04u.txt' % poseNumbers[i],'r')
		val = f.read()
		a_data = eval(val)
		b_data = []
		a_trans = []
		for p in a_data:
			a_trans.append(gen_icp.dispOffset(p, est1))
			
		gen_icp.draw( a_trans, b_data,"foo.png", fileWrite = False) 


		gnd1 = gndPoses[i]
		a_data = eval(val)
		b_data = []
		a_trans = []
		for p in a_data:
			a_trans.append(gen_icp.dispOffset(p, gnd1))
		gen_icp.draw(b_data, a_trans, "foo.png", fileWrite = False) 
		
		
	plotEnv()

	pylab.xlim(-3,6)
	pylab.ylim(-4,4)

	gen_icp.save("finalMap.png")

	#pylab.xlim(-4.5,4.5)
	#pylab.ylim(-4,4)
	#pylab.savefig("finalMap.png")
	#pylab.clf()


		
		
	
