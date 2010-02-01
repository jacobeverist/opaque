#!/usr/bin/python

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
	

if __name__ == '__main__':


	# 1. compute all the relative offsets between the odd-numbered poses
	# 2. compute the global position of all the poses wrt to first pose

	# TUNE ME:  threshold cost difference between iterations to determine if converged
	#costThresh = 0.004
	costThresh = 0.1

	# TUNE ME:   minimum match distance before point is discarded from consideration
	minMatchDist = 2.0

	# plot the best fit at each iteration of the algorithm?
	plotIteration = True
	#plotIteration = False

	# initial guess for x, y, theta parameters
	offset = [0.0,0.0,0.0]

	# sample data
	a_data = []
	b_data = []

	finalOffsets = []
	estPoses = []
	gndPoses = []
	poseNumbers = []

	for i in range(1,9):
	
		f = open('cntpose%04u.txt' % i,'r')
		val = f.read()
		val = val.rstrip('\r\n')
		b_cnt = eval(val)

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

	plotEnv()	

	count = 0
	pylab.savefig("uncorrectedMap_%04u.png" % count)

	count += 1

	for i in range(len(estPoses)):

		est1 = estPoses[i]

		#f = open('alpha_bound_%04u.txt' % poseNumbers[i],'r')
		f = open('est_occ_points_%04u.txt' % poseNumbers[i],'r')
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


	pylab.xlim(-3,6)
	pylab.ylim(-4,4)
	gen_icp.save("uncorrectedMap.png")

	#exit()

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

		#f = open('alpha_bound_%04u.txt' % poseNumbers[i],'r')
		f = open('est_occ_points_%04u.txt' % poseNumbers[i],'r')
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


		
		
	
