#!/usr/bin/python

import numpy
import scipy 
import scipy.linalg
import scipy.optimize
import math
import pylab
import pca_module

from copy import copy
from copy import deepcopy
from coord import Pose
from SplineFit import SplineFit


plotCount = 0

def normalizeAngle(angle):

	while angle>math.pi:
		angle=angle-2*math.pi

	while angle<=-math.pi:
		angle=angle+2*math.pi

	return angle 

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


def disp(ai, bi, T):

	temp = T*ai
	result = bi-temp	

	result[2] = 1.0

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
	newGuess[0] -= finalVec[0,0] * 1.0
	newGuess[1] -= finalVec[1,0] * 1.0
	
	return newGuess

def computeMatchError(offset, a, b, Ca, Cb):

	# NOTE: Ca needs to be pre-transformed apparently?

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
			[math.sin(theta), math.cos(theta), yd],
			[0.0, 0.0, 1.0]
		    ])

	R = numpy.matrix([	[math.cos(theta), -math.sin(theta)],
		[math.sin(theta), math.cos(theta)] ])

	# now rotate the covariance matrix by the rotation matrix
	#Cv = R * Ca * numpy.transpose(R)
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

	# NOTE:  seems to create opposite, not orthogonal, but colinear vectors if data is colinear

	try:
		#print cov_a
		#cov_a = numpy.array([[ 1.,  1.], [ 1., 1.0000001]])
		#print cov_a
			
		scores, loadings, E = pca_module.nipals_mat(cov_a, 2, 0.000001, False)
		#print loadings
	except:
		raise

	if len(loadings) < 2:
		raise

	#return scores, loadings, E
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

	#Ca = R * Cv * numpy.transpose(R)
	Ca = numpy.transpose(R) * Cv * R

	return Ca

def findRotationTangent(pnt):
	rotVec = [-pnt[1],pnt[0]]
	return rotVec

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


def cost_func(offset, match_pairs):

	#xd = offset[0]
	#yd = offset[1]
	#theta = offset[2]

	#T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
	#		[math.sin(theta), math.cos(theta), yd],
	#		[0.0, 0.0, 1.0]
	#	    ])

	sum = 0.0
	for pair in match_pairs:

		a = numpy.matrix([[pair[0][0]],[pair[0][1]]])
		b = numpy.matrix([[pair[1][0]],[pair[1][1]]])
		sum += computeMatchError(offset, a, b, pair[2], pair[3])

		# NOTE:  standard point-to-point ICP
		#a = numpy.matrix([[pair[0][0]],[pair[0][1]],[1.0]])
		#b = numpy.matrix([[pair[1][0]],[pair[1][1]],[1.0]])

		#distVec = disp(a, b, T)
		#mag = distVec[0,0]**2 + distVec[1,0]**2
		#sum += mag

	return sum

def precomputeCovariance(points, high_var=1.0, low_var=0.001):

	for p in points:
		C = numpy.matrix([	[high_var, 0.0],
				[0.0, high_var]
				])

		try:
			# the covariance matrix that enforces the point-to-plane constraint
			normVec = findLocalNormal(p,points)
			C = computeVectorCovariance(normVec,low_var,high_var)


			# the covariance matrix that represents increasing uncertainty as distance from root increases
			pVec = p
			dist = math.sqrt(pVec[0]**2 + pVec[1]**2)
			
			# 1. get the orthogonal vector, rotation 90 degrees
			rotVec = [0.,0.]
			rotVec[0] = -pVec[1]
			rotVec[1] = pVec[0]

			# 2. normalize vector	
			rotVec[0] /= dist
			rotVec[1] /= dist

			# 3. computeVectorCovariance() with x_var, y_var proportional to dist, and x_var larger
			Cr = computeVectorCovariance(rotVec,dist*0.1,dist*0.01)
			
			#print Cr + C

		except:
			pass

		#p.append(C)
		p.append(C + Cr)

def gen_ICP(offset, a_data, b_data):


	# DONE:  store the covariance for each point untransformed so we don't have to repeat
	# DONE:  transform each covariance by the appropriate rotation at each point

	# 1. compute the local line for each point
	# 2. compute the covariance for point-to-line constraint
	# 3. compute the covariance for rotational error constraint for A points only
	# 4. add the covariances together

	# 5. find closest points in A of B
	# 6. discard matches beyond threshold d_max < | ai - T *bi |
	# 7. optimize T for the sum of computeMatchError
	# 8. if converged, stop, else go to 5 with b offset by new T

	# TODO:  bias the points of the alpha boundary with custom covariances
	# 1. rotational error wrt root node
	# 2. amplitfy covariance for boundary points that are off the true boundary

	# TODO:  translate the offset parameters between the Probe code offset method

	# TODO:  port the code over to the Probe codebase

	# compute the local covariance matrix for a point-to-plane contrainst
	precomputeCovariance(a_data)
	precomputeCovariance(b_data)

	#costThresh = 0.4
	costThresh = 0.04
	lastCost = 1e100
	minMatchDist = 1.0
	numIterations = 0
	while True:

		a_trans = []

		# pre-transform the A points and their associated covariances
		for p in a_data:
			a_trans.append(dispPoint(p, offset))

		match_pairs = []
		for i in range(len(a_trans)):
			a_p = a_trans[i]

			# for every transformed point of A, find it's closest neighbor in B
			b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])

			if minDist <= minMatchDist:
	
				# add to the list of match pairs less than 1.0 distance apart 
				# keep A points and covariances untransformed
				Ca = a_data[i][2]
				
				#print "selection:"
				#print a_data[i][2]
				#print a_p[2]

				#Ca = a_p[2]
				Cb = b_p[2]

				# we store the untransformed point, but the transformed covariance of the A point
				match_pairs.append([a_data[i],b_p,Ca,Cb])

				#findRotationTangent(pnt)

		
		#offset[2] += math.pi/16
		#newCost = cost_func(offset, match_pairs)
		#exit()

		# optimize the match error for the current list of match pairs
		newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs])

		# get the current cost
		newCost = cost_func(newOffset, match_pairs)

		# check for convergence condition, different between last and current cost is below threshold
		if abs(lastCost - newCost) < costThresh:
			offset = newOffset
			lastCost = newCost
			break

		# save the current offset and cost
		offset = newOffset
		lastCost = newCost
		numIterations += 1

		# optionally draw the position of the points in current transform
		draw(a_trans,b_data, numIterations)

	return offset

def draw(a_pnts, b_pnts, index):

	# optionally draw the position of the points in current transform
	for a in a_pnts:
		pylab.scatter([a[0]],[a[1]],linewidth=1, color='b')
		
	for b in b_pnts:
		pylab.scatter([b[0]],[b[1]],linewidth=1, color='r')

	pylab.xlim(-4.5,4.5)
	pylab.ylim(-4,4)
	pylab.savefig("ICP_alpha_%04u.png" % index)
	pylab.clf()


if __name__ == '__main__':

	# DONE:  store the covariance for each point untransformed so we don't have to repeat
	# DONE:  transform each covariance by the appropriate rotation at each point

	# 1. compute the local line for each point
	# 2. compute the covariance for point-to-line constraint
	# 3. compute the covariance for rotational error constraint for A points only
	# 4. add the covariances together

	# 5. find closest points in A of B
	# 6. discard matches beyond threshold d_max < | ai - T *bi |
	# 7. optimize T for the sum of computeMatchError
	# 8. if converged, stop, else go to 5 with b offset by new T

	#offset = [0.0,0.0,-math.pi/2]
	offset = [0.0,0.0,0.0]
	a_data = [[1.0,1.0],[1.1,1.1],[1.2,1.2],[1.3,1.31],[1.4,1.4],[1.51,1.5],[1.6,1.6]]
	b_data = [[0.3,1.0],[0.3,1.1],[0.3,1.2],[0.31,1.3],[0.3,1.4],[0.3,1.5],[0.3,1.6]]

	f = open('alpha_bound_0007.txt','r')
	#f = open('alpha_bound_0001.txt','r')
	val = f.read()
	b_data = eval(val)
	f = open('alpha_bound_0009.txt','r')
	#f = open('alpha_bound_0003.txt','r')
	val = f.read()
	a_data = eval(val)

	f = open('estpose0007.txt','r')
	#f = open('estpose0001.txt','r')
	val = f.read()
	estPose1 = eval(val)

	f = open('estpose0009.txt','r')
	#f = open('estpose0003.txt','r')
	val = f.read()
	estPose2 = eval(val)


	# get the offset of node 2 with respect to node 1
	pose1 = Pose(estPose1)
	offset = pose1.convertGlobalPoseToLocal(estPose2)


	a_trans = []
	for p in a_data:
		a_trans.append(dispOffset(p, offset))
	draw(a_data, b_data, -1)
	draw(a_trans, b_data, 0)

	offset = gen_ICP(offset, a_data, b_data)

	a_trans = []
	for p in a_data:
		a_trans.append(dispPoint(p, offset))

	draw(a_trans, b_data, 100)

