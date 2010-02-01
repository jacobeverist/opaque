#!/usr/bin/python

import gen_icp
import math
import numpy 
import random
from coord import Pose



if __name__ == '__main__':
	x_var = 0.01
	y_var = 1.0

	cov1 = gen_icp.computeVectorCovariance([1.0,0.0],x_var,y_var)
	cov2 = gen_icp.computeVectorCovariance([1.0,1.0],x_var,y_var)
	cov3 = gen_icp.computeVectorCovariance([0.0,1.0],x_var,y_var)
	#print cov1
	#print cov2
	#print cov3

	a_data = []
	for i in range(20):
		#a_data.append([i*0.1 + random.gauss(0.0,0.01),i*0.1 + random.gauss(0.0,0.01)])
		#a_data.append([1.0 + i*0.1,1.0 + i*0.1 ])
		a_data.append([i*0.1,i*0.1 ])

	norm1 = gen_icp.findLocalNormal([1.5,1.5],a_data)
	#print norm1

	a_data = []
	for i in range(20):
		a_data.append([-i*0.1 + random.gauss(0.0,0.01),i*0.1 + random.gauss(0.0,0.01)])

	norm2 = gen_icp.findLocalNormal([-1.5,1.5],a_data)
	#print norm2


	offset = [0.0,0.0,0.0]
	a = numpy.matrix([[0.0],[0.0]])
	b = numpy.matrix([[0.1],[0.0]])
	c = numpy.matrix([[0.0],[0.1]])
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	offset = [-0.1,0.0,0.0]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	offset = [0.1,0.0,0.0]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	offset = [0.0,0.0,math.pi]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	offset = [0.0,0.0,math.pi/2]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	# colinear on the axis Y but with opposing variances on Y
	offset = [0.1,0.1,math.pi/2]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	offset = [0.1,0.2,math.pi/2]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	# colinear on the high variance axis Y
	offset = [0.1,0.1, 0.0]
	sum = gen_icp.computeMatchError(offset, a, b, cov1, cov1)
	#print sum

	offset = [0.0,0.0, 0.0]
	sum = gen_icp.computeMatchError(offset, a, c, cov1, cov1)
	#print sum


	f = open('alpha_bound_0007.txt','r')
	val = f.read()
	b_data = eval(val)

	f = open('alpha_bound_0009.txt','r')
	val = f.read()
	a_data = eval(val)

	f = open('estpose0007.txt','r')
	val = f.read()
	estPose1 = eval(val)

	f = open('estpose0009.txt','r')
	val = f.read()
	estPose2 = eval(val)


	pose1 = Pose(estPose1)
	pose2 = Pose(estPose2)

	global1 = pose2.convertLocalToGlobal(a_data[10])
	local1 = pose1.convertGlobalToLocal(global1)

	global2 = pose2.convertLocalToGlobal(a_data[300])
	local2 = pose1.convertGlobalToLocal(global2)

	relPose2 = pose1.convertGlobalPoseToLocal(estPose2)

	print "estPose1:", estPose1
	print "estPose2:", estPose2
	print "relPose2:", relPose2
	print
	print "local point of pose2:", a_data[10]
	print "global point of pose2:", global1
	print "point with respect to pose1:", local1
	print
	print "local point of pose2:", a_data[300]
	print "global point of pose2:", global2
	print "point with respect to pose1:", local2
	
	print
	print


	#offset = [estPose2[0]-estPose1[0],estPose2[1]-estPose2[1],-gen_icp.normalizeAngle(estPose2[2]-estPose1[2])]
	offset = relPose2
	print gen_icp.dispOffset(a_data[10], offset)
	print gen_icp.dispOffset(a_data[300], offset)

	# 1. problem is that the ICP uses homogeneous transformations of points in the local coordinate systmes,
	# 2. however, the legacy system uses a more complex transformation of coordinate systems in global coordinate systems

	# 3. so the offset in ICP means something different than the offset in legacy coordinates

	# Test I can do.... 
	# 1. get estpose of 2 nodes.
	# 2. compute relative position of point in node 2 with respect to node 1
	# 3. use legacy coordinate system
	# 4. try to get ICP to do the same thing
	# Hope... leads to a conversion

	control_points = [[0.0,0.0],[1.0,1.0],[2.0,2.0],[3.01,2.9],[4.0,4.0]]
	print offset
	#print gen_icp.modifyGuess(estPose1, offset, control_points)


	print [0.0,0.0]
	result = gen_icp.dispOffset([0.0,0.0], offset)	
	print result
	result = gen_icp.dispOffset(result, [-offset[0],-offset[1],0.0])	
	print result





	

