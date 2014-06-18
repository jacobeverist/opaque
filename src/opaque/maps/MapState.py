

from Pose import Pose
import math
import ctypes, os, sys
import multiprocessing as processing
from copy import copy, deepcopy
from functions import *
import random
from subprocess import Popen, PIPE
import Image
from medialaxis import computeMedialAxis
import graph
from LocalNode import getLongestPath
import gen_icp
from SplineFit import SplineFit
import pylab
import numpy
from operator import itemgetter
import hashlib
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath, getTipAngles
from MapProcess import selectLocalCommonOrigin, selectCommonOrigin
from ParticleFilter import multiParticleFitSplice, batchLocalizeParticle, batchDisplaceParticles, Particle
import time
import traceback

import alphamod

from shoots import computeShootSkeleton, spliceSkeletons, computeGlobalControlPoses

pylab.ioff()


" 1) probability that a path is the same as another path "
" 2) new ID for each new path "
" 3) test fitness of being its own path and being the same as other path "
" 4) collapse path to old path if pass threshold "

"""
Path Data Structure

parentID: int, the parent path ID
branchNodeID: int, the node from which the branching point is computed
localJunctionPose: [float*3], the pose of the branching point in local coordinates of branchNodeID
sameProb: dict of floats for each pathID indicating probability this path is the same as older path 

"""

" Likeness Test "
" 1) % overlap of two paths via ICP "
" 2) successful overlap of member nodes translated to old path "
" 2a) maintain separate set of intra-path constraints for each path "
" 2b) measure X_theta error for each path hypothesis "

" New Path Condition "
" 1) any new departure automatically creates a new path "
" 2) a new departure that does not overlap correctly creates a new path "


" Overlap State "
" 1) perform best overlap of each splice "
" 2) perform ICP of each splice "
" 3) maintain path-state machine to keep track of overlaps and directions and likely configurations "
" 4) extract path-state from particle filter motion estimation and subsequent ICP fit on splice "


" FIXME:  Branch direction does not seem to be used properly yet.  Only used in getOverlapDeparture() "

"""

Things that are happening with the data set with branches 
-----------

Second branch detection fails because is deemed "similar" in departure angle to path 1.
The difference is about 1.01 radians but is just under the difference threshold.
Paradoxically, this node is not added to the similar path and is instead added to path 0.
This is because the similar path 1 is not part of the orderedOverlap set.

When the next instance of branch detection occurs, it is determined to be successful because of both dist and ang diff.
However, this new branch eventually is constrained to the junction where it actually belongs and on top of the node 
with the previous branch detection fail.  The mechanism of how this happens is?   NOTE:  This happened because of a bug in my code
that allowed cross-path constraints between member nodes (node path membership always returned True).  Remarkably, there were no map errors from this permissiveness.

In the event that the first node is successfully branch detected, subsequent branch nodes are added to existing path,
but not properly merged or constrained.

"""

" 1) splice is signified by the paired terminals (t0,t1) "
" 2) splice can only change by one terminal.  (t0,t1) -> (t0,t2) "
" assume at most 3 paths for possible splice states "

" junctions indexed by their pathID "
" terminals indexed by pathID+1,  root is 0 and 1 "
" t%u % (pathID+1) dictionary returns pathID and index "
" j%u % (pathID) dictionary returns pathID and index "

" 3) splice state for node0 and node1 are same or separate by one terminal delta "
" 4) splice state for node2 and node3 are at most one delta away from node0 and node1 "
" 5) generate possible splice states for node2 given splice state of node0 "
" 6) apply displacement of node2 from node0's position for every given splice "
" 7) compute its fit according to getMultiDeparture(), classifying by matchCount, cost, and contiguity (or ICP?) (see selectSplice function)"
" 8) select its most likely splice state "
" 9) repeat for node1 and node3 "
" 10) difference of splice state between node2 and node3 is at most 1, but usually zero "

			
" FIXME: branch from set of ordered overlap pathIDs (0,2) selects 2 as the parent, but should really be 0 "
" How did it overlap with 2 in the first place?  It departed from 0, but selected 2 as parent since it's the terminating pathID "
" 2 was marked positive by overlapCondition since it was parallel to 2 within the minMatchDist of 1.0 "

renderGlobalPlotCount = 0
pathPlotCount = 0

qin_branch = None
qout_branch =  None
pool_branch = []


def __num_processors():
	
	return 4
	
	if os.name == 'nt': # Windows
		return int(os.getenv('NUMBER_OF_PROCESSORS'))
	else: # glibc (Linux, *BSD, Apple)
		get_nprocs = ctypes.cdll.libc.get_nprocs
		get_nprocs.restype = ctypes.c_int
		get_nprocs.argtypes = []
		return get_nprocs()

def __remote_prof_multiBranch(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_multiBranch(rank, qin, qout)', "branch_%d.prof" % pid )
		cProfile.runctx("__remote_multiBranch(rank, qin, qout)", globals(), locals(), "branch_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


def __remote_multiBranch(rank, qin, qout):

	sys.stdout = open("branch_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("branch_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

	print "started __remote_multiBranch"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		# process data
		results = []
		for job in args:

			#branchJobs.append((pathID, parentID, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs, arcDist))

			pathID = job[0]
			parentID = job[1]
			origGlobJuncPose = job[2]
			controlPose = job[3]
			childPath = job[4]
			parentPath = job[5]
			trimmedParent = job[6]
			localPathSegs = job[7]
			arcDist = job[8]
			localSkeletons = job[9]
			controlPoses = job[10]
			junctionPoses = job[11]
			parentPathIDs = job[12]
			numNodes = job[13]
			hypID = job[14]

			result = computeBranch(pathID, parentID, origGlobJuncPose, controlPose, childPath, parentPath, trimmedParent, localPathSegs, arcDist, localSkeletons, controlPoses, junctionPoses, parentPathIDs, numNodes, hypID)

			results.append(result)
						   
		# write to output queue
		qout.put((nc,results))

def batchBranch(branchJobs):

	global renderGlobalPlotCount
	global pool_branch
	global qin_branch
	global qout_branch

	ndata = len(branchJobs)
	
	args = []
	
	for k in range(ndata):
		arg = branchJobs[k]
		args.append(arg)
	
	nproc = __num_processors()
	
	nproc *= 2
	print "nproc =", nproc
	
	# compute chunk size
	chunk_size = ndata / nproc
	chunk_size = 2 if chunk_size < 2 else chunk_size
	
	print "chunk_size =", chunk_size
	print "max_size =", ndata/chunk_size
	
	# set up a pool of processes
	if len(pool_branch) == 0:

		qin_branch = processing.Queue(maxsize=ndata/chunk_size)
		qout_branch = processing.Queue(maxsize=ndata/chunk_size)
		pool_branch = [processing.Process(target=__remote_multiBranch,
		#pool_branch = [processing.Process(target=__remote_prof_multiBranch,
				args=(rank, qin_branch, qout_branch))
					for rank in range(nproc)]
		for p in pool_branch: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	#print "args =", args
	while 1:
		#_data = data[:,cur:cur+chunk_size]
		_data = args[cur:cur+chunk_size]
		#print "_data =", _data
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		#print "put(", (nc,_data), ")"
		qin_branch.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout_branch.get()]
	
	print "received output:", knn
		
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "sorted output:", knn
		

	#print "terminating pool"
	# terminate workers
	#for p in pool_branch:
	#	print "terminate"
	#	p.terminate()
	#	print "terminated"
		
	print "returning"
	return knn

def printStack():

	flist = traceback.format_stack()
	flist = flist[:-1]
	
	printStr = ""
	for line in flist:
		printStr += line
		
	print printStr

@logFunction
def getTangentIntersections(path1, path2, frontDepI, backDepI, path1FrontDepI, path1BackDepI, path2JuncI, path1JuncI, plotCount, hypothesisID = 0, nodeID = 0, plotIter = False):


	""" for tangent on path2, find the intersection point on path1 """

	interPoints = []
	indices1 = []
	foreEdges = []
	backEdges = []

	forePoints = []
	backPoints = []

	frontDepI = path2JuncI
	backDepI = path2JuncI

	frontBoundI = frontDepI - 40
	if frontBoundI < 0:
		frontBoundI = 0

	backBoundI = backDepI + 40
	if backBoundI >= len(path2):
		backBoundI = len(path2)-1

	frontDepI = frontDepI - 10
	if frontDepI < 0:
		frontDepI = 0

	backDepI = backDepI + 10
	if backDepI >= len(path2):
		backDepI = len(path2)-1

	print "frontBoundI, frontDepI:", frontBoundI, frontDepI
	print "backBoundI, backDepI:", backBoundI, backDepI

	for i in range(frontBoundI, frontDepI+1):

		p = path2[i]

		""" generate tangent segment """
		angle2 = p[2]

		pA = [p[0] + 3*cos(angle2), p[1] + 3*sin(angle2)]
		pB = [p[0] - 3*cos(angle2), p[1] - 3*sin(angle2)]

		edge2 = deepcopy([pA,pB])
		foreEdges.append(edge2)

		""" find the intersection with path1 """

		for k in range(len(path1)-1):
			pathEdge = [path1[k],path1[k+1]]
			isIntersect1, point1 = Intersect(edge2, pathEdge)
			if isIntersect1:

				vec = [p[0]-point1[0],p[1]-point1[1]]
				mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1])

				vec[0] /= mag
				vec[1] /= mag
				
				tanAng = acos(vec[0])
				if vec[1] < 0.0:
					tanAng = -tanAng

				""" save point and index on path1 """
				interPoints.append(point1 + [tanAng])
				indices1.append(k)

				if k <= path1FrontDepI:
					forePoints.append(point1 + [tanAng])

	for i in range(backDepI, backBoundI+1):

		p = path2[i]

		""" generate tangent segment """
		angle2 = p[2]

		print i, ":", angle2, p[:2]

		pA = [p[0] + 3*cos(angle2), p[1] + 3*sin(angle2)]
		pB = [p[0] - 3*cos(angle2), p[1] - 3*sin(angle2)]

		edge2 = deepcopy([pA,pB])
		#print "back edge:", p, edge2
		backEdges.append(edge2)

		""" find the intersection with path1 """

		for k in range(len(path1)-1):
			pathEdge = [path1[k],path1[k+1]]
			isIntersect1, point1 = Intersect(edge2, pathEdge)
			if isIntersect1:

				vec = [p[0]-point1[0],p[1]-point1[1]]
				mag = sqrt(vec[0]*vec[0] + vec[1]*vec[1])

				vec[0] /= mag
				vec[1] /= mag
				
				tanAng = acos(vec[0])
				if vec[1] < 0.0:
					tanAng = -tanAng

				""" save point and index on path1 """
				interPoints.append(point1 + [tanAng])
				indices1.append(k)

				if k >= path1BackDepI:
					backPoints.append(point1 + [tanAng])

	foreDists = [0.0 for k in range(len(forePoints))]
	if len(forePoints) > 1:
		for k in range(0,len(forePoints)):
			if k == 0:
				dist1 = sqrt((forePoints[k][0]-forePoints[k+1][0])**2 + (forePoints[k][1]-forePoints[k+1][1])**2)
				totalDist = 2*dist1

			elif k == len(forePoints)-1:
				dist2 = sqrt((forePoints[k][0]-forePoints[k-1][0])**2 + (forePoints[k][1]-forePoints[k-1][1])**2)
				totalDist = 2*dist2

			else:

				dist1 = sqrt((forePoints[k][0]-forePoints[k+1][0])**2 + (forePoints[k][1]-forePoints[k+1][1])**2)
				dist2 = sqrt((forePoints[k][0]-forePoints[k-1][0])**2 + (forePoints[k][1]-forePoints[k-1][1])**2)

				totalDist = dist1 + dist2


			foreDists[k] = totalDist

	maxForeDist = 0.0
	minForeDist = 1e100
	minForeI = 0
	newForeDists = [0.0 for k in range(len(forePoints))]
	for k in range(0,len(forePoints)):

		leftI = k-2
		rightI = k+2

		if leftI < 0:
			leftI = 0

		if rightI > len(forePoints)-1:
			rightI = len(forePoints)-1

		leftDist = foreDists[leftI]
		rightDist = foreDists[rightI]

		totalDist = leftDist + rightDist + foreDists[k]
		newForeDists[k] = totalDist

		if totalDist > maxForeDist:
			maxForeDist = totalDist

		if totalDist < minForeDist:
			minForeDist = totalDist
			minForeI = k


	foreDists = newForeDists

	backDists = [0.0 for k in range(len(backPoints))]
	if len(backPoints) > 1:
		for k in range(0,len(backPoints)):
			if k == 0:
				dist1 = sqrt((backPoints[k][0]-backPoints[k+1][0])**2 + (backPoints[k][1]-backPoints[k+1][1])**2)
				totalDist = 2*dist1

			elif k == len(backPoints)-1:
				dist2 = sqrt((backPoints[k][0]-backPoints[k-1][0])**2 + (backPoints[k][1]-backPoints[k-1][1])**2)
				totalDist = 2*dist2

			else:

				dist1 = sqrt((backPoints[k][0]-backPoints[k+1][0])**2 + (backPoints[k][1]-backPoints[k+1][1])**2)
				dist2 = sqrt((backPoints[k][0]-backPoints[k-1][0])**2 + (backPoints[k][1]-backPoints[k-1][1])**2)

				totalDist = dist1 + dist2


			backDists[k] = totalDist


	maxBackDist = 0.0
	minBackDist = 1e100
	minBackI = 0
	newBackDists = [0.0 for k in range(len(backPoints))]
	for k in range(0,len(backPoints)):

		leftI = k-2
		rightI = k+2

		if leftI < 0:
			leftI = 0

		if rightI > len(backPoints)-1:
			rightI = len(backPoints)-1

		leftDist = backDists[leftI]
		rightDist = backDists[rightI]

		totalDist = leftDist + rightDist + backDists[k]
		newBackDists[k] = totalDist

		if totalDist > maxBackDist:
			maxBackDist = totalDist

		if totalDist < minBackDist:
			minBackDist = totalDist
			minBackI = k

	backDists = newBackDists

	print "foreDists:", foreDists
	print "backDists:", backDists
	#print "interPoints:", interPoints



	foreAvg = [0.0,0.0,0.0]
	foreWeight = 0.0
	for k in range(len(forePoints)):
		if maxForeDist > 0.0:
			foreAvg[0] += forePoints[k][0] * (1.0 - foreDists[k]/maxForeDist)
			foreAvg[1] += forePoints[k][1] * (1.0 - foreDists[k]/maxForeDist)
			foreAvg[2] += forePoints[k][2] * (1.0 - foreDists[k]/maxForeDist)
			foreWeight += (1.0 - foreDists[k]/maxForeDist)
		else:
			foreAvg[0] += forePoints[k][0]
			foreAvg[1] += forePoints[k][1]
			foreAvg[2] += forePoints[k][2]
			foreWeight += 1.0 


	if len(forePoints) > 0:

		if foreWeight >= 0:
			foreAvg[0] /= foreWeight
			foreAvg[1] /= foreWeight
			foreAvg[2] /= foreWeight
			foreAvg[2] = normalizeAngle(foreAvg[2])

		tempAvg = copy(forePoints[minForeI])
		tempAvg[2] = foreAvg[2]
		foreAvg = tempAvg

	else:

		p_f = path1[path1FrontDepI] 
		if len(forePoints) > 0:
			a_f, i_af, minDist_af = gen_icp.findClosestPointInA(forePoints, p_f)
			juncForeAng = forePoints[i_af][2]
		elif len(interPoints) > 0:
			a_f, i_af, minDist_af = gen_icp.findClosestPointInA(interPoints, p_f)
			juncForeAng = interPoints[i_af][2]
		else:
			juncForeAng = 0.0

		foreAvg = copy(path1[path1FrontDepI])
		foreAvg[2] = juncForeAng

	backAvg = [0.0,0.0,0.0]
	backWeight = 0.0
	for k in range(len(backPoints)):
		if maxBackDist > 0.0:
			backAvg[0] += backPoints[k][0] * (1.0 - backDists[k]/maxBackDist)
			backAvg[1] += backPoints[k][1] * (1.0 - backDists[k]/maxBackDist)
			backAvg[2] += backPoints[k][2] * (1.0 - backDists[k]/maxBackDist)
			backWeight += (1.0 - backDists[k]/maxBackDist)
		else:
			backAvg[0] += backPoints[k][0]
			backAvg[1] += backPoints[k][1]
			backAvg[2] += backPoints[k][2]
			backWeight += 1.0 

	if len(backPoints) > 0:

		if backWeight > 0:
			backAvg[0] /= backWeight
			backAvg[1] /= backWeight
			backAvg[2] /= backWeight
			backAvg[2] = normalizeAngle(backAvg[2])

		tempAvg = copy(backPoints[minBackI])
		tempAvg[2] = backAvg[2]
		backAvg = tempAvg

	else:

		p_b = path1[path1BackDepI]

		if len(backPoints) > 0:
			a_b, i_ab, minDist_ab = gen_icp.findClosestPointInA(backPoints, p_b)
			juncBackAng = backPoints[i_ab][2]
		elif len(interPoints) > 0:
			a_b, i_ab, minDist_ab = gen_icp.findClosestPointInA(interPoints, p_b)
			juncBackAng = interPoints[i_ab][2]
		else:
			juncBackAng = 0.0

		backAvg = copy(path1[path1BackDepI])
		backAvg[2] = juncBackAng


	juncForeAng1 = foreAvg[2]
	juncBackAng1 = backAvg[2]


	p_f, i_f, minDist_f = gen_icp.findClosestPointInA(path1, foreAvg)
	p_b, i_b, minDist_b = gen_icp.findClosestPointInA(path1, backAvg)

	if len(forePoints) > 0:
		a_f, i_af, minDist_af = gen_icp.findClosestPointInA(forePoints, p_f)
		juncForeAng = forePoints[i_af][2]
	elif len(interPoints) > 0:
		a_f, i_af, minDist_af = gen_icp.findClosestPointInA(interPoints, p_f)
		juncForeAng = interPoints[i_af][2]
	else:
		juncForeAng = 0.0

	if len(backPoints) > 0:
		a_b, i_ab, minDist_ab = gen_icp.findClosestPointInA(backPoints, p_b)
		juncBackAng = backPoints[i_ab][2]
	elif len(interPoints) > 0:
		a_b, i_ab, minDist_ab = gen_icp.findClosestPointInA(interPoints, p_b)
		juncBackAng = interPoints[i_ab][2]
	else:
		juncBackAng = 0.0



	if plotIter:
		pylab.clf()
		xP = []
		yP = []
		for p in path2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.5,0.5,1.0),zorder=500)

		if True:	

			""" draw the tangents """
			for edge in foreEdges:
				xP = [edge[0][0], edge[1][0]]
				yP = [edge[0][1], edge[1][1]]
				pylab.plot(xP,yP, color=(0.5,1.0,0.5), linewidth=1, alpha=0.5,zorder=1)		   

			for edge in backEdges:
				xP = [edge[0][0], edge[1][0]]
				yP = [edge[0][1], edge[1][1]]
				pylab.plot(xP,yP, color=(1.0,0.5,0.5), linewidth=1, alpha=0.5,zorder=1)		   

			for pnt in forePoints + backPoints:
				xP = [pnt[0],]
				yP = [pnt[1],]
				pylab.scatter(xP,yP,color='k', linewidth=1, zorder=502)

			pylab.scatter([foreAvg[0],backAvg[0]], [foreAvg[1],backAvg[1]], color='m', linewidth=1, zorder=503)
			pylab.scatter([p_f[0],], [p_f[1],], color='r', linewidth=1, zorder=504)
			pylab.scatter([p_b[0],], [p_b[1],], color='g', linewidth=1, zorder=504)

			pylab.scatter([path2[path2JuncI][0],], [path2[path2JuncI][1],], color='y', linewidth=1, zorder=504)

		xP = []
		yP = []
		for p in path1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(1.0,0.5,0.5),zorder=500)

		print "intersectDeparture:", plotCount
		printStack()				 

		pylab.title("%d intersections, %d %d tangent segments, angles %1.2f %1.2f %1.2f %1.2f %d %d" % (len(interPoints), len(foreEdges), len(backEdges), juncForeAng1, juncBackAng1, juncForeAng, juncBackAng, hypothesisID, nodeID))
		pylab.savefig("intersectDeparture_%04u_%04u_%04u.png" % (nodeID, hypothesisID, plotCount))
		

		return i_f, i_b, juncForeAng, juncBackAng


@logFunction
def computePathAngleVariance(pathSamples):
	pathVar = []
	
	" compute the local variance of the angle "
	VAR_WIDTH = 40
	for i in range(len(pathSamples)):
		
		lowK = i - VAR_WIDTH/2
		if lowK < 0:
			lowK = 0
			
		highK = i + VAR_WIDTH/2
		if highK >= len(pathSamples):
			highK = len(pathSamples)-1
		
		localSamp = []
		for k in range(lowK, highK+1):
			localSamp.append(pathSamples[k][2])
		
		sum = 0.0
		for val in localSamp:
			sum += val
			
		meanSamp = sum / float(len(localSamp))
		

		sum = 0
		for k in range(len(localSamp)):
			sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
	
		varSamp = sum / float(len(localSamp))
		
		pathVar.append((meanSamp, varSamp))		 

	return pathVar

" get the trimmed version of child and parent paths that are overlapping in some fashion "
@logFunction
def getOverlapDeparture(globalJunctionPose, parentPathID, childPathID, path1, path2, plotIter = False):

	"Assumption:  one section of the medial axis is closely aligned with the path "		   
		
	print "getOverlapDeparture():"
	
	isExist1 = False
	isInterior1 = False
	departurePoint1 = 0
	isExist2 = False
	isInterior2 = False
	departurePoint2 = 0

	" orienting the medial axis of the branch node correctly "
	
	" return exception if we receive an invalid path "		  
	if len(path1) == 0:
		print "path1 has zero length"
		raise


	" make sure the overlap of both paths are oriented the same way "
	path1Spline = SplineFit(path1, smooth=0.1)
	path2Spline = SplineFit(path2, smooth=0.1)

	path2Reverse = deepcopy(path2)
	path2Reverse.reverse()
	path2SplineReverse = SplineFit(path2Reverse, smooth=0.1)
	
	orientedPath2 = orientPath(path2, path1, dist_thresh=0.1)
		
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
	
	"TODO:	why is there a 3-offset in this comparison? "
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


	" NOTE:  This guarantees detection of a departure.	We've assumed that there exists a"
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
	
	angleSum2 = 0.0
	overlapSum2 = 0.0
	matchCount2 = 0
	
	" path section for our back departure hypothesis "
	pathSec2 = pathPoints2[backDepI:]

	print "pathSec1 hypothesis angle and overlap sum and match count:", angleSum1, overlapSum1, matchCount1
	print "pathSec2 hypothesis angle and overlap sum and match count:", angleSum2, overlapSum2, matchCount2

	" distance of departure point from known junction point "
	p0 = pathSec1[0]
	juncDist1 = sqrt((globalJunctionPose[0]-p0[0])**2 + (globalJunctionPose[1]-p0[1])**2)

	p0 = pathSec2[0]
	juncDist2 = sqrt((globalJunctionPose[0]-p0[0])**2 + (globalJunctionPose[1]-p0[1])**2)

	" juncAngle "
	ang1 = pathSec1[0][2]
	ang2 = pathSec2[0][2]

	frontVec = [pathPoints2[frontDepI][0] - pathPoints2[frontDepI+1][0], pathPoints2[frontDepI][0] - pathPoints2[frontDepI+1][0]]
	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	frontVec[0] /= frontMag
	frontVec[1] /= frontMag

	frontJuncAng = acos(-frontVec[0])
	if -frontVec[1] < 0.0:
		frontJuncAng = -frontJuncAng

	backVec = [pathPoints2[backDepI][0] - pathPoints2[backDepI-1][0], pathPoints2[backDepI][0] - pathPoints2[backDepI-1][0]]
	#backVec = [pathSec2[1][0] - pathSec2[0][0], pathSec2[1][0] - pathSec2[0][0]]
	backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
	backVec[0] /= backMag
	backVec[1] /= backMag

	backJuncAng = acos(-backVec[0])
	if -backVec[1] < 0.0:
		backJuncAng = -backJuncAng

	#ang1 = diffAngle(pathSec1[0][2],math.pi)
	#ang2 = pathSec2[0][2]
	#secP2 = pathPoints2[frontDepI]
	#secP1 = pathPoints2[backDepI]

	print "hypothesis discrepancy distance:", juncDist1, juncDist2, globalJunctionPose[2], diffAngle(frontJuncAng, globalJunctionPose[2]), diffAngle(backJuncAng, globalJunctionPose[2])


	#if plotIter:
	if False:

		hypothesisID = 0
		numNodes = 0
		pathPlotCount = 0

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

			
		pylab.scatter([globalJunctionPose[0]],[globalJunctionPose[1]], color='r')		   

		xP = []
		yP = []
		for p in path1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(1.0,0.5,0.5))

		pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		pylab.title("hyp %d nodeID %d, %d %d %d %d %d %3.2f %3.2f %3.2f %d %3.2f %3.2f %3.2f" % ( hypothesisID, numNodes, isExist1, isExist2, isInterior1, isInterior2, matchCount1, overlapSum1, angleSum1, juncDist1, matchCount2, overlapSum2, angleSum2, juncDist2))
		pylab.savefig("trimDeparture_%04u_%04u.png" % (hypothesisID, pathPlotCount))

		print "saving trimDeparture_%04u_%04u.png" % (hypothesisID, pathPlotCount)
		
		pathPlotCount += 1

	secP1 = []
	secP2 = []

	"""
	cases
	1) near-zero length of path section, never departs, not the path, high juncDist
	2) near-zero length, no match count, low juncDist, correct path
	3) high path length, high match count, low juncDist, incorrect path 
	4) no match count, high junc dist, incorrect
	5) no match count, low junc dist, correct
	"""
	
	if juncDist1 < juncDist2:
		
		" FIXME: if [0,frontDepI] has high overlap compared to total length, then we reject it "
		
		secP1 = pathPoints2[0]
		secP2 = pathPoints2[frontDepI]
	else:
		" FIXME: if [backDepI,len(distances)-1] has high overlap compared to total length, then we reject it "
		secP1 = pathPoints2[backDepI]
		secP2 = pathPoints2[len(distances)-1]

	if len(secP1) == 0:
		print "no departures found"
		raise
		
	
	return secP1, secP2


@logFunction
def getOverlapCondition(medial2, estPose2, supportLine, nodeID, plotIter = False, overlapPlotCount = 0):

	if len(supportLine) == 0:
		return 1e100

	#medial2 = self.poseData.medialAxes[nodeID]

	#estPose2 = self.nodePoses[nodeID]
	
	minMatchDist2 = 0.2
		
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
	
	support_pairs = []
	for i in range(len(vecPoints2)):
		
		p_2 = vecPoints2[i]

		" for every transformed point of A, find it's closest neighbor in B "
		try:
			p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)

			if minDist <= minMatchDist2:
				C2 = points2[i][2]
				C1 = supportPoints[minI][2]

				" we store the untransformed point, but the transformed covariance of the A point "
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
		pylab.savefig("overlapCost_%04u.png" % overlapPlotCount)
		#overlapPlotCount += 1
	
	if len(support_pairs) == 0:
		return 1e100

	return cost


@logFunction
def getAngleDerivatives(pathPoints1):

	""" compute the angle derivative at each point of the curves """

	angDerivs1 = []
	for k in range(len(pathPoints1)):
		
		if (k-6) > 0 and (k+6) < len(pathPoints1):

			diffs = []

			for l in range(0,7):

				if (k-6+l) >= 0 and (k+l) < len(pathPoints1):
					ang1 = pathPoints1[k-6+l][2]
					ang2 = pathPoints1[k+l][2]
					diffs.append(fabs(ang1-ang2))

			angDeriv = sum(diffs) / float(len(diffs))

			angDerivs1.append(angDeriv)
		else:
			angDerivs1.append(0.0)

	return angDerivs1

@logFunction
def getClosestPairs(pathPoints2, pathPoints1):

	""" for each point on the child path, find its closest pair on the parent path """
	distances = []
	indices = []

	for i in range(0,len(pathPoints2)):
		p_2 = pathPoints2[i]
		p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
		
		""" keep the distance information """
		distances.append(minDist)
		""" and the associated index of the point on the parent path """
		indices.append(i_1)
		
	return distances, indices

@logFunction
def getPointDistances(pathPoints2, globalJunctionPose):

	juncDists = []
	minDist2 = 1e100
	juncI = 0		 

	for i in range(0,len(pathPoints2)):
		p_2 = pathPoints2[i]
		
		juncDist = sqrt((p_2[0]-globalJunctionPose[0])**2 + (p_2[1]-globalJunctionPose[1])**2)

		if juncDist < minDist2:
			minDist2 = juncDist
			juncI = i

		juncDists.append(juncDist)

	return juncDists, minDist2, juncI

@logFunction
def computeDivergencePoint(juncI, frontInd, distances, indices, pathPoints1, pathPoints2):

	""" initial distances of the indexed point on the child shoot, here classified as
		the first point within radius of the temporary junction point

		pathPoints2 are aligned with distances and indices
	"""		   
	maxFront = distances[frontInd]
	#maxBack = distances[backInd]
	print "match distances of tip points:", maxFront#, maxBack

	""" walk back from first point until we have a non-monotic increase in match distance """
	""" this becomes our departure point """
	
	INDEX_DELTA = 3
	""" FIXME:  index out of bounds case, occurs when the number of shoot points is smaller than various deltas """

	""" the 3-offset is so that the any small blips do not affect the distance increase threshold """
	currI = frontInd + 1
	try:
		while distances[currI+INDEX_DELTA] < maxFront or distances[currI+INDEX_DELTA] > 0.1:
			maxFront = distances[currI]
			currI += 1
	except:
		pass
	
	""" departure index on child shoot """
	frontDepI = currI

	if frontDepI+3 >= len(pathPoints2)-1:
		frontDepI = juncI

	
	""" departure index of child path and match distance """
	frontPoint = [frontDepI, distances[frontDepI]]


	""" starting from the frontDepI, reverse direction and step until the distance is greater than 0.1 """
	frontAngleRefI = frontDepI
	while distances[frontAngleRefI] < 0.1:
		frontAngleRefI -= 1
		
		if frontAngleRefI < 0:
			frontAngleRefI = 0
			break
					
	""" corresponding point on the parent shoot, the angle on parent becomes our reference angle """
	forePathIndex = indices[frontAngleRefI]
	forePathAngle = pathPoints1[forePathIndex][2]
	forePathAngle = normalizeAngle(forePathAngle)
	

	""" 40 steps away from initial departure point on child shoot, 
		step back towards the parent until the angle difference is
		< pi/3
	"""
	newFrontDepI = frontDepI - 40
	if newFrontDepI < 4:
		newFrontDepI = 4

	while fabs(diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]), forePathAngle)) > pi/3.0:
		
		if newFrontDepI < frontDepI:
			newFrontDepI += 1
		else:
			break
	
	""" compute the angle difference from the parent shoot reference angle """
	frontDiffAngles = []
	for k in range(len(pathPoints2)):
		if k == newFrontDepI:
			frontDiffAngles.append(None)
		else:
			frontDiffAngles.append(fabs(diffAngle(normalizeAngle(pathPoints2[k][2]), forePathAngle)))

	#print "frontDepI, backDepI:", frontDepI, backDepI
	#print "frontAngleRefI, backAngleRefI:", frontAngleRefI, backAngleRefI
	#print "forePathAngle, backPathAngle:", forePathAngle, backPathAngle
	#print "newFrontDepI, newBackDepI:", newFrontDepI, newBackDepI
	#print "foreDiff, backDiff:", diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]+pi), forePathAngle), diffAngle(normalizeAngle(pathPoints2[newBackDepI][2]), backPathAngle)
	#print "frontDiffAngles:", frontDiffAngles
	#print "backDiffAngles:", backDiffAngles
	#print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint

	return frontDepI, newFrontDepI, forePathIndex, maxFront, frontDiffAngles

	"""
	maxFront = distances[frontInd]
	maxBack = distances[backInd]
	print "match distances of tip points:", maxFront, maxBack



	INDEX_DELTA = 3

	currI = frontInd + 1
	try:
		while distances[currI+INDEX_DELTA] < maxFront or distances[currI+INDEX_DELTA] > 0.1:
			maxFront = distances[currI]
			currI += 1
	except:
		pass
	
	frontDepI = currI

	if frontDepI+3 >= len(pathPoints2)-1:
		frontDepI = juncI

	
	frontPoint = [frontDepI, distances[frontDepI]]


	frontAngleRefI = frontDepI
	while distances[frontAngleRefI] < 0.1:
		frontAngleRefI -= 1
		
		if frontAngleRefI < 0:
			frontAngleRefI = 0
			break
					
	forePathIndex = indices[frontAngleRefI]
	forePathAngle = pathPoints1[forePathIndex][2]
	forePathAngle = normalizeAngle(forePathAngle + pi)
	

	newFrontDepI = frontDepI - 40
	if newFrontDepI < 4:
		newFrontDepI = 4

	while fabs(diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]+pi), forePathAngle)) > pi/3.0:
		
		if newFrontDepI < frontDepI:
			newFrontDepI += 1
		else:
			break
	


	currI = 1 + len(pathPoints2) - backInd
	try:
		while distances[-currI-INDEX_DELTA] < maxBack or distances[-currI-INDEX_DELTA] > 0.1:
			maxBack = distances[-currI]
			currI += 1
	except:
		pass

	backDepI = len(distances) - currI

	if backDepI-3 <= 0: 
		backDepI = juncI

	backPoint = [backDepI, distances[backDepI]]


	backAngleRefI = backDepI
	while distances[backAngleRefI] < 0.1:
		backAngleRefI += 1
		
		if backAngleRefI >= len(distances):
			backAngleRefI = len(distances)-1
			break
	
	backPathIndex = indices[backAngleRefI]
	backPathAngle = pathPoints1[backPathIndex][2]
	


	newBackDepI = backDepI + 40
	if newBackDepI > len(distances)-5:
		newBackDepI = len(distances)-5

	while fabs(diffAngle(pathPoints2[newBackDepI][2], backPathAngle)) > pi/3.0:
		
		if newBackDepI > backDepI:
			newBackDepI -= 1
		else:
			break		 
	

	frontDiffAngles = []
	for k in range(len(pathPoints2)):
		if k == newFrontDepI:
			frontDiffAngles.append(None)
		else:
			frontDiffAngles.append(fabs(diffAngle(normalizeAngle(pathPoints2[k][2]+pi), forePathAngle)))

	backDiffAngles = []
	for k in range(len(pathPoints2)):
		if k == newBackDepI:
			backDiffAngles.append(None)
		else:
			backDiffAngles.append(fabs(diffAngle(normalizeAngle(pathPoints2[k][2]), backPathAngle)))
	"""

@logFunction
def ensureEnoughPoints(newPath2, max_spacing = 0.08, minPoints = 5):

	if len(newPath2) < 2:
		print len(newPath2), "not enough points to expand with", newPath2
		raise

	max_spacing = 0.08
	newPath3 = []

	""" make sure path is greater than 5 points """
	while len(newPath3) <= minPoints:

		max_spacing /= 2
		print "max_spacing =", max_spacing
		
		newPath3 = [copy(newPath2[0])]
							
		for i in range(len(newPath2)-1):
			p0 = newPath2[i]
			p1 = newPath2[(i+1)]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			vec[0] /= dist
			vec[1] /= dist
			
			
			if dist > max_spacing:
				""" cut into pieces max_spacing length or less """
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					newPath3.append(newP)

			newPath3.append(copy(p1))			 


	return newPath3
	

@logFunction
def getBranchPoint(globalJunctionPose, parentPathID, childPathID, path1, path2, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):
	""" get the trimmed version of child and parent paths that are overlapping in some fashion """

	"""Assumption:  one section of the child shoot is overlapped and closely aligned with the parent shoot """
	
	global pathPlotCount

	print "getBranchPoint():"

	print "lengths of parent and child paths:", len(path1), len(path2)
	
	""" return exception if we receive an invalid path """		  
	if len(path1) == 0:
		print "path1 has zero length"
		raise

	if len(path2) == 0:
		print "path2 has zero length"
		raise
	
	""" make sure the overlap of both shoots are oriented the same way """
	orientedPath2 = orientPath(path2, path1)
	
	""" compute the control point """
	commonU1, commonU2, commonP1, commonP2 = selectCommonOrigin(path1, orientedPath2)

	""" compute spline for each set of points """
	path1Spline = SplineFit(path1, smooth=0.1)			  
	path2Spline = SplineFit(orientedPath2, smooth=0.1)

	
	""" get uniform selection of points along curve segments with the tangent angle interpolated """
	pathPoints1 = path1Spline.getUniformSamples(interpAngle=True)
	pathPoints2 = path2Spline.getUniformSamples(interpAngle=True)


	""" compute the angle derivative at each point of the curves """
	angDerivs1 = getAngleDerivatives(pathPoints1)
	angDerivs2 = getAngleDerivatives(pathPoints2)


	""" for each point on the child path, find its closest pair on the parent path """
	distances, indices = getClosestPairs(pathPoints2, pathPoints1)

	pathPoints2_rvrs = deepcopy(pathPoints2)
	pathPoints2_rvrs.reverse()

	""" get the closest point to the current junction pose and each points' distance to it """
	juncDists, minDist2, juncI = getPointDistances(pathPoints2, globalJunctionPose)
	print "minDist2, juncI:", minDist2, juncI
	juncI_rvrs = len(pathPoints2) - juncI - 1
	#juncDists_rvrs, minDist2, juncI_rvrs = getPointDistances(pathPoints2_rvrs, globalJunctionPose)



	""" find the first point on the child shoot that is within radius 1.0 of junction point """
	frontInd = 0
	backInd = len(pathPoints2)-1

	path2Indices = range(0,len(pathPoints2))
		
	frontFound = False
	backFound = False
	for k in path2Indices:
		if not frontFound and juncDists[k] <= 1.0:
			frontInd = k
			frontFound = True

	path2Indices.reverse()
	for k in path2Indices:
		if not backFound and juncDists[k] <= 1.0:
			backInd = k
			backFound = True

	print "frontFound, backFound:", frontFound, backFound
	print "frontInd, backInd:", frontInd, backInd


	""" compute the divergence point of the child shoot from the parent in the front """
	frontDepI, newFrontDepI, forePathIndex, maxFront, frontDiffAngles = computeDivergencePoint(juncI, frontInd, distances, indices, pathPoints1, pathPoints2)


	""" reverse the direction of all of the input datas and compute the divergence point from the rear """
	distances_rvrs = copy(distances)
	distances_rvrs.reverse()
	indices_rvrs = copy(indices)
	indices_rvrs.reverse()

	backDepI, newBackDepI, backPathIndex, maxBack, backDiffAngles = computeDivergencePoint(juncI_rvrs, backInd, distances_rvrs, indices_rvrs, pathPoints1, pathPoints2_rvrs)

	""" unreverse output data """
	backDepI = len(pathPoints2) - backDepI - 1
	newBackDepI = len(pathPoints2) - newBackDepI - 1
	backDiffAngles.reverse()




	""" reset to the tip match distance """
	maxFront = distances[frontInd]
	maxBack = distances[backInd]
	
	
	""" sum of closest points on front and back """
	""" select the one with minimal cost """		

	""" path section for our front departure hypothesis """
	pathSec1 = pathPoints2[:frontDepI+1]
	pathSec1.reverse()

	""" path section for our back departure hypothesis """
	pathSec2 = pathPoints2[backDepI:]
	
	""" distance of departure point from known junction point """
	p0 = pathSec1[0]
	juncDist1 = sqrt((globalJunctionPose[0]-p0[0])**2 + (globalJunctionPose[1]-p0[1])**2)

	p0 = pathSec2[0]
	juncDist2 = sqrt((globalJunctionPose[0]-p0[0])**2 + (globalJunctionPose[1]-p0[1])**2)

	print "pathSec1 hypothesis discrepancy distance:", juncDist1
	print "pathSec2 hypothesis discrepancy distance:", juncDist2


	secP1 = []
	secP2 = []

	"""
	cases
	1) near-zero length of path section, never departs, not the path, high juncDist
	2) near-zero length, no match count, low juncDist, correct path
	3) high path length, high match count, low juncDist, incorrect path 
	4) no match count, high junc dist, incorrect
	5) no match count, low junc dist, correct
	"""

	""" FIXME:  Assumes that we are only diverging in one direction, selects only one point """
	""" retrieve terminator points of the diverging section of the child shoot """
	if juncDist1 <= juncDist2:
		""" FIXME: if [0,frontDepI] has high overlap compared to total length, then we reject it """
		secP1 = pathPoints2[0]
		secP2 = pathPoints2[frontDepI]
		
		newPath2 = pathPoints2[:newFrontDepI]
		newPath2.reverse()

	else:
		""" FIXME: if [backDepI,len(distances)-1] has high overlap compared to total length, then we reject it """
		secP1 = pathPoints2[backDepI]
		secP2 = pathPoints2[len(distances)-1]

		newPath2 = pathPoints2[newBackDepI:]


	""" convert path so that the points are uniformly distributed """
	""" make sure path is greater than 5 points """
	newPath3 = ensureEnoughPoints(newPath2, max_spacing = 0.08, minPoints = 5)
	


	leafPath = deepcopy(newPath3)

	
	""" compute direction from which the curve diverges from parent shoot """
	frontVec = [0.,0.]
	indic = range(3)
	indic.reverse()
	
	for i in indic:
		if i+2 < len(leafPath):
			p1 = leafPath[i+2]
			p2 = leafPath[i]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			frontVec[0] += vec[0]
			frontVec[1] += vec[1]

	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	frontVec[0] /= frontMag
	frontVec[1] /= frontMag


	""" extend edge back towards the parent shoot such that it most likely overlaps """
	newP1 = (leafPath[0][0] + frontVec[0]*2, leafPath[0][1] + frontVec[1]*2)

	leafPath.insert(0,newP1)
	
	medial2 = deepcopy(leafPath)

	""" take the short length segment that we hope overlaps the parent shoot """
	edge1 = medial2[0:2]
	
	""" make a smaller version of these edges, supposedly so that we are not
		intersecting the parent in more than one place
	"""
	#newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)

	#edge1 = [newP1, edge1[1]]

	""" find the intersection points with the parent shoot """
	interPoints = []
	for k in range(len(path1)-1):
		shootEdge = [path1[k],path1[k+1]]
		isIntersect1, point1 = Intersect(edge1, shootEdge)
		if isIntersect1:
			interPoints.append(point1)
			break

	
	""" replace the extended edges with a termination point on the parent shoot """
	medial2 = medial2[1:]
	
	""" take the closest intersection point """
	if isIntersect1:
		termPoint = newPath3[0]
		minDist = 1e100
		minIndex = 0

		for k in range(len(interPoints)):
			intPoint = interPoints[k]
			dist = sqrt((intPoint[0]-termPoint[0])**2 + (intPoint[1]-termPoint[1])**2 )

			if dist < minDist:
				minDist = dist
				minIndex = k

		intPoint = interPoints[minIndex]	
		medial2.insert(0, intPoint)


	""" angle inverted, starting from parent shoot, out in direction of divergence """
	juncAng = acos(-frontVec[0])
	if -frontVec[1] < 0.0:
		juncAng = -juncAng
	
	
	""" get divergence point by using tangent intersection method """
	""" This method instead of angle/dist threshold approaches """

	try:

		foreIntI, backIntI, juncForeAng, juncBackAng = getTangentIntersections(pathPoints1, pathPoints2, frontDepI, backDepI, indices[frontDepI], indices[backDepI], juncI, indices[juncI], pathPlotCount, hypothesisID = hypothesisID, nodeID = nodeID, plotIter = True)
		print "foreIntI, backIntI:", foreIntI, backIntI

		""" junction distances are equal if both of the indices selected juncI as the departure point on path2'
			occurs if path2 does not come close enough to path1 to get under distance 0.1
		"""
		foreGlobJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], globalJunctionPose[2]]
		foreControlPoint = pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1]
		foreAngDeriv = angDerivs1[forePathIndex]

		backGlobJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], globalJunctionPose[2]]
		backControlPoint = pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1]
		backAngDeriv = angDerivs1[backPathIndex]

		print "juncDist1, juncDist2 =", juncDist1, juncDist2
		if juncDist1 == juncDist2:

			juncDist2 = sqrt((globalJunctionPose[0]-p0[0])**2 + (globalJunctionPose[1]-p0[1])**2)

			foreDist = sqrt((pathPoints1[foreIntI][0]-globalJunctionPose[0])**2 +  (pathPoints1[foreIntI][1]-globalJunctionPose[1])**2)
			backDist = sqrt((pathPoints1[backIntI][0]-globalJunctionPose[0])**2 +  (pathPoints1[backIntI][1]-globalJunctionPose[1])**2)

			print "foreDist, backDist =", foreDist, backDist

			if foreDist < backDist:
				globJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], juncForeAng]
				controlPoint = pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1]
				angDeriv = angDerivs1[forePathIndex]

			else:
				globJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], juncBackAng]
				controlPoint = pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1]
				angDeriv = angDerivs1[backPathIndex]
	
		elif juncDist1 < juncDist2:
			globJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], juncForeAng]
			controlPoint = pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1]
			angDeriv = angDerivs1[forePathIndex]
		else:
			globJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], juncBackAng]
			controlPoint = pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1]
			angDeriv = angDerivs1[backPathIndex]

		print "globJuncPose =", globJuncPose

	except:
		print "getTangentIntersections() failed!"
		foreDist = sqrt((pathPoints1[forePathIndex][0]-globalJunctionPose[0])**2 +  (pathPoints1[forePathIndex][1]-globalJunctionPose[1])**2)
		backDist = sqrt((pathPoints1[backPathIndex][0]-globalJunctionPose[0])**2 +  (pathPoints1[backPathIndex][1]-globalJunctionPose[1])**2)

		foreGlobJuncPose = [pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1], globalJunctionPose[2]]
		foreControlPoint = pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1]
		foreAngDeriv = angDerivs1[forePathIndex]

		backGlobJuncPose = [pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1], globalJunctionPose[2]]
		backControlPoint = pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1]
		backAngDeriv = angDerivs1[backPathIndex]

		if foreDist < backDist:
			globJuncPose = [pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1], globalJunctionPose[2]]
			controlPoint = pathPoints1[forePathIndex][0], pathPoints1[forePathIndex][1]
			angDeriv = angDerivs1[forePathIndex]
		else:
			globJuncPose = [pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1], globalJunctionPose[2]]
			controlPoint = pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1]
			angDeriv = angDerivs1[backPathIndex]

		#globJuncPose = [medial2[0][0], medial2[0][1], juncAng]
		#controlPoint = pathPoints1[backPathIndex][0], pathPoints1[backPathIndex][1]
	
	""" last junction point is the intersection point """

	if plotIter:

		pylab.clf()
		xP = []
		yP = []
		for p in path2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.5,0.5,1.0))

		if False:	
			P1 = pathPoints2[frontInd]
			P2 = pathPoints2[frontDepI]

			
			P3 = pathPoints2[backDepI]
			P4 = pathPoints2[backInd]

			" 0, frontDepI "
			xP = [P1[0], P2[0]]
			yP = [P1[1], P2[1]]
			pylab.scatter(xP,yP, color='b')		   

			" backDepI, -1 "
			xP = [P3[0], P4[0]]
			yP = [P3[1], P4[1]]
			pylab.scatter(xP,yP, color='g')		   

			" angle references "
			PF = pathPoints1[forePathIndex]
			xP = [PF[0],]
			yP = [PF[1],]
			pylab.scatter(xP,yP, color='b', alpha=0.5, zorder=100)		  

			PB = pathPoints1[backPathIndex]
			xP = [PB[0],]
			yP = [PB[1],]
			pylab.scatter(xP,yP, color='g', alpha=0.5, zorder=100)		  


			
		pylab.scatter([globalJunctionPose[0]],[globalJunctionPose[1]], color='r')		   
		pylab.scatter([globJuncPose[0]],[globJuncPose[1]], color='k')		 
		pylab.scatter([controlPoint[0]],[controlPoint[1]], color='k', alpha = 0.5)		   

		pylab.scatter([commonP1[0],commonP2[0]], [commonP1[1],commonP2[1]], color='b')

		xP = []
		yP = []
		for p in path1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(1.0,0.5,0.5))


		xP = []
		yP = []
		cutMedial = medial2[1:]
		for p in cutMedial:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k')

		print "trimDeparture:", pathPlotCount
		printStack()				 

		pylab.title("%1.2f %1.2f %1.2f %1.2f" % (juncDist1, juncDist2, juncAng, angDeriv))
		#pylab.savefig("trimDeparture2_%04u_%04u_%04u.png" % (hypothesisID, nodeID, pathPlotCount))
		pylab.savefig("trimDeparture2_%04u_%04u_%1.1f_%04u.png" % (hypothesisID, nodeID, arcDist, pathPlotCount))
		
		pathPlotCount += 1


	
	return globJuncPose, controlPoint, angDeriv

@logFunction
def computeBranch(pathID, parentID, origGlobJuncPose, origControlPose, childPath, parentPath, trimmedParent, localPathSegs, arcDist, localSkeletons, controlPoses, junctionPoses, parentPathIDs, numNodes=0, hypothesisID=0):

	if parentID != None: 

		""" parent path information """
		pathSpline = SplineFit(parentPath)

		""" initial position of junction """
		origJuncPose = copy(origGlobJuncPose)
		origJuncPose[2] = 0.0

		""" f(arcDist) = position on the spline curve of the parent path """
		""" NOTE:  angle is the tangent angle, not branch angle """
		""" LOCATION OF MODIFIED CONTROL POINT ON PARENT """
		newArcDist = arcDist
		newJuncPose = pathSpline.getPointOfDist(newArcDist)
		modControlPose = copy(newJuncPose)
		modControlPose[2] = 0.0

		"""
		1) point that has the junction
		2) direction of junction
		3) set of long paths that include junction path  (only 2)

		"""

		""" convert smooth segments from global to local coordinates """
		origJuncOrigin = Pose(origJuncPose)

		origControlOrigin = Pose(origControlPose)
		localJuncPose = origControlOrigin.convertGlobalPoseToLocal(origGlobJuncPose)
		offsetOrigin1 = Pose(modControlPose)
		modJuncPose = offsetOrigin1.convertLocalOffsetToGlobal(localJuncPose)


		""" convert junction point from canonical to local to modified global coordinates """

		#junctionDetails = self.leaf2LeafPathJunctions[pathID]
		#localPathSegs = junctionDetails["localSegments"]

		#localPathSegs = []

		#for k in range(len(smoothPathSegs)):
		#	pathSeg = smoothPathSegs[k]
		#	localSeg = []
		#	for p in pathSeg:
		#		p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
		#		localSeg.append(p1)
		#	
		#	localPathSegs.append(localSeg)


		""" place the segments back into global coordinates with modified junction pose instead """
		offsetOrigin1 = Pose(modControlPose)
		placedPathSegs = []
		for k in range(len(localPathSegs)):
			localSeg = localPathSegs[k]
			placedSeg = []
			for p in localSeg:
				p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
				placedSeg.append(p1)
			placedPathSegs.append(placedSeg)


		""" curves of parent shoot """
		globalPath2 = parentPath

		""" point on parent shoot curve closest to junction pose """
		#globalSpline2 = SplineFit(globalPath2, smooth=0.1)
		globalSpline2 = pathSpline
		globalSamples2 = globalSpline2.getUniformSamples(spacing = 0.04)
		originU2 = globalSpline2.findU(modControlPose)
		minDist, originU2, uPoint = globalSpline2.findClosestPoint(modControlPose)
		
		""" u1 value is meaningless since data is not a spline curve """
		""" angGuess is offset by tangent angle on curve, resulting in original branch angle """
		u2 = originU2
		u1 = 0.5
		angGuess = -uPoint[2] + modControlPose[2]
		

		""" use ICP to align the branch segments to the parent shoot """
		initGuess = [u1,u2,angGuess]
		resultPose1, lastCost1, matchCount1 = gen_icp.branchEstimateCost(initGuess, modControlPose, placedPathSegs, globalPath2, plotIter = True, n1 = pathID, n2 = parentID, arcDist = arcDist)


		""" perform computations to get new aligned branch pose from ICP result """
		poseOrigin = Pose(resultPose1)

		inversePose = poseOrigin.doInverse(resultPose1)
		poseOrigin = Pose(inversePose)

		finalPose = poseOrigin.convertLocalOffsetToGlobal(modControlPose)
		#poseOrigin2 = Pose(finalPose)

		origJuncPose = copy(origGlobJuncPose)
		origJuncPose[2] = 0.0
		#distDisc = sqrt((finalPose[0]-origJuncPose[0])**2 + (finalPose[1]-origJuncPose[1])**2)
		distDisc = 0.0

		""" FIXME:  what does this angle difference represent?  is finalPose a unique value? """
		angDisc = fabs(normalizeAngle(finalPose[2]-origJuncPose[2]))


		""" get the trimmed child shoot at the new designated branch point from parent """
		newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths =  trimBranch(pathID, parentID, origControlPose, modControlPose, origGlobJuncPose, modJuncPose, childPath, parentPath, trimmedParent, localPathSegs, plotIter=True, arcDist = arcDist, nodeID=numNodes, hypothesisID = hypothesisID)

		junctionPoses[pathID] = newGlobJuncPose
		controlPoses[pathID] = modControlPose

		spliceSkeleton = spliceSkeletons(localSkeletons, controlPoses, junctionPoses, parentPathIDs)

		""" the terminals of the parent shoot """
		parentTerm1 = parentPath[0]
		parentTerm2 = parentPath[-1]

		""" the terminals of the child shoot """
		dist1 = sqrt((newGlobJuncPose[0]-newPath3[0][0])**2 + (newGlobJuncPose[1]-newPath3[0][1])**2)
		dist2 = sqrt((newGlobJuncPose[0]-newPath3[-1][0])**2 + (newGlobJuncPose[1]-newPath3[-1][1])**2)

		if dist1 < dist2:
			childTerm = newPath3[-1]
		else:
			childTerm = newPath3[0]


		""" two different splices between child and parent terminals """
		termPaths = [(parentTerm1, childTerm), (parentTerm2, childTerm)]

		splicedPaths = []
		for termPath in termPaths:

			""" use splice skeleton to find shortest path """
			startPose = termPath[0]
			endPose = termPath[-1]

			minStartDist = 1e100
			minStartNode = None
			minEndDist = 1e100
			minEndNode = None
			
			for edge in spliceSkeleton.edges():
			
				globalNodePoint1 = edge[0]
				globalNodePoint2 = edge[1]

				dist1 = sqrt((globalNodePoint1[0]-startPose[0])**2 + (globalNodePoint1[1]-startPose[1])**2)
				dist2 = sqrt((globalNodePoint2[0]-startPose[0])**2 + (globalNodePoint2[1]-startPose[1])**2)

				if dist1 < minStartDist:
					minStartDist = dist1
					minStartNode = globalNodePoint1

				if dist2 < minStartDist:
					minStartDist = dist2
					minStartNode = globalNodePoint2

				dist1 = sqrt((globalNodePoint1[0]-endPose[0])**2 + (globalNodePoint1[1]-endPose[1])**2)
				dist2 = sqrt((globalNodePoint2[0]-endPose[0])**2 + (globalNodePoint2[1]-endPose[1])**2)

				if dist1 < minEndDist:
					minEndDist = dist1
					minEndNode = globalNodePoint1

				if dist2 < minEndDist:
					minEndDist = dist2
					minEndNode = globalNodePoint2


			startNode = minStartNode
			endNode = minEndNode


			shortestSpliceTree, shortestSpliceDist = spliceSkeleton.shortest_path(endNode)
			currNode = shortestSpliceTree[startNode]					 
			splicedSkel = []
			while currNode != endNode:
				splicedSkel.append(currNode)
				nextNode = shortestSpliceTree[currNode]
				currNode = nextNode
			splicedSkel.append(currNode)

			splicedSkel = ensureEnoughPoints(splicedSkel, max_spacing = 0.08, minPoints = 5)
			spliceSpline1 = SplineFit(splicedSkel, smooth=0.1)
			splicePoints1 = spliceSpline1.getUniformSamples()

			splicedPaths.append(splicePoints1)
		

		""" cache the generated result and package it for later retrieval """
		newSplices = deepcopy(splicedPaths)
		# initProb = 0.0
		#initProb = matchCount1*(1000.0-lastCost1)
		initProb = matchCount1
		#part2 = (parentID, modJuncPose, modControlPose, newArcDist, pathID, matchCount1, lastCost1, distDisc, angDisc, initProb, newSplices, juncDiscAngle, juncDiscDist)

		part2 = {}
		part2["parentID"] = parentID
		part2["modJuncPose"] = newGlobJuncPose
		part2["modControlPose"] = modControlPose
		part2["newArcDist"] = newArcDist
		part2["pathID"] = pathID
		part2["matchCount"] = matchCount1
		part2["lastCost"] = lastCost1
		part2["distDisc"] = distDisc
		part2["angDisc"] = angDisc
		part2["initProb"] = initProb
		part2["newSplices"] = newSplices
		part2["juncDiscAngle"] = juncDiscAngle
		part2["juncDiscDist"] = juncDiscDist


		"""
		result values are as follows:
		0 = parent ID
		1 = pose of the branch point
		2 = arc distance on the parent curve of the branch point
		3 = child ID
		4 = match count from ICP algorithm
		5 = last cost from ICP algorithm
		6 = discrepancy distance of new pose from old pose, zeroed for this function
		7 = angle discrepancy difference of new pose from old pose
		8 = initial probability value of this branch position, initialized to 0.0 in this function
		9 = set of splices including the branch at this position
		10 = angle discrepancy 
		11 = dist discrepancy 
		"""

		return part2
 
@logFunction
def trimBranch2(pathID, parentPathID, origControlPose, modControlPose, origGlobJuncPose, modJuncPose, childPath, parentPath, trimmedParent, localPathSegs, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):

	global pathPlotCount 



	globJuncPose = origGlobJuncPose

	origJuncPose = copy(origGlobJuncPose)
	origJuncPose[2] = 0.0
	origControlOrigin = Pose(origControlPose)


	offsetOrigin1 = Pose(modControlPose)
	
	partPathSegs = []

	path1 = parentPath
	path2 = childPath

	localPath2 = []
	for p in path2:
		p1 = origControlOrigin.convertGlobalToLocal(p)
		localPath2.append(p1)

	particlePath2 = []
	for p in localPath2:
		p1 = offsetOrigin1.convertLocalToGlobal(p)
		particlePath2.append(p1)

	""" get departing sections of overlapped curves """
	""" branch point as input """
	secP1, secP2 = getOverlapDeparture(modJuncPose, parentPathID, pathID, path1, particlePath2, plotIter = False)				 

	minI_1 = 0		  
	minI_2 = 0
	minDist_1 = 1e100
	minDist_2 = 1e100		 
	for i in range(len(particlePath2)):
		pnt = particlePath2[i]
		dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
		dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
	
		if dist1 < minDist_1:
			minDist_1 = dist1
			minI_1 = i

		if dist2 < minDist_2:
			minDist_2 = dist2
			minI_2 = i

	term0 = particlePath2[0]
	termN = particlePath2[-1]

	""" smallest distance is the terminal point """
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
	
	if minDist_1 < minDist_2:
		junctionPoint_K = secP1
		juncI_K = minI_1
	else:
		junctionPoint_K = secP2
		juncI_K = minI_2
	
	minDist2 = 1e100
	juncI = 0		 
	for i in range(len(particlePath2)):
		pnt = particlePath2[i]
		dist = sqrt((pnt[0]-modJuncPose[0])**2 + (pnt[1]-modJuncPose[1])**2)
	
		if dist < minDist2:
			minDist2 = dist
			juncI = i
	 
	 
	
	print "len(path1):", len(path1)
	print "len(particlePath2):", len(particlePath2)
	print "juncI:", juncI
	print "minDist:", minDist_1, minDist_2
	
	""" now we have closest point to departure point. """
	""" Which side is the departing side? """	 

			
	if minI == 0:
		"""secP1 is terminal 0"""
		index = juncI-10
		if index < 1:
			index = 1

		
		newPath2 = particlePath2[:index+1]
		newPath2.reverse()
	
	elif minI == 1:
		"""secP1 is terminal N"""
		index = juncI+10
		if index >= len(particlePath2)-1:
			
			""" ensure at least 2 elements in path """
			index = len(particlePath2)-2

		newPath2 = particlePath2[index:]

		
	elif minI == 2:
		"""secP2 is terminal 0"""
		index = juncI-10
		if index < 1:
			index = 1

		newPath2 = particlePath2[:index+1]
		newPath2.reverse()
		
	elif minI == 3:
		"""secP2 is terminal N"""
		index = juncI+10
		if index >= len(particlePath2)-1:
			""" ensure at least 2 elements in path """
			index = len(particlePath2)-2
		
		newPath2 = particlePath2[index:]
	
	else:
		print """no terminal found"""
		raise
					
	""" convert path so that the points are uniformly distributed """
	newPath3 = ensureEnoughPoints(newPath2, max_spacing = 0.08, minPoints = 5)


	newGlobJuncPose1, controlPoint1, angDeriv1 = getBranchPoint(modJuncPose, parentPathID, pathID, trimmedParent, particlePath2, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)
	newGlobJuncPose2, controlPoint2, angDeriv2 = getBranchPoint(modJuncPose, pathID, parentPathID, particlePath2, trimmedParent, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)


	juncDiscDist1 = sqrt((newGlobJuncPose1[0]-modJuncPose[0])**2 + (newGlobJuncPose1[1]-modJuncPose[1])**2)
	juncDiscDist2 = sqrt((newGlobJuncPose2[0]-modJuncPose[0])**2 + (newGlobJuncPose2[1]-modJuncPose[1])**2)

	if angDeriv1 < angDeriv2:
		newGlobJuncPose = newGlobJuncPose1
		controlPoint = controlPoint1
		angDeriv = angDeriv1
	else:

		""" only diverge from the parent if the determined branch point is in the neighborhood of the original """
		if juncDiscDist2 < 1.0:
			newGlobJuncPose = newGlobJuncPose2
			controlPoint = controlPoint2
			angDeriv = angDeriv2
		else:
			newGlobJuncPose = newGlobJuncPose1
			controlPoint = controlPoint1
			angDeriv = angDeriv1

	""" junction discrepency distance """
	juncDiscDist = 0.0
	juncDiscAngle = 0.0



	""" for each path, attempt to join with its parent path """
	junctionPoint = [newGlobJuncPose[0],newGlobJuncPose[1]]






	localSkeletons = {}
	controlPoses = {}
	junctionPoses = {}
	parentPathIDs = {}
	for pathID in pathIDs:
		localSkeletons[pathID] = self.leaf2LeafPathJunctions[pathID]["skeletonGraph"]
		controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
		junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

		cPath = self.getPath(pathID)
		parentPathID = cPath["parentID"]
		parentPathIDs[pathID] = parentPathID

	#def trimBranch2(pathID, parentPathID, origControlPose, modControlPose, origGlobJuncPose, modJuncPose, childPath, parentPath, trimmedParent, localPathSegs, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):
	spliceSkeleton = spliceSkeletons(localSkeletons, controlPoses, junctionPoses, parentPathIDs)

	""" the terminals of the parent shoot """
	parentTerm1 = trimmedParent[0]
	parentTerm2 = trimmedParent[-1]

	""" the terminals of the child shoot """
	dist1 = sqrt((newGlobJuncPose[0]-newPath3[0][0])**2 + (newGlobJuncPose[1]-newPath3[0][1])**2)
	dist2 = sqrt((newGlobJuncPose[0]-newPath3[-1][0])**2 + (newGlobJuncPose[1]-newPath3[-1][1])**2)

	if dist1 < dist2:
		childTerm = newPath3[-1]
	else:
		childTerm = newPath3[0]


	""" two different splices between child and parent terminals """
	termPaths = [(parentTerm1, childTerm), (parentTerm2, childTerm)]

	splicedPaths = []
	for termPath in termPaths:

		""" use splice skeleton to find shortest path """
		startPose = termPath[0]
		endPose = termPath[-1]

		minStartDist = 1e100
		minStartNode = None
		minEndDist = 1e100
		minEndNode = None
		
		for edge in spliceSkeleton.edges():
		
			globalNodePoint1 = edge[0]
			globalNodePoint2 = edge[1]

			dist1 = sqrt((globalNodePoint1[0]-startPose[0])**2 + (globalNodePoint1[1]-startPose[1])**2)
			dist2 = sqrt((globalNodePoint2[0]-startPose[0])**2 + (globalNodePoint2[1]-startPose[1])**2)

			if dist1 < minStartDist:
				minStartDist = dist1
				minStartNode = globalNodePoint1

			if dist2 < minStartDist:
				minStartDist = dist2
				minStartNode = globalNodePoint2

			dist1 = sqrt((globalNodePoint1[0]-endPose[0])**2 + (globalNodePoint1[1]-endPose[1])**2)
			dist2 = sqrt((globalNodePoint2[0]-endPose[0])**2 + (globalNodePoint2[1]-endPose[1])**2)

			if dist1 < minEndDist:
				minEndDist = dist1
				minEndNode = globalNodePoint1

			if dist2 < minEndDist:
				minEndDist = dist2
				minEndNode = globalNodePoint2


		startNode = minStartNode
		endNode = minEndNode


		shortestSpliceTree, shortestSpliceDist = spliceSkeleton.shortest_path(endNode)
		currNode = shortestSpliceTree[startNode]					 
		splicedSkel = []
		while currNode != endNode:
			splicedSkel.append(currNode)
			nextNode = shortestSpliceTree[currNode]
			currNode = nextNode
		splicedSkel.append(currNode)

		spliceSpline1 = SplineFit(splicedSkel, smooth=0.1)
		splicePoints1 = spliceSpline1.getUniformSamples()

		#splicedPaths.append(splicedSkel)
		splicedPaths.append(splicedPoints1)
		

	if plotIter:
		#if False:
			
		#hypothesisID = 0
		numNodes = 0
		pathPlotCount = 0

		
		pylab.clf()
		xP = []
		yP = []
		for p in newPath3:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.5,0.5,1.0))

		xP = []
		yP = []
		for p in path1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(1.0,0.5,0.5))

		
		pylab.scatter([newGlobJuncPose[0],], [newGlobJuncPose[1],], color='r')
		pylab.scatter([modControlPose[0],], [modControlPose[1],], color='b')
		pylab.scatter([origControlPose[0],], [origControlPose[1],], color='g')
		pylab.scatter([modJuncPose[0],], [modJuncPose[1],], color='m')

		xP = []
		yP = []
		for p in splicePoints1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k', alpha=0.3)

		xP = []
		yP = []
		for p in splicePoints2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k', alpha=0.3)


		pylab.axis("equal")
		pylab.title("hyp %d nodeID %d %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f" % ( hypothesisID, nodeID, juncDiscDist, juncDiscAngle, origGlobJuncPose[2], newGlobJuncPose[2], angDeriv1, angDeriv2))
		pylab.savefig("trimDeparture_%04u_%04u_%1.1f.png" % (hypothesisID, nodeID, arcDist))

		print "saving trimDeparture_%04u_%04u_%1.1f.png" % (hypothesisID, nodeID, arcDist)
		
		pathPlotCount += 1


	return newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths


@logFunction
def trimBranch(pathID, parentPathID, origControlPose, modControlPose, origGlobJuncPose, modJuncPose, childPath, parentPath, trimmedParent, localPathSegs, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):

	global pathPlotCount 



	globJuncPose = origGlobJuncPose

	origJuncPose = copy(origGlobJuncPose)
	origJuncPose[2] = 0.0
	origControlOrigin = Pose(origControlPose)

	#localPathSegs = []

	#for k in range(len(smoothPathSegs)):
	#	pathSeg = smoothPathSegs[k]
	#	localSeg = []
	#	for p in pathSeg:
	#		p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
	#		localSeg.append(p1)

	#	localPathSegs.append(localSeg)


	offsetOrigin1 = Pose(modControlPose)
	
	partPathSegs = []

	#for k in range(len(localPathSegs)):
	#	localSeg = localPathSegs[k]
	#	for p in localSeg:
	#		p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)

	path1 = parentPath
	path2 = childPath

	localPath2 = []
	for p in path2:
		p1 = origControlOrigin.convertGlobalToLocal(p)
		localPath2.append(p1)

	particlePath2 = []
	for p in localPath2:
		p1 = offsetOrigin1.convertLocalToGlobal(p)
		particlePath2.append(p1)

	""" get departing sections of overlapped curves """
	""" branch point as input """
	secP1, secP2 = getOverlapDeparture(modJuncPose, parentPathID, pathID, path1, particlePath2, plotIter = False)				 

	minI_1 = 0		  
	minI_2 = 0
	minDist_1 = 1e100
	minDist_2 = 1e100		 
	for i in range(len(particlePath2)):
		pnt = particlePath2[i]
		dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
		dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
	
		if dist1 < minDist_1:
			minDist_1 = dist1
			minI_1 = i

		if dist2 < minDist_2:
			minDist_2 = dist2
			minI_2 = i

	term0 = particlePath2[0]
	termN = particlePath2[-1]

	""" smallest distance is the terminal point """
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
	
	if minDist_1 < minDist_2:
		junctionPoint_K = secP1
		juncI_K = minI_1
	else:
		junctionPoint_K = secP2
		juncI_K = minI_2
	
	minDist2 = 1e100
	juncI = 0		 
	for i in range(len(particlePath2)):
		pnt = particlePath2[i]
		dist = sqrt((pnt[0]-modJuncPose[0])**2 + (pnt[1]-modJuncPose[1])**2)
	
		if dist < minDist2:
			minDist2 = dist
			juncI = i
	 
	 
	
	print "len(path1):", len(path1)
	print "len(particlePath2):", len(particlePath2)
	print "juncI:", juncI
	print "minDist:", minDist_1, minDist_2
	
	""" now we have closest point to departure point. """
	""" Which side is the departing side? """	 

			
	if minI == 0:
		"""secP1 is terminal 0"""
		index = juncI-10
		if index < 1:
			index = 1

		
		newPath2 = particlePath2[:index+1]
		newPath2.reverse()
	
	elif minI == 1:
		"""secP1 is terminal N"""
		index = juncI+10
		if index >= len(particlePath2)-1:
			
			""" ensure at least 2 elements in path """
			index = len(particlePath2)-2

		newPath2 = particlePath2[index:]

		
	elif minI == 2:
		"""secP2 is terminal 0"""
		index = juncI-10
		if index < 1:
			index = 1

		newPath2 = particlePath2[:index+1]
		newPath2.reverse()
		
	elif minI == 3:
		"""secP2 is terminal N"""
		index = juncI+10
		if index >= len(particlePath2)-1:
			""" ensure at least 2 elements in path """
			index = len(particlePath2)-2
		
		newPath2 = particlePath2[index:]
	
	else:
		print """no terminal found"""
		raise
					
	""" convert path so that the points are uniformly distributed """
	newPath3 = ensureEnoughPoints(newPath2, max_spacing = 0.08, minPoints = 5)


	newGlobJuncPose1, controlPoint1, angDeriv1 = getBranchPoint(modJuncPose, parentPathID, pathID, trimmedParent, particlePath2, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)
	newGlobJuncPose2, controlPoint2, angDeriv2 = getBranchPoint(modJuncPose, pathID, parentPathID, particlePath2, trimmedParent, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)
	#newGlobJuncPose1, controlPoint1, angDeriv1 = getBranchPoint(modControlPose, parentPathID, pathID, path1, particlePath2, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)
	#newGlobJuncPose2, controlPoint2, angDeriv2 = getBranchPoint(modControlPose, pathID, parentPathID, particlePath2, path1, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)

	juncDiscDist1 = sqrt((newGlobJuncPose1[0]-modJuncPose[0])**2 + (newGlobJuncPose1[1]-modJuncPose[1])**2)
	juncDiscDist2 = sqrt((newGlobJuncPose2[0]-modJuncPose[0])**2 + (newGlobJuncPose2[1]-modJuncPose[1])**2)

	if angDeriv1 < angDeriv2:
		newGlobJuncPose = newGlobJuncPose1
		controlPoint = controlPoint1
		angDeriv = angDeriv1
	else:

		""" only diverge from the parent if the determined branch point is in the neighborhood of the original """
		if juncDiscDist2 < 1.0:
			newGlobJuncPose = newGlobJuncPose2
			controlPoint = controlPoint2
			angDeriv = angDeriv2
		else:
			newGlobJuncPose = newGlobJuncPose1
			controlPoint = controlPoint1
			angDeriv = angDeriv1

	""" junction discrepency distance """
	#origControlPose = copy(origGlobJuncPose)
	juncDiscDist = 0.0
	juncDiscAngle = 0.0
	#origControlPose = copy(origGlobJuncPose)
	#juncDiscDist = sqrt((newGlobJuncPose[0]-modControlPose[0])**2 + (newGlobJuncPose[1]-modControlPose[1])**2)
	#juncDiscAngle = normalizeAngle(origControlPose[2]-newGlobJuncPose[2])



	""" for each path, attempt to join with its parent path """
	junctionPoint = [newGlobJuncPose[0],newGlobJuncPose[1]]

	path1 = newPath3
	path2 = trimmedParent
	
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
	
	joins = []
	joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])


	""" get junctions """ 
	
	junctions = {}
	junctions[pathID] = [0, newGlobJuncPose, (parentPathID,minI2), path2[minI2], minI1]


	""" create a tree with the node IDs and then stitch them together with joins """
	pathGraph = graph.graph()

	for k in range(len(path1)):
		pathGraph.add_node((pathID, k), path1[k])
	for k in range(len(path1)-1):
		pathGraph.add_edge((pathID, k), (pathID, k+1))
	
	for k in range(len(path2)):
		pathGraph.add_node((parentPathID, k), path2[k])
	for k in range(len(path2)-1):
		pathGraph.add_edge((parentPathID, k), (parentPathID, k+1))

	""" join with the junction in between the join points """
	for k in range(len(joins)):
		join = joins[k]
		pathGraph.add_edge(join[0], join[1])

	""" get terminals """
	topDict = {}
	terminals = {}


	""" parent is either the root path or some sub path """
	""" we only care about the two terminals of our immediate parent """


	terminals[parentPathID] = [(parentPathID, 0), path2[0]]
	terminals[parentPathID+1] = [(parentPathID, len(path2)-1), path2[len(path2)-1]]
	topDict["t%u" % 0] = terminals[parentPathID][0]
	topDict["t%u" % 1] = terminals[parentPathID+1][0]
		
	minI1 = junctions[pathID][4]
	
	""" get the two terminals of the current pathID path """
	""" determine which side is junction and which side is terminal """
	if minI1 > len(path1)-1 - minI1:    
		terminals[pathID+1] = [(pathID,0), path1[0]]
		topDict["t%u" % (pathID+1)] = (pathID,0)
	else:
		terminals[pathID+1] = [(pathID,len(path1)-1), path1[len(path1)-1]]
		topDict["t%u" % (pathID+1)] = (pathID,len(path1)-1)



	""" determine which paths are leaves """
	isAParent = {}
	parents = {}
	pathIDGraph = graph.graph()

	isAParent[pathID] = False
	pathIDGraph.add_node(pathID, [])

	isAParent[parentPathID] = True
	pathIDGraph.add_node(parentPathID, [])
	pathIDGraph.add_edge(pathID, parentPathID)
	
	leafCount = 1

	""" build topological graph """
	topGraph = graph.graph()
	
	""" add topological nodes to the graph """
	topNodes = []
	for k, term in terminals.iteritems():
		topGraph.add_node(term[0], term[1])
		topNodes.append(term[0])
	for k, junc in junctions.iteritems():
		topGraph.add_node(junc[2], junc[3])
		topNodes.append(junc[2])

	topGraph.add_edge(terminals[parentPathID][0], junctions[pathID][2])
	topGraph.add_edge(terminals[parentPathID+1][0], junctions[pathID][2])
	topGraph.add_edge(terminals[pathID+1][0], junctions[pathID][2])



	startKey1 = terminals[parentPathID][0]
	endKey1 = terminals[pathID+1][0]

	startKey2 = terminals[parentPathID+1][0]
	endKey2 = terminals[pathID+1][0]

	joinPairs = []
	shortestPathSpanTree, shortestDist = pathGraph.shortest_path(endKey1)
	currNode = shortestPathSpanTree[startKey1]					 
	splicedPath = []
	while currNode != endKey1:
		splicedPath.append(pathGraph.get_node_attributes(currNode))
		nextNode = shortestPathSpanTree[currNode]
		
		for join in joins:
			if join[0] == currNode and join[1] == nextNode:
				joinPairs.append((len(splicedPath)-1,len(splicedPath)))
			elif join[1] == currNode and join[0] == nextNode:
				joinPairs.append((len(splicedPath)-1,len(splicedPath)))
			
		currNode = nextNode
		
	splicedPath.append(pathGraph.get_node_attributes(currNode))

	""" if there are any joins, we must interpolate a smooth transition """
	lastIndex = 0
	newSplicePath1 = []
	for pair in joinPairs:
		index1 = pair[0]
		index2 = pair[1]
		cPoints = [splicedPath[index1-1], splicedPath[index1], splicedPath[index2], splicedPath[index2+1]]				 
		spline = SplineFit(cPoints, smooth=0.0, kp=3)
		points = spline.getUniformSamples(spacing = 0.1)
		
		newSplicePath1 += splicedPath[lastIndex:index1]
		newSplicePath1 += points
		
		lastIndex = index2+2

	newSplicePath1 += splicedPath[lastIndex:]




	joinPairs = []
	shortestPathSpanTree, shortestDist = pathGraph.shortest_path(endKey2)
	currNode = shortestPathSpanTree[startKey2]					 
	splicedPath = []
	while currNode != endKey2:
		splicedPath.append(pathGraph.get_node_attributes(currNode))
		nextNode = shortestPathSpanTree[currNode]
		
		for join in joins:
			if join[0] == currNode and join[1] == nextNode:
				joinPairs.append((len(splicedPath)-1,len(splicedPath)))
			elif join[1] == currNode and join[0] == nextNode:
				joinPairs.append((len(splicedPath)-1,len(splicedPath)))
			
		currNode = nextNode
		
	splicedPath.append(pathGraph.get_node_attributes(currNode))

	""" if there are any joins, we must interpolate a smooth transition """
	lastIndex = 0
	newSplicePath2 = []
	for pair in joinPairs:
		index1 = pair[0]
		index2 = pair[1]
		cPoints = [splicedPath[index1-1], splicedPath[index1], splicedPath[index2], splicedPath[index2+1]]				 
		spline = SplineFit(cPoints, smooth=0.0, kp=3)
		points = spline.getUniformSamples(spacing = 0.1)
		
		newSplicePath2 += splicedPath[lastIndex:index1]
		newSplicePath2 += points
		
		lastIndex = index2+2

	newSplicePath2 += splicedPath[lastIndex:]


	spliceSpline1 = SplineFit(newSplicePath1, smooth=0.1)
	spliceSpline2 = SplineFit(newSplicePath2, smooth=0.1)

	splicePoints1 = spliceSpline1.getUniformSamples()
	splicePoints2 = spliceSpline2.getUniformSamples()

	splicedPaths = []
	splicedPaths.append(newSplicePath1)
	splicedPaths.append(newSplicePath2)


	if plotIter:
		#if False:
			
		#hypothesisID = 0
		numNodes = 0
		pathPlotCount = 0

		
		pylab.clf()
		xP = []
		yP = []
		for p in newPath3:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.5,0.5,1.0))

		xP = []
		yP = []
		for p in path1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(1.0,0.5,0.5))

		
		pylab.scatter([newGlobJuncPose[0],], [newGlobJuncPose[1],], color='r')
		pylab.scatter([modControlPose[0],], [modControlPose[1],], color='b')
		pylab.scatter([origControlPose[0],], [origControlPose[1],], color='g')
		pylab.scatter([modJuncPose[0],], [modJuncPose[1],], color='m')

		xP = []
		yP = []
		for p in splicePoints1:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k', alpha=0.3)

		xP = []
		yP = []
		for p in splicePoints2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k', alpha=0.3)


		pylab.axis("equal")
		pylab.title("hyp %d nodeID %d %1.2f %1.2f %1.2f %1.2f %1.2f %1.2f" % ( hypothesisID, nodeID, juncDiscDist, juncDiscAngle, origGlobJuncPose[2], newGlobJuncPose[2], angDeriv1, angDeriv2))
		pylab.savefig("trimDeparture_%04u_%04u_%1.1f.png" % (hypothesisID, nodeID, arcDist))

		print "saving trimDeparture_%04u_%04u_%1.1f.png" % (hypothesisID, nodeID, arcDist)
		
		pathPlotCount += 1


	return newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths


class MapState:
	
	def __init__(self, poseData, hypothesisID):
		
		""" raw sensor data products """
		self.poseData = deepcopy(poseData)
		
		""" pose particles for different combinations of poses and branch points """
		#self.numPoseParticles = 100
		self.numPoseParticles = 40
		self.poseParticles = {}
		self.poseParticles["numParticles"] = self.numPoseParticles
		self.poseParticles["updateCount"] = 0
		self.poseParticles["snapshots2"] = {0 : ()}


		""" pose information for each of the spatial map nodes """
		self.isNodeBranching = {}
		self.nodeRawPoses = {}
		self.nodePoses = {}
		self.gndPoses = {}
		self.gndRawPoses = {}
		self.origPoses = {}

		""" evaluation result of map consistency """
		self.utility = 0.0

		""" if a branching decision was made, and this is the no-branch state """
		self.isNotBranched = False

		""" criteria for rejecting this map state """
		self.mapOverlapSum = 0.0
		self.isNoLocalize = False

		""" unique ID for this MapState instantation """
		self.hypothesisID = hypothesisID

		""" the alpha hull and maximum length path of shoot skeleton """
		self.paths = {0 : []}
		self.hulls = {0 : []}

		""" intermediate data structure for computing the shoot skeleton """
		self.leaf2LeafPathJunctions = {}

		""" shoot spine only, separated from the skeleton at the branch point """
		self.trimmedPaths  = {}

		""" membership, parent-relationship, branch point, and control point for a shoot """
		self.pathClasses = {}
		self.pathClasses[0] = {"parentID" : None,
					"branchNodeID" : None,
					"localJunctionPose" : None, 
					"sameProb" : {},
					"nodeSet" : [],
					"globalJunctionPose" : None,
					"controlPose" : [0.0,0.0,0.0] }		

		""" number of shoots in the current map state """
		self.pathIDs = 1

		""" travel destinations on the shoot map """
		self.pathTermsVisited = {0: False}

		""" number of branch hypotheses and spacing between them for a branch point distribution """
		self.DIV_LEN = 0.2
		self.NUM_BRANCHES = 5
		#self.DIV_LEN = 0.1
		#self.NUM_BRANCHES = 1


		""" the entrance point of the environment """
		self.rootPoint = [-3.2, 0.0]
		
		""" data structures for computing splices on the shoot map """
		self.pathGraph = graph.graph()
		self.joins = []
		self.junctions = {}
		self.terminals = {}
		self.allSplices = {}
		self.isChanged = True



		""" plotting colors """
		self.colors = []

		self.colors.append([127, 127, 255])
		self.colors.append([127, 255, 127])
		self.colors.append([255, 127, 127])
		self.colors.append([255, 127, 255])
		self.colors.append([255, 255, 127])
		self.colors.append([127, 255, 255])

		for color in self.colors:
			color[0] = float(color[0])/256.0
			color[1] = float(color[1])/256.0
			color[2] = float(color[2])/256.0

		for i in range(1000):
			self.colors.append((random.random(),random.random(),random.random()))


		""" plot counts """
		self.alphaPlotCount = 0
		self.medialCount = 0
		self.topCount = 0
		self.spliceCount = 0
		self.multiDepCount = 0
		self.overlapPlotCount = 0
		self.tempCount = 0
		self.pathPlotCount = 0
		self.pathPlotCount2 = 0
		self.overlapPlotCount2 = 0
		self.medialSoupCount = 0
		self.posePlotCount = 0



	def getNodePose(self, nodeID):
		return copy(self.nodePoses[nodeID])
	
	def setNodePose(self, nodeID, newPose):

		print "setNodePose(", nodeID, ",", newPose, ")"

		oldGPACPose = self.getNodePose(nodeID)

		self.nodePoses[nodeID] = newPose

		gpacProfile = Pose(oldGPACPose)
		
		localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[nodeID])
		
		" go back and convert this from GPAC pose to estPose "
		newProfile = Pose(newPose)
		newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
		
		self.nodeRawPoses[nodeID] = newEstPose

	@logFunction
	def initializePoseParticles(self):

	
		pose0 = self.getNodePose(0)
		pose1 = self.getNodePose(1)

		path0 = deepcopy(self.paths[0])

		termPoint = self.terminals[0][1]

		termDist0 = sqrt((path0[0][0]-termPoint[0])**2 + (path0[0][1]-termPoint[1])**2)
		termDistN = sqrt((path0[-1][0]-termPoint[0])**2 + (path0[-1][1]-termPoint[1])**2)

		if termDist0 > termDistN: 
			path0.reverse()

		pathSpline = SplineFit(path0, smooth=0.1)

		minDist0, newU0, newP0 = pathSpline.findClosestPoint(pose0)
		minDist1, newU1, newP1 = pathSpline.findClosestPoint(pose1)

		dist0 = pathSpline.dist_u(newU0)
		dist1 = pathSpline.dist_u(newU1)

		#newAng0 = pathSpline.angle_u(newU0)
		#newAng1 = pathSpline.angle_u(newU1)
		newAng0 = pose0[2]
		newAng1 = pose1[2]

		totalLen = pathSpline.length()

		pathID = 0
		avgDist = (dist0 + dist1)/2.0
		spliceID = ('t0', 't1')

		#pID, pathI = self.topDict["t%u" % (pathID+1)]

		part = (pathID, avgDist, ('t0','t1'))

		#self.poseParticles = {}
		#self.poseParticles["numParticles"] = 40
		#self.poseParticles["updateCount"] = 0
		#self.poseParticles["snapshots"] = {0 : ()}

		numParticles = self.poseParticles["numParticles"]
		#initDist = []
		initDist2 = []

		" create the initial particle distibution "
		while len(initDist2) < numParticles:
			hypDist = random.gauss(avgDist, 0.5)

			if hypDist > 0.0 and hypDist <= totalLen:
				hypPoint = pathSpline.getPointOfDist(hypDist)
				hypPose0 = copy(hypPoint)
				hypPose0[2] = newAng0
				hypPose1 = copy(hypPoint)
				hypPose1[2] = newAng1
				#initDist.append((hypPose, pathID, hypDist, ('t0', 't1'), 0.0))
				#initDist.append((hypPose0, hypPose1, pathID, hypDist, ('t0', 't1'), 0.0))

				particleObj = Particle(hypPose0, hypPose1, pathID, hypDist, 0.0, self.hypothesisID)
				particleObj.spliceCurve = deepcopy(self.paths[0])
				particleObj.addNode(0,0)
				particleObj.addNode(1,0)
				initDist2.append(particleObj)

		#print "setting pose particles", len(initDist), initDist
		self.poseParticles["snapshots2"][0] = initDist2

		#self.isChanged = True
		#self.junctions = {}
		#self.terminals = {}
		#self.allSplices = {}

	@logFunction
	def batchDisplaceParticles(self, nodeID0, nodeID1):

		batchJobs = []

		particleDist2 = self.poseParticles["snapshots2"][0]
		for particleIndex in range(len(particleDist2)):

			part = particleDist2[particleIndex]

			hypPose0 = part.pose0
			hypPose1 = part.pose1
			prevPose0 = part.prevPose0
			prevPose1 = part.prevPose1

			staticSplicedPaths0, spliceTerms0, staticSplicePathIDs0 = self.getSplicesByNearJunction(hypPose0)
			staticSplicedPaths1, spliceTerms1, staticSplicePathIDs1 = self.getSplicesByNearJunction(hypPose1)

			batchJobs.append([self.poseData, part, particleIndex, nodeID1, prevPose0, prevPose1, hypPose0, hypPose1, self.paths[0], staticSplicedPaths0, staticSplicedPaths1])

		results = batchDisplaceParticles(batchJobs)

		for result in results:
			particleIndex = result[0]
			part = particleDist2[particleIndex]
			part.pose0 = result[1]
			part.pose1 = result[2]

			part.displacePose(result[1], result[2])

	@logFunction
	def displacePoseParticles(self, nodeID0, nodeID1, travelDist0, travelDist1):

		#self.isChanged = True

		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]

		updateCount += 1
		self.poseParticles["updateCount"] = updateCount

		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]


		medialSpline = SplineFit(medial0, smooth=0.1)

		newPartDist2 = []
		for part in particleDist2:


			hypPose0 = part.pose0
			hypPose1 = part.pose1
			pathID = part.memberPaths[0]
			hypDist = part.hypDist
			spliceID = part.spliceName
			corresp = part.nodeCorrespondence

			currSplice0 = part.spliceCurve

			poseOrigin0 = Pose(hypPose0)
			globalMedial0 = []
			for p in medial0:
				globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

			orientedSplicePath = orientPath(currSplice0, globalMedial0)				
			pathSpline = SplineFit(orientedSplicePath, smooth=0.1)

			#pathSpline = SplineFit(currSplice0, smooth=0.1)


			minDist0, oldU0, oldP0 = pathSpline.findClosestPoint(hypPose0)
			minDist1, oldU1, oldP1 = pathSpline.findClosestPoint(hypPose1)

			thisDist0 = travelDist0
			thisDist1 = travelDist1

			moveChance = random.random()
			" 80% chance we move forward, 20% we move backward "    
			if moveChance >= 0.1:
				thisDist0 = travelDist0
				thisDist1 = travelDist1
			else:
				thisDist0 = -travelDist0
				thisDist1 = -travelDist1

			thisDist0 = random.gauss(thisDist0,0.6)
			thisDist1 = random.gauss(thisDist1,0.6)


			#newU0 = pathSpline.getUOfDist(oldU0, travelDist0)
			newU0 = pathSpline.getUOfDist(oldU0, thisDist0)
			newU1 = pathSpline.getUOfDist(oldU1, thisDist1)

			newDist0 = pathSpline.dist_u(newU0)
			newP0 = pathSpline.point_u(newU0)
			newAng0 = pathSpline.angle_u(newU0)
			oldAng0 = hypPose0[2]
			newPose0 = copy(newP0)
			#newPose[2] = oldAng0

			newDist1 = pathSpline.dist_u(newU1)
			newP1 = pathSpline.point_u(newU1)
			newAng1 = pathSpline.angle_u(newU1)
			oldAng1 = hypPose1[2]
			newPose1 = copy(newP1)
			#newPose[2] = oldAng0

			#particleObj = Particle(newPose0, newPose1, 0, newDist0, ('t0', 't1'), 0.0)
			particleObj = part.copy()
			#particleObj.pose0 = newPose0
			#particleObj.pose1 = newPose1
			particleObj.displacePose(newPose0, newPose1)
			particleObj.newDist0 = newDist0

			#uPath0, uMedialOrigin0 = selectLocalCommonOrigin(currSplice0, medial0, newPose)

			newPartDist2.append(particleObj)


		self.poseParticles["snapshots2"][0] = newPartDist2

	#@logFunction
	def localizePoseParticles(self, nodeID0, nodeID1):

		#cProfile.run('runTest(probe)', 'test_prof')

		poseData = self.poseData

		""" if the node has a landmark feature, then we try to localize on the junction point """
		""" FIXME: only use the first long path, we don't try and figure out the longest yet """

		if False:
		
			for nodeID in [nodeID0, nodeID1]:

				sF0 = poseData.spatialFeatures[nodeID][0]
				if len(self.pathClasses) > 1 and (sF0["bloomPoint"] != None or sF0["archPoint"] != None):

					""" FIXME: assume one landmark feature """

					""" get the point landmark in local coordinates """
					if sF0["bloomPoint"] != None:
						print "DO bloom localize to landmark feature", nodeID
						localJunctionPoint = sF0["bloomPoint"]

					elif sF0["archPoint"] != None:
						print "DO arch localize to landmark feature", nodeID
						localJunctionPoint = sF0["archPoint"]

					""" starting pose of nodeID """
					nodePose0 = self.getNodePose(nodeID)
					
					""" convert landmark into global coordinates """
					poseOrigin0 = Pose(nodePose0)
					globalLandmarkPoint = poseOrigin0.convertLocalToGlobal(localJunctionPoint)

	
					""" identify which junction this landmark is by the closest in cartesian distance from canonical """
					minPathID = None
					minJuncDist = 1e100
					junctionPose0 = None
					for pathID, values in self.pathClasses.iteritems():

						if pathID != 0 and values["branchNodeID"] != nodeID:

							currJuncPose = values["globalJunctionPose"]

							currDist = sqrt((currJuncPose[0]-globalLandmarkPoint[0])**2 + (currJuncPose[1]-globalLandmarkPoint[1])**2)

							if currDist < minJuncDist:
								minJuncDist = currDist
								minPathID = pathID
								junctionPose0 = currJuncPose



					""" compute new node pose based on current location of junction landmark feature """
					if minPathID != None:

						print "nodePose0 =", nodePose0
						print "global landmark =", globalLandmarkPoint
						print "localJunctionPoint =", localJunctionPoint
						print "junction pose =", junctionPose0
						print "nodeID, pathID =", nodeID, minPathID

						modJuncPose0 = [junctionPose0[0], junctionPose0[1], nodePose0[2]]


						juncOrigin0 = Pose(modJuncPose0)
						invPoint0 = [-localJunctionPoint[0], -localJunctionPoint[1]]
						newNodePoint0 = juncOrigin0.convertLocalToGlobal(invPoint0)

						""" set pose of node0 to compute offset, but keep same angle """
						junctionNodePose0 = [newNodePoint0[0], newNodePoint0[1], nodePose0[2]]
						
						print "junctionNodePose0 =", junctionNodePose0

						self.setNodePose(nodeID, junctionNodePose0)

						
		self.isNoLocalize = False
		self.resetBranches()

		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]

		nodePose0 = self.getNodePose(nodeID0)
		staticSplicedPaths, spliceTerms, staticSplicePathIDs = self.getSplicesByNearJunction(nodePose0)

		print "numStaticSplicedPaths =", len(staticSplicedPaths)

		pathIDs = self.getPathIDs()

		medial0 = poseData.medialAxes[nodeID0]
		medial1 = poseData.medialAxes[nodeID1]

		if nodeID0-2 >= 0:
			prevMedial0 = poseData.medialAxes[nodeID0-2]
			prevMedial1 = poseData.medialAxes[nodeID1-2]
		else:
			prevMedial0 = None
			prevMedial1 = None


		medialSpline0 = SplineFit(medial0, smooth=0.1)
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		minMedialDist0, oldMedialU0, oldMedialP0 = medialSpline0.findClosestPoint([0.0,0.0,0.0])
		minMedialDist1, oldMedialU1, oldMedialP1 = medialSpline1.findClosestPoint([0.0,0.0,0.0])


		localizeJobs = []
		allSplicedPaths = []
		spliceCount = 0

		#parentPathIDs = self.getParentHash()
		#controlPoses = self.getControlPoses()
		#globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		probSets = []
		for particleIndex in range(len(particleDist2)):
			part = particleDist2[particleIndex]
			hypPose0 = part.pose0
			pathIDs = self.getPathIDs()

			for pathID in pathIDs:
				if pathID != 0:

					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]
					controlPose = part.junctionData[pathID]["controlPose"]

					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)
					if dist1 < 3.0:
						probSets.append((pathID, globJuncPose, controlPose))


		""" precompute the evaluation for branches and cache result """
		""" each branch point is the max likelihood of each pose particle """
		self.batchEvalBranches(probSets)


		for particleIndex in range(len(particleDist2)):

			time1 = time.time()

			part = particleDist2[particleIndex]

			hypPose0 = part.pose0
			hypPose1 = part.pose1
			#pathID = part.memberPaths[0]
			hypDist = part.hypDist
			prevHypPose0 = part.prevPose0
			prevHypPose1 = part.prevPose1

			resultsBySplice = []

			thisSplicedPaths = []

			for pathID in pathIDs:

				if pathID != 0:

					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]
					controlPose = part.junctionData[pathID]["controlPose"]

					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)
					if dist1 < 3.0:

						#time1 = time.time()

						samples = self.evaluateBranches(pathID, controlPose, globJuncPose, nodeID0, particleIndex)

						probDist = []
						branchPoseDist = []
						controlPoseDist = []
						branchSplices = []
						for branchSamp in samples:
							"""
							8 = initial probability value of this branch position, initialized to 0.0 in this function
							1 = pose of the branch point
							9 = set of splices including the branch at this position
							"""
							probDist.append(branchSamp["initProb"])
							branchPoseDist.append(branchSamp["modJuncPose"])
							controlPoseDist.append(branchSamp["modControlPose"])
							branchSplices.append(branchSamp["newSplices"])

							print "initProb, modJuncPose, modControlPose:", branchSamp["initProb"], branchSamp["modJuncPose"], branchSamp["modControlPose"]


						probStr = ""
						for probVal in probDist:
							probStr += "%1.2f " % probVal

						print self.hypothesisID, particleIndex, pathID, "branch probDist:", probStr

						#self.poseParticles["snapshots2"][updateCount][particleIndex].junctionData[pathID] = {}
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["probDist"] = probDist
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["branchPoseDist"] = branchPoseDist
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["controlPoseDist"] = controlPoseDist
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["branchSplices"] = branchSplices


						#time2 = time.time()
						#print "evaluateBranches time", self.hypothesisID, particleIndex, pathID, time2-time1

						for j in range(len(samples)):

							branchPose = samples[j]["modJuncPose"]
							probVal = samples[j]["initProb"]
							#branchSplices = samples[j][9]
							numSplices = len(branchSplices[j])

							dist2 = sqrt((branchPose[0]-hypPose0[0])**2 + (branchPose[1]-hypPose0[1])**2)
							if dist2 < 3.0 and probVal > 0.0:

								for k in range(numSplices):
									splice = branchSplices[j][k]
									thisSplicedPaths.append((j, probVal, splice, pathID))


			""" consider single path splices """
			""" uses hypothesize currPath to localize the trimmedPath """
			""" FIXME: grandparent child paths should be localized in a chain up to the root """
			""" FIXME: do we still need comparison to max likelihood?  should this only be the childless parent? """
			for j in range(len(staticSplicePathIDs)):
				if len(staticSplicePathIDs[j]) == 1:

					pathID = staticSplicePathIDs[j][0]
					canonicalPath = staticSplicedPaths[j]

					if pathID != 0:

						partBranchIndex = part.junctionData[pathID]["maxLikelihoodBranch"]

						""" use the max likelihood branch for our 'static' case """
						partControlPose = part.junctionData[pathID]["controlPoseDist"][partBranchIndex]
						
						origControlPose = copy(self.pathClasses[pathID]["controlPose"])
						origControlPose[2] = 0.0

						modControlPose = copy(partControlPose)
						modControlPose[2] = 0.0

						origControlOrigin = Pose(origControlPose)
						offsetOrigin1 = Pose(modControlPose)


						localPath = []
						for p in canonicalPath:
							p1 = origControlOrigin.convertGlobalToLocal(p)
							localPath.append(p1)


						placedPath = []
						for p in localPath:
							p1 = offsetOrigin1.convertLocalToGlobal(p)
							placedPath.append(p1)

						""" static splices have 1.0 probability """
						thisSplicedPaths.append((None, 1.0, placedPath, pathID))
							
					else:

						""" static splices have 1.0 probability """
						thisSplicedPaths.append((None, 1.0, canonicalPath, pathID))


			print "particle:", particleIndex, ",",  len(thisSplicedPaths), "localize jobs"
			for spliceIndex in range(len(thisSplicedPaths)):
				
				branchSampleIndex = thisSplicedPaths[spliceIndex][0]
				probVal = thisSplicedPaths[spliceIndex][1]
				path = thisSplicedPaths[spliceIndex][2]

				
				localizeJobs.append([oldMedialP0, oldMedialU0, 0.0, oldMedialP1, oldMedialU1, 0.0, branchSampleIndex, spliceCount, path, medial0, medial1, deepcopy(hypPose0), deepcopy(hypPose1), prevMedial0, prevMedial1, prevHypPose0, prevHypPose1, [], nodeID0, nodeID1, particleIndex, updateCount, self.hypothesisID])

				self.pathPlotCount2 += 1
				spliceCount += 1



			allSplicedPaths += thisSplicedPaths

			time2 = time.time()
			print "localize construct", particleIndex, time2-time1
			print

		print len(localizeJobs), "total localize jobs"

		time1 = time.time()

		results = batchLocalizeParticle(localizeJobs)

		time2 = time.time()

		print "batch localize for map hypothesis", self.hypothesisID, "=", time2-time1 
		print len(results[0]), "result arguments"
		#for k in range(len(results)):
		#	resultArgs = results[k]


		#results.sort
		#sortedResults = sorted(results, key=itemgetter(0,1), reverse=False)

		""" sort by pose particle index followed by utility value """
		sortedResults = sorted(results, key=itemgetter(0,45), reverse=False)

		print "particle sort"
		thisParticleID = -1
		distMin = 1e100
		filteredParticles = []
		for res in sortedResults:
			if res[0] != thisParticleID:
				thisParticleID = res[0]
				distMin = res[1]
				filteredParticles.append(res)

		print "len(filteredParticles) =", len(filteredParticles)

		"""
		result2["parentID"] = result["parentID"]
		result2["modJuncPose"] = result["modJuncPose"]
		result2["modControlPose"] = result["modControlPose"]
		result2["newArcDist"] = result["newArcDist"]
		result2["pathID"] = result["pathID"]
		result2["matchCount"] = result["matchCount"]
		result2["lastCost"] = result["lastCost"]
		result2["distDisc"] = result["distDisc"]
		result2["angDisc"] = result["angDisc"]
		result2["initProb"] = probVal
		result2["newSplices"] = result["newSplices"]
		result2["juncDiscAngle"] = result["juncDiscAngle"]
		result2["juncDiscDist"] = result["juncDiscDist"]
		"""

		newParticleDist2 = []

		print "particle evaluation:", nodeID0, self.hypothesisID, updateCount
		for particleIndex in range(len(filteredParticles)):

			#if duplicateMatches[particleIndex] == particleIndex:

			part = filteredParticles[particleIndex]

			utilVal = part[45]
			spliceIndex = part[47]
			branchIndex = allSplicedPaths[spliceIndex][0]
			branchProbVal = allSplicedPaths[spliceIndex][1]
			spliceCurve = allSplicedPaths[spliceIndex][2]
			pathID = allSplicedPaths[spliceIndex][3]

			initPose0 = part[48]
			initPose1 = part[49]

			overlapSum = part[50]


			newPose0 = part[2]
			newDist0 = 0.5

			isInterior1_0 = part[9]
			isExist1_0 = part[10]
			isInterior2_0 = part[15]
			isExist2_0 = part[16]
			contigFrac_0 = part[19]
			overlapSum_0 = part[20]

			# return departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2
			# return (particleIndex, icpDist0, resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs0 + (isExist1_0 or isExist2_0,) + (icpDist1, resultPose1, lastCost1, matchCount1, currAng1, currU1) + resultArgs1 + (isExist1_1 or isExist2_1,)

			newPose1 = part[24]
			isInterior1_1 = part[31]
			isExist1_1 = part[32]
			isInterior2_1 = part[37]
			isExist2_1 = part[38]
			contigFrac_1 = part[41]
			overlapSum_1 = part[42]


			isReject = False
			""" divergence is prohibited """
			if isInterior1_0 or isInterior2_0 or isInterior1_1 or isInterior2_1:
				isReject = True

			""" only a single extension is permitted, but not two """
			#if isExist1_0 and isExist2_0 or isExist1_1 and isExist2_1:
			if isExist1_0 and isExist2_0 or isExist1_1 and isExist2_1:
				isReject = True
			
			""" horrifically low contiguity is rejected out of hand """
			if contigFrac_0 <= 0.5 or contigFrac_1 <= 0.5:
				isReject = True

			""" contiguity between current and previous pose """
			if overlapSum > 1e10:
				isReject = True

			if fabs(diffAngle(initPose0[2],newPose0[2])) > 2.0*pi/3.0:
				isReject = True
				
			""" probability is the contigFrac squared """ 
			newProb = 0.0
			if not isReject:
				newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0 * contigFrac_0 / overlapSum_0
				newProb *= branchProbVal
				#utilVal0 = (1.0-contigFrac_0) + (isExist1_0 or isExist2_0) + (1.0-contigFrac_1) + (isExist1_1 or isExist2_1)
			
			print "%d %d %d branchProbVal, utilVal, poseProbValue, overlapSum:" % (self.hypothesisID, particleIndex, spliceIndex), branchProbVal, utilVal, newProb, overlapSum, contigFrac_0, contigFrac_1, overlapSum, initPose0[2], newPose0[2], int(isExist1_0), int(isExist2_0), int(isExist1_1), int(isExist2_1), int(isInterior1_0), int(isInterior2_0), int(isInterior1_1), int(isInterior2_1), isReject, branchIndex

			#print "%d %1.2f %1.2f %1.2f %d %d %d %d %1.2f %1.2f %d %d %d %d" % (particleIndex, newProb, contigFrac_0, overlapSum_0, isInterior1_0, isInterior2_0, isExist1_0, isExist2_0, contigFrac_1, overlapSum_1, isInterior1_1, isInterior2_1, isExist1_1, isExist2_1)

			" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 "
			" return (particleIndex, icpDist, resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs + (isExist1 or isExist2,) "

			particleObj = particleDist2[particleIndex].copy()

			#particleObj.updatePose(pose0, pose1)
			particleObj.pose0 = newPose0
			particleObj.pose1 = newPose1
			particleObj.hypDist = newDist0
			particleObj.weightVal = newProb
			particleObj.spliceCurve = spliceCurve

			""" Change branch if we are within a junction, otherwise, keep the most recent branch index """
			if branchIndex != None:
				particleObj.junctionData[pathID]["maxLikelihoodBranch"] = branchIndex


			newParticleDist2.append(particleObj)

			#else:
			#	dupIndex = duplicateMatches[particleIndex]
			#	newParticleDist2.append(newParticleDist2[dupIndex].copy())


		updateCount += 1
		self.poseParticles["updateCount"] = updateCount
		self.poseParticles["snapshots2"][0] = newParticleDist2
		numParticles = self.poseParticles["numParticles"]

		self.drawPoseParticles()

		""" now resample the particles """
		particleDist2 = self.poseParticles["snapshots2"][0]
		probSum = 0.0
		for part in particleDist2:
			probVal = part.weightVal
			probSum += probVal

		probParticles = []
		for k in range(len(particleDist2)):
			part = particleDist2[k]

			poseProbVal = part.weightVal

			""" if all 0.0 probabilities, make uniform distribution, set no localization flag """
			if probSum > 0.0:
				probParticles.append(poseProbVal/probSum)
			else:
				probParticles.append(float(1)/float(numParticles))
				self.isNoLocalize = True

		print "probParticles:", self.isNoLocalize, probParticles 

		""" apply branch probability """
		probSum2 = 0.0
		for k in range(len(probParticles)):
			#part = particleDist2[k]
			#part = filteredParticles[particleIndex]
			part = filteredParticles[k]

			""" branch probability """
			spliceIndex = part[47]
			branchProbVal = allSplicedPaths[spliceIndex][1]

			probParticles[k] *= branchProbVal

			probSum2 += probParticles[k]

		print "applied branch probability:", probParticles

		""" renormalize """
		for k in range(len(probParticles)):
			if probSum2 > 0.0:
				probParticles[k] /= probSum2
			else:
				probParticles[k] = float(1)/float(numParticles)

		print "normalized probParticles:", probParticles 

		""" now resample """
		resampledParticles2 = []
		numParticles = self.poseParticles["numParticles"]

		""" find the maximum likelihood particle.  Take average if tied for max """
		maxVal = 0.0
		maxIndex = 0
		numMax = 0
		for k in range(len(probParticles)):
			if probParticles[k] > maxVal:
				maxVal = probParticles[k]
				maxIndex = k 
				numMax = 1
			elif probParticles[k] == maxVal:
				numMax += 1

		print "maxIndex, numMax, maxVal:", maxIndex, numMax, maxVal

		""" if more than one max, than we don't change the node pose, we accept the localize on the canonical pose """
		if numMax == 1:

			self.setNodePose(nodeID0, deepcopy(particleDist2[maxIndex].pose0))
			self.setNodePose(nodeID1, deepcopy(particleDist2[maxIndex].pose1))

			print "setting max likelihood poses", nodeID0, nodeID1, self.getNodePose(nodeID0), self.getNodePose(nodeID1)

			#self.nodePoses[nodeID0] = deepcopy(particleDist2[maxIndex].pose0)
			#self.nodePoses[nodeID1] = deepcopy(particleDist2[maxIndex].pose1)


			""" change to the maximum likelihood branch position as well """

			allPathIDs = self.getPathIDs()
			for pathID in allPathIDs:
				if pathID != 0:

					""" maximum likelihood pose particle """
					part = particleDist2[maxIndex]
					#partBranchIndex = particleDist2[maxIndex].junctionData[pathID]["maxLikelihoodBranch"]


					""" maximum likelihood branch point within maximum likelihood pose particle """
					partBranchIndex = part.junctionData[pathID]["maxLikelihoodBranch"]
					partControlPose = part.junctionData[pathID]["controlPoseDist"][partBranchIndex]

					#origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
					origControlPose = copy(self.pathClasses[pathID]["controlPose"])
					print "changing path", pathID, "branch position from", origControlPose, "to", partControlPose
					origControlPose[2] = 0.0

					modControlPose = copy(partControlPose)
					modControlPose[2] = 0.0

					origControlOrigin = Pose(origControlPose)
					offsetOrigin1 = Pose(modControlPose)

					memberNodes = self.getNodes(pathID)

					for nodeID in memberNodes:
						if nodeID != nodeID0 and nodeID != nodeID1:
							localPose = origControlOrigin.convertGlobalPoseToLocal(self.getNodePose(nodeID))
							newGlobalPose = offsetOrigin1.convertLocalOffsetToGlobal(localPose)
							#self.nodePoses[nodeID] = newGlobalPose
							self.setNodePose(nodeID, newGlobalPose)

					""" update of canonical branch point from maximum likelihood branch point """
					partJuncPose = part.junctionData[pathID]["branchPoseDist"][partBranchIndex]
					newJuncPose = deepcopy(partJuncPose)
					newJuncPose[2] = self.pathClasses[pathID]["globalJunctionPose"][2]
					self.pathClasses[pathID]["globalJunctionPose"] = newJuncPose

					self.pathClasses[pathID]["controlPose"] = partControlPose


		elif numMax > 1:
			pass

		""" FIXME:  do we really do nothing here when we tie for maximum? """


		print "particle evaluation:", nodeID0, self.hypothesisID, updateCount
		for i in range(numParticles):
			sampVal = random.random()

			k = 0
			probSum = probParticles[k]

			while sampVal > probSum:
				k += 1
				probSum += probParticles[k]

			oldPart2 = particleDist2[k]
			
			myKeys = oldPart2.junctionData.keys()
			printStr = ""
			for key in myKeys:
				printStr += repr(oldPart2.junctionData[key]["globalJunctionPose"]) 
				printStr += " " + repr(self.pathClasses[key]["globalJunctionPose"]) 

			#print "resampling particle", k, probParticles[k], printStr
			print "resampling particle", k, probParticles[k]
			print oldPart2.prevPose0, oldPart2.prevPose1
			print oldPart2.dispPose0, oldPart2.dispPose1
			print oldPart2.pose0, oldPart2.pose1

			resampledParticles2.append(oldPart2.copy())

		updateCount += 1
		self.poseParticles["updateCount"] = updateCount
		self.poseParticles["snapshots2"][0] = resampledParticles2

		self.drawPoseParticles()

		return


	@logFunction
	def drawPoseParticles(self):

		updateCount = self.poseParticles["updateCount"] 
		#particleDist = self.poseParticles["snapshots"][updateCount]
		particleDist2 = self.poseParticles["snapshots2"][0]


		localPathSets = {}

		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:

			if pathID != 0:

				junctionDetails = self.leaf2LeafPathJunctions[pathID]

				origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
				origJuncPose[2] = 0.0
				origJuncOrigin = Pose(origJuncPose)

				#smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]
				localPathSegs = junctionDetails["localSegments"]
				#for k in range(len(smoothPathSegs)):
				#	pathSeg = smoothPathSegs[k]
				#	localSeg = []
				#	for p in pathSeg:
				#		p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
				#		localSeg.append(p1)
				#
				#	localPathSegs.append(localSeg)

				#print "path segs:", len(junctionDetails["leafSegments"]), "+", len(junctionDetails["internalSegments"]), "=", len(smoothPathSegs)

				localPathSets[pathID] = localPathSegs


		#print "particles:", particleDist2

		pylab.ioff()
		print pylab.isinteractive()
		pylab.clf()

		for k,path in  self.trimmedPaths.iteritems():
			print "path has", len(path), "points"
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color = self.colors[k], linewidth=4)


		#self.drawWalls()

		"""
		xP = []
		yP = []
		for p in orientedPoints2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k')

		xP = []
		yP = []
		for p in orientedPoints3:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='k')
		"""

		nodeID0 = self.poseData.numNodes-2
		nodeID1 = self.poseData.numNodes-1
		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]

		# newPart = (newPose0, newPose1, 0, newDist0, ('t0', 't1'), newProb)
		hypPointsX_0 = []
		hypPointsY_0 = []
		hypPointsX_1 = []
		hypPointsY_1 = []
		for part in particleDist2:
			hypPose0 = part.pose0
			hypPose1 = part.pose1
			#hypPose0 = part[0]
			#hypPose1 = part[1]
			#pathID = part[1]
			#hypDist = part[2]
			#spliceID = part[3]
			hypPointsX_0.append(hypPose0[0])
			hypPointsY_0.append(hypPose0[1])
			hypPointsX_1.append(hypPose1[0])
			hypPointsY_1.append(hypPose1[1])

			thisSplice = part.spliceCurve




			poseOrigin0 = Pose(hypPose0)

			xP = []
			yP = []
			for p in medial0:
				p1 = poseOrigin0.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'r', alpha = 0.2, zorder=9)


			allPathIDs = self.getPathIDs()
			for pathID in allPathIDs:


				xP = []
				yP = []
				if pathID != 0:


					if pathID in part.junctionData.keys():

						probDist = part.junctionData[pathID]["probDist"]
						branchPoseDist = part.junctionData[pathID]["branchPoseDist"]
						controlPoseDist = part.junctionData[pathID]["controlPoseDist"]


						maxProb = -1e100
						maxIndex = 0
						for k in range(len(probDist)):
							modPose0 = deepcopy(controlPoseDist[k])

							xP.append(modPose0[0])
							yP.append(modPose0[1])

							if probDist[k] > maxProb:
								maxProb = probDist[k]
								maxIndex = k

						modPose0 = deepcopy(controlPoseDist[maxIndex])

						print "probDist:", probDist
						print "branchPoseDist:", branchPoseDist
						print "controlPoseDist:", controlPoseDist

					else:
						#modPose0 = deepcopy(part.junctionData[pathID]["globalJunctionPose"])
						modPose0 = deepcopy(part.junctionData[pathID]["controlPose"])

					print "modPose0:", modPose0
					modPose0[2] = 0.0
					modOrigin0 = Pose(modPose0)

					"""
					if thisSplice != None:

						xP1 = []
						yP1 = []

						for p1 in thisSplice:
							xP1.append(p1[0])
							yP1.append(p1[1])

						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.2)
					"""

					"""
					localPathSegs = localPathSets[pathID] 

					for k in range(len(localPathSegs)):
						localSeg = localPathSegs[k]
						xP1 = []
						yP1 = []
						for p in localSeg:
							p1 = modOrigin0.convertLocalOffsetToGlobal(p)
							xP1.append(p1[0])
							yP1.append(p1[1])
						#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part[8]/maxProb)
						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.2)
					"""

					#pylab.scatter([modPose0[0]], [modPose0[1]], color='k', linewidth=1, zorder=10, alpha=0.4)
					pylab.scatter(xP, yP, color='k', linewidth=1, zorder=10, alpha=0.4)




			poseOrigin1 = Pose(hypPose1)

			xP = []
			yP = []
			for p in medial1:
				p1 = poseOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'b', alpha = 0.2, zorder=9)

		#pylab.scatter(hypPointsX_0, hypPointsY_0, color='r', linewidth=1, zorder=10, alpha=0.2)
		#pylab.scatter(hypPointsX_1, hypPointsY_1, color='b', linewidth=1, zorder=10, alpha=0.2)

		pose0 = self.getNodePose(nodeID0)
		pose1 = self.getNodePose(nodeID1)
		medial0 = self.poseData.medialAxes[nodeID0]
		medial1 = self.poseData.medialAxes[nodeID1]

		poseOrigin0 = Pose(pose0)
		poseOrigin1 = Pose(pose1)

		xP = []
		yP = []
		for p in medial0:
			p1 = poseOrigin0.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		pylab.plot(xP,yP,color='k', zorder=9, alpha=0.2)

		xP = []
		yP = []
		for p in medial1:
			p1 = poseOrigin1.convertLocalToGlobal(p)
			xP.append(p1[0])
			yP.append(p1[1])
		pylab.plot(xP,yP,color='k', zorder=9, alpha=0.2)

		
		pylab.xlim(-6,16.48)
		#pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		#pylab.ylim(-5.85,5.85)
		#pylab.axis("equal")
		#pylab.title("%d particles" % len(particleDist))
		pylab.title("nodeIDs %d %d hypID %d" % (nodeID0, nodeID1, self.hypothesisID))


		#pylab.savefig("moveEstimate2_%04u_%04u.png" % (self.posePlotCount, self.hypothesisID))
		pylab.savefig("moveEstimate2_%04u_%04u.png" % (self.hypothesisID, self.posePlotCount))
		print "moveEstimate2_%04u_%04u.png" % (self.hypothesisID, self.posePlotCount)

		self.posePlotCount += 1

	@logFunction
	def updatePoseData(self, poseData):
		self.poseData = deepcopy(poseData)
		print "updated poseData:", poseData.numNodes, self.poseData.numNodes

	@logFunction
	def copy(self, hypothesisID):

		print "creating hypothesis", hypothesisID

		#newObj = deepcopy(self)
		#newObj.hypothesisID = hypothesisID

		#return newObj

		newObj = MapState(self.poseData, hypothesisID)

		newObj.numPoseParticles = deepcopy(self.numPoseParticles)

		newObj.poseParticles = deepcopy(self.poseParticles)

		updateCount = newObj.poseParticles["updateCount"] 
		particleDist2 = deepcopy(newObj.poseParticles["snapshots2"][0])
		for part in particleDist2:
			part.mapStateID = hypothesisID


		newObj.pathClasses = deepcopy(self.pathClasses)
		newObj.pathTermsVisisted = deepcopy(self.pathTermsVisited)
		newObj.pathIDs = deepcopy(self.pathIDs)
		newObj.paths = deepcopy(self.paths)
		newObj.hulls = deepcopy(self.hulls)
		newObj.trimmedPaths = deepcopy(self.trimmedPaths)
		newObj.nodeRawPoses = deepcopy(self.nodeRawPoses)
		newObj.nodePoses = deepcopy(self.nodePoses)
		newObj.gndPoses = deepcopy(self.gndPoses)
		newObj.gndRawPoses = deepcopy(self.gndRawPoses)
		newObj.origPoses = deepcopy(self.origPoses)
		newObj.utility = deepcopy(self.utility)
		newObj.isNodeBranching = deepcopy(self.isNodeBranching)

		newObj.alphaPlotCount = self.alphaPlotCount
		newObj.medialCount = self.medialCount
		newObj.topCount = self.topCount
		newObj.spliceCount = self.spliceCount
		newObj.multiDepCount = self.multiDepCount
		newObj.overlapPlotCount = self.overlapPlotCount
			
		newObj.pathPlotCount = self.pathPlotCount
		newObj.pathPlotCount2 = self.pathPlotCount2
		newObj.overlapPlotCount2 = self.overlapPlotCount2
		newObj.medialSouptCount = self.medialSoupCount

		
		newObj.rootPoint = deepcopy(self.rootPoint)
		
		newObj.pathGraph = deepcopy(self.pathGraph)
		newObj.joins = deepcopy(self.joins)
		newObj.junctions = deepcopy(self.junctions)
		newObj.terminals = deepcopy(self.terminals)
		newObj.topDict = deepcopy(self.topDict)
		newObj.allSplices = deepcopy(self.allSplices)

		newObj.orderedPathIDs1 = deepcopy(self.orderedPathIDs1)
		newObj.orderedPathIDs2 = deepcopy(self.orderedPathIDs2)

		" departure events for node 1 "
		newObj.departures1 = deepcopy(self.departures1) 
		newObj.interiors1 = deepcopy(self.interiors1)  
		newObj.depPoints1 = deepcopy(self.depPoints1) 
		newObj.distances1 = deepcopy(self.distances1)  
		newObj.depAngles1 = deepcopy(self.depAngles1)
		newObj.contig1 = deepcopy(self.contig1)

		" departure events for node 2 "
		newObj.departures2 = deepcopy(self.departures2) 
		newObj.interiors2 = deepcopy(self.interiors2) 
		newObj.depPoints2 = deepcopy(self.depPoints2)  
		newObj.distances2 = deepcopy(self.distances2) 
		newObj.depAngles2 = deepcopy(self.depAngles2)  
		newObj.contig2 = deepcopy(self.contig2)  

		newObj.departureResultSet1 = deepcopy(self.departureResultSet1)  
		newObj.departureResultSet2 = deepcopy(self.departureResultSet2)

		newObj.leaf2LeafPathJunctions = deepcopy(self.leaf2LeafPathJunctions) 

		return newObj

	@logFunction
	def computeEval(self):

		keys = self.trimmedPaths.keys()
		keys.sort()

		totalSum = 0.0
		for j in keys:
			path1 = self.trimmedPaths[j]
			for k in keys:
				path2 = self.trimmedPaths[k]

				if j < k:
					resultSum = self.getPathOverlapSum(path1, path2, j, k, plotIter = True)
					print "computing sum of", j, "and", k, "=", resultSum
					totalSum += resultSum

		self.mapOverlapSum = totalSum
		return totalSum
	
		utilSum = 0.0

		for nodeID1 in self.nodePoses.keys():

			estPose1 = self.getNodePose(nodeID1)
	
			orderedPathIDs1 = self.getOrderedOverlappingPaths(nodeID1)
			print nodeID1, "recompute orderedPathIDs1:", orderedPathIDs1

			hull1 = self.poseData.aHulls[nodeID1]
			medial1 = self.poseData.medialAxes[nodeID1]
			origPose1 = self.origPoses[nodeID1]
			#estPose1 = self.nodePoses[nodeID1]
				
			splicedPaths1 = self.splicePathIDs(orderedPathIDs1)

			#print "received", len(splicedPaths1), "spliced paths from path IDs", orderedPathIDs1

			" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 |"

			results1 = []		
			for k in range(len(splicedPaths1)):
				path = splicedPaths1[k]			
				result = getMultiDeparturePoint(path, medial1, origPose1, estPose1, orderedPathIDs1, nodeID1)
				results1.append(result+(k,))

			results1 = sorted(results1, key=itemgetter(14))
			results1 = sorted(results1, key=itemgetter(12), reverse=True)				

			#kIndex = results1[0][15]

			angDiff2 = results1[0][14]
			contigFrac = results1[0][12]
			#dist = results1[0][17]
			dist = sqrt((origPose1[0]-estPose1[0])**2 + (origPose1[1] - estPose1[1])**2)
			
			#isDeparture = results1[0][15]
			isDeparture = results1[0][2] and results1[0][3] or results1[0][8] and results1[0][9]
			#utilVal = (0.5-angDiff2)/0.5 + (contigFrac-0.5)/0.5 + (3.0-dist)/3.0 + 1-isDeparture*0.5

			maxFront = results1[0][5]
			maxBack = results1[0][11]

			#utilVal = contigFrac - 2*isDeparture

			utilVal = maxFront + maxBack

			print "util:", nodeID1, maxFront, maxBack, angDiff2, contigFrac, dist, isDeparture, utilVal
				

			utilSum += utilVal
			#return splicedPaths1[kIndex]

		self.utility = utilSum

		return utilSum


	""" state save and restore methods """
	@logFunction
	def saveState(self, saveCount):
		
		saveFile = ""
		
		saveFile += "self.pathClasses = " + repr(self.pathClasses) + "\n"
		saveFile += "self.pathTermsVisited = " + repr(self.pathTermsVisited) + "\n"
		saveFile += "self.pathIDs = " + repr(self.pathIDs) + "\n"		 
		saveFile += "self.paths = " + repr(self.paths) + "\n"		 
		saveFile += "self.hulls = " + repr(self.hulls) + "\n"		 
		saveFile += "self.trimmedPaths = " + repr(self.trimmedPaths) + "\n"		   

		saveFile += "self.nodeRawPoses = " + repr(self.nodeRawPoses) + "\n"		 
		saveFile += "self.nodePoses = " + repr(self.nodePoses) + "\n"		 
		saveFile += "self.gndPoses = " + repr(self.gndPoses) + "\n"		 
		saveFile += "self.gndRawPoses = " + repr(self.gndRawPoses) + "\n"		 
		saveFile += "self.hypothesisID = " + repr(self.hypothesisID) + "\n"

		f = open("mapStateSave_%04u_%04u.txt" % (self.hypothesisID, saveCount), 'w')
		f.write(saveFile)
		f.close()
		
	@logFunction
	def restoreState(self, dirName, numNodes):
		
		print "loading" + dirName + "/mapStateSave_%04u_%04u.txt" % (self.hypothesisID, numNodes+1)
		f = open(dirName + "/mapStateSave_%04u_%04u.txt" % (self.hypothesisID, numNodes+1), 'r')		 
		saveStr = f.read()
		print saveStr
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')
		
		exec(saveStr)
		

	""" member access methods """

	@logFunction
	def getControlPose(self, pathID):

		controlPose = self.pathClasses[pathID]["controlPose"]

		return controlPose

	@logFunction
	def getGlobalJunctionPose(self, pathID):

		try:
			globJuncPose = self.pathClasses[pathID]["globalJunctionPose"]
		except:
			
			branchNodeID = self.pathClasses[pathID]["branchNodeID"]
			
			if branchNodeID == None:
				return None
			
			localJunctionPose = self.pathClasses[pathID]["localJunctionPose"]
			
			poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
			globJuncPose = poseOrigin.convertLocalOffsetToGlobal(localJunctionPose)

		return globJuncPose
		
	@logFunction
	def getNodes(self, pathID):
		return self.pathClasses[pathID]["nodeSet"]
	
	@logFunction
	def getPathIDs(self):
		return self.pathClasses.keys()

	@logFunction
	def getParentHash(self):

		parentPathIDs = {}
		pathIDs = self.getPathIDs()

		for pathID in pathIDs:
			cPath = self.getPath(pathID)
			parentPathID = cPath["parentID"]
			parentPathIDs[pathID] = parentPathID

		return parentPathIDs

	@logFunction
	def getControlPoses(self):
		controlPoses = {}
		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]

		return controlPoses
	

	@logFunction
	def getFinalControlPose(self, pathID):
		
		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)


		return globalControlPoses[pathID]

	@logFunction
	def getParentPathID(self, pathID):
		parentID = self.pathClasses[pathID]["parentID"]
		return parentID
	
	@logFunction
	def getPath(self, pathID):
		return self.pathClasses[pathID]
	
	def getAllSplices(self, plotIter = False):
		if not self.isChanged:
			print "returning without computing, not isChanged"
			return self.allSplices, self.terminals, self.junctions
		else:
			print "recomputing allSplices"
			self.generatePaths()
			return self.allSplices, self.terminals, self.junctions

	#@logFunction
	def computeAllSplices(self, plotIter = False):

		print "paths:", self.paths
		print "junctions:", self.junctions
		print "terminals:", self.terminals

		#if not self.isChanged:
		#	print "returning without computing, not isChanged"
		#	return self.allSplices, self.terminals, self.junctions
		#else:
		#	print "computing allSplices"



		" determine which paths are leaves "
		
		pathIDs = self.getPathIDs()
		isAParent = {}
		parents = {}
		pathIDGraph = graph.graph()

		for pathID in pathIDs:
			isAParent[pathID] = False
			pathIDGraph.add_node(pathID, [])
		
		for pathID in pathIDs:	  
			parentPathID = self.getPath(pathID)["parentID"]
			parents[pathID] = parentPathID
			if parentPathID != None:
				isAParent[parentPathID] = True
				pathIDGraph.add_edge(pathID, parentPathID)
		
		leafCount = 0
		for pathID in pathIDs:
			if not isAParent[pathID]:
				leafCount += 1


		" build topological graph "
		topGraph = graph.graph()
		
		topNodes = []
					   
		for k, term in self.terminals.iteritems():
			topGraph.add_node(term[0], term[1])
			print "add_node", term[0], term[1]
			topNodes.append(term[0])
		for k, junc in self.junctions.iteritems():
			topGraph.add_node(junc[2], junc[3])
			print "add_node", junc[2], junc[3]
			topNodes.append(junc[2])


		for pathID in pathIDs:
			pathIndex = []
			if pathID == 0:
				term = self.terminals[0]
				graphNodeKey = term[0]
				pathIndex.append(graphNodeKey[1])
				
				term = self.terminals[1]
				graphNodeKey = term[0]
				pathIndex.append(graphNodeKey[1])
			else:
				term = self.terminals[pathID+1]
				graphNodeKey = term[0]
				pathIndex.append(graphNodeKey[1])
				
			for k, junc in self.junctions.iteritems():
				graphNodeKey = junc[2]
				if graphNodeKey[0] == pathID:
					pathIndex.append(graphNodeKey[1])
						
			pathIndex.sort()
			
			print "pathIndex:", pathIndex
			
			for k in range(len(pathIndex)-1):
				print "add_edge", (pathID,pathIndex[k]),(pathID,pathIndex[k+1])
				topGraph.add_edge((pathID,pathIndex[k]),(pathID,pathIndex[k+1]))

			if pathID != 0:
				graphNodeKey = self.junctions[pathID][2]
				print "add_edge", graphNodeKey,(pathID,pathIndex[0])
				topGraph.add_edge(graphNodeKey,(pathID,pathIndex[0]))


		results = {}
		for i in range(len(pathIDs)):
			pathID1 = pathIDs[i]
			for j in range(i,len(pathIDs)):
				
				pathID2 = pathIDs[j]
				
				if pathID1 == pathID2:
					"term1 to term 2"

					termPaths = []
					if pathID1 == 0:
						termIDs = [self.topDict["t%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
					else:
						termIDs = [self.topDict["j%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
					termPaths.append(termIDs)  

					orderedPathIDs = [pathID1]
				else:
					
					startPathID = pathID1
					endPathID = pathID2
		
					shortestPathSpanTree, shortestDist = pathIDGraph.shortest_path(endPathID)
					nextPathID = shortestPathSpanTree[startPathID]
			
					orderedPathIDs = [startPathID, nextPathID]
					while nextPathID != endPathID:
						nextPathID = shortestPathSpanTree[nextPathID]
						orderedPathIDs.append(nextPathID)			 

					
					" start at first path ID "
					termPaths = []
					if parents[orderedPathIDs[0]] == orderedPathIDs[1]:
						" going to parent "
						termIDs = [self.topDict["t%u" % (orderedPathIDs[0]+1)]]
						termPaths.append(termIDs)
					
					elif parents[orderedPathIDs[1]] == orderedPathIDs[0]:
						" going to child "
						termIDs = [self.topDict["t%u" % (orderedPathIDs[0]+1)]]
						termPaths.append(termIDs)

						if orderedPathIDs[0] == 0:
							termIDs = [self.topDict["t%u" % 0]]
							termPaths.append(termIDs)
						else:
							termIDs = [self.topDict["j%u" % orderedPathIDs[0]]]
							termPaths.append(termIDs)
					else:
						print parents
						print orderedPathIDs
						print termPaths
						raise
						
					" transitions between middle path IDs "
					for k in range(len(orderedPathIDs)-1):
						if parents[orderedPathIDs[k]] == orderedPathIDs[k+1]:
							" child -> parent:	prev -> junction "
							for ids in termPaths:
								ids.append(self.topDict["j%u" % orderedPathIDs[k]])

						elif parents[orderedPathIDs[k+1]] == orderedPathIDs[k]:
							" parent -> child: prev -> junction "
							for ids in termPaths:
								ids.append(self.topDict["j%u" % orderedPathIDs[k+1]])

						else:
							print parents
							print orderedPathIDs
							print termPaths
							print k
							raise
					
					" last path ID "
					if parents[orderedPathIDs[-2]] == orderedPathIDs[-1]:
						" child -> parent:	prev -> junction "
						oldTermPaths = termPaths
						termPaths = []
						" split into 2 possibilities into parent path ID "
						for ids in oldTermPaths:
							ids1 = deepcopy(ids)
							ids2 = deepcopy(ids)
							if orderedPathIDs[-1] == 0:
								ids1.append(self.topDict["t%u" % orderedPathIDs[-1]])
								ids2.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])								  
							else:
								ids1.append(self.topDict["j%u" % orderedPathIDs[-1]])
								ids2.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])
							
							termPaths.append(ids1)
							termPaths.append(ids2)

					elif parents[orderedPathIDs[-1]] == orderedPathIDs[-2]:
						" parent -> child: prev -> junction "
						for ids in termPaths:
							ids.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])
										
					else:
						print parents
						print orderedPathIDs
						print termPaths
						raise
				
				finalResults = []
				
				for termPath in termPaths:

					joinPairs = []

					startKey = termPath[0]
					endKey = termPath[-1]

					""" use splice skeleton to find shortest path """
					startPose = self.pathGraph.get_node_attributes(startKey)
					endPose = self.pathGraph.get_node_attributes(endKey)

					minStartDist = 1e100
					minStartNode = None
					minEndDist = 1e100
					minEndNode = None
					
					for edge in self.spliceSkeleton.edges():
					
						globalNodePoint1 = edge[0]
						globalNodePoint2 = edge[1]

						dist1 = sqrt((globalNodePoint1[0]-startPose[0])**2 + (globalNodePoint1[1]-startPose[1])**2)
						dist2 = sqrt((globalNodePoint2[0]-startPose[0])**2 + (globalNodePoint2[1]-startPose[1])**2)

						if dist1 < minStartDist:
							minStartDist = dist1
							minStartNode = globalNodePoint1

						if dist2 < minStartDist:
							minStartDist = dist2
							minStartNode = globalNodePoint2

						dist1 = sqrt((globalNodePoint1[0]-endPose[0])**2 + (globalNodePoint1[1]-endPose[1])**2)
						dist2 = sqrt((globalNodePoint2[0]-endPose[0])**2 + (globalNodePoint2[1]-endPose[1])**2)

						if dist1 < minEndDist:
							minEndDist = dist1
							minEndNode = globalNodePoint1

						if dist2 < minEndDist:
							minEndDist = dist2
							minEndNode = globalNodePoint2

			
					startNode = minStartNode
					endNode = minEndNode


					shortestSpliceTree, shortestSpliceDist = self.spliceSkeleton.shortest_path(endNode)
					currNode = shortestSpliceTree[startNode]					 
					splicedSkel = []
					while currNode != endNode:
						splicedSkel.append(currNode)
						nextNode = shortestSpliceTree[currNode]
						currNode = nextNode
					splicedSkel.append(currNode)


					""" find splice using the old method """
					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					splicedPath = []
					while currNode != endKey:
						splicedPath.append(self.pathGraph.get_node_attributes(currNode))
						nextNode = shortestPathSpanTree[currNode]
						
						for join in self.joins:
							if join[0] == currNode and join[1] == nextNode:
								joinPairs.append((len(splicedPath)-1,len(splicedPath)))
							elif join[1] == currNode and join[0] == nextNode:
								joinPairs.append((len(splicedPath)-1,len(splicedPath)))
							
						currNode = nextNode

					splicedPath.append(self.pathGraph.get_node_attributes(currNode))


					" if there are any joins, we must interpolate a smooth transition "
					lastIndex = 0
					newPath = []
					for pair in joinPairs:
						index1 = pair[0]
						index2 = pair[1]
						cPoints = [splicedPath[index1-1], splicedPath[index1], splicedPath[index2], splicedPath[index2+1]]				 
						spline = SplineFit(cPoints, smooth=0.0, kp=3)
						points = spline.getUniformSamples(spacing = 0.1)
						
						newPath += splicedPath[lastIndex:index1]
						newPath += points
						
						lastIndex = index2+2

					newPath += splicedPath[lastIndex:]


					sPath = {}
					sPath['orderedPathIDs'] = orderedPathIDs
					sPath['path'] = splicedSkel
					sPath['oldPath'] = newPath
					sPath['termPath'] = termPath
					sPath['skelPath'] = splicedSkel
					
					finalResults.append(sPath)
					
				results[(pathID1,pathID2)] = finalResults
		
		print len(results), "results:"
		for k, result in results.iteritems():
			print k, ":"
			for sPath in result:
				print sPath['orderedPathIDs'], sPath['termPath'], len(sPath['path'])


		if plotIter:
			print "plotting splicedPath"
			pylab.clf()

			
			#for k,path in  self.trimmedPaths.iteritems():
			#	print "path has", len(path), "points"
			#	xP = []
			#	yP = []
			#	for p in path:
			#		xP.append(p[0])
			#		yP.append(p[1])
			#	
			#	pylab.plot(xP,yP, color = self.colors[k])
	
	
			for k, result in results.iteritems():
				for sPath in result:
					#path = sPath['path']
					path = sPath['skelPath']
	
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])
					pylab.plot(xP,yP, color='b', zorder=2)


			for k in pathIDs:
				xP = []
				yP = []
				
				for p in self.paths[k]:
					xP.append(p[0])
					yP.append(p[1])
				pylab.plot(xP,yP, color=self.colors[k], linewidth=1, zorder=1)

			xP = []
			yP = []
			for k,path in  self.pathClasses.iteritems():
				if k != 0:
					globalJunctionPose = path["globalJunctionPose"]
					xP.append(globalJunctionPose[0])
					yP.append(globalJunctionPose[1])

			pylab.scatter(xP,yP, color='k', zorder=3)

			#self.drawWalls()
	
			print "saving splicedPath"
			pylab.title("Spliced Paths, pathIDs = %s" %  self.getPathIDs())
			pylab.savefig("splicedPath_%04u_%04u.png" % (self.hypothesisID, self.spliceCount))
			self.spliceCount += 1

		self.allSplices = results

		#self.isChanged = False
		
		return self.allSplices, self.terminals, self.junctions


	@logFunction
	def isNodeExist(self, nodeID, pathID):
		return nodeID in self.pathClasses[pathID]["nodeSet"]

	
	@logFunction
	def localizeBranchPoint(self, orientedPathSoup, junctionPose1, parentPathID, pathID1, plotIter = False):

		parentPath = self.paths[parentPathID]

		juncOrigin1 = Pose(junctionPose1)
		
		" curves of both paths "

		globalPath2 = parentPath

		globalSpline2 = SplineFit(globalPath2, smooth=0.1)
		globalSamples2 = globalSpline2.getUniformSamples(spacing = 0.04)
		originU2 = globalSpline2.findU(junctionPose1)
		minDist, originU2, uPoint = globalSpline2.findClosestPoint(junctionPose1)
		
		u2 = originU2
		u1 = 0.5  # value is meaningless since data is not a spline curve
		#angGuess = 0.0
		#angGuess = -uPoint[2] + junctionPose1[2]
		angGuess = junctionPose1[2]
		
		resultPose1, lastCost1, matchCount1 = gen_icp.branchEstimateICP([u1,u2,angGuess], junctionPose1, orientedPathSoup, globalPath2, plotIter = plotIter, n1 = pathID1, n2 = parentPathID)

		print "matchCount,cost,result:", matchCount1, lastCost1, resultPose1
		
		return resultPose1, lastCost1, matchCount1
	 

	@logFunction
	def resetBranches(self):


		self.branchEvaluations = {}

		#self.DIV_LEN = 0.2

		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			if self.pathClasses[pathID]["parentID"] != None:
				parentPathID = self.pathClasses[pathID]["parentID"]
				pathSpline = SplineFit(self.paths[parentPathID])

				totalDist = pathSpline.dist_u(1.0)

				branchSpace = {}

				currDist = 0.0
				branchSpace[currDist] = None
				while currDist <= totalDist:
					currDist += self.DIV_LEN
					branchSpace[currDist] = None

				self.branchEvaluations[pathID] = branchSpace

	@logFunction
	def getPrecomputedBranch(self, pathID, arcDist):

		"""
		result values are as follows:
		0 = parent ID
		1 = pose of the branch point
		2 = arc distance on the parent curve of the branch point
		3 = child ID
		4 = match count from ICP algorithm
		5 = last cost from ICP algorithm
		6 = discrepancy distance of new pose from old pose, zeroed for this function
		7 = angle discrepancy difference of new pose from old pose
		8 = initial probability value of this branch position, initialized to 0.0 in this function
		9 = set of splices including the branch at this position
		10 = angle discrepancy 
		11 = dist discrepancy 
		"""

		if self.pathClasses[pathID]["parentID"] != None: 
			currBranchSpace = self.branchEvaluations[pathID]

			currKeys = currBranchSpace.keys()
			currKeys.sort()


			" find the bin for this arc distance "
			thisKey = currKeys[0]
			for k in range(len(currKeys)):
				val = currKeys[k]
				if val >= arcDist:
					break
				thisKey = val

			#print "precompute value:", arcDist, thisKey

			if currBranchSpace[thisKey] != None:
				#print "returned precomputation of branch", thisKey
				return deepcopy(currBranchSpace[thisKey])
			else:

				print "returned precomputation of branch", arcDist, thisKey, None

			raise

			return None


	@logFunction
	def batchEvalBranches(self, probSets):

		""" precompute the evaluation for branches and cache result """
		""" each branch point is the max likelihood of each pose particle """
		""" computes discretized samples around each maximum likelihood branch point """

		""" probSets:  [ (pathID, globJuncPose), ] """
		""" possible locations of branch point shoot pathID """

		probSets = sorted(probSets, key=itemgetter(0), reverse=True)

		allPathIDs = self.getPathIDs()

		""" get spline curves of each parent shoot """
		pathSplines = {}
		for pathID in allPathIDs:

			""" path information """
			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			if parentID != None:
				pathSpline = SplineFit(self.paths[parentID])
				pathSplines[parentID] = pathSpline



		""" convert branch points from poses into arc distances on parent shoot """
		arcDists = []
		for prob in probSets:

			pathID = prob[0]
			estJuncPose = prob[1]
			controlPose = prob[2]

			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]
			pathSpline = pathSplines[parentID]

			""" get the arc distance of the control point """
			minDist, controlUVal, newControlPose = pathSpline.findClosestPoint(controlPose)
			arcDist = pathSpline.dist_u(controlUVal)

			#minDist, uVal, splinePoint = pathSpline.findClosestPoint(estJuncPose)
			#arcDist = pathSpline.dist_u(uVal)


			arcHigh = arcDist + self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
			arcLow = arcDist - self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)

			" sample points around each probSet branch point "
			for k in range(self.NUM_BRANCHES):
				newArcDist = arcLow + k * self.DIV_LEN
				arcDists.append((pathID, newArcDist))


		""" binned arc distances, discretized state space """
		binnedProblems = []

		for prob in arcDists:
			pathID = prob[0]
			arcDist = prob[1]

			currBranchSpace = self.branchEvaluations[pathID]
			currKeys = currBranchSpace.keys()
			currKeys.sort()

			" find the bin for this arc distance "
			thisKey = currKeys[0]
			for k in range(len(currKeys)):
				val = currKeys[k]
				if val >= arcDist:
					break
				thisKey = val

			" converted arc distances to bin keys "
			binnedProblems.append((pathID, thisKey))

		""" discretized space creates duplicates, remove duplicates """
		numProbs1 = len(binnedProblems)
		binnedProblems = list(set(binnedProblems))
		numProbs2 = len(binnedProblems)

		print "binnedProblems:", numProbs1, numProbs2
		branchJobs = []


		parentIDs = {}
		localSkeletons = {}
		controlPoses = {}
		junctionPoses = {}
		parentPathIDs = {}
		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			localSkeletons[pathID] = self.leaf2LeafPathJunctions[pathID]["skeletonGraph"]
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

			cPath = self.getPath(pathID)
			parentPathID = cPath["parentID"]
			parentPathIDs[pathID] = parentPathID


		for prob in binnedProblems:
			pathID = prob[0]
			arcDist = prob[1]

			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			junctionDetails = self.leaf2LeafPathJunctions[pathID]
			#smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]
			localPathSegs = junctionDetails["localSegments"]

			origControlPose = copy(self.getControlPose(pathID))

			origGlobJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
			childPath = self.paths[pathID]
			parentPath = self.paths[parentID]
			trimmedParent = self.trimmedPaths[parentID]

			branchJobs.append((pathID, parentID, origGlobJuncPose, origControlPose, childPath, parentPath, trimmedParent, localPathSegs, arcDist, localSkeletons, controlPoses, junctionPoses, parentPathIDs, len(self.nodePoses)-1, self.hypothesisID ))


		"""
		result values are as follows:
		0 = parent ID
		1 = pose of the branch point
		2 = pose of the control point
		3 = arc distance on the parent curve of the branch point
		4 = child ID
		5 = match count from ICP algorithm
		6 = last cost from ICP algorithm
		7 = discrepancy distance of new pose from old pose, zeroed for this function
		8 = angle discrepancy difference of new pose from old pose
		9 = initial probability value of this branch position, initialized to 0.0 in this function
		10 = set of splices including the branch at this position
		11 = angle discrepancy 
		12 = dist discrepancy 
		"""

		results = batchBranch(branchJobs)

		""" get the maximum value for each of our features """
		maxCost = -1e100
		maxMatchCount = -1000
		maxAngDiff = -1e100
		maxDist = 10.0

		""" find maximum values for each metric """
		for k in range(len(results)):
			part = results[k]
			matchCount = part["matchCount"]
			lastCost = part["lastCost"]
			dist = part["distDisc"]
			angDiff = part["angDisc"]

			if matchCount > maxMatchCount:
				maxMatchCount = matchCount

			if lastCost > maxCost:
				maxCost = lastCost

			if dist > maxDist:
				maxDist = dist

			if angDiff > maxAngDiff:
				maxAngDiff = angDiff


		probSum = 0.0
		for result in results:
			probVal = result["initProb"]
			probSum += probVal

		for result in results:
			arcDist = result["newArcDist"]
			pathID = result["pathID"]
			probVal = result["initProb"] / probSum
			currBranchSpace = self.branchEvaluations[pathID]

			result2 = {}
			result2["parentID"] = result["parentID"]
			result2["modJuncPose"] = result["modJuncPose"]
			result2["modControlPose"] = result["modControlPose"]
			result2["newArcDist"] = result["newArcDist"]
			result2["pathID"] = result["pathID"]
			result2["matchCount"] = result["matchCount"]
			result2["lastCost"] = result["lastCost"]
			result2["distDisc"] = result["distDisc"]
			result2["angDisc"] = result["angDisc"]
			result2["initProb"] = probVal
			result2["newSplices"] = result["newSplices"]
			result2["juncDiscAngle"] = result["juncDiscAngle"]
			result2["juncDiscDist"] = result["juncDiscDist"]


			#result2 = (result[0], result[1], result[3], result[4], result[5], result[6], result[7], result[8], probVal, result[10], result[11], result[12])

			currBranchSpace[arcDist] = result2

		#part2["parentID"] = parentID
		#part2["modJuncPose"] = modJuncPose
		#part2["modControlPose"] = modControlPose
		#part2["newArcDist"] = newArcDist
		#part2["pathID"] = pathID
		#part2["matchCount"] = matchCount
		#part2["lastCost"] = lastCost1
		#part2["distDisc"] = distDisc
		#part2["angDisc"] = angDisc
		#part2["initProb"] = initProb
		#part2["newSplices"] = newSplices
		#part2["juncDiscAngle"] = juncDiscAngle
		#part2["juncDiscDist"] = juncDiscDist

		maxProb = 0.0
		for result in results:
			probVal = result["initProb"]
			if probVal > maxProb:
				maxProb = probVal
			#probSum += probVal

		if True:

			pylab.clf() 


			for result in results:

				modJuncPose = result["modJuncPose"]
				modControlPose = result["modControlPose"]
				result2["newArcDist"] = result["newArcDist"]
				pathID = result["pathID"]
				matchCount = result["matchCount"]
				lastCost = result["lastCost"]
				probVal = result["initProb"]

				newSplices = result["newSplices"]

				#origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
				#origJuncPose[2] = 0.0
				#origJuncOrigin = Pose(origJuncPose)

				junctionDetails = self.leaf2LeafPathJunctions[pathID]
				#smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

				localPathSegs = junctionDetails["localSegments"]

				#localPathSegs = []
				#for k in range(len(smoothPathSegs)):
				#	pathSeg = smoothPathSegs[k]
				#	localSeg = []
				#	for p in pathSeg:
				#		p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
				#		localSeg.append(p1)
				#
				#	localPathSegs.append(localSeg)

				for path in newSplices:
					print "path has", len(path), "points"
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])

					#pylab.plot(xP,yP, color = 'k', linewidth=4)
					pylab.plot(xP,yP,color='k', zorder=9, alpha=0.10)

				for k,path in self.trimmedPaths.iteritems():
					print "path has", len(path), "points"
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])

					pylab.plot(xP,yP, color = self.colors[k], linewidth=4)

				xP = []
				yP = []
				xP.append(modJuncPose[0])
				yP.append(modJuncPose[1])
				pylab.scatter(xP, yP, color='k', zorder=8)

				offsetOrigin1 = Pose(modControlPose)

				"""
				for k in range(len(localPathSegs)):
					localSeg = localPathSegs[k]
					xP1 = []
					yP1 = []
					for p in localSeg:
						p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
						xP1.append(p1[0])
						yP1.append(p1[1])

					#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.5)
					#if maxCost > 0:
					if maxMatchCount > 0:
						#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=probVal/maxProb)
						#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=lastCost/maxCost)
						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=float(matchCount)/float(maxMatchCount))
					else:
						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.10)
				"""



			pylab.title("nodeID: %d hyp: %d, pathID: %d" % (len(self.nodePoses), self.hypothesisID, pathID))
			pylab.savefig("all_branch_plot_%04u_%02u_%04u.png" % (len(self.nodePoses), self.hypothesisID, self.tempCount) )
			self.tempCount += 1


	@logFunction
	def evaluateBranches(self, pathID, modControlPose, estJuncPose, nodeID0, particleIndex):

		# pathID
		# parentID
		# curr path
		# parent path
		# origJuncPose
		# junctionDetails = self.leaf2LeafPathJunctions[pathID]
		# smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]
		# hypID


		if self.pathClasses[pathID]["parentID"] != None: 

			""" path information """
			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			pathSpline = SplineFit(self.paths[parentID])

			""" initial position of junction """
			origJuncPose = copy(pathDesc["globalJunctionPose"])
			origJuncPose[2] = 0.0

			#minDist, uVal, splinePoint = pathSpline.findClosestPoint(estJuncPose)
			minDist, uVal, splinePoint = pathSpline.findClosestPoint(modControlPose)
			arcDist = pathSpline.dist_u(uVal)

			arcHigh = arcDist + self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
			arcLow = arcDist - self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)

			totalDist = pathSpline.dist_u(1.0)
		
			"""
			1) point that has the junction
			2) direction of junction
			3) set of long paths that include junction path  (only 2)
			"""

			""" initialize if this is the first time """
			particles = []
			for k in range(self.NUM_BRANCHES):
				#newArcDist = arcLow + (arcHigh-arcLow)*k/(self.NUM_BRANCHES-1)
				newArcDist = arcLow + k * self.DIV_LEN
				particle = self.getPrecomputedBranch(pathID, newArcDist)
				particles.append(particle)

			""" get the maximum value for each of our features """
			maxCost = -1e100
			maxMatchCount = -1000
			#maxDist = -1e100
			maxAngDiff = -1e100
			maxDist = 10.0


			""" find maximum values for each metric """
			for k in range(len(particles)):
				part = particles[k]
				matchCount = part["matchCount"]
				lastCost = part["lastCost"]
				dist = part["distDisc"]
				angDiff = part["angDisc"]

				if matchCount > maxMatchCount:
					maxMatchCount = matchCount

				#if matchCount > 0 and lastCost/matchCount > maxCost:
				#	maxCost = lastCost/matchCount

				if lastCost > maxCost:
					maxCost = lastCost

				if dist > maxDist:
					maxDist = dist

				if angDiff > maxAngDiff:
					maxAngDiff = angDiff


			""" invert the distance, angDiff and cost features """
			totalProbSum = 0.0
			for k in range(len(particles)):

				part = particles[k]

				lastCost = part["lastCost"]
				dist = part["distDisc"]
				angDiff = part["angDisc"]
				matchCount = part["matchCount"]

				#matchCost = 0.0
				#if maxCost > 0.0 and matchCount > 0:
				#	matchCost = (maxCost-lastCost/matchCount)/maxCost
				#else:
				#	matchCost = 0.0

							
				""" makes the maximum cost non-zero """
				matchCost = 0.1 + maxCost-lastCost
				#matchCost = maxCost-lastCost

				nomDist = 0.0
				if maxDist > 0.0:
					nomDist = (maxDist-dist)/maxDist
				else:
					nomDist = 0.0


				""" angDiff feature times matchCount times dist """
				#probVal = matchCount*matchCost + nomDist/4.0
				#probVal = matchCount*matchCost
				probVal = matchCount
				totalProbSum += probVal

				part2 = deepcopy(part)
				part2["lastCost"] = matchCost
				part2["distDisc"] = nomDist
				part2["initProb"] = probVal

				#part2 = (part[0], part[1], part[2], part[3], part[4], matchCost, nomDist, part[7], probVal, part[9], part[10], part[11])
				particles[k] = part2

			""" normalize the probability values """
			maxProb = 0.0
			for k in range(len(particles)):
				part = particles[k]

				probVal = 0.0

				if totalProbSum > 0.0:

					probVal = part["initProb"] / totalProbSum

					if probVal > maxProb:
						maxProb = probVal
				else:
					probVal = 0.0

				part2 = deepcopy(part)
				part2["initProb"] = probVal

				#part2 = (part[0], part[1], part[2], part[3], part[4], part[5], part[6], part[7], probVal, part[9], part[10], part[11])
				particles[k] = part2

				#part2 = (parentID, modJuncPose, newArcDist, pathID, matchCount1, lastCost1, distDisc, angDisc, initProb, newSplices, juncDiscAngle, juncDist)
				print "particle %04u %02u %02u %1.4f %03u %1.2f %1.2f %1.2f %1.2f %d %1.2f %1.2f" % (nodeID0, particleIndex, k, part2["modJuncPose"][2], part2["matchCount"], part2["lastCost"], part2["distDisc"], part2["angDisc"], part2["initProb"], len(part2["newSplices"]), part2["juncDiscAngle"], part2["juncDiscDist"])
				#print "particle %04u %02u %02u %1.4f %03u %1.2f %1.2f %1.2f %1.2f %d %1.2f %1.2f" % (nodeID0, particleIndex, k, part2[1][2], part2[4], part2[5], part2[6], part2[7], part2[8], len(part2[9]), part2[10], part2[11])
			print "particle"

			""" cartesian distance """
			DISC_THRESH = 0.5

			""" 60 degree threshold """
			#ANG_THRESH = 0.523 # pi/6
			ANG_THRESH = 1.5708 # pi/2

			
			""" get splices for new branch position """
			""" reject branch locations whose branch angle discrepancy is too high """ 
			for k in range(len(particles)):
				part = particles[k]

				#newSplices = deepcopy(splicedPaths)
				newProbVal = particles[k]["initProb"] 

				#part2["initProb"] = initProb
				#part2["newSplices"] = newSplices
				#part2["juncDiscAngle"] = juncDiscAngle
				#part2["juncDiscDist"] = juncDiscDist

				juncDiscAngle = part["juncDiscAngle"]

				if fabs(juncDiscAngle) > ANG_THRESH:
					newProbVal = 0.0

				part2 = deepcopy(part)
				part2["initProb"] = newProbVal

				#part2 = (part[0], part[1], part[2], part[3], part[4], part[5], part[6], part[7], newProbVal, part[9], part[10], part[11])
				particles[k] = part2



			#origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
			#origJuncPose[2] = 0.0
			#origJuncOrigin = Pose(origJuncPose)

			junctionDetails = self.leaf2LeafPathJunctions[pathID]
			localPathSegs = junctionDetails["localSegments"]
			#smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

			#print "path segs:", len(junctionDetails["leafSegments"]), "+", len(junctionDetails["internalSegments"]), "=", len(smoothPathSegs)

			#localPathSegs = []
			#for k in range(len(smoothPathSegs)):
			#	pathSeg = smoothPathSegs[k]
			#	localSeg = []
			#	for p in pathSeg:
			#		p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
			#		localSeg.append(p1)

			#	localPathSegs.append(localSeg)


			if True:
				pylab.clf() 
				for k,path in self.trimmedPaths.iteritems():
					print "path has", len(path), "points"
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])

					pylab.plot(xP,yP, color = self.colors[k], linewidth=4)

			if True:
				xP = []
				yP = []
				for part in particles:

					modControlPose = part["modControlPose"]
					modJuncPose = part["modJuncPose"]
					xP.append(modJuncPose[0])
					yP.append(modJuncPose[1])

					newSplices = part["newSplices"]

					for path in newSplices:
						xP1 = []
						yP1 = []
						for p in path:
							xP1.append(p[0])
							yP1.append(p[1])
						if maxProb > 0:
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part["initProb"]/maxProb)
						else:
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.10)


					"""
					offsetOrigin1 = Pose(modControlPose)

					for k in range(len(localPathSegs)):
						localSeg = localPathSegs[k]
						xP1 = []
						yP1 = []
						for p in localSeg:
							p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
							xP1.append(p1[0])
							yP1.append(p1[1])

						if maxProb > 0:
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part["initProb"]/maxProb)
						else:
							#pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part[8])
							pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.10)
					"""


				pylab.scatter(xP, yP, color='k', zorder=8)
				pylab.title("nodeID: %d hyp: %d, particleID: %d pathID: %d, localPathSegs %d" % (nodeID0, self.hypothesisID, particleIndex, pathID, len(localPathSegs)))
				
				pylab.savefig("bayes_plot_%04u_%02u_%03u_%04u.png" % (nodeID0, self.hypothesisID, particleIndex, self.tempCount) )
				self.tempCount += 1

		return particles

	@logFunction
	def getJuncSample(self, pathID, sampleIndex):
		if self.pathClasses[pathID]["parentID"] == None:
			raise

		path2 = self.paths[pathID]

		origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
		origJuncPose[2] = 0.0

		origJuncOrigin = Pose(origJuncPose)

		globJuncPose = branchPart[1]
		offsetOrigin1 = Pose(globJuncPose)
		
		path2 = self.paths[pathID]

		localPath2 = []
		for p in path2:
			p1 = origJuncOrigin.convertGlobalToLocal(p)
			localPath2.append(p1)

		samplePath2 = []
		for p in localPath2:
			p1 = offsetOrigin1.convertLocalToGlobal(p)
			samplePath2.append(p1)

		#localJunctionPoint = self.pathClasses[childPathID]["localJunctionPose"]
		#poseOrigin2 = Pose(self.nodeRawPoses[branchNodeID])
		#globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)

		#print "generated globJuncPose:",  globJuncPose

		return samplePath2
				
	""" generate the paths from the node and pose information """
	@logFunction
	def generatePaths(self):
		
		" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
		self.paths = {}
		self.hulls = {}
		self.leaf2LeafPathJunctions = {}
		
		pathIDs = self.getPathIDs()
		for pathID in pathIDs:
			print "computing path for node set", pathID, ":", self.getNodes(pathID)
			#self.paths[pathID], self.hulls[pathID] = self.computeShootSkeleton(pathID)
			results = computeShootSkeleton(self.poseData,
										pathID,
										self.pathClasses[pathID]["branchNodeID"],
										self.getGlobalJunctionPose(pathID),
										self.getControlPose(pathID),
										self.getNodes(pathID),
										self.nodePoses,
										self.hypothesisID,
										self.colors[pathID],
										self.topCount)

			self.leaf2LeafPathJunctions[pathID] = results[0]
			self.paths[pathID] = results[1]
			self.hulls[pathID] = results[2]

			self.topCount += 1
			
			#return theoryMedialLongPaths, medialLongPaths, leaf2LeafPathJunctions, medialLongPaths[maxIndex], vertices

		" compute the junction points between parent paths and child branches "
		for pathID in pathIDs:
			
			if self.pathClasses[pathID]["parentID"] != None:

				"""
				parentPathID = self.pathClasses[pathID]["parentID"]
				childPathID = pathID

				path1 = self.paths[parentPathID]
				path2 = self.paths[childPathID]
				branchNodeID = self.pathClasses[childPathID]["branchNodeID"]
				localJunctionPoint = self.pathClasses[childPathID]["localJunctionPose"]
				poseOrigin2 = Pose(self.nodeRawPoses[branchNodeID])
				globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)

				globJuncPose = getBranchPoint(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False)

				print "generated globJuncPose:",  globJuncPose
				
				self.pathClasses[childPathID]["globalJunctionPose"] = globJuncPose
				"""



		self.trimmedPaths = self.trimPaths(self.paths)

		localSkeletons = {}
		controlPoses = {}
		junctionPoses = {}
		for pathID in pathIDs:
			localSkeletons[pathID] = self.leaf2LeafPathJunctions[pathID]["skeletonGraph"]
			controlPoses[pathID] = self.pathClasses[pathID]["controlPose"]
			junctionPoses[pathID] = self.pathClasses[pathID]["globalJunctionPose"]

		parentPathIDs = self.getParentHash()
		
		finalPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)
		print "controlPoses:", controlPoses
		print "finalPoses:", finalPoses


		self.spliceSkeleton = spliceSkeletons(localSkeletons, controlPoses, junctionPoses, parentPathIDs)

		" for each path, attempt to join with its parent path "
		self.joins = []
		self.junctions = {}
		for pathID in pathIDs:

			cPath = self.getPath(pathID)
			
			parentPathID = cPath["parentID"]
			
			" parent does not concern us "
			if parentPathID == None:
				continue

			junctionPose = self.getGlobalJunctionPose(pathID)
			junctionPoint = [junctionPose[0],junctionPose[1]]

			print "node globJuncPose:",  junctionPose

			path1 = self.trimmedPaths[pathID]
			path2 = self.trimmedPaths[parentPathID]

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
			
			self.joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])

			" get junctions " 
			branchNodeID = cPath["branchNodeID"]
			
			self.junctions[pathID] = [branchNodeID, junctionPose, (parentPathID,minI2), path2[minI2], minI1]


		" create a tree with the node IDs and then stitch them together with joins "
		self.pathGraph = graph.graph()
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for k in range(len(path)):
				self.pathGraph.add_node((pathID, k), path[k])

			for k in range(len(path)-1):
				p1 = path[k]
				p2 = path[k+1]
				weight = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

				self.pathGraph.add_edge((pathID, k), (pathID, k+1), wt=weight)
			
			" parent does not concern us "
			cPath = self.getPath(pathID)			
			parentPathID = cPath["parentID"]
			
		" join with the junction in between the join points "
		for k in range(len(self.joins)):
			join = self.joins[k]

			pID1, k1 = join[0]
			pID2, k2 = join[1]

			p1 = self.trimmedPaths[pID1][k1]
			p2 = self.trimmedPaths[pID2][k2]
			weight = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)

			self.pathGraph.add_edge(join[0], join[1], wt=weight)



		self.topDict = {}
		
		
		" get terminals "

		self.terminals = {}
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			
			if len(path) > 0:
				if pathID == 0:
					
					dist1 = sqrt((self.trimmedPaths[pathID][0][0] - self.rootPoint[0])**2 + (self.trimmedPaths[pathID][0][1] - self.rootPoint[1])**2)
					dist2 = sqrt((self.trimmedPaths[pathID][len(path)-1][0] - self.rootPoint[0])**2 + (self.trimmedPaths[pathID][len(path)-1][1] - self.rootPoint[1])**2)

					if dist1 < dist2:
						self.terminals[0] = [(pathID,0), path[0]]
						self.terminals[1] = [(pathID,len(path)-1), path[len(path)-1]]
						self.topDict["t%u" % 0] = (pathID,0)
						self.topDict["t%u" % 1] = (pathID,len(path)-1)
						self.rootPoint = self.trimmedPaths[pathID][0]

					else:
						self.terminals[1] = [(pathID,0), path[0]]
						self.terminals[0] = [(pathID,len(path)-1), path[len(path)-1]]
						self.topDict["t%u" % 1] = (pathID,0)
						self.topDict["t%u" % 0] = (pathID,len(path)-1)
						self.rootPoint = self.trimmedPaths[pathID][len(path)-1]
					
				else:
					
					self.topDict["j%u" % pathID] = self.junctions[pathID][2]
					minI1 = self.junctions[pathID][4]
					
					" determine which side is junction and which side is terminal"
					if minI1 > len(path)-1 - minI1:    
						self.terminals[pathID+1] = [(pathID,0), path[0]]
						self.topDict["t%u" % (pathID+1)] = (pathID,0)
					else:
						self.terminals[pathID+1] = [(pathID,len(path)-1), path[len(path)-1]]
						self.topDict["t%u" % (pathID+1)] = (pathID,len(path)-1)


		#if len(self.trimmedPaths[pathID]) > 0:
		""" don't care about return values, are stored as object member variable """
		self.computeAllSplices(plotIter = True)
		self.isChanged = False

	@logFunction
	def getOrderedOverlappingPaths(self, nodeID):

		nodePose = self.getNodePose(nodeID)
		splicePaths, spliceTerms, splicePathIDs = self.getSplicesByNearJunction(nodePose)
		
		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.getNodePose(nodeID)


		resultSet = []

		for k in range(len(splicePaths)):
			
			print "comparing node", nodeID, "to splice", splicePathIDs[k], spliceTerms[k],  "plot", self.multiDepCount
			
			currPath = splicePaths[k]
			pathIDs = splicePathIDs[k]
			results = getMultiDeparturePoint(currPath, medial2, estPose2, estPose2, pathIDs, nodeID, pathPlotCount = self.multiDepCount, hypID = self.hypothesisID, plotIter = True)

			self.multiDepCount += 1

			resultSet.append(results+(k,))

			"departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departu rePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2"

		for result in resultSet:
			print "result:", result
		
		resultSet = sorted(resultSet, key=itemgetter(12), reverse=True)
		
		print "getOrderedOverlappingPaths() sorted node", nodeID, ":"
		for result in resultSet:
			print result[15]


		for k in range(len(resultSet)):
			result = resultSet[k]	 
			spliceIndex = result[15]
			overlappedPathIDs = splicePathIDs[spliceIndex]
			
			print "checking overlap of", overlappedPathIDs
			isNotOverlapped = False
			for pathID in overlappedPathIDs:

				#medial2 = self.poseData.medialAxes[nodeID]
				#estPose2 = self.nodePoses[nodeID]
				sum1 = getOverlapCondition(self.poseData.medialAxes[nodeID], self.getNodePose(nodeID), self.trimmedPaths[pathID], nodeID, plotIter = False, overlapPlotCount = self.overlapPlotCount)
				self.overlapPlotCount += 1

				#sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID, plotIter = False)

				print "pathID", pathID, "returned", sum1, "for splice", spliceIndex
				if sum1 > 1e10:
					isNotOverlapped = True
					
			if not isNotOverlapped:
				print "selected", k, "'th result"
				break
			
		
		orderedPathIDs = self.getPathOrdering(nodeID, overlappedPathIDs)

		print "returned ordered", orderedPathIDs, "from unordered", overlappedPathIDs
		
		return orderedPathIDs
		

	#@logFunction
	def getSplicesByNearJunction(self, initPose2):
		
		allSplices, terminals, junctions = self.getAllSplices(plotIter = False)

		#initPose2 = self.nodeHash[nodeID].getGlobalGPACPose()
		#initPose2 = self.nodePoses[nodeID]

		
		print "junctions:", junctions
		print "initPose2:", initPose2
		
		junctions2 = []
		for pathID, params in junctions.iteritems():
			junctionPose = params[1]
			
			print "checking", pathID, params
			
			dist2 = sqrt((junctionPose[0]-initPose2[0])**2 + (junctionPose[1]-initPose2[1])**2)
			
			print "dist2 =", dist2
			
			if dist2 < 3.0:
				junctions2.append((pathID,params[2][0]))
			
		"self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]"
		
		closePathID2 = self.getClosestPath(initPose2)

		print "closePathID2:", closePathID2
		
		pathSet2 = [closePathID2]

		junctionSet = junctions2
		junctionSet = list(set(junctionSet))

		print "junctions2:", junctions2

		for junc in junctionSet:
			pathSet2.append(junc[0])
			pathSet2.append(junc[1])

		print "pathSet2:", pathSet2

		pathSet2 = list(set(pathSet2))
		
		pathSet2.sort()

		print "pathSet2:", pathSet2
		
		termSet2 = []
		for pathID in pathSet2:
			if pathID == 0:
				termSet2.append(terminals[0][0])
				termSet2.append(terminals[1][0])
			else:
				termSet2.append(terminals[pathID+1][0])
				
		for pathID in pathSet2:
			parentID = self.getParentPathID(pathID)

			#if parentID != None and not parentID in pathSet2:
			if parentID != None:
				termSet2.append(junctions[pathID][2])
			

		

		splicePaths = []
		spliceTerms = []
		splicePathIDs = []

		print "termSet2:", termSet2
		print "allSplices:"

		for k, result in allSplices.iteritems():
			print "id:", k
			if k[0] in pathSet2 and k[1] in pathSet2:
				
				for sPath in result:
					termPath = sPath['termPath']
					print "termPath:", termPath
												
					if termPath[0] in termSet2 and termPath[-1] in termSet2:
						splicePaths.append(sPath['path'])
						spliceTerms.append(termPath)
						splicePathIDs.append(sPath['orderedPathIDs'])

		return splicePaths, spliceTerms, splicePathIDs




	@logFunction
	def comparePaths(self):
 
		allSplices, terminals, junctions = self.getAllSplices(plotIter = False)
		
		toBeMerged = []
		
		allPathIDs = self.getPathIDs()
		allPathIDs.sort()
		
		pathIDPairs = []
		for i in range(len(allPathIDs)):
			for j in range(i+1, len(allPathIDs)):
				pathIDPairs.append((allPathIDs[i],allPathIDs[j]))

		print "pathIDPairs:", pathIDPairs

		for pathPair in pathIDPairs:
			pathID1 = pathPair[0]
			pathID2 = pathPair[1]

			path1 = self.getPath(pathID1)
			path2 = self.getPath(pathID2)
			
			parentID1 = path1["parentID"]
			parentID2 = path2["parentID"]
			
			
			if parentID1 != None:
				pathParent1 = self.getPath(parentID1)
				grandParentID1 = pathParent1["parentID"]
			else:
				pathParent1 = None
				grandParentID1 = None

			if parentID2 != None:
				pathParent2 = self.getPath(parentID2)
				grandParentID2 = pathParent2["parentID"]
			else:
				pathParent2 = None
				grandParentID2 = None
			
			
			" Only consider case when two paths share the same parent "
			pathPairs = []
			
			rootPaths = []

			print "pathIDPairs:", pathID1, pathID2, parentID1, parentID2
			
			if False and grandParentID1 == pathID2 and parentID1 != None:
				
				if parentID2 != None:
					sPaths1 = allSplices[(parentID2, pathID1)]
					sPaths2 = allSplices[(parentID2, pathID2)]
				else:
					sPaths1 = allSplices[(pathID2, pathID1)]
					sPaths2 = allSplices[(pathID2, pathID2)]

				for sPath1 in sPaths1:
					termPath1 = sPath1['termPath']
					print "termPath1:", termPath1
					startKey = termPath1[0]
					endKey = termPath1[-1]

					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					path1 = []
					while currNode != endKey:
						path1.append(self.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path1.append(self.pathGraph.get_node_attributes(currNode))

					rootPaths.append(path1)

					for sPath2 in sPaths2:
						termPath2 = sPath2['termPath']
						
						print "termPath2:", termPath2


						startKey = termPath2[0]
						endKey = termPath2[-1]
	
						shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
						currNode = shortestPathSpanTree[startKey]					 
						path2 = []
						while currNode != endKey:
							path2.append(self.pathGraph.get_node_attributes(currNode))
							currNode = shortestPathSpanTree[currNode]
						path2.append(self.pathGraph.get_node_attributes(currNode))

						pathPairs.append((path1,path2,pathID1,pathID2))

			elif False and grandParentID2 == pathID1 and parentID2 != None:
				
				if parentID1 != None:
					sPaths1 = allSplices[(parentID1, pathID1)]
					sPaths2 = allSplices[(parentID1, pathID2)]
				else:
					sPaths1 = allSplices[(pathID1, pathID1)]
					sPaths2 = allSplices[(pathID1, pathID2)]

				for sPath1 in sPaths1:
					termPath1 = sPath1['termPath']
					print "termPath1:", termPath1
					startKey = termPath1[0]
					endKey = termPath1[-1]

					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					path1 = []
					while currNode != endKey:
						path1.append(self.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path1.append(self.pathGraph.get_node_attributes(currNode))

					rootPaths.append(path1)
					
					for sPath2 in sPaths2:
						termPath2 = sPath2['termPath']
						
						print "termPath2:", termPath2


						startKey = termPath2[0]
						endKey = termPath2[-1]
	
						shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
						currNode = shortestPathSpanTree[startKey]					 
						path2 = []
						while currNode != endKey:
							path2.append(self.pathGraph.get_node_attributes(currNode))
							currNode = shortestPathSpanTree[currNode]
						path2.append(self.pathGraph.get_node_attributes(currNode))

						pathPairs.append((path1,path2,pathID1,pathID2))

			elif parentID1 == pathID2:
				
				if parentID2 != None:
					sPaths1 = allSplices[(parentID2, pathID1)] + allSplices[(pathID1, pathID2)]
					sPaths2 = allSplices[(parentID2, pathID2)]
				else:
					sPaths1 = allSplices[(pathID2, pathID1)]
					sPaths2 = allSplices[(pathID2, pathID2)]

				for sPath1 in sPaths1:
					termPath1 = sPath1['termPath']
					print "termPath1:", termPath1

					startKey = termPath1[0]
					endKey = termPath1[-1] 
					if len(termPath1) > 2:		
						midKey = termPath1[1]
					else:
						midKey = None

					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					path1 = []
					midIndex = 0
					while currNode != endKey:

						if currNode == midKey:
							midIndex = len(path1)

						path1.append(self.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
						
						
					path1.append(self.pathGraph.get_node_attributes(currNode))
					
					print "path1:", len(path1), midIndex, midKey
					
					if midIndex > 150:
						path1 = path1[midIndex-150:]

					rootPaths.append(path1)
					
					for sPath2 in sPaths2:
						termPath2 = sPath2['termPath']
						
						print "termPath2:", termPath2


						startKey = termPath2[0]
						endKey = termPath2[-1]
						if len(termPath2) > 2:		
							midKey = termPath2[1]
						else:
							midKey = None
							
						shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
						currNode = shortestPathSpanTree[startKey]					 
						path2 = []
						midIndex = 0
						while currNode != endKey:
							
							if currNode == midKey:
								midIndex = len(path2)
															
							path2.append(self.pathGraph.get_node_attributes(currNode))
							currNode = shortestPathSpanTree[currNode]

						path2.append(self.pathGraph.get_node_attributes(currNode))

						print "path2:", len(path2), midIndex, midKey

						if midIndex > 150:
							path2 = path2[midIndex-150:]
							
						pathPairs.append((path2,path1,pathID2,pathID1))

			elif parentID2 == pathID1:
				
				if parentID1 != None:
					sPaths1 = allSplices[(parentID1, pathID1)] 
					sPaths2 = allSplices[(parentID1, pathID2)] + allSplices[(pathID1, pathID2)]
				else:
					sPaths1 = allSplices[(pathID1, pathID1)]
					sPaths2 = allSplices[(pathID1, pathID2)]

				for sPath1 in sPaths1:

					termPath1 = sPath1['termPath']
					print "termPath1:", termPath1

					startKey = termPath1[0]
					endKey = termPath1[-1]
					if len(termPath1) > 2:		
						midKey = termPath1[1]
					else:
						midKey = None

					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					path1 = []
					midIndex = 0
					
					while currNode != endKey:
						if currNode == midKey:
							midIndex = len(path1)

						path1.append(self.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path1.append(self.pathGraph.get_node_attributes(currNode))

					print "path1:", len(path1), midIndex, midKey

					if midIndex > 150:
						path1 = path1[midIndex-150:]
										
					rootPaths.append(path1)			   
					
					for sPath2 in sPaths2:
						termPath2 = sPath2['termPath']
						
						print "termPath2:", termPath2

						
						startKey = termPath2[0]
						endKey = termPath2[-1]
						if len(termPath2) > 2:		
							midKey = termPath2[1]
						else:
							midKey = None
	
						shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
						currNode = shortestPathSpanTree[startKey]					 
						path2 = []
						midIndex = 0
						while currNode != endKey:
							if currNode == midKey:
								midIndex = len(path2)
								
							path2.append(self.pathGraph.get_node_attributes(currNode))
							currNode = shortestPathSpanTree[currNode]
																						
						path2.append(self.pathGraph.get_node_attributes(currNode))

						print "path2:", len(path2), midIndex, midKey

						if midIndex > 150:
							path2 = path2[midIndex-150:]
							
						pathPairs.append((path1,path2,pathID1,pathID2))

			elif parentID1 == parentID2:
				
				values1 = allSplices[(parentID1,pathID1)]
				values2 = allSplices[(parentID2,pathID2)]
				
				print "sib:", pathID1, pathID2

				
				terms1 = []
				for sPath in values1:
					orderedIDs = sPath['orderedPathIDs']
					if len(orderedIDs) == 2:
						termPath1 = sPath['termPath']
						terms1.append(termPath1)
						print termPath1
				
				terms2 = []
				for sPath in values2:
					orderedIDs = sPath['orderedPathIDs']
					if len(orderedIDs) == 2:
						termPath2 = sPath['termPath']
						terms2.append(termPath2)
						print termPath2
				
				for k in range(len(terms1)):
					
					termPath1 = terms1[k]
					termPath2 = terms2[k]						 
					
					startKey = termPath1[0]
					endKey = termPath1[-1]

					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					path1 = []
					while currNode != endKey:
						path1.append(self.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path1.append(self.pathGraph.get_node_attributes(currNode))

					startKey = termPath2[0]
					endKey = termPath2[-1]

					shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
					currNode = shortestPathSpanTree[startKey]					 
					path2 = []
					while currNode != endKey:
						path2.append(self.pathGraph.get_node_attributes(currNode))
						currNode = shortestPathSpanTree[currNode]
					path2.append(self.pathGraph.get_node_attributes(currNode))

					rootPaths.append(path1)			   


				pathPairs.append((None,None,pathID1,pathID2))
				pathPairs.append((None,None,pathID2,pathID1))
			
			print
			
			for pair in pathPairs:
				pathID1 = pair[2]
				pathID2 = pair[3]

				if True:
					
					if pair[0] == None:
						
						print "comparing sibling paths", pathID1, pathID2
						resultPose0, resultPose1, lastCost0, matchCount0, juncAngDiff = self.makeSiblingPathCompare(pathID1, pathID2, parentID1, plotIter = False)
						
						print "sibling compare result:", resultPose0, resultPose1, lastCost0, matchCount0, juncAngDiff
						
						if fabs(resultPose0[2]) < 0.5 and fabs(resultPose1[2]) < 0.5 and fabs(juncAngDiff) < pi/4.0:

							print "queuing paths", pathID1, pathID2, "to be merged"

							toBeMerged.append((pathID1, pathID2, resultPose0, rootPaths, resultPose1, lastCost0, matchCount0))
					
					else:
						print "comparing paths", pathID1, pathID2
						print pair[1][0], pair[0][0]
						print pair[1][-1], pair[0][-1]
						resultPose0, lastCost0, matchCount0 = self.makePathCompare(pair[1], pair[0], pathID2, pathID1, plotIter = False)
					
						if fabs(resultPose0[2]) < 0.5:
		
							poseOrigin = Pose(resultPose0)
							path1 = pair[1]
							path1_offset = []
							for p in path1:
								result = poseOrigin.convertLocalToGlobal(p)
								path1_offset.append(result)					   
							
							result0 = self.getPathDeparture(pair[0], path1_offset, pathID1, pathID2, plotIter = False)
							print "getPathDeparture() = ", result0
		
							isInterior1 = result0[2]
							isInterior2 = result0[7]
		
							if not isInterior1 and not isInterior2:
			
								vals = allSplices[(pathID1,pathID1)]
								sPath = vals[0]
								termPath0 = sPath['termPath']
			
								startKey = termPath0[0]
								endKey = termPath0[-1]
				
								shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
								currNode = shortestPathSpanTree[startKey]
								path0 = []
								while currNode != endKey:
									path0.append(self.pathGraph.get_node_attributes(currNode))
									currNode = shortestPathSpanTree[currNode]
								path0.append(self.pathGraph.get_node_attributes(currNode))
								
								cost = self.getPathOverlapCondition(path0, path1_offset, pathID1, pathID2, plotIter = False)				   
								print "getPathOverlapCondition() =", cost
								
								if cost < 1e10:
									print pathID1, "and", pathID2, "are similar!"
									toBeMerged.append((pathID1, pathID2, resultPose0, rootPaths))
		
		return toBeMerged
					
	
	@logFunction
	def makeSiblingPathCompare(self, pathID1, pathID2, parentPathID, plotIter = False):


		globalPath1 = self.paths[pathID1]
		globalPath2 = self.paths[pathID2]
		parentPath = self.paths[parentPathID]

		junctionPose1 = self.pathClasses[pathID1]["globalJunctionPose"]
		junctionPose2 = self.pathClasses[pathID2]["globalJunctionPose"]
		
		print "sibling pathIDs:", pathID1, pathID2
		print "junctionPose1:", junctionPose1
		print "junctionPose2:", junctionPose2
		
		""" get closest point index onto the parent path to the child junctions """
		minDist_1 = 1e100
		minK_1 = 0
		minDist_2 = 1e100
		minK_2 = 0
		for k in range(len(parentPath)):
			p = parentPath[k]
		
			dist1 = sqrt((p[0]-junctionPose1[0])**2 + (p[1]-junctionPose1[1])**2)
			dist2 = sqrt((p[0]-junctionPose2[0])**2 + (p[1]-junctionPose2[1])**2)
		
			if dist1 < minDist_1:
				minDist_1 = dist1
				minK_1 = k

			if dist2 < minDist_2:
				minDist_2 = dist2
				minK_2 = k
		
		
		"""
		select which path segments to use by
		checking if the other child path's junction point is contained within

		excise portion of parent path between the two junction points 
		"""
		
		if minK_1 < minK_2:
			pathSeg1 = parentPath[:minK_1+1]
			pathSeg2 = parentPath[minK_2:]

		else:		 
			pathSeg2 = parentPath[:minK_2+1]
			pathSeg1 = parentPath[minK_1:]
		
		juncOrigin1 = Pose(junctionPose1)
		juncOrigin2 = Pose(junctionPose2)
		
		"""
		align the two junctions
		junction pose 2, with the orientation of junction pose 1		
		"""
		alignedJunctionPose2 = [junctionPose2[0], junctionPose2[1], junctionPose1[2]]

		""" aligned junction 2 with respect to junction 1 coordinate system """
		localJunctionPose2 = juncOrigin1.convertGlobalPoseToLocal(alignedJunctionPose2)
		
		print "alignedJunctionPose2:", alignedJunctionPose2
		print "localJunctionPose2:", localJunctionPose2
		
		""" pose of aligned junction 2 """
		alignJuncOrigin2 = Pose(alignedJunctionPose2)
		
		""" realigned segment 2 to be stitched onto segment 1"""
		newPathSeg2 = []
		for k in range(len(pathSeg2)):
			p = pathSeg2[k]
			p2 = alignJuncOrigin2.convertGlobalToLocal(p)			 
			p1 = juncOrigin1.convertLocalToGlobal(p2)
			newPathSeg2.append(p1)
		
		
		""" stitched path """
		if minK_1 < minK_2:
			stitchedPath = pathSeg1 + newPathSeg2
		else:
			stitchedPath = newPathSeg2 + pathSeg1
		
		
		if plotIter:
			
			pylab.clf()
			
			xP = []
			yP = []
			for p in stitchedPath:
				xP.append(p[0])
				yP.append(p[1])
				
			pylab.plot(xP,yP, color='k', zorder=10)
			
			xP = [junctionPose1[0], junctionPose2[0]]
			yP = [junctionPose1[1], junctionPose2[1]]
		
			pylab.scatter(xP,yP, color='k', zorder=10)
		
			xP = []
			yP = []
			for p in globalPath1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='r')

			xP = []
			yP = []
			for p in globalPath2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='b')


			xP = []
			yP = []
			for p in parentPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='g')
		
		
			pylab.title("parent = %d, siblings = %d %d" % (parentPathID, pathID1, pathID2))
			pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
		
			self.pathPlotCount2 += 1
		
		" compute the medial axis for each pose "
		globalPath1 = stitchedPath
		globalPath2 = parentPath


		
		globalSpline1 = SplineFit(globalPath1, smooth=0.1)
		globalSpline2 = SplineFit(globalPath2, smooth=0.1)
			
		orientedGlobalPath = orientPath(globalPath2, globalPath1)

		globalSpline1 = SplineFit(globalPath1, smooth=0.1)
		orientedGlobalSpline2 = SplineFit(orientedGlobalPath, smooth=0.1)


		globalSamples2 = orientedGlobalSpline2.getUniformSamples(spacing = 0.04)
		globalSamples1 = globalSpline1.getUniformSamples(spacing = 0.04)

		
		
		" compute the local variance of the angle "
		globalVar = computePathAngleVariance(globalSamples2)
		medialVar = computePathAngleVariance(globalSamples1)


		" now lets find closest points and save their local variances "			   
		closestPairs = []
		TERM_DIST = 20
		
		TERM_DIST1 = len(globalSamples1)/4
		TERM_DIST2 = len(globalSamples2)/4
		TERM_DIST1 = 0
		TERM_DIST2 = 0
		
		for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
			pG = globalSamples2[i]
			minDist = 1e100
			minJ = -1
			for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
				pM = globalSamples1[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minJ = j
			
			juncDist = math.sqrt((junctionPose1[0]-pM[0])**2 + (junctionPose1[1]-pM[1])**2)
			#if minDist < 0.1:
			if juncDist < 1.0:
				closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1], juncDist))

		for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
			pM = globalSamples1[j]
			minDist = 1e100
			minI = -1
			for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
				pG = globalSamples2[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minI = i
					
			juncDist = math.sqrt((junctionPose1[0]-pM[0])**2 + (junctionPose1[1]-pM[1])**2)
			#if minDist < 0.1:
			if juncDist < 1.0:
				closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1],juncDist))

		""" remove duplicates """
		closestPairs = list(set(closestPairs))
		
		""" sort by minDist, then lowest angular variance, then juncDist """
		closestPairs = sorted(closestPairs, key=itemgetter(2,5,6,7))
		
		print len(closestPairs), "closest pairs"

		if plotIter:
			
			pylab.clf()
			
			xP = [junctionPose1[0]]
			yP = [junctionPose1[1]]
		
			pylab.scatter(xP,yP, color='k', zorder=10)

			xP = []
			yP = []
			for p in globalPath1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='k')

			xP = []
			yP = []
			for p in globalPath2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='g')

			for pair in closestPairs:
				p1 = globalSamples1[pair[1]]
				p2 = globalSamples2[pair[0]]		
				pylab.plot([p1[0],p2[0]], [p1[1],p2[1]])

			pylab.title("parent = %d, siblings = %d %d" % (parentPathID, pathID1, pathID2))
			pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
		
			self.pathPlotCount2 += 1




		pathDict1 = self.getPath(pathID1)
		pathDict2 = self.getPath(pathID2)
		
		parentID1 = pathDict1["parentID"]
		parentID2 = pathDict2["parentID"]
		
		if len(closestPairs) > 0:

			globalJunctionPose1 = pathDict1["globalJunctionPose"]
			globalJunctionPose2 = pathDict2["globalJunctionPose"]			 

			originU2 = globalSpline1.findU(globalJunctionPose1)    
			originU1 = orientedGlobalSpline2.findU(globalJunctionPose1)
			isFound = False
			
			for pair in closestPairs:
				juncDist = pair[7]
				
				""" always True because we already reject juncDist > 1.0 """
				if juncDist < 1.0:
			
					originU2 = globalSpline1.findU(globalSamples1[pair[1]])    
					originU1 = orientedGlobalSpline2.findU(globalSamples2[pair[0]])
					isFound = True
					break
				
			if not isFound:
				""" theoretically, this should never be executed """
				for pair in closestPairs:
					print pair
				raise
		else:
			""" no closest pairs were found with our criteria """
			raise
		
		u2 = originU2
		u1 = originU1
		angGuess = 0.0


		resultPose1, lastCost1, matchCount1 = gen_icp.pathOverlapICP([u1,u2,angGuess], orientedGlobalPath, globalPath1, plotIter = plotIter, n1 = pathID1, n2 = pathID2)

		print "matchCount,cost,result:", matchCount1, lastCost1, resultPose1

		globalPath1 = self.paths[pathID1]
		globalPath2 = self.paths[pathID2]
		parentPath = self.paths[parentPathID]

		if plotIter:
			
			pylab.clf()
			
			xP = [junctionPose1[0], junctionPose2[0]]
			yP = [junctionPose1[1], junctionPose2[1]]
		
			pylab.scatter(xP,yP, color='k', zorder=10)

			poseOrigin = Pose(resultPose1)
			
			xP = []
			yP = []
			for p in stitchedPath:
				result = poseOrigin.convertLocalToGlobal(p)
				xP.append(result[0])
				yP.append(result[1])

			pylab.plot(xP,yP, color='k', zorder=10)
			   
			xP = []
			yP = []
			for p in globalPath1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='r')

			xP = []
			yP = []
			for p in globalPath2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='b')


			xP = []
			yP = []
			for p in parentPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='g')
		
		
			pylab.title("parent = %d, siblings = %d %d" % (parentPathID, pathID1, pathID2))
			pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
		
			self.pathPlotCount2 += 1
		

		
		juncOrigin1 = Pose(junctionPose1)
		juncOrigin2 = Pose(junctionPose2)
		
		print "junctionPose1:", junctionPose1
		print "junctionPose2:", junctionPose2

		print "resultPose1:", resultPose1
		
		
		" align the two junctions "
		alignedJunctionPose2 = [junctionPose2[0], junctionPose2[1], junctionPose1[2]]
		
		print "alignedJunctionPose2:", alignedJunctionPose2
		
		offsetOrigin1 = Pose(resultPose1)
		newJuncPose1 = offsetOrigin1.convertLocalOffsetToGlobal(junctionPose1)

		print "newJuncPose1:", newJuncPose1

		newJuncOrigin1 = Pose([newJuncPose1[0], newJuncPose1[1], 0.0])
		resultOffset = newJuncOrigin1.convertGlobalPoseToLocal([junctionPose2[0],junctionPose2[1], 0.0])
		
		print "resultOffset:", resultOffset
		
		resultPose2 = newJuncOrigin1.doInverse(resultOffset)



		def computeOffset(pose1, pose2):
		
			" corner points and orientations "
			corner1Pose = pose1
			corner2Pose = pose2
			
			" convert the desired intersection point on curve 1 into global coordinates "
			poseProfile1 = Pose([0.0,0.0,0.0])
				
			" now convert this point into a pose, and perform the inverse transform using corner2Pose "
			desGlobalPose2 = Pose(corner1Pose)
			
			" perform inverse offset from the destination pose "
			negCurve2Pose = desGlobalPose2.doInverse(corner2Pose)
			
			" relative pose between pose 1 and pose 2 to make corners coincide and same angle "
			resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
			localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
			
			return [localOffset[0], localOffset[1], localOffset[2]]
		
		
		resultPose2 = computeOffset(newJuncPose1, alignedJunctionPose2)

		print "resultPose2:", resultPose2

		globalPath1 = self.paths[pathID1]
		globalPath2 = self.paths[pathID2]
		parentPath = self.paths[parentPathID]

		offsetOrigin1 = Pose(resultPose1)
		offsetOrigin2 = Pose(resultPose2)

		offJuncPose1 = offsetOrigin1.convertLocalOffsetToGlobal(junctionPose1)
		offJuncPose2 = offsetOrigin2.convertLocalOffsetToGlobal(junctionPose2)

		juncAngDiff = diffAngle(offJuncPose1[2],offJuncPose2[2])

		if plotIter:
			
			pylab.clf()
			
			xP = [junctionPose1[0], junctionPose2[0]]
			yP = [junctionPose1[1], junctionPose2[1]]
		
			pylab.scatter(xP,yP, color='k', zorder=9)
		
			xP = []
			yP = []

			for p in globalPath1:
				p1 = offsetOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])
			pylab.plot(xP,yP,color='r', zorder=8)

			offJuncPoint1 = offsetOrigin1.convertLocalToGlobal([junctionPose1[0],junctionPose1[1]])
			pylab.scatter([offJuncPoint1[0]],[offJuncPoint1[1]], color='r', zorder=10)
			print "offJuncPoint1:", offJuncPoint1

			xP = []
			yP = []
			for p in globalPath2:
				p1 = offsetOrigin2.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])
			pylab.plot(xP,yP,color='b', zorder = 8)
			
			offJuncPoint2 = offsetOrigin2.convertLocalToGlobal([junctionPose2[0],junctionPose2[1]])
			pylab.scatter([offJuncPoint2[0]],[offJuncPoint2[1]], color='b', zorder=10)
			
			
			xP = []
			yP = []
			for p in parentPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP,color='g')
		
			pylab.title("resultPose = [%1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f] %1.2f\n" % (resultPose1[0], resultPose1[1], resultPose1[2], resultPose2[0], resultPose2[1], resultPose2[2], juncAngDiff))
		
			pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
		
			self.pathPlotCount2 += 1
		

		


		return resultPose1, resultPose2, lastCost1, matchCount1, juncAngDiff
	 



	@logFunction
	def makePathCompare(self, globalPath1, globalPath2, pathID1, pathID2, plotIter = False):

		" compute the medial axis for each pose "
		
		globalPath1Reverse = deepcopy(globalPath1)
		globalPath1Reverse.reverse()
		
		globalPath2Reverse = deepcopy(globalPath2)
		globalPath2Reverse.reverse()

		globalSpline1 = SplineFit(globalPath1, smooth=0.1)
		globalSpline1Reverse = SplineFit(globalPath1Reverse, smooth=0.1)
		globalSpline2 = SplineFit(globalPath2, smooth=0.1)
		globalSpline2Reverse = SplineFit(globalPath2Reverse, smooth=0.1)


		""" ORIENT PATH """
		orientedGlobalPath = orientPath(globalPath2, globalPath1)
			
		orientedGlobalSpline2 = SplineFit(orientedGlobalPath, smooth=0.1)


		globalSamples2 = orientedGlobalSpline2.getUniformSamples(spacing = 0.04)
		globalSamples1 = globalSpline1.getUniformSamples(spacing = 0.04)

		
		""" COMPUTE ANGLE VARIANCES """
		globalVar = computePathAngleVariance(globalSamples2)
		medialVar = computePathAngleVariance(globalSamples1)


		""" now lets find closest points and save their local variances """			   
		closestPairs = []
		TERM_DIST = 20
		
		TERM_DIST1 = len(globalSamples1)/4
		TERM_DIST2 = len(globalSamples2)/4
		
		
		for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
			pG = globalSamples2[i]
			minDist = 1e100
			minJ = -1
			for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
				pM = globalSamples1[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minJ = j
					
			#if minDist < 0.1:
			if True:
				closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

		for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
			pM = globalSamples1[j]
			minDist = 1e100
			minI = -1
			for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
				pG = globalSamples2[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minI = i
					
			#if minDist < 0.1:
			if True:
				closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))

		" remove duplicates "
		closestPairs = list(set(closestPairs))
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(2,5,6))

		pathDict1 = self.getPath(pathID1)
		pathDict2 = self.getPath(pathID2)
		
		parentID1 = pathDict1["parentID"]
		parentID2 = pathDict2["parentID"]
		
		if parentID1 == parentID2:
			
			globalJunctionPose1 = pathDict1["globalJunctionPose"]
			globalJunctionPose2 = pathDict2["globalJunctionPose"]


			originU2 = globalSpline1.findU(globalJunctionPose1)    
			originU1 = orientedGlobalSpline2.findU(globalJunctionPose2)
		
		elif len(closestPairs) > 0:
			originU2 = globalSpline1.findU(globalSamples1[closestPairs[0][1]])	  
			originU1 = orientedGlobalSpline2.findU(globalSamples2[closestPairs[0][0]])

		else:
			raise
		
		u2 = originU2
		u1 = originU1
		angGuess = 0.0

		initGuess = [u1,u2,angGuess]
		resultPose1, lastCost1, matchCount1 = gen_icp.pathOverlapICP(initGuess, orientedGlobalPath, globalPath1, plotIter = plotIter, n1 = pathID1, n2 = pathID2)
		


		print "matchCount,cost,result:", matchCount1, lastCost1, resultPose1


		return resultPose1, lastCost1, matchCount1
	 
	 
	@logFunction
	def delNode(self, nodeID, pathID):
		print "deleting node", nodeID, "from path", pathID
		self.pathClasses[pathID]["nodeSet"].remove(nodeID)
		

	@logFunction
	def moveNode(self, nodeID, splicePathIDs):
		" change the shoot membership based on fitted splice "

		particleDist = self.poseParticles["snapshots2"][0]

		" determine which paths are leaves "
		pathIDs = self.getPathIDs()

		" remove node from other shoots first"
		for pathID in pathIDs:
			if nodeID in self.pathClasses[pathID]["nodeSet"]:
				print "removing node", nodeID, "from path", pathID
				self.pathClasses[pathID]["nodeSet"].remove(nodeID)

		isAParent = {}
		for k in pathIDs:
			isAParent[k] = False
		for k in splicePathIDs:
			currPath = self.getPath(k)
			currParent = currPath["parentID"]
			if currParent != None:
				isAParent[currParent] = True

		"add nodes to paths that are the leaves "
		for pathID in splicePathIDs:
			if not isAParent[pathID]:
				print "adding node", nodeID, "to path", pathID
				self.pathClasses[pathID]["nodeSet"].append(nodeID)
				#self.addNode(nodeID,pathID)

				for p in particleDist:
					p.addNode(nodeID, pathID)

		print "shoot", pathID, "has nodes", self.pathClasses[pathID]["nodeSet"]

		self.isChanged = True		

	@logFunction
	def addNode(self, nodeID, pathID):


		print "adding node", nodeID, "to path", pathID
		if not nodeID in self.pathClasses[pathID]["nodeSet"]:
			self.pathClasses[pathID]["nodeSet"].append(nodeID)
		else:
			print "node", nodeID, "already present in path", pathID

		updateCount = self.poseParticles["updateCount"] 
		particleDist = self.poseParticles["snapshots2"][0]

		for p in particleDist:
			p.addNode(nodeID, pathID)

		print pathID, "has nodes", self.pathClasses[pathID]["nodeSet"]

		self.isChanged = True
		#self.junctions = {}
		#self.terminals = {}
		#self.allSplices = {}

	@logFunction
	def delPath(self, pathID, mergeTargetID):
		
		print "deleting path", pathID, "merging to path", mergeTargetID
		try: 
			del self.pathClasses[pathID]
			del self.pathTermsVisited[pathID]
			
			
			for k, pathClass in self.pathClasses.iteritems():
				
				if pathClass["parentID"] == pathID:
					pathClass["parentID"] = mergeTargetID
					
		except:
			pass
	
	@logFunction
	def addPath(self, parentID, branchNodeID, localJunctionPose):

		"""

		Create a new shoot branch given the following starting information:
		
		1) the parent shoot ID
		2) the spatial node ID that begins the new shoot with a detected divergence
		3) the pose of the divergence point in the local coordinates of the spatial node


		From this compute the following information:

		1) the junction point, representing the center of a hypothetical junction
		2) a control point on the parent shoot that allows a distribution of possible branch locations
		3) the computed shoot branch parameters for each currently active pose particle

		"""


		print "addPath(", parentID, branchNodeID, localJunctionPose
		
		""" the raw pose and posture pose """
		estPose = self.nodePoses[branchNodeID]
		nodePose = self.nodeRawPoses[branchNodeID]

		""" the localJunctionPose computed from raw local coordinates to global coordinates """
		nodeOrigin = Pose(nodePose)
		globalJunctionPose = nodeOrigin.convertLocalOffsetToGlobal(localJunctionPose)
	
		
		allPathIDs = self.getPathIDs()

		""" new shoot ID """
		newPathID = self.pathIDs

		""" posture curve of local spatial map """
		medial0 = self.poseData.medialAxes[branchNodeID]

		""" medial axis converted to global coordinates from local posture coordinates """
		poseOrigin0 = Pose(estPose)
		globalMedial0 = []
		for p in medial0:
			globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))


		""" compute branch point of posture curve diverging from parent shoot """
		newGlobJuncPose1, controlPoint1, angDeriv1 = getBranchPoint(globalJunctionPose, parentID, newPathID, self.trimmedPaths[parentID], globalMedial0, plotIter = True, hypothesisID = self.hypothesisID, nodeID = branchNodeID)

		""" compute branch point of parent shoot diverging from posture curve """
		newGlobJuncPose2, controlPoint2, angDeriv2 = getBranchPoint(globalJunctionPose, newPathID, parentID, globalMedial0, self.trimmedPaths[parentID], plotIter = True, hypothesisID = self.hypothesisID, nodeID = branchNodeID)

		""" the control point that is the flattest determines which curve is the diverging one """
		""" the curve that is sharpest at the control point is the diverging curve """
		if angDeriv1 > angDeriv2:
			newGlobJuncPose = newGlobJuncPose1
			controlPoint = controlPoint1
			angDeriv = angDeriv1
		else:
			newGlobJuncPose = newGlobJuncPose2
			controlPoint = controlPoint2
			angDeriv = angDeriv2


		commonU1, commonU2, cPoint1, cPoint2 = selectCommonOrigin(self.trimmedPaths[parentID],globalMedial0)	

		#controlPose = [controlPoint[0], controlPoint[1], 0.0]
		controlPose = [cPoint1[0],cPoint1[1], 0.0]

		""" convert the controlPose to coordinates local to the parent frame """
		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)
		
		""" compute the control pose relative to the parent """
		parentControlPose = globalControlPoses[parentID]
		parentFrame = Pose(parentControlPose)
		localControlPose = parentFrame.convertGlobalPoseToLocal(controlPose)


		""" basic shoot data structure """
		self.pathClasses[newPathID] = {"parentID" : parentID,
						"branchNodeID" : branchNodeID,
						"localJunctionPose" : localJunctionPose, 
						"sameProb" : {},
						"nodeSet" : [branchNodeID,],
						"globalJunctionPose" : newGlobJuncPose,
						"controlPose" : localControlPose }		
						#"controlPose" : controlPose }		

		print "new pathClass:", self.pathClasses[newPathID]


		self.pathTermsVisited[newPathID] = False
		
		self.pathIDs += 1


		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]

		""" FIXME:  Full path spline instead of just relevant section.  Does this always work? """
		pathSpline = SplineFit(self.paths[parentID])


		""" add new path to particles """
		for part in particleDist2:

			origPose = self.getNodePose(branchNodeID)

			""" get local transform of GPAC to raw pose """
			gpacProfile = Pose(origPose)
			localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[branchNodeID])
			
			""" go back and convert this from GPAC pose to estPose """
			particlePose = deepcopy(part.pose0)
			particlePose[2] = origPose[2]
			newProfile = Pose(particlePose)
			newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)

			localControlPose = gpacProfile.convertGlobalPoseToLocal(controlPose)
			particleControlPose = newProfile.convertLocalOffsetToGlobal(localControlPose)

			#print "control poses:", origPose, particlePose, controlPose, particleControlPose
		
			#nodePose = part.pose0

			""" location of junction point from raw to global in the particle pose """
			nodeOrigin = Pose(newEstPose)
			newGlobalJunctionPose = nodeOrigin.convertLocalOffsetToGlobal(localJunctionPose)

			""" location on the parent shoot closest to the locally computed branch point """
			minDist, uVal, newJuncPose = pathSpline.findClosestPoint(newGlobalJunctionPose)
			#arcDist = pathSpline.dist_u(uVal)

			""" angle of junction corrected to the original orientation """
			modJuncPose = copy(newJuncPose)
			modJuncPose[2] = globalJunctionPose[2]

			#""" get the arc distance of the given branch point """
			#arcDist = pathSpline.dist_u(uVal)
			""" get the arc distance of the control point """
			minDist, controlUVal, newControlPose = pathSpline.findClosestPoint(particleControlPose)
			arcDist = pathSpline.dist_u(controlUVal)

			""" compute the rails of the branch point distribution """
			arcHigh = arcDist + self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)
			arcLow = arcDist - self.DIV_LEN * floor(self.NUM_BRANCHES/2.0)

			""" compute the sample points for the branch point distribution """
			arcDists = []
			for k in range(self.NUM_BRANCHES):
				newArcDist = arcLow + k * self.DIV_LEN
				arcDists.append((newPathID, newArcDist))

			""" convert the controlPose to coordinates local to the parent frame """
			""" compute the control pose relative to the parent """
			parentPathIDs = self.getParentHash()
			controlPoses = self.getControlPoses()
			globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)
			parentControlPose = globalControlPoses[parentID]
			parentFrame = Pose(parentControlPose)

			controlPoses = []
			for k in range(self.NUM_BRANCHES):
				globalControlPose = pathSpline.getPointOfDist(arcDists[k][1])
				localControlPose = parentFrame.convertGlobalPoseToLocal(globalControlPose)
				controlPoses.append(localControlPose)

			localParticleControlPose = parentFrame.convertGlobalPoseToLocal(particleControlPose)

			#controlPoses = [pathSpline.getPointOfDist(arcDists[k][1]) for k in range(self.NUM_BRANCHES)]

			""" add the details of this junction given particle pose is true """
			part.addPath(newPathID, parentID, branchNodeID, localJunctionPose, modJuncPose, localParticleControlPose, self.NUM_BRANCHES, arcDists, controlPoses)

		print "newPath", newPathID, "=", self.pathClasses[newPathID]

		self.isChanged = True

		return newPathID
 
	@logFunction
	def trimPaths(self, foo = None):

		paths = self.paths
		trimmedPaths = {}
		
		print "path lengths:"
		for k,v in paths.iteritems():
			print k, len(v)
		
		
		if len(paths) <= 1:
			for k, v in paths.iteritems():
				trimmedPaths[k] = v
				
			self.trimmedPaths = trimmedPaths
			return trimmedPaths

		" path0 is root, and path1 is a child of path0 "

		" information on parentage is required to determine which path belongs to which "

		" global junction point for parent 0 and child 1 "


		pathIDs = self.getPathIDs()

		childParents = {}
		isParentComputed = {}
		for pathID in pathIDs:
			isParentComputed[pathID] = False
			parentPathID = self.pathClasses[pathID]["parentID"]
			childParents[pathID] = parentPathID

		isParentComputed[0] = True

		parentPathIDs = self.getParentHash()
		controlPoses = self.getControlPoses()
		globalControlPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		while False in isParentComputed.values():
			print "trim parents:", isParentComputed, childParents

			for pathID in pathIDs:
		
				if isParentComputed[pathID]:
				
					if childParents[pathID] == None:
						trimmedPaths[pathID] = deepcopy(paths[pathID])

					else:
						parentPathID = childParents[pathID]
						childPathID = pathID

						path1 = paths[parentPathID]
						path2 = paths[childPathID]

								
						origControlPose = self.pathClasses[childPathID]["controlPose"]
						origGlobalControlPose = self.getFinalControlPose(childPathID)


						branchNodeID = self.pathClasses[childPathID]["branchNodeID"]
						branchNodePose = self.getNodePose(branchNodeID)
						globalJunctionPose = self.getGlobalJunctionPose(childPathID)

						modJuncPose = [globalJunctionPose[0], globalJunctionPose[1], 0.0]

						globalControlPose = globalControlPoses[pathID]

						newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths = trimBranch(childPathID, parentPathID, origGlobalControlPose, globalControlPose, globalJunctionPose, modJuncPose, path2, path1, trimmedPaths[parentPathID], [], plotIter=False, hypothesisID=self.hypothesisID, nodeID=(self.poseData.numNodes-1))
						#newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths = trimBranch(childPathID, parentPathID, origControlPose, modJuncPose, globalJunctionPoint, path2, path1, trimmedPaths[parentPathID], [], plotIter=False, hypothesisID=self.hypothesisID, nodeID=(self.poseData.numNodes-1))

						trimmedPaths[pathID] = deepcopy(newPath3)

					for childID in pathIDs:
						if self.pathClasses[childID]["parentID"] == pathID:
							isParentComputed[childID] = True

		self.trimmedPaths = trimmedPaths
			
		return trimmedPaths


	@logFunction
	def pathTermVisited(self, pathID):
		
		self.pathTermsVisited[pathID] = True

	@logFunction
	def getPathTermsVisited(self):
		return self.pathTermsVisited

	@logFunction
	def resetTerms(self):
		
		for k, v in self.pathTermsVisited.iteritems():
			self.pathTermsVisited[k] = False

	@logFunction
	def getPathTerms(self):

		terms = {}

		pathIDs = self.getPathIDs()

		for pathID in pathIDs:

			pID, pathI = self.topDict["t%u" % (pathID+1)]
			path = self.trimmedPaths[pID]
			terms[pathID] = [path[pathI][0],path[pathI][1], 0.0]
			
		print "returning path terms:", terms

		return terms

	@logFunction
	def getClosestPath(self, locPose):
	
		pathIDs = self.getPathIDs()
		
		minDist = 1e100
		minPathID = 0
		for pathID in pathIDs:
			
			pathSeg = self.trimmedPaths[pathID]
			
			for p in pathSeg:
				dist = sqrt((p[0]-locPose[0])**2+(p[1]-locPose[1])**2)
				
				if dist < minDist:
					minDist = dist
					minPathID = pathID
		
		return minPathID

	@logFunction
	def getPathPath(self, startPathID, endPathID):
		" given start pathID and end pathID, what is the series of paths I need to follow to get from start to end "

		if startPathID == endPathID:
			return [startPathID]
		
		pathGraph = graph.graph()

		pathIDs = self.getPathIDs()

		for pathID in pathIDs:
			pathGraph.add_node(pathID, [])
		
		for pathID in pathIDs:	  
			parentPathID = self.getPath(pathID)["parentID"]
			if parentPathID != None:
				pathGraph.add_edge(pathID, parentPathID)
		
		shortestPathSpanTree, shortestDist = pathGraph.shortest_path(endPathID)
		nextPathID = shortestPathSpanTree[startPathID]

		splicedPaths = [startPathID, nextPathID]
		while nextPathID != endPathID:
			nextPathID = shortestPathSpanTree[nextPathID]
			splicedPaths.append(nextPathID)

		return splicedPaths

	@logFunction
	def getPathOrdering(self, nodeID, pathIDs):

		plotIter = False

		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.getNodePose(nodeID)
			
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
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			 
		
		pathSelect = []
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

				
				except:
					pass

			if minPathID != -1:
				pathSelect.append(minPathID)
		
		print "pathSelect:", pathSelect
		
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

	@logFunction
	def getPathDeparture(self, path1, path2, pathID1, pathID2, termPathIDs1 = [], termPathIDs2 = [], plotIter = False):
		
		isExist1 = False
		isInterior1 = False
		departurePoint1 = 0
		angle1 = 0.0
		isExist2 = False
		isInterior2 = False
		departurePoint2 = 0
		angle2 = 0.0
		
		if len(path1) == 0:
			return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
				
				
				
				
		orientedPath1 = orientPath(path1, path2)		   


		"""Assumption:  one section of the medial axis is closely aligned with the path """
		
		pathSpline2 = SplineFit(path2, smooth=0.1)
		pathPoints2 = pathSpline2.getUniformSamples()
		print "path2:", path2[0], pathPoints2[0], path2[-1], pathPoints2[-1]

		pathSpline1 = SplineFit(orientedPath1, smooth=0.1)
		pathPoints1 = pathSpline1.getUniformSamples()
		print "path1:", orientedPath1[0], pathPoints1[0], orientedPath1[-1], pathPoints1[-1]


		angle1, angle2 = getTipAngles(pathPoints2)

		distSum = 0.0
		contigCount = 0
		maxContig = 0
		distances = []
		indices = []
		for i in range(0,len(pathPoints2)):
			p_2 = pathPoints2[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
			distances.append(minDist)
			indices.append(i_1)
			distSum += minDist
			
			if minDist < 0.2:
				contigCount += 1
				if contigCount > maxContig:
					maxContig = contigCount
			else:
				contigCount = 0
		
		overlapSum = distSum / float(len(pathPoints2))
		
		print "maxContig,overlapSum:", maxContig, overlapSum

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


		frontAngleRefI = frontDepI
		while distances[frontAngleRefI] < 0.1:
			frontAngleRefI -= 1
			
			if frontAngleRefI < 0:
				frontAngleRefI = 0
				break
						
		forePathIndex = indices[frontAngleRefI]
		forePathAngle = pathPoints1[forePathIndex][2]
		forePathAngle = normalizeAngle(forePathAngle + pi)

		foreDiffAngle = diffAngle(angle1, forePathAngle)
		
		newFrontDepI = 4
		while fabs(diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]+pi), forePathAngle)) > pi/6.0:
			
			if newFrontDepI < frontDepI:
				newFrontDepI += 1
			else:
				break
		
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

		backAngleRefI = backDepI
		while distances[backAngleRefI] < 0.1:
			backAngleRefI += 1
			
			if backAngleRefI >= len(distances):
				backAngleRefI = len(distances)-1
				break
		
		backPathIndex = indices[backAngleRefI]
		backPathAngle = pathPoints1[backPathIndex][2]
		
		backDiffAngle = diffAngle(angle2, backPathAngle)
		
		newBackDepI = len(distances)-5
		while fabs(diffAngle(pathPoints2[newBackDepI][2], backPathAngle)) > pi/6.0:
			
			if newBackDepI > backDepI:
				newBackDepI -= 1
			else:
				break		 


		foreAngDiffs = []
		for i in range(0,len(pathPoints2)):
			foreAngDiffs.append(fabs(diffAngle(normalizeAngle(pathPoints2[i][2] + pi), forePathAngle)))

		backAngDiffs = []
		for i in range(0,len(pathPoints2)):
			backAngDiffs.append(fabs(diffAngle(pathPoints2[i][2], backPathAngle)))
			
		"reset to the maximum distances "
		maxFront = distances[0]
		maxBack = distances[-1]
		
		
		
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



		arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
		arr2 = numpy.array(deepcopy(indices[backDepI:]))
		matchVar1 = arr1.var()
		matchVar2 = arr2.var()

		tipMatch1 = pathPoints1[indices[0]]

		" discrepancy distance between tip closest point and average departure point "
		dist1 = sqrt((tipMatch1[0]-pathPoints1[max1][0])**2 + (tipMatch1[1] - pathPoints1[max1][1])**2)

		tipMatch2 = pathPoints1[indices[-1]]

		" discrepancy distance between tip closest point and average departure point "
		dist2 = sqrt((tipMatch2[0]-pathPoints1[max2][0])**2 + (tipMatch2[1] - pathPoints1[max2][1])**2)


		DEP_THRESH = 0.2

		if maxFront > DEP_THRESH:
			departurePoint1 = pathPoints2[newFrontDepI]
			isExist1 = True

			if max1 == 0 or max1 == len(pathPoints1)-1:
				isInterior1 = False
			else:
				isInterior1 = True

		if maxBack > DEP_THRESH:
			departurePoint2 = pathPoints2[newBackDepI]
			isExist2 = True

			if max2 == 0 or max2 == len(pathPoints1)-1:
				isInterior2 = False
			else:
				isInterior2 = True

		" sum of closest points on front and back "
		" select the one with minimal cost "
		
		if False:
			pylab.clf()
			xP = range(len(pathPoints2))
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
			pylab.title("%1.2f %1.2f %d %d %d" % ( maxFront, maxBack, len(pathPoints1), max1, max2))
			pylab.savefig("distances_%04u.png" % self.pathPlotCount)

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in pathPoints2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')   

			if True:	
				xP = [pathPoints1[max1][0]]
				yP = [pathPoints1[max1][1]]
				pylab.scatter(xP,yP, color='b')		   
	
				xP = [tipMatch1[0]]
				yP = [tipMatch1[1]]
				pylab.scatter(xP,yP, color='r')		   

				xP = [pathPoints2[frontDepI][0]]
				yP = [pathPoints2[frontDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints1[forePathIndex][0]]
				yP = [pathPoints1[forePathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [pathPoints2[newFrontDepI][0]]
				yP = [pathPoints2[newFrontDepI][1]]
				pylab.scatter(xP,yP, color='m')		   
	
	
			if True:
				xP = [pathPoints1[max2][0]]
				yP = [pathPoints1[max2][1]]
				pylab.scatter(xP,yP, color='b')		   
	
				xP = [tipMatch2[0]]
				yP = [tipMatch2[1]]
				pylab.scatter(xP,yP, color='r')		   

				xP = [pathPoints2[backDepI][0]]
				yP = [pathPoints2[backDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints1[backPathIndex][0]]
				yP = [pathPoints1[backPathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [pathPoints2[newBackDepI][0]]
				yP = [pathPoints2[newBackDepI][1]]
				pylab.scatter(xP,yP, color='m')		   
	
	
			xP = []
			yP = []
			for p in pathPoints1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], (%d,%d)" % (maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, pathID1, pathID2))
			pylab.savefig("pathDeparture_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1
		
		
		return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2

	@logFunction
	def getDeparturePoint(self, currPath, nodeID, plotIter = False):
		
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
		
		
		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.getNodePose(nodeID)
		
		"Assumption:  one section of the medial axis is closely aligned with the path "
		poseOrigin = Pose(estPose2)
		
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		points2 = medialSpline2.getUniformSamples(interpAngle=True)

		points2_offset = []
		for p in points2:
			result = poseOrigin.convertLocalOffsetToGlobal(p)
			points2_offset.append(result)


		pathPoints = orientPath(currPath, points2_offset)

		angle1, angle2 = getTipAngles(points2_offset)


		print "ang1:", angle1
		print "ang2:", angle2
		print "diff:", diffAngle(angle1, angle2)
		
		distSum = 0.0
		contigCount = 0
		maxContig = 0
		distances = []
		indices = []
		for i in range(0,len(points2_offset)):
			p_2 = points2_offset[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints, p_2)
			distances.append(minDist)
			indices.append(i_1)
			distSum += minDist
			
			if minDist < 0.2:
				contigCount += 1
				if contigCount > maxContig:
					maxContig = contigCount
			else:
				contigCount = 0
		
		overlapSum = distSum / float(len(points2_offset))
		
		print "maxContig,overlapSum:", maxContig, overlapSum
		
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


		frontAngleRefI = frontDepI
		while distances[frontAngleRefI] < 0.1:
			frontAngleRefI -= 1
			
			if frontAngleRefI < 0:
				frontAngleRefI = 0
				break
		
		forePathIndex = indices[frontAngleRefI]
		forePathAngle = pathPoints[forePathIndex][2]
		forePathAngle = normalizeAngle(forePathAngle + pi)

		foreDiffAngle = diffAngle(angle1, forePathAngle)
		
		newFrontDepI = 4
		while fabs(diffAngle(normalizeAngle(points2_offset[newFrontDepI][2]+pi), forePathAngle)) > pi/6.0:
			
			if newFrontDepI < frontDepI:
				newFrontDepI += 1
			else:
				break		 
		
		print "front:", nodeID, frontDepI, distances[frontDepI], frontAngleRefI, forePathIndex, forePathAngle, angle1, foreDiffAngle, newFrontDepI, diffAngle(normalizeAngle(points2_offset[newFrontDepI][2]+pi), forePathAngle) 
		

		" FIXME:  index out of bounds case "
		currI = 2
		try:
			while distances[-currI-3] < maxBack:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass
			""" raw sensor data products """

		backDepI = len(distances) - currI
		backPoint = [backDepI, distances[backDepI]]

		backAngleRefI = backDepI
		while distances[backAngleRefI] < 0.1:
			backAngleRefI += 1
			
			if backAngleRefI >= len(distances):
				backAngleRefI = len(distances)-1
				break
		
		backPathIndex = indices[backAngleRefI]
		backPathAngle = pathPoints[backPathIndex][2]
		
		backDiffAngle = diffAngle(angle2, backPathAngle)
		
		newBackDepI = len(distances)-5
		while fabs(diffAngle(points2_offset[newBackDepI][2], backPathAngle)) > pi/6.0:
			
			if newBackDepI > backDepI:
				newBackDepI -= 1
			else:
				break		 

		print "back:", nodeID, backDepI, distances[backDepI], backAngleRefI, backPathIndex, backPathAngle, angle2, backDiffAngle, newBackDepI, diffAngle(points2_offset[newBackDepI][2], backPathAngle) 


		foreAngDiffs = []
		for i in range(0,len(points2_offset)):
			foreAngDiffs.append(fabs(diffAngle(normalizeAngle(points2_offset[i][2] + pi), forePathAngle)))

		backAngDiffs = []
		for i in range(0,len(points2_offset)):
			backAngDiffs.append(fabs(diffAngle(points2_offset[i][2], backPathAngle)))



		"reset to the maximum distances "
		maxFront = distances[0]
		maxBack = distances[-1]



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



		arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
		arr2 = numpy.array(deepcopy(indices[backDepI:]))
		matchVar1 = arr1.var()
		matchVar2 = arr2.var()

		tipMatch1 = pathPoints[indices[0]]

		" discrepancy distance between tip closest point and average departure point "
		dist1 = sqrt((tipMatch1[0]-pathPoints[max1][0])**2 + (tipMatch1[1] - pathPoints[max1][1])**2)

		tipMatch2 = pathPoints[indices[-1]]

		" discrepancy distance between tip closest point and average departure point "
		dist2 = sqrt((tipMatch2[0]-pathPoints[max2][0])**2 + (tipMatch2[1] - pathPoints[max2][1])**2)


		DEP_THRESH = 0.3

		if maxFront > DEP_THRESH:
			departurePoint1 = points2_offset[newFrontDepI]
			isExist1 = True

			if max1 == 0 or max1 == len(pathPoints)-1:
				isInterior1 = False
			else:
				isInterior1 = True

		if maxBack > DEP_THRESH:
			departurePoint2 = points2_offset[newBackDepI]
			isExist2 = True
			
			if max2 == 0 or max2 == len(pathPoints)-1:
				isInterior2 = False
			else:
				isInterior2 = True

		" sum of closest points on front and back "
		" select the one with minimal cost "
		
		if False:
			pylab.clf()
			xP = range(len(points2_offset))
			yP = distances
			pylab.plot(xP,yP, color ='b')
			
			yP = foreAngDiffs
			pylab.plot(xP,yP, color ='r')

			yP = backAngDiffs
			pylab.plot(xP,yP, color ='g')
			
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
			pylab.title("nodeID %d: %1.2f %1.2f %d %d %d %d %1.2f" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2, maxContig, overlapSum))
			pylab.savefig("distances_%04u.png" % self.pathPlotCount)

		if plotIter:
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

				xP = [points2_offset[frontDepI][0]]
				yP = [points2_offset[frontDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints[forePathIndex][0]]
				yP = [pathPoints[forePathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [points2_offset[newFrontDepI][0]]
				yP = [points2_offset[newFrontDepI][1]]
				pylab.scatter(xP,yP, color='m')		   
	
	
			if True:
				xP = [pathPoints[max2][0]]
				yP = [pathPoints[max2][1]]
				pylab.scatter(xP,yP, color='b')		   
	
				xP = [tipMatch2[0]]
				yP = [tipMatch2[1]]
				pylab.scatter(xP,yP, color='r')		   

				xP = [points2_offset[backDepI][0]]
				yP = [points2_offset[backDepI][1]]
				pylab.scatter(xP,yP, color='g')		   

				xP = [pathPoints[backPathIndex][0]]
				yP = [pathPoints[backPathIndex][1]]
				pylab.scatter(xP,yP, color='y')		   

				xP = [points2_offset[newBackDepI][0]]
				yP = [points2_offset[newBackDepI][1]]
				pylab.scatter(xP,yP, color='m')		   

			xP = [points2_offset[-1][0]]
			yP = [points2_offset[-1][1]]
			pylab.scatter(xP,yP, color=(0.5,0.5,0.5))		   

	
			xP = []
			yP = []
			for p in currPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2))
			pylab.savefig("departure_%04u_%04u.png" % (self.hypothesisID, self.pathPlotCount))
			
			self.pathPlotCount += 1
		
		print "departure %d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, frontDepI, backDepI)
		
		
		" if the medial axis does not overlap the path contiguously enough, mark a high discrepancy "
		
		maxContig, overlapSum
		contigFrac = float(maxContig)/float(len(points2_offset))
		print "returning:", contigFrac, overlapSum
		return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum
		

	@logFunction
	def getPathOverlapCondition(self, path1, path2, pathID1, pathID2, plotIter = False):


		if len(path1) == 0:
			return 1e100
	
		medial2 = path2
			   
		minMatchDist2 = 0.2
					
	
		supportSpline = SplineFit(path1, smooth=0.1)		
		vecPoints1 = supportSpline.getUniformSamples()
		supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			 
		
		support_pairs = []
		for i in range(len(vecPoints2)):
			
			p_2 = vecPoints2[i]
	
			try:
				p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
	
				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = supportPoints[minI][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
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
			pylab.title("%d %d cost = %f, count = %d" % (pathID1, pathID2, cost, len(support_pairs)))
			pylab.savefig("pathOverlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 1e100

		return cost


	@logFunction
	def getPathOverlapSum(self, path1, path2, pathID1, pathID2, plotIter = False):


		if len(path1) == 0:
			return 0.0
	
		medial2 = path2
			   
		minMatchDist2 = 0.1
		angThresh = math.pi/8.0
					
	
		supportSpline = SplineFit(path1, smooth=0.1)		
		vecPoints1 = supportSpline.getUniformSamples()
		supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
			
		medialSpline2 = SplineFit(medial2, smooth=0.1)
		vecPoints2 = medialSpline2.getUniformSamples()
		points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

		" transformed points without associated covariance "
		poly2 = []
		for p in points2:
			poly2.append([p[0],p[1]])			 
		
		support_pairs = []
		for i in range(len(vecPoints2)):
			
			p_2 = vecPoints2[i]
	
			try:
				p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, angThresh)
	
				if minDist <= minMatchDist2:
					C2 = points2[i][2]
					C1 = supportPoints[minI][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					support_pairs.append([points2[i],supportPoints[minI],C2,C1])
			except:
				pass

		cost = 0.0
		if len(support_pairs) == 0:
			cost = 0.0
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
				
			#cost = sum1 / len(support_pairs)
			cost = sum1
			

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
				pylab.plot(xP,yP, color='k', alpha=0.3)
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d %d %d cost = %f, count = %d" % (self.hypothesisID, pathID1, pathID2, cost, len(support_pairs)))
			pylab.savefig("pathOverlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 0.0

		return cost


	@logFunction
	def checkUniqueBranch(self, parentPathID, nodeID1, depAngle, depPoint):


		print "checkUniqueBranch(", parentPathID, ",", nodeID1, ",", depAngle, ",", depPoint, ")"

		#foreTerm1 = frontInterior1 and frontExist1
		#foreTerm2 = frontInterior2 and frontExist2

		if depPoint == 0:
			print "REJECT, no departure point!", depPoint
			return False, -1

		" check cartesian distance to similar junction points from parent path "

		" BEING MORE PERMISSIVE IN CREATING NEW BRANCHES BECAUSE WE CAN MERGE LATER "

		" cartesian distance "
		DISC_THRESH = 0.5

		" 60 degree threshold "
		ANG_THRESH = 0.523 # pi/6
		
		" maximum distance to the parent terminal point before creating a junction point "
		TERM_THRESH = 0.8
		
		pathIDs = self.getPathIDs()
		
		for pathID in pathIDs:

			if pathID == parentPathID:
				termPoint = self.terminals[pathID+1]
				
				print "termPoint =", termPoint
				termDist = sqrt((depPoint[0]-termPoint[1][0])**2 + (depPoint[1]-termPoint[1][1])**2)
				
				print "termDist =", termDist
				
				if termDist < TERM_THRESH:
					print "REJECT, proposed junction point is too close to parent pathID terminal", pathID, termDist, termPoint, depPoint
					return False, pathID

			
			path = self.getPath(pathID)
			
			if path["parentID"] == parentPathID:

				
				junctionNodeID = path["branchNodeID"]
				
				junctionPoint = self.getGlobalJunctionPose(pathID)
				
				dist = sqrt((depPoint[0]-junctionPoint[0])**2 + (depPoint[1]-junctionPoint[1])**2 )

				" check difference of tangential angle "
				angDiff = diffAngle(depAngle, junctionPoint[2])

				print "result = ", junctionNodeID, dist, angDiff
		
				#isNodeFeatureless = self.nodeHash[nodeID1].getIsFeatureless()
				isNodeFeatureless = self.poseData.isNodeFeatureless[nodeID1]
				print "node", nodeID1, "is featureless =", isNodeFeatureless
		
				if dist < DISC_THRESH:
					if fabs(angDiff) < ANG_THRESH:
						
						print "DUPLICATE junction of", junctionPoint, "rejecting with differences of", dist, angDiff
						return False, pathID
				
				" Explicit limitation that we wont be able to detect a T-Junction coming from left or right"
				#if isNodeFeatureless:
				#	print "REJECT new branch because the branching node is featureless"
				#	return False, -1
		
		return True, -1



	@logFunction
	def determineBranchSingle(self, nodeID1, frontExist1, frontInterior1, depAngle1, depPoint1, parentPathID1, dirFlag):
	   
		"""
		3) one node is departing only
			- is only one departure ( other node has no departure, or has external departure but not same angle )
			- is a false departure	( falseness happens if bad map data, ignore )
		
		Now, we have the departing angles, but we need to compare them with each other. 
		
		dirFlag specifies which node we want to check branches on.	We do not want to branch with the anchored end of a probe sweep   
		"""

		foreTerm1 = frontInterior1 and frontExist1

		pathBranchID = -1
		isBranch = False
		isNew = False
		
		if foreTerm1:

			isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				
			if isUnique1:
				   
				" foreTerm1 has unique departure "	  
				pathID = parentPathID1
				branchNodeID = nodeID1
				globalJunctionPoint = depPoint1
				depAng = depAngle1
				junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

				poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
				newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
				pathBranchID = newPathID
				isBranch = True
				isNew = True

					
			if duplicatePathID1 != -1:
				pathBranchID = duplicatePathID1
				isBranch = True

		else:
			" no new departure, just add nodes to the leaf paths "
			pathID = parentPathID1
			branchNodeID = nodeID1
			globalJunctionPoint = depPoint1
			
		
		return isBranch, pathBranchID, isNew
		

	@logFunction
	def determineBranchPair(self, nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag, isUnique1, isUnique2, duplicatePathID1, duplicatePathID2):
		
		"""
		1)	both nodes departing
			- are same departure  ( equal departing angle and are on top of each other)  1 new path
			- are different departures ( different departing angle and are not necessarily on top of each other) 2 new paths
			
		2) neither node is departing
			- are on same path
			
		3) one node is departing only
			- is only one departure ( other node has no departure, or has external departure but not same angle )
			- is a false departure	( falseness happens if bad map data, ignore )
			- are both departing ( should have at least external departure on both, compare departure angle )
		
		Now, we have the departing angles, but we need to compare them with each other. 
		
		dirFlag specifies which node we want to check branches on.	We do not want to branch with the anchored end of a probe sweep   
		"""

		print "determineBranchPair(", nodeID1, ",", nodeID2, ",", frontExist1, ",", frontExist2, ",", frontInterior1, ",", frontInterior2, ",", depAngle1, ",", depAngle2, ",", depPoint1, ",", depPoint2, ",", parentPathID1, ",", parentPathID2, ",", dirFlag, ")"
		
		foreTerm1 = frontInterior1 and frontExist1
		foreTerm2 = frontInterior2 and frontExist2

		frontAngDiff = diffAngle(depAngle1, depAngle2)

		pathBranchIDs = [-1, -1]
		isBranch = [False, False]
		isNew = [False, False]
		
		" 60 degree threshold "
		ANG_THRESH = 1.047


		if foreTerm1 and foreTerm2:
			" both departing case: check if same or different "
			if fabs(frontAngDiff) < ANG_THRESH:
				" are on top of each other and are branching to the same path "

				#isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				#isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

				if isUnique1 and isUnique2:
					
					pathID = parentPathID1
					branchNodeID = nodeID1
					globalJunctionPoint = depPoint1
					depAng = depAngle1
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					
					newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID

					isBranch[0] = True
					isBranch[1] = True
					isNew[0] = True
					isNew[1] = True
					
					print "foreA"

				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True
						isNew[0] = True

				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True
						isNew[1] = True
					
				
				if duplicatePathID1 != -1:
					pathBranchIDs[0] = duplicatePathID1
					isBranch[0] = True
					

				if duplicatePathID2 != -1:
					pathBranchIDs[1] = duplicatePathID2
					isBranch[1] = True
					
				
			else:
				" these are branching to separate new paths "

				#isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				#isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

				if isUnique1 and isUnique2:

					if dirFlag == 0:
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID1 = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID1
						isBranch[0] = True
						isNew[0] = True


					if dirFlag == 1:
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID2 = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID2
						isBranch[1] = True
						isNew[1] = True

					print "foreB"

				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True
						isNew[0] = True

				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True
						isNew[1] = True

				if duplicatePathID1 != -1:
					pathBranchIDs[0] = duplicatePathID1
					isBranch[0] = True
					
				if duplicatePathID2 != -1:
					pathBranchIDs[1] = duplicatePathID2
					isBranch[1] = True
					
		
		elif foreTerm1 and not foreTerm2:

			#isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
				
			if isUnique1:
				   
				if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
					" both have same departure "
					pathID = parentPathID1
					branchNodeID = nodeID1
					globalJunctionPoint = depPoint1
					depAng = depAngle1
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID
					isBranch[0] = True
					isBranch[1] = True
					isNew[0] = True
					isNew[1] = True

					print "foreC"
					
				else:
					if dirFlag == 0:
						" foreTerm1 has unique departure "	  
						pathID = parentPathID1
						branchNodeID = nodeID1
						globalJunctionPoint = depPoint1
						depAng = depAngle1
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[0] = newPathID
						isBranch[0] = True
						isNew[0] = True

					print "foreD"
					
			if duplicatePathID1 != -1:
				pathBranchIDs[0] = duplicatePathID1
				isBranch[0] = True

		elif foreTerm2 and not foreTerm1:

			#isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

			if isUnique2:
				if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

					" both have same departure "
					pathID = parentPathID2
					branchNodeID = nodeID2
					globalJunctionPoint = depPoint2
					depAng = depAngle2
					junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

					poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
					newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
					pathBranchIDs[0] = newPathID
					pathBranchIDs[1] = newPathID
					isBranch[0] = True
					isBranch[1] = True
					isNew[0] = True
					isNew[1] = True

					print "foreE"

				else:
					
					if dirFlag == 1:
						" foreTerm2 has unique departure "	  
						pathID = parentPathID2
						branchNodeID = nodeID2
						globalJunctionPoint = depPoint2
						depAng = depAngle2
						junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

						poseOrigin = Pose(self.nodeRawPoses[branchNodeID])
						newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
						pathBranchIDs[1] = newPathID
						isBranch[1] = True
						isNew[1] = True

					print "foreF"

			if duplicatePathID2 != -1:
				pathBranchIDs[1] = duplicatePathID2
				isBranch[1] = True

		"""
		" CHECKING FOR BAD PATH MEMBERSHIP HERE.  May reuse later "
		else:
			print "foreG"
			" no new departure, just add nodes to the leaf paths "
			pathID = parentPathID1
			branchNodeID = nodeID1
			globalJunctionPoint = depPoint1
			
			if parentPathID1 != parentPathID2:
				print "departing path IDs of two colocated nodes are not the same"
				raise
		"""	
		
		return isBranch, pathBranchIDs, isNew
		
		
	@logFunction
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
			parent = self.getPath(currPathID)
			thisParentID = parent["parentID"]
			
			if thisParentID == None:
				break
			
			if pathIDs.count(thisParentID) == 0:
				break
			
			currPathID = thisParentID

		" currPathID is the root path "
		rootPathID = currPathID
		
			
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


			for join in self.joins:
				if join[0][0] == childPathID:
					childJunctionPoint = join[2]
					childI = join[0][1]

			if childI < abs(childI-len(childPath)-1):
				childTermI = len(childPath)-1
			else:
				childTermI = 0

			"find path between: (childPathID, childTermI) to (rootPathID, 0) and (rootPathID, len(rootPath)-1)"
			shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path((childPathID, childTermI))
			startNode1 = shortestPathSpanTree[(rootPathID, 0)]
			startNode2 = shortestPathSpanTree[(rootPathID, len(rootPath)-1)]
			
			currNode = startNode1
			splicedPath1 = []
			while currNode != (childPathID, childTermI):
				splicedPath1.append(self.pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath1.append(self.pathGraph.get_node_attributes(currNode))

			currNode = startNode2
			splicedPath2 = []
			while currNode != (childPathID, childTermI):
				splicedPath2.append(self.pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath2.append(self.pathGraph.get_node_attributes(currNode))

			newPaths.append(splicedPath1)
			newPaths.append(splicedPath2)
			
		else:
			"find entire path from rightPathID to leftPathID "
			" start from point farthest from child's junction point "
			
			rightPath = self.trimmedPaths[rightPathID]
			leftPath = self.trimmedPaths[leftPathID]
			
			" find the junction of this child's to its parent "    
			for join in self.joins:
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
			shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path((rightPathID, rightTermI))
			currNode = shortestPathSpanTree[(leftPathID, leftTermI)]

			splicedPath = []
			while currNode != (rightPathID, rightTermI):
				splicedPath.append(self.pathGraph.get_node_attributes(currNode))
				currNode = shortestPathSpanTree[currNode]
			splicedPath.append(self.pathGraph.get_node_attributes(currNode))

			newPaths.append(splicedPath)





		return newPaths

	@logFunction
	def getNearestPathPoint(self, originPoint):

		" find closest point on each path "
		minI1 = 0
		minJ1 = 0
		minDist1 = 1e100
		
		pathIDs = self.getPathIDs()    
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for j in range(len(path)):
				p0 = path[j]
				dist = sqrt((originPoint[0]-p0[0])**2 + (originPoint[1]-p0[1])**2)
				
				if dist < minDist1:
					minDist1 = dist
					minI1 = pathID
					minJ1 = j
		
		newPoint = self.trimmedPaths[minI1][minJ1]
		
		return newPoint
		
	@logFunction
	def computeNavigationPath(self, startPose, endPose):
		
		" get the trimmed paths "
		
		" find closest point on each path "
		minI1 = 0
		minJ1 = 0
		minDist1 = 1e100
		minI2 = 0
		minJ2 = 0
		minDist2 = 1e100
		
		pathIDs = self.getPathIDs()    
		for pathID in pathIDs:
			path = self.trimmedPaths[pathID]
			for j in range(len(path)):
				p0 = path[j]
				dist = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
				
				if dist < minDist1:
					minDist1 = dist
					minI1 = pathID
					minJ1 = j

				dist = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

				if dist < minDist2:
					minDist2 = dist
					minI2 = pathID
					minJ2 = j
		
		" select shortest one "
		
		" get origin and target path "
		
		" get splice of origin and target path "
		" splice minJ1 and minJ2 "
		orderPathIDs = self.getPathPath(minI1, minI2)		 
		
		splicedPaths = self.splicePathIDs(orderPathIDs)
		
		print len(splicedPaths), "spliced paths for path finding "
		
		" again, find closest point of start and end pose "
		
		minTotalDist = 3e100
		minSplicePathID = 0
		spliceIndex = (0,1)
		for i in range(len(splicedPaths)):

			minDist1 = 1e100
			minDist2 = 1e100

			minJ1 = 0
			minJ2 = 0

			
			path = splicedPaths[i]

			for j in range(len(path)):
				p0 = path[j]
				dist1 = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
				
				if dist1 < minDist1:
					minDist1 = dist1
					minJ1 = j

				dist2 = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

				if dist2 < minDist2:
					minDist2 = dist2
					minJ2 = j

			
			totalDist = minDist1 + minDist2
			if totalDist < minTotalDist:
				minTotalDist = totalDist
				minSplicePathID = i
				spliceIndex = (minJ1,minJ2)

			print i, path[0], path[-1], minJ1, minJ2, minDist1, minDist2, totalDist, minTotalDist, minSplicePathID, spliceIndex
		
		
		" return splicePath[startIndex:endIndex+1]"
		if spliceIndex[0] < spliceIndex[1]:
			pathSpline = SplineFit(splicedPaths[minSplicePathID])
			newPath = pathSpline.getUniformSamples()
			return newPath
		else:
			path = splicedPaths[minSplicePathID]
			path.reverse()
			pathSpline = SplineFit(path)
			newPath = pathSpline.getUniformSamples()
			return newPath

			return path

			

