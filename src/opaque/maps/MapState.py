

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
from LocalNode import getLongestPath, computePathSegments
import gen_icp
from SplineFit import SplineFit
import pylab
import numpy
from operator import itemgetter
import hashlib
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath, getTipAngles
from MapProcess import selectLocalCommonOrigin
from ParticleFilter import multiParticleFitSplice, batchParticle, Particle
import time
import traceback

import alphamod


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
			childPath = job[3]
			parentPath = job[4]
			trimmedParent = job[5]
			smoothPathSegs = job[6]
			arcDist = job[7]

			result = computeBranch(pathID, parentID, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs, arcDist)

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
def getTangentIntersections(path1, path2, frontDepI, backDepI, path1FrontDepI, path1BackDepI, path2JuncI, path1JuncI, plotCount):


	""" for tangent on path2, find the intersection point on path1 """

	interPoints = []
	indices1 = []
	edges = []

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


	for i in range(frontBoundI, frontDepI):

		p = path2[i]

		""" generate tangent segment """
		angle2 = p[2]

		pA = [p[0] + 3*cos(angle2), p[1] + 3*sin(angle2)]
		pB = [p[0] - 3*cos(angle2), p[1] - 3*sin(angle2)]

		edge2 = [pA,pB]
		edges.append(edge2)

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

		pA = [p[0] + 3*cos(angle2), p[1] + 3*sin(angle2)]
		pB = [p[0] - 3*cos(angle2), p[1] - 3*sin(angle2)]

		edge2 = [pA,pB]
		edges.append(edge2)

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
	print "interPoints:", interPoints



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



	if False:
		pylab.clf()
		xP = []
		yP = []
		for p in path2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.5,0.5,1.0),zorder=500)

		if True:	

			""" draw the tangents """
			for edge in edges:
				xP = [edge[0][0], edge[1][0]]
				yP = [edge[0][1], edge[1][1]]
				pylab.plot(xP,yP, color=(0.5,1.0,0.5), linewidth=1, alpha=0.5,zorder=1)		   

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

		pylab.title("%d intersections, %d tangent segments, angles %1.2f %1.2f %1.2f %1.2f" % (len(interPoints), len(edges), juncForeAng1, juncBackAng1, juncForeAng, juncBackAng))
		pylab.savefig("intersectDeparture_%04u.png" % plotCount)
		

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
def getOverlapDeparture(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False):

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
	juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

	p0 = pathSec2[0]
	juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

	print "pathSec1 hypothesis discrepancy distance:", juncDist1
	print "pathSec2 hypothesis discrepancy distance:", juncDist2

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

			
		pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')		   

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
def getBranchPoint(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False):
	""" get the trimmed version of child and parent paths that are overlapping in some fashion """

	"Assumption:  one section of the medial axis is closely aligned with the path "		   
	
	pathPlotCount = 0
	print "getBranchPoint():"
	
	" return exception if we receive an invalid path "		  
	if len(path1) == 0:
		print "path1 has zero length"
		raise
	
	" make sure the overlap of both paths are oriented the same way "
	orientedPath2 = orientPath(path2, path1)
	
	path1Spline = SplineFit(path1, smooth=0.1)			  
	path2Spline = SplineFit(orientedPath2, smooth=0.1)


	" for each point on the child path, find its closest pair on the parent path "
	
	pathPoints1 = path1Spline.getUniformSamples(interpAngle=True)
	pathPoints2 = path2Spline.getUniformSamples(interpAngle=True)

	distances = []
	indices = []
	juncDists = []

	minDist2 = 1e100
	juncI = 0		 

	for i in range(0,len(pathPoints2)):
		p_2 = pathPoints2[i]
		p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
		
		" keep the distance information "
		distances.append(minDist)
		" and the associated index of the point on the parent path "
		indices.append(i_1)
		
		juncDist = sqrt((p_2[0]-globalJunctionPoint[0])**2 + (p_2[1]-globalJunctionPoint[1])**2)

		if juncDist < minDist2:
			minDist2 = juncDist
			juncI = i

		juncDists.append(juncDist)


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
			backFound = k
			backFound = True

	print "frontFound, backFound:", frontFound, backFound
	print "frontInd, backInd:", frontInd, backInd

	print "minDist2, juncI:", minDist2, juncI

	" match distances of the tip points of child path "		   
	maxFront = distances[frontInd]
	maxBack = distances[backInd]
	print "match distances of tip points:", maxFront, maxBack

	" walk back from tip point until we have a non-monotic increase in match distance "
	" this becomes our departure point "
	
	"TODO:	why is there a 3-offset in this comparison? "
	currI = frontInd + 1
	try:
		while distances[currI+3] < maxFront or distances[currI+3] > 0.1:
			maxFront = distances[currI]
			currI += 1
	except:
		pass
	
	" departure index on child path "
	frontDepI = currI

	if frontDepI+3 >= len(pathPoints2)-1:
		frontDepI = juncI

	
	" departure index of child path and match distance "
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


	" FIXME:  index out of bounds case "
	currI = 1 + len(pathPoints2) - backInd
	try:
		while distances[-currI-3] < maxBack or distances[-currI-3] > 0.1:
			maxBack = distances[-currI]
			currI += 1
	except:
		pass

	" departure index on child path "
	backDepI = len(distances) - currI

	if backDepI-3 <= 0: 
		backDepI = juncI

	" departure index of child path and match distance "
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


	print "frontDepI, backDepI:", frontDepI, backDepI
	print "frontAngleRefI, backAngleRefI:", frontAngleRefI, backAngleRefI
	print "forePathAngle, backPathAngle:", forePathAngle, backPathAngle
	print "newFrontDepI, newBackDepI:", newFrontDepI, newBackDepI
	print "foreDiff, backDiff:", diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]+pi), forePathAngle), diffAngle(normalizeAngle(pathPoints2[newBackDepI][2]), backPathAngle)

	print "frontDiffAngles:", frontDiffAngles
	print "backDiffAngles:", backDiffAngles


	print "lengths of parent and child paths:", len(pathPoints1), len(pathPoints2)
	print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint


	"reset to the tip match distance "
	maxFront = distances[frontInd]
	maxBack = distances[backInd]
	
	
	" sum of closest points on front and back "
	" select the one with minimal cost "		

	" path section for our front departure hypothesis "
	pathSec1 = pathPoints2[:frontDepI+1]
	pathSec1.reverse()

	" path section for our back departure hypothesis "
	pathSec2 = pathPoints2[backDepI:]
	
	" distance of departure point from known junction point "
	#juncDist1 = -1
	#juncDist2 = -1
	#if len(pathSec1) > 0:
	p0 = pathSec1[0]
	juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

	#if len(pathSec2) > 0:
	p0 = pathSec2[0]
	juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

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


	junctionPoint = globalJunctionPoint		 
	minDist2 = 1e100
	juncI = 0		 
	for i in range(len(pathPoints2)):
		pnt = pathPoints2[i]
		dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
	
		if dist < minDist2:
			minDist2 = dist
			juncI = i
	
	
	
	if juncDist1 < juncDist2:
		" FIXME: if [0,frontDepI] has high overlap compared to total length, then we reject it "
		secP1 = pathPoints2[0]
		secP2 = pathPoints2[frontDepI]
		
		newPath2 = pathPoints2[:newFrontDepI]
		newPath2.reverse()

	else:
		" FIXME: if [backDepI,len(distances)-1] has high overlap compared to total length, then we reject it "
		secP1 = pathPoints2[backDepI]
		secP2 = pathPoints2[len(distances)-1]

		newPath2 = pathPoints2[newBackDepI:]

	if len(secP1) == 0:
		print "no departures found"
		raise
	
	" convert path so that the points are uniformly distributed "
	
	
	max_spacing = 0.08
	newPath3 = []

	" make sure path is greater than 5 points "
	while len(newPath3) <= 5:

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
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					newPath3.append(newP)

			newPath3.append(copy(p1))			 
	


	leafPath = deepcopy(newPath3)
	
	frontVec = [0.,0.]
	backVec = [0.,0.]
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


	newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)

	leafPath.insert(0,newP1)
	
	medial2 = deepcopy(leafPath)

	" take the long length segments at tips of medial axis"
	edge1 = medial2[0:2]
	
	frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	
	frontVec[0] /= frontMag
	frontVec[1] /= frontMag
	
	" make a smaller version of these edges "
	newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)

	edge1 = [newP1, edge1[1]]

	" find the intersection points with the hull "
	interPoints = []
	for k in range(len(path1)-1):
		hullEdge = [path1[k],path1[k+1]]
		isIntersect1, point1 = Intersect(edge1, hullEdge)
		if isIntersect1:
			interPoints.append(point1)
			break

	
	" replace the extended edges with a termination point at the hull edge "			
	medial2 = medial2[1:]
	
	if isIntersect1:
		medial2.insert(0, point1)


	juncAng = acos(-frontVec[0])
	if -frontVec[1] < 0.0:
		juncAng = -juncAng
	
	

	try:
		foreIntI, backIntI, juncForeAng, juncBackAng = getTangentIntersections(pathPoints1, pathPoints2, frontDepI, backDepI, indices[frontDepI], indices[backDepI], juncI, indices[juncI], pathPlotCount)
		print "foreIntI, backIntII:", foreIntI, backIntI

		""" junction distances are equal if both of the indices selected juncI as the departure point on path2'
			occurs if path2 does not come close enough to path1 to get under distance 0.1
		"""
		print "juncDist1, juncDist2 =", juncDist1, juncDist2
		if juncDist1 == juncDist2:

			juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

			foreDist = sqrt((pathPoints1[foreIntI][0]-globalJunctionPoint[0])**2 +  (pathPoints1[foreIntI][1]-globalJunctionPoint[1])**2)
			backDist = sqrt((pathPoints1[backIntI][0]-globalJunctionPoint[0])**2 +  (pathPoints1[backIntI][1]-globalJunctionPoint[1])**2)

			print "foreDist, backDist =", foreDist, backDist

			if foreDist < backDist:
				globJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], juncForeAng]
			else:
				globJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], juncBackAng]
	
		elif juncDist1 < juncDist2:
			globJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], juncForeAng]
		else:
			globJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], juncBackAng]

		print "globJuncPose =", globJuncPose

	except:
		print "getTangentIntersections() failed!"
		globJuncPose = [medial2[0][0], medial2[0][1], juncAng]
	
	" last junction point is the intersection point "

	#if plotIter:
	if False:

		pylab.clf()
		xP = []
		yP = []
		for p in path2:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.5,0.5,1.0))

		if True:	
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


			
		pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')		   
		pylab.scatter([globJuncPose[0]],[globJuncPose[1]], color='k')		 

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

		pylab.title("%3.2f %3.2f %3.2f" % (juncDist1, juncDist2, juncAng))
		pylab.savefig("trimDeparture_%04u.png" % pathPlotCount)
		
		pathPlotCount += 1


	
	return globJuncPose

@logFunction
def computeBranch(pathID, parentID, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs, arcDist):
#def trimBranch(pathID, parentPathID, globJuncPose, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs):

	if parentID != None: 

		" path information "
		pathSpline = SplineFit(parentPath)

		" initial position of junction "
		#origJuncPose = copy(pathDesc["globalJunctionPose"])
		origJuncPose = copy(origGlobJuncPose)
		origJuncPose[2] = 0.0

		" position on the spline curve of the parent path "
		newArcDist = arcDist

		" NOTE:  angle is the tangent angle, not branch angle "
		newJuncPose = pathSpline.getPointOfDist(newArcDist)
		modJuncPose = copy(newJuncPose)
		modJuncPose[2] = origJuncPose[2]

		"""
		1) point that has the junction
		2) direction of junction
		3) set of long paths that include junction path  (only 2)

		"""


		#medialLongPath = self.medialLongPaths[pathID]

		#junctionDetails = self.longPathJunctions[pathID]


		globalPath1 = childPath
		origJuncPose = copy(origGlobJuncPose)
		origJuncPose[2] = 0.0
		origJuncOrigin = Pose(origJuncPose)

		localPathSegs = []
		#smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

		#print "path segs:", len(junctionDetails["leafSegments"]), "+", len(junctionDetails["internalSegments"]), "=", len(smoothPathSegs)

		for k in range(len(smoothPathSegs)):
			pathSeg = smoothPathSegs[k]
			localSeg = []
			for p in pathSeg:
				p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
				localSeg.append(p1)

			localPathSegs.append(localSeg)

		args = []
		xP = []
		yP = []

		offsetOrigin1 = Pose(modJuncPose)

		placedPathSegs = []
		for k in range(len(localPathSegs)):
			localSeg = localPathSegs[k]
			placedSeg = []
			for p in localSeg:
				p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
				placedSeg.append(p1)
			placedPathSegs.append(placedSeg)


		#parentPath = self.paths[parentID]
		juncOrigin1 = Pose(modJuncPose)

		" curves of both paths "
		globalPath2 = parentPath

		globalSpline2 = SplineFit(globalPath2, smooth=0.1)
		globalSamples2 = globalSpline2.getUniformSamples(spacing = 0.04)
		originU2 = globalSpline2.findU(modJuncPose)
		minDist, originU2, uPoint = globalSpline2.findClosestPoint(modJuncPose)
		
		u2 = originU2
		u1 = 0.5  # value is meaningless since data is not a spline curve
		angGuess = -uPoint[2] + modJuncPose[2]
		
		initGuess = [u1,u2,angGuess]
		resultPose, lastCost, matchCount = gen_icp.branchEstimateCost(initGuess, modJuncPose, placedPathSegs, globalPath2, plotIter = False, n1 = pathID, n2 = parentID)
		result = (resultPose, lastCost, matchCount)	

		resultPose1 = result[0]
		lastCost1 = result[1]
		matchCount1 = result[2]

		poseOrigin = Pose(resultPose1)

		inversePose = poseOrigin.doInverse(resultPose1)
		poseOrigin = Pose(inversePose)

		finalPose = poseOrigin.convertLocalOffsetToGlobal(modJuncPose)
		poseOrigin2 = Pose(finalPose)

		origJuncPose = copy(origGlobJuncPose)
		origJuncPose[2] = 0.0
		#distDisc = sqrt((finalPose[0]-origJuncPose[0])**2 + (finalPose[1]-origJuncPose[1])**2)
		distDisc = 0.0
		angDisc = fabs(normalizeAngle(finalPose[2]-origJuncPose[2]))

		#junctionDetails = self.longPathJunctions[pathID]
		#smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

		#newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths =  trimBranch(pathID, parentID, modJuncPose, pathDesc["globalJunctionPose"], self.paths[pathID], self.paths[parentID], self.trimmedPaths[parentID], smoothPathSegs)
		newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths =  trimBranch(pathID, parentID, modJuncPose, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs)

		newSplices = deepcopy(splicedPaths)

		initProb = 0.0

		part2 = (parentID, modJuncPose, newArcDist, pathID, matchCount1, lastCost1, distDisc, angDisc, initProb, newSplices, juncDiscAngle)

		return part2
		#currBranchSpace[thisKey] = part2

		#print "completed precomputation of branch", thisKey

		#return deepcopy(currBranchSpace[thisKey])

@logFunction
def trimBranch(pathID, parentPathID, globJuncPose, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs):

	origJuncPose = copy(origGlobJuncPose)
	origJuncPose[2] = 0.0

	origJuncOrigin = Pose(origJuncPose)

	localPathSegs = []

	for k in range(len(smoothPathSegs)):
		pathSeg = smoothPathSegs[k]
		localSeg = []
		for p in pathSeg:
			p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
			localSeg.append(p1)

		localPathSegs.append(localSeg)


	offsetOrigin1 = Pose(globJuncPose)
	
	partPathSegs = []

	for k in range(len(localPathSegs)):
		localSeg = localPathSegs[k]
		for p in localSeg:
			p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)

	path1 = parentPath
	path2 = childPath

	localPath2 = []
	for p in path2:
		p1 = origJuncOrigin.convertGlobalToLocal(p)
		localPath2.append(p1)

	particlePath2 = []
	for p in localPath2:
		p1 = offsetOrigin1.convertLocalToGlobal(p)
		particlePath2.append(p1)

	" get departing sections of overlapped curves "
	secP1, secP2 = getOverlapDeparture(globJuncPose, parentPathID, pathID, path1, particlePath2, plotIter = False)				 

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

	" smallest distance is the terminal point "
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
		dist = sqrt((pnt[0]-globJuncPose[0])**2 + (pnt[1]-globJuncPose[1])**2)
	
		if dist < minDist2:
			minDist2 = dist
			juncI = i
	 
	 
	
	print "len(path1):", len(path1)
	print "len(particlePath2):", len(particlePath2)
	print "juncI:", juncI
	print "minDist:", minDist_1, minDist_2
	
	" now we have closest point to departure point. "
	" Which side is the departing side? "	 

			
	if minI == 0:
		"secP1 is terminal 0"
		index = juncI-10
		if index < 1:
			index = 1

		
		newPath2 = particlePath2[:index+1]
		newPath2.reverse()
	
	elif minI == 1:
		"secP1 is terminal N"
		index = juncI+10
		if index >= len(particlePath2)-1:
			
			" ensure at least 2 elements in path "
			index = len(particlePath2)-2

		newPath2 = particlePath2[index:]

		
	elif minI == 2:
		"secP2 is terminal 0"
		index = juncI-10
		if index < 1:
			index = 1

		newPath2 = particlePath2[:index+1]
		newPath2.reverse()
		
	elif minI == 3:
		"secP2 is terminal N"
		index = juncI+10
		if index >= len(particlePath2)-1:
			" ensure at least 2 elements in path "
			index = len(particlePath2)-2
		
		newPath2 = particlePath2[index:]
	
	else:
		print "no terminal found"
		raise
					
	" convert path so that the points are uniformly distributed "
	
	
	max_spacing = 0.08
	newPath3 = []

	" make sure path is greater than 5 points "
	while len(newPath3) <= 5:

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
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					newPath3.append(newP)

			newPath3.append(copy(p1))			 

	deepcopy(newPath3)

	newGlobJuncPose = getBranchPoint(globJuncPose, parentPathID, pathID, path1, particlePath2, plotIter = False)

	" junction discrepency distance "
	origJuncPose = copy(origGlobJuncPose)
	juncDiscDist = sqrt((newGlobJuncPose[0]-globJuncPose[0])**2 + (newGlobJuncPose[1]-globJuncPose[1])**2)
	juncDiscAngle = normalizeAngle(origJuncPose[2]-newGlobJuncPose[2])




	" for each path, attempt to join with its parent path "
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


	" get junctions " 
	
	junctions = {}
	junctions[pathID] = [0, newGlobJuncPose, (parentPathID,minI2), path2[minI2], minI1]


	" create a tree with the node IDs and then stitch them together with joins "
	pathGraph = graph.graph()

	for k in range(len(path1)):
		pathGraph.add_node((pathID, k), path1[k])
	for k in range(len(path1)-1):
		pathGraph.add_edge((pathID, k), (pathID, k+1))
	
	for k in range(len(path2)):
		pathGraph.add_node((parentPathID, k), path2[k])
	for k in range(len(path2)-1):
		pathGraph.add_edge((parentPathID, k), (parentPathID, k+1))

	" join with the junction in between the join points "
	for k in range(len(joins)):
		join = joins[k]
		pathGraph.add_edge(join[0], join[1])

	" get terminals "
	topDict = {}
	terminals = {}


	" parent is either the root path or some sub path "
	" we only care about the two terminals of our immediate parent "


	terminals[parentPathID] = [(parentPathID, 0), path2[0]]
	terminals[parentPathID+1] = [(parentPathID, len(path2)-1), path2[len(path2)-1]]
	topDict["t%u" % 0] = terminals[parentPathID][0]
	topDict["t%u" % 1] = terminals[parentPathID+1][0]
		
	minI1 = junctions[pathID][4]
	
	" get the two terminals of the current pathID path "
	" determine which side is junction and which side is terminal"
	if minI1 > len(path1)-1 - minI1:    
		terminals[pathID+1] = [(pathID,0), path1[0]]
		topDict["t%u" % (pathID+1)] = (pathID,0)
	else:
		terminals[pathID+1] = [(pathID,len(path1)-1), path1[len(path1)-1]]
		topDict["t%u" % (pathID+1)] = (pathID,len(path1)-1)



	" determine which paths are leaves "
	isAParent = {}
	parents = {}
	pathIDGraph = graph.graph()

	isAParent[pathID] = False
	pathIDGraph.add_node(pathID, [])

	isAParent[parentPathID] = True
	pathIDGraph.add_node(parentPathID, [])
	pathIDGraph.add_edge(pathID, parentPathID)
	
	leafCount = 1

	" build topological graph "
	topGraph = graph.graph()
	
	" add topological nodes to the graph "
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

	" if there are any joins, we must interpolate a smooth transition "
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

	" if there are any joins, we must interpolate a smooth transition "
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


	#if plotIter:
	if False:
			
		hypothesisID = 0
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
		pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='b')

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
		pylab.title("hyp %d nodeID %d %1.3f %1.3f %1.3f %1.3f" % ( hypothesisID, numNodes, juncDiscDist, juncDiscAngle, origJuncPose[2], newGlobJuncPose[2]))
		pylab.savefig("trimDeparture_%04u_%04u.png" % (hypothesisID, pathPlotCount))

		print "saving trimDeparture_%04u_%04u.png" % (hypothesisID, pathPlotCount)
		
		pathPlotCount += 1


	return newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths

def multiParticleFitSplice(initGuess, orientedPath, medialAxis, initPose, pathIDs, nodeID, pathPlotCount = 0):
	
	u1 = initGuess[0]
	u2 = initGuess[1]
	angGuess = initGuess[2]

	#minMedialDist0, oldMedialU0, oldMedialP0 = medialSpline.findClosestPoint([0.0,0.0,0.0])

	#uPath1, uMedialOrigin1 = selectLocalCommonOrigin(orientedSplicePath, medial1, pose1)
	#u1 = uPath1
	#resultPose0, lastCost0, matchCount0, currAng0, currU0 = gen_icp.globalPathToNodeOverlapICP2([uPath0, uMedialOrigin0, 0.0], currSplice0, medial0, plotIter = False, n1 = nodeID0, n2 = -1, arcLimit = 0.01)

	#localizeJobs.append([originU2, pathU1, 0.0, spliceIndex, orientedPath, medial1, hypPose, [], nodeID, particleIndex])

	resultPose0, lastCost0, matchCount0, currAng0, currU0 = gen_icp.globalPathToNodeOverlapICP2([u1, u2, angGuess], orientedPath, medialAxis, plotIter = False, n1 = nodeID, n2 = -1, arcLimit = 0.01, origPose = initPose)


	icpDist = sqrt((resultPose0[0]-initPose[0])**2 + (resultPose0[1]-initPose[1])**2)
	#print "nodeID:", nodeID0, nodeID1
	#print "travelDist:", thisDist, travelDist0, travelDist1
	#print "icpDist:", nodeID, pathPlotCount, icpDist, currU0, u1, currAng0, angGuess, u2

	resultArgs = getMultiDeparturePoint(orientedPath, medialAxis, initPose, resultPose0, pathIDs, nodeID, pathPlotCount, plotIter = False)
	#resultArgs = ([0.0, 0.0],  0.0, False, False, 0.0, 0.0, [0.0, 0.0], 0.0, False, False, 0.0, 0.0, 1.0, 0.0, 0.0)

	isExist1 = resultArgs[3]
	isExist2 = resultArgs[9]

	" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 "

	return (resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs + (isExist1 or isExist2,)




class MapState:
	
	def __init__(self, poseData, hypothesisID):
		
		self.poseData = deepcopy(poseData)
		

		" branch point particles "
		self.numJuncSamples = 20
		self.juncSamples = {}

		self.numPoseParticles = 20
		self.poseParticles = {}
		self.poseParticles["numParticles"] = self.numPoseParticles
		self.poseParticles["updateCount"] = 0
		self.poseParticles["snapshots2"] = {0 : ()}


		self.nodeRawPoses = {}
		self.nodePoses = {}
		self.gndPoses = {}
		self.gndRawPoses = {}
		self.origPoses = {}

		self.utility = 0.0

		" criteria for rejecting this map state "
		self.mapOverlapSum = 0.0
		self.isNoLocalize = False

		self.pathIDs = 0
		self.hypothesisID = hypothesisID

		self.paths = {0 : []}
		self.hulls = {0 : []}
		self.pathTermsVisited = {0: False}
		self.trimmedPaths  = {}
		self.pathClasses = {}
		self.pathClasses[0] = {"parentID" : None, "branchNodeID" : None, "localJunctionPose" : None, 
							"sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : None}
		self.pathIDs += 1


		self.longPathJunctions = {}
		self.medialLongPaths = {}
		self.theoryMedialLongPaths = {}

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



		self.rootPoint = [-3.2, 0.0]
		
		self.pathGraph = graph.graph()
		self.joins = []
		self.junctions = {}
		self.terminals = {}
		self.allSplices = {}
		self.isChanged = True

		self.colors = []

		self.colors.append([127, 127, 255])
		self.colors.append([127, 255, 127])
		self.colors.append([255, 127, 127])
		self.colors.append([255, 127, 255])
		self.colors.append([255, 255, 127])
		self.colors.append([127, 255, 255])
		"""
		self.colors.append([215, 48, 39])
		self.colors.append([252, 141, 89])
		self.colors.append([254, 224, 144])
		self.colors.append([224, 243, 248])
		self.colors.append([145, 191, 219])
		self.colors.append([69, 117, 180])
		"""

		for color in self.colors:
			color[0] = float(color[0])/256.0
			color[1] = float(color[1])/256.0
			color[2] = float(color[2])/256.0

		for i in range(1000):
			self.colors.append((random.random(),random.random(),random.random()))
	
	@logFunction
	def initializePoseParticles(self):

	
		pose0 = self.nodePoses[0]
		pose1 = self.nodePoses[1]

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

			"""
			minMedialDist0, oldMedialU0, oldMedialP0 = medialSpline.findClosestPoint([0.0,0.0,0.0])

			#uPath1, uMedialOrigin1 = selectLocalCommonOrigin(orientedSplicePath, medial1, pose1)
			#u1 = uPath1
			#resultPose0, lastCost0, matchCount0, currAng0, currU0 = gen_icp.globalPathToNodeOverlapICP2([uPath0, uMedialOrigin0, 0.0], currSplice0, medial0, plotIter = False, n1 = nodeID0, n2 = -1, arcLimit = 0.01)
			resultPose0, lastCost0, matchCount0, currAng0, currU0 = gen_icp.globalPathToNodeOverlapICP2([newU0, oldMedialU0, 0.0], currSplice0, medial0, plotIter = False, n1 = nodeID0, n2 = -1, arcLimit = 0.01)

			icpDist = sqrt((resultPose0[0]-newPose[0])**2 + (resultPose0[1]-newPose[1])**2)
			print "nodeID:", nodeID0, nodeID1
			print "travelDist:", thisDist, travelDist0, travelDist1
			print "icpDist:", icpDist, currU0, newU0, currAng0

			#resultPose1, lastCost1, matchCount1, currAng1, currU1 = gen_icp.globalPathToNodeOverlapICP2([uPath1, uMedialOrigin1, 0.0], orientedSplicePath, medial1, plotIter = False, n1 = nodeID, n2 = -1, arcLimit = 0.1)

			#newDist0 = pathSpline.dist_u(currU0)
			#minDist0, oldU0, oldP0 = pathSpline.findClosestPoint(hypPose)

			newPart = (resultPose0, 0, newDist0, ('t0', 't1'))
			"""

			newPartDist2.append(particleObj)


		self.poseParticles["snapshots2"][0] = newPartDist2

	#@logFunction
	def localizePoseParticles(self, nodeID0, nodeID1):

		#cProfile.run('runTest(probe)', 'test_prof')


		self.isNoLocalize = False
		self.resetBranches()

		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]
		poseData = self.poseData

		#nodeID = nodeID0

		nodePose0 = self.nodePoses[nodeID0]
		staticSplicedPaths, spliceTerms, staticSplicePathIDs = self.getSplicesByNearJunction(nodePose0)

		pathIDs = self.getPathIDs()

		medial0 = poseData.medialAxes[nodeID0]
		medial1 = poseData.medialAxes[nodeID1]

		if nodeID0-2 >= 0:
			prevMedial0 = poseData.medialAxes[nodeID0-2]
			prevMedial1 = poseData.medialAxes[nodeID1-2]
		else:
			prevMedial0 = None
			prevMedial1 = None


		#medialSamples0 = medialSpline0.getUniformSamples(spacing = 0.04)
		#medialSamples1 = medialSpline1.getUniformSamples(spacing = 0.04)

		medialSpline0 = SplineFit(medial0, smooth=0.1)
		medialSpline1 = SplineFit(medial1, smooth=0.1)
		minMedialDist0, oldMedialU0, oldMedialP0 = medialSpline0.findClosestPoint([0.0,0.0,0.0])
		minMedialDist1, oldMedialU1, oldMedialP1 = medialSpline1.findClosestPoint([0.0,0.0,0.0])


		""" let's cache the trimBranch results so we don't compute more than once """
		"""
		trimBranchCache = {}
		for particleIndex in range(len(particleDist2)):

			part = particleDist2[particleIndex]
			for pathID in pathIDs:
				if pathID != 0:
					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]

					thisKeys = trimBranchCache.keys()
					if repr(globJuncPose) not in thisKeys:
					
						newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, branchSplices =  self.trimBranch(pathID, globJuncPose)

						trimBranchCache[repr(globJuncPose)] = branchSplices

					#for j in range(len(branchSplices)):
					#	splice = branchSplices[j]
					#	thisSplicedPaths.append((0, deepcopy(splice)))
		"""


		#orientedPaths = []
		#orientedPaths.append(orientedPath)

		localizeJobs = []
		allSplicedPaths = []
		spliceCount = 0


		"""
		duplicateMatches = {}
		for particleIndex in range(len(particleDist2)):

			isDuplicate = False

			for i in range(0,particleIndex):
				if particleDist2[i] == particleDist2[particleIndex]:
					duplicateMatches[particleIndex] = i
					isDuplicate = True
					break

			if not isDuplicate:
				duplicateMatches[particleIndex] = particleIndex

		#print "duplicateMatches:", duplicateMatches.values()
		"""
		probSets = []
		for particleIndex in range(len(particleDist2)):
			part = particleDist2[particleIndex]
			hypPose0 = part.pose0
			pathIDs = self.getPathIDs()

			for pathID in pathIDs:
				if pathID != 0:

					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]

					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)
					if dist1 < 3.0:
						probSets.append((pathID, globJuncPose))

		self.batchEvalBranches(probSets)

		for particleIndex in range(len(particleDist2)):

			time1 = time.time()

			part = particleDist2[particleIndex]

			hypPose0 = part.pose0
			hypPose1 = part.pose1
			pathID = part.memberPaths[0]
			hypDist = part.hypDist
			prevHypPose0 = part.prevPose0
			prevHypPose1 = part.prevPose1

			resultsBySplice = []

			thisSplicedPaths = []

			for pathID in pathIDs:
				if pathID != 0:
					globJuncPose = part.junctionData[pathID]["globalJunctionPose"]

					dist1 = sqrt((globJuncPose[0]-hypPose0[0])**2 + (globJuncPose[1]-hypPose0[1])**2)
					if dist1 < 3.0:

						#time1 = time.time()

						samples = self.evaluateBranches(pathID, globJuncPose)

						probDist = []
						branchPoseDist = []
						branchSplices = []
						for branchSamp in samples:
							probDist.append(branchSamp[8])
							branchPoseDist.append(branchSamp[1])
							branchSplices.append(branchSamp[9])

						probStr = ""
						for probVal in probDist:
							probStr += "%1.2f " % probVal

						print self.hypothesisID, particleIndex, pathID, "branch probDist:", probStr

						#self.poseParticles["snapshots2"][updateCount][particleIndex].junctionData[pathID] = {}
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["probDist"] = probDist
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["branchPoseDist"] = branchPoseDist
						self.poseParticles["snapshots2"][0][particleIndex].junctionData[pathID]["branchSplices"] = branchSplices


						#time2 = time.time()
						#print "evaluateBranches time", self.hypothesisID, particleIndex, pathID, time2-time1

						for j in range(len(samples)):

							branchPose = samples[j][1]
							probVal = samples[j][8]
							#branchSplices = samples[j][9]
							numSplices = len(branchSplices[j])

							dist2 = sqrt((branchPose[0]-hypPose0[0])**2 + (branchPose[1]-hypPose0[1])**2)
							if dist2 < 3.0:

								for k in range(numSplices):
									splice = branchSplices[j][k]
									thisSplicedPaths.append((j, probVal, splice))


			" consider single path splices "
			" FIXME: uses static currPath, should use hypothesize currPath instead "
			for k in range(len(staticSplicePathIDs)):
				if len(staticSplicePathIDs[k]) == 1:
					" static splices have 1.0 probability "
					thisSplicedPaths.append((None, 1.0, staticSplicedPaths[k]))


			print "particle:", particleIndex, ",",  len(thisSplicedPaths), "localize jobs"
			for spliceIndex in range(len(thisSplicedPaths)):
				
				branchSampleIndex = thisSplicedPaths[spliceIndex][0]
				probVal = thisSplicedPaths[spliceIndex][1]
				
				path = thisSplicedPaths[spliceIndex][2]
				
				#poseOrigin = Pose(hypPose0)
				#
				#globalMedial0 = []
				#for p in medial0:
				#	globalMedial0.append(poseOrigin.convertLocalToGlobal(p))
				#
				#globalMedial1 = []
				#for p in medial1:
				#	globalMedial1.append(poseOrigin.convertLocalToGlobal(p))
				#globalMedialP0 = poseOrigin.convertLocalOffsetToGlobal(oldMedialP0)	
				#globalMedialP1 = poseOrigin.convertLocalOffsetToGlobal(oldMedialP1)	
				#orientedPath0 = orientPath(path, globalMedial0)
				#orientedPathSpline0 = SplineFit(orientedPath0, smooth=0.1)
				#pathU0 = orientedPathSpline0.findU(globalMedialP0)
				#pathU1 = orientedPathSpline0.findU(globalMedialP1)


				" add guesses in the neighborhood of the closest pathU "
				" 5 forward, 5 backward "
				
			
				localizeJobs.append([oldMedialP0, oldMedialU0, 0.0, oldMedialP1, oldMedialU1, 0.0, branchSampleIndex, spliceCount, path, medial0, medial1, deepcopy(hypPose0), deepcopy(hypPose1), prevMedial0, prevMedial1, prevHypPose0, prevHypPose1, [], nodeID0, nodeID1, particleIndex, updateCount, self.hypothesisID])
				#localizeJobs.append([pathU0, oldMedialU0, 0.0, pathU1, oldMedialU1, 0.0, branchSampleIndex, spliceCount, orientedPath0, medial0, medial1, deepcopy(hypPose0), deepcopy(hypPose1), [], nodeID0, nodeID1, particleIndex, updateCount, self.hypothesisID])

				self.pathPlotCount2 += 1
				spliceCount += 1



			allSplicedPaths += thisSplicedPaths

			time2 = time.time()
			print "localize construct", particleIndex, time2-time1

		print len(localizeJobs), "total localize jobs"

		time1 = time.time()

		results = batchParticle(localizeJobs)

		time2 = time.time()

		print "batch localize for map hypothesis", self.hypothesisID, "=", time2-time1 
		print len(results[0]), "result arguments"
		#for k in range(len(results)):
		#	resultArgs = results[k]


		#results.sort
		#sortedResults = sorted(results, key=itemgetter(0,1), reverse=False)

		" sort by pose particle index followed by utility value "
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

		newParticleDist2 = []

		print "particle evaluation:", nodeID0, self.hypothesisID, updateCount
		for particleIndex in range(len(filteredParticles)):

			#if duplicateMatches[particleIndex] == particleIndex:

			part = filteredParticles[particleIndex]

			utilVal = part[45]
			spliceIndex = part[47]
			branchProbVal = allSplicedPaths[spliceIndex][1]
			spliceCurve = allSplicedPaths[spliceIndex][2]

			initPose0 = part[48]
			initPose0 = part[49]

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
			" divergence is prohibited "
			if isInterior1_0 or isInterior2_0 or isInterior1_1 or isInterior2_1:
				isReject = True

			" only a single extension is permitted, but not two "
			#if isExist1_0 and isExist2_0 or isExist1_1 and isExist2_1:
			if isExist1_0 and isExist2_0 or isExist1_1 and isExist2_1:
				isReject = True
			
			" horrifically low contiguity is rejected out of hand "
			if contigFrac_0 <= 0.5 or contigFrac_1 <= 0.5:
				isReject = True

			" contiguity between current and previous pose "
			if overlapSum > 1e10:
				isReject = True
				
			" probability is the contigFrac squared " 
			newProb = 0.0
			if not isReject:
				newProb = (pi-fabs(diffAngle(initPose0[2],newPose0[2])) ) * contigFrac_0 * contigFrac_0 / overlapSum_0
				newProb *= branchProbVal
				#utilVal0 = (1.0-contigFrac_0) + (isExist1_0 or isExist2_0) + (1.0-contigFrac_1) + (isExist1_1 or isExist2_1)
				
			print "%d %d %d branchProbVal, utilVal, poseProbValue, overlapSum:" % (self.hypothesisID, particleIndex, spliceIndex), branchProbVal, utilVal, newProb, overlapSum 

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

			newParticleDist2.append(particleObj)

			#else:
			#	dupIndex = duplicateMatches[particleIndex]
			#	newParticleDist2.append(newParticleDist2[dupIndex].copy())


		updateCount += 1
		self.poseParticles["updateCount"] = updateCount
		self.poseParticles["snapshots2"][0] = newParticleDist2
		numParticles = self.poseParticles["numParticles"]

		self.drawPoseParticles()

		" now resample the particles "
		particleDist2 = self.poseParticles["snapshots2"][0]
		probSum = 0.0
		for part in particleDist2:
			probVal = part.weightVal
			probSum += probVal

		probParticles = []
		for k in range(len(particleDist2)):
			part = particleDist2[k]

			poseProbVal = part.weightVal

			" if all 0.0 probabilities, make uniform distribution, set no localization flag "
			if probSum > 0.0:
				probParticles.append(poseProbVal/probSum)
			else:
				probParticles.append(float(1)/float(numParticles))
				self.isNoLocalize = True


		" apply branch probability "
		probSum2 = 0.0
		for k in range(len(probParticles)):
			#part = particleDist2[k]
			part = filteredParticles[particleIndex]

			" branch probability "
			spliceIndex = part[47]
			branchProbVal = allSplicedPaths[spliceIndex][1]

			probParticles[k] *= branchProbVal

			probSum2 += probParticles[k]

		" renormalize "
		for k in range(len(probParticles)):
			if probSum2 > 0.0:
				probParticles[k] /= probSum2
			else:
				probParticles[k] = float(1)/float(numParticles)

		" now resample "
		resampledParticles2 = []
		numParticles = self.poseParticles["numParticles"]

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

			print "resampling particle", k, probParticles[k], printStr
			resampledParticles2.append(oldPart2.copy())

		updateCount += 1
		self.poseParticles["updateCount"] = updateCount
		self.poseParticles["snapshots2"][0] = resampledParticles2

		self.drawPoseParticles()

		#self.isChanged = True
		#self.junctions = {}
		#self.terminals = {}
		#self.allSplices = {}

		return

		print len(orientedPaths), "orientedPaths", len(splicedPaths1), "splicedPaths", [ initGuesses[guessIndex][3] for guessIndex in range(len(initGuesses))]
		print len(initGuesses), "initGuesses"
		#results = batchGlobalMultiFit(initGuesses, orientedPaths, medial1, estPose1, [], nodeID)
		results = batchGlobalMultiFit(initGuesses, orientedPaths, medial1, estPose1, [], nodeID)

		time2 = time.time()
		print "time:", time2 - time1

		resultsBySplice += results
		
		print "consistentFit node", nodeID
		
		print len(resultsBySplice), "results"
		
		for k in range(len(resultsBySplice)):
			result = resultsBySplice[k]
			resultPose = result[0]
			dist = sqrt((estPose1[0]-resultPose[0])**2 + (estPose1[1] - resultPose[1])**2)
			resultsBySplice[k] = result + (dist, k,)


		print "unfilteredResults:"
		for result in resultsBySplice:
			resultPose = result[0]
			isInterior1 = result[5]
			isExist1 = result[6]
			isInterior2 = result[11]
			isExist2 = result[12]
			contigFrac = result[15]
			angDiff2 = result[17]
			spliceIndex = result[19]
			#spliceIndex = result[21]

			#result[18] is isDepartureExists True/False
			dist = result[20]
			#dist = sqrt((estPose1[0]-resultPose[0])**2 + (estPose1[1] - resultPose[1])**2)
			lastCost = result[1]
			matchCount = result[2]
			overlapSum = result[16]
			print spliceIndex, [int(isInterior1), int(isInterior2), int(isExist1), int(isExist2)], contigFrac, dist, angDiff2, lastCost, matchCount, overlapSum, spliceTerms[spliceIndex], splicePathIDs[spliceIndex], resultPose

		
		" resultPose, lastCost, matchCount, departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2, isDeparture, dist, spliceIndex "


		"""
		1) departures:	(0,0) good, (1,0) okay, (1,1) bad 
		2) contigFrac:	1.0 good, 0.9 okay, 0.5 bad, 0.0 reject
		3) dist:   0.0 great, 0.5 good, 1.0 okay, 2.0 bad, 3.0 reject
		4) angDiff2:  0.0 great, 0.2 good, 0.5 bad/reject
		5) splicePathIDs:  exID in REJECT, currPathID in GOOD
		"""
		
		#(0.5-angDiff2)/0.5 + (contigFrac-0.5)/0.5 + (3.0-dist)/3.0 + 1-isDeparture*0.5
		
		filteredResults = []
		for result in resultsBySplice:
			resultPose = result[0]
			isInterior1 = result[5]
			isExist1 = result[6]
			isInterior2 = result[11]
			isExist2 = result[12]
			contigFrac = result[15]
			angDiff2 = result[17]
			isDeparture = result[18]
			spliceIndex = result[19]
			dist = result[20]
			
			
			" rejection filter "
			isReject = False
			
			if isInterior1 or isInterior2:
				isReject = True
			

			if isExist1 and isExist2:
				isReject = True
			
			for exID in excludePathIDs:
				if exID in splicePathIDs[spliceIndex]:
					isReject = True
		
			if contigFrac <= 0.5:
				isReject = True
				
			if dist >= 3.0:
				isReject = True
			
			if angDiff2 > 0.7:
				isReject = True
				
			if not isReject:
				filteredResults.append(result)

			
	
	@logFunction
	def drawPoseParticles(self):

		updateCount = self.poseParticles["updateCount"] 
		#particleDist = self.poseParticles["snapshots"][updateCount]
		particleDist2 = self.poseParticles["snapshots2"][0]


		localPathSets = {}

		allPathIDs = self.getPathIDs()
		for pathID in allPathIDs:

			if pathID != 0:

				junctionDetails = self.longPathJunctions[1]

				origJuncPose = copy(self.pathClasses[1]["globalJunctionPose"])
				origJuncPose[2] = 0.0
				origJuncOrigin = Pose(origJuncPose)

				smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]
				localPathSegs = []
				for k in range(len(smoothPathSegs)):
					pathSeg = smoothPathSegs[k]
					localSeg = []
					for p in pathSeg:
						p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
						localSeg.append(p1)

					localPathSegs.append(localSeg)

				#print "path segs:", len(junctionDetails["leafSegments"]), "+", len(junctionDetails["internalSegments"]), "=", len(smoothPathSegs)

				localPathSets[pathID] = localPathSegs


		#print "particles:", particleDist2

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


				if pathID != 0:


					if pathID in part.junctionData.keys():

						probDist = part.junctionData[pathID]["probDist"]
						branchPoseDist = part.junctionData[pathID]["branchPoseDist"]


						maxProb = -1e100
						maxIndex = 0
						for k in range(len(probDist)):
							if probDist[k] > maxProb:
								maxProb = probDist[k]
								maxIndex = k


						modPose0 = deepcopy(branchPoseDist[maxIndex])
					else:
						modPose0 = deepcopy(part.junctionData[pathID]["globalJunctionPose"])

					modPose0[2] = 0.0
					modOrigin0 = Pose(modPose0)

					if thisSplice != None:

						xP1 = []
						yP1 = []

						for p1 in thisSplice:
							xP1.append(p1[0])
							yP1.append(p1[1])

						pylab.plot(xP1,yP1,color='k', zorder=9, alpha=0.2)

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

					pylab.scatter([modPose0[0]], [modPose0[1]], color='k', linewidth=1, zorder=10, alpha=0.4)




			poseOrigin1 = Pose(hypPose1)

			xP = []
			yP = []
			for p in medial1:
				p1 = poseOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])

			pylab.plot(xP,yP, color = 'b', alpha = 0.2, zorder=9)

		pylab.scatter(hypPointsX_0, hypPointsY_0, color='r', linewidth=1, zorder=10, alpha=0.2)
		pylab.scatter(hypPointsX_1, hypPointsY_1, color='b', linewidth=1, zorder=10, alpha=0.2)
		
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

		newObj.numJuncSamples = deepcopy(self.numJuncSamples)
		newObj.numPoseParticles = deepcopy(self.numPoseParticles)

		newObj.juncSamples = deepcopy(self.juncSamples)
		newObj.poseParticles = deepcopy(self.poseParticles)

		updateCount = newObj.poseParticles["updateCount"] 
		particleDist2 = newObj.poseParticles["snapshots2"][0]
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
		newObj.junctions = self.junctions
		newObj.terminals = self.terminals
		newObj.allSplices = self.allSplices

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

		newObj.medialLongPaths = deepcopy(self.medialLongPaths)
		newObj.theoryMedialLongPaths = deepcopy(self.theoryMedialLongPaths) 

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
					resultSum = self.getPathOverlapSum(path1, path2, j, k, plotIter = False)
					print "computing sum of", j, "and", k, "=", resultSum
					totalSum += resultSum

		self.mapOverlapSum = totalSum
		return totalSum
	
		utilSum = 0.0

		for nodeID1, estPose1 in self.nodePoses.iteritems():
	
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


	@logFunction
	def addBranch(self):
		pass

	@logFunction
	def computeConsistency(self):
		pass

	@logFunction
	def updatePose(self, nodeID, pose):
		pass

	
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
	def getParentPathID(self, pathID):
		parentID = self.pathClasses[pathID]["parentID"]
		return parentID
	
	@logFunction
	def getPath(self, pathID):
		return self.pathClasses[pathID]
	
	#@logFunction
	def getAllSplices(self, plotIter = False):

		if not self.isChanged:
			print "returning without computing, not isChanged"
			return self.allSplices, self.terminals, self.junctions

		print "paths:", self.paths

		print "junctions:", self.junctions
		print "terminals:", self.terminals

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
					sPath['path'] = newPath
					sPath['termPath'] = termPath
					
					finalResults.append(sPath)
					
				results[(pathID1,pathID2)] = finalResults
		
		print "results:"
		for k, result in results.iteritems():
			print k, ":"
			for sPath in result:
				print sPath['orderedPathIDs'], sPath['termPath'], len(sPath['path'])


		if plotIter:
			pylab.clf()
			
			for k,path in  self.trimmedPaths.iteritems():
				print "path has", len(path), "points"
				xP = []
				yP = []
				for p in path:
					xP.append(p[0])
					yP.append(p[1])
	
				pylab.plot(xP,yP, color = self.colors[k])
	
	
			for k, result in results.iteritems():
				for sPath in result:
					path = sPath['path']
	
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])
					pylab.plot(xP,yP, color='b')


			for k in pathIDs:
				xP = []
				yP = []
				
				for p in self.paths[k]:
					xP.append(p[0])
					yP.append(p[1])
				pylab.plot(xP,yP, color=self.colors[k], linewidth=1)


			#self.drawWalls()
	
			pylab.title("Spliced Paths, pathIDs = %s" %  self.getPathIDs())
			pylab.savefig("splicedPath_%04u.png" % self.spliceCount)
			self.spliceCount += 1

		self.allSplices = results
		
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

		DIV_LEN = 0.2

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
					currDist += DIV_LEN
					branchSpace[currDist] = None

				self.branchEvaluations[pathID] = branchSpace

	@logFunction
	def computeBranch(self, pathID, arcDist):

		if self.pathClasses[pathID]["parentID"] != None: 
			currBranchSpace = self.branchEvaluations[pathID]

			currKeys = currBranchSpace.keys()
			currKeys.sort()

			" find the bin for this arc distance "
			thisKey = currKeys[0]
			for val in currKeys:
				if val >= arcDist:
					break
				thisKey = val

			if currBranchSpace[thisKey] != None:
				print "returned precomputation of branch", thisKey
				return deepcopy(currBranchSpace[thisKey])

			return None

			" path information "
			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]
			pathSpline = SplineFit(self.paths[parentID])

			" initial position of junction "
			origJuncPose = copy(pathDesc["globalJunctionPose"])
			origJuncPose[2] = 0.0


			" position on the spline curve of the parent path "
			newArcDist = thisKey

			" NOTE:  angle is the tangent angle, not branch angle "
			newJuncPose = pathSpline.getPointOfDist(newArcDist)
			modJuncPose = copy(newJuncPose)
			modJuncPose[2] = origJuncPose[2]

			"""
			1) point that has the junction
			2) direction of junction
			3) set of long paths that include junction path  (only 2)

			"""


			#medialLongPath = self.medialLongPaths[pathID]

			junctionDetails = self.longPathJunctions[pathID]


			globalPath1 = self.paths[pathID]
			origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
			origJuncPose[2] = 0.0
			origJuncOrigin = Pose(origJuncPose)

			localPathSegs = []
			smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

			print "path segs:", len(junctionDetails["leafSegments"]), "+", len(junctionDetails["internalSegments"]), "=", len(smoothPathSegs)

			for k in range(len(smoothPathSegs)):
				pathSeg = smoothPathSegs[k]
				localSeg = []
				for p in pathSeg:
					p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
					localSeg.append(p1)

				localPathSegs.append(localSeg)

			args = []
			xP = []
			yP = []

			offsetOrigin1 = Pose(modJuncPose)

			placedPathSegs = []
			for k in range(len(localPathSegs)):
				localSeg = localPathSegs[k]
				placedSeg = []
				for p in localSeg:
					p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
					placedSeg.append(p1)
				placedPathSegs.append(placedSeg)


			parentPath = self.paths[parentID]
			juncOrigin1 = Pose(modJuncPose)

			" curves of both paths "
			globalPath2 = parentPath

			globalSpline2 = SplineFit(globalPath2, smooth=0.1)
			globalSamples2 = globalSpline2.getUniformSamples(spacing = 0.04)
			originU2 = globalSpline2.findU(modJuncPose)
			minDist, originU2, uPoint = globalSpline2.findClosestPoint(modJuncPose)
			
			u2 = originU2
			u1 = 0.5  # value is meaningless since data is not a spline curve
			angGuess = -uPoint[2] + modJuncPose[2]
			
			initGuess = [u1,u2,angGuess]
			resultPose, lastCost, matchCount = gen_icp.branchEstimateCost(initGuess, modJuncPose, placedPathSegs, globalPath2, plotIter = False, n1 = pathID, n2 = parentID)
			result = (resultPose, lastCost, matchCount)	
	
			resultPose1 = result[0]
			lastCost1 = result[1]
			matchCount1 = result[2]

			poseOrigin = Pose(resultPose1)

			inversePose = poseOrigin.doInverse(resultPose1)
			poseOrigin = Pose(inversePose)

			finalPose = poseOrigin.convertLocalOffsetToGlobal(modJuncPose)
			poseOrigin2 = Pose(finalPose)

			origJuncPose = copy(pathDesc["globalJunctionPose"])
			origJuncPose[2] = 0.0
			#distDisc = sqrt((finalPose[0]-origJuncPose[0])**2 + (finalPose[1]-origJuncPose[1])**2)
			distDisc = 0.0
			angDisc = fabs(normalizeAngle(finalPose[2]-origJuncPose[2]))

			junctionDetails = self.longPathJunctions[pathID]
			smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

			newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths =  trimBranch(pathID, parentID, modJuncPose, pathDesc["globalJunctionPose"], self.paths[pathID], self.paths[parentID], self.trimmedPaths[parentID], smoothPathSegs)

			newSplices = deepcopy(splicedPaths)

			initProb = 0.0

			part2 = (parentID, modJuncPose, newArcDist, pathID, matchCount1, lastCost1, distDisc, angDisc, initProb, newSplices, juncDiscAngle)

			currBranchSpace[thisKey] = part2

			print "completed precomputation of branch", thisKey

			return deepcopy(currBranchSpace[thisKey])


	@logFunction
	def batchEvalBranches(self, probSets):

		#probSets.append((pathID, globJuncPose))
		probSets = sorted(probSets, key=itemgetter(0), reverse=True)

		allPathIDs = self.getPathIDs()

		pathSplines = {}
		for pathID in allPathIDs:

			" path information "
			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			if parentID != None:

				pathSpline = SplineFit(self.paths[parentID])

				pathSplines[parentID] = pathSpline

		arcDists = []
		for prob in probSets:

			pathID = prob[0]
			estJuncPose = prob[1]
			minDist, uVal, splinePoint = pathSpline.findClosestPoint(estJuncPose)
			arcDist = pathSpline.dist_u(uVal)


			arcHigh = arcDist + 0.5
			arcLow = arcDist - 0.5

			" precompute the problem space "
			#for k in range(11):
			for k in range(6):
				newArcDist = arcLow + (arcHigh-arcLow)*k/10
				arcDists.append((pathID, newArcDist))


		" binned arc distances, discretized state space "
		binnedProblems = []

		for prob in arcDists:
			pathID = prob[0]
			arcDist = prob[1]

			currBranchSpace = self.branchEvaluations[pathID]
			currKeys = currBranchSpace.keys()
			currKeys.sort()

			" find the bin for this arc distance "
			thisKey = currKeys[0]
			for val in currKeys:
				if val >= arcDist:
					break
				thisKey = val

			binnedProblems.append((pathID, thisKey))

		numProbs1 = len(binnedProblems)
		binnedProblems = list(set(binnedProblems))
		numProbs2 = len(binnedProblems)

		print "binnedProblems:", numProbs1, numProbs2
		branchJobs = []

		for prob in binnedProblems:
			pathID = prob[0]
			arcDist = prob[1]

			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			junctionDetails = self.longPathJunctions[pathID]
			smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

			origGlobJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
			childPath = self.paths[pathID]
			parentPath = self.paths[parentID]
			trimmedParent = self.trimmedPaths[parentID]

			branchJobs.append((pathID, parentID, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs, arcDist))

			#result = computeBranch(pathID, parentID, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs, arcDist)
			#currBranchSpace = self.branchEvaluations[pathID]
			#currBranchSpace[arcDist] = result

		#part2 = (parentID, modJuncPose, newArcDist, pathID, matchCount1, lastCost1, distDisc, angDisc, initProb, newSplices, juncDiscAngle)

		results = batchBranch(branchJobs)
		for result in results:
			arcDist = result[2]
			pathID = result[3]
			currBranchSpace = self.branchEvaluations[pathID]
			currBranchSpace[arcDist] = result

		#def computeBranch(pathID, parentID, origGlobJuncPose, childPath, parentPath, trimmedParent, smoothPathSegs, arcDist):
		#particles.append(particle)

	@logFunction
	def evaluateBranches(self, pathID, estJuncPose):

		# pathID
		# parentID
		# curr path
		# parent path
		# origJuncPose
		# junctionDetails = self.longPathJunctions[pathID]
		# smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]
		# hypID

		if self.pathClasses[pathID]["parentID"] != None: 

			" path information "
			pathDesc = self.pathClasses[pathID]
			parentID = pathDesc["parentID"]

			pathSpline = SplineFit(self.paths[parentID])

			" initial position of junction "
			origJuncPose = copy(pathDesc["globalJunctionPose"])
			origJuncPose[2] = 0.0

			#minDist, uVal, splinePoint = pathSpline.findClosestPoint(origJuncPose)
			#arcDist = pathSpline.dist_u(uVal)
			minDist, uVal, splinePoint = pathSpline.findClosestPoint(estJuncPose)
			arcDist = pathSpline.dist_u(uVal)

			arcHigh = arcDist + 0.5
			arcLow = arcDist - 0.5

			totalDist = pathSpline.dist_u(1.0)

			#1) point that has the junction
			#2) direction of junction
			#3) set of long paths that include junction path  (only 2)

			" initialize if this is the first time "
			particles = []
			for k in range(6):
				newArcDist = arcLow + (arcHigh-arcLow)*k/10
				particle = self.computeBranch(pathID, newArcDist)
				particles.append(particle)

			" get the maximum value for each of our features "
			maxCost = -1e100
			maxMatchCount = -1000
			#maxDist = -1e100
			maxAngDiff = -1e100
			maxDist = 10.0

			" find maximum values for each metric "
			for k in range(len(particles)):
				part = particles[k]
				matchCount = part[4]
				lastCost = part[5]
				dist = part[6]
				angDiff = part[7]

				if matchCount > maxMatchCount:
					maxMatchCount = matchCount

				if lastCost > maxCost:
					maxCost = lastCost

				if dist > maxDist:
					maxDist = dist

				if angDiff > maxAngDiff:
					maxAngDiff = angDiff

			" invert the distance, angDiff and cost features "
			totalProbSum = 0.0
			for k in range(len(particles)):

				part = particles[k]

				lastCost = part[5]
				dist = part[6]
				angDiff = part[7]

				matchCost = 0.0
				if maxCost > 0.0:
					matchCost = (maxCost-lastCost)/maxCost
				else:
					matchCost = 0.0

				nomDist = 0.0
				if maxDist > 0.0:
					nomDist = (maxDist-dist)/maxDist
				else:
					nomDist = 0.0

				matchCount = part[4]

				" angDiff feature times matchCount times dist "
				#probVal = matchCount*matchCost + nomDist/4.0
				probVal = matchCount*matchCost
				totalProbSum += probVal

				part2 = (part[0], part[1], part[2], part[3], part[4], matchCost, nomDist, part[7], probVal, part[9], part[10])
				particles[k] = part2

			" normalize the probability values "
			maxProb = 0.0
			for k in range(len(particles)):
				part = particles[k]

				probVal = 0.0

				if totalProbSum > 0.0:

					probVal = part[8] / totalProbSum

					if probVal > maxProb:
						maxProb = probVal
				else:
					probVal = 0.0

				part2 = (part[0], part[1], part[2], part[3], part[4], part[5], part[6], part[7], probVal, part[9], part[10])
				particles[k] = part2

				print "particle %02u %1.4f %03u %1.4f %1.4f %1.4f %1.4f %d" % (k, part2[1][2], part2[4], part2[5], part2[6], part2[7], part2[8], len(part2[9]))
			print "particle"

			" cartesian distance "
			DISC_THRESH = 0.5

			" 60 degree threshold "
			ANG_THRESH = 0.523 # pi/6

			
			" get splices for new branch position " 
			" reject branch locations whose branch angle discrepancy is too high " 
			for k in range(len(particles)):
				part = particles[k]

				thisPathID = part[3]
				globJuncPose = part[1]

				#newSplices = deepcopy(splicedPaths)
				newProbVal = particles[k][8] 

				juncDiscAngle = part[10]

				if fabs(juncDiscAngle) > ANG_THRESH:
					newProbVal = 0.0

				part2 = (part[0], part[1], part[2], part[3], part[4], part[5], part[6], part[7], newProbVal, part[9], part[10])
				particles[k] = part2

			pylab.clf() 
			for k,path in self.trimmedPaths.iteritems():
				print "path has", len(path), "points"
				xP = []
				yP = []
				for p in path:
					xP.append(p[0])
					yP.append(p[1])

				pylab.plot(xP,yP, color = self.colors[k], linewidth=4)


			origJuncPose = copy(self.pathClasses[pathID]["globalJunctionPose"])
			origJuncPose[2] = 0.0
			origJuncOrigin = Pose(origJuncPose)

			junctionDetails = self.longPathJunctions[pathID]
			smoothPathSegs = junctionDetails["leafSegments"] + junctionDetails["internalSegments"]

			print "path segs:", len(junctionDetails["leafSegments"]), "+", len(junctionDetails["internalSegments"]), "=", len(smoothPathSegs)

			localPathSegs = []
			for k in range(len(smoothPathSegs)):
				pathSeg = smoothPathSegs[k]
				localSeg = []
				for p in pathSeg:
					p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
					localSeg.append(p1)

				localPathSegs.append(localSeg)

			xP = []
			yP = []
			for part in particles:
				globJuncPose = part[1]
				xP.append(globJuncPose[0])
				yP.append(globJuncPose[1])

				offsetOrigin1 = Pose(globJuncPose)

				for k in range(len(localPathSegs)):
					localSeg = localPathSegs[k]
					xP1 = []
					yP1 = []
					for p in localSeg:
						p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
						xP1.append(p1[0])
						yP1.append(p1[1])
					pylab.plot(xP1,yP1,color='k', zorder=9, alpha=part[8]/maxProb)


			pylab.scatter(xP, yP, color='k', zorder=8)
			pylab.title("hyp: %d, pathID: %d, localPathSegs %d" % (self.hypothesisID, pathID, len(localPathSegs)))
			
			#self.plotEnv()			
			
			pylab.savefig("bayes_plot_%04u_%04u.png" % (self.hypothesisID, self.tempCount) )
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

		#globJuncPose = getBranchPoint(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False)
		#print "generated globJuncPose:",  globJuncPose

		return samplePath2
				
	""" generate the paths from the node and pose information """
	@logFunction
	def generatePaths(self):
		
		" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
		self.paths = {}
		self.hulls = {}
		
		pathIDs = self.getPathIDs()
		for k in pathIDs:
			print "computing path for node set", k, ":", self.getNodes(k)
			self.paths[k], self.hulls[k] = self.getTopology(k)

		" compute the junction points between parent paths and child branches "
		for pathID in pathIDs:
			
			if self.pathClasses[pathID]["parentID"] != None:

				parentPathID = self.pathClasses[pathID]["parentID"]
				childPathID = pathID

				path1 = self.paths[parentPathID]
				path2 = self.paths[childPathID]
				branchNodeID = self.pathClasses[childPathID]["branchNodeID"]
				localJunctionPoint = self.pathClasses[childPathID]["localJunctionPose"]
				poseOrigin2 = Pose(self.nodeRawPoses[branchNodeID])
				globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)

				globJuncPose = getBranchPoint(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False)
				#globJuncPose = getBranchPoint(parentPathID, childPathID, self.paths, plotIter = False)
				print "generated globJuncPose:",  globJuncPose
				
				self.pathClasses[childPathID]["globalJunctionPose"] = globJuncPose



		self.trimmedPaths = self.trimPaths(self.paths)

		" for each path, attempt to join with its parent path "
		self.joins = []
		self.junctions = {}
		for pathID in pathIDs:

			cPath = self.getPath(pathID)
			
			parentPathID = cPath["parentID"]
			
			" parent does not concern us "
			if parentPathID == None:
				continue
			
			#junctionNodeID = cPath["branchNodeID"]
			#localJunctionPoint = cPath["localJunctionPose"]
			#poseOrigin = Pose(self.nodeRawPoses[junctionNodeID])
			#junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)
			#junctionPose = poseOrigin.convertLocalOffsetToGlobal(localJunctionPoint)

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
				self.pathGraph.add_edge((pathID, k), (pathID, k+1))
			
			" parent does not concern us "
			cPath = self.getPath(pathID)			
			parentPathID = cPath["parentID"]
			
		" join with the junction in between the join points "
		for k in range(len(self.joins)):
			join = self.joins[k]
			self.pathGraph.add_edge(join[0], join[1])

			pID1, k1 = join[0]
			pID2, k2 = join[1]


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


		if len(self.trimmedPaths[pathID]) > 0:
			allSplices, terminals, junctions = self.getAllSplices(plotIter = False)

		isChanged = False

	@logFunction
	def getOrderedOverlappingPaths(self, nodeID):

		nodePose = self.nodePoses[nodeID]
		splicePaths, spliceTerms, splicePathIDs = self.getSplicesByNearJunction(nodePose)
		
		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.nodePoses[nodeID]


		resultSet = []

		for k in range(len(splicePaths)):
			
			print "comparing node", nodeID, "to splice", splicePathIDs[k], spliceTerms[k],  "plot", self.multiDepCount
			
			currPath = splicePaths[k]
			pathIDs = splicePathIDs[k]
			results = getMultiDeparturePoint(currPath, medial2, estPose2, estPose2, pathIDs, nodeID, pathPlotCount = self.multiDepCount, hypID = self.hypothesisID, plotIter = False)

			self.multiDepCount += 1

			resultSet.append(results+(k,))

			"departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departu rePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2"

		for result in resultSet:
			print "result:", results+(k,)
		
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
				sum1 = getOverlapCondition(self.poseData.medialAxes[nodeID], self.nodePoses[nodeID], self.trimmedPaths[pathID], nodeID, plotIter = False, overlapPlotCount = self.overlapPlotCount)
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
		self.junctions = {}
		self.terminals = {}
		self.allSplices = {}

	@logFunction
	def delPath(self, pathID, mergeTargetID):
		
		print "deleting path", pathID, "merging to path", mergeTargetID
		try: 
			del self.pathClasses[pathID]
			del self.pathTermsVisited[pathID]
			del self.theoryMedialLongPaths[pathID]
			del self.medialLongPaths[pathID]
			
			
			for k, pathClass in self.pathClasses.iteritems():
				
				if pathClass["parentID"] == pathID:
					pathClass["parentID"] = mergeTargetID
					
		except:
			pass
	
	@logFunction
	def addPath(self, parentID, branchNodeID, localJunctionPose):


		print "addPath(", parentID, branchNodeID, localJunctionPose
		
		nodePose = self.nodeRawPoses[branchNodeID]
		nodeOrigin = Pose(nodePose)
		globalJunctionPose = nodeOrigin.convertLocalOffsetToGlobal(localJunctionPose)
		
		
		oldPaths = self.getPathIDs()
		
		newPathID = self.pathIDs
		self.pathClasses[newPathID] = {"parentID" : parentID, "branchNodeID" : branchNodeID, "localJunctionPose" : localJunctionPose, 
							"sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : globalJunctionPose }		
				
		self.pathTermsVisited[newPathID] = False
		
		self.pathIDs += 1


		juncSamps = {}
		juncSamps["numParticles"] = self.numJuncSamples
		juncSamps["updateCount"] = 0
		juncSamps["snapshots"] = {0 : []}
		self.juncSamples[newPathID] = juncSamps



		updateCount = self.poseParticles["updateCount"] 
		particleDist2 = self.poseParticles["snapshots2"][0]


		pathSpline = SplineFit(self.paths[parentID])

		" add new path to particles "
		for part in particleDist2:

			origPose = self.nodePoses[branchNodeID]

			gpacProfile = Pose(origPose)
			localOffset = gpacProfile.convertGlobalPoseToLocal(self.nodeRawPoses[branchNodeID])
			
			" go back and convert this from GPAC pose to estPose "
			particlePose = deepcopy(part.pose0)
			particlePose[2] = origPose[2]
			newProfile = Pose(particlePose)
			newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
		


			#nodePose = part.pose0
			nodeOrigin = Pose(newEstPose)
			newGlobalJunctionPose = nodeOrigin.convertLocalOffsetToGlobal(localJunctionPose)

			minDist, uVal, newJuncPose = pathSpline.findClosestPoint(newGlobalJunctionPose)
			#minDist, uVal, splinePoint = pathSpline.findClosestPoint(globalJunctionPose)
			#arcDist = pathSpline.dist_u(uVal)

			#newJuncPose = pathSpline.getPointOfDist(newArcDist)
			modJuncPose = copy(newJuncPose)
			modJuncPose[2] = globalJunctionPose[2]

			#part.addPath(newPathID, parentID, branchNodeID, localJunctionPose, globalJunctionPose)
			part.addPath(newPathID, parentID, branchNodeID, localJunctionPose, modJuncPose)

		print "newPath", newPathID, "=", self.pathClasses[newPathID]
		#self.drawPoseParticles()

		self.isChanged = True
		self.junctions = {}
		self.terminals = {}
		self.allSplices = {}

		return newPathID
 
	@logFunction
	def getTopology(self, pathID):
		
		random.seed(0)		  
 
		self.theoryMedialLongPaths[pathID] = []
		self.medialLongPaths[pathID] = []
		self.longPathJunctions[pathID] = {}

		print "caseA"
		sys.stdout.flush()
		nodes = self.getNodes(pathID)


		junctionNodeID = self.pathClasses[pathID]["branchNodeID"]
		globalJunctionPoint = None

		if junctionNodeID != None:
		
			print "junctionNodeID:", junctionNodeID
			#print "nodes:", self.nodeHash.keys()
			print "nodes:", [k for k in range(self.poseData.numNodes)]
			
			globalJunctionPoint = self.getGlobalJunctionPose(pathID) 

		print "caseB"
		sys.stdout.flush()
		
		def convertAlphaUniform(a_vert, max_spacing = 0.04):
			
			" make the vertices uniformly distributed "
			
			new_vert = []
		
			for i in range(len(a_vert)):
				p0 = a_vert[i]
				p1 = a_vert[(i+1) % len(a_vert)]
				dist = math.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
	
				vec = [p1[0]-p0[0], p1[1]-p0[1]]
				
				if dist == 0:
					pass
				else:
										
					vec[0] /= dist
					vec[1] /= dist
					
					new_vert.append(copy(p0))
					
					if dist > max_spacing:
						" cut into pieces max_spacing length or less "
						numCount = int(math.floor(dist / max_spacing))
						
						for j in range(1, numCount+1):
							newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
							new_vert.append(newP)
			
			return new_vert			   

		medialPointSoup = []


		if nodes == []:
			return [], []


		print "caseC"
		sys.stdout.flush()

		if True:	
			for nodeID in nodes:

				#estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()		
				estPose1 = self.nodePoses[nodeID]
		
				#if self.nodeHash[nodeID].isBowtie:			  
				#	hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
				#else:
				#	hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
		
				hull1 = self.poseData.aHulls[nodeID]
		
				m = hashlib.md5()
				m.update(repr(hull1))
				print nodeID, "hull1 =", int(m.digest().encode('hex'),16)
		
				" set the origin of pose 1 "
				poseOrigin = Pose(estPose1)
		
				xP = []
				yP = []    
				for k in range(len(hull1)):
	 
					p = hull1[k]
					
					m = hashlib.md5()
					m.update(repr(p))
					
					p1 = poseOrigin.convertLocalToGlobal(p)
					
					m = hashlib.md5()
					m.update(repr(p1))
					
					medialPointSoup.append(p1)

		print "caseD"
		sys.stdout.flush()

		radius = 0.2

		numPoints = len(medialPointSoup)

		isDone = False

		print "caseE"
		sys.stdout.flush()
		
		while not isDone:
	
			perturbPoints = []
			
			for p in medialPointSoup:
				p2 = copy(p)
				
				" add a little bit of noise to avoid degenerate conditions in CGAL "
				" catch exception in case of math domain error "
				isReturned = False
				while not isReturned:
					try:
						p2[0] += random.gauss(0.0,0.000001)
					except:
						pass
					else:
						isReturned = True

				isReturned = False
				while not isReturned:
					try:
						p2[1] += random.gauss(0.0,0.000001)
					except:
						pass
					else:
						isReturned = True
						
	
				perturbPoints.append(p2)
		
			try:			
		
				saveFile = ""	 
				saveFile += "radius = " + repr(radius) + "\n"
				saveFile += "perturbPoints = " + repr(perturbPoints) + "\n"

				#print saveFile
				#sys.stdout.flush()
		
				isWritten = False
				while not isWritten:
					try:
						f = open("doAlphaInput_%08u.txt" % (self.alphaPlotCount), 'w')
						f.write(saveFile)
						f.close()
					except:
						pass
					else:
						isWritten = True
				print "caseF"
				sys.stdout.flush()
	
				vertices = alphamod.doAlpha(radius,perturbPoints)
				numVert = len(vertices)

				print "caseG"
				sys.stdout.flush()

				os.remove("doAlphaInput_%08u.txt" % (self.alphaPlotCount))
				self.alphaPlotCount += 1
				
				
				if numVert <= 2:
					print "Failed, hull had only", numVert, "vertices"
					raise
				
				isDone = True
			except:
				print "hull has holes!	retrying..."
				#print sArr		
		
		

		" cut out the repeat vertex "
		vertices = vertices[:-1]
		
		vertices = convertAlphaUniform(vertices)

		vertices.append(vertices[0])
		
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in vertices:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]
	
	
		" SPECIFY SIZE OF GRID AND CONVERSION PARAMETERS "
		PIXELSIZE = 0.05
		mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
		pixelSize = PIXELSIZE
		numPixel = int(2.0*mapSize / pixelSize + 1.0)
		divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)
		
		def realToGrid(point):
			indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
			indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
			return indexX, indexY
	
		def gridToReal(indices):
			i = indices[0]
			j = indices[1]
			pX = (i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0
			pY = (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0
			#point = ((i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0)
			point = (pX,pY)
			return point
	
		" CONVERT HULL TO GRID COORDINATES "
		gridHull = []
		for i in range(len(vertices)):
			p = vertices[i]
			gridHull.append(realToGrid(p))

		" BOUNDING BOX OF GRID HULL "	 
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in gridHull:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]

		" COMPUTE MEDIAL AXIS OF HULL "
		resultImg = Image.new('L', (numPixel,numPixel))
		resultImg = computeMedialAxis(self.medialCount, numPixel,numPixel, 5, resultImg, len(gridHull[:-2]), gridHull[:-2])
		self.medialCount += 1

		" EXTRACT POINTS FROM GRID TO LIST "
		imgA = resultImg.load()    
		points = []
		for i in range(1, numPixel-1):
			for j in range(1, numPixel-1):
				if imgA[i,j] == 0:
					points.append((i,j))
	
		" CREATE GRAPH NODES FOR EACH POINT "
		medialGraph = graph.graph()
		for p in points:
			medialGraph.add_node(p, [])
	
		" ADD EDGES BETWEEN NEIGHBORS "
		builtGraph = {}
		for i in range(2, numPixel-2):
			for j in range(2, numPixel-2):
				if imgA[i,j] == 0:
					builtGraph[(i,j)] = []
					for k in range(i-1, i+2):
						for l in range(j-1, j+2):
							if imgA[k,l] == 0 and not (k == i and l == j):
								builtGraph[(i,j)].append((k,l))
								medialGraph.add_edge((i,j), (k,l))
								
		" COMPUTE MINIMUM SPANNING TREE "
		mst = medialGraph.minimal_spanning_tree()
		
		" INITIALIZE DATA DICTS FOR UNIDIRECTIONAL MST"
		uni_mst = {}
		for k, v in mst.items():
			uni_mst[k] = []

		
		" ADD EDGES TO DICT TREE REPRESENTATION "
		for k, v in mst.items():
			if v != None:
				uni_mst[k].append(v)
				uni_mst[v].append(k)

		" LOCATE ALL LEAVES "
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
		
		" DELETE ALL NODES THAT ARE LEAVES, TO REMOVE SINGLE NODE BRANCHES "
		for leaf in leaves:
			medialGraph.del_node(leaf)	  
			
		" RECOMPUTE MST "
		mst = medialGraph.minimal_spanning_tree()

		" AGAIN, CREATE OUR DATA STRUCTURES AND IDENTIFY LEAVES "
		uni_mst = {}
		isVisited = {}
		nodeSum = {}
		for k, v in mst.items():
			uni_mst[k] = []
			isVisited[k] = 0
			nodeSum[k] = 0
		
		for k, v in mst.items():
			if v != None:
				uni_mst[k].append(v)
				uni_mst[v].append(k)
						
		" RECORD THE LEAVES AND JUNCTIONS "
		leaves = []
		junctions = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)

			if len(v) > 2:
				junctions.append(k)

		print "junctionNodeID:", junctionNodeID
		print "junctions:", junctions
		
		" SAVE FOR LATER, IDENTIFIED BY THEIR INDEX NOW "
		numLeafs = len(leaves)
		allLeafs = []
		allJunctions = []
		for leaf in leaves:
			allLeafs.append(gridToReal(leaf))

		if junctionNodeID != None:
		
			print "finding theoretical junction:"
			minKey = None
			minCand = None
			minJuncDist = 1e100
			for k, v in uni_mst.items():
				gCand = gridToReal(k)
				juncDist = sqrt((globalJunctionPoint[0]-gCand[0])**2 + (globalJunctionPoint[1]-gCand[1])**2)
				
				if juncDist < minJuncDist:
					minKey = k
					minJuncDist = juncDist
					minCand = gCand
		
			theoryVec = [cos(globalJunctionPoint[2]), sin(globalJunctionPoint[2])]

			theoryJunc = (minKey, minCand, minJuncDist, theoryVec)
			print "theoryJunc:", theoryJunc
			
			theoryJuncPoint = theoryJunc[1]
			theoryLeaf = []
			leafMag = 0.02
			for i in range(5):
				newPoint = (theoryVec[0]*leafMag + theoryJuncPoint[0], theoryVec[1]*leafMag + theoryJuncPoint[1])
				theoryLeaf.append(newPoint)
				leafMag += 0.02


			print "theoryLeaf:", theoryLeaf
			

			"""
			If there is a node of degree > 2, use these nodes as the junction points.
			If there is no nodes of degree > 2, find the closest node to the theoretical junction point and use this as the junction index 
			"""


			if len(junctions) > 0:
				
				for junc in junctions:
					gJunc = gridToReal(junc)
					juncDist = sqrt((globalJunctionPoint[0]-gJunc[0])**2 + (globalJunctionPoint[1]-gJunc[1])**2)
					allJunctions.append((junc, gJunc,juncDist))

			else:
				
				print "adding theoretical junction:", (minKey, minCand, minJuncDist)
				allJunctions.append(theoryJunc)						  
		
		print "allJunctions:", allJunctions
		
		" FIND ALL THE PATHS BETWEEN LEAVES "		 
		nodePaths = {}
		for leaf in leaves:
			
			" INITIALIZE DATA DICTS FOR VISITED, PATH SO FAR, AND NODE HOP COUNT"
			isVisited = {}
			nodeSum = {}
			nodePath = {}
			for k, v in uni_mst.items():
				isVisited[k] = 0
				nodeSum[k] = 0
				nodePath[k] = []
	
			" PERFORM DFS FOR THIS LEAF "
			getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
	
			" SAVE THE RESULTING PATH DATA STRUCTURE "
			nodePaths[leaf] = nodePath

		#leafSegments = []
		#internalSegments = []
		#isVisited = {}
		#nodeSum = {}
		#nodePath = {}
		leafSegments, internalSegments = computePathSegments(junctions, leaves, uni_mst)

		print "computePathSegments:", len(leafSegments), len(internalSegments), "paths from", len(junctions), "junctions and", len(leaves), "leaves", [len(pMem) for pMem in leafSegments], [len(pMem) for pMem in internalSegments]


		"""
		1) convert to real
		2) get theory leaf as well
		3) extrapolate the leaf segments to the boundary of the hull
		4) spline smooth
		5) convert to oriented points
		6) return set of points based on return condition
		"""

		realLeafSegments = []
		realInternalSegments = []
		for seg in leafSegments:
			newSeg = []
			for p in seg:
				newSeg.append(gridToReal(p))

			realLeafSegments.append(newSeg)

		for seg in internalSegments:
			newSeg = []
			for p in seg:
				newSeg.append(gridToReal(p))
			realInternalSegments.append(newSeg)

		longLeafSegments = []
		for seg in realLeafSegments:

			backVec = [0.,0.]
			indic = range(3)
			indic.reverse()
			
			for i in indic:
				if i+2 < len(seg):
					p1 = seg[-i-3]
					p2 = seg[-i-1]
					vec = [p2[0]-p1[0], p2[1]-p1[1]]
					backVec[0] += vec[0]
					backVec[1] += vec[1]
		
			backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
			backVec[0] /= backMag
			backVec[1] /= backMag
		
			newP2 = (seg[-1][0] + backVec[0]*10, seg[-1][1] + backVec[1]*10)
		
			realSeg = deepcopy(seg)

			realSeg.append(newP2)
			
	
			" take the long length segments at tips of medial axis"
			edge2 = realSeg[-2:]
			
			backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
			backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
			
			backVec[0] /= backMag
			backVec[1] /= backMag
			
			" make a smaller version of these edges "
			newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)
	
			edge2 = [edge2[0], newP2]

			" find the intersection points with the hull "
			hull = vertices
			interPoints = []
			for k in range(len(hull)-1):
				hullEdge = [hull[k],hull[k+1]]
				isIntersect2, point2 = Intersect(edge2, hullEdge)
				if isIntersect2:
					interPoints.append(point2)
					break
			
			" replace the extended edges with a termination point at the hull edge "			
			realSeg = realSeg[:-2]
			
			if isIntersect2:
				realSeg.append(point2)
	
			longLeafSegments.append(realSeg)


		smoothLeafSegments = []
		smoothInternalSegments = []

		for seg in longLeafSegments:
			leafSpline = SplineFit(seg, smooth=0.1)
			leafPoints = leafSpline.getUniformSamples()
			smoothLeafSegments.append(leafPoints)

		for seg in realInternalSegments:
			internalSpline = SplineFit(seg, smooth=0.1)
			internalPoints = internalSpline.getUniformSamples()
			smoothInternalSegments.append(internalPoints)


		

		" FOR EVERY PAIR OF LEAVES, SAVE ITS PATH IF ITS LONG ENOUGH "
		" should have X choose 2 combinations"
		longPaths = []
		
		MAX_LEN = 2
		for leaf1 in leaves:
			for leaf2 in leaves:
				if leaf1 < leaf2:
					if len(nodePaths[leaf1][leaf2]) > MAX_LEN:
						nPath = deepcopy(nodePaths[leaf1][leaf2])
						juncIndices = []
						for junc in allJunctions:
							try:
								index = nPath.index(junc[0])
							except:
								juncIndices.append(None)
							else:
								juncIndices.append(index)
								
						longPaths.append((len(nPath),nPath, juncIndices))


		if junctionNodeID != None:
			#theoryMedialLongPaths = []
			theoryPaths = []
			print "theoryPaths:"
			for leaf1 in leaves:
				nPath = deepcopy(nodePaths[leaf1][theoryJunc[0]])					 
				theoryPaths.append((len(nPath), nPath))
				print leaf1, theoryJunc, len(nPath)
	
			theoryPaths.sort(reverse=True)
	
			for k in range(len(theoryPaths)):
				path = theoryPaths[k][1]
				realPath = []
				for p in path:
					realPath.append(gridToReal(p))
				#realPath.reverse()
				
				
				theoryPaths[k] = realPath + theoryLeaf
				
				print "theoryPath(" , k , ")", len(realPath), realPath[-1], theoryLeaf
				
				theoryPath = theoryPaths[k]
				print "theoryPath:", len(theoryPath)
				
				leafPath = deepcopy(theoryPath)
				
				frontVec = [0.,0.]
				backVec = [0.,0.]
				indic = range(3)
				indic.reverse()
				
				for i in indic:
					if i+2 < len(leafPath):
						p1 = leafPath[i+2]
						p2 = leafPath[i]
						vec = [p2[0]-p1[0], p2[1]-p1[1]]
						frontVec[0] += vec[0]
						frontVec[1] += vec[1]
			
						p1 = leafPath[-i-3]
						p2 = leafPath[-i-1]
						vec = [p2[0]-p1[0], p2[1]-p1[1]]
						backVec[0] += vec[0]
						backVec[1] += vec[1]
			
				frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
				backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
			
				frontVec[0] /= frontMag
				frontVec[1] /= frontMag
				backVec[0] /= backMag
				backVec[1] /= backMag
			
				newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
				newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
			
				leafPath.insert(0,newP1)
				leafPath.append(newP2)
				
				medial2 = deepcopy(leafPath)
		
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
				hull = vertices
				interPoints = []
				for k in range(len(hull)-1):
					hullEdge = [hull[k],hull[k+1]]
					isIntersect1, point1 = Intersect(edge1, hullEdge)
					if isIntersect1:
						interPoints.append(point1)
						break
		
				for k in range(len(hull)-1):
					hullEdge = [hull[k],hull[k+1]]
					isIntersect2, point2 = Intersect(edge2, hullEdge)
					if isIntersect2:
						interPoints.append(point2)
						break
				
				" replace the extended edges with a termination point at the hull edge "			
				medial2 = medial2[1:-1]
				
				if isIntersect1:
					medial2.insert(0, point1)
				if isIntersect2:
					medial2.append(point2)
	
				self.theoryMedialLongPaths[pathID].append(deepcopy(medial2))








		print "longPaths:"
		for longPath in longPaths:
			print longPath[2]
		
		" SORT FOR THE LONGEST TO SHORTEST "
		longPaths.sort(reverse=True)
		
		" REMOVE SIZE FROM TUPLE  "
		juncGrids = []
		juncIndices = []
		longPathLengths = []
		for k in range(len(longPaths)):
			juncIndices.append(longPaths[k][2])

			juncInds = longPaths[k][2]
			jGrids = []
			for index in juncInds:
				if index != None:
					jGrids.append(longPaths[k][1][index])
			juncGrids.append(jGrids)

			longPathLengths.append(longPaths[k][0])
			longPaths[k] = longPaths[k][1]
		
		print "juncIndices:", juncIndices
		
		" GET THE LEAF INDEXES TO EACH PATH "
		leafPairs = []
		for path in longPaths:
			leaf1 = path[0]
			leaf2 = path[-1]	
			
			leafID1 = leaves.index(leaf1)
			leafID2 = leaves.index(leaf2)

			leafPairs.append((leafID1,leafID2))
		
		print "leafPairs:", leafPairs

		juncReals = []
		for k in range(len(longPaths)):
			path = longPaths[k]
			realPath = []
			for p in path:
				realPath.append(gridToReal(p))
			
			longPaths[k] = realPath
			#juncReals.append(longPaths[k][juncIndices[k]])

			juncInds = juncIndices[k]
			jReals = []
			for index in juncInds:
				if index != None:
					jReals.append(copy(longPaths[k][index]))
			juncReals.append(jReals)
					

		#self.medialLongPaths = []
		juncAngSet = []

		print len(longPaths), "long paths"
		for n in range(len(longPaths)):
			
			longPath = longPaths[n]
			print "longPath:", len(longPath)
			print longPath
			
			leafPath = deepcopy(longPath)
			
			frontVec = [0.,0.]
			backVec = [0.,0.]
			indic = range(3)
			indic.reverse()
			
			for i in indic:
				if i+2 < len(leafPath):
					p1 = leafPath[i+2]
					p2 = leafPath[i]
					vec = [p2[0]-p1[0], p2[1]-p1[1]]
					frontVec[0] += vec[0]
					frontVec[1] += vec[1]
		
					p1 = leafPath[-i-3]
					p2 = leafPath[-i-1]
					vec = [p2[0]-p1[0], p2[1]-p1[1]]
					backVec[0] += vec[0]
					backVec[1] += vec[1]
		
			frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
			backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
			frontVec[0] /= frontMag
			frontVec[1] /= frontMag
			backVec[0] /= backMag
			backVec[1] /= backMag
		
			newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
			newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
		
			leafPath.insert(0,newP1)
			leafPath.append(newP2)
			print "leafPath:", leafPath
			
			medial2 = deepcopy(leafPath)
	
			" take the long length segments at tips of medial axis"
			edge1 = medial2[0:2]
			edge2 = medial2[-2:]
			print "edge1, edge2 =", edge1, edge2
			
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

			print "edge1, edge2 =", edge1, edge2

			" find the intersection points with the hull "
			hull = vertices
			interPoints = []
			for k in range(len(hull)-1):
				hullEdge = [hull[k],hull[k+1]]
				isIntersect1, point1 = Intersect(edge1, hullEdge)
				if isIntersect1:
					interPoints.append(point1)
					break
	
			for k in range(len(hull)-1):
				hullEdge = [hull[k],hull[k+1]]
				isIntersect2, point2 = Intersect(edge2, hullEdge)
				if isIntersect2:
					interPoints.append(point2)
					break
			
			print isIntersect1, isIntersect2, interPoints
			" replace the extended edges with a termination point at the hull edge "			
			medial2 = medial2[1:-1]

			
			if isIntersect1:
				medial2.insert(0, point1)
			if isIntersect2:
				medial2.append(point2)

			print "medial2:", medial2

			self.medialLongPaths[pathID].append(deepcopy(medial2))

			print "globalJunctionPoint:", globalJunctionPoint
			print "juncIndices:", juncIndices[n]

			
			jReals = juncReals[n]
			mLongPath = self.medialLongPaths[pathID][n]
			jIndices = []

			for jPoint in jReals:

				try:
					index = mLongPath.index(jPoint)
				except:
					jIndices.append(None)
				else:
					jIndices.append(index)
					
			print "jIndices:", jIndices


			juncAngs = []
			if globalJunctionPoint != None:
				for juncInd in jIndices:
					if juncInd != None:

						frontVec = [0.,0.]
						backVec = [0.,0.]
						indic = range(3)
						indic.reverse()
						
						print "len(mLongPath):", len(mLongPath)
						print "juncInd:", juncInd

						highIndex = juncInd+4
						highMod = highIndex
						if highIndex+4 >= len(mLongPath):
							highMod = len(mLongPath) - 5
						
						lowIndex = juncInd-4
						lowMod = lowIndex
						if lowIndex-4 < 0:
							lowMod = 4
						
						print "highIndex, highMod:", highIndex, highMod
						print "lowIndex, lowMod:", lowIndex, lowMod
						
						
						for i in indic:
							p1 = mLongPath[i+highMod]
							p2 = mLongPath[i+2+highMod]
							vec = [p2[0]-p1[0], p2[1]-p1[1]]
							frontVec[0] += vec[0]
							frontVec[1] += vec[1]
					
							p1 = mLongPath[-i+lowMod]
							p2 = mLongPath[-i+lowMod-2]
							vec = [p2[0]-p1[0], p2[1]-p1[1]]
							backVec[0] += vec[0]
							backVec[1] += vec[1]
						
						frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
						backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
					
						frontVec[0] /= frontMag
						frontVec[1] /= frontMag
						backVec[0] /= backMag
						backVec[1] /= backMag
						
						print "frontVec:", frontVec
						print "backVec:", backVec
						
						foreAng = acos(frontVec[0])
						if frontVec[1] < 0.0:
							foreAng = -foreAng
		
						backAng = acos(backVec[0])
						if backVec[1] < 0.0:
							backAng = -backAng
		
						print "foreAng:", foreAng
						print "backAng:", backAng
		
						frontError = normalizeAngle(globalJunctionPoint[2]-foreAng)
						backError = normalizeAngle(globalJunctionPoint[2]-backAng)
						juncAngs.append((frontError,backError))
					else:
						juncAngs.append(None)
			juncAngSet.append(juncAngs)


		
		juncDists = []
		for junc in allJunctions:
			juncDists.append(junc[2])

		juncLongIndices = []
		for k in range(len(self.medialLongPaths[pathID])):
			
			jReals = juncReals[k]
			mLongPath = self.medialLongPaths[pathID][k]

			jIndices = []

			for jPoint in jReals:
				try:
					index = mLongPath.index(jPoint)
				except:
					print "index failure!"
					print "jPoint:", jPoint
					print "jReals:", jReals
					print "mLongPath:", mLongPath
					print "juncReals:", juncReals
					print "juncAngSet:", juncAngSet
					print "allJunctions:", allJunctions
					print "juncIndices:", juncIndices
					print "juncDists:", juncDists
					raise
				else:
					jIndices.append(index)

				"""
				try:
				except:
					jIndices.append(None)
				else:
				"""
					
			juncLongIndices.append(jIndices)

		juncMedialReals = []
		for k in range(len(self.medialLongPaths[pathID])):
			juncInds = juncLongIndices[k]
			jMedialReals = []
			for index in juncInds:
				jMedialReals.append(copy(self.medialLongPaths[pathID][k][index]))
			juncMedialReals.append(jMedialReals)

		self.longPathJunctions[pathID]["juncIndices"] = juncIndices
		self.longPathJunctions[pathID]["juncAngSet"] = juncAngSet
		self.longPathJunctions[pathID]["juncDists"] = juncDists
		self.longPathJunctions[pathID]["longPathLengths"] = longPathLengths
		self.longPathJunctions[pathID]["juncMedialReals"] = juncMedialReals
		self.longPathJunctions[pathID]["juncReals"] = juncReals
		self.longPathJunctions[pathID]["juncGrids"] = juncGrids


		self.longPathJunctions[pathID]["juncLongIndices"] = juncLongIndices

		self.longPathJunctions[pathID]["leafSegments"] = smoothLeafSegments
		self.longPathJunctions[pathID]["internalSegments"] = smoothInternalSegments


		juncPoints = []
		juncArmPoints = []
		juncDesc = {}
		for k in range(len(juncLongIndices)):
			juncInds = juncLongIndices[k]
			juncPnts = []
			juncArmPnts = []
			for index in juncInds:
				jPnt = copy(self.medialLongPaths[pathID][k][index])

				jaPnt1 = self.medialLongPaths[pathID][k][index+1]
				jaPnt2 = self.medialLongPaths[pathID][k][index-1]

				jPnt = tuple(jPnt)
				jaPnt1 = tuple(jaPnt1)
				jaPnt2 = tuple(jaPnt2)

				juncPnts.append(jPnt)
				juncArmPnts.append((jaPnt1, jaPnt2))

				try:
					juncDesc[jPnt].append(jaPnt1)
					juncDesc[jPnt].append(jaPnt2)
				except:
					juncDesc[jPnt] = []
					juncDesc[jPnt].append(jaPnt1)
					juncDesc[jPnt].append(jaPnt2)


			juncPoints.append(juncPnts)
			juncArmPoints.append(juncArmPnts)

		for k, v in juncDesc.iteritems():
			v1 = set(v)
			juncDesc[k] = list(v1)

					

		self.longPathJunctions[pathID]["juncPoints"] = juncPoints
		self.longPathJunctions[pathID]["juncArmPoints"] = juncArmPoints
		self.longPathJunctions[pathID]["juncDesc"] = juncDesc

		
		if False:
			pylab.clf()
	
			for path in self.medialLongPaths[pathID]:
				xP = []
				yP = []
				for p in path:
					xP.append(p[0])
					yP.append(p[1])
	
				pylab.plot(xP,yP)

	
			if junctionNodeID != None:
	
				for path in self.theoryMedialLongPaths[pathID]:
					xP = []
					yP = []
					for p in path:
						xP.append(p[0])
						yP.append(p[1])
	
					pylab.plot(xP,yP, color='k')

			globJuncPose = self.getGlobalJunctionPose(pathID)
			if globJuncPose != None:
				pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='k')

			
			xP = []
			yP = []
			for p in vertices:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			sizes = []
			for path in longPaths:
				sizes.append(len(path))
			
			bufStr1 = ""
			for dist in juncDists:
				if dist != None:
					bufStr1 += "%1.2f " % dist


			bufStr2 = ""
			for juncAngs in juncAngSet:
				for angs in juncAngs:
					if angs != None:
						bufStr2 += "%1.2f %1.2f " % (angs[0],angs[1])
			
			pylab.axis("equal")
			pylab.title("Path %d %s %s %s" % (pathID, sizes,bufStr1,bufStr2))
			pylab.savefig("medialOut2_%04u.png" % self.topCount)
			print "saving medialOut2_%04u.png" % self.topCount
			self.topCount += 1


		
		print "juncAngSet:", juncAngSet

		
		" searches for the branch from the parent junction "
		" selects splice of topology that has aligning junction "
		" does not select for distance or best medial axis representation of path "
		

		" sort by longest first "
		pathCands = []
		for k in range(len(longPaths)):
			pathCands.append((len(longPaths[k]), k))

		pathCands.sort(reverse=True)
		
		print "pathCands:", pathCands
		

		" TODO: find longest path going to right, longest path going to left through the junction "


		
		" FIXME:  Verify that leaf path exists for the junction point "
		
		" select longest path that has the best angle fit "
		branchArm = None
		if globalJunctionPoint != None:
			bestFit = -1
			minDiff = 1e100
			for cand in pathCands:
				k = cand[1]
				juncAngs = juncAngSet[k]
				jlIndices = juncLongIndices[k]
				
				#for angs in juncAngs:
				for l in range(len(juncAngs)):

					angs = juncAngs[l]
					jIndex = jlIndices[l]

					if angs != None:
						angDiff1 = angs[0]
						angDiff2 = angs[1]
						
						print "minDiff:", minDiff, angDiff1, angDiff2 
						
						if fabs(angDiff1) < minDiff:
							minDiff = fabs(angDiff1)
							bestFit = k
							branchArm = self.medialLongPaths[pathID][bestFit][jIndex+1]
	
						if fabs(angDiff2) < minDiff:
							minDiff = fabs(angDiff2)
							bestFit = k
							branchArm = self.medialLongPaths[pathID][bestFit][jIndex-1]

				pass



			if bestFit != -1:
				
				if fabs(minDiff) < 1.047:

					print "returning bestFit:", bestFit, minDiff
					self.longPathJunctions[pathID]["bestFit"] = bestFit
					self.longPathJunctions[pathID]["branchArm"] = branchArm
					return self.medialLongPaths[pathID][bestFit], vertices

				else:
					" FIXME:  return theory junction information if this is selected "
					print "returning bestFit theory:", theoryJunc
					self.longPathJunctions[pathID]["bestFit"] = 0
					self.longPathJunctions[pathID]["branchArm"] = branchArm
					return self.theoryMedialLongPaths[pathID][0], vertices
			
			else:
				print "not returning bestFit"

		print "returning longest fit"

		maxIndex = 0
		maxLen = 0
		for k in range(len(longPathLengths)):
			if longPathLengths[k] > maxLen:
				maxIndex = k
				maxLen = longPathLengths[k]

		self.longPathJunctions[pathID]["bestFit"] = maxIndex


		" FIXME: change from selecting any old junction point since "
		#jIndex = juncLongIndices[maxIndex][0]
		self.longPathJunctions[pathID]["branchArm"] = self.medialLongPaths[pathID][maxIndex][1]

		return self.medialLongPaths[pathID][maxIndex], vertices

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

		for pathID in pathIDs:
			
			if self.pathClasses[pathID]["parentID"] == None:
				trimmedPaths[pathID] = deepcopy(paths[pathID])

			else:
				parentPathID = self.pathClasses[pathID]["parentID"]
				childPathID = pathID

				path1 = paths[parentPathID]
				path2 = paths[childPathID]
				


				#secP1, secP2 = getOverlapDeparture(parentPathID, childPathID, paths, plotIter = False)				 

				branchNodeID = self.pathClasses[childPathID]["branchNodeID"]
				branchNodePose = self.nodePoses[branchNodeID]
				globalJunctionPoint = self.getGlobalJunctionPose(childPathID)

				#secP1, secP2 = getOverlapDeparture(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False)				 
				secP1, secP2 = getOverlapDeparture(globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False)				 

				minI_1 = 0		  
				minI_2 = 0
				minDist_1 = 1e100
				minDist_2 = 1e100		 
				for i in range(len(path2)):
					pnt = path2[i]
					dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
					dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
				
					if dist1 < minDist_1:
						minDist_1 = dist1
						minI_1 = i

					if dist2 < minDist_2:
						minDist_2 = dist2
						minI_2 = i

				term0 = path2[0]
				termN = path2[-1]

				" smallest distance is the terminal point "
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
				
				junctionPoint = self.getGlobalJunctionPose(pathID)
				
				if minDist_1 < minDist_2:
					junctionPoint_K = secP1
					juncI_K = minI_1
				else:
					junctionPoint_K = secP2
					juncI_K = minI_2
				
				minDist2 = 1e100
				juncI = 0		 
				for i in range(len(path2)):
					pnt = path2[i]
					dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
				
					if dist < minDist2:
						minDist2 = dist
						juncI = i
				 
				 
				
				print "len(path1):", len(path1)
				print "len(path2):", len(path2)
				print "juncI:", juncI
				print "minDist:", minDist_1, minDist_2
				
				" now we have closest point to departure point. "
				" Which side is the departing side? "	 

						
				if minI == 0:
					"secP1 is terminal 0"
					index = juncI-10
					if index < 1:
						index = 1

					
					newPath2 = path2[:index+1]
					newPath2.reverse()
				
				elif minI == 1:
					"secP1 is terminal N"
					index = juncI+10
					if index >= len(path2)-1:
						
						" ensure at least 2 elements in path "
						index = len(path2)-2

					newPath2 = path2[index:]

					
				elif minI == 2:
					"secP2 is terminal 0"
					index = juncI-10
					if index < 1:
						index = 1

					newPath2 = path2[:index+1]
					newPath2.reverse()
					
				elif minI == 3:
					"secP2 is terminal N"
					index = juncI+10
					if index >= len(path2)-1:
						" ensure at least 2 elements in path "
						index = len(path2)-2
					
					newPath2 = path2[index:]
				
				else:
					print "no terminal found"
					raise
								
				" convert path so that the points are uniformly distributed "
				
				
				max_spacing = 0.08
				newPath3 = []

				" make sure path is greater than 5 points "
				while len(newPath3) <= 5:
	
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
							" cut into pieces max_spacing length or less "
							numCount = int(floor(dist / max_spacing))
							
							for j in range(1, numCount+1):
								newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
								newPath3.append(newP)
	
						newPath3.append(copy(p1))			 
				
				trimmedPaths[pathID] = deepcopy(newPath3)

		self.trimmedPaths = trimmedPaths
		
		return trimmedPaths


	@logFunction
	def getBranchPoint(self, globalJunctionPoint, parentPathID, childPathID, path1, path2, plotIter = False):
		""" get the trimmed version of child and parent paths that are overlapping in some fashion """

		"Assumption:  one section of the medial axis is closely aligned with the path "		   
		
		print "getBranchPoint():"
		
		" return exception if we receive an invalid path "		  
		if len(path1) == 0:
			print "path1 has zero length"
			raise
		
		" make sure the overlap of both paths are oriented the same way "
		orientedPath2 = orientPath(path2, path1)
		
		path1Spline = SplineFit(path1, smooth=0.1)			  
		path2Spline = SplineFit(orientedPath2, smooth=0.1)


		" for each point on the child path, find its closest pair on the parent path "
		
		pathPoints1 = path1Spline.getUniformSamples(interpAngle=True)
		pathPoints2 = path2Spline.getUniformSamples(interpAngle=True)

		distances = []
		indices = []
		juncDists = []

		minDist2 = 1e100
		juncI = 0		 

		for i in range(0,len(pathPoints2)):
			p_2 = pathPoints2[i]
			p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
			
			" keep the distance information "
			distances.append(minDist)
			" and the associated index of the point on the parent path "
			indices.append(i_1)
			
			juncDist = sqrt((p_2[0]-globalJunctionPoint[0])**2 + (p_2[1]-globalJunctionPoint[1])**2)

			if juncDist < minDist2:
				minDist2 = juncDist
				juncI = i

			juncDists.append(juncDist)


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
				backFound = k
				backFound = True

		print "frontFound, backFound:", frontFound, backFound
		print "frontInd, backInd:", frontInd, backInd

		print "minDist2, juncI:", minDist2, juncI

		" match distances of the tip points of child path "		   
		maxFront = distances[frontInd]
		maxBack = distances[backInd]
		print "match distances of tip points:", maxFront, maxBack

		" walk back from tip point until we have a non-monotic increase in match distance "
		" this becomes our departure point "
		
		"TODO:	why is there a 3-offset in this comparison? "
		currI = frontInd + 1
		try:
			while distances[currI+3] < maxFront or distances[currI+3] > 0.1:
				maxFront = distances[currI]
				currI += 1
		except:
			pass
		
		" departure index on child path "
		frontDepI = currI

		if frontDepI+3 >= len(pathPoints2)-1:
			frontDepI = juncI

		
		" departure index of child path and match distance "
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


		" FIXME:  index out of bounds case "
		currI = 1 + len(pathPoints2) - backInd
		try:
			while distances[-currI-3] < maxBack or distances[-currI-3] > 0.1:
				maxBack = distances[-currI]
				currI += 1
		except:
			pass

		" departure index on child path "
		backDepI = len(distances) - currI

		if backDepI-3 <= 0: 
			backDepI = juncI

		" departure index of child path and match distance "
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


		print "frontDepI, backDepI:", frontDepI, backDepI
		print "frontAngleRefI, backAngleRefI:", frontAngleRefI, backAngleRefI
		print "forePathAngle, backPathAngle:", forePathAngle, backPathAngle
		print "newFrontDepI, newBackDepI:", newFrontDepI, newBackDepI
		print "foreDiff, backDiff:", diffAngle(normalizeAngle(pathPoints2[newFrontDepI][2]+pi), forePathAngle), diffAngle(normalizeAngle(pathPoints2[newBackDepI][2]), backPathAngle)

		print "frontDiffAngles:", frontDiffAngles
		print "backDiffAngles:", backDiffAngles


		print "lengths of parent and child paths:", len(pathPoints1), len(pathPoints2)
		print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint


		"reset to the tip match distance "
		maxFront = distances[frontInd]
		maxBack = distances[backInd]
		
		
		" sum of closest points on front and back "
		" select the one with minimal cost "		

		" path section for our front departure hypothesis "
		pathSec1 = pathPoints2[:frontDepI+1]
		pathSec1.reverse()

		" path section for our back departure hypothesis "
		pathSec2 = pathPoints2[backDepI:]
		
		" distance of departure point from known junction point "
		#juncDist1 = -1
		#juncDist2 = -1
		#if len(pathSec1) > 0:
		p0 = pathSec1[0]
		juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

		#if len(pathSec2) > 0:
		p0 = pathSec2[0]
		juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

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


		junctionPoint = globalJunctionPoint		 
		minDist2 = 1e100
		juncI = 0		 
		for i in range(len(pathPoints2)):
			pnt = pathPoints2[i]
			dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
		
			if dist < minDist2:
				minDist2 = dist
				juncI = i
		
		
		
		if juncDist1 < juncDist2:
			" FIXME: if [0,frontDepI] has high overlap compared to total length, then we reject it "
			secP1 = pathPoints2[0]
			secP2 = pathPoints2[frontDepI]
			
			newPath2 = pathPoints2[:newFrontDepI]
			newPath2.reverse()

		else:
			" FIXME: if [backDepI,len(distances)-1] has high overlap compared to total length, then we reject it "
			secP1 = pathPoints2[backDepI]
			secP2 = pathPoints2[len(distances)-1]

			newPath2 = pathPoints2[newBackDepI:]

		if len(secP1) == 0:
			print "no departures found"
			raise
		
		" convert path so that the points are uniformly distributed "
		
		
		max_spacing = 0.08
		newPath3 = []

		" make sure path is greater than 5 points "
		while len(newPath3) <= 5:

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
					" cut into pieces max_spacing length or less "
					numCount = int(floor(dist / max_spacing))
					
					for j in range(1, numCount+1):
						newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
						newPath3.append(newP)

				newPath3.append(copy(p1))			 
		


		leafPath = deepcopy(newPath3)
		
		frontVec = [0.,0.]
		backVec = [0.,0.]
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

	
		newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
	
		leafPath.insert(0,newP1)
		
		medial2 = deepcopy(leafPath)

		" take the long length segments at tips of medial axis"
		edge1 = medial2[0:2]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)

		edge1 = [newP1, edge1[1]]

		" find the intersection points with the hull "
		interPoints = []
		for k in range(len(path1)-1):
			hullEdge = [path1[k],path1[k+1]]
			isIntersect1, point1 = Intersect(edge1, hullEdge)
			if isIntersect1:
				interPoints.append(point1)
				break

		
		" replace the extended edges with a termination point at the hull edge "			
		medial2 = medial2[1:]
		
		if isIntersect1:
			medial2.insert(0, point1)


		juncAng = acos(-frontVec[0])
		if -frontVec[1] < 0.0:
			juncAng = -juncAng
		
		

		try:
			foreIntI, backIntI, juncForeAng, juncBackAng = getTangentIntersections(pathPoints1, pathPoints2, frontDepI, backDepI, indices[frontDepI], indices[backDepI], juncI, indices[juncI], self.pathPlotCount)
			print "foreIntI, backIntII:", foreIntI, backIntI

			""" junction distances are equal if both of the indices selected juncI as the departure point on path2'
				occurs if path2 does not come close enough to path1 to get under distance 0.1
			"""
			print "juncDist1, juncDist2 =", juncDist1, juncDist2
			if juncDist1 == juncDist2:

				juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

				foreDist = sqrt((pathPoints1[foreIntI][0]-globalJunctionPoint[0])**2 +  (pathPoints1[foreIntI][1]-globalJunctionPoint[1])**2)
				backDist = sqrt((pathPoints1[backIntI][0]-globalJunctionPoint[0])**2 +  (pathPoints1[backIntI][1]-globalJunctionPoint[1])**2)

				print "foreDist, backDist =", foreDist, backDist

				if foreDist < backDist:
					globJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], juncForeAng]
				else:
					globJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], juncBackAng]
		
			elif juncDist1 < juncDist2:
				globJuncPose = [pathPoints1[foreIntI][0], pathPoints1[foreIntI][1], juncForeAng]
			else:
				globJuncPose = [pathPoints1[backIntI][0], pathPoints1[backIntI][1], juncBackAng]

			print "globJuncPose =", globJuncPose

		except:
			print "getTangentIntersections() failed!"
			globJuncPose = [medial2[0][0], medial2[0][1], juncAng]
		
		" last junction point is the intersection point "

		if plotIter:
			pylab.clf()
			xP = []
			yP = []
			for p in path2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(0.5,0.5,1.0))
	
			if True:	
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

	
				
			pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')		   
			pylab.scatter([globJuncPose[0]],[globJuncPose[1]], color='k')		 
	
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

			print "trimDeparture:", self.pathPlotCount
			printStack()				 

			pylab.title("%3.2f %3.2f %3.2f" % (juncDist1, juncDist2, juncAng))
			pylab.savefig("trimDeparture_%04u.png" % self.pathPlotCount)
			
			self.pathPlotCount += 1


		
		return globJuncPose

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
		estPose2 = self.nodePoses[nodeID]
			
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
		estPose2 = self.nodePoses[nodeID]
		
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
				pylab.plot(xP,yP)
			
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.title("%d %d cost = %f, count = %d" % (pathID1, pathID2, cost, len(support_pairs)))
			pylab.savefig("pathOverlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 0.0

		return cost


	@logFunction
	def getOverlapCondition(self, supportLine, nodeID, plotIter = False):

		if len(supportLine) == 0:
			return 1e100

		#node2 = self.nodeHash[nodeID]

		#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
		hull2 = self.poseData.aHulls[nodeID]
		medial2 = self.poseData.medialAxes[nodeID]

		#estPose2 = node2.getGlobalGPACPose()		
		estPose2 = self.nodePoses[nodeID]
		
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
			pylab.savefig("overlapCost_%04u.png" % self.overlapPlotCount)
			self.overlapPlotCount += 1
		
		if len(support_pairs) == 0:
			return 1e100

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
				if isNodeFeatureless:
					print "REJECT new branch because the branching node is featureless"
					return False, -1
		
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

			

