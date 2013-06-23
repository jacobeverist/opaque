import pylab
import numpy
import multiprocessing as processing
import ctypes, os
from functions import normalizeAngle, diffAngle
from copy import copy, deepcopy
from SplineFit import SplineFit
from Pose import Pose
from math import pi, sqrt, acos, fabs, cos, sin
import nelminICP
from icp import matchPairs

renderGlobalPlotCount = 0

def __num_processors():
	
	return 4
	
	if os.name == 'nt': # Windows
		return int(os.getenv('NUMBER_OF_PROCESSORS'))
	else: # glibc (Linux, *BSD, Apple)
		get_nprocs = ctypes.cdll.libc.get_nprocs
		get_nprocs.restype = ctypes.c_int
		get_nprocs.argtypes = []
		return get_nprocs()

def __remote_multiFit(rank, qin, qout, globalPath, medial, initPose, pathIDs, nodeID):

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		#print rank, "received", nc, args
		# process data
		#knn = __do_nothing(data, nc, someArg2, someArg3)
		results = []
		for arg in args:
			results.append(multiFitSplice(arg[0], globalPath, medial, initPose, pathIDs, nodeID, arg[1]))
						   
			#resultPose, lastCost, matchCount = globalOverlapICP_GPU2(arg, globalPath, medial)
			#results.append([resultPose, lastCost, matchCount])
		# write to output queue
		qout.put((nc,results))


def multiFitSplice(initGuess, orientedPath, medialAxis, initPose, pathIDs, nodeID, pathPlotCount = 0):
	
	u1 = initGuess[0]
	u2 = initGuess[1]
	angGuess = initGuess[2]
	
	#u2 = originU2
	#u1 = currPathU
	#angGuess = 0.0
	
	resultPose, lastCost, matchCount = globalOverlapICP_GPU2([u1,u2,angGuess], orientedPath, medialAxis, globalPlotCount = pathPlotCount, plotIter = False, n1 = nodeID, n2 = -1)
		
	#self.nodeHash[nodeID].setGPACPose(resultPose)
	
	#pathPlotCount = 0
	resultArgs = getMultiDeparturePoint(orientedPath, medialAxis, initPose, resultPose, pathIDs, nodeID, pathPlotCount, plotIter = True)

	isExist1 = resultArgs[3]


	return (resultPose, lastCost, matchCount) + resultArgs


def batchGlobalMultiFit(initGuesses, globalPath, medial, initPose, pathIDs, nodeID):

	global renderGlobalPlotCount


	#initGuess

	#arg = [initGuess, initPose, pathID, nodeID]

	ndata = len(initGuesses)
	
	args = []
	
	for k in range(ndata):
		arg = [initGuesses[k], k + renderGlobalPlotCount]
		args.append(arg)
	
	renderGlobalPlotCount += len(initGuesses)
	
	
	nproc = __num_processors()
	
	print "nproc =", nproc
	
	# compute chunk size
	#chunk_size = ceil(float(ndata) / float(nproc))
	#chunk_size = int(chunk_size)
	chunk_size = ndata / nproc
	chunk_size = 2 if chunk_size < 2 else chunk_size
	
	print "chunk_size =", chunk_size
	print "max_size =", ndata/chunk_size
	
	# set up a pool of processes
	qin = processing.Queue(maxsize=ndata/chunk_size)
	qout = processing.Queue(maxsize=ndata/chunk_size)
	pool = [processing.Process(target=__remote_multiFit,
				args=(rank, qin, qout, globalPath, medial, initPose, pathIDs, nodeID))
					for rank in range(nproc)]
	for p in pool: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	while 1:
		#_data = data[:,cur:cur+chunk_size]
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


		




def getMultiDeparturePoint(currPath, medial2, initPose2, estPose2, pathIDs, nodeID, pathPlotCount = 0, plotIter = False):
	
	isExist1 = False
	isInterior1 = False
	departurePoint1 = 0
	angle1 = 0.0
	isExist2 = False
	isInterior2 = False
	departurePoint2 = 0
	angle2 = 0.0
	
	if len(currPath) == 0:
		return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0, 0.0, 0.0, 0.0
	
	#node2 = self.nodeHash[nodeID]
	
	#hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
	
	#estPose2 = node2.getGlobalGPACPose()		
	
	"Assumption:  one section of the medial axis is closely aligned with the path "
	poseOrigin = Pose(estPose2)
	
	medialSpline2 = SplineFit(medial2, smooth=0.1)
	points2 = medialSpline2.getUniformSamples(interpAngle=True)

	points2_offset = []
	for p in points2:
		result = poseOrigin.convertLocalOffsetToGlobal(p)
		points2_offset.append(result)


	globalMedialSpline = SplineFit(points2_offset, smooth=0.1)
	globalMedialPoints = globalMedialSpline.getUniformSamples()



	pathSpline = SplineFit(currPath, smooth=0.1)
	pathPoints = pathSpline.getUniformSamples()

	currPathReverse = deepcopy(currPath)
	currPathReverse.reverse()
	pathSplineReverse = SplineFit(currPathReverse, smooth=0.1)
	pathPointsReverse = pathSplineReverse.getUniformSamples()

	overlapMatch = []
	angleSum1 = 0.0
	angleSum2 = 0.0
	for i in range(0,len(points2_offset)):
		p_1 = points2_offset[i]
		p_2, j, minDist = findClosestPointInA(pathPoints, p_1)
		if minDist < 0.5:
			overlapMatch.append((i,j,minDist))

			pathU1 = globalMedialSpline.findU(p_1)	
			pathU2 = pathSpline.findU(p_2)	
			pathU2_R = pathSplineReverse.findU(p_2)	

			pathVec1 = globalMedialSpline.getUVector(pathU1)
			pathVec2 = pathSpline.getUVector(pathU2)
			pathVec2_R = pathSplineReverse.getUVector(pathU2_R)

			val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
			if val1 > 1.0:
				val1 = 1.0
			elif val1 < -1.0:
				val1 = -1.0
			ang1 = acos(val1)
			
			val2 = pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1]
			if val2 > 1.0:
				val2 = 1.0
			elif val2 < -1.0:
				val2 = -1.0
			ang2 = acos(val2)

			angleSum1 += ang1
			angleSum2 += ang2

	
	" select global path orientation based on which has the smallest angle between tangent vectors "
	print "getDeparturePoint:", i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
	if angleSum1 > angleSum2:
		pathPoints = pathPointsReverse
		pathSpline = pathSplineReverse

	" tip angles "
	angSum1 = 0.0
	angSum2 = 0.0
	angs1 = []
	angs2 = []
	phi1 = normalizeAngle(points2_offset[0][2])
	phi2 = normalizeAngle(points2_offset[-1][2])
	for i in range(10):
		ang1 = normalizeAngle(points2_offset[i][2]-phi1)
		ang2 = normalizeAngle(points2_offset[-i-1][2]-phi2)
		angSum1 += ang1
		angSum2 += ang2
		
		angs1.append(ang1+phi1)
		angs2.append(ang2+phi2)

	angle1 = angSum1 / 10.0 + phi1
	angle2 = angSum2 / 10.0 + phi2

	" invert one angle so opposite tips have opposite angles "
	angle1 = normalizeAngle(angle1 + pi)

	print "ang1:", angle1, angs1
	print "ang2:", angle2, angs2
	print "diff:", diffAngle(angle1, angle2)
	
	distSum = 0.0
	contigCount = 0
	maxContig = 0
	distances = []
	indices = []
	for i in range(0,len(points2_offset)):
		p_2 = points2_offset[i]
		p_1, i_1, minDist = findClosestPointInA(pathPoints, p_2)
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
	#for i in range(0,frontDepI+1):
	for i in range(0,len(points2_offset)):
		foreAngDiffs.append(fabs(diffAngle(normalizeAngle(points2_offset[i][2] + pi), forePathAngle)))

	backAngDiffs = []
	#for i in range(backDepI,len(points2_offset)):
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
		departurePoint1 = pathPoints[max1]
		
		departurePoint1 = points2_offset[newFrontDepI]
		
		isExist1 = True

		if max1 == 0 or max1 == len(pathPoints)-1:
			isInterior1 = False
		else:
			isInterior1 = True

	if maxBack > DEP_THRESH:
		departurePoint2 = pathPoints[max2]

		departurePoint2 = points2_offset[newBackDepI]

		isExist2 = True
		
		if max2 == 0 or max2 == len(pathPoints)-1:
			isInterior2 = False
		else:
			isInterior2 = True

	maxContig, overlapSum
	contigFrac = float(maxContig)/float(len(points2_offset))

	angDiff2 = abs(diffAngle(initPose2[2], estPose2[2]))
	
	print "returning:", angDiff2, contigFrac, overlapSum
		
	" sum of closest points on front and back "
	" select the one with minimal cost "
	
	if plotIter:
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
		#pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
		pylab.title("nodeID %d: %1.2f %1.2f %d %d %d %d %1.2f" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2, maxContig, overlapSum))
		pylab.savefig("multi_distances_%04u.png" % pathPlotCount)

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


		xP = []
		yP = []
		for p in currPath:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='r')
		
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		#pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
		#pylab.title("%d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2))
		pylab.title("%d %s: [%d,%d] [%d,%d] %1.2f %1.2f" % (nodeID, repr(pathIDs), isExist1, isExist2, isInterior1, isInterior2, contigFrac, angDiff2))
		pylab.savefig("multi_departure_%04u.png" % pathPlotCount)
			
	print "multi_departure %d: %1.2f %1.2f %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], [%d,%d]" % (nodeID, angDiff2, contigFrac, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, frontDepI, backDepI)
	
	
	" if the medial axis does not overlap the path contiguously enough, mark a high discrepancy "
	#if float(maxContig)/float(len(points2_offset)) < 0.4:
	#	dist1 = 1e100
	#	dist2 = 1e100
	
	#return departurePoint, isInterior, isExist
	#return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2
	return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum, angDiff2


def globalOverlapICP_GPU2(initGuess, globalPath, medialPoints, globalPlotCount = 0, plotIter = False, n1 = 0, n2 = 0, minMatchDist = 1.0):
	
	numIterations = 0
	#globalPlotCount = 0
	
	def computeOffset(point1, point2, ang1, ang2):
	
		" corner points and orientations "
		corner1Pose = [point1[0], point1[1], ang1]
		corner2Pose = [point2[0], point2[1], ang2]
		
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
	
	#globalSpline = SplineFit(globalPath, smooth=0.1)
	#medialSpline = SplineFit(medialPoints, smooth=0.1)

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

	" transform the new pose "
	poseOrigin = Pose(currPose)
	
	costThresh = 0.004
	#minMatchDist = 2.0
	#minMatchDist = 1.0
	lastCost = 1e100
	
	startIteration = numIterations

	" set the initial guess "
	poseOrigin = Pose(currPose)
	
	" sample a series of points from the medial curves "
	" augment points with point-to-line covariances "
	" treat the points with the point-to-line constraint "
	globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
	#globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)
	#localPoints = addGPACVectorCovariance(localSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
	points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
	#points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)

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
		
		" transformed points without associated covariance "
		#posePoly = []
		#for p in points_offset:
		#	posePoly.append([p[0],p[1]])
		
		#match_pairs = matchPairs2(points, points_offset, globalPoints, minMatchDist)

		#print points_offset

		#cProfile.run('match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)', 'match_prof')
		#cProfile.runctx('match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)', globals(), locals(), 'match_prof')

		match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)	   
		
		"""
		" get the circles and radii "
		radius2, center2 = computeEnclosingCircle(points_offset)

		for i in range(len(globalPoints)):
			p_1 = globalPoly[i]
			#p_1 = localPoly[i]
	
			if isInCircle(p_1, radius2, center2):
		
				" for every point of B, find it's closest neighbor in transformed A "
				p_2, i_2, minDist = findClosestPointInA(points_offset, p_1)
	
				if minDist <= minMatchDist:
		
					" add to the list of match pairs less than 1.0 distance apart "
					" keep A points and covariances untransformed "
					C2 = points[i_2][2]
					
					C1 = globalPoints[i][2]
	
					" we store the untransformed point, but the transformed covariance of the A point "
					match_pairs.append([points[i_2],globalPoints[i],C1,C2])
		"""
		#input_GPU = convertToGPU(match_pairs)

		#oldCost = medialOverlapCostFunc([currU, currAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)
		#oldCost = medialOverlapCostFunc_GPU([currU, currAng], input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1)
		

		#newParam = scipy.optimize.fmin(medialOverlapCostFunc_GPU, [currU,currAng], [input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1], disp = 0)
		
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
		#newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), initGuess, c_poses_1, c_poses_2, len(poses_1))

		newU = newParam[0]
		newAng = newParam[1]
		
		#newCost = medialOverlapCostFunc([newU, newAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)
		#newCost = medialOverlapCostFunc_GPU([newU, newAng], input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1)


		#saveStr = ""
		#saveStr += "initGuess = " + repr(initGuess) + '\n'
		#saveStr += "match_pairs = " + repr(match_pairs) + '\n'
		#saveStr += "numPairs = " + repr(len(match_pairs)) + '\n'
		#saveStr += "poses_1 = " + repr(poses_1) + '\n'
		#saveStr += "poses_2 = " + repr(poses_2) + '\n'

		#f = open("costTest.txt", 'w')
		#f.write(saveStr)
		#f.close()

		#exit()
		" set the current parameters "
		currAng = normalizeAngle(newAng)
		currU = newU


		#if currU+0.02 > 1.0:
		#	pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
		#elif currU < 0.0:
		#	pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
		#else:
		#	pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]
			
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
		
		#if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
		#	break



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


		#plotEnv()		
		pylab.title("(%u,%u) u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, n2, u1, currU, currAng, newCost))

		pylab.xlim(currPose[0]-4, currPose[0]+4)					
		pylab.ylim(currPose[1]-3, currPose[1]+3)
		pylab.savefig("spline_plot_%04u.png" % globalPlotCount)
		pylab.clf()
		
		" save inputs "
		#saveFile = ""
		#saveFile += "initGuess = " + repr(initGuess) + "\n"
		#saveFile += "globalPath = " + repr(globalPath) + "\n"
		#saveFile += "medialPoints = " + repr(medialPoints) + "\n"

		#f = open("icpInputSave_%04u.txt" % globalPlotCount, 'w')
		#globalPlotCount += 1
		#f.write(saveFile)
		#f.close()		
			
		#numIterations += 1


	return currPose, newCost, len(match_pairs)
	

def computeVectorCovariance(vec,x_var,y_var):

	c11 = x_var
	c12 = c21 = 0.0
	c22 = y_var

	mag = sqrt(vec[0]**2 + vec[1]**2)
	normVec = [vec[0]/mag, vec[1]/mag]

	if normVec[1] == 0:
		r22 = r11 = 1.0
		r12 = r21 = 0.0

	else:
		B = -1 / (normVec[1] + normVec[0]**2/normVec[1])
		A = -normVec[0]*B/normVec[1]
		#R = [[A, -B], [B, A]]

		r22 = r11 = A
		r12 = -B
		r21 = B

	res11 = r11*(c11*r11 + c12*r21) + r21*(c21*r11 + c22*r21)
	res12 = r11*(c11*r12 + c12*r22) + r21*(c21*r12 + c22*r22)
	res21 = r12*(c11*r11 + c12*r21) + r22*(c21*r11 + c22*r21)
	res22 = r12*(c11*r12 + c12*r22) + r22*(c21*r12 + c22*r22)
	
	Ca = [[res11, res12], [res21, res22]]
	
	return Ca


def addGPACVectorCovariance(points, high_var=1.0, low_var=0.001):

	newPoints = []

	for p in points:
		#C = numpy.matrix([	[high_var, 0.0],
		#		[0.0, high_var]
		#		])
		C = [ [high_var, 0.], [0., high_var]]

		# the covariance matrix that enforces the point-to-plane constraint
		normVec = [cos(p[2]+pi/2.0), sin(p[2]+pi/2.0)]		
		C = computeVectorCovariance(normVec,low_var,high_var)


		newPoints.append([p[0], p[1], C])

	return newPoints


def findClosestPointInA(a_trans, b):

	minDist = 1e100
	minPoint = None
	min_i = 0

	for i in range(len(a_trans)):
		p = a_trans[i]

		dist = sqrt((p[0]-b[0])**2 + (p[1]-b[1])**2)
		if dist < minDist:
			minPoint = copy(p)
			minDist = dist
			min_i = i

	if minPoint != None:
		return minPoint, min_i, minDist
	else:
		raise

" displace the point by the offset only.  No covariance "
def dispOffset(p, offset):
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	px = p[0]
	py = p[1]

	p_off = [px*cos(theta) - py*sin(theta) + xd, px*sin(theta) + py*cos(theta) + yd]
	
	return p_off

" displace the point by the offset plus modify it's covariance "
def dispPoint(p, offset):
	
	xd = offset[0]
	yd = offset[1]
	theta = offset[2]

	#T = numpy.matrix([	[math.cos(theta), -math.sin(theta), xd],
	#		[math.sin(theta), math.cos(theta), yd],
	#		[0.0, 0.0, 1.0]
	#		])
	
	px = p[0]
	py = p[1]
	
	tx = px*cos(theta) - py*sin(theta) + xd
	ty = px*sin(theta) + py*cos(theta) + yd
	
	#p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
	#temp = T*p_hom
	#p_off = [temp[0,0],temp[1,0]]
	p_off = [tx, ty]

	Cv = p[2]
	
	r11 = cos(theta)
	r12 = -sin(theta)
	r21 = sin(theta)
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
