import multiprocessing as processing
import ctypes, os, sys
from SplineFit import SplineFit
import gen_icp
from Pose import Pose
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath, orientPathLean
from functions import *
from StableCurve import StableCurve
import math
from operator import itemgetter
import time
import pylab
from landmarks import *
from shoots import *

import traceback 
import cProfile

renderGlobalPlotCount = 0

qin_move = None
qout_move =  None
pool_move = []

qin_eval = None
qout_eval =  None
pool_eval = []

qin_localize = None
qout_localize = None
pool_localize = []

qin_generate = None
qout_generate =  None
pool_generate = []

def __num_processors():
	
	return 4
	
	if os.name == 'nt': # Windows
		return int(os.getenv('NUMBER_OF_PROCESSORS'))
	else: # glibc (Linux, *BSD, Apple)
		get_nprocs = ctypes.cdll.libc.get_nprocs
		get_nprocs.restype = ctypes.c_int
		get_nprocs.argtypes = []
		return get_nprocs()


@logFunction
def getInPlaceGuess(poseData, nodeID1, nodeID2, estPose1, estPose2, supportLine, direction):
	
	#poseData = mapHyp.poseData

	" PERFORM INPLACE CONSTRAINT BETWEEN PAIR "
	#supportLine = mapHyp.paths[0]
	
	#posture1 = poseData.correctedPostures[nodeID1]
	#posture2 = poseData.correctedPostures[nodeID2]
	#medial1 = poseData.medialLongPaths[nodeID1][0]
	#medial2 = poseData.medialLongPaths[nodeID2][0]

	#estPose1 = mapHyp.nodePoses[nodeID1]		
	#estPose2 = mapHyp.nodePoses[nodeID2]


	transform = makeInPlaceConstraint(poseData, nodeID1, nodeID2)

	transform1 = transform
	offset1 = [transform[0,0], transform[1,0], transform[2,0]]			
	
	if len(supportLine) == 0:
		resultSum1 = 1e100
	else:
		#def checkSupport(estPose1, medial2, nodeID1, nodeID2, offset, supportLine):
		medial2 = poseData.medialAxes[nodeID2]
		#estPose1 = mapHyp.nodePoses[nodeID1]		
		resultSum1 = checkSupport(estPose1, medial2, nodeID1, nodeID2, offset1, supportLine)
	

	if poseData.numLeafs[nodeID1] > 2 or poseData.numLeafs[nodeID2] > 2:
		transform, overHist = makeMultiJunctionMedialOverlapConstraint(poseData, nodeID1, nodeID2, estPose1, estPose2, isMove = False, inPlace = False, isForward = direction )
	else:
		transform, overHist = makeMedialOverlapConstraint(poseData, nodeID1, nodeID2, estPose1, estPose2, isMove = False, inPlace = True, isForward = direction)
	
	
	transform3 = transform
	offset3 = [transform[0,0], transform[1,0], transform[2,0]]			
	if len(supportLine) == 0:
		resultSum3 = 1e100
	else:
		#def checkSupport(estPose1, medial2, nodeID1, nodeID2, offset, supportLine):
		medial2 = poseData.medialAxes[nodeID2]
		#estPose1 = mapHyp.nodePoses[nodeID1]		
		resultSum3 = checkSupport(estPose1, medial2, nodeID1, nodeID2, offset3, supportLine)

	print "INPLACE sums:", resultSum1, resultSum3

	#poseOrigin = Pose(mapHyp.nodePoses[nodeID1])
	poseOrigin = Pose(estPose1)

	if len(supportLine) > 0 and resultSum1 < resultSum3 and resultSum3 - resultSum1 > 100.0:
		newEstPose2 = poseOrigin.convertLocalOffsetToGlobal(offset1)
		#mapHyp.nodePoses[nodeID2] = estPose2
	else:
		newEstPose2 = poseOrigin.convertLocalOffsetToGlobal(offset3)
		#mapHyp.nodePoses[nodeID2] = estPose2

	return newEstPose2

@logFunction
def getStepGuess(poseData, nodeID0, nodeID2, estPose0, estPose2, direction, distEst = 0.5):
	
	#poseData = mapHyp.poseData

	" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
	if nodeID2 >= 2:

		#estPose0 = mapHyp.nodePoses[nodeID0]		
		#estPose2 = mapHyp.nodePoses[nodeID2]
		
		if poseData.numLeafs[nodeID0] > 2 or poseData.numLeafs[nodeID2] > 2:
			transform, hist1 = makeMultiJunctionMedialOverlapConstraint(poseData, nodeID0, nodeID2, estPose0, estPose2, isMove = True, isForward = direction )
		else:			
			transform1, hist1 = makeMedialOverlapConstraint(poseData, nodeID0, nodeID2, estPose0, estPose2, isMove = True, isForward = direction, distEst = distEst )
			if hist1[2] > 0 or hist1[1] > 10 and hist1[2] == 0:
				transform2, hist2 = makeMedialOverlapConstraint(poseData, nodeID0, nodeID2, estPose0, estPose2, isMove = True, isForward = not direction , distEst = distEst)

				if hist1[2] < hist2[2]:
					transform = transform1
				else:
					if hist1[1] <= hist2[1]:
						transform = transform1
					else:
						transform = transform2
			else:
				transform = transform1

		offset = [transform[0,0], transform[1,0], transform[2,0]]			
		poseOrigin = Pose(estPose0)
		newEstPose2 = poseOrigin.convertLocalOffsetToGlobal(offset)
		return newEstPose2

	#mapHyp.nodePoses[nodeID2] = estPose2
	
	return estPose2


@logFunction
def makeInPlaceConstraint(poseData, nodeID1, nodeID2):

	#poseData = mapHyp.poseData
	originPosture = poseData.correctedPostures[nodeID1]
	newPosture = poseData.correctedPostures[nodeID2]
	medial1 = poseData.medialLongPaths[nodeID1][0]
	medial2 = poseData.medialLongPaths[nodeID2][0]

	" compute the new posture "
	" 1. compute the GPAC pose "
	" 2. perform correction of GPAC pose "
	" 3. compute location of rootPose from corrected GPAC pose "
	
	" node1 is the front poke node "
	" nodes is the back poke node "


	#originPosture = posture1
	#newPosture = posture2

	originCurve = StableCurve(originPosture)
	originForePose, originBackPose = originCurve.getPoses()
	
	originForeProfile = Pose(originForePose)
	originBackProfile = Pose(originBackPose)
	
	originForePosture = []
	originBackPosture = []
	for i in range(len(originPosture)):
		originForePosture.append(originForeProfile.convertGlobalPoseToLocal(originPosture[i]))
	for i in range(len(originPosture)):
		originBackPosture.append(originBackProfile.convertGlobalPoseToLocal(originPosture[i]))


	localPosture = newPosture
	localCurve = StableCurve(localPosture)

	localForePose, localBackPose = localCurve.getPoses()
	 
	localForeProfile = Pose(localForePose)
	localBackProfile = Pose(localBackPose)
	
	localForePosture = []
	localBackPosture = []
	for i in range(len(localPosture)):
		localForePosture.append(localForeProfile.convertGlobalPoseToLocal(localPosture[i]))
	for i in range(len(localPosture)):
		localBackPosture.append(localBackProfile.convertGlobalPoseToLocal(localPosture[i]))

	foreCost = 0.0
	for j in range(0,20):
		p1 = originForePosture[j]
		p2 = localForePosture[j]
		
		foreCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)		   
	
	backCost = 0.0
	for j in range(20,39):
		p1 = originBackPosture[j]
		p2 = localBackPosture[j]
		
		backCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)		   

	angle = 0.0

	correctedGPACPose1 = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
	correctedGPACPose2 = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
	correctedProfile1 = Pose(correctedGPACPose1)
	correctedProfile2 = Pose(correctedGPACPose2)

	localRootOffset3_1 = correctedProfile1.convertGlobalPoseToLocal([0.0,0.0,0.0])
	localRootOffset3_2 = correctedProfile2.convertGlobalPoseToLocal([0.0,0.0,0.0])

	offset1 = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3_1)
	offset2 = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3_2)
	
	" except the case where one node has not updated "
	cost1, matchCount1 = computeMedialError(nodeID1, nodeID2, offset1, medial1, medial2)
	cost2, matchCount2 = computeMedialError(nodeID1, nodeID2, offset2, medial1, medial2)
				
	if cost1 < cost2:
		correctedGPACPose = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
	else:
		correctedGPACPose = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])

	correctedProfile = Pose(correctedGPACPose)

	" offset from 0,0 to newGPACPose "
	" rootPose from neg. offset from correctedGPACPose "
	localRootOffset3 = correctedProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
	
	if cost1 < cost2:
		offset = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3)
	else:
		offset = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3)
	
	transform = matrix([[offset[0]], [offset[1]], [offset[2]]])
	
	return transform


@logFunction
def computeMedialError(nodeID1, nodeID2, offset, medial1, medial2, minMatchDist = 2.0):

	#poseData = mapHyp.poseData
	#medial1 = poseData.medialLongPaths[nodeID1][tail1]
	#medial2 = poseData.medialLongPaths[nodeID2][tail2]
	
	medialSpline1 = SplineFit(medial1, smooth=0.1)
	medialSpline2 = SplineFit(medial2, smooth=0.1)

	points1 = gen_icp.addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
	points2 = gen_icp.addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)

	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	points2_offset = []
	for p in points2:
		result = gen_icp.dispPoint(p, offset)		
		points2_offset.append(result)


	poly1 = []		
	for p in points1:
		poly1.append([p[0],p[1]])	
	
	" find the matching pairs "
	match_pairs = []
	
	" transformed points without associated covariance "
	poly2 = []
	for p in points2_offset:
		poly2.append([p[0],p[1]])	
	
	" get the circles and radii "
	radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
	radius1, center1 = gen_icp.computeEnclosingCircle(points1)

	for i in range(len(points2_offset)):
		p_2 = poly2[i]

		if gen_icp.isInCircle(p_2, radius1, center1):

			" for every transformed point of A, find it's closest neighbor in B "
			p_1, minDist = gen_icp.findClosestPointInB(points1, p_2, [0.0,0.0,0.0])

			if minDist <= minMatchDist:
	
				" add to the list of match pairs less than 1.0 distance apart "
				" keep A points and covariances untransformed "
				C2 = points2_offset[i][2]
				C1 = p_1[2]

				" we store the untransformed point, but the transformed covariance of the A point "
				match_pairs.append([points2_offset[i],p_1,C2,C1])

	for i in range(len(points1)):
		p_1 = poly1[i]

		if gen_icp.isInCircle(p_1, radius2, center2):
	
			" for every point of B, find it's closest neighbor in transformed A "
			p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)

			if minDist <= minMatchDist:
	
				" add to the list of match pairs less than 1.0 distance apart "
				" keep A points and covariances untransformed "
				C2 = points2_offset[i_2][2]
				
				C1 = points1[i][2]

				" we store the untransformed point, but the transformed covariance of the A point "
				match_pairs.append([points2_offset[i_2],points1[i],C2,C1])

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
	
		val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
					
		vals.append(val)
		sum1 += val


	matchCount = len(vals)

	return sum1, matchCount



@logFunction
def checkSupport(estPose1, medial2, nodeID1, nodeID2, offset, supportLine):

	#poseData = mapHyp.poseData

	#medial2 = poseData.medialAxes[nodeID2]
	#estPose1 = mapHyp.nodePoses[nodeID1]		
	
	minMatchDist2 = 2.0
		
	" set the initial guess "
	poseOrigin = Pose(estPose1)
	
	localSupport = []
	for pnt in supportLine:
		localSupport.append(poseOrigin.convertGlobalToLocal(pnt))

	supportSpline = SplineFit(localSupport, smooth=0.1)		
	supportPoints = gen_icp.addGPACVectorCovariance(supportSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
		
	medialSpline2 = SplineFit(medial2, smooth=0.1)
	points2 = gen_icp.addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)


	" transform pose 2 by initial offset guess "	
	" transform the new pose "
	points2_offset = []
	for p in points2:
		result = gen_icp.dispPoint(p, offset)		
		points2_offset.append(result)

	" transformed points without associated covariance "
	poly2 = []
	for p in points2_offset:
		poly2.append([p[0],p[1]])			
	
	" get the circles and radii "
	radius2, center2 = gen_icp.computeEnclosingCircle(points2_offset)
			
	support_pairs = []
	for i in range(len(poly2)):
		
		p_2 = poly2[i]

		" for every transformed point of A, find it's closest neighbor in B "
		p_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

		if gen_icp.isInCircle(p_1, radius2, center2):

			if minDist <= minMatchDist2:
				C2 = points2[i][2]
				C1 = p_1[2]

				" we store the untransformed point, but the transformed covariance of the A point "
				support_pairs.append([points2[i],p_1,C2,C1])

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
	
		val = gen_icp.computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
		
		vals.append(val)
		sum1 += val
		
	return sum1



@logFunction
def makeMultiJunctionMedialOverlapConstraint(poseData, nodeID1, nodeID2, estPose1, estPose2, isMove = True, isForward = True, inPlace = False, uRange = 1.5):

	#poseData = mapHyp.poseData

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


	hull1 = poseData.aHulls[nodeID1]
	hull2 = poseData.aHulls[nodeID2]

	#estPose1 = mapHyp.nodePoses[nodeID1]		
	#estPose2 = mapHyp.nodePoses[nodeID2]

	originProfile = Pose(estPose1)
	diffOffset = originProfile.convertGlobalPoseToLocal(estPose2)
	
	initDiff1 = diffAngle(estPose2[2], estPose1[2])

	results = []

	for k in range(len(poseData.medialLongPaths[nodeID1])):
		medial1 = poseData.medialLongPaths[nodeID1][k]
		for l in range(len(poseData.medialLongPaths[nodeID2])):
			medial2 = poseData.medialLongPaths[nodeID2][l]
			

			if isMove:
				
				medialSpline1 = SplineFit(medial1, smooth=0.1)
				medialSpline2 = SplineFit(medial2, smooth=0.1)

				pose1 = medialSpline1.getUVecSet([0.5, 0.5+0.02])[0]	
				pose2 = medialSpline2.getUVecSet([0.5, 0.5+0.02])[0]
			
				point1 = [pose1[0],pose1[1]]
				point2 = [pose2[0],pose2[1]]
				ang1 = pose1[2]
				ang2 = pose2[2]
			
				offset = computeOffset(point1, point2, ang1, ang2)

				points1 = medialSpline1.getUniformSamples()
				points2 = medialSpline2.getUniformSamples()
			
			
				" transform the past poses "
				points2_offset = []
				medialProfile = Pose(offset)
				for p in points2:
					result = medialProfile.convertLocalOffsetToGlobal(p)
					points2_offset.append(result)

				" transform the past poses "
				medial2_offset = []
				for p in medial2:
					result = medialProfile.convertLocalToGlobal(p)
					medial2_offset.append(result)

				" if isMove is True, then we align medial2 with medial1 "
				points2R_offset = deepcopy(points2_offset)
				points2R_offset.reverse()
						
				overlapMatch = []
				angleSum1 = 0.0
				angleSum2 = 0.0
				for i in range(0,len(points1)):
					p_1 = points1[i]
					p_2, j, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
					if minDist < 0.5:
						p_2R, j, minDist = gen_icp.findClosestPointInA(points2R_offset, p_2)
						
						ang1 = abs(diffAngle(p_1[2], p_2[2]))
						ang2 = abs(diffAngle(p_1[2], p_2R[2]))
						
						angleSum1 += ang1
						angleSum2 += ang2
				
				" select global path orientation based on which has the smallest angle between tangent vectors "
				print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
				if angleSum1 > angleSum2:
					medial2Reverse = deepcopy(medial2)
					medial2Reverse.reverse()
					orientedMedial2 = deepcopy(medial2Reverse)
				else:
					orientedMedial2 = deepcopy(medial2)
				
				#orientedMedial2 = orientPath(medial2, medial1)
				medialSpline1 = SplineFit(medial1, smooth=0.1)
				medialSpline2 = SplineFit(orientedMedial2, smooth=0.1)

		
				originU1 = medialSpline1.findU([0.0,0.0])	
				originU2 = medialSpline2.findU([0.0,0.0])	
				
				if isForward:
					
					moveDist = 0.3
						

					if len(poseData.frontProbeError[nodeID2]) > 0:
		
						frontSum = 0.0
						frontProbeError = poseData.frontProbeError[nodeID2]
						for n in frontProbeError:
							frontSum += n
						foreAvg = frontSum / len(frontProbeError)
										
						if foreAvg >= 1.4:
							u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
						else:	
							u2 = medialSpline2.getUOfDist(originU2, moveDist, distIter = 0.001)
					else:	
						u2 = medialSpline2.getUOfDist(originU2, moveDist, distIter = 0.001)
						
				else:
					moveDist = -0.3
					u2 = medialSpline2.getUOfDist(originU2, moveDist, distIter = 0.001)
				
				print "computed u2 =", u2, "from originU2 =", originU2

				u1 = originU1
				angGuess = 0.0
			
				result, hist = gen_icp.overlapICP(estPose1, [u1, u2, angGuess], medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)

				transform = matrix([[result[0]], [result[1]], [result[2]]])
				
				print "making overlap constraint:", result[0], result[1], result[2]
				
				angDiff = abs(diffAngle(diffOffset[2], transform[2,0]))			
				#totalGuesses.append((angDiff, result[2], result[3], result[4]))

				points1 = medialSpline1.getUniformSamples()
				points2 = medialSpline2.getUniformSamples()
				p_1 = medialSpline1.getU(0.5)
				
				
				offset = [transform[0,0], transform[1,0], transform[2,0]]
				fitProfile = Pose(offset)
				
				points2_offset = []
				for p in points2:
					result = fitProfile.convertLocalOffsetToGlobal(p)
					points2_offset.append(result)		
				
				overlapMatch = []
				matchCount = 0
				overlapSum = 0.0
				angleSum = 0.0
				for m in range(0,len(points2_offset)):
					p_2 = points2_offset[m]
					p_1, n, minDist = gen_icp.findClosestPointInA(points1, p_2)
	
					if minDist < 0.1:
						overlapMatch.append((m,n,minDist))
						matchCount += 1
						
						overlapSum += minDist
						
						ang1 = p_1[2]
						ang2 = p_2[2]

						angleSum += abs(diffAngle(ang1,ang2))
							
				if matchCount > 0:
					angleSum /= matchCount
					overlapSum /= matchCount
		
		
				#medialError, matchCount = computeMedialError(nodeID1, nodeID2, offset, minMatchDist = 0.5, tail1=k, tail2=l)

				#poseData = mapHyp.poseData
				medial1 = poseData.medialLongPaths[nodeID1][k]
				medial2 = poseData.medialLongPaths[nodeID2][l]
				medialError, matchCount = computeMedialError(nodeID1, nodeID2, offset, medial1, medial2, minMatchDist = 0.5)

		
				" None result used to be covariance matrix "
				results.append((hist[0], angleSum, medialError, matchCount, angDiff, k, l, transform, hist, matchCount, overlapSum, angleSum))


			else:
				
				medialSpline1 = SplineFit(medial1, smooth=0.1)
				medialSpline2 = SplineFit(medial2, smooth=0.1)

				pose1 = medialSpline1.getUVecSet([0.5, 0.5+0.02])[0]	
				pose2 = medialSpline2.getUVecSet([0.5, 0.5+0.02])[0]
			
				point1 = [pose1[0],pose1[1]]
				point2 = [pose2[0],pose2[1]]
				ang1 = pose1[2]
				ang2 = pose2[2]
			
				offset = computeOffset(point1, point2, ang1, ang2)

				points1 = medialSpline1.getUniformSamples()
				points2 = medialSpline2.getUniformSamples()
			
			
				" transform the past poses "
				points2_offset = []
				medialProfile = Pose(offset)
				for p in points2:
					result = medialProfile.convertLocalOffsetToGlobal(p)
					points2_offset.append(result)

				" transform the past poses "
				medial2_offset = []
				for p in medial2:
					result = medialProfile.convertLocalToGlobal(p)
					medial2_offset.append(result)

			
				" if isMove is True, then we align medial2 with medial1 "
				points2R_offset = deepcopy(points2_offset)
				points2R_offset.reverse()
						
				overlapMatch = []
				angleSum1 = 0.0
				angleSum2 = 0.0
				for i in range(0,len(points1)):
					p_1 = points1[i]
					p_2, j, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
					if minDist < 0.5:
						p_2R, j, minDist = gen_icp.findClosestPointInA(points2R_offset, p_2)
						
						ang1 = abs(diffAngle(p_1[2], p_2[2]))
						ang2 = abs(diffAngle(p_1[2], p_2R[2]))
						
						angleSum1 += ang1
						angleSum2 += ang2
				
				" select global path orientation based on which has the smallest angle between tangent vectors "
				print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
				if angleSum1 > angleSum2:
					medial2Reverse = deepcopy(medial2)
					medial2Reverse.reverse()
					orientedMedial2 = deepcopy(medial2Reverse)
				else:
					orientedMedial2 = deepcopy(medial2)


				#medialSpline1 = SplineFit(medial1, smooth=0.1)
				medialSpline2 = SplineFit(orientedMedial2, smooth=0.1)
		
				originU1 = medialSpline1.findU([0.0,0.0])	
				originU2 = medialSpline2.findU([0.0,0.0])	
		
				if inPlace:
					originU1 = 0.6
					originU2 = 0.4
					u2 = originU2
					print "computed u2 =", u2, "from originU2 =", originU2

				else:
					poseOrigin = Pose(estPose1)
					offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
					
					points2 = medialSpline2.getUniformSamples()
					p_1 = medialSpline1.getU(0.5)
					
					points2_offset = []
					for p in points2:
						result = gen_icp.dispOffset(p, offset)		
						points2_offset.append(result)
			
					p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)
			
					u2 = medialSpline2.findU(points2[i_2])	
					
				u1 = originU1
				angGuess = 0.0
			
				" create the ground constraints "

				result, hist = gen_icp.overlapICP(estPose1, [u1, u2, angGuess], medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)

				transform = matrix([[result[0]], [result[1]], [result[2]]])
				
				print "making overlap constraint:", result[0], result[1], result[2]
				
				angDiff = abs(diffAngle(diffOffset[2], transform[2,0]))			

				points1 = medialSpline1.getUniformSamples()
				points2 = medialSpline2.getUniformSamples()
				p_1 = medialSpline1.getU(0.5)
				
				
				offset = [transform[0,0], transform[1,0], transform[2,0]]
				fitProfile = Pose(offset)
				
				points2_offset = []
				for p in points2:
					result = fitProfile.convertLocalOffsetToGlobal(p)
					points2_offset.append(result)		
				
				overlapMatch = []
				matchCount = 0
				overlapSum = 0.0
				angleSum = 0.0
				for m in range(0,len(points2_offset)):
					p_2 = points2_offset[m]
					p_1, n, minDist = gen_icp.findClosestPointInA(points1, p_2)
	
					if minDist < 0.1:
						overlapMatch.append((m,n,minDist))
						matchCount += 1
						
						overlapSum += minDist
						
						ang1 = p_1[2]
						ang2 = p_2[2]

						angleSum += abs(diffAngle(ang1,ang2))
							
				if matchCount > 0:
					angleSum /= matchCount
					overlapSum /= matchCount
		
		
				#poseData = mapHyp.poseData
				medial1 = poseData.medialLongPaths[nodeID1][k]
				medial2 = poseData.medialLongPaths[nodeID2][l]
				medialError, matchCount = computeMedialError(nodeID1, nodeID2, offset, medial1, medial2, minMatchDist = 0.5)

		
				" None result used to be covariance matrix "
				results.append((hist[0], angleSum, medialError, matchCount, angDiff, k, l, transform, hist, matchCount, overlapSum, angleSum))
	
	results.sort(reverse=True)
	
	print "Multi-Junction Overlap of", nodeID1, "and", nodeID2
	selectedIndex = 0
	for k in range(len(results)):
		print results[k]
		

	for k in range(len(results)):
		if results[k][8][1] < 10 and results[k][8][2] == 0:
			selectedIndex = k
			break

	transform = results[selectedIndex][7]
	hist = results[selectedIndex][8]
	
	return transform, hist


@logFunction
def makeMedialOverlapConstraint(poseData, nodeID1, nodeID2, estPose1, estPose2, isMove = True, isForward = True, inPlace = False, uRange = 0.1, distEst = 0.5 ):
#def makeMedialOverlapConstraint(mapHyp, nodeID1, nodeID2, isMove = True, isForward = True, inPlace = False, uRange = 0.1 ):

	#poseData = mapHyp.poseData

	#print "recomputing hulls and medial axis"
	" compute the medial axis for each pose "
	
	posture1 = poseData.correctedPostures[nodeID1]
	posture2 = poseData.correctedPostures[nodeID2]

	medial1 = poseData.medialAxes[nodeID1]
	medial2 = poseData.medialAxes[nodeID2]

	#estPose1 = mapHyp.nodePoses[nodeID1]
	#estPose2 = mapHyp.nodePoses[nodeID2]
	
	medialSpline1 = SplineFit(medial1, smooth=0.1)
	medialSpline2 = SplineFit(medial2, smooth=0.1)

	originU1 = medialSpline1.findU([0.0,0.0])	
	originU2 = medialSpline2.findU([0.0,0.0])	

	#distEst = 0.5

	if inPlace:
		" FULL LENGTH MEDIAL AXIS "
		originU1 = 0.5
		originU2 = 0.5
		
		u2 = originU2
		print "computed u2 =", u2, "from originU2 =", originU2

	elif isMove:
		
		if isForward:
			
			if len(poseData.frontProbeError[nodeID2]) > 0:

				frontSum = 0.0
				frontProbeError = poseData.frontProbeError[nodeID2]
				for n in frontProbeError:
					frontSum += n
				foreAvg = frontSum / len(frontProbeError)
								
				if foreAvg >= 1.4:
					u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
				else:	
					u2 = medialSpline2.getUOfDist(originU2, distEst, distIter = 0.001)
			else:	
				u2 = medialSpline2.getUOfDist(originU2, distEst, distIter = 0.001)
				
		else:

			#u2 = medialSpline2.getUOfDist(originU2, -0.3, distIter = 0.001)

			if len(poseData.backProbeError[nodeID2]) > 0:

				backSum = 0.0
				backProbeError = poseData.backProbeError[nodeID2]
				for n in backProbeError:
					backSum += n
				backAvg = backSum / len(backProbeError)
								
				if backAvg >= 1.4:
					u2 = medialSpline2.getUOfDist(originU2, 0.0, distIter = 0.001)
				else:	
					u2 = medialSpline2.getUOfDist(originU2, -distEst, distIter = 0.001)
			else:	
				u2 = medialSpline2.getUOfDist(originU2, -distEst, distIter = 0.001)
			
			
			#u2 = 0.4
		print "computed u2 =", u2, "from originU2 =", originU2
		 
	else:
		poseOrigin = Pose(estPose1)
		offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
		
		points2 = medialSpline2.getUniformSamples()
		p_1 = medialSpline1.getU(0.5)
		
		points2_offset = []
		for p in points2:
			result = gen_icp.dispOffset(p, offset)		
			points2_offset.append(result)

		p_2, i_2, minDist = gen_icp.findClosestPointInA(points2_offset, p_1)

		u2 = medialSpline2.findU(points2[i_2])	
		
		if u2 > 0.9 or u2 < 0.1:
			raise

	u1 = originU1
	angGuess = 0.0

	#result, hist = gen_icp.overlapICP_GPU2(estPose1, [u1, u2, angGuess], medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = nodeID1, n2 = nodeID2, uRange = uRange)
	result, hist = gen_icp.overlapICP_GPU2(estPose1, [u1, u2, angGuess], medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = True, n1 = nodeID1, n2 = nodeID2, uRange = uRange)

	transform = matrix([[result[0]], [result[1]], [result[2]]])
	
	print "making overlap constraint:", result[0], result[1], result[2]

	return transform, hist

@logFunction
def generateAll(hypSet):
	for pID, mapHyp in hypSet.iteritems():
		generateMap(mapHyp)

@logFunction
def generateMap(mapHyp):

	mapHyp.generatePaths()
	" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "


@logFunction
def addToPaths2(shootIDs, particleIDs, hypSet, nodeID1, nodeID2):

	poseData = hypSet.values()[0].poseData

	newHyps = {}

	" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
	isFirst = False
	for pID, mapHyp in hypSet.iteritems():
		if len(mapHyp.paths[0]) == 0:
			print "this is first nodes, just add them already!"
			isFirst = True

	if isFirst:
		for pID, mapHyp in hypSet.iteritems():
			mapHyp.addNode(nodeID1,0)
			mapHyp.addNode(nodeID2,0)
	
			mapHyp.generatePaths()

			mapHyp.initializePoseParticles()

		return shootIDs, particleIDs, hypSet


	for pID, mapHyp in hypSet.iteritems():
		generateMap(mapHyp)

	""" check branching of local spline from member path """
	computeLocalDivergence2(hypSet, nodeID1, nodeID2)

	for pID, mapHyp in hypSet.iteritems():
		mapHyp.isNotBranched = True
		mapHyp.isSplit = False


	""" check for branching conditions from path, spawn new hypotheses """
	hypSet, shootIDs, particleIDs = checkForeBranch2(hypSet, nodeID1, nodeID2, shootIDs, particleIDs)

	hypSet, shootIDs, particleIDs = checkBackBranch2(hypSet, nodeID1, nodeID2, shootIDs, particleIDs)

	addNodesToShoot(hypSet, nodeID1, nodeID2)


	""" if map is not branching, don't regenerate the maps with the new poses yet
		Keep them exactly how they were with two less poses """

	genSet = {}
	nonSet = {}
	for pID, currHyp in hypSet.iteritems():
		if not currHyp.isNotBranched:
			genSet[pID] = currHyp
		else:
			nonSet[pID] = currHyp


	for pID, mapHyp in genSet.iteritems():
		generateMap(mapHyp)

	nonSet.update(genSet)
	hypSet = nonSet

	return shootIDs, particleIDs, hypSet

@logFunction
def addNodesToShoot(hypSet, nodeID1, nodeID2):

	poseData = hypSet.values()[0].poseData

	for pID, mapHyp in hypSet.iteritems():


		medialAxis1 = mapHyp.poseData.medialAxes[nodeID1]
		medialSpline1 = SplineFit(medialAxis1, smooth=0.1)
		medial1_vec = medialSpline1.getUniformSamples()
		estPose1 = mapHyp.getNodePose(nodeID1)
		poseFrame1 = Pose(estPose1)
		globalMedial1 = []
		for p in medial1_vec:
			globalMedial1.append(poseFrame1.convertLocalOffsetToGlobal(p))

		medialAxis2 = mapHyp.poseData.medialAxes[nodeID2]
		medialSpline2 = SplineFit(medialAxis2, smooth=0.1)
		medial2_vec = medialSpline2.getUniformSamples()
		estPose2 = mapHyp.getNodePose(nodeID2)
		poseFrame2 = Pose(estPose2)
		globalMedial2 = []

		for p in medial2_vec:
			globalMedial2.append(poseFrame2.convertLocalOffsetToGlobal(p))
		
		
		" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "

		memberShootIDs1 = mapHyp.memberShootIDs1
		print "memberShootIDs1:", memberShootIDs1

		isNewID = False
		for shootID in memberShootIDs1:
			if shootID not in mapHyp.globalLongPaths.keys():
				isNewID = True
				mapHyp.addNode(nodeID1,shootID)

		if not isNewID:


			resultSet1 = []

			for shootID in memberShootIDs1:
				shootLongSplices = mapHyp.globalLongPaths[shootID]

				for longSplice in shootLongSplices:

					orientedSplice = orientPathLean(longSplice, globalMedial1)

					resultArgs1 = getMultiDeparturePoint(orientedSplice, medial1_vec, estPose1, estPose1, [shootID], nodeID1, pathPlotCount=0, hypID=mapHyp.hypothesisID, plotIter=False)

					resultSet1.append(resultArgs1+(shootID,))

			""" overlapSum secondary, contigFrac primary """
			resultSet1 = sorted(resultSet1, key=itemgetter(13), reverse=False)
			resultSet1 = sorted(resultSet1, key=itemgetter(12), reverse=True)


			isAdded1 = False
			for result in resultSet1:
				contigFrac0 = result[12]
				shootID = result[15]
				isInterior1 = result[2]
				isInterior2 = result[8]

				if not isInterior1 and not isInterior2:
					""" add nodeID1 to this shoot """
					mapHyp.addNode(nodeID1,shootID)
					isAdded1 = True
					break

			if not isAdded1:
				shootID = resultSet1[0][15]
				mapHyp.addNode(nodeID1,shootID)



		memberShootIDs2 = mapHyp.memberShootIDs2
		print "memberShootIDs2:", memberShootIDs2

		isNewID = False
		for shootID in memberShootIDs2:
			if shootID not in mapHyp.globalLongPaths.keys():
				isNewID = True
				mapHyp.addNode(nodeID2,shootID)

		if not isNewID:

			resultSet2 = []

			for shootID in memberShootIDs2:
				shootLongSplices = mapHyp.globalLongPaths[shootID]

				for longSplice in shootLongSplices:

					orientedSplice = orientPathLean(longSplice, globalMedial2)

					resultArgs2 = getMultiDeparturePoint(orientedSplice, medial2_vec, estPose2, estPose2, [shootID], nodeID2, pathPlotCount=0, hypID=mapHyp.hypothesisID, plotIter=False)

					#departurePoint1, depAngle1, isInterior1, isExist1, discDist1, maxFront, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, maxBack, contigFrac, overlapSum, angDiff2 = resultArgs2


					resultSet2.append(resultArgs2+(shootID,))


			""" overlapSum secondary, contigFrac primary """
			resultSet2 = sorted(resultSet2, key=itemgetter(13), reverse=False)
			resultSet2 = sorted(resultSet2, key=itemgetter(12), reverse=True)
			
			isAdded2 = False
			for result in resultSet2:
				contigFrac0 = result[12]
				shootID = result[15]
				isInterior1 = result[2]
				isInterior2 = result[8]

				if not isInterior1 and not isInterior2:
					""" add nodeID2 to this shoot """
					mapHyp.addNode(nodeID2,shootID)
					isAdded2 = True
					break

			if not isAdded2:
				shootID = resultSet2[0][15]
				mapHyp.addNode(nodeID2,shootID)
		
		#if isAdded1 or isAdded2:
		#	print isAdded1, isAdded2, "node not found non-diverging shoot splice"
		#	#raise


@logFunction
def computeLocalDivergence2(hypSet, nodeID1, nodeID2):

	poseData = hypSet.values()[0].poseData

	for pID, mapHyp in hypSet.iteritems():


		medialAxis1 = mapHyp.poseData.medialAxes[nodeID1]
		medialSpline1 = SplineFit(medialAxis1, smooth=0.1)
		medial1_vec = medialSpline1.getUniformSamples()
		estPose1 = mapHyp.getNodePose(nodeID1)
		poseFrame1 = Pose(estPose1)
		globalMedial1 = []
		for p in medial1_vec:
			globalMedial1.append(poseFrame1.convertLocalOffsetToGlobal(p))

		medialAxis2 = mapHyp.poseData.medialAxes[nodeID2]
		medialSpline2 = SplineFit(medialAxis2, smooth=0.1)
		medial2_vec = medialSpline2.getUniformSamples()
		estPose2 = mapHyp.getNodePose(nodeID2)
		poseFrame2 = Pose(estPose2)
		globalMedial2 = []

		for p in medial2_vec:
			globalMedial2.append(poseFrame2.convertLocalOffsetToGlobal(p))
		
		
		" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "

		allSplices = mapHyp.getAllSplices2()
		resultSet1 = []
		resultSet2 = []

		#for sPath in allSplices:

		for k in range(len(allSplices)):

			sPath = allSplices[k]
			path = sPath['skelPath']
			termPoints = sPath['termPath']
			memberShootIDs = sPath['memberShootIDs']

			orientedSplice = orientPathLean(path, globalMedial1)

			resultArgs1 = getMultiDeparturePoint(orientedSplice, medial1_vec, estPose1, estPose1, memberShootIDs, nodeID1, pathPlotCount=k, hypID=mapHyp.hypothesisID, plotIter=False)




			resultArgs2 = getMultiDeparturePoint(orientedSplice, medial2_vec, estPose2, estPose2, memberShootIDs, nodeID2, pathPlotCount=k, hypID=mapHyp.hypothesisID, plotIter=False)

			#mapHyp.departureResultSet2 = resultArgs2

			departurePoint1, depAngle1, isInterior1, isExist1, discDist1, maxFront, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, maxBack, contigFrac, overlapSum, angDiff2 = resultArgs2


			resultSet1.append(resultArgs1+(k,))
			resultSet2.append(resultArgs2+(k,))

			"departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2"



		""" overlapSum secondary, contigFrac primary """
		resultSet1 = sorted(resultSet1, key=itemgetter(13), reverse=False)
		resultSet1 = sorted(resultSet1, key=itemgetter(12), reverse=True)
		resultSet2 = sorted(resultSet2, key=itemgetter(13), reverse=False)
		resultSet2 = sorted(resultSet2, key=itemgetter(12), reverse=True)
		
		#print "getOrderedOverlappingPaths() sorted node", nodeID1, ":"
		#for result in resultSet1:
		#	print result[15], result[12], result[13]

		# FIXME:  temporary forced result 
		if len(resultSet1) > 1:
			contigFrac0 = resultSet1[0][12]
			contigFrac1 = resultSet1[1][12]
			index0 = resultSet1[0][15]
			index1 = resultSet1[1][15]

			sPath0 = allSplices[index0]
			sPath1 = allSplices[index1]
			memberShootIDs0 = sPath0['memberShootIDs']
			memberShootIDs1 = sPath1['memberShootIDs']

			diff = fabs(contigFrac0-contigFrac1)

			if diff < 0.05:

				if len(memberShootIDs0) < len(memberShootIDs1): 
					result = resultSet1[0]	 
				else:
					result = resultSet1[1]	 
			else:
				result = resultSet1[0]	 

		else:
			result = resultSet1[0]	 


		result1 = result
		spliceIndex1 = result1[15]
		sPath1 = allSplices[spliceIndex1]

		mapHyp.departureResultSet1 = result1
		mapHyp.overlapSplice1 = sPath1
		mapHyp.memberShootIDs1 = sPath1['memberShootIDs']

		if len(resultSet2) > 1:
			contigFrac0 = resultSet2[0][12]
			contigFrac1 = resultSet2[1][12]
			index0 = resultSet2[0][15]
			index1 = resultSet2[1][15]

			sPath0 = allSplices[index0]
			sPath1 = allSplices[index1]
			memberShootIDs0 = sPath0['memberShootIDs']
			memberShootIDs1 = sPath1['memberShootIDs']

			diff = fabs(contigFrac0-contigFrac1)

			if diff < 0.05:

				if len(memberShootIDs0) < len(memberShootIDs1): 
					result = resultSet2[0]	 
				else:
					result = resultSet2[1]	 
			else:
				result = resultSet2[0]	 

		else:
			result = resultSet2[0]	 


		result2 = result
		spliceIndex2 = result2[15]
		sPath2 = allSplices[spliceIndex2]

		mapHyp.departureResultSet2 = result2
		mapHyp.overlapSplice2 = sPath2
		mapHyp.memberShootIDs2 = sPath2['memberShootIDs']


@logFunction
def checkForeBranch2(hypSet, nodeID1, nodeID2, shootIDs, particleIDs):
	
	" 60 degree threshold "
	ANG_THRESH = 1.047

	newHyps = {}
	poseData = hypSet.values()[0].poseData


	for pID, mapHyp in hypSet.iteritems():
		
		departurePoint1, depAngle1, isInterior1, isExist1, discDist1, maxFront, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, maxBack, contigFrac, overlapSum, angDiff2, spliceIndex = mapHyp.departureResultSet1

		departures1 = [isExist1,isExist2]
		interiors1 = [isInterior1, isInterior2]
		depPoints1 = [departurePoint1, departurePoint2]
		distances1 = [discDist1, discDist2]
		depAngles1 = [depAngle1, depAngle2]
		contig1 = (contigFrac, overlapSum)

		departurePoint1, depAngle1, isInterior1, isExist1, discDist1, maxFront, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, maxBack, contigFrac, overlapSum, angDiff2, spliceIndex = mapHyp.departureResultSet2

		departures2 = [isExist1,isExist2]
		interiors2 = [isInterior1, isInterior2]
		depPoints2 = [departurePoint1, departurePoint2]
		distances2 = [discDist1, discDist2]
		depAngles2 = [depAngle1, depAngle2]
		contig2 = (contigFrac, overlapSum)

		print "F node departures", nodeID1, ":", departures1
		print "F node departures", nodeID2, ":", departures2
		print "F node  interiors", nodeID1, ":", interiors1
		print "F node  interiors", nodeID2, ":", interiors2
		print "F node contiguity", nodeID1, ":", contig1
		print "F node contiguity", nodeID2, ":", contig2
		print "F node depPoints", nodeID1, ":", depPoints1
		print "F node depPoints", nodeID2, ":", depPoints2
		print "F node distances", nodeID1, ":", distances1
		print "F node distances", nodeID2, ":", distances2
		print "F node depAngles", nodeID1, ":", depAngles1
		print "F node depAngles", nodeID2, ":", depAngles2
		print "F node contig", nodeID1, ":", contig1
		print "F node contig", nodeID2, ":", contig2

		" new junction finding logic "
		" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
		" if a terminal departure exists that is internal, than we have a new junction "

		" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
		frontExist1 = departures1[0]
		frontInterior1 = interiors1[0]
		foreTerm1 = frontExist1 and frontInterior1

		" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
		
		" NODE2: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
		frontExist2 = departures2[0]
		frontInterior2 = interiors2[0]
		foreTerm2 = frontExist2 and frontInterior2

		depPoint1 = depPoints1[0]
		depPoint2 = depPoints2[0]

		depAngle1 = depAngles1[0]
		depAngle2 = depAngles2[0]

		frontAngDiff = diffAngle(depAngle1, depAngle2)
		
		

		isFront1 = poseData.faceDirs[nodeID1]
		isFront2 = poseData.faceDirs[nodeID2]
		if isFront1 and not isFront2:
			dirFlag = 0
		elif not isFront1 and isFront2:
			dirFlag = 1
		else:
			print isFront1, isFront2
			raise	


		sF1 = poseData.spatialFeatures[nodeID1][0]
		sF2 = poseData.spatialFeatures[nodeID2][0]

		hasLandmark1 = True
		hasLandmark2 = True

		if sF1["bloomPoint"] == None and sF1["archPoint"] == None and sF1["inflectionPoint"] == None:
			hasLandmark1 = False

		if sF2["bloomPoint"] == None and sF2["archPoint"] == None and sF2["inflectionPoint"] == None:
			hasLandmark2 = False


		if hasLandmark1 or hasLandmark2:
			hasSpatialFeature = True
			print "current pair", nodeID1, nodeID2, "has spatial landmark feature(s)", depPoint1, depPoint2
		else:
			hasSpatialFeature = False
			print "current pair", nodeID1, nodeID2, "has no spatial landmark feature and cannot branch", depPoint1, depPoint2

			#print "REJECT, diverging pose has no spatial landmark feature", nodeID1, pathID, neighborCount, depPoint
			#return False, parentPathID

		""" call with depAngle instead of depPoint[2] angle
			since depPoint[2] angle refers to orientation of medial axis curve """
		""" depAngle refers to tip angle, which is not necessarily what we want """
		""" depPoint[2] is now inverted for forward angle """
		parentPathID1 = 0
		parentPathID2 = 0
		if depPoint1 != 0:
			#isUnique1, duplicatePathID1 = mapHyp.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
			isUnique1, duplicatePathID1 = mapHyp.checkUniqueBranch2(parentPathID1, nodeID1, depPoint1[2], depPoint1)
		else:
			isUnique1 = False
			duplicatePathID1 = -1
			
		if depPoint2 != 0:
			#isUnique2, duplicatePathID2 = mapHyp.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)
			isUnique2, duplicatePathID2 = mapHyp.checkUniqueBranch2(parentPathID2, nodeID2, depPoint2[2], depPoint2)
		else:
			isUnique2 = False
			duplicatePathID2 = -1


		""" temporarily disable this branch requirement """
		hasSpatialFeature = True

		if isUnique1 and duplicatePathID1 != -1:
			if not hasSpatialFeature:
				isUnique1 = False
				duplicatePathID1 = parentPathID1

		if isUnique2 and duplicatePathID2 != -1:
			if not hasSpatialFeature:
				isUnique2 = False
				duplicatePathID2 = parentPathID2


		isBranched = False
		if foreTerm1 and foreTerm2:
			if fabs(frontAngDiff) < ANG_THRESH:
				if isUnique1 and isUnique2:
					isBranched = True
				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						isBranched = True
				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						isBranched = True
			else:
				if isUnique1 and isUnique2:
					if dirFlag == 0:
						isBranched = True
					if dirFlag == 1:
						isBranched = True
				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						isBranched = True
				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						isBranched = True
		elif foreTerm1 and not foreTerm2:
			if isUnique1:
				if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
					isBranched = True
				else:
					if dirFlag == 0:
						isBranched = True
		elif foreTerm2 and not foreTerm1:
			if isUnique2:
				if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:
					isBranched = True
				else:
					if dirFlag == 1:		
						isBranched = True


		#if isUnique1 or isUnique2:
		branchHyp = mapHyp
		if isBranched:

			" create a new map state where a branch decision is not made "
			mapHyp.isNodeBranching[nodeID1] = True
			mapHyp.isNodeBranching[nodeID2] = True

			#newMapHyp = mapHyp.copy(particleIDs)
			#newHyps[particleIDs] = newMapHyp
			#newHyps[particleIDs].isNotBranched = False 

			" in case current hypothesis has already branched, this ensures that topology is regenerated when branching "
			hypSet[pID].isNotBranched = False

			#print "creating hyp", particleIDs, "from hyp", mapHyp.hypothesisID, ", len(paths) =", len(mapHyp.pathClasses)
			#particleIDs += 1

			#branchHyp = newMapHyp 

		shootIDs, pathBranchIDs = branchHyp.determineBranchPair2(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, dirFlag, isUnique1, isUnique2, shootIDs)

		if pathBranchIDs[0] != None:
			branchHyp.memberShootIDs1.append(pathBranchIDs[0])

		if pathBranchIDs[1] != None:
			branchHyp.memberShootIDs2.append(pathBranchIDs[1])


	hypSet.update(newHyps)

	return hypSet, shootIDs, particleIDs

@logFunction
def checkBackBranch2(hypSet, nodeID1, nodeID2, shootIDs, particleIDs):

	" 60 degree threshold "
	ANG_THRESH = 1.047

	newHyps = {}
	poseData = hypSet.values()[0].poseData

	for pID, mapHyp in hypSet.iteritems():

		departurePoint1, depAngle1, isInterior1, isExist1, discDist1, maxFront, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, maxBack, contigFrac, overlapSum, angDiff2, spliceIndex = mapHyp.departureResultSet1

		departures1 = [isExist1,isExist2]
		interiors1 = [isInterior1, isInterior2]
		depPoints1 = [departurePoint1, departurePoint2]
		distances1 = [discDist1, discDist2]
		depAngles1 = [depAngle1, depAngle2]
		contig1 = (contigFrac, overlapSum)

		departurePoint1, depAngle1, isInterior1, isExist1, discDist1, maxFront, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, maxBack, contigFrac, overlapSum, angDiff2, spliceIndex = mapHyp.departureResultSet2

		departures2 = [isExist1,isExist2]
		interiors2 = [isInterior1, isInterior2]
		depPoints2 = [departurePoint1, departurePoint2]
		distances2 = [discDist1, discDist2]
		depAngles2 = [depAngle1, depAngle2]
		contig2 = (contigFrac, overlapSum)


		print "B node departures", nodeID1, ":", departures1
		print "B node departures", nodeID2, ":", departures2
		print "B node  interiors", nodeID1, ":", interiors1
		print "B node  interiors", nodeID2, ":", interiors2
		print "B node contiguity", nodeID1, ":", contig1
		print "B node contiguity", nodeID2, ":", contig2
		print "B node depPoints", nodeID1, ":", depPoints1
		print "B node depPoints", nodeID2, ":", depPoints2
		print "B node distances", nodeID1, ":", distances1
		print "B node distances", nodeID2, ":", distances2
		print "B node depAngles", nodeID1, ":", depAngles1
		print "B node depAngles", nodeID2, ":", depAngles2
		print "B node contig", nodeID1, ":", contig1
		print "B node contig", nodeID2, ":", contig2


		backExist1 = departures1[1]
		backInterior1 = interiors1[1]
		backTerm1 = backExist1 and backInterior1

		backExist2 = departures2[1]
		backInterior2 = interiors2[1]
		backTerm2 = backExist2 and backInterior2

		depPoint1 = depPoints1[1]
		depPoint2 = depPoints2[1]
		
		depAngle1 = depAngles1[1]
		depAngle2 = depAngles2[1]
		
		backAngDiff = diffAngle(depAngle1, depAngle2)



		isFront1 = poseData.faceDirs[nodeID1]
		isFront2 = poseData.faceDirs[nodeID2]
		if isFront1 and not isFront2:
			dirFlag = 1
		elif not isFront1 and isFront2:
			dirFlag = 0
		else:
			print isFront1, isFront2
			raise	

		sF1 = poseData.spatialFeatures[nodeID1][0]
		sF2 = poseData.spatialFeatures[nodeID2][0]

		hasLandmark1 = True
		hasLandmark2 = True

		if sF1["bloomPoint"] == None and sF1["archPoint"] == None and sF1["inflectionPoint"] == None:
			hasLandmark1 = False

		if sF2["bloomPoint"] == None and sF2["archPoint"] == None and sF2["inflectionPoint"] == None:
			hasLandmark2 = False

		if hasLandmark1 or hasLandmark2:
			hasSpatialFeature = True
			print "current pair", nodeID1, nodeID2, "has spatial landmark feature(s)", depPoint1, depPoint2
		else:
			hasSpatialFeature = False
			print "current pair", nodeID1, nodeID2, "has no spatial landmark feature and cannot branch", depPoint1, depPoint2


		parentPathID1 = 0
		parentPathID2 = 0
		if depPoint1 != 0:
			#isUnique1, duplicatePathID1 = mapHyp.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
			isUnique1, duplicatePathID1 = mapHyp.checkUniqueBranch2(parentPathID1, nodeID1, depPoint1[2], depPoint1)
		else:
			isUnique1 = False
			duplicatePathID1 = -1
			
		if depPoint2 != 0:
			#isUnique2, duplicatePathID2 = mapHyp.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)
			isUnique2, duplicatePathID2 = mapHyp.checkUniqueBranch2(parentPathID2, nodeID2, depPoint2[2], depPoint2)
		else:
			isUnique2 = False
			duplicatePathID2 = -1

		""" temporarily disable this branch requirement """
		hasSpatialFeature = True

		if isUnique1 and duplicatePathID1 != -1:
			if not hasSpatialFeature:
				isUnique1 = False
				duplicatePathID1 = parentPathID1

		if isUnique2 and duplicatePathID2 != -1:
			if not hasSpatialFeature:
				isUnique2 = False
				duplicatePathID2 = parentPathID2


		isBranched = False
		if backTerm1 and backTerm2:
			if fabs(backAngDiff) < ANG_THRESH:
				if isUnique1 and isUnique2:
					isBranched = True
				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						isBranched = True
				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						isBranched = True
			else:
				if isUnique1 and isUnique2:
					if dirFlag == 0:
						isBranched = True
					if dirFlag == 1:
						isBranched = True
				if isUnique1 and not isUnique2:
					if dirFlag == 0:
						isBranched = True
				if not isUnique1 and isUnique2:
					if dirFlag == 1:
						isBranched = True
		elif backTerm1 and not backTerm2:
			if isUnique1:
				if backExist2 and fabs(backAngDiff) < ANG_THRESH:
					isBranched = True
				else:
					if dirFlag == 0:
						isBranched = True
		elif backTerm2 and not backTerm1:
			if isUnique2:
				if backExist1 and fabs(backAngDiff) < ANG_THRESH:
					isBranched = True
				else:
					if dirFlag == 1:		
						isBranched = True

		#if isUnique1 or isUnique2:
		branchHyp = mapHyp
		if isBranched:

			" create a new map state where a branch decision is not made "
			mapHyp.isNodeBranching[nodeID1] = True
			mapHyp.isNodeBranching[nodeID2] = True

			#newMapHyp = mapHyp.copy(particleIDs)
			#newHyps[particleIDs] = newMapHyp
			#newHyps[particleIDs].isNotBranched = False 

			" in case current hypothesis has already branched, this ensures that topology is regenerated when branching "
			hypSet[pID].isNotBranched = False
			#print "creating hyp", particleIDs, "from hyp", mapHyp.hypothesisID, ", len(paths) =", len(mapHyp.pathClasses)
			
			#particleIDs += 1

			#branchHyp = newMapHyp


		shootIDs, pathBranchIDs = branchHyp.determineBranchPair2(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, dirFlag, isUnique1, isUnique2, shootIDs)

		if pathBranchIDs[0] != None:
			branchHyp.memberShootIDs1.append(pathBranchIDs[0])

		if pathBranchIDs[1] != None:
			branchHyp.memberShootIDs2.append(pathBranchIDs[1])


	hypSet.update(newHyps)
	
	return hypSet, shootIDs, particleIDs


