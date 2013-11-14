import multiprocessing as processing
import ctypes, os
from SplineFit import SplineFit
import gen_icp
from Pose import Pose
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath
from functions import *
from StableCurve import StableCurve
from Paths import computePathAngleVariance
import math
from operator import itemgetter

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

def __remote_multiMovePath(rank, qin, qout, nodeID, direction, distEst):

	print "started __remote_multiMovePath"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		print rank, "received", nc, args
		# process data
		#knn = __do_nothing(data, nc, someArg2, someArg3)
		results = []
		for arg in args:
			mapHyp = arg[0]
			movePath(mapHyp, nodeID, direction, distEst = distEst)
			results.append(mapHyp)
						   
		# write to output queue
		qout.put((nc,results))



#def batchMovePath(, splices, medial, initPose, pathIDs, nodeID):
def batchMovePath(mapHyps, nodeID, direction, distEst = 1.0):


	global renderGlobalPlotCount


	#initGuess

	#arg = [initGuess, initPose, pathID, nodeID]

	ndata = len(mapHyps)
	
	args = []
	
	for k in range(ndata):
		arg = [mapHyps[k], k + renderGlobalPlotCount]
		args.append(arg)
	
	renderGlobalPlotCount += len(mapHyps)
	
	
	nproc = __num_processors()
	
	nproc *= 2
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
	#__remote_multiMovePath(rank, qin, qout, nodeID, direction, distEst):
	pool = [processing.Process(target=__remote_multiMovePath,
				args=(rank, qin, qout, nodeID, direction, distEst))
					for rank in range(nproc)]
	for p in pool: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	print "args =", args
	while 1:
		#_data = data[:,cur:cur+chunk_size]
		_data = args[cur:cur+chunk_size]
		print "_data =", _data
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		print "put(", (nc,_data), ")"
		qin.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout.get()]
	
	print "received output:", knn
		
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "sorted output:", knn
		

	print "terminating pool"
	# terminate workers
	for p in pool:
		print "terminate"
		p.terminate()
		print "terminated"
		
	print "returning"
	return knn


	

def movePath(mapHyp, nodeID, direction, distEst = 1.0):
	
	print "movePath(", nodeID, ",", direction, ",", distEst, ")"

	poseData = mapHyp.poseData
	
	if nodeID > 0:
		
		if nodeID % 2 == 1:
			getStepGuess(mapHyp, nodeID-1, direction)
			getInPlaceGuess(mapHyp, nodeID-1, nodeID, direction)
		else:
			print "movePath:  DO NOTHING"
		
		" guess path state from these guesses "
			
		" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
		if poseData.numNodes >= 4 and nodeID % 2 == 1:
			
			" 1) get spliced path that the previous node is on "
			hull2 = poseData.aHulls[nodeID-1]
			medial2 = poseData.medialAxes[nodeID-1]
			hull3 = poseData.aHulls[nodeID]
			medial3 = poseData.medialAxes[nodeID]
			
			print "hypothesis", mapHyp.hypothesisID, len(mapHyp.paths[0])
			allSplices, terminals, junctions = mapHyp.getAllSplices(plotIter = False)
			print "hypothesis", mapHyp.hypothesisID, terminals, junctions


			initPose2 = mapHyp.nodePoses[nodeID-1]
			initPose3 = mapHyp.nodePoses[nodeID]
			
			print "junctions:", junctions
			print "initPose2:", initPose2
			print "initPose3:", initPose3
			
			junctions2 = []
			junctions3 = []
			for pathID, params in junctions.iteritems():
				junctionPose = params[1]
				
				print "checking", pathID, params
				
				dist2 = sqrt((junctionPose[0]-initPose2[0])**2 + (junctionPose[1]-initPose2[1])**2)
				dist3 = sqrt((junctionPose[0]-initPose3[0])**2 + (junctionPose[1]-initPose3[1])**2)
				
				print "dist2 =", dist2
				print "dist3 =", dist3
				
				if dist2 < 3.0:
					junctions2.append((pathID,params[2][0]))

				if dist3 < 3.0:
					junctions3.append((pathID,params[2][0]))
						
				
			closePathID2 = mapHyp.getClosestPath(initPose2)

			print "closePathID2:", closePathID2
			
			pathSet2 = [closePathID2]

			junctionSet = junctions2 + junctions3
			junctionSet = list(set(junctionSet))

			print "junctions2:", junctions2
			print "junctions3:", junctions3

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
				parentID = mapHyp.getParentPathID(pathID)

				if parentID != None and not parentID in pathSet2:
					termSet2.append(junctions[pathID][2])
				

			

			splicePaths = []

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
							
									
			resultMoves2 = []
			resultMoves3 = []

			for path in splicePaths:		

				" 2) get pose of previous node, get point on path curve "
				pose0 = mapHyp.nodePoses[nodeID-3]
				pose1 = mapHyp.nodePoses[nodeID-2]

				" FIXME:  make sure path is oriented correctly wrt node medial axis "
				hull0 = poseData.aHulls[nodeID-3]
				medial0 = poseData.medialAxes[nodeID-3]
				
				poseOrigin0 = Pose(pose0)
				globalMedial0 = []
				for p in medial0:
					globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))
	
				orientedSplicePath = orientPath(path, globalMedial0)				
				currPathSpline = SplineFit(orientedSplicePath, smooth=0.1)

				
				print "pose0,pose1:", pose0, pose1
				
				minDist0, u0, p0 = currPathSpline.findClosestPoint(pose0[0:2])
				minDist1, u1, p1 = currPathSpline.findClosestPoint(pose1[0:2])
				
				print "u0,u1:", u0, u1
				print "len(distPoints):", len(currPathSpline.distPoints)
				
				arcDist0 = currPathSpline.distPoints[int(u0*1000)]
				arcDist1 = currPathSpline.distPoints[int(u1*1000)]

				
				" 3) step distance along the path in direction of travel "
				
				if direction:
					arcDist2 = arcDist0 - distEst
					arcDist3 = arcDist1 - distEst
				else:
					arcDist2 = arcDist0 + distEst
					arcDist3 = arcDist1 + distEst
				
				print "arcDist0, arcDist1, arcDist2, arcDist3:", arcDist0, arcDist1, arcDist2, arcDist3
				
				pose2 = initPose2
				pose3 = initPose3
				
				
				print "pose2,pose3:", pose2, pose3
				
				" 4) set as pose of new node "
				mapHyp.nodePoses[nodeID-1] = pose2
				mapHyp.nodePoses[nodeID] = pose3

				uPath2, uMedialOrigin2 = selectLocalCommonOrigin(orientedSplicePath, medial2, pose2)

				uPath3, uMedialOrigin3 = selectLocalCommonOrigin(orientedSplicePath, medial3, pose3)
				
				u2 = uPath2
				u3 = uPath3

				print "input: uMedialOrigin2, u2, pose2:", uMedialOrigin2, u2, pose2
				print "input: uMedialOrigin3, u3, pose3:", uMedialOrigin3, u3, pose3

				
				resultPose2, lastCost2, matchCount2, currAng2 = gen_icp.globalPathToNodeOverlapICP2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial2, plotIter = False, n1 = nodeID-1, n2 = -1, arcLimit = 0.1)
				resultPose3, lastCost3, matchCount3, currAng3 = gen_icp.globalPathToNodeOverlapICP2([u3, uMedialOrigin3, 0.0], orientedSplicePath, medial3, plotIter = False, n1 = nodeID, n2 = -1, arcLimit = 0.1)
				
				print "resultPoses:", resultPose2, resultPose3

				multiDepCount = 0

				result2 = getMultiDeparturePoint(orientedSplicePath, medial2, pose2, resultPose2, [], nodeID-1, pathPlotCount = multiDepCount, plotIter = False)
				multiDepCount += 1
				result3 = getMultiDeparturePoint(orientedSplicePath, medial3, pose3, resultPose3, [], nodeID, pathPlotCount = multiDepCount, plotIter = False)
				multiDepCount += 1
				
				#results1.append(result+(k,))
				" (departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 )"
				
				" (resultPose2,lastCost2,matchCount2,fabs(currAng2)) "
				
				angDiff2 = abs(diffAngle(pose2[2],resultPose2[2]))
				angDiff3 = abs(diffAngle(pose3[2],resultPose3[2]))
				
				print "angDiff:", angDiff2, angDiff3
				
				" NOTE:  added minimum threshold to angle difference "
				" NOTE:  guess pose is now the inter-nodal estimate instead of path-based estimate "
				if angDiff2 < 0.5:
					#resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(currAng2)) + result2)
					resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(0.0)) + result2)
				if angDiff3 < 0.5:
					#resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(currAng3)) + result3)				
					resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(0.0)) + result3)				
			
			
			
			resultMoves2 = sorted(resultMoves2, key=itemgetter(18))
			resultMoves2 = sorted(resultMoves2, key=itemgetter(16), reverse=True)

			resultMoves3 = sorted(resultMoves3, key=itemgetter(18))
			resultMoves3 = sorted(resultMoves3, key=itemgetter(16), reverse=True)
			
			print "resultMoves2:"
			for res in resultMoves2:
				print res
			print "resultMoves3:"
			for res in resultMoves3:
				print res

			if len(resultMoves2) > 0:
				mapHyp.nodePoses[nodeID-1] = resultMoves2[0][0]
			else:
				print "node", nodeID-1, "not movePathed because no valid pose"
			if len(resultMoves3) > 0:
				mapHyp.nodePoses[nodeID] = resultMoves3[0][0]
			else:
				print "node", nodeID, "not movePathed because no valid pose"




def getInPlaceGuess(mapHyp, nodeID1, nodeID2, direction):
	
	poseData = mapHyp.poseData

	" PERFORM INPLACE CONSTRAINT BETWEEN PAIR "
	supportLine = mapHyp.paths[0]
	
	transform = makeInPlaceConstraint(mapHyp, nodeID1, nodeID2)
	transform1 = transform
	offset1 = [transform[0,0], transform[1,0], transform[2,0]]			
	
	if len(supportLine) == 0:
		resultSum1 = 1e100
	else:
		resultSum1 = checkSupport(mapHyp, nodeID1, nodeID2, offset1, supportLine)
	

	if poseData.numLeafs[nodeID1] > 2 or poseData.numLeafs[nodeID2] > 2:
		transform, overHist = makeMultiJunctionMedialOverlapConstraint(mapHyp, nodeID1, nodeID2, isMove = False, inPlace = False, isForward = direction )
	else:
		transform, overHist = makeMedialOverlapConstraint(mapHyp, nodeID1, nodeID2, isMove = False, inPlace = True, isForward = direction)
	
	
	transform3 = transform
	offset3 = [transform[0,0], transform[1,0], transform[2,0]]			
	if len(supportLine) == 0:
		resultSum3 = 1e100
	else:
		resultSum3 = checkSupport(mapHyp, nodeID1, nodeID2, offset3, supportLine)

	print "INPLACE sums:", resultSum1, resultSum3

	poseOrigin = Pose(mapHyp.nodePoses[nodeID1])

	if len(supportLine) > 0 and resultSum1 < resultSum3 and resultSum3 - resultSum1 > 100.0:
		estPose2 = poseOrigin.convertLocalOffsetToGlobal(offset1)
		mapHyp.nodePoses[nodeID2] = estPose2
	else:
		estPose2 = poseOrigin.convertLocalOffsetToGlobal(offset3)
		mapHyp.nodePoses[nodeID2] = estPose2


def getStepGuess(mapHyp, nodeID, direction):
	
	poseData = mapHyp.poseData

	" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
	if nodeID >= 2:
		
		if poseData.numLeafs[nodeID-2] > 2 or poseData.numLeafs[nodeID] > 2:
			transform, hist1 = makeMultiJunctionMedialOverlapConstraint(mapHyp, nodeID-2, nodeID, isMove = True, isForward = direction )
		else:			
			transform1, hist1 = makeMedialOverlapConstraint(mapHyp, nodeID-2, nodeID, isMove = True, isForward = direction )
			if hist1[2] > 0 or hist1[1] > 10 and hist1[2] == 0:
				transform2, hist2 = makeMedialOverlapConstraint(mapHyp, nodeID-2, nodeID, isMove = True, isForward = not direction )

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
		poseOrigin = Pose(mapHyp.nodePoses[nodeID-2])
		estPose2 = poseOrigin.convertLocalOffsetToGlobal(offset)
		mapHyp.nodePoses[nodeID] = estPose2


@logFunction
def makeInPlaceConstraint(mapHyp, nodeID1, nodeID2):

	poseData = mapHyp.poseData

	" compute the new posture "
	" 1. compute the GPAC pose "
	" 2. perform correction of GPAC pose "
	" 3. compute location of rootPose from corrected GPAC pose "
	
	" node1 is the front poke node "
	" nodes is the back poke node "

	originPosture = poseData.correctedPostures[nodeID1]
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

	newPosture = poseData.correctedPostures[nodeID2]

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
	cost1, matchCount1 = computeMedialError(mapHyp, nodeID1, nodeID2, offset1)
	cost2, matchCount2 = computeMedialError(mapHyp, nodeID1, nodeID2, offset2)
				
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
def computeMedialError(mapHyp, i, j, offset, minMatchDist = 2.0, tail1=0, tail2=0):

	poseData = mapHyp.poseData

	medial1 = poseData.medialLongPaths[i][tail1]
	medial2 = poseData.medialLongPaths[j][tail2]
	
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
def checkSupport(mapHyp, nodeID1, nodeID2, offset, supportLine):

	poseData = mapHyp.poseData

	hull2 = poseData.aHulls[nodeID2]
	medial2 = poseData.medialAxes[nodeID2]
	
	estPose1 = mapHyp.nodePoses[nodeID1]		
	
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
def makeMultiJunctionMedialOverlapConstraint(mapHyp, nodeID1, nodeID2, isMove = True, isForward = True, inPlace = False, uRange = 1.5):

	poseData = mapHyp.poseData

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

	estPose1 = mapHyp.nodePoses[nodeID1]		
	estPose2 = mapHyp.nodePoses[nodeID2]

	originProfile = Pose(estPose1)
	diffOffset = originProfile.convertGlobalPoseToLocal(estPose2)
	
	initDiff1 = diffAngle(estPose2[2], estPose1[2])

	results = []

	for k in range(len(poseData.medialLongPaths[nodeID1])):
		medial1 = mapHyp.medialLongPaths[nodeID1][k]
		for l in range(len(poseData.medialLongPaths[nodeID2])):
			medial2 = mapHyp.medialLongPaths[nodeID2][l]
			

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
			
				" create the ground constraints "
				gndGPAC1Pose = mapHyp.gndPoses[nodeID1]
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = mapHyp.gndPoses[nodeID2]
				gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)
				
				result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)

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
		
		
				medialError, matchCount = computeMedialError(nodeID1, nodeID2, offset, minMatchDist = 0.5, tail1=k, tail2=l)

		
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

				gndGPAC1Pose = mapHyp.gndPoses[nodeID1]
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = mapHyp.gndPoses[nodeID2]
				gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)

				
				result, hist = gen_icp.overlapICP(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, orientedMedial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = (nodeID1,k), n2 = (nodeID2,l), uRange = uRange)

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
		
		
				medialError, matchCount = computeMedialError(mapHyp, nodeID1, nodeID2, offset, minMatchDist = 0.5, tail1=k, tail2=l)

		
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
def makeMedialOverlapConstraint(mapHyp, nodeID1, nodeID2, isMove = True, isForward = True, inPlace = False, uRange = 0.1 ):

	poseData = mapHyp.poseData

	#print "recomputing hulls and medial axis"
	" compute the medial axis for each pose "
	
	posture1 = poseData.correctedPostures[nodeID1]
	posture2 = poseData.correctedPostures[nodeID2]

	hull1 = poseData.aHulls[nodeID1]
	medial1 = poseData.medialAxes[nodeID1]
	hull2 = poseData.aHulls[nodeID2]
	medial2 = poseData.medialAxes[nodeID2]

	estPose1 = mapHyp.nodePoses[nodeID1]
	estPose2 = mapHyp.nodePoses[nodeID2]
	
	medialSpline1 = SplineFit(medial1, smooth=0.1)
	medialSpline2 = SplineFit(medial2, smooth=0.1)

	originU1 = medialSpline1.findU([0.0,0.0])	
	originU2 = medialSpline2.findU([0.0,0.0])	

	distEst = 0.5

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

	" create the ground constraints "
	gndGPAC1Pose = mapHyp.gndPoses[nodeID1]
	currProfile = Pose(gndGPAC1Pose)

	gndGPAC2Pose = mapHyp.gndPoses[nodeID2]
	gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)



	result, hist = gen_icp.overlapICP_GPU2(estPose1, gndOffset, [u1, u2, angGuess], hull1, hull2, medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = nodeID1, n2 = nodeID2, uRange = uRange)

	transform = matrix([[result[0]], [result[1]], [result[2]]])
	
	print "making overlap constraint:", result[0], result[1], result[2]

	return transform, hist


@logFunction
def selectLocalCommonOrigin(globalPath, medial1, estPose1):
	poseOrigin = Pose(estPose1)
	
	" FIXME:  change function so that it always returns a common pair, biasing towards it's current location "
	" alternately, make a new function that gives the option to enforce locality to common pair "
	
	
	#globalMedial = []
	#for p in medial1:
	#	globalMedial.append(poseOrigin.convertLocalToGlobal(p))
	
	#medialSpline1 = SplineFit(globalMedial, smooth=0.1)

	globalSpline = SplineFit(globalPath, smooth=0.1)
	medialSpline1 = SplineFit(medial1, smooth=0.1)


	globalSamples = globalSpline.getUniformSamples(spacing = 0.04)
	medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
	
	globalMedialSamples = []
	for p in medialSamples:
		result = poseOrigin.convertLocalOffsetToGlobal(p)	
		globalMedialSamples.append(result)
	
	" compute the local variance of the angle "
	globalVar = computePathAngleVariance(globalSamples)
	medialVar = computePathAngleVariance(globalMedialSamples)

  
	" now lets find closest points and save their local variances "			
	closestPairs = []
	allPairs = []
	TERM_DIST = 20

	pathRail = int(len(globalSamples) / 20.0)
	medialRail = int(len(globalMedialSamples) / 20.0)

	print "rails:", len(globalSamples), len(globalMedialSamples), pathRail, medialRail

	
	for i in range(pathRail, len(globalSamples)-pathRail):
		pG = globalSamples[i]
		minDist = 1e100
		minJ = -1
		for j in range(medialRail, len(globalMedialSamples)-medialRail):
			pM = globalMedialSamples[j]
			dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
			
			
			
			if dist < minDist:
				minDist = dist
				minJ = j
		
		#poseDist = math.sqrt((pG[0]-estPose1[0])**2 + (pG[1]-estPose1[1])**2)
		#if poseDist < 0.3:
		if True:
			allPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))
				
		#if minDist < 0.1:
		#	closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

	for j in range(medialRail, len(globalMedialSamples)-medialRail):
		pM = globalMedialSamples[j]
		minDist = 1e100
		minI = -1
		for i in range(pathRail, len(globalSamples)-pathRail):
			pG = globalSamples[i]
			dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
			
			
			
			if dist < minDist:
				minDist = dist
				minI = i

		#pG = globalSamples[minI]
		#poseDist = math.sqrt((pG[0]-estPose1[0])**2 + (pG[1]-estPose1[1])**2)
		#if poseDist < 0.3:
		if True:		
			allPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))
		
		#if minDist < 0.1:
		#	closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))



	" remove duplicates "
	allPairs = list(set(allPairs))



	
	allPairs = sorted(allPairs, key=itemgetter(2))
	
	maxDistThresh = allPairs[-1][2]
	minDistThresh = 0.1
	
	print "minDistThresh,maxDistThresh =", allPairs[0][2], allPairs[-1][2]


	if len(allPairs) == 0:
		raise

	
	originU2 = 0.5
	originU1 = 0.5
	
	while minDistThresh <= maxDistThresh:
		closestPairs = []
		for pair in allPairs:			
			if pair[2] < minDistThresh:				
				closestPairs.append(pair)
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(5,6))

		print len(closestPairs), "closest pairs for dist", minDistThresh

		if len(closestPairs) > 0:
			originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
			originU1 = globalSpline.findU(globalSamples[closestPairs[0][0]])
			
			break
		
		minDistThresh += 0.1
	
	u2 = originU2
	u1 = originU1
	angGuess = 0.0
	
	return u1, u2
		


