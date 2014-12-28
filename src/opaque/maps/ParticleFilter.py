import multiprocessing as processing
import random
from copy import deepcopy
from SplineFit import SplineFit
import pylab
from math import sqrt, pi
from Pose import Pose
import gen_icp
from math import acos, asin, fabs
from numpy import array
import ctypes, os, sys
from operator import itemgetter

from functions import diffAngle
from Splices import getMultiDeparturePoint, orientPath, getCurveOverlap, orientPathLean, getTipAngles
from MapProcess import getStepGuess, getInPlaceGuess
from shoots import selectLocalCommonOrigin

from scipy.spatial import KDTree, cKDTree
from scipy.cluster.hierarchy import fclusterdata, fcluster

from matplotlib.pyplot import show, savefig, clf
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from numpy.random import rand

import traceback 
import cProfile

renderGlobalPlotCount = 0

qin_posePart = None
qout_posePart =  None
pool_posePart = []

qin_dispPosePart = None
qout_dispPosePart =  None
pool_dispPosePart = []

qin_dispPosePart2 = None
qout_dispPosePart2 =  None
pool_dispPosePart2 = []

qin_localLandmark = None
qout_localLandmark = None
pool_localLandmark = []

qin_branch = None
qout_branch =  None
pool_branch = []

def printStack():

	flist = traceback.format_stack()
	flist = flist[:-1]
	
	printStr = ""
	for line in flist:
		printStr += line
		
	print printStr


def __num_processors():

	return processing.cpu_count()
	
	return 4
	
	if os.name == 'nt': # Windows
		return int(os.getenv('NUMBER_OF_PROCESSORS'))
	else: # glibc (Linux, *BSD, Apple)
		get_nprocs = ctypes.cdll.libc.get_nprocs
		get_nprocs.restype = ctypes.c_int
		get_nprocs.argtypes = []
		return get_nprocs()

def batchDisplaceParticles2(displaceJobs):

	global renderGlobalPlotCount
	global pool_dispPosePart2
	global qin_dispPosePart2
	global qout_dispPosePart2

	ndata = len(displaceJobs)
	
	args = []
	
	for k in range(ndata):
		arg = displaceJobs[k]
		args.append(arg)
	
	nproc = __num_processors()
	
	#nproc *= 2
	print "nproc =", nproc
	
	# compute chunk size
	chunk_size = ndata / nproc
	chunk_size = 2 if chunk_size < 2 else chunk_size
	
	print "chunk_size =", chunk_size
	print "max_size =", ndata/chunk_size
	
	# set up a pool of processes
	if len(pool_dispPosePart2) == 0:

		qin_dispPosePart2 = processing.Queue(maxsize=ndata/chunk_size)
		qout_dispPosePart2 = processing.Queue(maxsize=ndata/chunk_size)
		pool_dispPosePart2 = [processing.Process(target=__remote_displaceParticle2,
		#pool_dispPosePart = [processing.Process(target=__remote_prof_displaceParticle2,
				#args=(rank, qin_dispPosePart2, qout_dispPosePart2, splices, medial, initPose, pathIDs, nodeID))
				args=(rank, qin_dispPosePart2, qout_dispPosePart2))
				#args=(rank, qin_dispPosePart2, qout_dispPosePart2))
					for rank in range(nproc)]
		for p in pool_dispPosePart2: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	#print "args =", args
	while 1:
		_data = args[cur:cur+chunk_size]
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		qin_dispPosePart2.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	isFail = False
	while len(knn) < nc:
		thisKnn = qout_dispPosePart2.get()
		if thisKnn[0] != None:
			knn += [thisKnn]
		else:
			isFail = True
			break
	
	if isFail:

		print "isFail =", isFail
		print "knn =", knn


		qin_dispPosePart2.close()
		qin_dispPosePart2.join_thread()

		qout_dispPosePart2.close()
		qout_dispPosePart2.cancel_join_thread()

		for p in pool_dispPosePart2:
			p.terminate()

		raise


	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	for p in pool_dispPosePart2:
		p.terminate()
	pool_dispPosePart2 = []

	print "returning"
	return knn

def __remote_prof_displaceParticle2(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_displaceParticle(rank, qin, qout)', "particle_%d.prof" % pid )
		cProfile.runctx("__remote_displaceParticle2(rank, qin, qout)", globals(), locals(), "displaceParticle2_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


def __remote_displaceParticle2(rank, qin, qout):


	try:

		#sys.stdout = open("displaceParticle_" + str(os.getpid()) + ".out", "w")
		#sys.stderr = open("displaceParticle_" + str(os.getpid()) + ".err", "w")
		sys.stdout = open("displaceParticle2_" + str(rank) + ".out", "a")
		sys.stderr = open("displaceParticle2_" + str(rank) + ".err", "a")
		print 'module name:', __name__
		print 'parent process:', os.getppid()
		print 'process id:', os.getpid()

		print "started __remote_displaceParticle2"

		while 1:
			# read input queue (block until data arrives)
			results = []
			nc, args = qin.get()
			
			print "__remote_displaceParticle(", nc, len(args), args
			sys.stdout.flush()

			#foo = None
			#badVal = foo[0] 

			for job in args:
				
				#[particleIndex, nodeID3, prevPose0, prevPose1, initPose2, initPose3, supportLine, pathSplices2, pathSplices3]
				poseData = job[0]
				particleIndex = job[1]
				nodeID3 = job[2]
				prevPose0 = job[3]
				prevPose1 = job[4]
				initPose2 = job[5]
				initPose3 = job[6]
				supportLine = job[7]
				pathSplices2 = job[8]
				pathSplices3 = job[9]
				landmarks_G = job[10]
				landmarks_N = job[11]

				result = displaceParticle2( poseData, pathSplices2, pathSplices3, supportLine, nodeID3, initPose2, initPose3, prevPose0, prevPose1, particleIndex, landmarks_G, landmarks_N)
				results.append((particleIndex,) + result)

				print "result:", (particleIndex,) + result
				sys.stdout.flush()


			print "qout.put(", nc, results
			sys.stdout.flush()
							   
			# write to output queue
			qout.put((nc,results))

	except:
		print "Worker process failed. Exiting"
		#printStack()
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]
		
		sys.stdout.flush()
		sys.stderr.flush()
		qout.put((None,None))
		raise
	
	print "process exited incorrectly"



def batchLocalizeLandmark(localizeJobs):

	global renderGlobalPlotCount
	global pool_localLandmark
	global qin_localLandmark
	global qout_localLandmark

	ndata = len(localizeJobs)
	
	args = []
	
	for k in range(ndata):
		arg = localizeJobs[k]
		args.append(arg)
	
	nproc = __num_processors()
	
	#nproc *= 2
	print "nproc =", nproc
	
	# compute chunk size
	chunk_size = ndata / nproc
	chunk_size = 2 if chunk_size < 2 else chunk_size
	
	print "chunk_size =", chunk_size
	print "max_size =", ndata/chunk_size
	
	# set up a pool of processes
	if len(pool_localLandmark) == 0:

		qin_localLandmark = processing.Queue(maxsize=ndata/chunk_size)
		qout_localLandmark = processing.Queue(maxsize=ndata/chunk_size)
		pool_localLandmark = [processing.Process(target=__remote_localizeLandmark,
		#pool_dispPosePart = [processing.Process(target=__remote_prof_localizeLandmark,
				#args=(rank, qin_dispPosePart2, qout_dispPosePart2, splices, medial, initPose, pathIDs, nodeID))
				args=(rank, qin_localLandmark, qout_localLandmark))
				#args=(rank, qin_dispPosePart2, qout_dispPosePart2))
					for rank in range(nproc)]
		for p in pool_localLandmark: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	#print "args =", args
	while 1:
		_data = args[cur:cur+chunk_size]
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		qin_localLandmark.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	isFail = False
	while len(knn) < nc:
		thisKnn = qout_localLandmark.get()
		if thisKnn[0] != None:
			knn += [thisKnn]
		else:
			isFail = True
			break
	
	if isFail:

		print "isFail =", isFail
		print "knn =", knn


		qin_localLandmark.close()
		qin_localLandmark.join_thread()

		qout_localLandmark.close()
		qout_localLandmark.cancel_join_thread()

		for p in pool_localLandmark:
			p.terminate()

		raise


	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	for p in pool_localLandmark:
		p.terminate()
	pool_localLandmark = []

	print "returning"
	return knn

def __remote_prof_localizeLandmark(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_localizeLandmark(rank, qin, qout)', "particle_%d.prof" % pid )
		cProfile.runctx("__remote_localizeLandmark(rank, qin, qout)", globals(), locals(), "localizeLandmark_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


def __remote_localizeLandmark(rank, qin, qout):


	try:

		sys.stdout = open("localizeLandmark_" + str(rank) + ".out", "a")
		sys.stderr = open("localizeLandmark_" + str(rank) + ".err", "a")
		print 'module name:', __name__
		print 'parent process:', os.getppid()
		print 'process id:', os.getpid()

		print "started __remote_localizeLandmark"

		while 1:
			# read input queue (block until data arrives)
			results = []
			nc, args = qin.get()
			
			print "__remote_localizeLandmark(", nc, len(args), args
			sys.stdout.flush()

			#foo = None
			#badVal = foo[0] 

			for job in args:
				
				poseData = job[0]
				spliceIndex = job[1]
				pathID = job[2]
				nodeID3 = job[3]
				sampDist = job[4]
				estPose0 = job[5]
				pathSplices2 = job[6]
				landmarks_G = job[7]
				landmarks_N = job[8]

				result = localizeLandmarkPose(poseData, pathSplices2, pathID, nodeID3, spliceIndex, sampDist, estPose0, landmarks_G, landmarks_N)
				results.append((spliceIndex,nodeID3, pathID) + result)

				print "result:", (spliceIndex,nodeID3,pathID) + result
				sys.stdout.flush()


			print "qout.put(", nc, results
			sys.stdout.flush()
							   
			# write to output queue
			qout.put((nc,results))

	except:
		print "Worker process failed. Exiting"
		#printStack()
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]
		
		sys.stdout.flush()
		sys.stderr.flush()
		qout.put((None,None))
		raise
	
	print "process exited incorrectly"


def localizeLandmarkPose( poseData, pathSplices2, pathID, nodeID, sIndex, sampDist, estPose, landmarks_G, landmarks_N):

	sys.stdout.flush()

	medial0 = poseData.medialAxes[nodeID]

	medialSpline0 = SplineFit(medial0, smooth=0.1)
	medial0_vec = medialSpline0.getUniformSamples()

	direction = poseData.travelDirs[nodeID]

	currPose2 = estPose
	currProb2 = 0.0

	foreAngle0, backAngle0 = getTipAngles(medial0_vec, estPose)

	splicePaths = pathSplices2

	""" if the node has a spatial feature, localize it to the closest landmark """

	""" if the node has a landmark feature, then we try to localize on the junction point """
	""" FIXME: only use the first long path, we don't try and figure out the longest yet """
	currPoses = {nodeID: currPose2}

	""" if the node has been localized to a landmark, we need the updated version """
	currPose2 = currPoses[nodeID]

	currProb2 = 0.0

	" now we fit the guessed pose to a splice of the paths "
	" this constrains it to the previously explored environment "
	" all possible splices are used, and they are filtered to determine the best result "

	resultMoves2 = []

	for spliceIndex in range(len(splicePaths)):		

		path = splicePaths[spliceIndex]

		" 2) get pose of previous node, get point on path curve "
		pose0 = estPose

		poseOrigin0 = Pose(pose0)
		globalMedial0 = []
		for p in medial0:
			globalMedial0.append(poseOrigin0.convertLocalToGlobal(p))

		orientedSplicePath = orientPath(path, globalMedial0)				
		currPathSpline = SplineFit(orientedSplicePath, smooth=0.1)

		minDist0, u0, p0 = currPathSpline.findClosestPoint(pose0[0:2])
		
		arcDist0 = currPathSpline.distPoints[int(u0*1000)]
		
		" 3) step distance along the path in direction of travel "
		
		pose2 = currPose2
		
		uPath2, uMedialOrigin2 = selectLocalCommonOrigin(orientedSplicePath, medial0, pose2)
		
		u2 = uPath2

		print "input: uMedialOrigin2, u2, pose2:", uMedialOrigin2, u2, pose2
		
		resultPose2, lastCost2, matchCount2, currAng2, currU2 = gen_icp.globalPathToNodeOverlapICP2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial0, plotIter = False, n1 = nodeID-1, n2 = -1, arcLimit = 0.01)
		
		print "resultPoses:", resultPose2,

		multiDepCount = 0

		result2 = getMultiDeparturePoint(orientedSplicePath, medial0_vec, pose2, resultPose2, [], nodeID-1, hypID = 0, pathPlotCount = 0, particleIndex=spliceIndex, spliceIndex=spliceIndex, plotIter = False)
		multiDepCount += 1
		
		" (departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 )"
		
		" (resultPose2,lastCost2,matchCount2,fabs(currAng2)) "

		contigFrac_2 = result2[12]
		overlapSum2 = result2[13]

		foreAngle2, backAngle2 = getTipAngles(medial0_vec, resultPose2)

		angDiff2 = min(abs(diffAngle(foreAngle2,foreAngle0)), abs(diffAngle(backAngle2,backAngle0)))



		frame0 = Pose(resultPose2)
		landmark0_N = landmarks_N[nodeID]
		landmark0_G = None
		thresh0 = None
		poseSum = 0.0

		if landmark0_N != None:
			landmark0_G = (frame0.convertLocalToGlobal(landmark0_N[0]), landmark0_N[1], landmark0_N[2])
			thresh0 = landmark0_G[1]

		print nodeID, sIndex, "landmarks:", landmark0_G, landmarks_G

		LANDMARK_THRESH = 7.0
		for i in range(len(landmarks_G)):
			p1 = landmarks_G[i][0]
			threshG = landmarks_G[i][1]

			if landmark0_G != None:

				maxThresh = thresh0
				if threshG > thresh0:
					maxThresh = threshG

				dist0 = sqrt((p1[0]-landmark0_G[0][0])**2 + (p1[1]-landmark0_G[0][1])**2)
				if dist0 < LANDMARK_THRESH:
					poseSum += sqrt(dist0*dist0/(threshG*threshG + thresh0*thresh0))


		" NOTE:  added minimum threshold to angle difference "
		" NOTE:  guess pose is now the inter-nodal estimate instead of path-based estimate "

		if overlapSum2 > 1e10:
			newProb2 = 0.0	
		else:
			#newProb2 = (pi-angDiff2) * contigFrac_2
			if poseSum > 0.0:
				newProb2 = 1.0/poseSum
			else:
				newProb2 = 0.0

		print "angDiff/contigFrac:", angDiff2, contigFrac_2
		print "overlapSum:", overlapSum2
		print "newProb:", newProb2

		if angDiff2 < 0.5 and contigFrac_2 > 0.7:
			resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(0.0)) + result2 + (orientedSplicePath, newProb2))
	
	
	""" sort by angDiff follwed by contigFrac """
	#resultMoves2 = sorted(resultMoves2, key=itemgetter(18))
	#resultMoves2 = sorted(resultMoves2, key=itemgetter(16), reverse=True)

	resultMoves2 = sorted(resultMoves2, key=itemgetter(20))

	currSplice2 = []

	" select the best pose for each node "
	" record the splice used for the fit "
	if len(resultMoves2) > 0:
		currPose2 = resultMoves2[0][0]
		currSplice2 = resultMoves2[0][19]
		currProb2 = resultMoves2[0][20]

	else:
		print "node", nodeID-1, "not movePathed because no valid pose"
		currProb2 = 0.0

	" move the pose particles along their paths "	
	print "fitted poses:", currPose2
	print "currProbs:", currProb2

	#if True:
	if False:

		pylab.clf()


		for splicePath in pathSplices2:

			xP = []
			yP = []
			for p in splicePath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r', zorder=10, alpha=0.3)

		for result in resultMoves2:
			resultPose = result[0]

			poseOrigin = Pose(resultPose)
			
			points_offset = []
			for p in medial0_vec:
				p1 = poseOrigin.convertLocalOffsetToGlobal(p)
				points_offset.append(p1)

			xP = []
			yP = []
			for p in points_offset:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color='k', zorder=10, alpha=1.0)




		xP = []
		yP = []
		for point_G, pointThresh, pointName in landmarks_G:
			xP.append(point_G[0])
			yP.append(point_G[1])

		pylab.scatter(xP, yP, color='k', zorder=9, alpha=0.3)

		pylab.scatter([landmark0_G[0][0],], [landmark0_G[0][1],], color='k', zorder=9, alpha=1.0)


		pylab.axis("equal")
		pylab.title("pathID %d nodeID %d currProb %1.2f" % ( pathID, nodeID, currProb2 ))

		nameStr = "landmarkLocalize_%04u_%04u_%04u_%1.1f.png" % (nodeID, pathID, sIndex, sampDist)

		pylab.savefig(nameStr)
		print "saving", nameStr


	return currPose2, currProb2


def displaceParticle2( poseData, pathSplices2, pathSplices3, supportLine, nodeID3, initPose2, initPose3, prevPose0, prevPose1, particleIndex, landmarks_G, landmarks_N):

	#print "movePath(", nodeID, ",", direction, ",", distEst, ")"
	print "displaceParticle2()"
	sys.stdout.flush()

	estPose0 = prevPose0
	estPose1 = prevPose1
	estPose2 = initPose2
	estPose3 = initPose3

	nodeID2 = nodeID3-1
	nodeID1 = nodeID2-1
	nodeID0 = nodeID1-1

	nodeID = nodeID3

	medial0 = poseData.medialAxes[nodeID0]
	medial1 = poseData.medialAxes[nodeID1]
	medial2 = poseData.medialAxes[nodeID2]
	medial3 = poseData.medialAxes[nodeID3]

	medialSpline0 = SplineFit(medial0, smooth=0.1)
	medial0_vec = medialSpline0.getUniformSamples()
	medialSpline1 = SplineFit(medial1, smooth=0.1)
	medial1_vec = medialSpline1.getUniformSamples()

	medialSpline2 = SplineFit(medial2, smooth=0.1)
	medial2_vec = medialSpline2.getUniformSamples()
	medialSpline3 = SplineFit(medial3, smooth=0.1)
	medial3_vec = medialSpline3.getUniformSamples()

	direction = poseData.travelDirs[nodeID]

	currPose2 = estPose2
	currPose3 = estPose3
	currProb2 = 0.0
	currProb3 = 0.0

	foreAngle0, backAngle0 = getTipAngles(medial0_vec, estPose0)
	foreAngle1, backAngle1 = getTipAngles(medial1_vec, estPose1)



	#splicePaths = list(set(pathSplices2 + pathSplices3))
	splicePaths = pathSplices2 + pathSplices3

	if nodeID > 0:
		
		currPose3 = getInPlaceGuess(poseData, nodeID-1, nodeID, currPose2, estPose3, [], direction)

		if False:
			if nodeID % 2 == 1:

				distEst = 0.5
				
				" the guess process gives a meaningful guess for these node's poses "
				" before this, it is meaningless "
				currPose2 = getStepGuess(poseData, nodeID-3, nodeID-1, estPose0, estPose2, direction, distEst = distEst)

				#supportLine = mapHyp.paths[0]

				#currPose3 = getInPlaceGuess(poseData, nodeID-1, nodeID, currPose2, estPose3, supportLine, direction)
				currPose3 = getInPlaceGuess(poseData, nodeID-1, nodeID, currPose2, estPose3, [], direction)

				print "displace old to new poses:", estPose0, estPose2, estPose3, currPose2, currPose3

			else:
				print "movePath:  DO NOTHING"

		#return currPose2, currPose3


		""" if the node has a spatial feature, localize it to the closest landmark """

		""" if the node has a landmark feature, then we try to localize on the junction point """
		""" FIXME: only use the first long path, we don't try and figure out the longest yet """
		currPoses = {nodeID2: currPose2, nodeID3: currPose3}
		poseOffsets = {nodeID2: [0.0,0.0,0.0], nodeID3: [0.0,0.0,0.0]}
	
		""" if the node has been localized to a landmark, we need the updated version """
		currPose2 = currPoses[nodeID2]
		currPose3 = currPoses[nodeID3]

		print "final poses:", currPose2, currPose3
		currProb2 = 0.0
		currProb3 = 0.0


		
		" now we fit the guessed pose to a splice of the paths "
		" this constrains it to the previously explored environment "
		" all possible splices are used, and they are filtered to determine the best result "

		if poseData.numNodes >= 4 and nodeID % 2 == 1:
			
			" 1) get spliced path that the previous node is on "
			#medial2 = medialAxis0
			#medial3 = medialAxis1
			
			#initPose2 = initPose0
			#initPose3 = initPose1
			
			#print "junctions:", junctions
			print "initPose2:", initPose2
			print "initPose3:", initPose3
									
			#splicePaths = []

			resultMoves2 = []
			resultMoves3 = []

			#for path in splicePaths:		
			for spliceIndex in range(len(splicePaths)):		

				path = splicePaths[spliceIndex]

				" 2) get pose of previous node, get point on path curve "
				#pose0 = mapHyp.nodePoses[nodeID-3]
				#pose1 = mapHyp.nodePoses[nodeID-2]
				pose0 = prevPose0
				pose1 = prevPose1

				" FIXME:  make sure path is oriented correctly wrt node medial axis "
				#medial0 = prevMedialAxis0
				
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
				
				#pose2 = initPose2
				#pose3 = initPose3
				pose2 = currPose2
				pose3 = currPose3
				
				print "pose2,pose3:", pose2, pose3
				
				" 4) set as pose of new node "
				#mapHyp.nodePoses[nodeID-1] = pose2
				#mapHyp.nodePoses[nodeID] = pose3

				uPath2, uMedialOrigin2 = selectLocalCommonOrigin(orientedSplicePath, medial2, pose2)

				uPath3, uMedialOrigin3 = selectLocalCommonOrigin(orientedSplicePath, medial3, pose3)
				
				u2 = uPath2
				u3 = uPath3

				print "input: uMedialOrigin2, u2, pose2:", uMedialOrigin2, u2, pose2
				print "input: uMedialOrigin3, u3, pose3:", uMedialOrigin3, u3, pose3

				
				resultPose2, lastCost2, matchCount2, currAng2, currU2 = gen_icp.globalPathToNodeOverlapICP2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial2, plotIter = False, n1 = nodeID-1, n2 = -1, arcLimit = 0.01)
				resultPose3, lastCost3, matchCount3, currAng3, currU3 = gen_icp.globalPathToNodeOverlapICP2([u3, uMedialOrigin3, 0.0], orientedSplicePath, medial3, plotIter = False, n1 = nodeID, n2 = -1, arcLimit = 0.01)
				
				print "resultPoses:", resultPose2, resultPose3

				multiDepCount = 0

				#self.mapStateID = mapStateID
				result2 = getMultiDeparturePoint(orientedSplicePath, medial2_vec, pose2, resultPose2, [], nodeID-1, hypID = 0, pathPlotCount = 0, particleIndex=particleIndex, spliceIndex=spliceIndex, plotIter = False)
				multiDepCount += 1
				result3 = getMultiDeparturePoint(orientedSplicePath, medial3_vec, pose3, resultPose3, [], nodeID, hypID = 0, pathPlotCount = 0, particleIndex=particleIndex, spliceIndex=spliceIndex, plotIter = False)
				multiDepCount += 1
				
				#results1.append(result+(k,))
				" (departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 )"
				
				" (resultPose2,lastCost2,matchCount2,fabs(currAng2)) "

				contigFrac_2 = result2[12]
				contigFrac_3 = result3[12]
				overlapSum2 = result2[13]
				overlapSum3 = result3[13]

				foreAngle2, backAngle2 = getTipAngles(medial2_vec, resultPose2)
				foreAngle3, backAngle3 = getTipAngles(medial3_vec, resultPose3)

				angDiff2 = min(abs(diffAngle(foreAngle2,foreAngle0)), abs(diffAngle(backAngle2,backAngle0)))
				angDiff3 = min(abs(diffAngle(foreAngle3,foreAngle1)), abs(diffAngle(backAngle3,backAngle1)))
				
				#angDiff2 = abs(diffAngle(pose2[2],resultPose2[2]))
				#angDiff3 = abs(diffAngle(pose3[2],resultPose3[2]))

				if overlapSum2 > 1e10:
					newProb2 = 0.0	
				else:
					newProb2 = (pi-angDiff2) * contigFrac_2

				if overlapSum3 > 1e10:
					newProb3 = 0.0	
				else:
					newProb3 = (pi-angDiff3) * contigFrac_3
				
				print "angDiff:", angDiff2, angDiff3
				print "overlapSum:", overlapSum2, overlapSum3
				print "newProb:", newProb2, newProb3
				


				" NOTE:  added minimum threshold to angle difference "
				" NOTE:  guess pose is now the inter-nodal estimate instead of path-based estimate "
				if angDiff2 < 0.5:
					#resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(currAng2)) + result2)
					resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(0.0)) + result2 + (orientedSplicePath, newProb2))
				if angDiff3 < 0.5:
					#resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(currAng3)) + result3)				
					resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(0.0)) + result3 + (orientedSplicePath, newProb3))				
			
			
			
			""" sort by angDiff follwed by contigFrac """
			resultMoves2 = sorted(resultMoves2, key=itemgetter(18))
			resultMoves2 = sorted(resultMoves2, key=itemgetter(16), reverse=True)

			resultMoves3 = sorted(resultMoves3, key=itemgetter(18))
			resultMoves3 = sorted(resultMoves3, key=itemgetter(16), reverse=True)
			
			#print "resultMoves2:"
			#for res in resultMoves2:
			#	print res
			#print "resultMoves3:"
			#for res in resultMoves3:
			#	print res

			currSplice2 = []
			currSplice3 = []

			" select the best pose for each node "
			" record the splice used for the fit "
			if len(resultMoves2) > 0:
				currPose2 = resultMoves2[0][0]
				currSplice2 = resultMoves2[0][19]
				currProb2 = resultMoves2[0][20]

			else:
				print "node", nodeID-1, "not movePathed because no valid pose"
				currProb2 = 0.0

			if len(resultMoves3) > 0:
				currPose3 = resultMoves3[0][0]
				currSplice3 = resultMoves3[0][19]
				currProb3 = resultMoves3[0][20]

			else:
				print "node", nodeID, "not movePathed because no valid pose"
				currProb3 = 0.0

			if len(currSplice2) == 0 and len(currSplice3) > 0:
				currSplice2 = currSplice3

			if len(currSplice3) == 0 and len(currSplice2) > 0:
				currSplice3 = currSplice2

			" default to the 1st splice, no particular meaning "
			" FIXME: default to previous pose's splice configuration "
			if len(currSplice2) == 0:
				currSplice2 = splicePaths[0]
			if len(currSplice3) == 0:
				currSplice3 = splicePaths[0]

			" move the pose particles along their paths "	


			print "fitted poses:", currPose2, currPose3
			print "currProbs:", currProb2, currProb3


		return currPose2, currPose3, currProb2, currProb3


def __remote_prof_multiParticle(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_multiParticle(rank, qin, qout)', "particle_%d.prof" % pid )
		cProfile.runctx("__remote_multiParticle(rank, qin, qout)", globals(), locals(), "particle_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


def __remote_multiParticle(rank, qin, qout):


	try:
		#sys.stdout = open("particle_" + str(os.getpid()) + ".out", "w")
		#sys.stderr = open("particle_" + str(os.getpid()) + ".err", "w")
		sys.stdout = open("particle_" + str(rank) + ".out", "a")
		sys.stderr = open("particle_" + str(rank) + ".err", "a")
		print 'module name:', __name__
		print 'parent process:', os.getppid()
		print 'process id:', os.getpid()

		print "started __remote_multiParticle"

		while 1:
			# read input queue (block until data arrives)
			results = []
			nc, args = qin.get()
			#print rank, "received", nc, args
			# process data
			#knn = __do_nothing(data, nc, someArg2, someArg3)

			#foo = None
			#badVal = foo[0] 

			for job in args:

				#localizeJobs.append([pathU0, oldMedialU0, 0.0, pathU1, oldMedialU1, 0.0, spliceIndex, deepcopy(orientedPath0), deepcopy(medial0), deepcopy(medial1), deepcopy(hypPose0), deepcopy(hypPose1), [], nodeID0, particleIndex, updateCount, self.hypothesisID])

				branchIndex = job[6]
				spliceIndex = job[7]
				#globalPath = job[8]
				orientedPath0 = job[8]
				medial0 = job[9]
				medial1 = job[10]
				initPose0 = job[11]
				initPose1 = job[12]

				prevMedial0 = job[13]
				prevMedial1 = job[14]
				prevPose0 = job[15]
				prevPose1 = job[16]

				pathIDs = job[17]
				nodeID0 = job[18]
				nodeID1 = job[19]
				particleIndex = job[20]
				updateCount = job[21]
				hypID = job[22]
				branchProbVal = job[23]
				landmarks_G = job[24]
				landmark0_N = job[25]
				landmark1_N = job[26]

				oldMedialP0 = job[0]
				oldMedialU0 = job[1]
				angGuess0 = job[2]
				oldMedialP1 = job[3]
				oldMedialU1 = job[4]
				angGuess1 = job[5]

				poseOrigin = Pose(initPose0)
				
				globalMedial0 = []
				for p in medial0:
					globalMedial0.append(poseOrigin.convertLocalOffsetToGlobal(p))
				
				#globalMedial1 = []
				#for p in medial1:
				#	globalMedial1.append(poseOrigin.convertLocalOffsetToGlobal(p))


				globalMedialP0 = poseOrigin.convertLocalOffsetToGlobal(oldMedialP0)	
				globalMedialP1 = poseOrigin.convertLocalOffsetToGlobal(oldMedialP1)	

				#orientedPath0 = orientPathLean(globalPath, globalMedial0)
				orientedPathSpline0 = SplineFit(orientedPath0, smooth=0.1)
				pathU0 = orientedPathSpline0.findU(globalMedialP0)
				pathU1 = orientedPathSpline0.findU(globalMedialP1)


				params0 = [pathU0, oldMedialU0, angGuess0]
				params1 = [pathU1, oldMedialU1, angGuess1]

				result = multiParticleFitSplice(params0, params1, orientedPath0, medial0, medial1, initPose0, initPose1, prevMedial0, prevMedial1, prevPose0, prevPose1, pathIDs, nodeID0, nodeID1, landmarks_G, landmark0_N, landmark1_N, particleIndex, hypID = hypID, pathPlotCount = updateCount, branchIndex = branchIndex, spliceIndex = spliceIndex, branchProbVal = branchProbVal)
				results.append(result)
							   
			# write to output queue
			qout.put((nc,results))
	except:
		print "Worker process failed. Exiting"
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]
		
		sys.stdout.flush()
		sys.stderr.flush()

		qout.put((None,None))
		raise


def batchLocalizeParticle(localizeJobs):

	global renderGlobalPlotCount
	global pool_posePart
	global qin_posePart
	global qout_posePart

	ndata = len(localizeJobs)
	
	args = []
	
	for k in range(ndata):
		arg = localizeJobs[k]
		args.append(arg)
	
	nproc = __num_processors()
	
	#nproc *= 2
	print "nproc =", nproc
	
	# compute chunk size
	chunk_size = ndata / nproc
	chunk_size = 2 if chunk_size < 2 else chunk_size
	
	print "chunk_size =", chunk_size
	print "max_size =", ndata/chunk_size
	
	# set up a pool of processes
	if len(pool_posePart) == 0:

		qin_posePart = processing.Queue(maxsize=ndata/chunk_size)
		qout_posePart = processing.Queue(maxsize=ndata/chunk_size)
		pool_posePart = [processing.Process(target=__remote_multiParticle,
		#pool_posePart = [processing.Process(target=__remote_prof_multiParticle,
				#args=(rank, qin_posePart, qout_posePart, splices, medial, initPose, pathIDs, nodeID))
				args=(rank, qin_posePart, qout_posePart))
					for rank in range(nproc)]
		for p in pool_posePart: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	#print "args =", args
	while 1:
		_data = args[cur:cur+chunk_size]
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		qin_posePart.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	isFail = False
	while len(knn) < nc:
		thisKnn = qout_posePart.get()
		if thisKnn[0] != None:
			knn += [thisKnn]
		else:
			isFail = True
			break
	
	if isFail:
		qin_posePart.close()
		qin_posePart.join_thread()

		qout_posePart.close()
		qout_posePart.cancel_join_thread()

		for p in pool_posePart:
			p.terminate()

		raise

	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	for p in pool_posePart:
		p.terminate()
	pool_posePart = []

	print "returning"
	return knn



def multiParticleFitSplice(initGuess0, initGuess1, orientedPath, medialAxis0, medialAxis1, initPose0, initPose1, prevMedialAxis0, prevMedialAxis1, prevPose0, prevPose1, pathIDs, nodeID0, nodeID1, landmarks_G, landmark0_N, landmark1_N, particleIndex, hypID = 0, pathPlotCount = 0, branchIndex = None, spliceIndex = 0, branchProbVal = 1.0):

	print "multiParticleFitSplice()"
	sys.stdout.flush()


	u1 = initGuess0[0]
	u2 = initGuess0[1]
	angGuess = initGuess0[2]

	resultPose0, lastCost0, matchCount0, currAng0, currU0 = gen_icp.globalPathToNodeOverlapICP2([u1, u2, angGuess], orientedPath, medialAxis0, plotIter = False, n1 = nodeID0, n2 = -1, arcLimit = 0.5, origPose = initPose0)

	u1 = initGuess1[0]
	u2 = initGuess1[1]
	angGuess = initGuess1[2]
	resultPose1, lastCost1, matchCount1, currAng1, currU1 = gen_icp.globalPathToNodeOverlapICP2([u1, u2, angGuess], orientedPath, medialAxis1, plotIter = False, n1 = nodeID1, n2 = -1, arcLimit = 0.5, origPose = initPose1)


	resultArgs0 = getMultiDeparturePoint(orientedPath, medialAxis0, initPose0, resultPose0, pathIDs, nodeID0, pathPlotCount = pathPlotCount, hypID = hypID, particleIndex = particleIndex, spliceIndex = spliceIndex, plotIter = False)

	resultArgs1 = getMultiDeparturePoint(orientedPath, medialAxis1, initPose1, resultPose1, pathIDs, nodeID1, pathPlotCount = pathPlotCount, hypID = hypID, particleIndex = particleIndex, spliceIndex = spliceIndex, plotIter = False)
	#resultArgs = ([0.0, 0.0],  0.0, False, False, 0.0, 0.0, [0.0, 0.0], 0.0, False, False, 0.0, 0.0, 1.0, 0.0, 0.0)

	isExist1_0 = resultArgs0[3]
	isExist2_0 = resultArgs0[9]
	contigFrac0 = resultArgs0[12]

	isExist1_1 = resultArgs1[3]
	isExist2_1 = resultArgs1[9]
	contigFrac1 = resultArgs1[12]

	currPoseOrigin = Pose(resultPose0)
	prevPoseOrigin = Pose(prevPose0)

	print "curr/prev pose:", resultPose0, prevPose0

	currGlobalMedial0 = []
	for p in medialAxis0:
		currGlobalMedial0.append(currPoseOrigin.convertLocalOffsetToGlobal(p))

	if prevMedialAxis0 != None:
		prevGlobalMedial0 = []
		for p in prevMedialAxis0:
			prevGlobalMedial0.append(prevPoseOrigin.convertLocalOffsetToGlobal(p))

		#currSum = getCurveOverlap(currGlobalMedial0, prevGlobalMedial0, plotIter = True, overlapPlotCount = pathPlotCount*20+particleIndex)
		currSum = getCurveOverlap(currGlobalMedial0, prevGlobalMedial0, plotIter = True, overlapPlotCount=nodeID0)
	else:
		currSum = None

	" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 "
	icpDist0 = sqrt((resultPose0[0]-initPose0[0])**2 + (resultPose0[1]-initPose0[1])**2)
	icpDist1 = sqrt((resultPose1[0]-initPose1[0])**2 + (resultPose1[1]-initPose1[1])**2)

	utilVal0 = (1.0-contigFrac0) + (isExist1_0 or isExist2_0) + (1.0-contigFrac1) + (isExist1_1 or isExist2_1)
	#utilVal0 = (1.0-contigFrac0) + (1.0-contigFrac1)

	frame0 = Pose(resultPose0)
	frame1 = Pose(resultPose1)

	if prevMedialAxis0 != None:
		foreAngle0, backAngle0 = getTipAngles(prevMedialAxis0, prevPose0)
		foreAngle2, backAngle2 = getTipAngles(medialAxis0, resultPose0)
		foreAngle1, backAngle1 = getTipAngles(prevMedialAxis1, prevPose1)
		foreAngle3, backAngle3 = getTipAngles(medialAxis1, resultPose1)

		angDiff0 = min(abs(diffAngle(foreAngle2,foreAngle0)), abs(diffAngle(backAngle2,backAngle0)))
		angDiff1 = min(abs(diffAngle(foreAngle3,foreAngle1)), abs(diffAngle(backAngle3,backAngle1)))
	else:
		angDiff0 = 0.0
		angDiff1 = 0.0
		
	#LANDMARK_THRESH = 1e100
	#LANDMARK_THRESH = 3.0
	LANDMARK_THRESH = 7.0
	#LANDMARK_THRESH = 4.5
	CLOSE_THRESH = 0.3
	poseSum = 0.0
	landmark0_G = None
	landmark1_G = None
	thresh0 = None
	thresh1 = None


	if landmark0_N != None:
		landmark0_G = (frame0.convertLocalToGlobal(landmark0_N[0]), landmark0_N[1], landmark0_N[2])
		thresh0 = landmark0_G[1]

	if landmark1_N != None:
		landmark1_G = (frame1.convertLocalToGlobal(landmark1_N[0]), landmark1_N[1], landmark1_N[2])
		thresh1 = landmark1_G[1]

	print hypID, nodeID0, particleIndex, spliceIndex, "landmarks:", landmark0_G, landmark1_G, landmarks_G

	for i in range(len(landmarks_G)):
		p1 = landmarks_G[i][0]
		threshG = landmarks_G[i][1]

		if landmark0_G != None:

			maxThresh = thresh0
			if threshG > thresh0:
				maxThresh = threshG

			dist0 = sqrt((p1[0]-landmark0_G[0][0])**2 + (p1[1]-landmark0_G[0][1])**2)
			if dist0 < LANDMARK_THRESH:

				#poseSum += sqrt(dist0*dist0/(maxThresh*maxThresh))
				poseSum += sqrt(dist0*dist0/(threshG*threshG + thresh0*thresh0))

				#if dist0 > maxThresh:
				#	poseSum += 10.0*dist0*dist0
				#else:
				#	poseSum += dist0*dist0

		if landmark1_G != None:

			maxThresh = thresh1
			if threshG > thresh1:
				maxThresh = threshG

			dist1 = sqrt((p1[0]-landmark1_G[0][0])**2 + (p1[1]-landmark1_G[0][1])**2)
			if dist1 < LANDMARK_THRESH:
				#poseSum += sqrt(dist1*dist1/(maxThresh*maxThresh))
				
				poseSum += sqrt(dist1*dist1/(threshG*threshG + thresh1*thresh1))

				#if dist1 > maxThresh:
				#	poseSum += 10.0*dist1*dist1
				#else:
				#	poseSum += dist1*dist1

	#utilVal0 = utilVal0 * (1.0-branchProbVal)

	#return (particleIndex, icpDist0, resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs0 + (isExist1_0 or isExist2_0,)
	return (particleIndex, icpDist0, resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs0 + (isExist1_0 or isExist2_0,) + (icpDist1, resultPose1, lastCost1, matchCount1, currAng1, currU1) + resultArgs1 + (isExist1_1 or isExist2_1,) + (utilVal0, branchIndex, spliceIndex, initPose0, initPose1, currSum, poseSum, angDiff0, angDiff1)


class Particle:

	def __init__(self, pose0, pose1, pathID, hypDist, weightVal, mapStateID):
		
		" converted to arc distance on-the-fly as needed "
		" remains invariant to regenerating paths "
		#self.pose0 = (0.0,0.0,0.0)
		#self.pose1 = (0.0,0.0,0.0)
		self.pose0 = deepcopy(pose0)
		self.pose1 = deepcopy(pose1)

		self.prevPose0 = deepcopy(pose0)
		self.prevPose1 = deepcopy(pose1)

		self.dispPose0 = deepcopy(pose0)
		self.dispPose1 = deepcopy(pose1)

		" classification or correspondence as member of path "
		self.memberPaths = []
		self.memberPaths.append(pathID)

		" junction description data "
		self.junctionData = {}

		" a currently undiagnosed necessary evil. "
		" do I make splice a correspondence problem as well? "
		" pose0 and pose1 splices may be different "
		self.spliceName = ('t0', 't1')
		#self.spliceName = spliceName
		self.spliceCurve = None

		self.junctionData[0] = {
									"parentID" : None,
									"branchNodeID" : None,
									"localJunctionPose" : None, 
									"globalJunctionPose" : None,
									"nodeSet" : [],
									"localNodePoses" : {},
									"controlPose" : [0.0,0.0,0.0],
									"probDist" : [ 1.0 ],
									"branchPoseDist" : None,
									"controlPoseDist" : None,
									"arcDists" : None,
									"maxLikelihoodBranch" : 0
									}		

		self.addNode(0,0, pose0)
		self.addNode(1,0, pose1)

		" temporary values, computed on-the-fly for each evaluation "

		self.displaceProb0 = 0.0
		self.displaceProb1 = 0.0

		self.landmarkSum = 0.0

		self.hypDist = hypDist
		self.weightVal = weightVal
		self.mapStateID = mapStateID
		self.branchArcDists = []
		self.branchControls = []

		self.maxLikelihoodBranch = None

			

	def displacePose(self, pose0, pose1):

		self.prevPose0 = self.pose0
		self.prevPose1 = self.pose1

		self.dispPose0 = deepcopy(pose0)
		self.dispPose1 = deepcopy(pose1)

		self.pose0 = deepcopy(pose0)
		self.pose1 = deepcopy(pose1)

	def getControlPoses(self):

		controlPoses = {}
		allPathIDs = self.junctionData.keys()
		for pathID in allPathIDs:
			controlPoses[pathID] = self.junctionData[pathID]["controlPose"]

		return controlPoses

	def addPath(self, pathID, parentID, branchNodeID, localNodePose, localDivergencePose, globalJunctionPose, localJunctionPose, controlPose, numBranches, arcDists, controlPoses):

		self.junctionData[pathID] = {
									"parentID" : parentID,
									"branchNodeID" : branchNodeID,
									"localDivergencePose" : localDivergencePose, 
									"localJunctionPose" : localJunctionPose, 
									"globalJunctionPose" : globalJunctionPose,
									"nodeSet" : [branchNodeID,],
									"localNodePoses" : {branchNodeID : localNodePose},
									"controlPose" : controlPose,
									"probDist" : [ 1.0 / float(numBranches) for k in range(numBranches) ],
									"branchPoseDist" : [deepcopy(globalJunctionPose) for k in range(numBranches)],
									"controlPoseDist" : controlPoses,
									"arcDists" : [deepcopy(arcDists)],
									"maxLikelihoodBranch" : 0
									}		

	def delNode(self, nodeID, pathID):

		self.junctionData[pathID]["nodeSet"].remove(nodeID)
		del self.junctionData[pathID]["localNodePoses"][nodeID]
		
	def addNode(self, nodeID, pathID, localNodePose):

		self.junctionData[pathID]["nodeSet"].append(nodeID)
		self.junctionData[pathID]["localNodePoses"][nodeID] = localNodePose

	def copy(self):
		newParticle = Particle(deepcopy(self.pose0), deepcopy(self.pose1), 0, self.hypDist, self.weightVal, self.mapStateID)
		newParticle.memberPaths = deepcopy(self.memberPaths)
		newParticle.junctionData = deepcopy(self.junctionData)

		newParticle.spliceName = deepcopy(self.spliceName)
		newParticle.spliceCurve = deepcopy(self.spliceCurve)

		newParticle.prevPose0 = deepcopy(self.prevPose0)
		newParticle.prevPose1 = deepcopy(self.prevPose1)

		newParticle.dispPose0 = deepcopy(self.dispPose0)
		newParticle.dispPose1 = deepcopy(self.dispPose1)

		newParticle.branchArcDists = deepcopy(self.branchArcDists)
		newParticle.branchControls = deepcopy(self.branchControls)
		newParticle.maxLikelihoodBranch = deepcopy(self.maxLikelihoodBranch)

		return newParticle

	def __str__(self):
		return str(self.__dict__)

	def __eq__(self, other): 
		resultVal = self.__dict__ == other.__dict__
		#print "compare equals", resultVal, self.__dict__, other.__dict__
		return resultVal
		#return self.__dict__ == other.__dict__

