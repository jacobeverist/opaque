import multiprocessing as processing
import random
from copy import deepcopy
from SplineFit import SplineFit
import pylab
from math import sqrt
from Pose import Pose
import gen_icp
from math import acos, asin, fabs
from numpy import array
import ctypes, os, sys
from operator import itemgetter

from functions import diffAngle
from Splices import getMultiDeparturePoint, orientPath, getCurveOverlap
from MapProcess import getStepGuess, getInPlaceGuess, selectLocalCommonOrigin

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

def batchDisplaceParticles(displaceJobs):

	global renderGlobalPlotCount
	global pool_dispPosePart
	global qin_dispPosePart
	global qout_dispPosePart

	ndata = len(displaceJobs)
	
	args = []
	
	for k in range(ndata):
		arg = displaceJobs[k]
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
	if len(pool_dispPosePart) == 0:

		qin_dispPosePart = processing.Queue(maxsize=ndata/chunk_size)
		qout_dispPosePart = processing.Queue(maxsize=ndata/chunk_size)
		pool_dispPosePart = [processing.Process(target=__remote_displaceParticle,
		#pool_dispPosePart = [processing.Process(target=__remote_prof_displaceParticle,
				#args=(rank, qin_dispPosePart, qout_dispPosePart, splices, medial, initPose, pathIDs, nodeID))
				args=(rank, qin_dispPosePart, qout_dispPosePart))
				#args=(rank, qin_dispPosePart, qout_dispPosePart))
					for rank in range(nproc)]
		for p in pool_dispPosePart: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	#print "args =", args
	while 1:
		_data = args[cur:cur+chunk_size]
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		qin_dispPosePart.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout_dispPosePart.get()]
	
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "returning"
	return knn

def __remote_prof_displaceParticle(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_displaceParticle(rank, qin, qout)', "particle_%d.prof" % pid )
		cProfile.runctx("__remote_displaceParticle(rank, qin, qout)", globals(), locals(), "displaceParticle_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


def __remote_displaceParticle(rank, qin, qout):

	sys.stdout = open("displaceParticle_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("displaceParticle_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

	print "started __remote_displaceParticle"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()

		results = []
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

			result = displaceParticle( poseData, pathSplices2, pathSplices3, supportLine, nodeID3, initPose2, initPose3, prevPose0, prevPose1)
			results.append((particleIndex,) + result)
						   
		# write to output queue
		qout.put((nc,results))

def displaceParticle( poseData, pathSplices2, pathSplices3, supportLine, nodeID3, initPose2, initPose3, prevPose0, prevPose1):

	#print "movePath(", nodeID, ",", direction, ",", distEst, ")"

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


	direction = poseData.travelDirs[nodeID]

	currPose2 = estPose2
	currPose3 = estPose3


	#splicePaths = list(set(pathSplices2 + pathSplices3))
	splicePaths = pathSplices2 + pathSplices3

	if nodeID > 0:
		
		if nodeID % 2 == 1:

			rState = random.getstate()
			random.seed()	
			moveChance = random.random()
			moveDist = random.random()
			
			if moveChance >= 0.1:
				#distEst = 0.5 + (0.5-moveDist)/2.0
				distEst = random.gauss(0.7,0.6)
			else:
				#distEst = -0.5 - (0.5-moveDist)/2.0
				distEst = random.gauss(-0.5,0.6)

			random.setstate(rState)

			#distEst = moveChance * 2.0
			#distEst = 2.0*moveChance - 0.4

			print "displacing", moveChance, distEst

			" the guess process gives a meaningful guess for these node's poses "
			" before this, it is meaningless "
			currPose2 = getStepGuess(poseData, nodeID-3, nodeID-1, estPose0, estPose2, direction, distEst = distEst)

			#supportLine = mapHyp.paths[0]

			currPose3 = getInPlaceGuess(poseData, nodeID-1, nodeID, currPose2, estPose3, supportLine, direction)



		else:
			print "movePath:  DO NOTHING"

		#return currPose2, currPose3
		
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

			for path in splicePaths:		

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
					resultMoves2.append((resultPose2,lastCost2,matchCount2,fabs(0.0)) + result2 + (orientedSplicePath,))
				if angDiff3 < 0.5:
					#resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(currAng3)) + result3)				
					resultMoves3.append((resultPose3,lastCost3,matchCount3,fabs(0.0)) + result3 + (orientedSplicePath,))				
			
			
			
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

			currSplice2 = []
			currSplice3 = []

			" select the best pose for each node "
			" record the splice used for the fit "
			if len(resultMoves2) > 0:
				currPose2 = resultMoves2[0][0]
				currSplice2 = resultMoves2[0][19]
			else:
				print "node", nodeID-1, "not movePathed because no valid pose"

			if len(resultMoves3) > 0:
				currPose3 = resultMoves3[0][0]
				currSplice3 = resultMoves3[0][19]
			else:
				print "node", nodeID, "not movePathed because no valid pose"

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


			"FIXME:  choose a splice that is common to the paired old and new poses "
			orientedPathSpline2 = SplineFit(currSplice2, smooth=0.1)
			orientedPathSpline3 = SplineFit(currSplice3, smooth=0.1)

			#estPose0 = prevPose0
			#estPose1 = prevPose1
			#estPose2 = initPose2
			#estPose3 = initPose3

			oldPose0 = estPose0
			minDist0, oldU0, oldP0 = orientedPathSpline2.findClosestPoint(oldPose0)

			oldPose1 = estPose1
			minDist1, oldU1, oldP1 = orientedPathSpline3.findClosestPoint(oldPose1)

			newPose2 = currPose2
			minDist2, newU2, newP2 = orientedPathSpline2.findClosestPoint(newPose2)

			newPose3 = currPose3
			minDist3, newU3, newP3 = orientedPathSpline3.findClosestPoint(newPose3)

			travelDist2 = orientedPathSpline2.dist(oldU0, newU2)
			if oldU0 > newU2:
				travelDist2 = -travelDist2
			travelDist3 = orientedPathSpline3.dist(oldU1, newU3)
			if oldU1 > newU3:
				travelDist3 = -travelDist3
			
			orientedPoints2 = orientedPathSpline2.getUniformSamples()
			orientedPoints3 = orientedPathSpline3.getUniformSamples()

			#mapHyp.displacePoseParticles(nodeID-1, nodeID, travelDist2, travelDist3)

			medialSpline = SplineFit(medial2, smooth=0.1)
			#medial3

			currSplice2

			poseOrigin2 = Pose(currPose2)
			globalMedial2 = []
			for p in medial2:
				globalMedial2.append(poseOrigin2.convertLocalToGlobal(p))

			orientedSplicePath = orientPath(currSplice2, globalMedial2)				
			pathSpline = SplineFit(orientedSplicePath, smooth=0.1)


			" original pose, before displacement "
			minDist2, oldU2, oldP2 = pathSpline.findClosestPoint(estPose0)
			minDist3, oldU3, oldP3 = pathSpline.findClosestPoint(estPose1)

			thisDist2 = travelDist2
			thisDist3 = travelDist3

			#moveChance = random.random()
			" 80% chance we move forward, 20% we move backward "    
			#if moveChance >= 0.1:
			#	thisDist2 = travelDist2
			#	thisDist3 = travelDist3
			#else:
			#	thisDist2 = -travelDist2
			#	thisDist3 = -travelDist3

			thisDist2 = travelDist2
			thisDist3 = travelDist3
			#thisDist2 = random.gauss(thisDist2,0.6)
			#thisDist3 = random.gauss(thisDist3,0.6)


			newU2 = pathSpline.getUOfDist(oldU2, thisDist2)
			newU3 = pathSpline.getUOfDist(oldU3, thisDist3)

			newDist2 = pathSpline.dist_u(newU2)
			newP2 = pathSpline.point_u(newU2)
			currPose2 = deepcopy(newP2)

			" choose angle of fit to selected splice "
			currPose2[2] = newPose2[2]

			newDist3 = pathSpline.dist_u(newU3)
			newP3 = pathSpline.point_u(newU3)
			currPose3 = deepcopy(newP3)

			" choose angle of fit to selected splice "
			currPose3[2] = newPose3[2]


		return currPose2, currPose3

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

	sys.stdout = open("particle_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("particle_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

	print "started __remote_multiParticle"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		#print rank, "received", nc, args
		# process data
		#knn = __do_nothing(data, nc, someArg2, someArg3)
		results = []
		for job in args:

			#localizeJobs.append([pathU0, oldMedialU0, 0.0, pathU1, oldMedialU1, 0.0, spliceIndex, deepcopy(orientedPath0), deepcopy(medial0), deepcopy(medial1), deepcopy(hypPose0), deepcopy(hypPose1), [], nodeID0, particleIndex, updateCount, self.hypothesisID])

			branchIndex = job[6]
			spliceIndex = job[7]
			globalPath = job[8]
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

			oldMedialP0 = job[0]
			oldMedialU0 = job[1]
			angGuess0 = job[2]
			oldMedialP1 = job[3]
			oldMedialU1 = job[4]
			angGuess1 = job[5]

			poseOrigin = Pose(initPose0)
			
			globalMedial0 = []
			for p in medial0:
				globalMedial0.append(poseOrigin.convertLocalToGlobal(p))
			
			#globalMedial1 = []
			#for p in medial1:
			#@	globalMedial1.append(poseOrigin.convertLocalToGlobal(p))


			globalMedialP0 = poseOrigin.convertLocalOffsetToGlobal(oldMedialP0)	
			globalMedialP1 = poseOrigin.convertLocalOffsetToGlobal(oldMedialP1)	

			orientedPath0 = orientPath(globalPath, globalMedial0)
			orientedPathSpline0 = SplineFit(orientedPath0, smooth=0.1)
			pathU0 = orientedPathSpline0.findU(globalMedialP0)
			pathU1 = orientedPathSpline0.findU(globalMedialP1)


			params0 = [pathU0, oldMedialU0, angGuess0]
			params1 = [pathU1, oldMedialU1, angGuess1]

			result = multiParticleFitSplice(params0, params1, orientedPath0, medial0, medial1, initPose0, initPose1, prevMedial0, prevMedial1, prevPose0, prevPose1, pathIDs, nodeID0, nodeID1, particleIndex, hypID = hypID, pathPlotCount = updateCount, branchIndex = branchIndex, spliceIndex = spliceIndex)
			results.append(result)
						   
		# write to output queue
		qout.put((nc,results))

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
	
	nproc *= 2
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
	while len(knn) < nc:
		knn += [qout_posePart.get()]
	
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "returning"
	return knn


def multiParticleFitSplice(initGuess0, initGuess1, orientedPath, medialAxis0, medialAxis1, initPose0, initPose1, prevMedialAxis0, prevMedialAxis1, prevPose0, prevPose1, pathIDs, nodeID0, nodeID1, particleIndex, hypID = 0, pathPlotCount = 0, branchIndex =0, spliceIndex = 0):


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
		currGlobalMedial0.append(currPoseOrigin.convertLocalToGlobal(p))

	if prevMedialAxis0 != None:
		prevGlobalMedial0 = []
		for p in prevMedialAxis0:
			prevGlobalMedial0.append(prevPoseOrigin.convertLocalToGlobal(p))

		#currSum = getCurveOverlap(currGlobalMedial0, prevGlobalMedial0, plotIter = True, overlapPlotCount = pathPlotCount*20+particleIndex)
		currSum = getCurveOverlap(currGlobalMedial0, prevGlobalMedial0)
	else:
		currSum = None

	" departurePoint1, angle1, isInterior1, isExist1, dist1, maxFront, departurePoint2, angle2, isInterior2, isExist2, dist2, maxBack, contigFrac, overlapSum, angDiff2 "
	icpDist0 = sqrt((resultPose0[0]-initPose0[0])**2 + (resultPose0[1]-initPose0[1])**2)
	icpDist1 = sqrt((resultPose1[0]-initPose1[0])**2 + (resultPose1[1]-initPose1[1])**2)

	utilVal0 = (1.0-contigFrac0) + (isExist1_0 or isExist2_0) + (1.0-contigFrac1) + (isExist1_1 or isExist2_1)

	#return (particleIndex, icpDist0, resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs0 + (isExist1_0 or isExist2_0,)
	return (particleIndex, icpDist0, resultPose0, lastCost0, matchCount0, currAng0, currU0) + resultArgs0 + (isExist1_0 or isExist2_0,) + (icpDist1, resultPose1, lastCost1, matchCount1, currAng1, currU1) + resultArgs1 + (isExist1_1 or isExist2_1,) + (utilVal0, branchIndex, spliceIndex, initPose0, initPose1, currSum)


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

		self.nodeCorrespondence = {}

		" temporary values, computed on-the-fly for each evaluation "

		self.hypDist = hypDist

		self.weightVal = weightVal

		self.mapStateID = mapStateID
	
	def displacePose(self, pose0, pose1):

		self.prevPose0 = self.pose0
		self.prevPose1 = self.pose1

		self.dispPose0 = deepcopy(pose0)
		self.dispPose1 = deepcopy(pose1)

		self.pose0 = deepcopy(pose0)
		self.pose1 = deepcopy(pose1)

	def addPath(self, pathID, parentID, branchNodeID, localJunctionPose, globalJunctionPose, numBranches):

		self.junctionData[pathID] = {
									"parentID" : parentID,
									"branchNodeID" : branchNodeID,
									"localJunctionPose" : localJunctionPose, 
									"globalJunctionPose" : globalJunctionPose,
									"probDist" : [ 1.0 / float(numBranches) for k in range(numBranches) ],
									"branchPoseDist" : [deepcopy(globalJunctionPose) for k in range(numBranches)],
									"maxLikelihoodBranch" : 0
									}		


	def addNode(self, nodeID, pathID):
		self.nodeCorrespondence[nodeID] = pathID

	def copy(self):
		newParticle = Particle(deepcopy(self.pose0), deepcopy(self.pose1), 0, self.hypDist, self.weightVal, self.mapStateID)
		newParticle.nodeCorrespondence = deepcopy(self.nodeCorrespondence)
		newParticle.memberPaths = deepcopy(self.memberPaths)
		newParticle.junctionData = deepcopy(self.junctionData)

		newParticle.spliceName = deepcopy(self.spliceName)
		newParticle.spliceCurve = deepcopy(self.spliceCurve)

		newParticle.prevPose0 = deepcopy(self.prevPose0)
		newParticle.prevPose1 = deepcopy(self.prevPose1)

		newParticle.dispPose0 = deepcopy(self.dispPose0)
		newParticle.dispPose1 = deepcopy(self.dispPose1)

		return newParticle

	def __str__(self):
		return str(self.__dict__)

	def __eq__(self, other): 
		resultVal = self.__dict__ == other.__dict__
		#print "compare equals", resultVal, self.__dict__, other.__dict__
		return resultVal
		#return self.__dict__ == other.__dict__

class ParticleFilter:

    def __init__(self, initPath, nodeHash):
    
        self.numParticles = 20

        " pathID, distMean "
        self.initPose = [0, 2.0]
        
        self.rootPose = [-2.0,0.0]
 
        self.nodeHash = nodeHash
        #self.distMove = 0.628
        #self.distVar = 0.0964       
        self.distMove = 1.0
        self.distVar = 0.5
        
        self.updateCount = 0
        
        self.paths = {}
        self.paths[self.updateCount] = deepcopy(initPath)


        p0 = self.paths[self.updateCount][0]
        pN = self.paths[self.updateCount][-1]
        
        dist0 = sqrt((self.rootPose[0]-p0[0])*(self.rootPose[0]-p0[0]) + (self.rootPose[1]-p0[1])*(self.rootPose[1]-p0[1]))
        distN = sqrt((self.rootPose[0]-pN[0])*(self.rootPose[0]-pN[0]) + (self.rootPose[1]-pN[1])*(self.rootPose[1]-pN[1]))

        if dist0 > distN:
            self.paths[self.updateCount].reverse()
        
        " each particle:  [pathID, distMean, distVar] "
        pathID = self.initPose[0]
        distMean = self.initPose[1]
        
        particles = []
        for i in range(self.numParticles):
            particles.append([pathID, distMean])

        self.particleSnapshots = {}
        self.particleSnapshots[self.updateCount] = particles

        self.mostLikely = 0

        self.mostLikely1 = 0
        self.mostLikely2 = 0

        self.lossCount = 0
        self.moveForward = 0
        self.moveBackward = 0
        self.avgDist = 0.0

        #self.drawState()
        
        self.resultVector = []
        
        #self.updateCount += 1
        self.colors = []
        for i in range(1000):
            self.colors.append((random.random(),random.random(),random.random()))
                    
    def update(self, newPath0, nodeID, isForward = True):

        self.lossCount = 0
        self.moveForward = 0
        self.moveBackward = 0
        self.avgDist = 0.0


        " old path + old positions => x,y "
        " x,y => new path + old positions "
        " new path + old positions + move => new path + new positions "
        
        oldPath = self.paths[self.updateCount]
        oldSpline = SplineFit(oldPath, smooth=0.1)
        

        " convert particles to x,y positions "
        xyParticles = []
        partPoints = []
        oldParticles = self.particleSnapshots[self.updateCount]
        for part in oldParticles:
            pathID = part[0]
            distMean = part[1]

            point = oldSpline.getPointOfDist(distMean)
            #print "dist:", distMean, "point:", point

            newPart = [pathID, point, distMean]
            partPoints.append(point[:2])
            xyParticles.append(newPart)

        newPath = deepcopy(newPath0)
        p0 = newPath[0]
        pN = newPath[-1]
        
        dist0 = sqrt((self.rootPose[0]-p0[0])*(self.rootPose[0]-p0[0]) + (self.rootPose[1]-p0[1])*(self.rootPose[1]-p0[1]))
        distN = sqrt((self.rootPose[0]-pN[0])*(self.rootPose[0]-pN[0]) + (self.rootPose[1]-pN[1])*(self.rootPose[1]-pN[1]))

        if dist0 > distN:
            newPath.reverse()
        
        
        " new path "
        newSpline = SplineFit(newPath, smooth = 0.1)
        
        newSpline.precompute()
        splinePoints = newSpline.densePoints
        cartPoints = []
        for p in splinePoints:
            cartPoints.append(p[:2])
        tree = cKDTree(cartPoints)
        
        qResults = tree.query(partPoints)
        distances = qResults[0]
        indices = qResults[1]

        " fit old particles onto new path and compute their distance "
        newParticles = []        
        for k in range(len(xyParticles)):
            
            #distances[k]
            #indices[k]
            
            arcDist = newSpline.distPoints[indices[k]]

            part = xyParticles[k]
            newPart = [part[0], arcDist]

            " move forward or backward "
            moveChance = random.random()
            " 80% chance we move forward, 20% we move backward "    
            if isForward:
                if moveChance >= 0.2:
                    newPart[1] += self.distMove
                    self.moveForward += 1
                else:
                    newPart[1] += -self.distMove
                    self.moveBackward += 1
            else:
                if moveChance >= 0.2:
                    newPart[1] += -self.distMove
                    self.moveBackward += 1
                else:
                    newPart[1] += self.distMove
                    self.moveForward += 1

            newPart[1] = random.gauss(newPart[1],self.distVar)

            if newPart[1] <= newSpline.dist_u(1.0) and newPart[1] >= newSpline.dist_u(0.0):  
                newParticles.append(newPart)
            else:
                self.lossCount += 1
                pass
                " add in random sample "


        " compute the medial axis for each pose "
        node1 = self.nodeHash[nodeID]

        medial1 = node1.getBestMedialAxis()
        medial1 = node1.medialLongPaths[0]
        
        #medial1 = computeHullAxis(nodeID, node1, tailCutOff = False)
        globalPath = newPath
        globalPathReverse = deepcopy(globalPath)
        globalPathReverse.reverse()
     

     
        estPose1 = node1.getGlobalGPACPose()        
        poseOrigin = Pose(estPose1)
 

        medialSpline2 = SplineFit(medial1, smooth=0.1)
        totalMedialDist = medialSpline2.dist_u(1.0)
                
        minDist, uVal, minPiont = medialSpline2.findClosestPoint([0.0,0.0])
        midMedialDist = medialSpline2.dist_u(uVal)
        
        print "node medial dist:", totalMedialDist, midMedialDist, totalMedialDist-midMedialDist, uVal, minDist
        
                
        globalMedial = []
        for p in medial1:
            globalMedial.append(poseOrigin.convertLocalToGlobal(p))
        
        

        
        
        medialSpline1 = SplineFit(globalMedial, smooth=0.1)
        globalSpline = SplineFit(globalPath, smooth=0.1)
        globalSplineReverse = SplineFit(globalPathReverse, smooth=0.1)
        medialReverse = deepcopy(medial1)
        medialReverse.reverse()




        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(globalMedial)):
            p_1 = globalMedial[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(globalPath, p_1)
            if minDist < 0.5:
                overlapMatch.append((i,j,minDist))

                pathU1 = medialSpline1.findU(p_1)    
                pathU2 = globalSpline.findU(p_2)    
                pathU2_R = globalSplineReverse.findU(p_2)    

                pathVec1 = medialSpline1.getUVector(pathU1)
                pathVec2 = globalSpline.getUVector(pathU2)
                pathVec2_R = globalSplineReverse.getUVector(pathU2_R)

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
        print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
        if angleSum1 > angleSum2:
            orientedGlobalPath = globalPathReverse
            isReverse = True
        else:
            orientedGlobalPath = globalPath
            isReverse = False
        

        #globalSplineReverse = SplineFit(globalPathReverse, smooth=0.1)

        orientedSpline = SplineFit(orientedGlobalPath, smooth=0.1)
        orientedSpline.precompute()
        cartPoints1 = []
        for p in orientedSpline.densePoints:
            cartPoints1.append(p[:2])
        treeO = cKDTree(cartPoints1)
        
        
        " sample the entire path for fitting likelihood "
        
        totalDist = orientedSpline.dist_u(1.0)
        
        """
        #totalDist = newSpline.dist_u(1.0)
        distInc = 0.1
        numInc = int(totalDist / distInc) + 1
        #uSet = [newSpline.dist_search(i*distInc)*0.001 for i in range(numInc)]
        uSet = [orientedSpline.dist_search(i*distInc)*0.001 for i in range(numInc)]

        args = []
        for i in range(len(uSet)):
    
            if isReverse:
                u1 = uSet[i]
            else:
                u1 = uSet[i]
                

            args.append([u1, 0.5, 0.0])
            
            # args.append([u1,u2,angGuess])

        print "ICP", len(args), "times"
        print args
        #exit()
        fitResults = gen_icp.batchGlobalICP(orientedGlobalPath, medial1, args)
        #fitResults = gen_icp.batchGlobalICP(newPath, medial1, args)
        print "RESULTS_SAMPLED:", len(fitResults)
 
        measResults1 = []
        queryPoints1 = []
        for k in range(len(fitResults)):
            resultPose, lastCost = fitResults[k]
            
            measResults1.append([resultPose,lastCost])
            queryPoints1.append(resultPose[:2])
        
        qResults1 = treeO.query(queryPoints1)
        distances1 = qResults1[0]
        indices1 = qResults1[1]
            
        measParticles1 = []
        for k in range(len(measResults1)):
            pos = measResults1[k][1]
            cost = measResults1[k][1]
            #arcDist = newSpline.distPoints[indices1[k]]
            arcDist = orientedSpline.distPoints[indices1[k]]
            
            if isReverse:
                arcDist = totalDist - arcDist

            print k*distInc, uSet[k], arcDist, cost

            #if measResults1[k][1] < 500:
            #   measParticles1.append([newPart[0],arcDist])
        """




                
        
        " move adjustment with ICP "       
        measResults = []
        queryPoints = []

        args = []
        for k in range(len(newParticles)):
            newPart = newParticles[k]
            
            dist = newPart[1]
            #uVal = newSpline.dist_search(dist) / 1000.

            if isReverse:
                dist = totalDist - dist
            
            uVal = orientedSpline.dist_search(dist) / 1000.
            #print dist, uVal2
            
            u2 = 0.5
            u1 = uVal
            angGuess = 0.0
            
            args.append([u1,u2,angGuess])
        
                
        results = gen_icp.batchGlobalICP(orientedGlobalPath, medial1, args)
        #print "len(results) =", len(results), "len(newParticles) =", len(newParticles)
        
        
        for k in range(len(newParticles)):

            #resultPose, lastCost = gen_icp.globalOverlapICP_GPU2(args[k], orientedGlobalPath, medial1)
            #print "estimating branch pose", nodeID, "at",  resultPose[0], resultPose[1], resultPose[2]
            
            resultPose, lastCost, matchCount = results[k]

            
            measResults.append([resultPose,lastCost, matchCount])
            queryPoints.append(resultPose[:2])
        

        qResults2 = treeO.query(queryPoints)
        distances2 = qResults2[0]
        indices2 = qResults2[1]
        
        hypParticles = []
        measParticles = []
        maxCost = 0.0
        for k in range(len(newParticles)):
            newPart = newParticles[k]
            #arcDist = newSpline.distPoints[indices2[k]]
            arcDist = orientedSpline.distPoints[indices2[k]]
            
            if isReverse:
                arcDist = totalDist - arcDist


            #print queryPoints[k], arcDist, newPart[1], measResults[k][1]
            
            hypParticles.append([measResults[k][1], measResults[k][2], distances2[k], newPart[0], arcDist])

            if measResults[k][1] > maxCost:
                maxCost = measResults[k][1]
            #if measResults[k][1] < 500:
            #    measParticles.append([newPart[0],arcDist])
        
        totalSum = 0.0
        for hyp in hypParticles:
            hyp[0] = (maxCost-hyp[0])/maxCost
            #hyp[0] = 1.0/hyp[0]
            totalSum += hyp[0]
         
        for hyp in hypParticles:
            hyp[0] = hyp[0] / totalSum
        
        hypParticles.sort(reverse=True)
        
        print "HYP particles"
        for hyp in hypParticles:
            print hyp
        
        
        measParticles = []
        weights = []
        for k in range(self.numParticles):
            
            roll = random.random()
            probTotal = 0.0
            
            partFound = False
            
            for hyp in hypParticles:
                probTotal += hyp[0]
                
                if probTotal > roll:
                    measParticles.append([hyp[3],hyp[4]])
                    weights.append(hyp[0])
                    partFound = True
                
                if partFound:
                    break
            if not partFound:
                measParticles.append([hypParticles[-1][3],hypParticles[-1][4]])
                weights.append(hyp[0])
                

        distances = []
        for hyp in measParticles:
            distances.append(hyp[1])


        #X = rand( 5, 3 )
        #X[0:5, :] *= 2

        features  = array(distances)
        features = features.reshape(len(distances),1)
        
        X = features
        Y = pdist( X )
        Z = linkage( Y )
        clf()
        labels = ['' for i in range(len(distances))]
        #labels[0] = 'label1'
        #labels[1] = 'label2'
        
        dendrogram( Z, labels=labels )
        savefig("part_dendro_%04u.png" % (self.updateCount+1))

        #resultVector = fcluster(Z,t=2, criterion='maxclust')
        self.resultVector = fcluster(Z,t=0.5, criterion='distance')
        
        #show()

        
        features = features.reshape(len(distances),1)
        #resultVector = fclusterdata(features,t=2, criterion='maxclust')
  
        print "resultCluster:", self.resultVector
 
        blueVals = []
        redVals = []
        blueAvg = 0.0
        redAvg = 0.0
        blueWAvg = 0.0
        redWAvg = 0.0
        blueTotal = 0.0
        redTotal = 0.0
        for k in range(len(distances)):
            val = features[k,0]
            clustID = self.resultVector[k]
            if clustID == 1:
                blueVals.append(val)
                blueTotal += weights[k]
                blueAvg += val
                blueWAvg += val * weights[k]
            if clustID == 2:
                redVals.append(val)
                redTotal += weights[k]
                redAvg += val
                redWAvg += val * weights[k]

        if len(redVals) > 0:
            redAvg /= len(redVals)
            redWAvg /= redTotal
            
        if len(blueVals) > 0:
            blueAvg /= len(blueVals)
            blueWAvg /= blueTotal

        " most likely "
        if len(blueVals) > len(redVals):
            self.mostLikely = [0,blueWAvg]
        else:
            self.mostLikely = [0,redWAvg]
            
        self.mostLikely1 = [0,blueWAvg]
        self.mostLikely2 = [0,redWAvg]
 
        self.avgDist = self.mostLikely[1]
        
        " generate random particles with weighted probabilities around existing particles "
        numGenParticles = self.numParticles - len(measParticles)
        numExistingParticles = len(measParticles)
        genParticles = []
        for k in range(numGenParticles):
            
            " randomly select particle out of remaining "
            selectParticleInd = random.randrange(0, numExistingParticles)
            " select its position with a variance "
            #newDist = random.gauss(measParticles[selectParticleInd][1], self.distVar)
            newDist = measParticles[selectParticleInd][1]
            
            " add particle "
            genParticles.append([measParticles[selectParticleInd][0], newDist])
        

        totalParticles = measParticles + genParticles
        
        #newSpline.findClosestPoint(resultPose)
            
            #self.makeGlobalMedialOverlapConstraint(nodeID, newPath, uVal, originU2 = 0.5)
        
        
        
        """
        " fit old particles onto new path and compute their distance "
        newParticles = []        
        for part in xyParticles:
            #print "part:", part
            pos = part[1]
            minDist, minU, closestPoint = newSpline.findClosestPoint(pos)
            arcDist = newSpline.dist(0.0, minU)

            #print "oldDist,newDist,pos,minU,closestPoint:", part[2], arcDist, pos, minU, closestPoint

            newPart = [part[0], arcDist]

            " move forward or backward "        
            if isForward:
                newPart[1] += self.distMove
            else:
                newPart[1] += -self.distMove

            newPart[1] = random.gauss(newPart[1],self.distVar)

            if newPart[1] <= newSpline.dist_u(1.0) and newPart[1] >= newSpline.dist_u(0.0):  
                newParticles.append(newPart)
        """
        
        self.updateCount += 1
        #self.particleSnapshots[self.updateCount] = newParticles
        self.particleSnapshots[self.updateCount] = totalParticles
        self.paths[self.updateCount]  = deepcopy(newPath)    


        self.drawState(nodeID)


    def drawState(self, nodeID = 0):
    
        pylab.clf()
    
        path = self.paths[self.updateCount]
        currSpline = SplineFit(path, smooth=0.1)

        splinePoints = currSpline.getUniformSamples()

        xP = []
        yP = []
        for p in splinePoints:
            xP.append(p[0])
            yP.append(p[1])
        
        pylab.plot(xP,yP, color='r')
        
        particles = self.particleSnapshots[self.updateCount]
        
        
        xP = []
        yP = []
        clusters = {}
        for k in range(len(particles)):
            p = particles[k]
            distVal = p[1]
            
            splinePoint = currSpline.getPointOfDist(distVal)
            
            if k < len(self.resultVector):
                try:
                    clusters[self.resultVector[k]]
                except:
                    clusters[self.resultVector[k]] = []
                
                clusters[self.resultVector[k]].append(splinePoint)
            else:
                xP.append(splinePoint[0])
                yP.append(splinePoint[1])
         
        if len(xP) > 0:       
            pylab.scatter(xP,yP, color='k')
        
        for k, v in clusters.iteritems():
            
            xP = []
            yP = []
            for p in v:
                xP.append(p[0])
                yP.append(p[1])
                
            pylab.scatter(xP,yP, color=self.colors[k], linewidth=1, zorder=10)    
 


        if False and self.mostLikely1 != 0:
            likelyPoint = currSpline.getPointOfDist(self.mostLikely1[1])
            xP = [likelyPoint[0]]
            yP = [likelyPoint[1]]
            pylab.scatter(xP,yP, color='k')

        if False and self.mostLikely2 != 0:
            likelyPoint = currSpline.getPointOfDist(self.mostLikely2[1])
            xP = [likelyPoint[0]]
            yP = [likelyPoint[1]]
            pylab.scatter(xP,yP, color='k')

        #if self.mostLikely != 0:
        #    likelyPoint = currSpline.getPointOfDist(self.mostLikely[1])
        #    xP = [likelyPoint[0]]
        #    yP = [likelyPoint[1]]
        #    pylab.scatter(xP,yP, color='k')

        
        originPoint = currSpline.getPointOfDist(0.0)
        pylab.scatter([originPoint[0]], [originPoint[1]], color='g')
        
        
        node = self.nodeHash[nodeID]
        
        gndPose1 = node.getGlobalGPACPose()
        medial1 = node.getBestMedialAxis()
        
        " set the origin of pose 1 "
        poseGnd = Pose(gndPose1)

        xP = []
        yP = []
        for p in medial1:
            p1 = poseGnd.convertLocalToGlobal(p)
            xP.append(p1[0])
            yP.append(p1[1])
        
        pylab.plot(xP,yP, color='k')
            

        for k in range(nodeID+1):
            print self.nodeHash[k].getGndGlobalGPACPose()
        
        
        pylab.xlim(-5,10)
        pylab.ylim(-8,8)
        
        pylab.title("%d %d %d %3.2f" % (self.lossCount, self.moveForward, self.moveBackward, self.avgDist))
        
        
        pylab.savefig("particles_%04u.png" % self.updateCount)
    
