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
			partObj = job[1]
			particleIndex = job[2]
			nodeID3 = job[3]
			prevPose0 = job[4]
			prevPose1 = job[5]
			initPose2 = job[6]
			initPose3 = job[7]
			supportLine = job[8]
			pathSplices2 = job[9]
			pathSplices3 = job[10]

			result = displaceParticle( poseData, partObj, pathSplices2, pathSplices3, supportLine, nodeID3, initPose2, initPose3, prevPose0, prevPose1)
			results.append((particleIndex,) + result)
						   
		# write to output queue
		qout.put((nc,results))

def displaceParticle( poseData, partObj, pathSplices2, pathSplices3, supportLine, nodeID3, initPose2, initPose3, prevPose0, prevPose1):

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


		""" if the node has a spatial feature, localize it to the closest landmark """

		""" if the node has a landmark feature, then we try to localize on the junction point """
		""" FIXME: only use the first long path, we don't try and figure out the longest yet """
		currPoses = {nodeID2: currPose2, nodeID3: currPose3}
		poseOffsets = {nodeID2: [0.0,0.0,0.0], nodeID3: [0.0,0.0,0.0]}
		isPoseChanged = {nodeID2: False, nodeID3: False}

		for nodeID in [nodeID2, nodeID3]:

			"""
			self.junctionData[pathID] = {
										"parentID" : parentID,
										"branchNodeID" : branchNodeID,
										"localJunctionPose" : localJunctionPose, 
										"globalJunctionPose" : globalJunctionPose,
										"probDist" : [ 1.0 / float(numBranches) for k in range(numBranches) ],
										"branchPoseDist" : [deepcopy(globalJunctionPose) for k in range(numBranches)],
										"controlPoseDist" : [deepcopy(controlPose) for k in range(numBranches)],
										"arcDists" : [deepcopy(arcDists)],
										"maxLikelihoodBranch" : 0
										}		
			"""



			sF0 = poseData.spatialFeatures[nodeID][0]
			if len(partObj.junctionData) > 0 and (sF0["bloomPoint"] != None or sF0["archPoint"] != None):

				""" FIXME: assume one landmark feature """

				if sF0["bloomPoint"] != None:
					print "DO bloom localize to landmark feature", nodeID
					localJunctionPoint = sF0["bloomPoint"]

				elif sF0["archPoint"] != None:
					print "DO arch localize to landmark feature", nodeID
					localJunctionPoint = sF0["archPoint"]

				""" starting pose of nodeID """
				nodePose0 = currPoses[nodeID]
				poseOrigin0 = Pose(nodePose0)
				globalLandmarkPoint = poseOrigin0.convertLocalToGlobal(localJunctionPoint)

				# self.pathClasses[0] = {"parentID" : None, "branchNodeID" : None, "localJunctionPose" : None, 
								#"sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : None}

				minPathID = None
				minJuncDist = 1e100
				junctionPose0 = None
				#for pathID, values in pathClasses.iteritems():
				for pathID, values in partObj.junctionData.iteritems():

					if pathID != 0 and values["branchNodeID"] != nodeID:

						currJuncPose = values["globalJunctionPose"]

						currDist = sqrt((currJuncPose[0]-globalLandmarkPoint[0])**2 + (currJuncPose[1]-globalLandmarkPoint[1])**2)

						if currDist < minJuncDist:
							minJuncDist = currDist
							minPathID = pathID
							junctionPose0 = currJuncPose

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
					
					poseOffset = [newNodePoint0[0]-nodePose0[0], newNodePoint0[1]-nodePose0[1]]
					print "poseOffset =", poseOffset
					#poseOffset = poseOrigin0.convertGlobalToLocal(newNodePoint0)
					#print "poseOffset =", poseOffset
					
					print "junctionNodePose0 =", junctionNodePose0

			
					print "changing nodeID", nodeID, "from ", nodePose0, "to", junctionNodePose0, "with offset", poseOffset
					poseOffsets[nodeID] = poseOffset
					currPoses[nodeID] = junctionNodePose0
					isPoseChanged[nodeID] = True

					#self.setNodePose(nodeID, junctionNodePose0)

					
				"""
				global pose of node0
				local pose of landmark feature  (localOffset0)
				global pose of existing junction (junction0)
				assuming that existing junction is landmark, convert back to get global pose of node
				"""
	
		if isPoseChanged[nodeID2] and not isPoseChanged[nodeID3]:
			poseOffset2 = poseOffsets[nodeID2]
			tempPose3 = deepcopy(currPose3)
			tempPose3[0] = poseOffset2[0] + currPose3[0]
			tempPose3[1] = poseOffset2[1] + currPose3[1]
			print "changing nodeID3", nodeID3, "from ", currPose3, "to", tempPose3, "with offset", poseOffset2
			currPoses[nodeID3] = tempPose3

		elif isPoseChanged[nodeID3] and not isPoseChanged[nodeID2]:
			poseOffset3 = poseOffsets[nodeID3]
			tempPose2 = deepcopy(currPose2)
			tempPose2[0] = poseOffset3[0] + currPose2[0]
			tempPose2[1] = poseOffset3[1] + currPose2[1]
			print "changing nodeID2", nodeID2, "from ", currPose2, "to", tempPose2, "with offset", poseOffset3
			currPoses[nodeID2] = tempPose2
	
		""" if the node has been localized to a landmark, we need the updated version """
		currPose2 = currPoses[nodeID2]
		currPose3 = currPoses[nodeID3]

		print "final poses:", currPose2, currPose3

		
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
			branchProbVal = job[23]

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

			result = multiParticleFitSplice(params0, params1, orientedPath0, medial0, medial1, initPose0, initPose1, prevMedial0, prevMedial1, prevPose0, prevPose1, pathIDs, nodeID0, nodeID1, particleIndex, hypID = hypID, pathPlotCount = updateCount, branchIndex = branchIndex, spliceIndex = spliceIndex, branchProbVal = branchProbVal)
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


def multiParticleFitSplice(initGuess0, initGuess1, orientedPath, medialAxis0, medialAxis1, initPose0, initPose1, prevMedialAxis0, prevMedialAxis1, prevPose0, prevPose1, pathIDs, nodeID0, nodeID1, particleIndex, hypID = 0, pathPlotCount = 0, branchIndex =0, spliceIndex = 0, branchProbVal = 1.0):


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

	utilVal0 = utilVal0 * branchProbVal

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

		return newParticle

	def __str__(self):
		return str(self.__dict__)

	def __eq__(self, other): 
		resultVal = self.__dict__ == other.__dict__
		#print "compare equals", resultVal, self.__dict__, other.__dict__
		return resultVal
		#return self.__dict__ == other.__dict__

