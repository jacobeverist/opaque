import multiprocessing as processing
import ctypes, os, sys
from SplineFit import SplineFit
import gen_icp
from Pose import Pose
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath
from functions import *
from StableCurve import StableCurve
from Paths import computePathAngleVariance
import math
from operator import itemgetter
import time
import pylab

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


def __remote_multiGenerate(rank, qin, qout):

	sys.stdout = open("gen_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("gen_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

	print "started __remote_multiGenerate"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		print rank, "received", nc, args
		# process data
		#knn = __do_nothing(data, nc, someArg2, someArg3)
		results = []
		for arg in args:
			mapHyp = arg[0]
			generateMap(mapHyp)
			results.append(mapHyp)
						   
		# write to output queue
		qout.put((nc,results))


def __remote_prof_multiGenerate(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_multiGenerate(rank, qin, qout)', "gen_%d.prof" % pid )
		cProfile.runctx("__remote_multiGenerate(rank, qin, qout)", globals(), locals(), "gen_%d.prof" % pid)	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


@logFunction
def batchGenerate(mapHyps):


	global renderGlobalPlotCount
	global pool_generate
	global qin_generate
	global qout_generate



	ndata = len(mapHyps)
	
	args = []
	keys = mapHyps.keys()
	
	for k in range(ndata):
		arg = [mapHyps[keys[k]], k + renderGlobalPlotCount]
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
	

	if len(pool_generate) == 0:
		# set up a pool of processes and data queues

		" make infinite size queues so we don't have to worry about specifying max size "
		qin_generate = processing.Queue(maxsize=0)   
		qout_generate = processing.Queue(maxsize=0)
		pool_generate = [processing.Process(target=__remote_multiGenerate,
		#pool_generate = [processing.Process(target=__remote_prof_multiGenerate,
					args=(rank, qin_generate, qout_generate))
						for rank in range(nproc)]
		for p in pool_generate: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	print "args =", args
	while 1:
		#_data = data[:,cur:cur+chunk_size]
		_data = args[cur:cur+chunk_size]
		#print "_data =", _data
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		#print "put(", (nc,_data), ")"
		qin_generate.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout_generate.get()]
	
	print "received output:", knn
		
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "sorted output:", knn
		

	#print "terminating pool"
	# terminate workers
	#for p in pool:
	#	print "terminate"
	#	p.terminate()
	#	print "terminated"

	hypSet = {}
	for hyp in knn:
		hypSet[hyp.hypothesisID] = hyp
		
	print "returning"
	return hypSet

def __remote_prof_multiEval(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_multiEval(rank, qin, qout)', "eval_%d.prof" % pid )
		cProfile.runctx("__remote_multiEval(rank, qin, qout)", globals(), locals(), "eval_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]



def __remote_multiEval(rank, qin, qout):

	sys.stdout = open("eval_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("eval_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

	print "started __remote_multiEval"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		print rank, "received", nc, args
		# process data
		results = []
		for arg in args:
			mapHyp = arg[0]
			nodeID = arg[1]
			utilVal = computeEvalNode(mapHyp, nodeID)
			results.append(utilVal)
						   
		# write to output queue
		qout.put((nc,results))



@logFunction
def batchEval(mapHyps):


	global renderGlobalPlotCount
	global pool_eval
	global qin_eval
	global qout_eval

	for hid, mapHyp in mapHyps.iteritems():

		poseData = mapHyps[0].poseData

		ndata = len(mapHyp.nodePoses)

		args = []
		keys = mapHyp.nodePoses.keys()

		for k in range(ndata):
			nodeID = keys[k]
			arg = [mapHyp, nodeID, k + renderGlobalPlotCount]
			args.append(arg)
		
		renderGlobalPlotCount += len(mapHyp.nodePoses)

		nproc = __num_processors()
		
		nproc *= 2
		
		# compute chunk size
		chunk_size = ndata / nproc
		chunk_size = 2 if chunk_size < 2 else chunk_size
		
		print "chunk_size =", chunk_size
		print "max_size =", ndata/chunk_size
		
		if len(pool_eval) == 0:
			# set up a pool of processes and data queues

			" make infinite size queues so we don't have to worry about specifying max size "
			qin_eval = processing.Queue(maxsize=0)   
			qout_eval = processing.Queue(maxsize=0)
			pool_eval = [processing.Process(target=__remote_multiEval,
			#pool_eval = [processing.Process(target=__remote_prof_multiEval,
						args=(rank, qin_eval, qout_eval))
							for rank in range(nproc)]
			for p in pool_eval: p.start()
		
		# put data chunks in input queue
		cur, nc = 0, 0
		print "args =", args
		while 1:
			#_data = data[:,cur:cur+chunk_size]
			_data = args[cur:cur+chunk_size]
			#print "_data =", _data
			print "nc = ", nc
			print "cur =", cur
			if len(_data) == 0: break
			#print "put(", (nc,_data), ")"
			qin_eval.put((nc,_data))
			print "DONE"
			cur += chunk_size
			nc += 1
		
		print "BATCH FINISHED"
		
		
		# read output queue
		knn = []
		while len(knn) < nc:
			knn += [qout_eval.get()]
		
		print "received output:", knn
			
		# avoid race condition
		_knn = [n for i,n in sorted(knn)]
		knn = []
		for tmp in _knn:
			knn += tmp

		print "sorted output:", knn
			

		#print "terminating pool"
		# terminate workers
		#for p in pool:
		#	print "terminate"
		#	p.terminate()
		#	print "terminated"

		utilSum = 0.0
		for utilVal in knn:
			utilSum += utilVal

		mapHyp.utility = utilSum

	print "returning"
	#return hypSet


@logFunction
def computeEvalNode(mapHyp, nodeID1):

	poseData = mapHyp.poseData

	print "mapHyp:", mapHyp
	print "nodeID1:", nodeID1
	print "mapHyp.nodePoses:", mapHyp.nodePoses

	estPose1 = mapHyp.nodePoses[nodeID1]

	utilVal = 0.0

	orderedPathIDs1 = mapHyp.getOrderedOverlappingPaths(nodeID1)
	print nodeID1, "recompute orderedPathIDs1:", orderedPathIDs1

	hull1 = poseData.aHulls[nodeID1]
	medial1 = poseData.medialAxes[nodeID1]
	origPose1 = mapHyp.origPoses[nodeID1]

	splicedPaths1 = mapHyp.splicePathIDs(orderedPathIDs1)

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
		

	return utilVal


def __remote_prof_multiMovePath(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_multiMovePath(rank, qin, qout)', "move_%d.prof" % pid )
		cProfile.runctx("__remote_multiMovePath(rank, qin, qout)", globals(), locals(), "move_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]




def __remote_multiMovePath(rank, qin, qout):

	sys.stdout = open("move_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("move_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

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
			nodeID = arg[1]
			direction = arg[2]
			distEst = arg[3]
			movePath(mapHyp, nodeID, direction, distEst = distEst)
			results.append(mapHyp)
						   
		# write to output queue
		qout.put((nc,results))



@logFunction
def batchMovePath(mapHyps, nodeID, direction, distEst = 1.0):


	global renderGlobalPlotCount
	global pool_move 
	global qin_move 
	global qout_move 


	#initGuess

	#arg = [initGuess, initPose, pathID, nodeID]

	ndata = len(mapHyps)
	
	args = []
	keys = mapHyps.keys()
	
	for k in range(ndata):
		arg = [mapHyps[keys[k]], nodeID, direction, distEst, k + renderGlobalPlotCount]
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
	
	#__remote_multiMovePath(rank, qin, qout, nodeID, direction, distEst):

	if len(pool_move) == 0:
		# set up a pool of processes and data queues
		#qin = processing.Queue(maxsize=ndata/chunk_size)
		#qout = processing.Queue(maxsize=ndata/chunk_size)

		" make infinite size queues so we don't have to worry about specifying max size "
		qin_move = processing.Queue(maxsize=0)   
		qout_move = processing.Queue(maxsize=0)
		pool_move = [processing.Process(target=__remote_multiMovePath,
		#pool_move = [processing.Process(target=__remote_prof_multiMovePath,
					args=(rank, qin_move, qout_move))
						for rank in range(nproc)]
		for p in pool_move: p.start()
	
	# put data chunks in input queue
	cur, nc = 0, 0
	print "args =", args
	while 1:
		#_data = data[:,cur:cur+chunk_size]
		_data = args[cur:cur+chunk_size]
		#print "_data =", _data
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		#print "put(", (nc,_data), ")"
		qin_move.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout_move.get()]
	
	print "received output:", knn
		
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "sorted output:", knn
		

	#print "terminating pool"
	# terminate workers
	#for p in pool:
	#	print "terminate"
	#	p.terminate()
	#	print "terminated"

	hypSet = {}
	for hyp in knn:
		hypSet[hyp.hypothesisID] = hyp
		
	print "returning"
	return hypSet

def __remote_prof_multiLocalize(rank, qin, qout):

	try:
		pid = os.getpid()
		" while loop "		
		#cProfile.run('__remote_multiLocalize(rank, qin, qout)', "localize_%d.prof" % pid )
		cProfile.runctx("__remote_multiLocalize(rank, qin, qout)", globals(), locals(), "localize_%d.prof" % pid)	
	
	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]


def __remote_multiLocalize(rank, qin, qout):

	sys.stdout = open("localize_" + str(os.getpid()) + ".out", "w")
	sys.stderr = open("localize_" + str(os.getpid()) + ".err", "w")
	print 'module name:', __name__
	print 'parent process:', os.getppid()
	print 'process id:', os.getpid()

	print "started __remote_multiLocalize"

	while 1:
		# read input queue (block until data arrives)
		nc, args = qin.get()
		print rank, "received", nc, args
		# process data
		#knn = __do_nothing(data, nc, someArg2, someArg3)
		results = []
		for arg in args:
			mapHyp = arg[0]
			nodeID1 = arg[1]
			nodeID2 = arg[2]
			localizePair(mapHyp, nodeID1, nodeID2)
			results.append(mapHyp)
						   
		# write to output queue
		qout.put((nc,results))


@logFunction
def batchLocalizePair(mapHyps, nodeID1, nodeID2):


	global renderGlobalPlotCount
	global pool_localize
	global qin_localize
	global qout_localize


	#initGuess

	#arg = [initGuess, initPose, pathID, nodeID]

	ndata = len(mapHyps)
	
	args = []
	keys = mapHyps.keys()
	
	for k in range(ndata):
		arg = [mapHyps[keys[k]], nodeID1, nodeID2, k + renderGlobalPlotCount]
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
	
	if len(pool_localize) == 0:
		# set up a pool of processes
		#qin = processing.Queue(maxsize=ndata/chunk_size)
		#qout = processing.Queue(maxsize=ndata/chunk_size)

		" make infinite size queues so we don't have to worry about specifying max size "
		qin_localize = processing.Queue(maxsize=0)   
		qout_localize = processing.Queue(maxsize=0)

		pool_localize = [processing.Process(target=__remote_multiLocalize,
		#pool_localize = [processing.Process(target=__remote_prof_multiLocalize,
					args=(rank, qin_localize, qout_localize))
						for rank in range(nproc)]

		for p in pool_localize: p.start()
		
	# put data chunks in input queue
	cur, nc = 0, 0
	print "args =", args
	while 1:
		#_data = data[:,cur:cur+chunk_size]
		_data = args[cur:cur+chunk_size]
		#print "_data =", _data
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		#print "put(", (nc,_data), ")"
		qin_localize.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout_localize.get()]
	
	print "received output:", knn
		
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "sorted output:", knn
		

	#print "terminating pool"
	# terminate workers
	#for p in pool:
	#	print "terminate"
	#	p.terminate()
	#	print "terminated"

	hypSet = {}
	for hyp in knn:
		hypSet[hyp.hypothesisID] = hyp
		
	print "returning"
	return hypSet

	

@logFunction
def movePath(mapHyp, nodeID, direction, distEst = 1.0):
	
	print "movePath(", nodeID, ",", direction, ",", distEst, ")"

	poseData = mapHyp.poseData
	
	if nodeID > 0:
		
		if nodeID % 2 == 1:
			" the guess process gives a meaningful guess for these node's poses "
			" before this, it is meaningless "
			#getStepGuess(mapHyp, nodeID-1, direction)
			estPose0 = mapHyp.nodePoses[nodeID-3]		
			estPose2 = mapHyp.nodePoses[nodeID-1]
			mapHyp.nodePoses[nodeID-1] = getStepGuess(poseData, nodeID-3, nodeID-1, estPose0, estPose2, direction)

			estPose3 = mapHyp.nodePoses[nodeID]
			#getInPlaceGuess(mapHyp, nodeID-1, nodeID, direction)
			supportLine = mapHyp.paths[0]
			mapHyp.nodePoses[nodeID] = getInPlaceGuess(poseData, nodeID-1, nodeID, estPose2, estPose3, supportLine, direction)


		else:
			print "movePath:  DO NOTHING"
		
		" now we fit the guessed pose to a splice of the paths "
		" this constrains it to the previously explored environment "
		" all possible splices are used, and they are filtered to determine the best result "

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
					pID = params[2][0]
					junctions2.append((pathID,pID))

				if dist3 < 3.0:
					pID = params[2][0]
					junctions3.append((pathID,pID))
						
				
			closePathID2 = mapHyp.getClosestPath(initPose2)

			print "closePathID2:", closePathID2
			
			pathSet2 = [closePathID2]

			junctionSet = junctions2 + junctions3
			junctionSet = list(set(junctionSet))

			print "junctions2:", junctions2
			print "junctions3:", junctions3

			" add the parent and the child path ID connected to this junction "
			for junc in junctionSet:
				pathSet2.append(junc[0])
				pathSet2.append(junc[1])

			print "pathSet2:", pathSet2

			pathSet2 = list(set(pathSet2))
			
			pathSet2.sort()

			print "pathSet2:", pathSet2
			

			" find the terminals for each path ID that we have selected "
			termSet2 = []
			for pathID in pathSet2:
				if pathID == 0:
					termSet2.append(terminals[0][0])
					termSet2.append(terminals[1][0])
				else:
					termSet2.append(terminals[pathID+1][0])
					
			" add the parent path ID for each path ID that we have selected "
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

				
				resultPose2, lastCost2, matchCount2, currAng2, currU2 = gen_icp.globalPathToNodeOverlapICP2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial2, plotIter = False, n1 = nodeID-1, n2 = -1, arcLimit = 0.1)
				resultPose3, lastCost3, matchCount3, currAng3, currU3 = gen_icp.globalPathToNodeOverlapICP2([u3, uMedialOrigin3, 0.0], orientedSplicePath, medial3, plotIter = False, n1 = nodeID, n2 = -1, arcLimit = 0.1)
				
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



			" select the best pose for each node "
			" record the splice used for the fit "
			currSplice2 = []
			currSplice3 = []
			if len(resultMoves2) > 0:
				mapHyp.nodePoses[nodeID-1] = resultMoves2[0][0]
				currSplice2 = resultMoves2[0][19]
			else:
				print "node", nodeID-1, "not movePathed because no valid pose"

			if len(resultMoves3) > 0:
				mapHyp.nodePoses[nodeID] = resultMoves3[0][0]
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

			oldPose0 = mapHyp.nodePoses[nodeID-3]
			minDist0, oldU0, oldP0 = orientedPathSpline2.findClosestPoint(oldPose0)
			#oldP0, oldI0, minDist0 = gen_icp.findClosestPointInA(currSplice2, oldPose0)

			oldPose1 = mapHyp.nodePoses[nodeID-2]
			minDist1, oldU1, oldP1 = orientedPathSpline3.findClosestPoint(oldPose1)
			#oldP1, oldI1, minDist1 = gen_icp.findClosestPointInA(currSplice3, oldPose1)


			newPose2 = mapHyp.nodePoses[nodeID-1]
			minDist2, newU2, newP2 = orientedPathSpline2.findClosestPoint(newPose2)
			#newP2, newI2, minDist2 = gen_icp.findClosestPointInA(currSplice2, newPose2)

			newPose3 = mapHyp.nodePoses[nodeID]
			minDist3, newU3, newP3 = orientedPathSpline3.findClosestPoint(newPose3)
			#newP3, newI3, minDist3 = gen_icp.findClosestPointInA(currSplice3, newPose3)

			travelDist2 = orientedPathSpline2.dist(oldU0, newU2)
			if oldU0 > newU2:
				travelDist2 = -travelDist2
			travelDist3 = orientedPathSpline3.dist(oldU1, newU3)
			if oldU1 > newU3:
				travelDist3 = -travelDist3
			
			orientedPoints2 = orientedPathSpline2.getUniformSamples()
			orientedPoints3 = orientedPathSpline3.getUniformSamples()

			#mapHyp.displacePoseParticles(nodeID-1, nodeID, travelDist2, travelDist3)

			if False:

				pylab.clf()

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

				pylab.scatter([oldPose0[0],oldPose1[0]], [oldPose0[1],oldPose1[1]], color='r')
				pylab.scatter([newPose2[0],newPose3[0]], [newPose2[1],newPose3[1]], color='b')

				medial2 = poseData.medialAxes[nodeID-1]
				poseOrigin2 = Pose(newPose2)
				xP = []
				yP = []
				for p in medial2:
					p1 = poseOrigin2.convertLocalToGlobal(p)
					xP.append(p1[0])
					yP.append(p1[1])

				pylab.plot(xP,yP, color = 'b', alpha = 0.2, zorder=9)	

				medial3 = poseData.medialAxes[nodeID]
				poseOrigin3 = Pose(newPose3)
				xP = []
				yP = []
				for p in medial3:
					p1 = poseOrigin3.convertLocalToGlobal(p)
					xP.append(p1[0])
					yP.append(p1[1])

				pylab.plot(xP,yP, color = 'b', alpha = 0.2, zorder=9)	

				pylab.xlim(-5,10)
				pylab.ylim(-8,8)
				pylab.title("%1.3f %1.3f" % (travelDist2, travelDist3))
				pylab.savefig("moveEstimate_%04u_%04u.png" % (nodeID, mapHyp.hypothesisID))
				

				#mapHyp.drawPoseParticles()
				

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

	result, hist = gen_icp.overlapICP_GPU2(estPose1, [u1, u2, angGuess], medial1, medial2, [0.0,0.0], [0.0,0.0], inPlace = inPlace, plotIter = False, n1 = nodeID1, n2 = nodeID2, uRange = uRange)

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
		


@logFunction
def generateAll(hypSet):
	for pID, mapHyp in hypSet.iteritems():
		generateMap(mapHyp)

@logFunction
def generateMap(mapHyp):

	mapHyp.generatePaths()
	" TRIM THE RAW PATHS TO ONLY EXTRUDE FROM THEIR BRANCHING POINT "
	# redundant
	#mapHyp.trimPaths()


@logFunction
def computeLocalDivergence(hypSet, nodeID1, nodeID2):

	poseData = hypSet.values()[0].poseData

	for pID, mapHyp in hypSet.iteritems():

		
		" GET THE ORDERED LIST OF OVERLAPPING PATHS FOR EACH NODE "
		
		" the overlapping paths are computed from the initial guess of position "
		orderedPathIDs1 = mapHyp.getOrderedOverlappingPaths(nodeID1)
		orderedPathIDs2 = mapHyp.getOrderedOverlappingPaths(nodeID2)
		
		mapHyp.orderedPathIDs1 = orderedPathIDs1
		mapHyp.orderedPathIDs2 = orderedPathIDs2

		print nodeID1, "orderedPathIDs1:", orderedPathIDs1
		print nodeID2, "orderedPathIDs2:", orderedPathIDs2

		" departure events for node 1 "
		mapHyp.departures1 = []
		mapHyp.interiors1 = []
		mapHyp.depPoints1 = []
		mapHyp.distances1 = []
		mapHyp.depAngles1 = []
		mapHyp.contig1 = []

		" departure events for node 2 "
		mapHyp.departures2 = []
		mapHyp.interiors2 = []
		mapHyp.depPoints2 = []
		mapHyp.distances2 = []
		mapHyp.depAngles2 = []
		mapHyp.contig2 = []

		
		" COMPUTE DEPARTURE EVENTS FOR EACH OVERLAPPING PATH SECTION "
		for pathID in orderedPathIDs1:
			resultSet = mapHyp.getDeparturePoint(mapHyp.trimmedPaths[pathID], nodeID1, plotIter = False)
			mapHyp.departureResultSet1 = resultSet

			departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = resultSet
			mapHyp.departures1.append([isExist1,isExist2])
			mapHyp.interiors1.append([isInterior1, isInterior2])
			mapHyp.depPoints1.append([departurePoint1, departurePoint2])
			mapHyp.distances1.append([discDist1, discDist2])
			mapHyp.depAngles1.append([depAngle1, depAngle2])
			mapHyp.contig1.append((contigFrac, overlapSum))

		for pathID in orderedPathIDs2:
			resultSet = mapHyp.getDeparturePoint(mapHyp.trimmedPaths[pathID], nodeID2, plotIter = False)
			mapHyp.departureResultSet2 = resultSet

			departurePoint1, depAngle1, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = resultSet
			mapHyp.departures2.append([isExist1,isExist2])
			mapHyp.interiors2.append([isInterior1, isInterior2])
			mapHyp.depPoints2.append([departurePoint1, departurePoint2])
			mapHyp.distances2.append([discDist1, discDist2])
			mapHyp.depAngles2.append([depAngle1, depAngle2])
			mapHyp.contig2.append((contigFrac, overlapSum))


		print "node departures", nodeID1, ":", mapHyp.departures1
		print "node  interiors", nodeID1, ":", mapHyp.interiors1
		print "node departures", nodeID2, ":", mapHyp.departures2
		print "node  interiors", nodeID2, ":", mapHyp.interiors2
		print "node contiguity", nodeID1, ":", mapHyp.contig1
		print "node contiguity", nodeID2, ":", mapHyp.contig2



@logFunction
def checkForeBranch(hypSet, nodeID1, nodeID2, particleIDs):
	
	" 60 degree threshold "
	ANG_THRESH = 1.047

	newHyps = {}
	poseData = hypSet.values()[0].poseData


	for pID, mapHyp in hypSet.iteritems():
		" new junction finding logic "
		" if terminal departures for each medial axis are None or exterior, than we stay on existing paths "
		" if a terminal departure exists that is internal, than we have a new junction "

		" NODE1: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
		frontExist1 = mapHyp.departures1[0][0]
		frontInterior1 = mapHyp.interiors1[0][0]
		foreTerm1 = frontExist1 and frontInterior1

		" DISCREPANCY BETWEEN TIP-CLOSEST and DEPARTURE-CLOSEST POINT ON PATH "				
		
		" NODE2: CHECK FRONT AND BACK FOR BRANCHING EVENTS "
		frontExist2 = mapHyp.departures2[0][0]
		frontInterior2 = mapHyp.interiors2[0][0]
		foreTerm2 = frontExist2 and frontInterior2

		depAngle1 = mapHyp.depAngles1[0][0]
		depAngle2 = mapHyp.depAngles2[0][0]

		frontAngDiff = diffAngle(depAngle1, depAngle2)
		
		depPoint1 = mapHyp.depPoints1[0][0]
		depPoint2 = mapHyp.depPoints2[0][0]
		
		parentPathID1 = mapHyp.orderedPathIDs1[0]
		parentPathID2 = mapHyp.orderedPathIDs2[0]

		isFront1 = poseData.faceDirs[nodeID1]
		isFront2 = poseData.faceDirs[nodeID2]
		if isFront1 and not isFront2:
			dirFlag = 0
		elif not isFront1 and isFront2:
			dirFlag = 1
		else:
			print isFront1, isFront2
			raise	


		isUnique1, duplicatePathID1 = mapHyp.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
		isUnique2, duplicatePathID2 = mapHyp.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

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
		if isBranched:

			" create a new map state where a branch decision is not made "
			newHyps[particleIDs] = mapHyp.copy(particleIDs)
			newHyps[particleIDs].isNotBranched = True
			newHyps[particleIDs].isNodeBranching[nodeID1] = True
			newHyps[particleIDs].isNodeBranching[nodeID2] = True
			print "creating hyp", particleIDs, "from hyp", mapHyp.hypothesisID, ", len(paths) =", len(mapHyp.pathClasses)
			particleIDs += 1


		isBranch, pathBranchIDs, isNew = mapHyp.determineBranchPair(nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag, isUnique1, isUnique2, duplicatePathID1, duplicatePathID2)


		print "determineBranchPair:"
		print isBranch
		print pathBranchIDs

		if isBranch[0]:
			mapHyp.orderedPathIDs1.insert(0,pathBranchIDs[0])
			mapHyp.departures1.insert(0, [False, False])
			mapHyp.interiors1.insert(0, [False, False])
			mapHyp.depPoints1.insert(0, [None, None])

		if isBranch[1]:
			mapHyp.orderedPathIDs2.insert(0,pathBranchIDs[1])
			mapHyp.departures2.insert(0, [False, False])
			mapHyp.interiors2.insert(0, [False, False])
			mapHyp.depPoints2.insert(0, [None, None])


	hypSet.update(newHyps)

	return hypSet, particleIDs

@logFunction
def checkBackBranch(hypSet, nodeID1, nodeID2, particleIDs):

	" 60 degree threshold "
	ANG_THRESH = 1.047

	newHyps = {}
	poseData = hypSet.values()[0].poseData

	for pID, mapHyp in hypSet.iteritems():

		backExist1 = mapHyp.departures1[-1][1]
		backInterior1 = mapHyp.interiors1[-1][1]
		backTerm1 = backExist1 and backInterior1

		backExist2 = mapHyp.departures2[-1][1]
		backInterior2 = mapHyp.interiors2[-1][1]
		backTerm2 = backExist2 and backInterior2
		
		depAngle1 = mapHyp.depAngles1[-1][1]
		depAngle2 = mapHyp.depAngles2[-1][1]
		
		backAngDiff = diffAngle(depAngle1, depAngle2)

		depPoint1 = mapHyp.depPoints1[-1][1]
		depPoint2 = mapHyp.depPoints2[-1][1]

		parentPathID1 = mapHyp.orderedPathIDs1[-1]
		parentPathID2 = mapHyp.orderedPathIDs2[-1]

		isFront1 = poseData.faceDirs[nodeID1]
		isFront2 = poseData.faceDirs[nodeID2]
		if isFront1 and not isFront2:
			dirFlag = 1
		elif not isFront1 and isFront2:
			dirFlag = 0
		else:
			print isFront1, isFront2
			raise	

		isUnique1, duplicatePathID1 = mapHyp.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
		isUnique2, duplicatePathID2 = mapHyp.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

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
		if isBranched:

			" create a new map state where a branch decision is not made "
			newHyps[particleIDs] = mapHyp.copy(particleIDs)
			newHyps[particleIDs].isNotBranched = True
			newHyps[particleIDs].isNodeBranching[nodeID1] = True
			newHyps[particleIDs].isNodeBranching[nodeID2] = True
			print "creating hyp", particleIDs, "from hyp", mapHyp.hypothesisID, ", len(paths) =", len(mapHyp.pathClasses)
			particleIDs += 1

		isBranch, pathBranchIDs, isNew = mapHyp.determineBranchPair(nodeID1, nodeID2, backExist1, backExist2, backInterior1, backInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag, isUnique1, isUnique2, duplicatePathID1, duplicatePathID2)
		print "determineBranchPair:"
		print isBranch
		print pathBranchIDs

		#isNew1 = isNew1 or isNew[0]
		#isNew2 = isNew2 or isNew[1]

		if isBranch[0]:
			mapHyp.orderedPathIDs1.append(pathBranchIDs[0])
			mapHyp.departures1.append([False, False])
			mapHyp.interiors1.append([False, False])
			mapHyp.depPoints1.append([None, None])

		if isBranch[1]:
			mapHyp.orderedPathIDs2.append(pathBranchIDs[1])
			mapHyp.departures2.append([False, False])
			mapHyp.interiors2.append([False, False])
			mapHyp.depPoints2.append([None, None])


	hypSet.update(newHyps)
	
	return hypSet, particleIDs

@logFunction
def addToPaths(particleIDs, hypSet, nodeID1, nodeID2):

	poseData = hypSet.values()[0].poseData

	newHyps = {}

	#generateAll(hypSet)

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
			mapHyp.trimPaths()

			mapHyp.initializePoseParticles()


		return particleIDs, hypSet

	hypSet = batchGenerate(hypSet)

	""" check branching of local spline from member path """
	computeLocalDivergence(hypSet, nodeID1, nodeID2)

	for pID, mapHyp in hypSet.iteritems():
		mapHyp.isNotBranched = False


	""" check for branching conditions from path, spawn new hypotheses """
	hypSet, particleIDs = checkForeBranch(hypSet, nodeID1, nodeID2, particleIDs)

	hypSet, particleIDs = checkBackBranch(hypSet, nodeID1, nodeID2, particleIDs)


	for pID, currHyp in hypSet.iteritems():

		orderedPathIDs1 = currHyp.orderedPathIDs1
		orderedPathIDs2 = currHyp.orderedPathIDs2

		" determine which paths are leaves "
		pathIDs = currHyp.getPathIDs()
		isAParent = {}
		for k in pathIDs:
			isAParent[k] = False
		for k in orderedPathIDs1:
			print "index:", k
			currPath = currHyp.getPath(k)
			currParent = currPath["parentID"]
			if currParent != None:
				isAParent[currParent] = True

		
		"add nodes to paths that are the leaves "
		for pathID in orderedPathIDs1:
			if not isAParent[pathID]:				
				currHyp.addNode(nodeID1,pathID)

		pathIDs = currHyp.getPathIDs()
		isAParent = {}
		for k in pathIDs:
			isAParent[k] = False
		for k in orderedPathIDs2:
			print "index:", k
			currPath = currHyp.getPath(k)
			currParent = currPath["parentID"]
			if currParent != None:
				isAParent[currParent] = True

		for pathID in orderedPathIDs2:
			if not isAParent[pathID]:				
				currHyp.addNode(nodeID2,pathID)

	#generateAll(hypSet)
	genSet = {}
	nonSet = {}
	for pID, currHyp in hypSet.iteritems():
		if not currHyp.isNotBranched:
			genSet[pID] = currHyp
		else:
			nonSet[pID] = currHyp


	genSet = batchGenerate(genSet)

	#hypSet = genSet + nonSet
	nonSet.update(genSet)
	hypSet = nonSet
	#hypSet.update(genSet)

	#hypSet = batchGenerate(hypSet)

	#for pID, currHyp in hypSet.iteritems():
	#	currHyp.drawPoseParticles()
					
	return particleIDs, hypSet




@logFunction
def localizePair(mapHyp, nodeID1, nodeID2):

	poseData = mapHyp.poseData

	""" LOCALIZE NODE PAIR """
	if len(mapHyp.paths[0]) > 0:
		
		"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
		"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
		"3)  select choice with the lowest cost "
		"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
		
		" nodeID1:	is it in a junction or a single path? "
		
		
		try:

			mapHyp.origPoses[nodeID1] = mapHyp.nodePoses[nodeID1]
			result = consistentFit(mapHyp, nodeID1, mapHyp.nodePoses[nodeID1], numGuesses = 11)

		except IndexError:
			print "failed to consistentFit node", nodeID1
			pass

			
		#self.drawConstraints(mapHyp, self.statePlotCount)
		#self.statePlotCount += 1
		#self.drawPathAndHull(mapHyp)

		try:
			mapHyp.origPoses[nodeID2] = mapHyp.nodePoses[nodeID2]
			consistentFit(mapHyp, nodeID2, mapHyp.nodePoses[nodeID2], numGuesses = 11)

		except IndexError:
			print "failed to consistentFit node", nodeID2
			pass


		#mapHyp.computeEval()

		mapHyp.generatePaths()

		#self.drawConstraints(mapHyp, self.statePlotCount)
		#self.statePlotCount += 1
		#self.drawPathAndHull(mapHyp)


@logFunction
def consistentFit(mapHyp, nodeID, estPose, numGuesses = 11, excludePathIDs = []):

	poseData = mapHyp.poseData

	splicedPaths1, spliceTerms, splicePathIDs = mapHyp.getSplicesByNearJunction(estPose)


	print "consistentFit(", nodeID

	initGuesses = []
	orientedPaths = []
	
	hull1 = poseData.aHulls[nodeID]
	medial1 = poseData.medialAxes[nodeID]
	medialSpline1 = SplineFit(medial1, smooth=0.1)

	estPose1 = mapHyp.nodePoses[nodeID]		
	poseOrigin = Pose(estPose1)
	
	globalMedial = []
	for p in medial1:
		globalMedial.append(poseOrigin.convertLocalToGlobal(p))

	globalMedialSpline1 = SplineFit(globalMedial, smooth=0.1)
				
	originU1 = medialSpline1.findU([0.0,0.0])	

	resultsBySplice = []

	for spliceIndex in range(len(splicedPaths1)):
		
		path = splicedPaths1[spliceIndex]


		pathSpline = SplineFit(path, smooth=0.1)
		pathU1 = pathSpline.findU(estPose1[:2])
		pathForeU = pathSpline.getUOfDist(originU1, 1.0, distIter = 0.001)
		pathBackU = pathSpline.getUOfDist(originU1, -1.0, distIter = 0.001)
		
		orientedPath = orientPath(path, globalMedial)
		orientedPathSpline = SplineFit(orientedPath, smooth=0.1)

		#medialSpline1 = SplineFit(medial1, smooth=0.1)


		globalSamples = orientedPathSpline.getUniformSamples(spacing = 0.04)
		medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
		
		globalMedialSamples = []
		for p in medialSamples:
			result = poseOrigin.convertLocalOffsetToGlobal(p)	
			globalMedialSamples.append(result)
		
		
		globalVar = computePathAngleVariance(globalSamples)
		medialVar = computePathAngleVariance(globalMedialSamples)

		" now lets find closest points and save their local variances "			
		closestPairs = []
		TERM_DIST = 20
		
		for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
			pG = globalSamples[i]
			minDist = 1e100
			minJ = -1
			for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
				pM = globalMedialSamples[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minJ = j
					
			if minDist < 0.1:
				closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

		for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
			pM = globalMedialSamples[j]
			minDist = 1e100
			minI = -1
			for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
				pG = globalSamples[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minI = i
					
			if minDist < 0.1:
				closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))

		" remove duplicates "
		closestPairs = list(set(closestPairs))
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(5,6))
		print len(closestPairs), "closest pairs"


		if len(closestPairs) > 0:
			originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
			pathU1 = orientedPathSpline.findU(globalSamples[closestPairs[0][0]])

		else:
		
			originU2 = medialSpline1.findU([0.0,0.0])
			pathU1 = orientedPathSpline.findU(estPose1[:2])
		

		time1 = time.time()

		" add guesses in the neighborhood of the closest pathU "
		" 5 forward, 5 backward "
		

		" while loop "
		halfNum = numGuesses/2
		dists = [(k-halfNum)*0.2 for k in range(numGuesses)]	
		for dist in dists:
			
			
			#pathUGuess = i*0.05
			pathUGuess = orientedPathSpline.getUOfDist(pathU1, dist)
			#initGuesses.append([pathUGuess, originU2, 0.0, spliceIndex])
			initGuesses.append([pathUGuess, originU2, 0.0, spliceIndex, orientedPath, medial1, estPose1, [], nodeID])
			
		orientedPaths.append(orientedPath)


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

		"""			
		if isExist1 or isExist2:
			depVal = 50.0
		else:
			depVal = 100.0
		
		
		if not isInterior1 and not isInterior2 and not isExist1 and not isExist2:
			exclude = False
			for exID in excludePathIDs:
				if exID in splicePathIDs[spliceIndex]:
					exclude = True
			
			if dist < 3.0 and contigFrac > 0.7 and angDiff2 < 0.5 and not exclude:
				filteredResults.append(result)
		"""
		

	print "filteredResults:"
	for result in filteredResults:
		resultPose = result[0]
		isInterior1 = result[5]
		isExist1 = result[6]
		isInterior2 = result[11]
		isExist2 = result[12]
		contigFrac = result[15]
		angDiff2 = result[17]
		spliceIndex = result[19]
		dist = result[20]
		lastCost = result[1]
		matchCount = result[2]
		overlapSum = result[16]
		resultIndex = result[21]
		print resultIndex, spliceIndex, [int(isInterior1), int(isInterior2), int(isExist1), int(isExist2)], contigFrac, dist, angDiff2, lastCost, matchCount, overlapSum, spliceTerms[spliceIndex], splicePathIDs[spliceIndex], resultPose



	utilResults = []
	for k in range(len(filteredResults)):
		result = filteredResults[k]
		angDiff2 = result[17]
		contigFrac = result[15]
		dist = result[20]
		isDeparture = result[18]
		utilVal = (0.5-angDiff2)/0.5 + (contigFrac-0.5)/0.5 + (3.0-dist)/3.0 + 1-isDeparture*0.5
		
		filteredResults[k] = result + (utilVal,)
		

	sortedResults = sorted(filteredResults, key=itemgetter(15), reverse=True)
	contigSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		contigSort.append((result[15],result[21]))
	
	sortedResults = sorted(filteredResults, key=itemgetter(20))
	distSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		distSort.append((result[20],result[21]))
	
	sortedResults = sorted(filteredResults, key=itemgetter(17))
	angSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		angSort.append((result[17],result[21]))
	
	sortedResults = sorted(filteredResults, key=itemgetter(18))
	depSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		depSort.append((result[18],result[21]))


	sortedResults = sorted(filteredResults, key=itemgetter(22), reverse = True)
	utilSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		utilSort.append((result[22],result[21]))
	

	
	print "contigSort:", contigSort
	print "distSort:", distSort
	print "angSort:", angSort
	print "depSort:", depSort
	print "utilSort:", utilSort


	#filteredResults = sorted(filteredResults, key=itemgetter(16))
	filteredResults = sorted(filteredResults, key=itemgetter(22), reverse=True)

	spliceIndex = filteredResults[0][19]
	guessPose = filteredResults[0][0]
	mapHyp.nodePoses[nodeID] = guessPose

	if mapHyp.isNodeBranching[nodeID]:
		mapHyp.moveNode(nodeID, splicePathIDs[spliceIndex])


	return filteredResults[0], splicePathIDs[spliceIndex]


@logFunction
def consistentParticleFit(mapHyp, nodeID, estPose, excludePathIDs = []):

	poseData = mapHyp.poseData

	splicedPaths1, spliceTerms, splicePathIDs = mapHyp.getSplicesByNearJunction(estPose)


	print "consistentFit(", nodeID

	initGuesses = []
	orientedPaths = []
	
	hull1 = poseData.aHulls[nodeID]
	medial1 = poseData.medialAxes[nodeID]
	medialSpline1 = SplineFit(medial1, smooth=0.1)

	estPose1 = mapHyp.nodePoses[nodeID]		
	poseOrigin = Pose(estPose1)
	
	globalMedial = []
	for p in medial1:
		globalMedial.append(poseOrigin.convertLocalToGlobal(p))

	globalMedialSpline1 = SplineFit(globalMedial, smooth=0.1)
				
	originU1 = medialSpline1.findU([0.0,0.0])	

	resultsBySplice = []

	for spliceIndex in range(len(splicedPaths1)):
		
		path = splicedPaths1[spliceIndex]


		pathSpline = SplineFit(path, smooth=0.1)
		pathU1 = pathSpline.findU(estPose1[:2])
		pathForeU = pathSpline.getUOfDist(originU1, 1.0, distIter = 0.001)
		pathBackU = pathSpline.getUOfDist(originU1, -1.0, distIter = 0.001)
		
		orientedPath = orientPath(path, globalMedial)
		orientedPathSpline = SplineFit(orientedPath, smooth=0.1)

		#medialSpline1 = SplineFit(medial1, smooth=0.1)


		globalSamples = orientedPathSpline.getUniformSamples(spacing = 0.04)
		medialSamples = medialSpline1.getUniformSamples(spacing = 0.04)
		
		globalMedialSamples = []
		for p in medialSamples:
			result = poseOrigin.convertLocalOffsetToGlobal(p)	
			globalMedialSamples.append(result)
		
		
		globalVar = computePathAngleVariance(globalSamples)
		medialVar = computePathAngleVariance(globalMedialSamples)

		" now lets find closest points and save their local variances "			
		closestPairs = []
		TERM_DIST = 20
		
		for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
			pG = globalSamples[i]
			minDist = 1e100
			minJ = -1
			for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
				pM = globalMedialSamples[j]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minJ = j
					
			if minDist < 0.1:
				closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

		for j in range(TERM_DIST, len(globalMedialSamples)-TERM_DIST):
			pM = globalMedialSamples[j]
			minDist = 1e100
			minI = -1
			for i in range(TERM_DIST, len(globalSamples)-TERM_DIST):
				pG = globalSamples[i]
				dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
				
				if dist < minDist:
					minDist = dist
					minI = i
					
			if minDist < 0.1:
				closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))

		" remove duplicates "
		closestPairs = list(set(closestPairs))
		
		" sort by lowest angular variance"
		closestPairs = sorted(closestPairs, key=itemgetter(5,6))
		print len(closestPairs), "closest pairs"


		if len(closestPairs) > 0:
			originU2 = medialSpline1.findU(medialSamples[closestPairs[0][1]])	
			pathU1 = orientedPathSpline.findU(globalSamples[closestPairs[0][0]])

		else:
		
			originU2 = medialSpline1.findU([0.0,0.0])
			pathU1 = orientedPathSpline.findU(estPose1[:2])
		

		time1 = time.time()

		" add guesses in the neighborhood of the closest pathU "
		" 5 forward, 5 backward "
		

		" while loop "
		halfNum = numGuesses/2
		dists = [(k-halfNum)*0.2 for k in range(numGuesses)]	
		for dist in dists:
			
			
			#pathUGuess = i*0.05
			pathUGuess = orientedPathSpline.getUOfDist(pathU1, dist)
			#initGuesses.append([pathUGuess, originU2, 0.0, spliceIndex])
			initGuesses.append([pathUGuess, originU2, 0.0, spliceIndex, orientedPath, medial1, estPose1, [], nodeID])
			
		orientedPaths.append(orientedPath)


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

		"""			
		if isExist1 or isExist2:
			depVal = 50.0
		else:
			depVal = 100.0
		
		
		if not isInterior1 and not isInterior2 and not isExist1 and not isExist2:
			exclude = False
			for exID in excludePathIDs:
				if exID in splicePathIDs[spliceIndex]:
					exclude = True
			
			if dist < 3.0 and contigFrac > 0.7 and angDiff2 < 0.5 and not exclude:
				filteredResults.append(result)
		"""
		

	print "filteredResults:"
	for result in filteredResults:
		resultPose = result[0]
		isInterior1 = result[5]
		isExist1 = result[6]
		isInterior2 = result[11]
		isExist2 = result[12]
		contigFrac = result[15]
		angDiff2 = result[17]
		spliceIndex = result[19]
		dist = result[20]
		lastCost = result[1]
		matchCount = result[2]
		overlapSum = result[16]
		resultIndex = result[21]
		print resultIndex, spliceIndex, [int(isInterior1), int(isInterior2), int(isExist1), int(isExist2)], contigFrac, dist, angDiff2, lastCost, matchCount, overlapSum, spliceTerms[spliceIndex], splicePathIDs[spliceIndex], resultPose



	utilResults = []
	for k in range(len(filteredResults)):
		result = filteredResults[k]
		angDiff2 = result[17]
		contigFrac = result[15]
		dist = result[20]
		isDeparture = result[18]
		utilVal = (0.5-angDiff2)/0.5 + (contigFrac-0.5)/0.5 + (3.0-dist)/3.0 + 1-isDeparture*0.5
		
		filteredResults[k] = result + (utilVal,)
		

	sortedResults = sorted(filteredResults, key=itemgetter(15), reverse=True)
	contigSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		contigSort.append((result[15],result[21]))
	
	sortedResults = sorted(filteredResults, key=itemgetter(20))
	distSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		distSort.append((result[20],result[21]))
	
	sortedResults = sorted(filteredResults, key=itemgetter(17))
	angSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		angSort.append((result[17],result[21]))
	
	sortedResults = sorted(filteredResults, key=itemgetter(18))
	depSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		depSort.append((result[18],result[21]))


	sortedResults = sorted(filteredResults, key=itemgetter(22), reverse = True)
	utilSort = []
	for k in range(len(sortedResults)):
		result = sortedResults[k]
		utilSort.append((result[22],result[21]))
	

	
	print "contigSort:", contigSort
	print "distSort:", distSort
	print "angSort:", angSort
	print "depSort:", depSort
	print "utilSort:", utilSort


	#filteredResults = sorted(filteredResults, key=itemgetter(16))
	filteredResults = sorted(filteredResults, key=itemgetter(22), reverse=True)

	guessPose = filteredResults[0][0]
	mapHyp.nodePoses[nodeID] = guessPose

	return filteredResults[0]


