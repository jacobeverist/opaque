import multiprocessing as processing
import random
from copy import deepcopy
from SplineFit import SplineFit
import pylab
from math import sqrt
from Pose import Pose
import gen_icp
from math import acos, asin
from numpy import array
import ctypes, os, sys

from Splices import getMultiDeparturePoint, orientPath, getCurveOverlap

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

"""
def mysum(lvals):
    s2 = 0
    s = 0
    for e in l:
        s += e
        s2 += e * e
    return (s, s2)

N = len(lvals)

s, s2 = mysum(lvals)

var = (s2 - (s*s) / N ) / N
mean = s / N

"""

def batchParticle(localizeJobs):

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
		#_data = data[:,cur:cur+chunk_size]
		_data = args[cur:cur+chunk_size]
		#print "_data =", _data
		print "nc = ", nc
		print "cur =", cur
		if len(_data) == 0: break
		#print "put(", (nc,_data), ")"
		qin_posePart.put((nc,_data))
		print "DONE"
		cur += chunk_size
		nc += 1
	
	print "BATCH FINISHED"
	
	
	# read output queue
	knn = []
	while len(knn) < nc:
		knn += [qout_posePart.get()]
	
	print "received output:", knn
		
	# avoid race condition
	_knn = [n for i,n in sorted(knn)]
	knn = []
	for tmp in _knn:
		knn += tmp

	print "sorted output:", knn
		

	#print "terminating pool"
	# terminate workers
	#for p in pool_posePart:
	#	print "terminate"
	#	p.terminate()
	#	print "terminated"
		
	print "returning"
	return knn


#multiParticleFitSplice(params0, params1, globalPath, medial0, medial1, initPose0, initPose1, pathIDs, nodeID0, nodeID1, particleIndex, hypID = hypID, pathPlotCount = updateCount, spliceIndex = spliceIndex)
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

		self.pose0 = deepcopy(pose0)
		self.pose1 = deepcopy(pose1)

	def addPath(self, pathID, parentID, branchNodeID, localJunctionPose, globalJunctionPose):

		self.junctionData[pathID] = {
									"parentID" : parentID,
									"branchNodeID" : branchNodeID,
									"localJunctionPose" : localJunctionPose, 
									"globalJunctionPose" : globalJunctionPose,
									#"probDist" : [ 1.0 / 11.0 for k in range(11) ],
									"probDist" : [ 1.0 / 6.0 for k in range(6) ],
									#"branchPoseDist" : [deepcopy(globalJunctionPose) for k in range(11)],
									"branchPoseDist" : [deepcopy(globalJunctionPose) for k in range(6)],
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
    
