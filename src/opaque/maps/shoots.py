import random
import math
import alphamod
import os
import sys
import graph
from functions import *
import Image
import hashlib
from medialaxis import computeMedialAxis
from SplineFit import SplineFit
from Pose import Pose
from LocalNode import getLongestPath
import pylab
import gen_icp
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath, getTipAngles
from MapProcess import selectLocalCommonOrigin, selectCommonOrigin
import traceback
import ctypes
import multiprocessing as processing
from scipy.spatial import cKDTree
from numpy import array



renderGlobalPlotCount = 0
pathPlotCount = 0

qin_branch = None
qout_branch =  None
pool_branch = []

alphaPlotCount = 0
medialCount = 0

def extendBackToHull(path, hull):

	currPath = deepcopy(path)

	backVec = [0.,0.]
	indic = range(3)
	indic.reverse()
	
	""" use the last few edges of the leaf to average the leaf direction """
	for i in indic:
		if i+2 < len(currPath):
			p1 = currPath[-i-3]
			p2 = currPath[-i-1]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			backVec[0] += vec[0]
			backVec[1] += vec[1]

	backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])

	backVec[0] /= backMag
	backVec[1] /= backMag

	""" create an edge with the computed direction extending over the hull boundary """
	newP2 = (currPath[-1][0] + backVec[0]*2, currPath[-1][1] + backVec[1]*2)
	edge2 = [currPath[-1], newP2]
	

	""" find the intersection points with the hull """
	interPoints = []
	for k in range(len(hull)-1):
		hullEdge = [hull[k],hull[k+1]]
		isIntersect2, point2 = Intersect(edge2, hullEdge)
		if isIntersect2:
			interPoints.append(point2)
			break
	
	""" add an extended edge to the termination point at the hull boundary if it intersects """			
	if isIntersect2:
		currPath.append(point2)

	return currPath

def extendToHull(path, hull):

	currPath = deepcopy(path)

	frontVec = [0.,0.]
	backVec = [0.,0.]
	indic = range(3)
	indic.reverse()
	
	for i in indic:
		if i+2 < len(currPath):
			p1 = currPath[i+2]
			p2 = currPath[i]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			frontVec[0] += vec[0]
			frontVec[1] += vec[1]

			p1 = currPath[-i-3]
			p2 = currPath[-i-1]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			backVec[0] += vec[0]
			backVec[1] += vec[1]

	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])

	frontVec[0] /= frontMag
	frontVec[1] /= frontMag
	backVec[0] /= backMag
	backVec[1] /= backMag

	newP1 = (currPath[0][0] + frontVec[0]*2, currPath[0][1] + frontVec[1]*2)
	newP2 = (currPath[-1][0] + backVec[0]*2, currPath[-1][1] + backVec[1]*2)
	edge1 = [newP1, currPath[0]]
	edge2 = [currPath[-1], newP2]

	""" find the intersection points with the hull """
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

	""" replace the extended edges with a termination point at the hull edge """			
	if isIntersect1:
		currPath.insert(0, point1)
	if isIntersect2:
		currPath.append(point2)

	return currPath

def computeSkeletonFromImage(medialPointSoup):

	global alphaPlotCount
	global medialCount

	""" radius constant for alpha shape algorithm """
	ALPHA_RADIUS = 0.2


	""" attempt to compute the alpha hull multiple times until it is successful """
	isDone = False
	while not isDone:

		""" add a little bit of noise to each of the points to avoid degenerate
			conditions in CGAL.  Catch the exception in case of a math domain error """
		perturbPoints = []
		for p in medialPointSoup:

			p2 = copy(p)
			
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
	
			""" the data to be input in text form """
			saveFile = ""	 
			saveFile += "ALPHA_RADIUS = " + repr(ALPHA_RADIUS) + "\n"
			saveFile += "perturbPoints = " + repr(perturbPoints) + "\n"

			""" save the input data for debug purposes in case we crash """
			isWritten = False
			while not isWritten:
				try:
					f = open("doAlphaInput_%08u.txt" % (alphaPlotCount), 'w')
					f.write(saveFile)
					f.close()
				except:
					pass
				else:
					isWritten = True

			""" call CGAL alpha shape 2D function with noise-added points """
			vertices = alphamod.doAlpha(ALPHA_RADIUS,perturbPoints)

			""" delete the temporary file since we were successful """
			os.remove("doAlphaInput_%08u.txt" % (alphaPlotCount))
			alphaPlotCount += 1
			
			
			""" if the hull doesn't have enough vertices, then we throw an exception, retry again """
			numVert = len(vertices)
			if numVert <= 2:
				print "Failed, hull had only", numVert, "vertices"
				raise
			
			""" success!  Now exit the infinite loop """
			isDone = True

		except:
			print "hull has holes!	retrying..."
	
	""" cut out the repeat vertex """
	vertices = vertices[:-1]
	
	""" make the edges of the hull uniform with maximum length """
	vertices = makePolygonUniform(vertices)

	""" add in the repeat vertex again """
	vertices.append(vertices[0])
	

	""" find the bounding box of the points """
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


	""" SPECIFY SIZE OF GRID AND CONVERSION PARAMETERS from the bounding box """
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
		point = (pX,pY)
		return point

	""" CONVERT HULL TO GRID COORDINATES """
	gridHull = []
	for i in range(len(vertices)):
		p = vertices[i]
		gridHull.append(realToGrid(p))

	""" BOUNDING BOX OF GRID HULL """	 
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

	""" COMPUTE MEDIAL AXIS OF HULL, the shoot skeleton """
	resultImg = Image.new('L', (numPixel,numPixel))
	resultImg = computeMedialAxis(medialCount, numPixel,numPixel, 5, resultImg, len(gridHull[:-2]), gridHull[:-2])
	medialCount += 1

	""" EXTRACT POINTS FROM GRID TO LIST """
	imgA = resultImg.load()    
	points = []
	for i in range(1, numPixel-1):
		for j in range(1, numPixel-1):
			if imgA[i,j] == 0:
				points.append((i,j))

	""" CREATE GRAPH NODES FOR EACH POINT """
	medialGraph = graph.graph()
	for p in points:
		medialGraph.add_node(p, [])

	""" ADD EDGES BETWEEN NEIGHBORS """
	for i in range(2, numPixel-2):
		for j in range(2, numPixel-2):
			if imgA[i,j] == 0:
				for k in range(i-1, i+2):
					for l in range(j-1, j+2):
						if imgA[k,l] == 0 and not (k == i and l == j):
							medialGraph.add_edge((i,j), (k,l))
							
	""" COMPUTE MINIMUM SPANNING TREE """
	mst = medialGraph.minimal_spanning_tree()
	
	""" INITIALIZE DATA DICTS FOR UNIDIRECTIONAL MST """
	uni_mst = {}
	for k, v in mst.items():
		uni_mst[k] = []

	
	""" ADD EDGES TO DICT TREE REPRESENTATION """
	for k, v in mst.items():
		if v != None:
			uni_mst[k].append(v)
			uni_mst[v].append(k)

	""" LOCATE ALL LEAVES """
	leaves = []
	for k, v in uni_mst.items():
		if len(v) == 1:
			leaves.append(k)
	
	""" DELETE ALL NODES THAT ARE LEAVES, TO REMOVE SINGLE NODE BRANCHES """
	for leaf in leaves:
		medialGraph.del_node(leaf)	  
		
	""" RECOMPUTE MST """
	mst = medialGraph.minimal_spanning_tree()

	""" AGAIN, CREATE OUR DATA STRUCTURES AND IDENTIFY LEAVES """
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

	""" create dictionary hash of converting grid indexes to real points """
	gridHash = {}
	for k, v in mst.items():
		gridHash[k] = gridToReal(k)
					
	""" RECORD THE LEAVES AND JUNCTIONS """
	leaves = []
	junctions = []
	for k, v in uni_mst.items():
		if len(v) == 1:
			leaves.append(k)

		if len(v) > 2:
			junctions.append(k)

	print "junctions:", junctions
	print "leaves:", leaves

	return vertices, junctions, leaves, uni_mst, gridHash


def computePathSegments(juncIDs, leafIDs, tree, gridHash, vertices):

	leafSegments = []
	internalSegments = []

	retPaths = []

	if len(juncIDs) > 0:

		for i in range(len(juncIDs)): 
			isVisited = {}
			for k, v in tree.items():
				isVisited[k] = 0

			jID = juncIDs[i]
			otherJuncIDs = copy(juncIDs)
			otherJuncIDs.remove(jID)
			isVisited[jID] = 1

			for childNode in tree[jID]:
				resultPath = getSegPaths(childNode, otherJuncIDs, leafIDs, [jID], tree, isVisited)
				if len(resultPath) > 2:
					retPaths.append(resultPath)
	else:
		isVisited = {}
		for k, v in tree.items():
			isVisited[k] = 0

		lID = leafIDs[0]
		otherLeafIDs = copy(leafIDs)
		otherLeafIDs.remove(lID)
		isVisited[lID] = 1

		for childNode in tree[lID]:
			resultPath = getSegPaths(childNode, juncIDs, otherLeafIDs, [lID], tree, isVisited)
			if len(resultPath) > 2:
				retPaths.append(resultPath)
	
	for path in retPaths:

		""" save leaf segments leaf node at end of the list """
		if path[0] in leafIDs:
			path.reverse()
			leafSegments.append(path)

		elif path[-1] in leafIDs:
			leafSegments.append(path)

			""" save internal segment since both terminals are junction nodes """
		elif path[0] in juncIDs and path[-1] in juncIDs:
			internalSegments.append(path)

		else:
			print len(path), path[0], path[-1], juncIDs, leafIDs
			raise


	" remove duplicate internal path segments "
	savedSegments = []
	for j in range(len(internalSegments)):
		thisSegment = internalSegments[j]	

		thisP0 = thisSegment[0]
		thisP1 = thisSegment[-1]

		isUnique = True
		for k in range(len(savedSegments)):
			savedSeg = savedSegments[k]
			
			p0 = savedSeg[0]
			p1 = savedSeg[-1]

			if thisP0 == p0 and thisP1 == p1 or thisP0 == p1 and thisP1 == p0:
				isUnique = False
				
		if isUnique:
			savedSegments.append(thisSegment)

	internalSegments = savedSegments

	
	""" convert grid segments to real segments """
	realLeafSegments = []
	realInternalSegments = []
	for seg in leafSegments:
		newSeg = []
		for p in seg:
			newSeg.append(gridHash[p])

		realLeafSegments.append(newSeg)

	for seg in internalSegments:
		newSeg = []
		for p in seg:
			newSeg.append(gridHash[p])
		realInternalSegments.append(newSeg)



	""" EXTEND LEAF SEGMENTS TO BOUNDARY OF HULL """
	longLeafSegments = []

	if len(realLeafSegments) > 1:

		for seg in realLeafSegments:

			realSeg = extendBackToHull(seg, vertices)

			""" the extended version of the leaf segment """
			longLeafSegments.append(realSeg)

	else:

		seg = realLeafSegments[0]

		realSeg = extendBackToHull(seg, vertices)
		realSeg.reverse()
		realSeg2 = extendBackToHull(realSeg, vertices)
		realSeg2.reverse()

		""" the extended version of the leaf segment """
		longLeafSegments.append(realSeg2)



	""" create smoothed versions of the segments with spline curves """
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
	
	allRealSegments = longLeafSegments + realInternalSegments

	smoothPathSegs = smoothLeafSegments + smoothInternalSegments



	skeletonGraph = graph.graph()

	""" add the leaf and their attached junction points,
		add internal segments, inclusive of junction terminals """
	for seg in allRealSegments:
		for k in range(0,len(seg)):
			skeletonGraph.add_node(tuple(seg[k]))
	
	""" now add edges weighted with cartesian distance """
	for seg in allRealSegments:
		for k in range(0,len(seg)-1):
			p1 = tuple(seg[k])
			p2 = tuple(seg[k+1])

			dist = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

			skeletonGraph.add_edge(p1,p2,wt=dist)
	
	"""
	origControlOrigin = Pose(controlPose)

	localSkeletonGraph = graph.graph()
	
	for edge in skeletonGraph.edges():
	
		nodePoint1 = edge[0]
		nodePoint2 = edge[1]

		localNodePoint1 = tuple(origControlOrigin.convertGlobalToLocal(nodePoint1))
		localNodePoint2 = tuple(origControlOrigin.convertGlobalToLocal(nodePoint2))

		dist = sqrt((localNodePoint1[0]-localNodePoint2[0])**2 + (localNodePoint1[1]-localNodePoint2[1])**2)

		localSkeletonGraph.add_node(localNodePoint1)
		localSkeletonGraph.add_node(localNodePoint2)
		localSkeletonGraph.add_edge(localNodePoint1, localNodePoint2, wt=dist)
	"""

	return smoothLeafSegments, smoothInternalSegments, skeletonGraph
	#return smoothLeafSegments, smoothInternalSegments, localSkeletonGraph

def getSegPaths(node, juncIDs, leafIDs, currPath, tree, isVisited):

	if node in juncIDs:
		isVisited[node] = 1
		nextPath = copy(currPath) + [node]
		return nextPath

	if node in leafIDs:
		isVisited[node] = 1
		nextPath = copy(currPath) + [node]
		return nextPath

	isVisited[node] = 1
	nextPath = copy(currPath) + [node]

	for childNode in tree[node]:
		" if we're not going backwards "
		if not isVisited[childNode]:
			return getSegPaths(childNode, juncIDs, leafIDs, nextPath, tree, isVisited)

	#allPaths = []
	#for childNode in tree[node]:
	#	thesePaths = getSegPaths(childNode, juncIDs, leafIDs, nextPath, tree, isVisited)
	#	allPaths += thesePaths
	
	" should never reach, defective data or algorithm"

	print "nextPath:", nextPath
	print "juncIDs:", juncIDs
	print "leafIDs:", leafIDs
	raise
	return nextPath

@logFunction
def computeGlobalControlPoses(controlPoses, parentPathIDs):

	pathIDs = parentPathIDs.keys()

	isComputed = {}

	finalPoses = {}

	for pathID in pathIDs:
		isComputed[pathID] = False

		if pathID == 0:
			finalPoses[pathID] = controlPoses[pathID]
			isComputed[pathID] = True


	while False in isComputed.values():

		for pathID in pathIDs:
			if pathID != 0:
				parentPathID = parentPathIDs[pathID]

				if isComputed[parentPathID]:
					localFrame = Pose(finalPoses[parentPathID])
					localPose = controlPoses[pathID]

					finalPoses[pathID] = localFrame.convertLocalOffsetToGlobal(localPose)
					isComputed[pathID] = True

	return finalPoses

@logFunction
def computeShootSkeleton(poseData, pathID, globalJunctionPose, nodeSet, nodePoses, hypothesisID, color, topCount, plotIter = False):
	
	print "computeShootSkeleton(", pathID, globalJunctionPose, hypothesisID, topCount, ")"

	""" establish recomputability for debugging """
	random.seed(0)		  

	""" flush the print buffer for debugging """
	sys.stdout.flush()

	""" delete the previous computed shoot skeleton """
	theoryMedialLongPaths = []
	""" Select the leaf-to-leaf path that will be our defacto shoot """
	medialLongPaths = []
	leaf2LeafPathJunctions = {}


	""" get the junction point of this shoot if it is not the root """
	#branchNodeID = self.pathClasses[pathID]["branchNodeID"]
	#globalJunctionPose = None

	#if branchNodeID != None:
	#	globalJunctionPose = self.getGlobalJunctionPose(pathID) 
	

	""" if there is no data, return an empty skeleton """
	#nodeSet = self.getNodes(pathID)
	if nodeSet == []:
		return [], [], {}, [], []


	""" collect the local hulls of each the spatial images attached to this shoot """
	medialPointSoup = []
	for nodeID in nodeSet:

		""" global GPAC pose and alpha hull of the local spatial image """
		#estPose1 = self.getNodePose(nodeID)
		#return copy(self.nodePoses[nodeID])
		estPose1 = copy(nodePoses[nodeID])
		hull1 = poseData.aHulls[nodeID]

		""" print out a hash value for this unique hull for debugging """
		m = hashlib.md5()
		m.update(repr(hull1))
		print nodeID, "hull1 =", int(m.digest().encode('hex'),16)

		""" convert hull points to global coordinates and save """
		poseOrigin = Pose(estPose1)
		for k in range(len(hull1)):

			p = hull1[k]
			
			m = hashlib.md5()
			m.update(repr(p))
			
			p1 = poseOrigin.convertLocalToGlobal(p)
			
			m = hashlib.md5()
			m.update(repr(p1))
			
			medialPointSoup.append(p1)

	""" compute the alpha shape and then the medial axis skeleton from the result """
	vertices, junctions, leaves, uni_mst, gridHash = computeSkeletonFromImage(medialPointSoup)

	allJunctions = []
	if globalJunctionPose != None:

		""" If this is a branching shoot, then we want to find the point on the skeleton
			that is closest to the originally located junction point.
			This is so that we can add an edge that roughly corresponds to point 
			and orientation if one does not already exist """
	
		""" go through all grid points in the unidirectional MST and find the closest point in real coordinates """
		minKey = None
		minCand = None
		minJuncDist = 1e100
		for k, v in uni_mst.items():
			gCand = gridHash[k]
			juncDist = sqrt((globalJunctionPose[0]-gCand[0])**2 + (globalJunctionPose[1]-gCand[1])**2)
			
			if juncDist < minJuncDist:
				minKey = k
				minJuncDist = juncDist
				minCand = gCand
	
		""" compute the direction vector based on globalJunctionPose angle """
		theoryVec = [cos(globalJunctionPose[2]), sin(globalJunctionPose[2])]

		""" junction description:
			0: grid coordinates
			1: real coordinates
			2: closest distance
			3: direction vector """
		theoryJunc = (minKey, minCand, minJuncDist, theoryVec)
		print "theoryJunc:", theoryJunc
		
		""" create a series of edges that extends off of closest point mimicking the theoretical branch """
		theoryJuncPoint = theoryJunc[1]
		theoryLeaf = []
		leafMag = 0.02
		for i in range(5):
			newPoint = (theoryVec[0]*leafMag + theoryJuncPoint[0], theoryVec[1]*leafMag + theoryJuncPoint[1])
			theoryLeaf.append(newPoint)
			leafMag += 0.02


		print "theoryLeaf:", theoryLeaf
		

		"""  Find the junctions and add them to allJunctions:

			If there is a node of degree > 2, use these nodes as the junction points.
			If there is no nodes of degree > 2, find the closest node to the theoretical junction point and use this as the junction index 
		"""
		if len(junctions) > 0:
			
			for junc in junctions:
				gJunc = gridHash[junc]
				juncDist = sqrt((globalJunctionPose[0]-gJunc[0])**2 + (globalJunctionPose[1]-gJunc[1])**2)

				""" 0: grid coordinates
					1: real coordinates
					2: distance from point to real junction
					3: branch direction (None) """
				allJunctions.append((junc, gJunc, juncDist, None))

		else:
			
			print "adding theoretical junction:", theoryJunc
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

	""" get the shoot skeleton as sets of order points representing curve segments
		with junctions and leaves as terminals """
	smoothLeafSegments, smoothInternalSegments, localSkeletonGraph = computePathSegments(junctions, leaves, uni_mst, gridHash, vertices)

	print "computePathSegments:", len(smoothLeafSegments), len(smoothInternalSegments), "paths from", len(junctions), "junctions and", len(leaves), "leaves", [len(pMem) for pMem in smoothLeafSegments], [len(pMem) for pMem in smoothInternalSegments]


	"""
	1) convert to real
	2) get theory leaf as well
	3) extrapolate the leaf segments to the boundary of the hull
	4) spline smooth
	5) convert to oriented points
	6) return set of points based on return condition
	"""
	

	""" FOR EVERY PAIR OF LEAVES, SAVE ITS PATH IF ITS LONG ENOUGH """
	leafToLeafPathTuples = []
	MAX_LEN = 2

	for leaf1 in leaves:
		for leaf2 in leaves:
			if leaf1 < leaf2:
				""" make sure selected route is greater than 2 points """
				if len(nodePaths[leaf1][leaf2]) > MAX_LEN:
					nPath = deepcopy(nodePaths[leaf1][leaf2])

					""" find which junctions are on this path, if any """
					juncIndices = []
					for junc in allJunctions:
						try:
							index = nPath.index(junc[0])
						except:
							juncIndices.append(None)
						else:
							juncIndices.append(index)
							
					""" save the path length, the path, and the junction indexes """
					leafToLeafPathTuples.append((len(nPath),nPath, juncIndices))


	if globalJunctionPose != None:


		""" every path from a leaf to the theoretical junction point """ 
		theoryPaths = []
		print "theoryPaths:"
		for leaf1 in leaves:
			nPath = deepcopy(nodePaths[leaf1][theoryJunc[0]])					 
			theoryPaths.append((len(nPath), nPath))
			print leaf1, theoryJunc, len(nPath)

		""" sort by longest path first """
		theoryPaths.sort(reverse=True)

		""" extend these theoretical leaves to the boundary of the hull """
		for k in range(len(theoryPaths)):
			path = theoryPaths[k][1]
			realPath = []
			for p in path:
				realPath.append(gridHash[p])
			#realPath.reverse()
			
			theoryPaths[k] = realPath + theoryLeaf
			
			print "theoryPath(" , k , ")", len(realPath), realPath[-1], theoryLeaf
			
			theoryPath = theoryPaths[k]
			print "theoryPath:", len(theoryPath)
			
			leafPath = deepcopy(theoryPath)
			
			resultPath = extendToHull(theoryPath, vertices)
			theoryMedialLongPaths.append(resultPath)



	""" SORT FROM THE LONGEST TO SHORTEST """
	leafToLeafPathTuples.sort(reverse=True)
	
	""" Separate the tuple information into individual lists """
	juncGrids = []
	juncReals = []
	juncIndices = []
	leaf2LeafPathLengths = []
	leafToLeafPaths = []
	leafToLeafGridPaths = []
	for k in range(len(leafToLeafPathTuples)):

		juncInds = leafToLeafPathTuples[k][2]

		""" separate list of junc indices for each leaf2leaf path """
		juncIndices.append(juncInds)


		""" only the list of paths itself in real coordinates """
		gridPath = leafToLeafPathTuples[k][1]
		realPath = []
		for p in gridPath:
			realPath.append(gridHash[p])
		leafToLeafPaths.append(realPath)

		leafToLeafGridPaths.append(deepcopy(gridPath))


		""" only the grid and real values of each junction point """
		jGrids = []
		jReals = []
		for index in juncInds:
			if index != None:
				jGrids.append(copy(gridPath[index]))
				jReals.append(copy(leafToLeafPaths[k][index]))
		juncGrids.append(jGrids)
		juncReals.append(jReals)

		""" only the list of path lengths """
		leaf2LeafPathLengths.append(leafToLeafPathTuples[k][0])

		#leafToLeafPaths.append(leafToLeafPathTuples[k][1])
	
	print "juncIndices:", juncIndices

	""" extend these leaf2leaf paths to the boundary of the hull """
	juncAngSet = []
	juncLongIndices = []
	print len(leafToLeafPaths), "long paths"
	for leafIndex in range(len(leafToLeafPaths)):
		
		leaf2LeafPath = leafToLeafPaths[leafIndex]
		print "leaf2LeafPath:", len(leaf2LeafPath)
		print leaf2LeafPath
		
		leafPath = extendToHull(leaf2LeafPath, vertices)
		medialLongPaths.append(leafPath)

		print "globalJunctionPose:", globalJunctionPose
		print "juncIndices:", juncIndices[leafIndex]

		
		""" recompute the junction indices based on the concatenated edges """
		jReals = juncReals[leafIndex]
		mLongPath = medialLongPaths[leafIndex]
		jIndices = []

		""" junction indices in real path """
		for jPoint in jReals:
			try:
				index = mLongPath.index(jPoint)
			except:
				jIndices.append(None)
			else:
				jIndices.append(index)
		
		print "jIndices:", jIndices

		juncLongIndices.append(jIndices)

		""" compute difference of junction index point angle with declared junction angle """
		juncAngs = []
		if globalJunctionPose != None:
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
	
					frontError = normalizeAngle(globalJunctionPose[2]-foreAng)
					backError = normalizeAngle(globalJunctionPose[2]-backAng)
					juncAngs.append((frontError,backError))
				else:
					juncAngs.append(None)
		juncAngSet.append(juncAngs)


	
	""" distance of points to declared junction """
	juncDists = []
	for junc in allJunctions:
		juncDists.append(junc[2])

	""" record all the data into data structure for leaf-to-leaf paths """
	leaf2LeafPathJunctions["juncIndices"] = juncIndices
	leaf2LeafPathJunctions["juncAngSet"] = juncAngSet
	leaf2LeafPathJunctions["juncDists"] = juncDists
	leaf2LeafPathJunctions["leaf2LeafPathLengths"] = leaf2LeafPathLengths
	leaf2LeafPathJunctions["juncReals"] = juncReals
	leaf2LeafPathJunctions["juncGrids"] = juncGrids
	leaf2LeafPathJunctions["juncLongIndices"] = juncLongIndices
	leaf2LeafPathJunctions["leafSegments"] = smoothLeafSegments
	leaf2LeafPathJunctions["internalSegments"] = smoothInternalSegments
	leaf2LeafPathJunctions["skeletonGraph"] = localSkeletonGraph

	localPathSegs = smoothLeafSegments + smoothInternalSegments

	"""
	if globalJunctionPose != None:
		origControlOrigin = Pose(controlPose)
		smoothPathSegs = smoothLeafSegments + smoothInternalSegments

		localPathSegs = []
		for k in range(len(smoothPathSegs)):
			pathSeg = smoothPathSegs[k]
			localSeg = []
			for p in pathSeg:
				p1 = origControlOrigin.convertGlobalPoseToLocal(p)
				localSeg.append(p1)

			localPathSegs.append(localSeg)
	else:
		localPathSegs = smoothLeafSegments + smoothInternalSegments
	"""

	leaf2LeafPathJunctions["localSegments"] = localPathSegs


	""" get the neighbor points along the path at the junction so we can describe the occupancy state """
	juncPoints = []
	juncArmPoints = []
	juncDesc = {}
	#juncGridPoints = []
	#juncGridArmPoints = []
	#juncGridDesc = {}
	for k in range(len(juncLongIndices)):

		print "junc indexes:", juncLongIndices[k], juncIndices[k]

		juncInds = juncLongIndices[k]
		juncPnts = []
		juncArmPnts = []
		for index in juncInds:
			jPnt = copy(medialLongPaths[k][index])

			jaPnt1 = medialLongPaths[k][index+1]
			jaPnt2 = medialLongPaths[k][index-1]

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

		#juncGridInds = juncIndices[k]
		#juncGridPnts = []
		#juncGridArmPnts = []
		#for index in juncGridInds:
		#	jPnt = copy(leafToLeafGridPaths[k][index])
		#
		#	jaPnt1 = leafToLeafGridPaths[k][index+1]
		#	jaPnt2 = leafToLeafGridPaths[k][index-1]
		#
		#	jPnt = tuple(jPnt)
		#	jaPnt1 = tuple(jaPnt1)
		#	jaPnt2 = tuple(jaPnt2)
		#
		#	juncGridPnts.append(jPnt)
		#	juncGridArmPnts.append((jaPnt1, jaPnt2))
		#
		#	try:
		#		juncGridDesc[jPnt].append(jaPnt1)
		#		juncGridDesc[jPnt].append(jaPnt2)
		#	except:
		#		juncGridDesc[jPnt] = []
		#		juncGridDesc[jPnt].append(jaPnt1)
		#		juncGridDesc[jPnt].append(jaPnt2)


		#juncGridPoints.append(juncGridPnts)
		#juncGridArmPoints.append(juncGridArmPnts)


	""" remove duplicates """
	for k, v in juncDesc.iteritems():
		v1 = set(v)
		juncDesc[k] = list(v1)

	#for k, v in juncGridDesc.iteritems():
	#	v1 = set(v)
	#	juncGridDesc[k] = list(v1)

				
	""" update to data structure """
	leaf2LeafPathJunctions["juncPoints"] = juncPoints
	leaf2LeafPathJunctions["juncArmPoints"] = juncArmPoints
	leaf2LeafPathJunctions["juncDesc"] = juncDesc

	
	if plotIter:
		pylab.clf()

		for path in medialLongPaths:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP)


		if globalJunctionPose != None:

			for path in theoryMedialLongPaths:
				xP = []
				yP = []
				for p in path:
					xP.append(p[0])
					yP.append(p[1])

				pylab.plot(xP,yP, color='k')

		#globalJunctionPose = self.getGlobalJunctionPose(pathID)
		if globalJunctionPose != None:
			pylab.scatter([globalJunctionPose[0],], [globalJunctionPose[1],], color='k')

		
		xP = []
		yP = []
		for p in vertices:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color='r')
		
		sizes = []
		for path in leafToLeafPaths:
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
		pylab.title("Path %d %d %s %s %s" % (hypothesisID, pathID, sizes,bufStr1,bufStr2))
		pylab.savefig("medialOut2_%02u_%04u.png" % (hypothesisID, topCount))
		print "saving medialOut2_%02u_%04u.png" % (hypothesisID, topCount)

		" 1) plot the pose of local splines and postures "
		" 2) plot the alpha shape of the union of pose alpha shapes and medial axis tree "
		" 3) plot the alpha shape of the union of pose alpha shapes and the long medial axis "
		" 4) plot the trimmed path "

		#fig = plt.figure()
		#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True)
		#ax1.plot(x)
		
		pylab.clf()
		#pathIDs = self.getPathIDs()

		#allNodes = []
		#for k in pathIDs:
		#	nodeSet = self.getNodes(k)
		#	allNodes += copy(nodeSet)
		
		#allNodes.sort() 
		
		#if len(allNodes) > 0:
		#	highestNodeID = allNodes[-1]
		#else:
		#	highestNodeID = 1e100		
			
			
		xP = []
		yP = []
		
		for path in medialLongPaths:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color=color, linewidth=4)

		if globalJunctionPose != None:

			for path in theoryMedialLongPaths:
				xP = []
				yP = []
				for p in path:
					xP.append(p[0])
					yP.append(p[1])

				pylab.plot(xP,yP, color='k')

		#nodeSet = self.getNodes(pathID)

		print "drawing pathID", pathID, "for nodes:", nodeSet
		for nodeID in nodeSet:
			xP = []
			yP = []

			estPose1 = nodePoses[nodeID]
	
			if poseData.isBowties[nodeID]:			
				hull1 = poseData.aHulls[nodeID]
				#medial1 = poseData.medialAxes[nodeID]
				#hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
			else:
				#hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
				hull1 = poseData.aHulls[nodeID]
	
			" set the origin of pose 1 "
			poseOrigin = Pose(estPose1)
	
			points = []
			for p in hull1:
				p1 = poseOrigin.convertLocalToGlobal(p)
				points.append(p1)
			
			for p in points:
				xP.append(p[0])
				yP.append(p[1])
			
			#if nodeID == highestNodeID:
			#	pylab.plot(xP,yP, color=(0,0,0))
			#elif nodeID == highestNodeID-1:
			#	pylab.plot(xP,yP, color=(0.5,0.5,0.5))
			#else:
			#	pylab.plot(xP,yP, color=color)
			pylab.plot(xP,yP, color=color)
				
		
		xP = []
		yP = []
		for p in vertices:
			xP.append(p[0])
			yP.append(p[1])
			
		pylab.plot(xP,yP, '--', color=color, linewidth=4)

		#globalJunctionPose = self.getGlobalJunctionPose(pathID)
		if globalJunctionPose != None:
			pylab.scatter([globalJunctionPose[0],], [globalJunctionPose[1],], color='k', zorder=10)

		#self.plotEnv()
		
		print "pathAndHull:", topCount

		pylab.axis("equal")
		pylab.title("Path %d %d" % (hypothesisID, pathID))
		#pylab.title("paths: %s numNodes: %d %d, hyp %d %3.2f" % (repr(mapHyp.getPathIDs()), poseData.numNodes, highestNodeID, mapHyp.hypothesisID, mapHyp.utility))
		pylab.savefig("pathAndHull_%02u_%04u.png" % (hypothesisID, topCount))

		#self.topCount += 1

	print "juncAngSet:", juncAngSet

	
	""" Select the leaf-to-leaf path that will be our defacto shoot """

	" searches for the branch from the parent junction "
	" selects splice of topology that has aligning junction "
	" does not select for distance or best medial axis representation of path "
	

	" sort by longest first "
	pathCands = []
	for k in range(len(leafToLeafPaths)):
		pathCands.append((len(leafToLeafPaths[k]), k))

	pathCands.sort(reverse=True)
	
	print "pathCands:", pathCands
	

	" TODO: find longest path going to right, longest path going to left through the junction "


	
	" FIXME:  Verify that leaf path exists for the junction point "
	
	" select longest path that has the best angle fit "
	branchArm = None
	if globalJunctionPose != None:
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
						branchArm = medialLongPaths[bestFit][jIndex+1]

					if fabs(angDiff2) < minDiff:
						minDiff = fabs(angDiff2)
						bestFit = k
						branchArm = medialLongPaths[bestFit][jIndex-1]

			pass



		if bestFit != -1:
			
			if fabs(minDiff) < 1.047:

				print "returning bestFit:", bestFit, minDiff
				leaf2LeafPathJunctions["bestFit"] = bestFit
				leaf2LeafPathJunctions["branchArm"] = branchArm
				leaf2LeafPathJunctions["longPaths"] = medialLongPaths
				leaf2LeafPathJunctions["theoryPaths"] = theoryMedialLongPaths
				return leaf2LeafPathJunctions, medialLongPaths[bestFit], vertices
			else:
				print "not returning bestFit"

			#else:
			#	" FIXME:  return theory junction information if this is selected "
			#	print "returning bestFit theory:", theoryJunc
			#	leaf2LeafPathJunctions["bestFit"] = 0
			#	leaf2LeafPathJunctions["branchArm"] = branchArm
			#	leaf2LeafPathJunctions["longPaths"] = medialLongPaths
			#	leaf2LeafPathJunctions["theoryPaths"] = theoryMedialLongPaths
			#	return leaf2LeafPathJunctions, theoryMedialLongPaths[0], vertices
		
		else:
			print "not returning bestFit"

	print "returning longest fit"

	maxIndex = 0
	maxLen = 0
	for k in range(len(leaf2LeafPathLengths)):
		if leaf2LeafPathLengths[k] > maxLen:
			maxIndex = k
			maxLen = leaf2LeafPathLengths[k]

	leaf2LeafPathJunctions["bestFit"] = maxIndex


	" FIXME: change from selecting any old junction point since "
	#jIndex = juncLongIndices[maxIndex][0]
	leaf2LeafPathJunctions["branchArm"] = medialLongPaths[maxIndex][1]

	leaf2LeafPathJunctions["longPaths"] = medialLongPaths
	leaf2LeafPathJunctions["theoryPaths"] = theoryMedialLongPaths

	return leaf2LeafPathJunctions, medialLongPaths[maxIndex], vertices
	#return medialLongPaths[maxIndex], vertices




@logFunction
def spliceSkeletons(localSkeletons, controlPoses, junctionPoses, parentPathIDs):

	SPLICE_DIST = 0.05

	globalSkeletons = {}
	junctionNodes = {}

	pathIDs = localSkeletons.keys()
	
	for pathID in pathIDs: 

		skel = localSkeletons[pathID]
		controlPose = controlPoses[pathID]
		juncPose = junctionPoses[pathID]
		parentPathID = parentPathIDs[pathID]

		origControlOrigin = Pose(controlPose)

		if parentPathID != None:
			globalJuncPose = origControlOrigin.convertLocalOffsetToGlobal(juncPose)
		else:
			globalJuncPose = None

		globalSkeletonGraph = graph.graph()

		minChildNode = None
		minChildDist = 1e100
	
		for edge in skel.edges():
		
			nodePoint1 = edge[0]
			nodePoint2 = edge[1]

			globalNodePoint1 = tuple(origControlOrigin.convertLocalToGlobal(nodePoint1))
			globalNodePoint2 = tuple(origControlOrigin.convertLocalToGlobal(nodePoint2))

			dist = sqrt((globalNodePoint1[0]-globalNodePoint2[0])**2 + (globalNodePoint1[1]-globalNodePoint2[1])**2)

			if parentPathID != None:
				dist1 = sqrt((globalNodePoint1[0]-globalJuncPose[0])**2 + (globalNodePoint1[1]-globalJuncPose[1])**2)
				dist2 = sqrt((globalNodePoint2[0]-globalJuncPose[0])**2 + (globalNodePoint2[1]-globalJuncPose[1])**2)

				if dist1 < minChildDist:
					minChildDist = dist1
					minChildNode = globalNodePoint1

				if dist2 < minChildDist:
					minChildDist = dist2
					minChildNode = globalNodePoint2


			globalSkeletonGraph.add_node(globalNodePoint1)
			globalSkeletonGraph.add_node(globalNodePoint2)
			globalSkeletonGraph.add_edge(globalNodePoint1, globalNodePoint2, wt=dist)

		if parentPathID != None:

			""" add junction node and edge """
			juncNode = (globalJuncPose[0], globalJuncPose[1])

			junctionNodes[pathID] = juncNode

			print "add child junc edge", juncNode, minChildNode, minChildDist
			globalSkeletonGraph.add_node(juncNode)
			globalSkeletonGraph.add_edge(juncNode, minChildNode, wt=minChildDist)


		else:
			junctionNodes[pathID] = None

		globalSkeletons[pathID] = globalSkeletonGraph


	spliceSkeleton = graph.graph()
	for k in range(len(globalSkeletons)):
		#nodes = globalSkeletons[k].nodes()
		edges = globalSkeletons[k].edges()
		#spliceSkeleton.add_nodes(nodes)

		for edge in edges:
			globalNodePoint1 = edge[0]
			globalNodePoint2 = edge[1]

			weight = globalSkeletons[k].get_edge_weight(globalNodePoint1, globalNodePoint2)

			#dist = sqrt((globalNodePoint1[0]-globalNodePoint2[0])**2 + (globalNodePoint1[1]-globalNodePoint2[1])**2)
			spliceSkeleton.add_node(globalNodePoint1)
			spliceSkeleton.add_node(globalNodePoint2)
			spliceSkeleton.add_edge(globalNodePoint1, globalNodePoint2, wt=weight)


	""" connect the children to their parents via the junction node """
	for pathID in pathIDs:
		juncNode = junctionNodes[pathID]
		parentPathID = parentPathIDs[pathID]

		if parentPathID != None and juncNode != None:
			parentSkel = globalSkeletons[parentPathID]

			minParentdNode = None
			minParentDist = 1e100
		
			for edge in parentSkel.edges():
			
				globalNodePoint1 = edge[0]
				globalNodePoint2 = edge[1]

				dist1 = sqrt((globalNodePoint1[0]-juncNode[0])**2 + (globalNodePoint1[1]-juncNode[1])**2)
				dist2 = sqrt((globalNodePoint2[0]-juncNode[0])**2 + (globalNodePoint2[1]-juncNode[1])**2)

				if dist1 < minParentDist:
					minParentDist = dist1
					minParentNode = globalNodePoint1

				if dist2 < minParentDist:
					minParentDist = dist2
					minParentNode = globalNodePoint2

			print "add parent junc edge", juncNode, minParentNode, minParentDist
			spliceSkeleton.add_edge(juncNode, minParentNode, wt=minParentDist)
	

	for j in range(len(globalSkeletons)):
		for k in range(j, len(globalSkeletons)):
			nodes1 = globalSkeletons[j].nodes()
			nodes2 = globalSkeletons[k].nodes()

			distances1, indices1 = getClosestPairs(nodes1, nodes2)
			distances2, indices2 = getClosestPairs(nodes2, nodes1)

			for m in range(len(distances1)):
				dist1 = distances1[m]

				node1 = nodes1[m]
				node2 = nodes2[indices1[m]]

				if dist1 <= SPLICE_DIST:
					#print "add_edge(", node1, node2, dist1
					spliceSkeleton.add_edge(node1, node2, wt=dist1)
	
			for m in range(len(distances2)):
				dist2 = distances2[m]

				node2 = nodes2[m]
				node1 = nodes1[indices2[m]]

				if dist2 <= SPLICE_DIST:
					#print "add_edge(", node1, node2, dist2
					spliceSkeleton.add_edge(node2, node1, wt=dist2)
		
	return spliceSkeleton

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
def getSimpleDeparture(globalJunctionPose, path1, path2, plotIter = False):

	" return exception if we receive an invalid path "		  
	if len(path1) == 0:
		print "path1 has zero length"
		raise

	" return exception if we receive an invalid path "		  
	if len(path2) == 0:
		print "path2 has zero length"
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

	DIST_THRESH = 0.1

	""" find the point that first crosses distance threshold """
	currI = 0
	try:
		while distances[currI] > DIST_THRESH:
			currI += 1
	except:
		""" index out of bounds, back track """
		currI -= 1
	
	frontPoint = pathPoints1[indices[currI]]

	currI = len(distances)-1
	try:
		while distances[currI] > DIST_THRESH:
			currI -= 1
	except:
		""" index out of bounds, back track """
		currI += 1

	backPoint = pathPoints1[indices[currI]]
	
	juncDist1 = sqrt((globalJunctionPose[0]-frontPoint[0])**2 + (globalJunctionPose[1]-frontPoint[1])**2)
	juncDist2 = sqrt((globalJunctionPose[0]-backPoint[0])**2 + (globalJunctionPose[1]-backPoint[1])**2)


	if juncDist1 < juncDist2:
		return frontPoint
	else:
		return backPoint


@logFunction
def getSimpleSoupDeparture(globalJunctionPose, pointSoup, path2, plotIter = False):

	" return exception if we receive an invalid path "		  
	if len(path2) == 0:
		print "path2 has zero length"
		raise

	pointSoupTree = cKDTree(array(pointSoup))

	" for each point on the child path, find its closest pair on the parent path "
	path2Spline = SplineFit(path2, smooth=0.1)
	pathPoints2 = path2Spline.getUniformSamples()

	distances = []
	indices = []
	for i in range(0,len(pathPoints2)):
		p_2 = pathPoints2[i]
		queryPoint = p_2[0:2]
		minDist, i_1 = pointSoupTree.query(array(queryPoint))
		#p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
		
		" keep the distance information "
		distances.append(minDist)
		" and the associated index of the point on the parent path "
		indices.append(i_1)

	DIST_THRESH = 0.1

	""" find the point that first crosses distance threshold """
	currI = 0
	try:
		while distances[currI] > DIST_THRESH:
			currI += 1
	except:
		""" index out of bounds, back track """
		currI -= 1
	
	frontPoint = pointSoup[indices[currI]]
	frontIndex = currI

	currI = len(distances)-1
	try:
		while distances[currI] > DIST_THRESH:
			currI -= 1
	except:
		""" index out of bounds, back track """
		currI += 1

	backPoint = pointSoup[indices[currI]]
	backIndex = currI
	
	juncDist1 = sqrt((globalJunctionPose[0]-frontPoint[0])**2 + (globalJunctionPose[1]-frontPoint[1])**2)
	juncDist2 = sqrt((globalJunctionPose[0]-backPoint[0])**2 + (globalJunctionPose[1]-backPoint[1])**2)

	#pathPoints2[0:frontIndex+1]
	frontIndexEnd = frontIndex-20
	if frontIndexEnd < 0:
		frontIndexEnd = 0
	if frontIndex-frontIndexEnd > 0:
		frontVecSum = [0.0,0.0]

		revRange = range(frontIndexEnd+1, frontIndex+1)
		revRange.reverse()
		for k in revRange:
			p1 = pathPoints2[k]
			p2 = pathPoints2[k-1]
			frontVecSum[0] += p2[0]-p1[0]
			frontVecSum[1] += p2[1]-p1[1]
		
		vecMag = sqrt(frontVecSum[0]**2 + frontVecSum[1]**2)

		frontVec = [frontVecSum[0]/vecMag, frontVecSum[1]/vecMag]

		frontVecAngle = acos(frontVec[0])
		if frontVec[1] < 0.0:
			frontVecAngle = -frontVecAngle

		frontPose = [frontPoint[0], frontPoint[1], frontVecAngle]
	else:
		frontPose = None


	#pathPoints2[backIndex:len(distances)]

	backIndexEnd = backIndex+20
	if backIndexEnd >= len(distances):
		backIndexEnd = len(distances)-1
	if backIndexEnd-backIndex > 0:
		backVecSum = [0.0,0.0]
		for k in range(backIndex, backIndexEnd):
			p1 = pathPoints2[k]
			p2 = pathPoints2[k+1]
			backVecSum[0] += p2[0]-p1[0]
			backVecSum[1] += p2[1]-p1[1]
		
		vecMag = sqrt(backVecSum[0]**2 + backVecSum[1]**2)

		backVec = [backVecSum[0]/vecMag, backVecSum[1]/vecMag]

		backVecAngle = acos(backVec[0])
		if backVec[1] < 0.0:
			backVecAngle = -backVecAngle

		backPose = [backPoint[0], backPoint[1], backVecAngle]
	else:
		backPose = None


	print "getSimpleSoupDeparture:", frontPoint, backPoint, frontPose, backPose, juncDist1, juncDist2, len(pathPoints2), frontIndex, backIndex
	if juncDist1 < juncDist2 and frontPose != None:
		return frontPose
	elif juncDist1 > juncDist2 and backPose != None:
		return backPose
	else:
		raise


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
def getInitSkeletonBranchPoint(globalJunctionPose, currShootID, globalMedial_G, parentShootIDs, localPathSegsByID, globalControlPoses, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):

	"""
	1)  get the control point, common length between current and most senior parent 
	2)  branch point is the part of current skeleton that diverges from all non-descendent skeletons

	"""

	allPointSoup_G = []

	CONTIG_DIST = 0.2

	allPathIDs = localPathSegsByID.keys()

	pathSegsByID_G = []
	for shootID in allPathIDs:

		localSegs_L = localPathSegsByID[shootID]
		controlPose_G = globalControlPoses[shootID]

		localSegs_G = []
		descenFrame = Pose(controlPose_G)
		for seg in localSegs_L:
			newSeg = []
			for p in seg:
				newPose = descenFrame.convertLocalOffsetToGlobal(p)
				allPointSoup_G.append(newPose[0:2])
				newSeg.append(newPose)
			localSegs_G.append(newSeg)

		pathSegsByID_G.append(localSegs_G)



	maxSkeletonID = 0
	maxSkeletonContigFrac = 0.0
	maxSkeletonSegID = 0

	for shootID in allPathIDs:
		shootSegs_G = pathSegsByID_G[shootID]

		maxContigFrac = 0.0
		maxContigIndex = 0

		for segIndex in range(len(shootSegs_G)):

			seg = shootSegs_G[segIndex]

			""" make sure the overlap of both shoots are oriented the same way """
			orientedMedial_G = orientPath(globalMedial_G, seg)
	
			kdInput = []
			for j in range(0,len(seg)):
				p_1 = seg[j]
				kdInput.append(p_1[0:2])

			kdTree = cKDTree(array(kdInput))
			
			contigCount = 0
			maxContig = 0
			for i in range(0,len(orientedMedial_G)):

				p_2 = orientedMedial_G[i]

				queryPoint = p_2[0:2]
				minDist, i_1 = kdTree.query(array(queryPoint))
				p_1 = seg[i_1]

				if minDist < CONTIG_DIST:
					contigCount += 1
					if contigCount > maxContig:
						maxContig = contigCount
				else:
					contigCount = 0

			contigFrac = float(maxContig)/float(len(orientedMedial_G))

			if contigFrac > maxContigFrac:
				maxContigFrac = contigFrac
				maxContigIndex = segIndex
			
		if maxContigFrac > maxSkeletonContigFrac:
			maxSkeletonID = shootID
			maxSkeletonSegID = maxContigIndex
			maxSkeletonContigFrac = maxContigFrac


	controlSeg = pathSegsByID_G[maxSkeletonID][maxSkeletonSegID]
	orientedMedial_G = orientPath(globalMedial_G, controlSeg)
	
	""" compute the control point """
	commonU1, commonU2, commonP1, commonP2 = selectCommonOrigin(controlSeg, orientedMedial_G)



	""" get branch point from all skeletons """
	branchPose_G = getSimpleSoupDeparture(globalJunctionPose, allPointSoup_G, globalMedial_G)

	""" convert control pose from global to parent frame """
	controlPose_G = commonP1
	controlPose_G[2] = 0.0
	controlParentID = maxSkeletonID
	#parentControlPose_G = globalControlPoses[controlParentID]
	#shootFrame = Pose(parentControlPose_G)
	#controlPose_P = shootFrame.convertGlobalPoseToLocal(controlPose_G)

	#currFrame = Pose(controlPose_G)
	#branchPose_C = currFrame.convertGlobalPoseToLocal(branchPose_G)

	""" placeholder values, not true branch point """
	#branchParentID = controlParentID

	print "skeletonBranchPoint:", maxSkeletonID, maxSkeletonSegID, maxSkeletonContigFrac, controlPose_G, branchPose_G

	return controlPose_G, controlParentID, branchPose_G
	


@logFunction
def getSkeletonBranchPoint(globalJunctionPose, currShootID, parentShootIDs, localPathSegsByID, globalControlPoses, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):

	"""
	2)  branch point is the part of current skeleton that diverges from all non-descendent skeletons
	"""

	""" bookkeeping, find all descendants of current shoot """
	isDescendant = {}

	for key in parentShootIDs.keys():
		isDescendant[key] = False

	currID = currShootID
	currDescen = [currShootID]
	while True:
		
		for val in currDescen:
			isDescendant[val] = True

		newDescen = [currShootID]

		for key, val in parentShootIDs.iteritems():
			if val in currDescen:
				newDescen.append(key)

		if len(currDescen) == len(newDescen):
			break

		currDescen = list(set(currDescen + newDescen))


	""" we are potentially branching or controlling from any of the non-descendant skeletons """
	nonDescentIDs = []
	for shootID, isDescen in isDescendant.iteritems():
		if not isDescen:
			nonDescentIDs.append(shootID)
	

	currSegs_L = localPathSegsByID[currShootID]
	currControlPose_G = globalControlPoses[currShootID]

	""" convert to global frame """
	currSegs_G = []
	currFrame = Pose(currControlPose_G)
	for seg in currSegs_L:
		newSeg = []
		for p in seg:
			newSeg.append(currFrame.convertLocalOffsetToGlobal(p))
		currSegs_G.append(newSeg)


	allPointSoup_G = []
	CONTIG_DIST = 0.2

	pathSegsByID_G = []
	for shootID in nonDescentIDs:

		localSegs_L = localPathSegsByID[shootID]
		controlPose_G = globalControlPoses[shootID]

		localSegs_G = []
		descenFrame = Pose(controlPose_G)
		for seg in localSegs_L:
			newSeg = []
			for p in seg:
				newPose = descenFrame.convertLocalOffsetToGlobal(p)
				allPointSoup_G.append(newPose[0:2])
				newSeg.append(newPose)
			localSegs_G.append(newSeg)

		pathSegsByID_G.append(localSegs_G)

	"""
	2)  branch point is the part of current skeleton that diverges from all non-descendent skeletons
	"""

	minDist = 1e100
	minBranchPose = None
	for k in range(len(currSegs_G)):

		seg_G = currSegs_G[k]

		try:
			""" get branch point from all skeletons """
			branchPose_G = getSimpleSoupDeparture(globalJunctionPose, allPointSoup_G, seg_G)

			dist = sqrt((branchPose_G[0]-globalJunctionPose[0])**2 + (branchPose_G[1]-globalJunctionPose[1])**2)

			if dist < minDist:
				minDist = dist
				minBranchPose = branchPose_G

		except:
			branchPose_G = None



	print "skeletonBranchPoint:", minBranchPose

	if minBranchPose == None:
		raise

	return minBranchPose
	
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
def computeBranch(pathID, parentID, oldLocalJuncPose_C, origControlPose, childPath_C, parentPath_P, trimmedParent_P, localPathSegsByID, arcDist, localSkeletons, controlPoses, junctionPoses, parentPathIDs, numNodes=0, hypothesisID=0):

	if parentID != None: 

		""" parent path information """
		pathSpline_P = SplineFit(parentPath_P)

		""" f(arcDist) = position on the spline curve of the parent path """
		""" NOTE:  angle is the tangent angle, not branch angle """
		""" LOCATION OF MODIFIED CONTROL POINT ON PARENT """
		newArcDist = arcDist
		newControlPose_P = pathSpline_P.getPointOfDist(newArcDist)
		modControlPose_P = copy(newControlPose_P)
		modControlPose_P[2] = 0.0


		""" old junction point in parent frame """
		offsetOrigin1 = Pose(modControlPose_P)
		modJuncPose_P = offsetOrigin1.convertLocalOffsetToGlobal(oldLocalJuncPose_C)


		""" place the segments back into parent frame with control point from arc distance """

		localPathSegs = localPathSegsByID[pathID]
		offsetOrigin1 = Pose(modControlPose_P)
		placedPathSegs = []
		for k in range(len(localPathSegs)):
			localSeg = localPathSegs[k]
			placedSeg = []
			for p in localSeg:
				p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)
				placedSeg.append(p1)
			placedPathSegs.append(placedSeg)

		# FIXME:  add in parent skeleton, all other non-descendant skeletons

		""" evaluate overlap of parent and child skeleton """
		lastCost1, matchCount1 = gen_icp.branchEstimateCost2(modControlPose_P, placedPathSegs, parentPath_P, plotIter = True, n1 = pathID, n2 = parentID, arcDist = arcDist)

		""" originally designed to compute the discrepancy between initial and final poses of ICP,
			now there is no ICP so there is no discrepancy """
		distDisc = 0.0
		angDisc = 0.0

		finalPoses_G = computeGlobalControlPoses(controlPoses, parentPathIDs)

		""" get the trimmed child shoot at the new designated branch point from parent """
		newPath3, localJuncPose_C, juncDiscDist, juncDiscAngle, splicedPaths, particlePath =  trimBranch(pathID, parentID, modControlPose_P, oldLocalJuncPose_C, childPath_C, parentPath_P, trimmedParent_P, localPathSegsByID, parentPathIDs, finalPoses_G, plotIter=True, arcDist = arcDist, nodeID=numNodes, hypothesisID = hypothesisID)


		newGlobJuncPose = offsetOrigin1.convertLocalOffsetToGlobal(localJuncPose_C)

		#junctionPoses[pathID] = newGlobJuncPose
		junctionPoses[pathID] = localJuncPose_C
		controlPoses[pathID] = modControlPose_P

		finalPoses = computeGlobalControlPoses(controlPoses, parentPathIDs)

		spliceSkeleton_G = spliceSkeletons(localSkeletons, finalPoses, junctionPoses, parentPathIDs)


		""" the terminals of the parent shoot """
		parentTerm1 = parentPath_P[0]
		parentTerm2 = parentPath_P[-1]

		""" following terminal comparisons are in local coordinates """

		""" the terminals of the child shoot """
		dist1 = sqrt((localJuncPose_C[0]-newPath3[0][0])**2 + (localJuncPose_C[1]-newPath3[0][1])**2)
		dist2 = sqrt((localJuncPose_C[0]-newPath3[-1][0])**2 + (localJuncPose_C[1]-newPath3[-1][1])**2)

		if dist1 < dist2:
			childTerm = newPath3[-1]
		else:
			childTerm = newPath3[0]

		dist3 = sqrt((childTerm[0]-particlePath[0][0])**2 + (childTerm[1]-particlePath[0][1])**2)
		dist4 = sqrt((childTerm[0]-particlePath[-1][0])**2 + (childTerm[1]-particlePath[-1][1])**2)

		if dist3 < dist4:
			childTerm2 = particlePath[-1]
		else:
			childTerm2 = particlePath[0]

		""" two different splices between child and parent terminals """
		globalChildTerm1 = offsetOrigin1.convertLocalToGlobal(childTerm)
		globalChildTerm2 = offsetOrigin1.convertLocalToGlobal(childTerm2)
		termPaths = [(parentTerm1, globalChildTerm1), (parentTerm2, globalChildTerm1), (globalChildTerm2, globalChildTerm1)]

		splicedPaths = []
		for termPath in termPaths:

			""" use splice skeleton to find shortest path """
			startPose = termPath[0]
			endPose = termPath[-1]

			minStartDist = 1e100
			minStartNode = None
			minEndDist = 1e100
			minEndNode = None
			
			for edge in spliceSkeleton_G.edges():
			
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


			shortestSpliceTree, shortestSpliceDist = spliceSkeleton_G.shortest_path(endNode)
			currNode = shortestSpliceTree[startNode]					 
			splicedSkel = [startNode]
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
		#part2["modJuncPose"] = newGlobJuncPose
		part2["modJuncPose"] = localJuncPose_C
		part2["modControlPose"] = modControlPose_P
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
		part2["termPaths"] = termPaths

		#if plotIter:
		if True:

			#numNodes = 0
			#pathPlotCount = 0

			
			pylab.clf()

			#for splicePath in newSplices:

			#	xP = []
			#	yP = []
			#	for p in splicePath:
			#		xP.append(p[0])
			#		yP.append(p[1])
			#	pylab.plot(xP,yP, color='k', alpha=0.8)


			xP = []
			yP = []
			for p in newPath3:
				p1 = offsetOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])
				pylab.plot(xP,yP, color='k', alpha=0.8)

			xP = []
			yP = []
			for p in particlePath:
				p1 = offsetOrigin1.convertLocalToGlobal(p)
				xP.append(p1[0])
				yP.append(p1[1])
				pylab.plot(xP,yP, color='m', alpha=0.5)

			for edge in spliceSkeleton_G.edges():
			
				globalNodePoint1 = edge[0]
				globalNodePoint2 = edge[1]

				xP = [globalNodePoint1[0], globalNodePoint2[0]]
				yP = [globalNodePoint1[1], globalNodePoint2[1]]

				pylab.plot(xP,yP, color='k', alpha=0.2)

			for termPath in termPaths:
				startPose = termPath[0]
				endPose = termPath[-1]

				xP = [startPose[0], endPose[0]]
				yP = [startPose[1], endPose[1]]

				pylab.scatter(xP,yP, color='r')

			pylab.scatter([newGlobJuncPose[0],], [newGlobJuncPose[1],], color='m')
			pylab.scatter([modControlPose_P[0],], [modControlPose_P[1],], color='b')

			pylab.axis("equal")
			pylab.title("hyp %d nodeID %d %1.2f %1.2f %1.2f" % ( hypothesisID, numNodes, juncDiscDist, juncDiscAngle, newGlobJuncPose[2] ))
			pylab.savefig("computeBranch_%04u_%04u_%1.1f.png" % (hypothesisID, numNodes, arcDist))

			print "saving computeBranch_%04u_%04u_%1.1f.png" % (hypothesisID, numNodes, arcDist)

			
		#pathPlotCount += 1

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
def trimBranch(pathID, parentPathID, modControlPose, origGlobJuncPose, childPath, parentPath, trimmedParent, localPathSegsByID, parentPathIDs, controlPoses_G, plotIter = False, hypothesisID = 0, nodeID = 0, arcDist = 0.0):

	global pathPlotCount 





	#origControlOrigin = Pose(origControlPose)

	#localPathSegs = []

	#for k in range(len(smoothPathSegs)):
	#	pathSeg = smoothPathSegs[k]
	#	localSeg = []
	#	for p in pathSeg:
	#		p1 = origJuncOrigin.convertGlobalPoseToLocal(p)
	#		localSeg.append(p1)

	#	localPathSegs.append(localSeg)


	offsetOrigin1 = Pose(modControlPose)

	globJuncPose = offsetOrigin1.convertLocalOffsetToGlobal(origGlobJuncPose)
	
	partPathSegs = []

	#for k in range(len(localPathSegs)):
	#	localSeg = localPathSegs[k]
	#	for p in localSeg:
	#		p1 = offsetOrigin1.convertLocalOffsetToGlobal(p)

	path1 = parentPath
	path2 = childPath

	#localPath2 = []
	#for p in path2:
	#	p1 = origControlOrigin.convertGlobalToLocal(p)
	#	localPath2.append(p1)

	particlePath2 = []
	#for p in localPath2:
	for p in path2:
		p1 = offsetOrigin1.convertLocalToGlobal(p)
		particlePath2.append(p1)

	""" get departing sections of overlapped curves """
	""" branch point as input """
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
		dist = sqrt((pnt[0]-globJuncPose[0])**2 + (pnt[1]-globJuncPose[1])**2)
	
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


	branchPose_G = getSkeletonBranchPoint(globJuncPose, pathID, parentPathIDs, localPathSegsByID, controlPoses_G)

	newGlobJuncPose1, controlPoint1, angDeriv1 = getBranchPoint(globJuncPose, parentPathID, pathID, trimmedParent, particlePath2, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)
	newGlobJuncPose2, controlPoint2, angDeriv2 = getBranchPoint(globJuncPose, pathID, parentPathID, particlePath2, trimmedParent, plotIter = plotIter, hypothesisID = hypothesisID, nodeID = nodeID, arcDist = arcDist)


	juncDiscDist1 = sqrt((newGlobJuncPose1[0]-globJuncPose[0])**2 + (newGlobJuncPose1[1]-globJuncPose[1])**2)
	juncDiscDist2 = sqrt((newGlobJuncPose2[0]-globJuncPose[0])**2 + (newGlobJuncPose2[1]-globJuncPose[1])**2)

	if angDeriv1 < angDeriv2:
		newGlobJuncPose = newGlobJuncPose1
		angDeriv = angDeriv1
	else:

		""" only diverge from the parent if the determined branch point is in the neighborhood of the original """
		if juncDiscDist2 < 1.0:
			newGlobJuncPose = newGlobJuncPose2
			angDeriv = angDeriv2
		else:
			newGlobJuncPose = newGlobJuncPose1
			angDeriv = angDeriv1

	""" override this with simple getSimpleDeparture() """
	#newGlobJuncPose = junctionPoint_K

	newGlobJuncPose = getSimpleDeparture(globJuncPose, path1, particlePath2)

	""" regardless, take the branch angle of the child from parent shoot branch point """
	newGlobJuncPose[2] = newGlobJuncPose1[2]


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

		xP = []
		yP = []
		for p in particlePath2:
			xP.append(p[0])
			yP.append(p[1])
			pylab.plot(xP,yP, color='m', alpha=0.5)

		
		pylab.scatter([newGlobJuncPose[0],], [newGlobJuncPose[1],], color='r')
		pylab.scatter([modControlPose[0],], [modControlPose[1],], color='b')
		#pylab.scatter([origControlPose[0],], [origControlPose[1],], color='g')
		pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='m')

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

	""" convert back to local coordinate system for each shoot """

	localNewPath3 = []
	for p in newPath3:
		p1 = offsetOrigin1.convertGlobalToLocal(p)
		localNewPath3.append(p1)

	localParticlePath2 = []
	for p in particlePath2:
		p1 = offsetOrigin1.convertGlobalToLocal(p)
		localParticlePath2.append(p1)

	localJuncPose = offsetOrigin1.convertGlobalPoseToLocal(newGlobJuncPose)

	localSplicedPaths = []
	for splicePath in splicedPaths:
		localSplicePath = []
		for p in splicePath:
			p1 = offsetOrigin1.convertGlobalToLocal(p)
			localSplicePath.append(p1)
		localSplicedPaths.append(localSplicePath)

	#return newPath3, newGlobJuncPose, juncDiscDist, juncDiscAngle, splicedPaths, particlePath2
	return localNewPath3, localJuncPose, juncDiscDist, juncDiscAngle, localSplicedPaths, localParticlePath2





