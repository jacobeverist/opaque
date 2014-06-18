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

	#print "branchNodeID:", branchNodeID
	print "junctions:", junctions
	print "leaves:", leaves

	return vertices, junctions, leaves, uni_mst, gridHash


def computePathSegments(juncIDs, leafIDs, tree, gridHash, vertices, controlPose):

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
	for seg in realLeafSegments:

		realSeg = extendBackToHull(seg, vertices)

		""" the extended version of the leaf segment """
		longLeafSegments.append(realSeg)


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

	origControlOrigin = Pose(controlPose)
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

	return smoothLeafSegments, smoothInternalSegments, localSkeletonGraph

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
def computeShootSkeleton(poseData, pathID, branchNodeID, globalJunctionPose, controlPose, nodeSet, nodePoses, hypothesisID, color, topCount):
	
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
	if branchNodeID != None:

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
	smoothLeafSegments, smoothInternalSegments, localSkeletonGraph = computePathSegments(junctions, leaves, uni_mst, gridHash, vertices, controlPose)

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


	if branchNodeID != None:


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

	if globalJunctionPose != None:
		#origJuncPose = copy(controlPose)
		#origJuncPose = copy(globalJunctionPose)
		#origJuncPose[2] = 0.0
		#origJuncOrigin = Pose(origJuncPose)
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

	leaf2LeafPathJunctions["localSegments"] = localPathSegs


	""" get the neighbor points along the path at the junction so we can describe the occupancy state """
	juncPoints = []
	juncArmPoints = []
	juncDesc = {}
	juncGridPoints = []
	juncGridArmPoints = []
	juncGridDesc = {}
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

		#juncGridPoints = []
		#juncGridArmPoints = []
		#juncGridDesc = {}
		#leafToLeafGridPaths = []
		juncGridInds = juncIndices[k]
		juncGridPnts = []
		juncGridArmPnts = []
		for index in juncGridInds:
			jPnt = copy(leafToLeafGridPaths[k][index])

			jaPnt1 = leafToLeafGridPaths[k][index+1]
			jaPnt2 = leafToLeafGridPaths[k][index-1]

			jPnt = tuple(jPnt)
			jaPnt1 = tuple(jaPnt1)
			jaPnt2 = tuple(jaPnt2)

			juncGridPnts.append(jPnt)
			juncGridArmPnts.append((jaPnt1, jaPnt2))

			try:
				juncGridDesc[jPnt].append(jaPnt1)
				juncGridDesc[jPnt].append(jaPnt2)
			except:
				juncGridDesc[jPnt] = []
				juncGridDesc[jPnt].append(jaPnt1)
				juncGridDesc[jPnt].append(jaPnt2)


		juncGridPoints.append(juncGridPnts)
		juncGridArmPoints.append(juncGridArmPnts)


	""" remove duplicates """
	for k, v in juncDesc.iteritems():
		v1 = set(v)
		juncDesc[k] = list(v1)

	for k, v in juncGridDesc.iteritems():
		v1 = set(v)
		juncGridDesc[k] = list(v1)

				
	""" update to data structure """
	leaf2LeafPathJunctions["juncPoints"] = juncPoints
	leaf2LeafPathJunctions["juncArmPoints"] = juncArmPoints
	leaf2LeafPathJunctions["juncDesc"] = juncDesc

	
	if True:
		pylab.clf()

		for path in medialLongPaths:
			xP = []
			yP = []
			for p in path:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP)


		if branchNodeID != None:

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

		if branchNodeID != None:

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
				#return medialLongPaths[bestFit], vertices
				#return theoryMedialLongPaths, medialLongPaths, leaf2LeafPathJunctions, medialLongPaths[bestFit], vertices

			else:
				" FIXME:  return theory junction information if this is selected "
				print "returning bestFit theory:", theoryJunc
				leaf2LeafPathJunctions["bestFit"] = 0
				leaf2LeafPathJunctions["branchArm"] = branchArm
				leaf2LeafPathJunctions["longPaths"] = medialLongPaths
				leaf2LeafPathJunctions["theoryPaths"] = theoryMedialLongPaths
				return leaf2LeafPathJunctions, theoryMedialLongPaths[0], vertices
				#return theoryMedialLongPaths[0], vertices
				#return theoryMedialLongPaths, medialLongPaths, leaf2LeafPathJunctions, theoryMedialLongPaths[0], vertices
		
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
				dist1 = sqrt((globalNodePoint1[0]-juncPose[0])**2 + (globalNodePoint1[1]-juncPose[1])**2)
				dist2 = sqrt((globalNodePoint2[0]-juncPose[0])**2 + (globalNodePoint2[1]-juncPose[1])**2)

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
			juncNode = (juncPose[0], juncPose[1])

			junctionNodes[pathID] = juncNode

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

		if parentPathID != None:
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
					spliceSkeleton.add_edge(node1, node2, wt=dist1)
	
			for m in range(len(distances2)):
				dist2 = distances2[m]

				node2 = nodes2[m]
				node1 = nodes1[indices2[m]]

				if dist2 <= SPLICE_DIST:
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






