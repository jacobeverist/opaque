import sys
import os

from LocalOccMap import *
from LocalBoundaryMap import *
from LocalObstacleMap import *
from SplineFit import *
from Pose import Pose
#from cornerDetection import extractCornerCandidates
from voronoi import Site, computeVoronoiDiagram

import estimateMotion
from GPACurve import GPACurve
from StableCurve import StableCurve
#from gen_icp import computeMedialAxis, getLongestPath
from medialaxis import computeMedialAxis

import graph
import alphamod

import random
import functions
from functions import decimatePoints
#import Image
from numpy import array, dot, transpose

estPlotCount = 0
alphaPlotCount = 0
splineCount = 0

def computeBareHull(node1, sweep = False, static = False):
	
	if static:
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		a_data = decimatePoints(a_data)

		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
	else:
				
		" Read in data of Alpha-Shapes without their associated covariances "
		node1.computeAlphaBoundary(sweep = sweep)
		a_data = node1.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
	
	return a_data_GPAC


def computeHullAxis(nodeID, node2, tailCutOff = False):

	medial2 = node2.getBestMedialAxis()

	if tailCutOff:
		medial2 = node2.medialTailCuts[0]
	else:
		medial2 = node2.medialLongPaths[0]


	if node2.isBowtie:			
		hull2 = computeBareHull(node2, sweep = False, static = True)
		hull2.append(hull2[0])

	else:
		hull2 = computeBareHull(node2, sweep = False)
		hull2.append(hull2[0])
	
	
	return hull2, medial2

"""
def computeHullAxis(nodeID, node2, tailCutOff = False):
	
	if node2.isBowtie:			
		hull2 = computeBareHull(node2, sweep = False, static = True)
		hull2.append(hull2[0])
		medial2 = node2.getStaticMedialAxis()

	else:
		hull2 = computeBareHull(node2, sweep = False)
		hull2.append(hull2[0])
		medial2 = node2.getMedialAxis(sweep = False)

	" take the long length segments at tips of medial axis"
	edge1 = medial2[0:2]
	edge2 = medial2[-2:]
	
	frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
	backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
	frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
	
	frontVec[0] /= frontMag
	frontVec[1] /= frontMag
	backVec[0] /= backMag
	backVec[1] /= backMag
	
	" make a smaller version of these edges "
	newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
	newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

	edge1 = [newP1, edge1[1]]
	edge2 = [edge2[0], newP2]

	
	" find the intersection points with the hull "
	interPoints = []
	for k in range(len(hull2)-1):
		hullEdge = [hull2[k],hull2[k+1]]
		isIntersect1, point1 = Intersect(edge1, hullEdge)
		if isIntersect1:
			interPoints.append(point1)
			break

	for k in range(len(hull2)-1):
		hullEdge = [hull2[k],hull2[k+1]]
		isIntersect2, point2 = Intersect(edge2, hullEdge)
		if isIntersect2:
			interPoints.append(point2)
			break
	
	" replace the extended edges with a termination point at the hull edge "			
	medial2 = medial2[1:-2]
	if isIntersect1:
		medial2.insert(0, point1)
	if isIntersect2:
		medial2.append(point2)
	

	" cut off the tail of the non-sweeping side "
	TAILDIST = 0.5

	if tailCutOff:
		
		if nodeID % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]
	
		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
			
	return hull2, medial2
"""

def getLongestPath(node, currSum, currPath, tree, isVisited, nodePath, nodeSum):

	if isVisited[node]:
		return

	isVisited[node] = 1
	nodePath[node] = copy(currPath) + [node]

	if nodeSum[node] < currSum:
		nodeSum[node] = currSum

	for childNode in tree[node]:
		getLongestPath(childNode, currSum + 1, nodePath[node], tree, isVisited, nodePath, nodeSum)

	isVisited[node] = 0


		

	def computeVoronoi(boundaryPoints):
		
		'''
		Compute the voronoi diagram using a 3rd party library.
		
		Voronoi diagram data structures
		----
		
		self.vertices, 2-tuples
		self.lines, 3-tuples: a,b,c for a*x + b*y = c
		self.edges, 3-tuples, (l, v1, v2)  line, vertex 1, vertex 2
		vertex -1 means edge extends to infinity
		'''

		sites = []
		for p in boundaryPoints:
			sites.append(Site(p[0],p[1]))

		vertices = []
		lines = []
		edges = []			
		if sites != []:
			vertices, lines, edges = computeVoronoiDiagram(sites)

		return vertices, lines, edges

	def pruneEdges(occMap, vertices, lines, edges):
		
		'''remove edges and vertices that are outside the explored area'''

		#NEIGHBOR_WIDTH = 3
		NEIGHBOR_WIDTH = 0

		roadGraph = graph.graph()	
		
		xSize, ySize = numPixel, numPixel
		#occPix = self.occMapImage.load()

		#occPix = self.occMap.image
		occPix = occMap.load()

		newVertices = {}
		newEdges = []

		# checking all edges to check within the free space
		for edge in edges:
			v1 = edge[1]
			v2 = edge[2]

			if v1 != -1 and v2 != -1:

				# get grid value of real point
				i1, j1 = realToGrid(vertices[v1])
				i2, j2 = realToGrid(vertices[v2])
				
				if i1 < xSize and i1 >= 0 and i2 < xSize and i2 >= 0 and j1 < ySize and j1 >= 0 and j2 < ySize and j2 >= 0 :
					seg = [(i1,j1),(i2,j2)]

					indX1 = seg[0][0]
					indY1 = seg[0][1]
					indX2 = seg[1][0]
					indY2 = seg[1][1]

					isOK = True

					for i in range(-NEIGHBOR_WIDTH,NEIGHBOR_WIDTH+1,1):
						for j in range(-NEIGHBOR_WIDTH,NEIGHBOR_WIDTH+1,1):

							if indX1+i > numPixel or indX1+i < 0 or indY1+j > numPixel or indY1+j < 0:
								isOK = False
							else:
								val1 = occPix[indX1 + i, indY1 + j]
								if val1 <= 127: # 0 or 127
									isOK = False

							if indX2+i > numPixel or indX2+i < 0 or indY2+j > numPixel or indY2+j < 0:
								isOK = False
							else:
								val2 = occPix[indX2 + i, indY2 + j]
								if val2 <= 127: # 0 or 127
									isOK = False

					if isOK:
						# 1. add node if not already added
						# 2. add edge

						newVertices[v1] = vertices[v1]
						newVertices[v2] = vertices[v2]

						newEdges.append([v1,v2])

		for key, attr in newVertices.iteritems():
			roadGraph.add_node(key, attr)

		for edge in newEdges:
			roadGraph.add_edge(edge[0], edge[1])

			# assign weights to the edges according to their length
			v1 = roadGraph.get_node_attributes(edge[0])
			v2 = roadGraph.get_node_attributes(edge[1])
			length = sqrt((v2[0]-v1[0])**2+(v2[1]-v1[1])**2)
			roadGraph.set_edge_weight(edge[0],edge[1],length)
	
		return roadGraph
	



class LocalNode:

	def __init__(self, probe, contacts, nodeID, rootNode, pixelSize, stabilizePose = False, faceDir = True, travelDir = True):

		self.pixelSize = pixelSize

		self.stabilizePose = stabilizePose
		self.robotParam = probe.robotParam
		self.numSegs = self.robotParam['numSegs']
		self.segLength = self.robotParam['segLength']
		self.mapSize = self.segLength*self.numSegs + 2.0 + 2.0
		
		self.faceDir = faceDir
		self.travelDir = travelDir
		self.partnerNodeID = 0

		self.nodeID = nodeID
		self.probe = probe
		self.contacts = contacts
		self.rootNode = rootNode
		
		
		self.frontProbeError = []
		self.backProbeError = []

		self._refnodes = []
		self._refent = []
		
		self.dirty = False
		
		"""
		initial information for the local node
		
		1. estimated pose
		2. profile of the anchor points 
		"""
		
		self.setEstPose(self.contacts.getAveragePose(self.rootNode))			
		self.setGndPose(self.probe.getActualJointPose(self.rootNode))

		# MAPS
		self.sweepMap = LocalOccMap(self, sweepDir = self.faceDir)
		self.occMap = LocalOccMap(self)
		self.boundaryMap = LocalBoundaryMap(self)
		self.obstacleMap = LocalObstacleMap(self)
		
		self.saveCount = 0
		
		self.localPosture = []
		self.correctedPosture = []
		self.resetPosture()
		
		" compute a fitted curve of a center point "
		self.centerCurve = GPACurve(self.localPosture, rotated=True)
		upoints = arange(0,1.0,0.01)
		self.splPoints = self.centerCurve.getUSet(upoints)

		self.rootPose = [0.0,0.0,0.0]
		self.gndRootPose = [0.0,0.0,0.0]

		if self.stabilizePose:
			self.initialGPACPose = self.getLocalGPACPose()
	
			print "rootPose:", self.rootPose


		" alpha hull boundary "
		self.a_vert = []   # full hull
		self.b_vert = []   # sweep hull
		self.c_vert = []   # static hull
		self.medialPathA = []
		self.medialPathB = []
		self.medialPathC = []
		self.hullAComputed = False
		self.hullBComputed = False		
		self.hullCComputed = False		
		self.medialAComputed = False
		self.medialBComputed = False		
		self.medialCComputed = False		
		self.isBowtie = False
		
		self.medialPathCut = []
		
		self.cornerCandidates = []

		self.longPaths = []
		self.medialLongPaths = []
		self.medialTailCuts = []
		self.longMedialWidths = []
		self.bowtieValues = []
		
		self.hasDeparture = False
		
		self.isFeatureless = None
				
		#pylab.clf()

		#xP = []
		#yP = []
		#for p in self.localPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='b')
		
		
	def setPartnerNodeID(self, nodeID):
		print "setting", self.nodeID, "partner to", nodeID
		self.partnerNodeID = nodeID	
		
	def getIsFeatureless(self):
		
		if self.isFeatureless != None:
			return self.isFeatureless
		
		print "isBowtie:", self.isBowtie
		
		#if self.isBowtie:			
		#	medialAxis = self.getStaticMedialAxis()
		
		#else:
		#	medialAxis = self.getMedialAxis(sweep = False)
		
		hull, medialAxis = computeHullAxis(self.nodeID, self, tailCutOff = False)
		#medialAxis = medialAxis[1:-2]

		print "len(medialAxis) =", len(medialAxis)
		print medialAxis[0], medialAxis[-1]

		medialSpline = SplineFit(medialAxis,smooth=0.1)
		
		#medialPoints = medialSpline.getUniformSamples(spacing = 0.05)
		medialPoints = medialSpline.getUniformSamples(spacing = 0.1)

		print "len(medialPoints) =", len(medialPoints)
		print medialPoints[0], medialPoints[-1]

		xP = []
		yP = []
		for p in medialPoints:
			xP.append(p[0])
			yP.append(p[1])
		
		(ar1,br1)= scipy.polyfit(xP,yP,1)
		ar1 = round(ar1,5)
		br1 = round(br1,5)
		
		print "(ar1,br1) =", (ar1,br1)
		
		" compute point distances from the fitted line "
		xf1 = [medialPoints[0][0], medialPoints[-1][0]]
		yf1 = scipy.polyval([ar1,br1],xf1)			
		linePoints1 = [[xf1[0],yf1[0]], [xf1[-1],yf1[-1]]]
		#linePoints2 = functions.makePointsUniform( linePoints1, max_spacing = 0.04)
		linePoints2 = functions.makePointsUniform( linePoints1, max_spacing = 0.1)
		
		print "linePoints1 =", linePoints1
		print "len(linePoints2) =", len(linePoints2)
		print linePoints2[0], linePoints2[-1]
		
		
		print self.nodeID, "isFeatureless() with", len(medialPoints), "medial points and spline length =", medialSpline.length()
		
		distances = []
		for p in medialPoints:
			p_1, i_1, minDist = functions.findClosestPoint(linePoints2, p)
			distances.append(minDist)
			
		sum = 0.0
		for dist in distances:
			sum += dist
		
		distAvg = sum / len(distances)
		
		print "node", self.nodeID, "has featurelessness of", distAvg
		if distAvg > 0.04:
			self.isFeatureless = False
			return False
		
		else:
			self.isFeatureless = True
			return True
		
		
	def setDeparture(self, val):
		self.hasDeparture = val
		
	def getDeparture(self):
		return self.hasDeparture
		
	def resetPosture(self):
		
		" FIXME:  while reseting the reference posture, if the root node is offset by angle, "
		" the posture will be at an angle and the occ map will have error, two poses overlapped "
		" new posture should angle correct from previous reference posture just like a stabilizeRootPose() call "
		
		#self.localPosture = []
		
		if len(self.localPosture) == 0:
			for j in range(self.numSegs-1):
				self.localPosture.append(self.probe.getJointPose([0.0,0.0,0.0], self.rootNode, j))
		else:
			self.stabilizeRootPose()
			correctedPosture = []
			for j in range(self.numSegs-1):
				correctedPosture.append(self.probe.getJointPose(self.rootPose, self.rootNode, j))

			self.localPosture = correctedPosture
	
		self.correctedPosture = deepcopy(self.localPosture)
	
	def stabilizeRootPose(self):
		
		global estPlotCount 
		
		" compute the new posture "
		" 1. compute the GPAC pose "
		" 2. perform correction of GPAC pose "
		" 3. compute location of rootPose from corrected GPAC pose "



		originPosture = self.localPosture
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

		newPosture = []
		for j in range(self.numSegs-1):
			newPosture.append(self.probe.getJointPose([0.0,0.0,0.0], self.rootNode, j))


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
		
		if foreCost > backCost:
			#plotPosture(originBackPosture, localBackPosture)
			correctedGPACPose = localBackProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		else:
			#plotPosture(originForePosture, localForePosture)					
			correctedGPACPose = localForeProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])

		correctedProfile = Pose(correctedGPACPose)
	
		" offset from 0,0 to newGPACPose "
		" rootPose from neg. offset from correctedGPACPose "
		localRootOffset3 = correctedProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		
		if foreCost > backCost:
			self.rootPose = originBackProfile.convertLocalOffsetToGlobal(localRootOffset3)
		else:
			self.rootPose = originForeProfile.convertLocalOffsetToGlobal(localRootOffset3)


		if isnan(self.rootPose[0]) or isnan(self.rootPose[1]) or isnan(self.rootPose[2]): 
			print "isNaN:  self.rootPose =", self.rootPose
			print "self.localPosture =", self.localPosture
			print "originForePose, originBackPose =", originForePose, originBackPose
			print "originForePosture =", originForePosture
			print "originBackPosture =", originBackPosture
			print "newPosture = ", newPosture
			print "localForePose, localBackPose =", localForePose, localBackPose
			print "localForePosture =", localForePosture 
			print "localBackPosture =", localBackPosture
			print "foreCost, backCost, angle =", foreCost, backCost, angle
			print "correctedGPACPose = ", correctedGPACPose
			print "localRootOffset3 = ", localRootOffset3

		self.correctedPosture = []
		for j in range(self.numSegs-1):
			self.correctedPosture.append(self.probe.getJointPose(self.rootPose, self.rootNode, j))
	
		gndProf = Pose(self.gndPose)
		currGndPose = self.probe.getActualJointPose(self.rootNode)
		self.gndRootPose = gndProf.convertGlobalPoseToLocal(currGndPose)



		#xP = []
		#yP = []
		#for p in self.localPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='b')

		#xP = []
		#yP = []
		#for p in correctedPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='r')    

		#pylab.xlim(-2,2)
		#pylab.ylim(-2,2)
		
		#pylab.savefig("plotPosture%04u.png" % estPlotCount)
		#pylab.clf()		
		
		#estPlotCount += 1
	
		return
		
		






		newPosture = []
		for j in range(self.numSegs-1):
			newPosture.append(self.probe.getJointPose([0.0,0.0,0.0], self.rootNode, j))
					
		newGPAC = GPACurve(newPosture, rotated=True)

		newGPACPose = newGPAC.getPose()
		#print "newGPACPose:", newGPACPose
		
		localGPACPose = self.centerCurve.getPose()
		#print "localGPACPose:", localGPACPose
		
		
		localGPACProfile = Pose(localGPACPose)
		newGPACProfile = Pose(newGPACPose)
				
		gpacPosture1 = []
		gpacPosture2 = []
		for i in range(len(self.localPosture)):
			gpacPosture1.append(localGPACProfile.convertGlobalPoseToLocal(self.localPosture[i]))
			gpacPosture2.append(newGPACProfile.convertGlobalPoseToLocal(newPosture[i]))	

		#angle, cost = estimateMotion.correctOrientation(gpacPosture1, gpacPosture2)
		angle, cost = estimateMotion.correctOrientation2(gpacPosture1, gpacPosture2)

		correctedGPACPose = newGPACProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
		#print "correctedGPACPose:", correctedGPACPose
		correctedProfile = Pose(correctedGPACPose)
	
		" offset from 0,0 to newGPACPose "
		" rootPose from neg. offset from correctedGPACPose "
		localRootOffset1 = localGPACProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		localRootOffset2 = newGPACProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		localRootOffset3 = correctedProfile.convertGlobalPoseToLocal([0.0,0.0,0.0])
		#print "localRootOffset1:", localRootOffset1
		#print "localRootOffset2:", localRootOffset2
		#print "localRootOffset3:", localRootOffset3
		
		#self.rootPose = correctedProfile.convertLocalOffsetToGlobal(localRootOffset)
		#self.rootPose = newGPACProfile.convertLocalOffsetToGlobal(localRootOffset2)
		#print "rootPose:", self.rootPose
		
		
		self.rootPose = localGPACProfile.convertLocalOffsetToGlobal(localRootOffset3)
		
		#return
		
		correctedPosture = []
		for j in range(self.numSegs-1):
			correctedPosture.append(self.probe.getJointPose(self.rootPose, self.rootNode, j))
		
		" plot the initial posture with the stabilized posture "
		#estimateMotion.plotOffset([0.0,0.0,0.0], [0.0,0.0,angle], gpacPosture1, gpacPosture2, cost)
		#estimateMotion.plotOffset([0.0,0.0,0.0], [0.0,0.0,0.0], self.localPosture, correctedPosture, cost)

		#xP = []
		#yP = []
		#for p in self.localPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='b')

		#xP = []
		#yP = []
		#for p in correctedPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='r')
		
		#pylab.scatter([0.0,self.rootPose[0], correctedGPACPose[0], localGPACPose[0]], [0.0,self.rootPose[1], correctedGPACPose[1], localGPACPose[1]])
		
		#pylab.title("Cost = %2.5f, angle = %2.5f" % (cost,angle))

		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)
		
		#pylab.savefig("plotPosture%04u.png" % estPlotCount)
		#pylab.clf()
		#estPlotCount += 1				
		
		#return
		#pose2 = poseProfile1.convertLocalOffsetToGlobal(offset)
		
		poseProfile2 = Pose([0.0,0.0,0.0])
		
		posture2_offset = []
		for p in newPosture:
			nP = poseProfile2.convertLocalOffsetToGlobal(p)
			posture2_offset.append(nP)
		
		xP = []
		yP = []
		for p in posture2_offset:
			#nP = poseProfile2.convertLocalToGlobal(p)
			xP.append(p[0])
			yP.append(p[1])
		
		pylab.plot(xP,yP, color='r')    
		
		rootPose1 = localGPACProfile.convertLocalOffsetToGlobal(localRootOffset1)
		rootPose2 = correctedProfile.convertLocalOffsetToGlobal(localRootOffset2)

		pylab.scatter([0.0, newGPACPose[0], localGPACPose[0]],[0.0, newGPACPose[1], localGPACPose[1]],color='k')
		
		pylab.xlim(-2,2)
		pylab.ylim(-2,2)
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)
		
		pylab.savefig("plotPosture%04u.png" % estPlotCount)
		pylab.clf()		
		
		estPlotCount += 1
	

		xP = []
		yP = []
		for p in gpacPosture1:
			xP.append(p[0])
			yP.append(p[1])
		
		pylab.plot(xP,yP, color='b')
		
		#pose2 = poseProfile1.convertLocalOffsetToGlobal(offset)
		
		poseProfile2 = Pose([0.0,0.0,0.0])
		#poseProfile2 = Pose([0.0,0.0,angle])
		
		posture2_offset = []
		for p in gpacPosture2:
			nP = poseProfile2.convertLocalOffsetToGlobal(p)
			posture2_offset.append(nP)
		
		xP = []
		yP = []
		for p in posture2_offset:
			#nP = poseProfile2.convertLocalToGlobal(p)
			xP.append(p[0])
			yP.append(p[1])
		
		pylab.plot(xP,yP, color='r')    
		
		rootPose1 = localGPACProfile.convertLocalOffsetToGlobal(localRootOffset1)
		rootPose2 = correctedProfile.convertLocalOffsetToGlobal(localRootOffset2)

		print rootPose1, localRootOffset1
		print rootPose2, localRootOffset2

		pylab.scatter([0.0, localRootOffset1[0], localRootOffset2[0]],[0.0, localRootOffset1[1], localRootOffset2[1]],color='k')
		
		#pylab.scatter([0.0, rootPose1[0], rootPose2[0]],[0.0, rootPose1[1], rootPose2[1]],color='k')
		
		pylab.xlim(-2,2)
		pylab.ylim(-2,2)
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)
		
		pylab.savefig("plotPosture%04u.png" % estPlotCount)
		pylab.clf()		
		
		estPlotCount += 1
		
		
		#correctedPosture2 = []
		#for i in range(len(gpacPosture2)):
		#	correctedPosture2.append(localGPACProfile.convertLocalOffsetToGlobal(gpacPosture2[i]))
		
		rootPose3 = localGPACProfile.convertLocalOffsetToGlobal(localRootOffset3)

		correctedPosture2 = []
		for j in range(self.numSegs-1):
			correctedPosture2.append(self.probe.getJointPose(rootPose3, self.rootNode, j))


		xP = []
		yP = []
		for p in self.localPosture:
			xP.append(p[0])
			yP.append(p[1])
		
		pylab.plot(xP,yP, color='b')

		xP = []
		yP = []
		for p in correctedPosture2:
			#nP = poseProfile2.convertLocalToGlobal(p)
			xP.append(p[0])
			yP.append(p[1])
		
		pylab.plot(xP,yP, color='r')    

		pylab.xlim(-2,2)
		pylab.ylim(-2,2)
		#pylab.xlim(-4,4)
		#pylab.ylim(-4,4)
		
		pylab.savefig("plotPosture%04u.png" % estPlotCount)
		pylab.clf()		
		
		estPlotCount += 1
	

	def getStableGPACPosture(self):
		
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		posture_GPAC = []
		for pose in self.correctedPosture:
			posture_GPAC.append(localGPACProfile.convertGlobalPoseToLocal(pose))
		
		return posture_GPAC		

	def getStableGPACurve(self):
		
		posture_GPAC = self.getStableGPACPosture()
		
		GPA_curve = SplineFit(posture_GPAC, smooth = 0.5, kp = 2)

		return GPA_curve
	
	def getGPACPosture(self):

		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		posture_GPAC = []
		for pose in self.localPosture:
			posture_GPAC.append(localGPACProfile.convertGlobalPoseToLocal(pose))
		
		return posture_GPAC
	
	def getGPACurve(self):
		posture_GPAC = self.getGPACPosture()
		
		GPA_curve = SplineFit(posture_GPAC, smooth = 0.5, kp = 2)

		return GPA_curve

	def getLocalGPACPose(self):
		
		return self.centerCurve.getPose()
	
	
	def getGlobalGPACPose(self):

		localGPACPose = self.getLocalGPACPose()	
		#globalGPACPose = self.convertLocalOffsetToGlobal(localGPACPose)
		
		globalEst = [0.0,0.0,0.0]
		
		
		finalVec = array([[localGPACPose[0]], [localGPACPose[1]]])
		transVec = dot(transpose(self.R), finalVec)
		resVec = dot(self.backR, transVec)
		resVec[0, 0] += self.dist
		tempVec = dot(self.foreR, resVec)

		#transVec = dot(transpose(self.gndR), finalVec)
		#resVec = dot(self.gndBackR, transVec)
		#resVec[0, 0] += self.gndDist
		#tempVec = dot(self.gndForeR, resVec)
		
		globalEst[0] = tempVec[0, 0]
		globalEst[1] = tempVec[1, 0]
		
		
		globalEst[2] = self.normalizeAngle(self.estPose[2] + localGPACPose[2])
		#globalEst[2] = self.normalizeAngle(self.gndPose[2] + localGPACPose[2])

		globalGPACPose = globalEst

		return globalGPACPose
		
	def setGPACPose(self, newPose):
		
		gpacProfile = Pose(self.getGlobalGPACPose())
		
		localOffset = gpacProfile.convertGlobalPoseToLocal(self.estPose)
		
		" go back and convert this from GPAC pose to estPose "
		newProfile = Pose(newPose)
		newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
		
		self.setEstPose(newEstPose)

	def getGndGlobalGPACPose(self):

		localGPACPose = self.getLocalGPACPose()	
		#globalGPACPose = self.convertLocalOffsetToGlobal(localGPACPose)
		
		globalEst = [0.0,0.0,0.0]
		
		
		finalVec = array([[localGPACPose[0]], [localGPACPose[1]]])
		transVec = dot(transpose(self.gndR), finalVec)
		resVec = dot(self.gndBackR, transVec)
		resVec[0, 0] += self.gndDist
		tempVec = dot(self.gndForeR, resVec)
		
		globalEst[0] = tempVec[0, 0]
		globalEst[1] = tempVec[1, 0]		
		globalEst[2] = self.normalizeAngle(self.gndPose[2] + localGPACPose[2])

		globalGndGPACPose = globalEst

		return globalGndGPACPose
	
	def getSweepCentroid(self):
		" first get the centroid of the sweepMap "
		
		" 1. pick out the points "
		numPixel = self.sweepMap.numPixel
		mapImage = self.sweepMap.getMap()
		image = mapImage.load()
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					pnt = self.sweepMap.gridToReal([j,k])
					points.append(pnt)
		
		" average the points "
		numPoints = len(points)
		xAvg = 0.0
		yAvg = 0.0
		for p in points:
			xAvg += p[0]
			yAvg += p[1]
			
		xAvg /= numPoints
		yAvg /= numPoints
		

		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		return localGPACProfile.convertGlobalToLocal([xAvg,yAvg])	
	
	def setEstPose(self, newPose):

			
		self.estPose = copy(newPose)
		self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)

		" avoid numerical errors "
		if fabs(self.estPose[0]) > self.dist:
			self.dist = fabs(self.estPose[0])
		elif fabs(self.estPose[1]) > self.dist:
			self.dist = fabs(self.estPose[1])

		
		self.vecAng = acos(self.estPose[0]/self.dist)
		if asin(self.estPose[1]/self.dist) < 0:
			self.vecAng = -self.vecAng
		
		self.backR = array([[cos(self.vecAng), sin(self.vecAng)],[-sin(self.vecAng),cos(self.vecAng)]])
		self.foreR = array([[cos(self.vecAng), -sin(self.vecAng)],[sin(self.vecAng),cos(self.vecAng)]])
		
		self.R = array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])
		
	def getEstPose(self):
		return copy(self.estPose)

	def setGndPose(self, newPose):
		
		self.gndPose = copy(newPose)
		self.gndDist = sqrt(self.gndPose[0]**2 + self.gndPose[1]**2)

		" avoid numerical errors "
		if fabs(self.gndPose[0]) > self.gndDist:
			self.gndDist = fabs(self.gndPose[0])
		elif fabs(self.gndPose[1]) > self.gndDist:
			self.gndDist = fabs(self.gndPose[1])		
		
		
		self.gndVecAng = acos(self.gndPose[0]/self.gndDist)
		if asin(self.gndPose[1]/self.gndDist) < 0:
			self.gndVecAng = -self.gndVecAng
		
		self.gndBackR = array([[cos(self.gndVecAng), sin(self.gndVecAng)],[-sin(self.gndVecAng),cos(self.gndVecAng)]])
		self.gndForeR = array([[cos(self.gndVecAng), -sin(self.gndVecAng)],[sin(self.gndVecAng),cos(self.gndVecAng)]])
		
		self.gndR = array([[cos(self.gndPose[2]), sin(self.gndPose[2])],[-sin(self.gndPose[2]),cos(self.gndPose[2])]])

	def getGndPose(self):
		return copy(self.gndPose)
		
	# this function converts the angle to its equivalent # in the range [-pi,pi]
	def normalizeAngle(self, angle):
	
		while angle>pi:
			angle=angle-2*pi
	
		while angle<=-pi:
			angle=angle+2*pi
	
		return angle 
	
	def convertLocalOffsetToGlobal(self, offset):

		globalEst = [0.0,0.0,0.0]

		finalVec = array([[offset[0]], [offset[1]]])
		transVec = dot(transpose(self.R), finalVec)
		resVec = dot(self.backR, transVec)
		resVec[0, 0] += self.dist
		tempVec = dot(self.foreR, resVec)
		globalEst[0] = tempVec[0, 0]
		globalEst[1] = tempVec[1, 0]
		globalEst[2] = self.normalizeAngle(self.estPose[2] + offset[2])

		return globalEst

	def convertGlobalPoseToLocal(self, pose):

		" transform pnt to local coordinates"
		globalVec = array([[pose[0]],[pose[1]]])

		" perform translation correction "
		tempVec = dot(self.backR, globalVec)
		tempVec[0,0] -= self.dist
		transVec = dot(self.foreR, tempVec)

		" now, apply rotation correction with respect to origin "
		localVec = dot(self.R, transVec)

		localPose = [localVec[0,0], localVec[1,0], self.normalizeAngle(pose[2] - self.estPose[2])]

		return localPose

	def convertLocalToGndGlobal(self, pnt):

		finalVec = array([[pnt[0]], [pnt[1]]])
		transVec = dot(transpose(self.gndR), finalVec)
		resVec = dot(self.gndBackR, transVec)
		resVec[0, 0] += self.gndDist
		tempVec = dot(self.gndForeR, resVec)

		newPoint = [tempVec[0,0],tempVec[1,0]]

		return newPoint

	def convertLocalToGlobal(self, pnt):

		finalVec = array([[pnt[0]], [pnt[1]]])
		transVec = dot(transpose(self.R), finalVec)
		resVec = dot(self.backR, transVec)
		resVec[0, 0] += self.dist
		tempVec = dot(self.foreR, resVec)

		newPoint = [tempVec[0,0],tempVec[1,0]]

		return newPoint

	def convertGlobalToLocal(self, pnt):

		" transform pnt to local coordinates"
		globalVec = array([[pnt[0]],[pnt[1]]])

		" perform translation correction "
		tempVec = dot(self.backR, globalVec)
		tempVec[0,0] -= self.dist
		transVec = dot(self.foreR, tempVec)

		" now, apply rotation correction with respect to origin "
		localVec = dot(self.R, transVec)

		newPoint = [localVec[0,0], localVec[1,0]]
		return newPoint
		
	def getForeTip(self):
		
		segLength = self.probe.segLength
		
		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0		

		joints = range(-1,self.rootNode)
		joints.reverse()

		for i in joints:
			totalAngle = totalAngle + self.probe.getServo(i+1)
			totalAngle = self.normalizeAngle(totalAngle)
			
			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

		return [xTotal, zTotal, totalAngle]

	def getAftTip(self):
		segLength = self.probe.segLength
		
		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0		

		joints = range(self.rootNode, 39)

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)

			if i < 38:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = self.normalizeAngle(totalAngle)
						
		return [xTotal, zTotal, totalAngle]

	def getJointPose(self, jointI):
		
		segLength = self.probe.segLength

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0
	
		if jointI > self.rootNode:
			#joints = range(self.rootNode, self.probe.numSegs-1)
			joints = range(self.rootNode, jointI)
	
			for i in joints:
				xTotal = xTotal + segLength*cos(totalAngle)
				zTotal = zTotal + segLength*sin(totalAngle)
	
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = self.normalizeAngle(totalAngle)
		
		elif jointI < self.rootNode:
			
			joints = range(jointI,self.rootNode)
			joints.reverse()
	
			for i in joints:
				totalAngle = totalAngle + self.probe.getServo(i+1)
				totalAngle = self.normalizeAngle(totalAngle)
				
				xTotal = xTotal - segLength*cos(totalAngle)
				zTotal = zTotal - segLength*sin(totalAngle)
		
		elif jointI == self.rootNode:
			xTotal = 0.0
			zTotal = 0.0

		else:
			print "rootNode =", self.rootNode, "jointI = ", jointI
			raise
	
		return [xTotal, zTotal, totalAngle]
			
	def getProbe(self):
		return self.probe
	
	def getOccMap(self):
		return self.occMap
	
	def getBoundMap(self):
		return self.boundaryMap
	
	def getObstacleMap(self):
		return self.obstacleMap
	
	def getContacts(self):
		return self.contacts
	
	def getNodeID(self):
		return self.nodeID
	
	def getRootNode(self):
		return self.rootNode
	
	#def getRootPose(self):
	#	return self.rootPose
	
	def saveToFile(self):
		"""
		1. save the rectangles
		2. save the pose profile
		3. save the estimated position
		"""
		
		#self.occMap.saveToFile("occ%04u.txt" % self.nodeID)
		#self.occMap.saveMap("occ%04u.png" % self.nodeID)
		#self.poseProfile.saveToFile("prof%04u.txt" % self.nodeID)
		
		f = open("estpose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.estPose))
		f.write("\n")
		f.close()

		f = open("gndpose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.gndPose))
		f.write("\n")
		f.close()

		f = open("posture%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.localPosture))
		f.write("\n")
		f.close()

		f = open("direction%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.faceDir))
		f.write("\n")
		f.close()

		f = open("travelDir%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.travelDir))
		f.write("\n")
		f.close()

		f = open("frontProbeError%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.frontProbeError))
		f.write("\n")
		f.close()

		f = open("backProbeError%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.backProbeError))
		f.write("\n")
		f.close()

		f = open("correctedPosture%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.correctedPosture))
		f.write("\n")
		f.close()

		f = open("rootPose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.rootPose))
		f.write("\n")
		f.close()

		f = open("gndRootPose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.gndRootPose))
		f.write("\n")
		f.close()


	def saveToFile2(self):

		saveFile = ""
		saveFile += "self.localPosture = " + repr(self.localPosture) + "\n"
		saveFile += "self.faceDir = " + repr(self.faceDir) + "\n"
		saveFile += "self.travelDir = " + repr(self.travelDir) + "\n"
		saveFile += "self.frontProbeError = " + repr(self.frontProbeError) + "\n"
		saveFile += "self.backProbeError = " + repr(self.backProbeError) + "\n"
		saveFile += "self.correctedPosture = " + repr(self.correctedPosture) + "\n"
		saveFile += "self.rootPose = " + repr(self.rootPose) + "\n"
		saveFile += "self.gndRootPose = " + repr(self.gndRootPose) + "\n"
		saveFile += "self.partnerNodeID = " + repr(self.partnerNodeID) + "\n"


		saveFile += "self.estPose = " + repr(self.estPose) + "\n"
		saveFile += "self.dist = " + repr(self.dist) + "\n"
		saveFile += "self.vecAng = " + repr(self.vecAng) + "\n"
		saveFile += "self.backR = " + repr(self.backR) + "\n"
		saveFile += "self.foreR = " + repr(self.foreR) + "\n"
		saveFile += "self.R = " + repr(self.R) + "\n"

		saveFile += "self.gndPose = " + repr(self.gndPose) + "\n"
		saveFile += "self.gndDist = " + repr(self.gndDist) + "\n"
		saveFile += "self.gndVecAng = " + repr(self.gndVecAng) + "\n"
		saveFile += "self.gndBackR = " + repr(self.gndBackR) + "\n"
		saveFile += "self.gndForeR = " + repr(self.gndForeR) + "\n"
		saveFile += "self.gndR = " + repr(self.gndR) + "\n"



		f = open("localStateSave_%04u.txt" % (self.nodeID), 'w')
		f.write(saveFile)
		f.close()		


	def readFromFile2(self, dirName, nodeID, forcedPose = []):

		print "loading " + dirName + "/localStateSave_%04u.txt" % nodeID

		self.nodeID = nodeID

		" occupancy map "
		self.occMap.readFromFile(dirName)
	
		" obstacle map "
		self.obstacleMap.readFromFile(dirName)

		" occupancy map "
		self.sweepMap.readFromFile(dirName)


		f = open(dirName + "/localStateSave_%04u.txt" % nodeID, 'r')		
		saveStr = f.read()
		f.close()
		
		saveStr = saveStr.replace('\r\n','\n')		
		exec(saveStr)
		
		if len(forcedPose) != 0:
			self.setEstPose(forcedPose)
			

		" compute a fitted curve of a center point "
		self.centerCurve = GPACurve(self.localPosture, rotated=True)
		upoints = arange(0,1.0,0.01)
		self.splPoints = self.centerCurve.getUSet(upoints)
		
		self.synch()
		self.getBestMedialAxis()
		
	
	def readFromFile(self, dirName, nodeID, forcedPose = []):
		
		self.nodeID = nodeID
		#self.poseProfile.readFromFile("prof%04u.txt" % self.nodeID)
		
		" occupancy map "
		self.occMap.readFromFile(dirName)
	
		" obstacle map "
		self.obstacleMap.readFromFile(dirName)

		" occupancy map "
		self.sweepMap.readFromFile(dirName)

		if len(forcedPose) == 0:
			f = open(dirName + "/estpose%04u.txt" % self.nodeID, 'r')
			estPose = eval(f.read().rstrip())
			f.close()
		else:
			estPose = forcedPose

		#self.setEstPose(estPose)

		f = open(dirName + "/gndpose%04u.txt" % self.nodeID, 'r')
		gndPose = eval(f.read().rstrip())
		f.close()


		self.setEstPose(estPose)
		self.setGndPose(gndPose)

		f = open(dirName + "/posture%04u.txt" % self.nodeID, 'r')
		self.localPosture = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/travelDir%04u.txt" % self.nodeID, 'r')
		self.travelDir = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/direction%04u.txt" % self.nodeID, 'r')
		self.faceDir = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/correctedPosture%04u.txt" % self.nodeID, 'r')
		self.correctedPosture = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/rootPose%04u.txt" % self.nodeID, 'r')
		self.rootPose = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/gndRootPose%04u.txt" % self.nodeID, 'r')
		self.gndRootPose = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/frontProbeError%04u.txt" % self.nodeID, 'r')
		self.frontProbeError = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/backProbeError%04u.txt" % self.nodeID, 'r')
		self.backProbeError = eval(f.read().rstrip())
		f.close()


		#gndProf = Pose(self.gndPose)
		#self.setGndPose(gndProf.convertLocalOffsetToGlobal(self.gndRootPose))
		
		#f = open(dirName + "/posture%04u.txt" % self.nodeID, 'r')
		#self.localPosture = eval(f.read().rstrip())
		#f.close()
		#self.centerCurve = GPACurve(self.localPosture, rotated=True)
		
		" compute a fitted curve of a center point "
		self.centerCurve = GPACurve(self.localPosture, rotated=True)
		upoints = arange(0,1.0,0.01)
		self.splPoints = self.centerCurve.getUSet(upoints)
		
		self.synch()
		
		
		" Plot the GPAC curve and origin "
		"""
		pylab.clf()

		xP = []
		yP = []
		for p in self.localPosture:
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP)	

		
		xP = []
		yP = []
		for p in self.splPoints:
			xP.append(p[0])
			yP.append(p[1])
		
		pylab.plot(xP,yP)


		pose = self.getLocalGPACPose()

		pylab.scatter([pose[0]],[pose[1]])


		pylab.xlim(-3,3)
		pylab.ylim(-2,2)
		pylab.savefig("plotGPAC%04u.png" % self.nodeID)

		"""


		self.getBestMedialAxis()

		
	def obstCallBack(self, direction):
		#self.currNode.obstCallBack(direction)		
		self.obstacleMap.updateContact(direction)

	def updateCorrectPosture(self):
		self.correctedPosture = []
		for j in range(self.numSegs-1):
			self.correctedPosture.append(self.probe.getJointPose(self.rootPose, self.rootNode, j))
		
	def update(self):
		
		if self.stabilizePose:
			self.stabilizeRootPose()
			
		self.occMap.update()
		self.sweepMap.update()
		
		self.dirty = True
		self.hullAComputed = False
		self.hullBComputed = False
		self.hullCComputed = False
		self.medialAComputed = False
		self.medialBComputed = False		
		self.medialCComputed = False		

		self.longPaths = []
		self.medialLongPaths = []
		self.medialTailCuts = []
		self.longMedialWidths = []
		self.bowtieValues = []
	
	def isDirty(self):
		return self.dirty
	
	def synch(self):
		
		self.occMap.buildMap()
		self.sweepMap.buildMap()
		
		#self.computeAlphaBoundary()

		self.boundaryMap.update()
		self.obstacleMap.update()
		
		"""
		self.cornerCandidates = []
		cornerCandidates = extractCornerCandidates(self.sweepMap.getMap(), estPose = self.getGndPose())


		" convert to real coordinates from pixel "
		for cand in cornerCandidates:
			pnt = self.sweepMap.gridToReal(cand[0])
		
			inwardVec = cand[2]
			ang = acos(inwardVec[0])
			if inwardVec[1] < 0.0:
				ang = -ang
		
			localGPACPose = self.getLocalGPACPose()
			localGPACProfile = Pose(localGPACPose)
	
			pnt2 = localGPACProfile.convertGlobalToLocal([pnt[0],pnt[1]])	
			pose2 = localGPACProfile.convertGlobalPoseToLocal([0.0,0.0,ang])


			distFromCenter = sqrt(pnt2[0]**2 + pnt2[1]**2)
			#print "candidate:", cand
			#print "distFromCenter =", distFromCenter, self.nodeID


			#finalCandidates.append((point, cornerAngle, inwardVec))
			
			#self.cornerCandidates.append((pnt2, cand[1], cand[2]))
			if distFromCenter >= 1.35:
				self.cornerCandidates.append((pnt2, cand[1], pose2[2]))
		"""

		self.dirty = False
		
	def saveMap(self):
		
		" FIXME: uncomment these to build the other maps "
		
		" synchronize maps first "
		if self.isDirty():
			self.synch()

		# save a copy of each version of the map
		self.occMap.saveMap()
		self.sweepMap.saveMap()
		self.obstacleMap.saveMap()
		
		self.saveToFile2()

	def getOccPoints(self, sweep = False):



		" 1. pick out the points "
		if sweep:
			numPixel = self.sweepMap.numPixel
			mapImage = self.sweepMap.getMap()
		else:
			numPixel = self.occMap.numPixel
			mapImage = self.occMap.getMap()
		image = mapImage.load()
		
		#print numPixel

		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
				
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					if sweep:
						pnt = self.sweepMap.gridToReal([j,k])
					else:
						pnt = self.occMap.gridToReal([j,k])
						
					points.append(localGPACProfile.convertGlobalToLocal(pnt))
	

		return points

	def getNumLeafs(self):
		
		if len(self.longPaths) > 0:
			return self.numLeafs
		
		self.getBestMedialAxis()
		
		return self.numLeafs

	def computeFullSkeleton(self, alphaHull):
		
		global splineCount

		" DO NOT RECOMPUTE, RETURN PREVIOUS "
		if len(self.longPaths) > 0:
			print "returning caseE"
			return self.longPaths, self.medialLongPaths, self.medialTailCuts, self.longMedialWidths, self.bowtieValues



		
		" FIND BOUNDING BOX "		
		hull = alphaHull
		hull.append(hull[0])
		minX1 = 1e100
		maxX1 = -1e100
		minY1 = 1e100
		maxY1 = -1e100
		for p in hull:
			if p[0] > maxX1:
				maxX1 = p[0]
			if p[0] < minX1:
				minX1 = p[0]
			if p[1] > maxY1:
				maxY1 = p[1]
			if p[1] < minY1:
				minY1 = p[1]
	
	
		" SPECIFY SIZE OF GRID AND CONVERSION PARAMETERS "
		PIXELSIZE = 0.05
		mapSize = 2*max(max(maxX1,math.fabs(minX1)),max(maxY1,math.fabs(minY1))) + 1
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
			point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
			return point


		" CONVERT HULL TO GRID COORDINATES "

		gridHull = []
		for i in range(len(hull)):
			p = hull[i]
			gridHull.append(realToGrid(p))

		" BOUNDING BOX OF GRID HULL "
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

		" FIND INTERIOR POINTS OF GRID HULL "
		#xRange = range(minX, maxX + 1)
		#yRange = range(minY, maxY + 1)	

		#interior = []
			
		#for i in xRange:
		#	for j in yRange:
		#		if point_inside_polygon(i, j, gridHull):
		#			interior.append((i,j))
		
		" POPULATE AN IMAGE WITH GRID HULL AND INTERIOR "
		#inputImg = Image.new('L', (numPixel,numPixel), 0)
		#imga = inputImg.load()
		
		#for i in range(len(gridHull)):
		#	p = gridHull[i]
		#	imga[p[0],p[1]] = 255
		
		#for p in interior:
		#	imga[p[0],p[1]] = 255


		" COMPUTE MEDIAL AXIS OF HULL "
		resultImg = Image.new('L', (numPixel,numPixel))
		resultImg = computeMedialAxis(self.nodeID, numPixel,numPixel, 0, resultImg, len(gridHull[:-2]), gridHull[:-2])
		#resultImg.save("medialOut_%04u.png" % self.nodeID)
		
		
		" EXTRACT POINTS FROM GRID TO LIST "
		imgA = resultImg.load()
	
		points = []
		for i in range(1, numPixel-1):
			for j in range(1, numPixel-1):
				if imgA[i,j] == 0:
					points.append((i,j))
	
	
		" CREATE GRAPH NODES FOR EACH POINT "
		medialGraph = graph.graph()
		for p in points:
			medialGraph.add_node(p, [])
	
		" ADD EDGES BETWEEN NEIGHBORS "
		builtGraph = {}
		for i in range(2, numPixel-2):
			for j in range(2, numPixel-2):
				if imgA[i,j] == 0:
					builtGraph[(i,j)] = []
					for k in range(i-1, i+2):
						for l in range(j-1, j+2):
							if imgA[k,l] == 0 and not (k == i and l == j):
								builtGraph[(i,j)].append((k,l))
								medialGraph.add_edge((i,j), (k,l))
				
		" COMPUTE MINIMUM SPANNING TREE "
		mst = medialGraph.minimal_spanning_tree()
	
		" INITIALIZE DATA DICTS FOR UNIDIRECTIONAL MST"
		uni_mst = {}
		for k, v in mst.items():
			uni_mst[k] = []
	
		" ADD EDGES TO DICT TREE REPRESENTATION "
		for k, v in mst.items():
			if v != None:
				uni_mst[k].append(v)
				uni_mst[v].append(k)
	
		
		" LOCATE ALL LEAVES "
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
		
		" DELETE ALL NODES THAT ARE LEAVES, TO REMOVE SINGLE NODE BRANCHES "
		for leaf in leaves:
			medialGraph.del_node(leaf)	
			
		" RECOMPUTE MST "
		mst = medialGraph.minimal_spanning_tree()
	
		" AGAIN, CREATE OUR DATA STRUCTURES AND IDENTIFY LEAVES "
		uni_mst = {}
		for k, v in mst.items():
			uni_mst[k] = []
	
		for k, v in mst.items():
			if v != None:
				uni_mst[k].append(v)
				uni_mst[v].append(k)
	
	
		" RECORD THE LEAVES AND JUNCTIONS "
		leaves = []
		junctions = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)

			if len(v) > 2:
				junctions.append(k)

		" SAVE FOR LATER, IDENTIFIED BY THEIR INDEX NOW "
		self.numLeafs = len(leaves)
		self.allLeafs = []
		self.allJunctions = []
		for leaf in leaves:
			self.allLeafs.append(gridToReal(leaf))
		
		for junc in junctions:
			self.allJunctions.append(gridToReal(junc))
				
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

		" FOR EVERY PAIR OF LEAVES, SAVE ITS PATH IF ITS LONG ENOUGH "
		" should have X choose 2 combinations"
		
		MAX_LEN = 2
		for leaf1 in leaves:
			for leaf2 in leaves:
				if leaf1 < leaf2:
					if len(nodePaths[leaf1][leaf2]) > MAX_LEN:
						self.longPaths.append((len(nodePaths[leaf1][leaf2]),deepcopy(nodePaths[leaf1][leaf2])))
		
		" SORT FOR THE LONGEST TO SHORTEST "
		self.longPaths.sort(reverse=True)
		
		" REMOVE SIZE FROM TUPLE  "
		for k in range(len(self.longPaths)):
			self.longPaths[k] = self.longPaths[k][1]
			
		
		" GET THE LEAF INDEXES TO EACH PATH "
		self.leafPairs = []
		for path in self.longPaths:
			leaf1 = path[0]
			leaf2 = path[-1]	
			
			leafID1 = leaves.index(leaf1)
			leafID2 = leaves.index(leaf2)

			self.leafPairs.append((leafID1,leafID2))
			
		print "leafPairs:", self.leafPairs

		for k in range(len(self.longPaths)):
			path = self.longPaths[k]
			realPath = []
			for p in path:
				realPath.append(gridToReal(p))
			
			self.longPaths[k] = realPath
					
		if False:
			pylab.clf()
	
			for path in self.longPaths:
				xP = []
				yP = []
				for p in path:
					#p2 = gridToReal(p)
					xP.append(p[0])
					yP.append(p[1])
	
				pylab.plot(xP,yP)
			#pylab.scatter(xP,yP, color='b')
	
			xP = []
			yP = []
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			sizes = []
			for path in self.longPaths:
				sizes.append(len(path))
			
			pylab.title("Paths %s" % sizes)
			pylab.savefig("medialOut2_%04u.png" % self.nodeID)
	

		#longMedialWidths = []
		#medialLongPaths = []
		#medialTailCuts = []
		print len(self.longPaths), "long paths"
		for longPath in self.longPaths:
			print "longPath:", len(longPath)
			
			leafPath = deepcopy(longPath)
			
			frontVec = [0.,0.]
			backVec = [0.,0.]
			indic = range(3)
			#indic = range(16)
			indic.reverse()
		
			for i in indic:
				p1 = leafPath[i+2]
				p2 = leafPath[i]
				vec = [p2[0]-p1[0], p2[1]-p1[1]]
				frontVec[0] += vec[0]
				frontVec[1] += vec[1]
		
				p1 = leafPath[-i-3]
				p2 = leafPath[-i-1]
				vec = [p2[0]-p1[0], p2[1]-p1[1]]
				backVec[0] += vec[0]
				backVec[1] += vec[1]
		
			frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
			backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
			frontVec[0] /= frontMag
			frontVec[1] /= frontMag
			backVec[0] /= backMag
			backVec[1] /= backMag
		
			newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
			newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
			#newP1 = (longPath[0][0] + frontVec[0]*500, longPath[0][1] + frontVec[1]*500)
			#newP2 = (longPath[-1][0] + backVec[0]*500, longPath[-1][1] + backVec[1]*500)

		
			leafPath.insert(0,newP1)
			leafPath.append(newP2)
			
			


			" convert path points to real "
			#realPath = []
			#for p in longPath:
			#	realPath.append(gridToReal(p))

			"""
			pylab.clf()

			xP = []
			yP = []
			for p in longPath:
				p1 = gridToReal(p)
				xP.append(p1[0])
				yP.append(p1[1])
			pylab.plot(xP,yP, color='b')
			
			#pylab.scatter([linePoint[0]], [linePoint[1]])
			
			pylab.xlim(-25,25)
			pylab.ylim(-25,25)		
			
			pylab.title("%d" % self.nodeID)
			pylab.savefig("spline_%06u.png" % splineCount)
			splineCount += 1
			"""
		
		
			" convert path points to real "
			realPath = []
			for p in leafPath:
				realPath.append(p)
				#realPath.append(gridToReal(p))

			"""
			pylab.clf()

			xP = []
			yP = []
			for p in realPath:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='b')
			
			#pylab.scatter([linePoint[0]], [linePoint[1]])
			
			pylab.xlim(-25,25)
			pylab.ylim(-25,25)		
			
			pylab.title("%d" % self.nodeID)
			pylab.savefig("spline_%06u.png" % splineCount)
			splineCount += 1
			"""
			
			posture1 = self.getStableGPACPosture()
			curve1 = StableCurve(posture1)
			uniform1 = curve1.getPlot()
		
			" check if the medial axis is ordered correctly "
			" ORDER NOT GUARANTEED IF MORE THAN ONE AXIS "
			distFore1 = (realPath[1][0]-uniform1[1][0])**2 + (realPath[1][1]-uniform1[1][1])**2
			distBack1 = (realPath[-2][0]-uniform1[-2][0])**2 + (realPath[-2][1]-uniform1[-2][1])**2
		
			if distFore1 > 1 and distBack1 > 1:
				realPath.reverse()
	
	
			medial2 = deepcopy(realPath)
					
			TAILDIST = 0.5
	
			" take the last segments at tips of medial axis"
			edge1 = medial2[0:2]
			edge2 = medial2[-2:]
			
			frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
			backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
			frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
			backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
			
			frontVec[0] /= frontMag
			frontVec[1] /= frontMag
			backVec[0] /= backMag
			backVec[1] /= backMag
			
			" make a longer version of these edges "
			newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
			newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)
	
			edge1 = [newP1, edge1[1]]
			edge2 = [edge2[0], newP2]
		
			" find the intersection points with the hull "
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
			
			" replace the extended edges with a termination point at the hull edge "			
			medial2 = medial2[1:-2]


			
			if isIntersect1:
				medial2.insert(0, point1)
			if isIntersect2:
				medial2.append(point2)
			

			self.medialLongPaths.append(deepcopy(medial2))
			
			" TAIL CUTOFF MEDIAL AXIS "
			medial3 = deepcopy(medial2)
			
			" cut off the tail of the non-sweeping side "
			if self.nodeID % 2 == 0:
				termPoint = medial3[-1]
				for k in range(len(medial3)):
					candPoint = medial3[-k-1]
					dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
					if dist > TAILDIST:
						break
				medial3 = medial3[:-k-1]
	
			else:
				termPoint = medial3[0]
				for k in range(len(medial3)):
					candPoint = medial3[k]
					dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
					if dist > TAILDIST:
						break
				medial3 = medial3[k:]
	
			self.medialTailCuts.append(medial3)
	
			
			print "len(medial2) =", len(medial2)
			medialSpline2 = SplineFit(medial2, smooth=0.1)
	
	
			mPoints = medialSpline2.getUniformSamples()
			tPoints = medialSpline2.getTransformCurve()

			
			if False:
				pylab.clf()
				
				#fig1.clf()
				fig1 = pylab.figure(2)
				#fig1 = pylab.gcf()
				axes2 = fig1.add_subplot(122)
				#fig1.subplot(122)
				xP = []
				yP = []
				for p in mPoints:
					xP.append(p[0])
					yP.append(p[1])
				axes2.plot(xP,yP, color='r')
	
				xP = []
				yP = []
				for p in hull:
					xP.append(p[0])
					yP.append(p[1])
				axes2.plot(xP,yP, color='g')
	
				xP = []
				yP = []
				for p in medial2:
					xP.append(p[0])
					yP.append(p[1])
				axes2.plot(xP,yP, color='b')
	
				#pylab.scatter(xP,yP, color='k')
				
				#pylab.scatter([linePoint[0]], [linePoint[1]])
				
				axes2.set_xlim(-4,4)
				axes2.set_ylim(-4,4)		
				
				axes2.set_title("%d" % self.nodeID)
	
				axes1 = fig1.add_subplot(121)
				axes1.grid(True)
				xP = []
				yP = []
				for p in tPoints:
					xP.append(p[0])
					yP.append(p[1])
				axes1.plot(xP,yP, color='k')
				axes1.set_xlabel("distance")
				axes1.set_ylabel("angle (radians)")
				axes1.set_xlim(0,10)
				axes1.set_ylim(-4,4)		
				
				fig1.set_size_inches(12,6)
				fig1.savefig("spline_%06u.png" % splineCount)
				pylab.clf()			
				pylab.figure(1)
				splineCount += 1
				
			initU = 0.0
			currU = initU
			#termU = 0.9
			termU = 1.0
			
			medialWidth= []
			while True:


								
				linePoint = medialSpline2.getU(currU)
				vec = medialSpline2.getUVector(currU)


										
				rightVec = [vec[0]*cos(pi/2.0) + vec[1]*sin(pi/2.0), -vec[0]*sin(pi/2.0) + vec[1]*cos(pi/2.0)]
				leftVec = [vec[0]*cos(pi/2.0) - vec[1]*sin(pi/2.0), vec[0]*sin(pi/2.0) + vec[1]*cos(pi/2.0)]
	
				edgeR = [linePoint, [linePoint[0]+rightVec[0]*1.0, linePoint[1]+rightVec[1]*1.0]]
				edgeL = [linePoint, [linePoint[0]+leftVec[0]*1.0, linePoint[1]+leftVec[1]*1.0]]
	
				rightPoint = []
				for k in range(len(hull)-1):
					hullEdge = [hull[k],hull[k+1]]
					isIntersect1, point1 = Intersect(edgeR, hullEdge)
					if isIntersect1:
						rightPoint = point1
						break
	
				leftPoint = []
				for k in range(len(hull)-1):
					hullEdge = [hull[k],hull[k+1]]
					isIntersect1, point1 = Intersect(edgeL, hullEdge)
					if isIntersect1:
						leftPoint = point1
						break
	
				if len(rightPoint) > 0:
					distR = sqrt((rightPoint[0]-linePoint[0])**2 + (rightPoint[1]-linePoint[1])**2)
				else:
					distR = 0.0
				
				if len(leftPoint) > 0:
					distL = sqrt((leftPoint[0]-linePoint[0])**2 + (leftPoint[1]-linePoint[1])**2)
				else:
					distL = 0.0
				
				if distL == 0.0:
					distL = distR
					
				if distR == 0.0:
					distR = distL				
				
				medialWidth.append([linePoint[0],linePoint[1],currU, distR, distL])

				if currU >= termU:
					break
	
				nextU = medialSpline2.getUOfDist(currU, 0.04, distIter = 0.001)			
				currU = nextU
					
			self.longMedialWidths.append(medialWidth)
			
			numSamples = len(medialWidth)
			
			" divide into 5 histograms "
			histDiv = numSamples / 5
			currDiv = 0
			divSums = [0.0 for k in range(5)]
			divIndexes = [0 for k in range(5)]
			
			for k in range(5):
				divIndexes[k] = (k+1) * histDiv
			
			" set last bin to infinite boundary "
			divIndexes[-1] = 1e100
			
			" pad the middle with the remainder "
			remainderTotal = numSamples % 5
			divIndexes[2] += remainderTotal
			divIndexes[3] += remainderTotal
				
			print "numSamples:", numSamples
			print "histDiv:", histDiv
			print "divIndexes:", divIndexes
			for k in range(numSamples):
				
				width = medialWidth[k][3] + medialWidth[k][4]
				if k > divIndexes[currDiv]:
					currDiv += 1
				
				divSums[currDiv] += width
				
	
			greatCount = 0
			for k in range(5):
				if divSums[k] > 12.0:
					greatCount += 1
			
			self.bowtieValues.append((greatCount,divSums))


		print "returning caseF"
		return self.longPaths, self.medialLongPaths, self.medialTailCuts, self.longMedialWidths, self.bowtieValues
				

	def getBestMedialAxis(self):
		
		global splineCount
		
		print "getBestMedialAxis(", self.nodeID, "):", self.isBowtie, self.medialAComputed, self.medialCComputed
		
		
		if not self.isBowtie and self.medialAComputed:
			print "returning caseA"
			return deepcopy(self.medialPathA)


		elif self.isBowtie and self.medialCComputed:
			print "returning caseB"
			return deepcopy(self.medialPathC)	

		" make sure alpha boundary is built "
		self.computeAlphaBoundary()
		
		a_data = self.getAlphaBoundary()
		a_data = decimatePoints(a_data)

		" CONVERT HULL TO GPAC COORDINATES "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)

		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))

		" FIND BOUNDING BOX "		
		hull = a_data_GPAC
		hull.append(hull[0])

		longPaths, medialLongPaths, medialTailCuts, longMedialWidths, bowtieValues = self.computeFullSkeleton(hull)

		if len(longPaths) > 1:
			print "MULTI-JUNCTION SITUATION"
			
			" record the edge points, the intersection points, each of splice paths "
			self.splicePaths = deepcopy(longPaths)
			self.splicedMedialPaths = deepcopy(medialLongPaths)
			
			self.numLeafs = len(longPaths)
		
		
		
		" FIRST MEDIAL AXIS IS THE LONGEST "
		realPath = medialLongPaths[0]

		self.medialAComputed = True
		self.medialPathA = realPath
		
		medial2 = medialTailCuts[0]

		
		medialSpline2 = SplineFit(medial2, smooth=0.1)


		if False:
			pylab.clf()
	
			xP = []
			yP = []
			for p in medial2:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP)
	
			xP = []
			yP = []
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
						
			pylab.title("Medial %d" % self.nodeID)
			pylab.savefig("medialOut_%04u.png" % self.nodeID)


		greatCount, divSums = bowtieValues[0]
		print "divSums:", greatCount, divSums
		
		if divSums[2] < divSums[0] and divSums[2] < divSums[4]:
			self.isBowtie = True
			print "BOWTIE"
		elif greatCount >= 4:
			self.isBowtie = True
			print "BOWTIE"
		else:
			self.medialPathCut = medial2
			print "returning caseC"
			return deepcopy(self.medialPathA)


		self.longPaths = []
		self.medialLongPaths = []
		self.medialTailCuts = []
		self.longMedialWidths = []
		self.bowtieValues = []


		print "computing medial axis"

		" make sure alpha boundary is built "
		self.computeStaticAlphaBoundary()
		
		a_data = self.getAlphaBoundary(static = True)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
				
		hull = a_data_GPAC
		hull.append(hull[0])

		longPaths, medialLongPaths, medialTailCuts, longMedialWidths, bowtieValues = self.computeFullSkeleton(hull)

		maxPath = deepcopy(longPaths[0])

		self.medialCComputed = True		
		self.medialPathC = maxPath
		self.medialPathCut = maxPath

		print "returning caseD"
		return deepcopy(self.medialPathC)
		
	
		raise


	def getMedialAxis(self, sweep = False):
		
		global splineCount
		
		if sweep and self.medialBComputed:
			return deepcopy(self.medialPathB)

		if not sweep and self.medialAComputed:
			return deepcopy(self.medialPathA)

		" make sure alpha boundary is built "
		self.computeAlphaBoundary(sweep)
		
		a_data = self.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)

		" CONVERT HULL TO GPAC COORDINATES "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)

		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))

		" FIND BOUNDING BOX "		
		hull = a_data_GPAC
		hull.append(hull[0])

		longPaths, medialLongPaths, medialTailCuts, longMedialWidths, bowtieValues = self.computeFullSkeleton(hull)

		if len(longPaths) > 1:
			print "MULTI-JUNCTION SITUATION"
			
			" record the edge points, the intersection points, each of splice paths "
			self.splicePaths = deepcopy(longPaths)
			self.splicedMedialPaths = deepcopy(medialLongPaths)
			
			self.numLeafs = len(longPaths)
		
		
		
		" FIRST MEDIAL AXIS IS THE LONGEST "
		realPath = medialLongPaths[0]


		if sweep:
			self.medialBComputed = True		
			self.medialPathB = realPath

		if not sweep:
			self.medialAComputed = True
			self.medialPathA = realPath
		
		medial2 = medialTailCuts[0]

		
		medialSpline2 = SplineFit(medial2, smooth=0.1)



		if False:
			pylab.clf()
	
			xP = []
			yP = []
			for p in medial2:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP)
	
			xP = []
			yP = []
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
						
			pylab.title("Medial %d" % self.nodeID)
			pylab.savefig("medialOut_%04u.png" % self.nodeID)


		greatCount, divSums = bowtieValues[0]
		print "divSums:", greatCount, divSums
		
		if divSums[2] < divSums[0] and divSums[2] < divSums[4]:
			self.isBowtie = True
			print "BOWTIE"
		elif greatCount > 4:
			self.isBowtie = True
			print "BOWTIE"
		else:
			self.medialPathCut = medial2


		if sweep:
			return deepcopy(self.medialPathB)

		if not sweep:
			return deepcopy(self.medialPathA)
		
	
		raise


	def getMedialAxis_Old(self, sweep = False):

		global splineCount
		
		if sweep and self.medialBComputed:
			#print "sweep medial axis already computed "
			return deepcopy(self.medialPathB)

		if not sweep and self.medialAComputed:
			#print "non-sweep medial axis already computed "
			return deepcopy(self.medialPathA)

		print "computing medial axis"

		" make sure alpha boundary is built "
		self.computeAlphaBoundary(sweep)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)

		a_data = self.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)

		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))

		occPoints = self.estOccPoints
		occPoints_GPAC = []
		for pnt in occPoints:
			occPoints_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		self.computeAlphaBoundary(sweep = True)
		b_data = self.getAlphaBoundary(sweep = True)
		b_data = decimatePoints(b_data)
				
		b_data_GPAC = []
		for pnt in b_data:
			b_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		


		#print "len(occPoints) =", len(occPoints_GPAC)

		#return a_data_GPAC		
		
		hull = a_data_GPAC
		hull.append(hull[0])

		sweepHull = b_data_GPAC
		sweepHull.append(sweepHull[0])
		
		minX1 = 1e100
		maxX1 = -1e100
		minY1 = 1e100
		maxY1 = -1e100
		for p in hull:
			if p[0] > maxX1:
				maxX1 = p[0]
			if p[0] < minX1:
				minX1 = p[0]
			if p[1] > maxY1:
				maxY1 = p[1]
			if p[1] < minY1:
				minY1 = p[1]
	
	
		PIXELSIZE = 0.05
		mapSize = 2*max(max(maxX1,math.fabs(minX1)),max(maxY1,math.fabs(minY1))) + 1
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
			point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
			return point

		"""
		def computeVoronoi(boundaryPoints):
			
			'''
			Compute the voronoi diagram using a 3rd party library.
			
			Voronoi diagram data structures
			----
			
			self.vertices, 2-tuples
			self.lines, 3-tuples: a,b,c for a*x + b*y = c
			self.edges, 3-tuples, (l, v1, v2)  line, vertex 1, vertex 2
			vertex -1 means edge extends to infinity
			'''
	
			sites = []
			for p in boundaryPoints:
				sites.append(Site(p[0],p[1]))

			vertices = []
			lines = []
			edges = []			
			if sites != []:
				vertices, lines, edges = computeVoronoiDiagram(sites)

			return vertices, lines, edges
	
		def pruneEdges(occMap, vertices, lines, edges):
			
			'''remove edges and vertices that are outside the explored area'''
	
			#NEIGHBOR_WIDTH = 3
			NEIGHBOR_WIDTH = 0

			roadGraph = graph.graph()	
			
			xSize, ySize = numPixel, numPixel
			#occPix = self.occMapImage.load()
	
			#occPix = self.occMap.image
			occPix = occMap.load()
	
			newVertices = {}
			newEdges = []
	
			# checking all edges to check within the free space
			for edge in edges:
				v1 = edge[1]
				v2 = edge[2]
	
				if v1 != -1 and v2 != -1:
	
					# get grid value of real point
					i1, j1 = realToGrid(vertices[v1])
					i2, j2 = realToGrid(vertices[v2])
					
					if i1 < xSize and i1 >= 0 and i2 < xSize and i2 >= 0 and j1 < ySize and j1 >= 0 and j2 < ySize and j2 >= 0 :
						seg = [(i1,j1),(i2,j2)]
	
						indX1 = seg[0][0]
						indY1 = seg[0][1]
						indX2 = seg[1][0]
						indY2 = seg[1][1]
	
						isOK = True
	
						for i in range(-NEIGHBOR_WIDTH,NEIGHBOR_WIDTH+1,1):
							for j in range(-NEIGHBOR_WIDTH,NEIGHBOR_WIDTH+1,1):
	
								if indX1+i > numPixel or indX1+i < 0 or indY1+j > numPixel or indY1+j < 0:
									isOK = False
								else:
									val1 = occPix[indX1 + i, indY1 + j]
									if val1 <= 127: # 0 or 127
										isOK = False
	
								if indX2+i > numPixel or indX2+i < 0 or indY2+j > numPixel or indY2+j < 0:
									isOK = False
								else:
									val2 = occPix[indX2 + i, indY2 + j]
									if val2 <= 127: # 0 or 127
										isOK = False
	
						if isOK:
							# 1. add node if not already added
							# 2. add edge
	
							newVertices[v1] = vertices[v1]
							newVertices[v2] = vertices[v2]
	
							newEdges.append([v1,v2])
	
			for key, attr in newVertices.iteritems():
				roadGraph.add_node(key, attr)
	
			for edge in newEdges:
				roadGraph.add_edge(edge[0], edge[1])
	
				# assign weights to the edges according to their length
				v1 = roadGraph.get_node_attributes(edge[0])
				v2 = roadGraph.get_node_attributes(edge[1])
				length = sqrt((v2[0]-v1[0])**2+(v2[1]-v1[1])**2)
				roadGraph.set_edge_weight(edge[0],edge[1],length)
		
			return roadGraph
		
		"""

		gridOccPoints = []
		for pnt in occPoints_GPAC:
			gridOccPoints.append(realToGrid(pnt))
			#gridOccPoints.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		gridHull = []
		for i in range(len(hull)):
			p = hull[i]
			gridHull.append(realToGrid(p))

		gridSweep = []
		for i in range(len(sweepHull)):
			p = sweepHull[i]
			gridSweep.append(realToGrid(p))
	
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

		xRange = range(minX, maxX + 1)
		yRange = range(minY, maxY + 1)	

		interior = []
		interior2 = []
	
		for i in xRange:
			for j in yRange:
				if point_inside_polygon(i, j, gridHull):
					interior.append((i,j))

		for i in xRange:
			for j in yRange:
				if point_inside_polygon(i, j, gridSweep):
					interior2.append((i,j))
		
		
		
		inputImg = Image.new('L', (numPixel,numPixel), 0)
		imga = inputImg.load()
		
		for i in range(len(gridHull)):
			p = gridHull[i]
			imga[p[0],p[1]] = 255
		
		for p in interior:
			imga[p[0],p[1]] = 255
			
		
		inputImg2 = Image.new('L', (numPixel,numPixel), 0)
		imga = inputImg2.load()
		
		for i in range(len(gridSweep)):
			p = gridSweep[i]
			imga[p[0],p[1]] = 255
		
		for p in interior2:
			imga[p[0],p[1]] = 255

		
		#inputImg.save("medialOut_%04u_1.png" % self.nodeID)

		
		"""
		vertices, lines, edges = computeVoronoi(hull[:-2])
		roadGraph = pruneEdges(inputImg, vertices, lines, edges)
		voronoiImg = Image.new('L', (numPixel,numPixel), 255)
		
		
		voronoiDraw = ImageDraw.Draw(voronoiImg)

		edges = roadGraph.edges()
		for edge in edges:
			v1 = roadGraph.get_node_attributes(edge[0])
			v2 = roadGraph.get_node_attributes(edge[1])
			
			# get grid value of real point
			i1, j1 = realToGrid(v1)
			i2, j2 = realToGrid(v2)

			voronoiDraw.line([(i1,j1),(i2,j2)], fill = 0)


		vertices, lines, edges = computeVoronoi(sweepHull[:-2])
		roadGraph = pruneEdges(inputImg2, vertices, lines, edges)
		voronoiImg2 = Image.new('L', (numPixel,numPixel), 255)
				
		voronoiDraw = ImageDraw.Draw(voronoiImg2)

		edges = roadGraph.edges()
		for edge in edges:
			v1 = roadGraph.get_node_attributes(edge[0])
			v2 = roadGraph.get_node_attributes(edge[1])
			
			# get grid value of real point
			i1, j1 = realToGrid(v1)
			i2, j2 = realToGrid(v2)

			voronoiDraw.line([(i1,j1),(i2,j2)], fill = 0)
		"""


	
		resultImg = Image.new('L', (numPixel,numPixel))
		resultImg2 = Image.new('L', (numPixel,numPixel))
		resultImg = computeMedialAxis(self.nodeID, numPixel,numPixel, 0, resultImg, len(gridHull[:-2]), gridHull[:-2])
		resultImg2 = computeMedialAxis(self.nodeID, numPixel,numPixel, 0, resultImg2, len(gridSweep[:-2]), gridSweep[:-2])


		occImg = Image.new('L', (numPixel,numPixel), 0)
		imga = occImg.load()
		for i in range(len(gridOccPoints)):
			p = gridOccPoints[i]
			imga[p[0],p[1]] = 255
			

		"""
		imgOut = Image.new('L', (4*numPixel,2*numPixel))
		imgOut.paste(occImg, (0,0))
		imgOut.paste(inputImg, (numPixel,0))
		imgOut.paste(resultImg, (2*numPixel, 0))
		imgOut.paste(voronoiImg, (3*numPixel, 0))
		imgOut.paste(inputImg2, (numPixel, numPixel))
		imgOut.paste(resultImg2, (2*numPixel, numPixel))
		imgOut.paste(voronoiImg2, (3*numPixel, numPixel))
		imgOut.save("medialOut_%04u.png" % self.nodeID)
		"""
		
		#resultImg.save("medialOut_%04u.png" % self.nodeID)
		imgA = resultImg.load()
	
		points = []
		for i in range(1, numPixel-1):
			for j in range(1, numPixel-1):
				if imgA[i,j] == 0:
					points.append((i,j))
	
		medialGraph = graph.graph()
		for p in points:
			medialGraph.add_node(p, [])
	
		builtGraph = {}
		for i in range(2, numPixel-2):
			for j in range(2, numPixel-2):
				if imgA[i,j] == 0:
					builtGraph[(i,j)] = []
					for k in range(i-1, i+2):
						for l in range(j-1, j+2):
							if imgA[k,l] == 0 and not (k == i and l == j):
								builtGraph[(i,j)].append((k,l))
								medialGraph.add_edge((i,j), (k,l))
								
		mst = medialGraph.minimal_spanning_tree()
	
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
	
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
				
				
		#print len(leaves), "leaves"

		for leaf in leaves:
			#print "delete", leaf
			medialGraph.del_node(leaf)	
		mst = medialGraph.minimal_spanning_tree()
	
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
	
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
				
				
		print len(leaves), "leaves"

		if True:
			imgA = resultImg2.load()
			points2 = []
			for i in range(1, numPixel-1):
				for j in range(1, numPixel-1):
					if imgA[i,j] == 0:
						points2.append((i,j))
		
			medialGraph2 = graph.graph()
			for p in points2:
				medialGraph2.add_node(p, [])
		
			builtGraph2 = {}
			for i in range(2, numPixel-2):
				for j in range(2, numPixel-2):
					if imgA[i,j] == 0:
						builtGraph2[(i,j)] = []
						for k in range(i-1, i+2):
							for l in range(j-1, j+2):
								if imgA[k,l] == 0 and not (k == i and l == j):
									builtGraph2[(i,j)].append((k,l))
									medialGraph2.add_edge((i,j), (k,l))
									
			mst2 = medialGraph2.minimal_spanning_tree()
		
			uni_mst2 = {}
			isVisited2 = {}
			nodeSum2 = {}
			for k, v in mst2.items():
				uni_mst2[k] = []
				isVisited2[k] = 0
				nodeSum2[k] = 0
		
			for k, v in mst2.items():
				if v != None:
					uni_mst2[k].append(v)
					uni_mst2[v].append(k)
		
			leaves2 = []
			for k, v in uni_mst2.items():
				if len(v) == 1:
					leaves2.append(k)
					
					
			#print len(leaves), "leaves"
	
			for leaf in leaves2:
				#print "delete", leaf
				medialGraph2.del_node(leaf)	
			mst2 = medialGraph2.minimal_spanning_tree()
		
			uni_mst2 = {}
			isVisited2 = {}
			nodeSum2 = {}
			for k, v in mst2.items():
				uni_mst2[k] = []
				isVisited2[k] = 0
				nodeSum2[k] = 0
		
			for k, v in mst2.items():
				if v != None:
					uni_mst2[k].append(v)
					uni_mst2[v].append(k)
		
			leaves2 = []
			for k, v in uni_mst2.items():
				if len(v) == 1:
					leaves2.append(k)
					
					
			#print len(leaves2), "leaves"
			

	
		" find the longest path from each leaf"
		maxPairDist = 0
		maxPair = None
		maxPath = []
		nodePaths = {}
		for leaf in leaves:
	
			isVisited = {}
			nodeSum = {}
			nodePath = {}
			for k, v in uni_mst.items():
				isVisited[k] = 0
				nodeSum[k] = 0
				nodePath[k] = []
	
			getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
	
			nodePaths[leaf] = nodePath
	
			maxDist = 0
			maxNode = None
			for k, v in nodeSum.items():
				#print k, v
				if v > maxDist:
					maxNode = k
					maxDist = v
	
			#print leaf, "-", maxNode, maxDist
	
			if maxDist > maxPairDist:
				maxPairDist = maxDist
				maxPair = [leaf, maxNode]
				maxPath = nodePath[maxNode]

		"""
		self.longPaths = []

		MAX_LEN = 2
		for leaf1 in leaves:
			for leaf2 in leaves:
				if leaf1 < leaf2:
					if len(nodePaths[leaf1][leaf2]) > MAX_LEN:
						self.longPaths.append((len(nodePaths[leaf1][leaf2]),deepcopy(nodePaths[leaf1][leaf2])))
					
		self.longPaths.sort(reverse=True)
		
		" remove the size from tuple "
		for k in range(len(self.longPaths)):
			self.longPaths[k] = self.longPaths[k][1]
		"""
					
		if False:
			pylab.clf()
	
			for path in self.longPaths:
				xP = []
				yP = []
				for p in path:
					p2 = gridToReal(p)
					xP.append(p2[0])
					yP.append(p2[1])
	
				pylab.plot(xP,yP)
			#pylab.scatter(xP,yP, color='b')
	
			xP = []
			yP = []
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color='r')
			
			sizes = []
			for path in self.longPaths:
				sizes.append(len(path))
			
			pylab.title("Paths %s" % sizes)
			pylab.savefig("medialOut_%04u.png" % self.nodeID)
	
	


	
	
		
		frontVec = [0.,0.]
		backVec = [0.,0.]
		#indic = range(6)
		indic = range(16)
		indic.reverse()
	
		for i in indic:
			p1 = maxPath[i+2]
			p2 = maxPath[i]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			frontVec[0] += vec[0]
			frontVec[1] += vec[1]
	
			p1 = maxPath[-i-3]
			p2 = maxPath[-i-1]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			backVec[0] += vec[0]
			backVec[1] += vec[1]
	
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
	
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
	
		newP1 = (maxPath[0][0] + frontVec[0]*500, maxPath[0][1] + frontVec[1]*500)
		newP2 = (maxPath[-1][0] + backVec[0]*500, maxPath[-1][1] + backVec[1]*500)
	
		maxPath.insert(0,newP1)
		maxPath.append(newP2)
	
	
		" convert path points to real "
		realPath = []
		for p in maxPath:
			realPath.append(gridToReal(p))

		posture1 = self.getStableGPACPosture()
		curve1 = StableCurve(posture1)
		uniform1 = curve1.getPlot()
	
		" check if the medial axis is ordered correctly "
		distFore1 = (realPath[1][0]-uniform1[1][0])**2 + (realPath[1][1]-uniform1[1][1])**2
		distBack1 = (realPath[-2][0]-uniform1[-2][0])**2 + (realPath[-2][1]-uniform1[-2][1])**2
	
		if distFore1 > 1 and distBack1 > 1:
			realPath.reverse()



		if sweep:
			self.medialBComputed = True		
			self.medialPathB = realPath

		if not sweep:
			self.medialAComputed = True
			self.medialPathA = realPath


		if sweep:
			medial2 = self.medialPathB
		else:
			medial2 = self.medialPathA
		
		TAILDIST = 0.5

		" take the long length segments at tips of medial axis"
		edge1 = medial2[0:2]
		edge2 = medial2[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]
	
		" find the intersection points with the hull "
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
		
		" replace the extended edges with a termination point at the hull edge "			
		medial2 = medial2[1:-2]
		if isIntersect1:
			medial2.insert(0, point1)
		if isIntersect2:
			medial2.append(point2)
		
		" cut off the tail of the non-sweeping side "
		if self.nodeID % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]

		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]


		
		medialSpline2 = SplineFit(medial2, smooth=0.1)


		initU = 0.0
		currU = initU
		#termU = 0.9
		termU = 1.0
		self.medialWidth = []
		self.lineEdges = []
		while True:
			
			linePoint = medialSpline2.getU(currU)
			vec = medialSpline2.getUVector(currU)
			
			angle = acos(vec[0])
			if asin(vec[1]) < 0:
				angle = -angle


			rightAng = angle + pi/2.0
			leftAng = angle - pi/2.0
			
			rightVec = [vec[0]*cos(pi/2.0) + vec[1]*sin(pi/2.0), -vec[0]*sin(pi/2.0) + vec[1]*cos(pi/2.0)]
			leftVec = [vec[0]*cos(pi/2.0) - vec[1]*sin(pi/2.0), vec[0]*sin(pi/2.0) + vec[1]*cos(pi/2.0)]

			edgeR = [linePoint, [linePoint[0]+rightVec[0]*1.0, linePoint[1]+rightVec[1]*1.0]]
			edgeL = [linePoint, [linePoint[0]+leftVec[0]*1.0, linePoint[1]+leftVec[1]*1.0]]

			rightPoint = []
			for k in range(len(hull)-1):
				hullEdge = [hull[k],hull[k+1]]
				isIntersect1, point1 = Intersect(edgeR, hullEdge)
				if isIntersect1:
					rightPoint = point1
					break

			leftPoint = []
			for k in range(len(hull)-1):
				hullEdge = [hull[k],hull[k+1]]
				isIntersect1, point1 = Intersect(edgeL, hullEdge)
				if isIntersect1:
					leftPoint = point1
					break

			#self.lineEdges.append(edgeR)
			#self.lineEdges.append(edgeL)

			#print leftPoint, rightPoint, linePoint
			if len(rightPoint) > 0:
				distR = sqrt((rightPoint[0]-linePoint[0])**2 + (rightPoint[1]-linePoint[1])**2)
				self.lineEdges.append([rightPoint,linePoint])
			else:
				distR = 0.0
			
			if len(leftPoint) > 0:
				distL = sqrt((leftPoint[0]-linePoint[0])**2 + (leftPoint[1]-linePoint[1])**2)
				self.lineEdges.append([leftPoint,linePoint])
			else:
				distL = 0.0
			
			if distL == 0.0:
				distL = distR
				
			if distR == 0.0:
				distR = distL				
			
			#print "distR,distL =", distR, distL
			#print "currU =", currU, "avgDist =", (distR+distL)/2.0, "diff =", distR-distL
			self.medialWidth.append([linePoint[0],linePoint[1],currU, distR, distL])

			#print "currU =", currU

			if currU >= termU:
				break

			nextU = medialSpline2.getUOfDist(currU, 0.04, distIter = 0.001)			
			currU = nextU
				
		
		numSamples = len(self.medialWidth)
		
		" divide into 5 histograms "
		histDiv = numSamples / 5
		currDiv = 0
		divSums = [0.0 for k in range(5)]
		divIndexes = [0 for k in range(5)]
		
		for k in range(5):
			divIndexes[k] = (k+1) * histDiv
		
		" set last bin to infinite boundary "
		divIndexes[-1] = 1e100
		
		" pad the middle with the remainder "
		remainderTotal = numSamples % 5
		divIndexes[2] += remainderTotal
		divIndexes[3] += remainderTotal
			
		print "numSamples:", numSamples
		print "histDiv:", histDiv
		print "divIndexes:", divIndexes
		for k in range(numSamples):
			
			width = self.medialWidth[k][3] + self.medialWidth[k][4]
			#print k, ">", divIndexes[currDiv]
			if k > divIndexes[currDiv]:
				currDiv += 1
			
			divSums[currDiv] += width
			

		greatCount = 0
		for k in range(5):
			if divSums[k] > 12.0:
				greatCount += 1

		print "divSums:", greatCount, divSums
			
		if divSums[2] < divSums[0] and divSums[2] < divSums[4]:
			self.isBowtie = True
			print "BOWTIE"
		elif greatCount > 4:
			self.isBowtie = True
			print "BOWTIE"
		else:
			self.medialPathCut = medial2




		if sweep:
			return deepcopy(self.medialPathB)

		if not sweep:
			return deepcopy(self.medialPathA)
		
	
		raise


	def getStaticMedialAxis(self):

		if self.medialCComputed:
			return deepcopy(self.medialPathC)

		print "computing medial axis"

		" make sure alpha boundary is built "
		self.computeStaticAlphaBoundary()
		
		a_data = self.getAlphaBoundary(static = True)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
				
		hull = a_data_GPAC
		hull.append(hull[0])

		longPaths, medialLongPaths, medialTailCuts, longMedialWidths, bowtieValues = self.computeFullSkeleton(hull)

		maxPath = deepcopy(longPaths[0])

		self.medialCComputed = True		
		self.medialPathC = maxPath

		return deepcopy(self.medialPathC)

	def getStaticMedialAxis_Old(self):

		if self.medialCComputed:
			return deepcopy(self.medialPathC)

		print "computing medial axis"

		" make sure alpha boundary is built "
		self.computeStaticAlphaBoundary()

		a_data = self.getAlphaBoundary(static = True)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
				
		hull = a_data_GPAC
		hull.append(hull[0])

		
		minX = 1e100
		maxX = -1e100
		minY = 1e100
		maxY = -1e100
		for p in hull:
			if p[0] > maxX:
				maxX = p[0]
			if p[0] < minX:
				minX = p[0]
			if p[1] > maxY:
				maxY = p[1]
			if p[1] < minY:
				minY = p[1]
	
	
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
			point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
			return point
	
		gridHull = []
		for i in range(len(hull)):
			p = hull[i]
			gridHull.append(realToGrid(p))
	
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

		resultImg = Image.new('L', (numPixel,numPixel))
		resultImg = computeMedialAxis(self.nodeID, numPixel,numPixel, 0, resultImg, len(gridHull[:-2]), gridHull[:-2])
		
		
		#resultImg.save("medialOut_%04u.png" % self.nodeID)
		imgA = resultImg.load()
	
		points = []
		for i in range(1, numPixel-1):
			for j in range(1, numPixel-1):
				if imgA[i,j] == 0:
					points.append((i,j))
	
		medialGraph = graph.graph()
		for p in points:
			medialGraph.add_node(p, [])
	
		builtGraph = {}
		for i in range(2, numPixel-2):
			for j in range(2, numPixel-2):
				if imgA[i,j] == 0:
					builtGraph[(i,j)] = []
					for k in range(i-1, i+2):
						for l in range(j-1, j+2):
							if imgA[k,l] == 0:
								builtGraph[(i,j)].append((k,l))
								medialGraph.add_edge((i,j), (k,l))
								
		mst = medialGraph.minimal_spanning_tree()
	
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
	
		leaves = []
		for k, v in uni_mst.items():
			if len(v) == 1:
				leaves.append(k)
	
	
		" find the longest path from each leaf"
		maxPairDist = 0
		maxPair = None
		maxPath = []
		for leaf in leaves:
	
			isVisited = {}
			nodeSum = {}
			nodePath = {}
			for k, v in uni_mst.items():
				isVisited[k] = 0
				nodeSum[k] = 0
				nodePath[k] = []
	
			getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
	
			maxDist = 0
			maxNode = None
			for k, v in nodeSum.items():
				#print k, v
				if v > maxDist:
					maxNode = k
					maxDist = v
	
			#print leaf, "-", maxNode, maxDist
	
			if maxDist > maxPairDist:
				maxPairDist = maxDist
				maxPair = [leaf, maxNode]
				maxPath = nodePath[maxNode]
	
	
		frontVec = [0.,0.]
		backVec = [0.,0.]
		#indic = range(6)
		indic = range(16)
		indic.reverse()
	
		for i in indic:
			p1 = maxPath[i+2]
			p2 = maxPath[i]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			frontVec[0] += vec[0]
			frontVec[1] += vec[1]
	
			p1 = maxPath[-i-3]
			p2 = maxPath[-i-1]
			vec = [p2[0]-p1[0], p2[1]-p1[1]]
			backVec[0] += vec[0]
			backVec[1] += vec[1]
	
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
	
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
	
		newP1 = (maxPath[0][0] + frontVec[0]*500, maxPath[0][1] + frontVec[1]*500)
		newP2 = (maxPath[-1][0] + backVec[0]*500, maxPath[-1][1] + backVec[1]*500)
	
		maxPath.insert(0,newP1)
		maxPath.append(newP2)
	
	
		" convert path points to real "
		realPath = []
		for p in maxPath:
			realPath.append(gridToReal(p))

		posture1 = self.getStableGPACPosture()
		curve1 = StableCurve(posture1)
		uniform1 = curve1.getPlot()
	
		" check if the medial axis is ordered correctly "
		distFore1 = (realPath[1][0]-uniform1[1][0])**2 + (realPath[1][1]-uniform1[1][1])**2
		distBack1 = (realPath[-2][0]-uniform1[-2][0])**2 + (realPath[-2][1]-uniform1[-2][1])**2
	
		if distFore1 > 1 and distBack1 > 1:
			realPath.reverse()



		self.medialCComputed = True		
		self.medialPathC = realPath

		medial2 = self.medialPathC

		
		TAILDIST = 0.5

		" take the long length segments at tips of medial axis"
		edge1 = medial2[0:2]
		edge2 = medial2[-2:]
		
		frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
		backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
		frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
		backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
		
		frontVec[0] /= frontMag
		frontVec[1] /= frontMag
		backVec[0] /= backMag
		backVec[1] /= backMag
		
		" make a smaller version of these edges "
		newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
		newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

		edge1 = [newP1, edge1[1]]
		edge2 = [edge2[0], newP2]
	
		" find the intersection points with the hull "
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
		
		" replace the extended edges with a termination point at the hull edge "			
		medial2 = medial2[1:-2]
		if isIntersect1:
			medial2.insert(0, point1)
		if isIntersect2:
			medial2.append(point2)
		
		" cut off the tail of the non-sweeping side "
		if self.nodeID % 2 == 0:
			termPoint = medial2[-1]
			for k in range(len(medial2)):
				candPoint = medial2[-k-1]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[:-k-1]

		else:
			termPoint = medial2[0]
			for k in range(len(medial2)):
				candPoint = medial2[k]
				dist = sqrt((termPoint[0]-candPoint[0])**2 + (termPoint[1]-candPoint[1])**2)
				if dist > TAILDIST:
					break
			medial2 = medial2[k:]
			
		self.medialPathCut = medial2
			


		return deepcopy(self.medialPathC)
	
		raise


	def computeStaticAlphaBoundary(self):

		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		if self.hullCComputed:
			return
		
		points = []

		segLength = self.probe.segLength
		segWidth = self.probe.segWidth

		for pose in self.correctedPosture:
			
			xTotal = pose[0]
			zTotal = pose[1]
			totalAngle = pose[2]
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			points.append(p4)
			points.append(p3)
			points.append(p2)
			points.append(p1)
		
		if len(points) > 0:
				
			self.c_vert = self.computeAlpha2(points, radius = 0.12)
			" cut out the repeat vertex "
			self.c_vert = self.c_vert[:-1]
			
			self.c_vert = self.convertAlphaUniform(self.c_vert)
		
		else:
			self.c_vert = []
			
		self.hullCComputed = True
					
	def computeAlphaBoundary(self, sweep = False):


		" synchronize maps first "
		if self.isDirty():
			self.synch()

		if sweep and self.hullBComputed:
			#print "sweep hull already computed "
			return

		if not sweep and self.hullAComputed:
			#print "non-sweep hull already computed "
			return
		
		" 1. pick out the points "
		if sweep:
			numPixel = self.sweepMap.numPixel
			mapImage = self.sweepMap.getMap()
		else:
			numPixel = self.occMap.numPixel
			mapImage = self.occMap.getMap()
		image = mapImage.load()
		
		#print numPixel
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					if sweep:
						pnt = self.sweepMap.gridToReal([j,k])
					else:
						pnt = self.occMap.gridToReal([j,k])
					points.append(pnt)
		
		self.estOccPoints = points
		#print len(points)
		
		if len(points) > 0:
			#print points
		
			if sweep:	
				self.b_vert = self.computeAlpha2(points, radius = 0.001)
				" cut out the repeat vertex "
				self.b_vert = self.b_vert[:-1]
				
				self.b_vert = self.convertAlphaUniform(self.b_vert)

			else:
				self.a_vert = self.computeAlpha2(points, radius = 0.2)
				" cut out the repeat vertex "
				self.a_vert = self.a_vert[:-1]
				
				self.a_vert = self.convertAlphaUniform(self.a_vert)
		
		else:
			if sweep:
				self.b_vert = []
			else:
				self.a_vert = []
		
		if sweep:		
			self.hullBComputed = True
		else:
			self.hullAComputed = True			
			
	def getAlphaBoundary(self, sweep = False, static = False):
		if static:
			return self.c_vert
		if sweep:
			return self.b_vert
		else:
			return self.a_vert
			

	def getBareHull(self):
		
		if self.isBowtie:
			self.computeStaticAlphaBoundary()
	
			a_data = self.getAlphaBoundary(static=True)
			a_data = decimatePoints(a_data)
	
			" convert hull points to GPAC coordinates before adding covariances "
			localGPACPose = self.getLocalGPACPose()
			localGPACProfile = Pose(localGPACPose)
			
			a_data_GPAC = []
			for pnt in a_data:
				a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
			
		else:
					
			" Read in data of Alpha-Shapes without their associated covariances "
			self.computeAlphaBoundary()
			a_data = self.getAlphaBoundary()
			a_data = decimatePoints(a_data)
			
			" convert hull points to GPAC coordinates "
			localGPACPose = self.getLocalGPACPose()
			localGPACProfile = Pose(localGPACPose)
			
			a_data_GPAC = []
			for pnt in a_data:
				a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		return a_data_GPAC
	
	def convertAlphaUniform(self, a_vert, max_spacing = 0.04):
		
		" make the vertices uniformly distributed "
		
		new_vert = []
		
		#max_spacing = 0.04
		#max_spacing = 0.1
		
		for i in range(len(a_vert)):
			p0 = a_vert[i]
			p1 = a_vert[(i+1) % len(a_vert)]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			vec[0] /= dist
			vec[1] /= dist
			
			new_vert.append(copy(p0))
			
			if dist > max_spacing:
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					new_vert.append(newP)
		
		return new_vert			
		#self.a_vert = new_vert
	
	def computeAlpha2(self, points, radius = 0.2):
		
		global alphaPlotCount
		plotIter = True
		
		random.seed(0)		
		
		isDone = False
		
		while not isDone:
	
			perturbPoints = []
			
			for p in points:
				p2 = copy(p)
				" add a little bit of noise to avoid degenerate conditions in CGAL "
				p2[0] += random.gauss(0.0,0.000001)
				p2[1] += random.gauss(0.0,0.000001)
	
				perturbPoints.append(p2)
		
			try:			

				saveFile = ""	
				saveFile += "radius = " + repr(radius) + "\n"
				saveFile += "perturbPoints = " + repr(perturbPoints) + "\n"
				
				isWritten = False
				while not isWritten:
					try:
						f = open("doLocalAlphaInput_%08u.txt" % (alphaPlotCount), 'w')
						f.write(saveFile)
						f.close()
					except:
						pass
					else:
						isWritten = True
	
				vertices = alphamod.doAlpha(radius,perturbPoints)
				numVert = len(vertices)
	
				os.remove("doLocalAlphaInput_%08u.txt" % (alphaPlotCount))
	
				if plotIter:
                    pylab.clf()
					xP = []
					yP = []
					for p in vertices:
						xP.append(p[0])
						yP.append(p[1])
						
					pylab.plot(xP,yP, color='r')
	
					pylab.xlim(-3,3)
					pylab.ylim(-3,3)
					pylab.title("nodeID = %d, radius = %f, numPoints = %d" % (self.nodeID, radius, numVert))
					pylab.savefig("alphaResult_%04d.png" % alphaPlotCount)
	
				alphaPlotCount += 1
									
				if numVert <= 2:
					print "Failed, hull had only", numVert, "vertices"
					raise
				
				isDone = True
			except:
				print "hull has holes!  retrying..."
				#print sArr	

		return vertices
			
		plotIter = False
		
		numPoints = len(points)

		isDone = False
		
		while not isDone:

			inputStr = str(numPoints) + " "
	
			" alpha shape circle radius "
			inputStr += str(radius) + " "
			
			xP = []
			yP = []
			for p in points:
				p2 = copy(p)
				" add a little bit of noise to avoid degenerate conditions in CGAL "
				p2[0] += random.gauss(0.0,0.0001)
				p2[1] += random.gauss(0.0,0.0001)
	
				xP.append(p2[0])
				yP.append(p2[1])
	
				inputStr += str(p2[0]) + " " + str(p2[1]) + " "
			
			inputStr += "\n"
			
			if plotIter:
				pylab.clf()
				pylab.scatter(xP,yP, color='b')
			try:			
				" start the subprocess "
				if sys.platform == "win32":
					subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
				else:
					subProc = Popen(["./alpha2"], stdin=PIPE, stdout=PIPE)
					
				
				" send input and receive output "
				sout, serr = subProc.communicate(inputStr)
		
				#print numPoints
				#print sout
				
				" convert string output to typed data "
				sArr = sout.split(" ")
				
		
				#print sArr[0]
				numVert = int(sArr[0])
				
				sArr = sArr[1:]
				
				
				vertices = []
				for i in range(len(sArr)/2):
					vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
					
				if plotIter:
                    pylab.clf()
					xP = []
					yP = []
					for p in vertices:
						xP.append(p[0])
						yP.append(p[1])
						
					pylab.plot(xP,yP, color='r')
	
					pylab.xlim(-3,3)
					pylab.ylim(-3,3)
					pylab.title("nodeID = %d, radius = %f, numPoints = %d" % (self.nodeID, radius, numPoints))
					pylab.savefig("alphaResult_%04d.png" % alphaCount)
					alphaCount += 1
								
				isDone = True
			except:
				if plotIter:
                    pylab.clf()
					pylab.xlim(-3,3)
					pylab.ylim(-3,3)
					pylab.title("FAIL: nodeID = %d, radius = %f, numPoints = %d" % (self.nodeID, radius, numPoints))
					pylab.savefig("alphaResult_%04d.png" % alphaCount)
					alphaCount += 1
				print "hull has holes!  retrying..."
				#print sArr
		"""
		else:
			maxX = 0
			minX = 1e10
			maxY = 0
			minY = 1e10
			
			for p in points:
				if p[0] > maxX:
					maxX = p[0]
				
				if p[0] < minX:
					minX = p[0]
				
				if p[1] > maxY:
					maxY = p[1]
					
				if p[1] < minY:
					minY = p[1]
			
			vertices = []
			vertices.append([maxX,maxY])
			vertices.append([maxX,minY])
			vertices.append([minX,minY])
			vertices.append([minX,maxY])
		"""
			
		return vertices
							
