import sys

from LocalOccMap import *
from LocalBoundaryMap import *
from LocalObstacleMap import *
from SplineFit import *
from Pose import Pose
from cornerDetection import extractCornerCandidates

import estimateMotion
from GPACurve import GPACurve
from StableCurve import StableCurve
#from gen_icp import computeMedialAxis, getLongestPath
from medialaxis import computeMedialAxis

import graph

import random
import functions
#import Image
from numpy import array, dot, transpose

estPlotCount = 0

def decimatePoints(points):
	result = []

	for i in range(len(points)):
		#if i%2 == 0:
		if i%4 == 0:
			result.append(points[i])

	return result

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
		self.a_vert = []
		self.b_vert = []
		self.c_vert = []
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
		
		
		self.hasDeparture = False
		
		self.isFeatureless = None
		
		#pylab.clf()

		#xP = []
		#yP = []
		#for p in self.localPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='b')
		
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
		self.medialAComputed = False
		self.medialBComputed = False		
	
	def isDirty(self):
		return self.dirty
	
	def synch(self):
		
		self.occMap.buildMap()
		self.sweepMap.buildMap()
		
		#self.computeAlphaBoundary()

		self.boundaryMap.update()
		self.obstacleMap.update()
		
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
		
		self.saveToFile()

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

	def getMedialAxis(self, sweep = False):

		if sweep and self.medialBComputed:
			#print "sweep medial axis already computed "
			return deepcopy(self.medialPathB)

		if not sweep and self.medialAComputed:
			#print "non-sweep medial axis already computed "
			return deepcopy(self.medialPathA)

		print "computing medial axis"

		" make sure alpha boundary is built "
		self.computeAlphaBoundary(sweep)
		
		def decimatePoints(points):
			result = []
		
			for i in range(len(points)):
				#if i%2 == 0:
				if i%4 == 0:
					result.append(points[i])
		
			return result

		a_data = self.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = self.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
		#return a_data_GPAC		
		
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

		"""		
		xRange = range(minX, maxX + 1)
		yRange = range(minY, maxY + 1)	

		interior = []
	
		for i in xRange:
			for j in yRange:
				if point_inside_polygon(i, j, gridHull):
					interior.append((i,j))
		
		inputImg = Image.new('L', (numPixel,numPixel), 255)
		imga = inputImg.load()
		
		for i in range(len(gridHull)):
			p = gridHull[i]
			imga[p[0],p[1]] = 0
		
		for p in interior:
			imga[p[0],p[1]] = 0
		"""
		
		#inputImg.save("medialIn_%04u.png" % self.nodeID)
	
		#inputImg = Image.new('L', (numPixel,numPixel), 255)
		resultImg = Image.new('L', (numPixel,numPixel))
		#resultImg = Image.new('L', (numPixel,numPixel))
		#resultImg = computeMedialAxis(inputImg, resultImg)
		#resultImg = computeMedialAxis(self.nodeID, numPixel,numPixel, inputImg, resultImg, len(gridHull[:-2]), gridHull[:-2])
		resultImg = computeMedialAxis(self.nodeID, numPixel,numPixel, resultImg, len(gridHull[:-2]), gridHull[:-2])
		
		
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
			#print sqrt(rightVec[0]*rightVec[0] + rightVec[1]*rightVec[1]), sqrt(leftVec[0]*leftVec[0] + leftVec[1]*leftVec[1])
			#rightVec = [cos(rightAng), sin(rightAng)]
			#leftVec = [cos(leftAng), sin(leftAng)]

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

			nextU = medialSpline2.getUOfDist(currU, 0.04)			
			currU = nextU
				
		

		totalWidth = 0.0
		for samp in self.medialWidth:
			width = samp[3] + samp[4]
			totalWidth += width

		#if totalWidth > 65.0:
		#if totalWidth > 60.0:
		if totalWidth > 50.0:
			self.isBowtie = True
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
		
		def decimatePoints(points):
			result = []
		
			for i in range(len(points)):
				#if i%2 == 0:
				if i%4 == 0:
					result.append(points[i])
		
			return result

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
		resultImg = computeMedialAxis(self.nodeID, numPixel,numPixel, resultImg, len(gridHull[:-2]), gridHull[:-2])
		
		
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
				self.b_vert = self.computeAlpha2(points, radius = 0.02)
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
		
		numPoints = len(points)

		isDone = False
		
		while not isDone:

			inputStr = str(numPoints) + " "
	
			" alpha shape circle radius "
			inputStr += str(radius) + " "
			
			for p in points:
				p2 = copy(p)
				" add a little bit of noise to avoid degenerate conditions in CGAL "
				p2[0] += random.gauss(0.0,0.0001)
				p2[1] += random.gauss(0.0,0.0001)
	
				inputStr += str(p2[0]) + " " + str(p2[1]) + " "
			
			inputStr += "\n"
			
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
				isDone = True
			except:
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
							