import sys
import os

from LocalOccMap import *
from LocalBoundaryMap import *
from LocalObstacleMap import *
from SplineFit import *
from Pose import Pose
#from cornerDetection import extractCornerCandidates
#from voronoi import Site, computeVoronoiDiagram

from GPACurve import GPACurve
from StableCurve import StableCurve
#from gen_icp import computeMedialAxis, getLongestPath
from medialaxis import computeMedialAxis
import matplotlib.pyplot as plt

import graph
import alphamod

import random
import functions
from functions import decimatePoints
#from PIL import Image
from numpy import array, dot, transpose

estPlotCount = 0
alphaPlotCount = 0
splineCount = 0
THRESH_WIDTH = 0.45

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

	def __init__(self, probe, contacts, nodeID, rootNode, pixelSize, stabilizePose = False, faceDir = True, travelDir = True, useBloom = True, useBend = True):

		self.pixelSize = pixelSize

		self.stabilizePose = stabilizePose
		self.robotParam = probe.robotParam
		self.numSegs = self.robotParam['numSegs']
		self.segLength = self.robotParam['segLength']
		self.mapSize = self.segLength*self.numSegs + 2.0 + 2.0
		
		self.faceDir = faceDir
		self.travelDir = travelDir
		self.partnerNodeID = 0

		self.useBloom = useBloom
		self.useBend = useBend

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
		self.spatialFeatures = []
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
			
		sumVal = 0.0
		for dist in distances:
			sumVal += dist
		
		distAvg = sumVal / len(distances)
		
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

	def convertGPACtoRawPose(self, newPose):

		gpacProfile = Pose(self.getGlobalGPACPose())
		
		localOffset = gpacProfile.convertGlobalPoseToLocal(self.estPose)
		
		" go back and convert this from GPAC pose to estPose "
		newProfile = Pose(newPose)
		newEstPose = newProfile.convertLocalOffsetToGlobal(localOffset)
		
		return newEstPose

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

		#self.nodeID = nodeID

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
		self.spatialFeatures = []
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

	@logFunction
	def computeFullSkeleton(self, alphaHull):
		
		global splineCount

		print "computeFullSkeleton(", self.nodeID, ")"

		" DO NOT RECOMPUTE, RETURN PREVIOUS "
		if len(self.longPaths) > 0:
			print "returning caseE"
			return self.longPaths, self.medialLongPaths, self.medialTailCuts, self.longMedialWidths, self.bowtieValues

			self.spatialFeatures.append(results)


		
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
					
		if True:
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
		for pathIndex in range(len(self.longPaths)):
			longPath = self.longPaths[pathIndex]
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

			"""
			Are we comparing the correct coordinate systems here?
			realPath is in localGPAC because input alphaHull is in localGPAC
			posture1 looks to be in localGPAC
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


			
			if self.faceDir and isIntersect1:
				medial2.insert(0, point1)
			if not self.faceDir and isIntersect2:
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
			#medialSpline2 = SplineFit(medial2, smooth=0.05)
	
	
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

				curveAngle = acos(vec[0])
				if asin(vec[1]) < 0:
					curveAngle = -curveAngle
										
				linePoint.append(curveAngle)

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
				
				widthEntry = {}
				widthEntry["linePoint"] = copy(linePoint)
				widthEntry["currU"] = currU
				widthEntry["distR"] = distR
				widthEntry["distL"] = distL
				widthEntry["widthSum"] = distR+distL
				widthEntry["asymm"] = distR-distL
				widthEntry["rightPoint"] = copy(rightPoint)
				widthEntry["leftPoint"] = copy(leftPoint)
				widthEntry["distU"] = medialSpline2.dist_u(currU)

				medialWidth.append(widthEntry)
				#medialWidth.append([copy(linePoint),currU, distR, distL, copy(rightPoint), copy(leftPoint)])

				if currU >= termU:
					break
	
				nextU = medialSpline2.getUOfDist(currU, 0.04, distIter = 0.001)			
				currU = nextU

	
			#DERIV_WIDTH = 6
			DERIV_WIDTH = 16


			for k in range(len(medialWidth)):
				
				if (k-DERIV_WIDTH) > 0 and (k+DERIV_WIDTH) < len(medialWidth):

					diffs = []

					for l in range(0,DERIV_WIDTH+1):

						if (k-DERIV_WIDTH+l) >= 0 and (k+l) < len(medialWidth):
							ang1 = medialWidth[k-DERIV_WIDTH+l]["linePoint"][2]
							ang2 = medialWidth[k+l]["linePoint"][2]
							diffs.append(fabs(diffAngle(normalizeAngle(ang1),normalizeAngle(ang2))))


					angDeriv = sum(diffs) / float(len(diffs))

					medialWidth[k]["angDeriv"] = angDeriv
				else:
					medialWidth[k]["angDeriv"] = 0.0


				if k > 1:

					diffs = []
					ang1 = medialWidth[k-2]["linePoint"][2]
					ang2 = medialWidth[k]["linePoint"][2]

					distU1 = medialWidth[k-2]["distU"]
					distU2 = medialWidth[k]["distU"]

					angDeriv = fabs(diffAngle(normalizeAngle(ang1),normalizeAngle(ang2))/(distU2-distU1))

					medialWidth[k]["angDeriv2"] = angDeriv
				else:
					medialWidth[k]["angDeriv2"] = 0.0


			self.longMedialWidths.append(medialWidth)


			"""
			1) bloom feature detector ( contigCount, near edge )
			2) arch feature detector
			3) bending feature detector
			4) bowtie feature detector

			{linePoint, currU, distR, distL, rightPoint, leftPoint, angDeriv, widthSum}
			"""



			if True:

				longPathWidth = medialWidth

				widthX = []
				widthY = []

				contigCount = 0
				contigPoints = []
				contigWeights = []

				contigAsymmSum = []
				contigBoundaries = []
				contigPointIndex = []
				contigCenterMass = []
				contigLen = []
				contigArea = []
				contigCounts = []
				contigDensity = []
				contigAsymm = []
				contigNegArea = []
				contigSlopes = []
				contigInflection = 0
				maxDeriv = 0.0

				contigStartIndex = 0
				contigFinalIndex = 0


				TURN_WIDTH = 15

				#for val in longPathWidth:
				for kIndex in range(len(longPathWidth)):

					val = longPathWidth[kIndex]
					distR = val["distR"]
					distL = val["distL"]
					widthSum = val["widthSum"]
					angDeriv = val["angDeriv2"]
					distU = val["distU"]
					asymm = val["asymm"]

					if angDeriv > maxDeriv and kIndex > TURN_WIDTH and kIndex < (len(longPathWidth)-TURN_WIDTH):
						maxDeriv = angDeriv
						contigInflection = kIndex

					if widthSum > THRESH_WIDTH:

						if contigCount == 0:
							contigStartIndex = kIndex

						contigCount += 1
						contigPoints.append(distU)
						contigWeights.append(widthSum)
						contigAsymm.append(asymm)

					else:

						if contigCount > 0:

							contigFinalIndex = kIndex-1

							weightSum = sum(contigWeights)
							centerMass = 0.0
							threshArea = 0.0
							maxWidth = 0.0
							maxIndex = None
							totalArea = 0.0

							for k in range(len(contigPoints)):
								centerMass += contigPoints[k] * contigWeights[k]
								threshArea += contigWeights[k] - THRESH_WIDTH
								totalArea += contigWeights[k]

								if contigWeights[k] > maxWidth:
									maxWidth = contigWeights[k]
									maxIndex = k


							maxArea = (maxWidth-THRESH_WIDTH) * len(contigPoints)
							#contigNegArea.append((maxArea-threshArea)/maxArea)
							#contigNegArea.append(maxArea-threshArea)
							#contigNegArea.append((maxArea-threshArea)/float(len(contigPoints)))
							#contigNegArea.append(maxWidth - THRESH_WIDTH)
							contigNegArea.append((maxArea-threshArea)/(maxWidth*float(len(contigPoints))))

							" search for local valleys "
							minIndex1 = contigStartIndex
							minWidth1 = longPathWidth[minIndex1]["widthSum"]

							minIndex2 = contigFinalIndex
							minWidth2 = longPathWidth[minIndex2]["widthSum"]

							while True:

								nextIndex1 = minIndex1 - 1
								nextIndex2 = minIndex2 + 1

								if nextIndex1 < 0:
									break

								if nextIndex2 >= len(longPathWidth):
									break

								nextWidth1 = longPathWidth[nextIndex1]["widthSum"]
								nextWidth2 = longPathWidth[nextIndex2]["widthSum"]

								minReached = False

								if nextWidth1 < minWidth1:
									minWidth1 = nextWidth1
									minIndex1 = nextIndex1
								else:
									minReached = True

								if nextWidth2 < minWidth2:
									minWidth2 = nextWidth2
									minIndex2 = nextIndex2
								else:
									minReached = True

								if minReached:
									break

							maxMinWidth = minWidth1
							if minWidth2 > minWidth1:
								maxMinWidth = minWidth2

							x1 = longPathWidth[minIndex1]["distU"]
							x2 = longPathWidth[contigStartIndex+maxIndex]["distU"]
							x3 = longPathWidth[minIndex2]["distU"]

							p1 = [x1, maxMinWidth]
							p2 = [x2, longPathWidth[contigStartIndex+maxIndex]["widthSum"]]
							p3 = [x3, maxMinWidth]
							#p1 = [contigStartIndex, longPathWidth[contigStartIndex]["widthSum"]]
							#p2 = [contigStartIndex + maxIndex, contigWeights[maxIndex]]
							#p3 = [contigFinalIndex, longPathWidth[contigFinalIndex]["widthSum"]]

							print "indices", contigStartIndex, maxIndex, contigFinalIndex, p1[0], p2[0], p3[0]

							if p2[0]-p1[0] != 0.0:
								slope1 = (p2[1]-p1[1])/(p2[0]-p1[0])
							else:
								slope1 = 1e100

							if p2[0]-p3[0] != 0.0:
								slope2 = (p2[1]-p3[1])/(p2[0]-p3[0])
							else:
								slope2 = 1e100

							#contigSlopes.append((slope1,slope2))
	

							vec1 = [p1[0]-p2[0],p1[1]-p2[1]]
							mag1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1])
							if mag1 > 0.0:
								vec1 = [vec1[0]/mag1, vec1[1]/mag1]
								angle1 = acos(vec1[0])
								if asin(vec1[1]) < 0:
									angle1 = -angle1

								#dAng1 = -angle1
								dAng1 = diffAngle(angle1,pi)
							else:
								dAng1 = 0.0


							vec2 = [p3[0]-p2[0],p3[1]-p2[1]]
							mag2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1])
							if mag2 > 0.0:
								vec2 = [vec2[0]/mag2, vec2[1]/mag2]
								angle2 = acos(vec2[0])
								if asin(vec2[1]) < 0:
									angle2 = -angle2

								dAng2 = -angle2
								#dAng2 = diffAngle(angle2,pi)
							else:
								dAng2 = 0.0

							contigSlopes.append((dAng1,dAng2,p1,p2,p3))
							



							centerMass = centerMass / weightSum
							contigCenterMass.append(centerMass)

							contigArea.append(threshArea)
							contigCounts.append(contigCount)
							contigDensity.append(threshArea/float(contigCount))
							contigAsymmSum.append(sum(contigAsymm))

							for k in range(len(longPathWidth)-1):
								distU1 = longPathWidth[k]["distU"]
								distU2 = longPathWidth[k+1]["distU"]

								if distU1 <= centerMass and distU2 > centerMass:
									if fabs(distU1-centerMass) > fabs(distU2-centerMass):
										contigPointIndex.append(k+1)
									else:
										contigPointIndex.append(k)


							contigLen.append(len(contigPoints))
							contigBoundaries.append((contigStartIndex, contigFinalIndex))

							contigCount = 0
							contigPoints = []
							contigWeights = []
							contigStartIndex = 0
							contigFinalIndex = 0
							contigAsymm = []


				#contigBoundaries.append((contigStartIndex, contigFinalIndex))

				print "spatial", self.nodeID, contigInflection, maxDeriv, len(longPathWidth), contigPointIndex, contigBoundaries
				print "mass deriv:", [(k, contigLen[k], contigArea[k], contigDensity[k], longPathWidth[contigPointIndex[k]]["angDeriv2"]) for k in range(len(contigPointIndex))]
				print "contigArea:", self.nodeID, contigArea
				print "contigLen:", self.nodeID, contigLen
				print "contigDensity:", self.nodeID, contigDensity

				#if maxDeriv >= 0.35:
				#if maxDeriv >= 0.6 and longPathWidth[contigInflection-1]["angDeriv2"] < maxDeriv and longPathWidth[contigInflection+1]["angDeriv2"] < maxDeriv:
				#if maxDeriv >= 0.6:
				#if False:
				if self.useBend and maxDeriv >= 0.6:

					inflectionPoint = None

					for k in range(len(contigBoundaries)):
						lowIndex, highIndex = contigBoundaries[k]
						if contigInflection >= lowIndex and contigInflection <= highIndex:

							centerMassIndex = contigPointIndex[k]

							angDeriv = longPathWidth[centerMassIndex]["angDeriv2"]
							#if angDeriv >= maxDeriv - 0.2:
							if angDeriv >= 0.6:
								inflectionPoint = copy(longPathWidth[centerMassIndex]["linePoint"])
								#print self.nodeID, "inflection boundaries:", inflectionPoint, contigInflection, angDeriv, maxDeriv, len(longPathWidth)
								print self.nodeID, "inflection boundaries:", inflectionPoint, centerMassIndex, angDeriv, maxDeriv, len(longPathWidth)

					if inflectionPoint == None:
						if longPathWidth[contigInflection-1]["angDeriv2"] < maxDeriv and longPathWidth[contigInflection+1]["angDeriv2"] < maxDeriv:
							angDeriv = longPathWidth[contigInflection]["angDeriv2"]
							if angDeriv >= maxDeriv - 0.2:
								inflectionPoint = copy(longPathWidth[contigInflection]["linePoint"])
								print self.nodeID, "inflection boundaries:", inflectionPoint, contigInflection, angDeriv, maxDeriv, len(longPathWidth)


				else:
					inflectionPoint = None

				" conditions for a bloom detection "
				bloomPoint = None
				if self.useBloom and len(contigPointIndex) == 1:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					frontBoundIndex, backBoundIndex = contigBoundaries[0]
					pointIndex = contigPointIndex[0]
					angDeriv = longPathWidth[pointIndex]["angDeriv"]
					asymm = contigAsymmSum[0]
					negArea = contigNegArea[0]
					slopes = contigSlopes[0]

					print self.nodeID, "bloom boundaries:", frontBoundIndex, backBoundIndex, len(longPathWidth), cLen, dens, angDeriv, asymm, area, negArea, slopes[0]+slopes[1]

					print frontBoundIndex <= 6, backBoundIndex >= len(longPathWidth)-7, cLen < 25, dens >= 0.05

					if (frontBoundIndex <= 6 or backBoundIndex >= len(longPathWidth)-7) and cLen < 25 and dens >= 0.05:
						bloomPoint = longPathWidth[pointIndex]["linePoint"]

						#ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')

						print "bloom ACCEPT"
					else:
						print "bloom REJECT"


				" conditions for a arch detection "
				archPoint = None
				if self.useBloom and len(contigPointIndex) == 1 and bloomPoint == None:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					frontBoundIndex, backBoundIndex = contigBoundaries[0]
					pointIndex = contigPointIndex[0]
					angDeriv = longPathWidth[pointIndex]["angDeriv"]
					asymm = contigAsymmSum[0]
					negArea = contigNegArea[0]
					slopes = contigSlopes[0]

					print self.nodeID, "arch boundaries:", frontBoundIndex, backBoundIndex, len(longPathWidth), cLen, dens, angDeriv, asymm, area, negArea, slopes[0]+slopes[1]
					#print frontBoundIndex > 6, backBoundIndex < len(longPathWidth)-7, angDeriv < 0.03, cLen > 5, cLen <= 20, area > 0.1, (slopes[0]+slopes[1]) > 0.37
					print frontBoundIndex > 6, backBoundIndex < len(longPathWidth)-7, angDeriv < 0.1, cLen > 5, cLen <= 20, area > 0.1, (slopes[0]+slopes[1]) > 0.37


					#if (frontBoundIndex > 6 and backBoundIndex < len(longPathWidth)-7) and angDeriv < 0.03 and cLen > 5 and cLen <= 20 and area > 0.2 and dens >= 0.03:
					#if (frontBoundIndex > 6 and backBoundIndex < len(longPathWidth)-7) and angDeriv < 0.03 and cLen > 5 and cLen <= 20 and area > 0.2 and (slopes[0]+slopes[1]) > 0.37:
					if (frontBoundIndex > 6 and backBoundIndex < len(longPathWidth)-7) and angDeriv < 0.1 and cLen > 5 and cLen <= 20 and area > 0.1 and (slopes[0]+slopes[1]) > 0.37:
						# and angDeriv < 0.1:
						#pass and cLen < 25 and dens >= 0.03:
						archPoint = longPathWidth[pointIndex]["linePoint"]

						print "arch ACCEPT"
					else:
						print "arch REJECT"


				for k in range(len(contigLen)):

					if contigLen[k] > 40 and contigDensity[k] > 0.1:
						print self.nodeID, "BOWTIE from spatial features"
						self.isBowtie = True

					#print "contigLen:", self.nodeID, contigLen
					#print "contigDensity:", self.nodeID, contigDensity



				results = {}
				results["inflectionPoint"] = inflectionPoint
				results["maxDeriv"] = maxDeriv
				results["maxDerivU"] = longPathWidth[contigInflection]["distU"]
				results["centerMasses"] = contigCenterMass
				results["contigCounts"] = contigCounts
				results["contigDensities"] = contigDensity
				results["contigAreas"] = contigArea
				results["contigLens"] = contigLen
				results["contigPointIndex"] = contigPointIndex
				results["contigBoundaries"] = contigBoundaries
				results["contigAsymmSum"] = contigAsymmSum
				results["contigNegArea"] = contigNegArea
				results["contigSlopes"] = contigSlopes
				results["medialWidths"] = longPathWidth
				results["bloomPoint"] = bloomPoint
				results["archPoint"] = archPoint


				self.spatialFeatures.append(results)


			widthStr = ""
			widthStr += "medialWidth: %d %d %d " % (self.nodeID, pathIndex, len(longPath))
			for width in medialWidth:
				widthSum = width["widthSum"]
				#widthStr += "%0.2f " % widthSum
				if widthSum > 0.4:
					widthStr += "1 "
				else:
					widthStr += "0 "

			print widthStr
					
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
				
				width = medialWidth[k]["widthSum"]
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


		if True:

			fig, (ax1, ax2, ax3) = plt.subplots(3,  sharex=False, sharey=False)

			ax1.set_xlim(-3, 3)
			ax1.set_aspect("equal")

			ax2.set_ylim(0.0, 1.0)
			ax2.set_aspect("equal")

			ax3.set_ylim(-pi, pi)
			#ax3.set_aspect("equal")


			xP = []
			yP = []
			#for p in medial2:
			for p in realPath:
				xP.append(p[0])
				yP.append(p[1])

			ax1.plot(xP,yP, color='b')
	
			xP = []
			yP = []
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			ax1.plot(xP,yP, color='r')

			if True:

				longPathWidth = longMedialWidths[0]

				results = self.spatialFeatures[0]

				maxDeriv = results["maxDeriv"]
				contigCenterMass = results["centerMasses"]
				contigCounts = results["contigCounts"]
				contigDensity = results["contigDensities"]
				contigArea = results["contigAreas"]
				contigLen = results["contigLens"]
				contigPointIndex = results["contigPointIndex"]
				contigBoundaries = results["contigBoundaries"]
				longPathWidth = results["medialWidths"]
				contigAsymmSum = results["contigAsymmSum"]
				contigNegArea = results["contigNegArea"]
				contigSlopes = results["contigSlopes"]

				numPoints = len(longPathWidth)

				for kIndex in range(len(contigCenterMass)):

					area = contigArea[kIndex]
					cLen = contigLen[kIndex]
					frontBoundIndex, backBoundIndex = contigBoundaries[kIndex]

				pointX = []
				pointY = []
				widthX = []
				widthY = []
				angDerivX = []
				angDerivY = []
				angDeriv2X = []
				angDeriv2Y = []
				angX = []
				angY = []

				#for val in longPathWidth:
				for kIndex in range(len(longPathWidth)):

					val = longPathWidth[kIndex]
					distR = val["distR"]
					distL = val["distL"]
					linePoint = val["linePoint"]
					rightPoint = val["rightPoint"]
					leftPoint = val["leftPoint"]
					widthSum = val["widthSum"]
					angDeriv = val["angDeriv"]
					angDeriv2 = val["angDeriv2"]
					distU = val["distU"]

					pointX.append(linePoint[0])
					pointY.append(linePoint[1])


					widthX.append(distU)
					widthY.append(widthSum)
					angDerivX.append(distU)
					angDerivY.append(angDeriv)
					angDeriv2X.append(distU)
					angDeriv2Y.append(angDeriv2)
					angX.append(distU)
					angY.append(linePoint[2])

					#if len(widthX) > 0:
					#	widthX.append(widthX[-1] + 0.04)
					#else:
					#	widthX.append(0.0)

					#print "distR:", distR
					#print "distL:", distL
					#print "linePoint:", linePoint
					#print "leftPoint:", leftPoint
					#print "rightPoint:", rightPoint

					if len(rightPoint) > 0.0:
						xP = [linePoint[0], rightPoint[0]]
						yP = [linePoint[1], rightPoint[1]]
						if widthSum > THRESH_WIDTH:
							ax1.plot(xP,yP, color='k')
						else:
							ax1.plot(xP,yP, color='b')

					if len(leftPoint) > 0.0:
						xP = [linePoint[0], leftPoint[0]]
						yP = [linePoint[1], leftPoint[1]]
						if widthSum > THRESH_WIDTH:
							ax1.plot(xP,yP, color='k')
						else:
							ax1.plot(xP,yP, color='b')

				ax1.plot(pointX,pointY, color='k')
				ax3.plot(angDerivX,angDerivY, color='k')
				ax3.plot(angDeriv2X,angDeriv2Y, color='b')
				ax3.plot(angX,angY, color='r')
				ax2.plot(widthX,widthY, color='k')
				ax2.plot([0.0,widthX[-1]], [THRESH_WIDTH,THRESH_WIDTH], color='r')
				#contigY = [THRESH_WIDTH for k in range(len(contigCenterMass))]
				#ax2.scatter(contigCenterMass,contigY, color='k')

				#centerPointsX = [longPathWidth[k]["linePoint"][0] for k in contigPointIndex]
				#centerPointsY = [longPathWidth[k]["linePoint"][1] for k in contigPointIndex]
				#angDerivs = [longPathWidth[k]["linePoint"] for k in contigPointIndex]

				#if len(centerPointsX) > 0:
				#	ax1.scatter(centerPointsX,centerPointsY, color='k')

				#for k in range(len(centerPointsX)):
				#	ax2.annotate("%1.2f" % contigDensity[k], xy=(contigCenterMass[k], 0.3), xytext=(contigCenterMass[k], 0.8))
					#ax2.annotate("%1.2f" % contigArea[k], xy=(contigCenterMass[k], 0.3), xytext=(contigCenterMass[k], 0.8))
					#ax2.annotate("%1.2f" % angDerivs[k], xy=(contigCenterMass[k], 0.3), xytext=(contigCenterMass[k], 0.8))


				bloomPoint = results["bloomPoint"]
				if bloomPoint != None:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					asymm = contigAsymmSum[0]
					negArea = contigNegArea[0]
					slopes = contigSlopes[0]
					p1 = slopes[2]
					p2 = slopes[3]
					p3 = slopes[4]
					ax1.scatter([bloomPoint[0],],[bloomPoint[1],], color='k', zorder=4)
					ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='k')
					ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')
					ax2.plot([p1[0],p2[0],p3[0]],[p1[1],p2[1],p3[1]], color='y')
			
				archPoint = results["archPoint"]
				if archPoint != None:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					asymm = contigAsymmSum[0]
					negArea = contigNegArea[0]
					slopes = contigSlopes[0]
					p1 = slopes[2]
					p2 = slopes[3]
					p3 = slopes[4]
					ax1.scatter([archPoint[0],],[archPoint[1],], color='b', zorder=4)
					ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='b')
					ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='b')
					ax2.plot([p1[0],p2[0],p3[0]],[p1[1],p2[1],p3[1]], color='y')

				"""
				print "num contigs:", len(contigPointIndex), len(contigCenterMass), len(contigLen), len(contigArea)
				if len(contigPointIndex) > 1:
					contigPointIndex = []
					contigCenterMass = []
					contigLen = []
					contigArea = []
					contigCounts = []
					contigDensity = []


				" conditions for a bloom detection "
				if len(contigPointIndex) == 1:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					frontBoundIndex, backBoundIndex = contigBoundaries[0]


					if (frontBoundIndex <= 6 or backBoundIndex >= len(longPathWidth)-7) and cLen < 25 and dens >= 0.03:
						pointIndex = contigPointIndex[0]
						linePoint = longPathWidth[pointIndex]["linePoint"]

						bloomPoint = results["bloomPoint"]
						ax1.scatter([bloomPoint[0],],[bloomPoint[1],], color='k')
						ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='k')
						ax2.annotate("%1.2f %1.2f" % (area, dens), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')
						print self.nodeID, "boundaries:", frontBoundIndex, backBoundIndex, len(longPathWidth), cLen, area, dens

				"""

				#if len(contigPointIndex) > 0 and (contigCounts[0] <= 5 or contigCounts[0] >= 20 or contigArea[0] < 0.20):
				#if len(contigPointIndex) > 0 and (contigCounts[0] <= 5 or contigCounts[0] >= 30 or contigDensity[0] < 0.05):
				#	contigPointIndex = []
				#	contigCenterMass = []
				#	contigLen = []
				#	contigArea = []
				#	contigCounts = []
				#	contigDensity = []

				#if maxDeriv >= 0.25:
				#	ax1.scatter([inflectionPoint[0]],[inflectionPoint[1]], color='r')
				#	ax1.annotate("%1.2f" % maxDeriv, xy=(inflectionPoint[0], 0.3), xytext=(inflectionPoint[0], 0.8))

				maxDeriv = results["maxDeriv"]
				inflectionU = results["maxDerivU"]
				inflectionPoint = results["inflectionPoint"]
				if inflectionPoint != None:
					ax1.scatter([inflectionPoint[0]],[inflectionPoint[1]], color='r')
					ax1.annotate("%1.2f" % maxDeriv, xy=(inflectionPoint[0], 0.3), xytext=(inflectionPoint[0], 0.8), color='r', zorder=4)

					ax3.scatter([inflectionU,],[maxDeriv,], color='k')

					#ax1.scatter([bloomPoint[0],],[bloomPoint[1],], color='k')
					#ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='k')
					#ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')
					#ax2.plot([p1[0],p2[0],p3[0]],[p1[1],p2[1],p3[1]], color='y')


			ax1.set_title("Medial %d, %s, %d" % (self.nodeID, repr(contigBoundaries), numPoints))
			plt.savefig("medialOut_%04u.png" % self.nodeID)
			plt.clf()
			plt.close()


		greatCount, divSums = bowtieValues[0]
		print "divSums:", greatCount, divSums

		if self.isBowtie:
			pass
			print "BOWTIE"
		
		elif divSums[2] < divSums[0] and divSums[2] < divSums[4]:
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
		self.spatialFeatures = []
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

		#maxPath = deepcopy(longPaths[0])
		#maxPath = deepcopy(medialLongPaths[0])
		maxPath = deepcopy(medialTailCuts[0])

		self.medialCComputed = True		
		self.medialPathC = maxPath
		self.medialPathCut = maxPath


		if True:

			fig, (ax1, ax2, ax3) = plt.subplots(3,  sharex=False, sharey=False)

			ax1.set_xlim(-3, 3)
			ax1.set_aspect("equal")

			ax2.set_ylim(0.0, 1.0)
			ax2.set_aspect("equal")

			ax3.set_ylim(-pi, pi)
			#ax3.set_aspect("equal")


			xP = []
			yP = []
			#for p in medial2:
			#for p in realPath:
			for p in maxPath:
				xP.append(p[0])
				yP.append(p[1])

			ax1.plot(xP,yP, color='b')
	
			xP = []
			yP = []
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			ax1.plot(xP,yP, color='r')

			if True:

				longPathWidth = longMedialWidths[0]

				results = self.spatialFeatures[0]

				maxDeriv = results["maxDeriv"]
				contigCenterMass = results["centerMasses"]
				contigCounts = results["contigCounts"]
				contigDensity = results["contigDensities"]
				contigArea = results["contigAreas"]
				contigLen = results["contigLens"]
				contigPointIndex = results["contigPointIndex"]
				contigBoundaries = results["contigBoundaries"]
				longPathWidth = results["medialWidths"]
				contigAsymmSum = results["contigAsymmSum"]
				contigNegArea = results["contigNegArea"]
				contigSlopes = results["contigSlopes"]

				numPoints = len(longPathWidth)

				for kIndex in range(len(contigCenterMass)):

					area = contigArea[kIndex]
					cLen = contigLen[kIndex]
					frontBoundIndex, backBoundIndex = contigBoundaries[kIndex]

				pointX = []
				pointY = []
				widthX = []
				widthY = []
				angDerivX = []
				angDerivY = []
				angDeriv2X = []
				angDeriv2Y = []
				angX = []
				angY = []

				#for val in longPathWidth:
				for kIndex in range(len(longPathWidth)):

					val = longPathWidth[kIndex]
					distR = val["distR"]
					distL = val["distL"]
					linePoint = val["linePoint"]
					rightPoint = val["rightPoint"]
					leftPoint = val["leftPoint"]
					widthSum = val["widthSum"]
					angDeriv = val["angDeriv"]
					angDeriv2 = val["angDeriv2"]
					distU = val["distU"]

					pointX.append(linePoint[0])
					pointY.append(linePoint[1])


					widthX.append(distU)
					widthY.append(widthSum)
					angDerivX.append(distU)
					angDerivY.append(angDeriv)
					angDeriv2X.append(distU)
					angDeriv2Y.append(angDeriv2)
					angX.append(distU)
					angY.append(linePoint[2])

					#if len(widthX) > 0:
					#	widthX.append(widthX[-1] + 0.04)
					#else:
					#	widthX.append(0.0)

					#print "distR:", distR
					#print "distL:", distL
					#print "linePoint:", linePoint
					#print "leftPoint:", leftPoint
					#print "rightPoint:", rightPoint

					if len(rightPoint) > 0.0:
						xP = [linePoint[0], rightPoint[0]]
						yP = [linePoint[1], rightPoint[1]]
						if widthSum > THRESH_WIDTH:
							ax1.plot(xP,yP, color='k')
						else:
							ax1.plot(xP,yP, color='b')

					if len(leftPoint) > 0.0:
						xP = [linePoint[0], leftPoint[0]]
						yP = [linePoint[1], leftPoint[1]]
						if widthSum > THRESH_WIDTH:
							ax1.plot(xP,yP, color='k')
						else:
							ax1.plot(xP,yP, color='b')

				ax1.plot(pointX,pointY, color='k')
				ax3.plot(angDerivX,angDerivY, color='k')
				ax3.plot(angDeriv2X,angDeriv2Y, color='b')
				ax3.plot(angX,angY, color='r')
				ax2.plot(widthX,widthY, color='k')
				ax2.plot([0.0,widthX[-1]], [THRESH_WIDTH,THRESH_WIDTH], color='r')
				#contigY = [THRESH_WIDTH for k in range(len(contigCenterMass))]
				#ax2.scatter(contigCenterMass,contigY, color='k')

				#centerPointsX = [longPathWidth[k]["linePoint"][0] for k in contigPointIndex]
				#centerPointsY = [longPathWidth[k]["linePoint"][1] for k in contigPointIndex]
				#angDerivs = [longPathWidth[k]["linePoint"] for k in contigPointIndex]

				#if len(centerPointsX) > 0:
				#	ax1.scatter(centerPointsX,centerPointsY, color='k')

				#for k in range(len(centerPointsX)):
				#	ax2.annotate("%1.2f" % contigDensity[k], xy=(contigCenterMass[k], 0.3), xytext=(contigCenterMass[k], 0.8))
					#ax2.annotate("%1.2f" % contigArea[k], xy=(contigCenterMass[k], 0.3), xytext=(contigCenterMass[k], 0.8))
					#ax2.annotate("%1.2f" % angDerivs[k], xy=(contigCenterMass[k], 0.3), xytext=(contigCenterMass[k], 0.8))


				bloomPoint = results["bloomPoint"]
				if bloomPoint != None:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					asymm = contigAsymmSum[0]
					negArea = contigNegArea[0]
					slopes = contigSlopes[0]
					p1 = slopes[2]
					p2 = slopes[3]
					p3 = slopes[4]
					ax1.scatter([bloomPoint[0],],[bloomPoint[1],], color='k', zorder=4)
					ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='k')
					ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')
					ax2.plot([p1[0],p2[0],p3[0]],[p1[1],p2[1],p3[1]], color='y')
			
				archPoint = results["archPoint"]
				if archPoint != None:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					asymm = contigAsymmSum[0]
					negArea = contigNegArea[0]
					slopes = contigSlopes[0]
					p1 = slopes[2]
					p2 = slopes[3]
					p3 = slopes[4]
					ax1.scatter([archPoint[0],],[archPoint[1],], color='b', zorder=4)
					ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='b')
					ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='b')
					ax2.plot([p1[0],p2[0],p3[0]],[p1[1],p2[1],p3[1]], color='y')

				"""
				print "num contigs:", len(contigPointIndex), len(contigCenterMass), len(contigLen), len(contigArea)
				if len(contigPointIndex) > 1:
					contigPointIndex = []
					contigCenterMass = []
					contigLen = []
					contigArea = []
					contigCounts = []
					contigDensity = []


				" conditions for a bloom detection "
				if len(contigPointIndex) == 1:
					dens = contigDensity[0]
					area = contigArea[0]
					cLen = contigLen[0]
					frontBoundIndex, backBoundIndex = contigBoundaries[0]


					if (frontBoundIndex <= 6 or backBoundIndex >= len(longPathWidth)-7) and cLen < 25 and dens >= 0.03:
						pointIndex = contigPointIndex[0]
						linePoint = longPathWidth[pointIndex]["linePoint"]

						bloomPoint = results["bloomPoint"]
						ax1.scatter([bloomPoint[0],],[bloomPoint[1],], color='k')
						ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='k')
						ax2.annotate("%1.2f %1.2f" % (area, dens), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')
						print self.nodeID, "boundaries:", frontBoundIndex, backBoundIndex, len(longPathWidth), cLen, area, dens

				"""

				#if len(contigPointIndex) > 0 and (contigCounts[0] <= 5 or contigCounts[0] >= 20 or contigArea[0] < 0.20):
				#if len(contigPointIndex) > 0 and (contigCounts[0] <= 5 or contigCounts[0] >= 30 or contigDensity[0] < 0.05):
				#	contigPointIndex = []
				#	contigCenterMass = []
				#	contigLen = []
				#	contigArea = []
				#	contigCounts = []
				#	contigDensity = []

				#if maxDeriv >= 0.25:
				#	ax1.scatter([inflectionPoint[0]],[inflectionPoint[1]], color='r')
				#	ax1.annotate("%1.2f" % maxDeriv, xy=(inflectionPoint[0], 0.3), xytext=(inflectionPoint[0], 0.8))

				maxDeriv = results["maxDeriv"]
				inflectionU = results["maxDerivU"]
				inflectionPoint = results["inflectionPoint"]
				if inflectionPoint != None:
					ax1.scatter([inflectionPoint[0]],[inflectionPoint[1]], color='r')
					ax1.annotate("%1.2f" % maxDeriv, xy=(inflectionPoint[0], 0.3), xytext=(inflectionPoint[0], 0.8), color='r', zorder=4)

					ax3.scatter([inflectionU,],[maxDeriv,], color='k')

					#ax1.scatter([bloomPoint[0],],[bloomPoint[1],], color='k')
					#ax2.scatter([contigCenterMass[0],],[THRESH_WIDTH,], color='k')
					#ax2.annotate("%1.2f %1.2f %1.2f %1.3f %s" % (dens, asymm, area, negArea, repr(slopes[0]+slopes[1])), xy=(contigCenterMass[0], 0.3), xytext=(contigCenterMass[0], 0.8), color='k')
					#ax2.plot([p1[0],p2[0],p3[0]],[p1[1],p2[1],p3[1]], color='y')


			ax1.set_title("Medial %d, %s, %d" % (self.nodeID, repr(contigBoundaries), numPoints))
			plt.savefig("medialOut_%04u.png" % self.nodeID)
			plt.clf()
			plt.close()

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
			#pylab.savefig("medialOut_%04u.png" % self.nodeID)


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


	def computeStaticAlphaBoundary(self):

		"""
		1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		2. for desired joint configuration, set occupancy to obstacle if not already free space
		3. for actual joint configuration, set to free space, even if previously set as obstacle
		"""

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
		plotIter = False
		
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

