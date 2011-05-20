
from LocalOccMap import *
from LocalBoundaryMap import *
from LocalObstacleMap import *
from SplineFit import *
from Pose import Pose

import estimateMotion
from GPACurve import GPACurve
from StableCurve import StableCurve

import random
#import Image
from numpy import array, dot, transpose

estPlotCount = 0



class LocalNode:

	def __init__(self, probe, contacts, nodeID, rootNode, pixelSize, stabilizePose = False, direction = True):

		self.pixelSize = pixelSize

		self.stabilizePose = stabilizePose
		self.robotParam = probe.robotParam
		self.numSegs = self.robotParam['numSegs']
		self.segLength = self.robotParam['segLength']
		self.mapSize = self.segLength*self.numSegs + 2.0 + 2.0
		
		self.direction = direction

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
		self.sweepMap = LocalOccMap(self, direction = self.direction)
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
		
		if self.stabilizePose:
			self.initialGPACPose = self.getLocalGPACPose()
	
			print "rootPose:", self.rootPose


		" alpha hull boundary "
		self.a_vert = []
		self.hullComputed = False
		self.computedSweep = False
		
		#pylab.clf()

		#xP = []
		#yP = []
		#for p in self.localPosture:
		#	xP.append(p[0])
		#	yP.append(p[1])
		
		#pylab.plot(xP,yP, color='b')
		
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

		#f = open("frontProbeError%04u.txt" % self.nodeID, 'w')
		#f.write(repr(self.frontProbeError))
		#f.write("\n")
		#f.close()

		#f = open("backProbeError%04u.txt" % self.nodeID, 'w')
		#f.write(repr(self.backProbeError))
		#f.write("\n")
		#f.close()

		f = open("correctedPosture%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.correctedPosture))
		f.write("\n")
		f.close()

		f = open("rootPose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.rootPose))
		f.write("\n")
		f.close()

	
	def readFromFile(self, dirName, nodeID):
		
		self.nodeID = nodeID
		#self.poseProfile.readFromFile("prof%04u.txt" % self.nodeID)
		
		" occupancy map "
		self.occMap.readFromFile(dirName)
	
		" obstacle map "
		self.obstacleMap.readFromFile(dirName)

		" occupancy map "
		self.sweepMap.readFromFile(dirName)

		
		f = open(dirName + "/estpose%04u.txt" % self.nodeID, 'r')
		estPose = eval(f.read().rstrip())
		f.close()

		#self.setEstPose(estPose)

		f = open(dirName + "/gndpose%04u.txt" % self.nodeID, 'r')
		gndPose = eval(f.read().rstrip())
		f.close()

		self.setEstPose(estPose)
		self.setGndPose(gndPose)

		f = open(dirName + "/posture%04u.txt" % self.nodeID, 'r')
		self.localPosture = eval(f.read().rstrip())
		f.close()


		f = open(dirName + "/correctedPosture%04u.txt" % self.nodeID, 'r')
		self.correctedPosture = eval(f.read().rstrip())
		f.close()

		f = open(dirName + "/rootPose%04u.txt" % self.nodeID, 'r')
		self.rootPose = eval(f.read().rstrip())
		f.close()
		
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
		
	def update(self, isForward = True):
		
		if self.stabilizePose:
			self.stabilizeRootPose()
			
		self.occMap.update()
		self.sweepMap.update()
		
		self.dirty = True
		self.hullComputed = False
	
	def isDirty(self):
		return self.dirty
	
	def synch(self):
		
		self.occMap.buildMap()
		self.sweepMap.buildMap()
		
		self.computeAlphaBoundary()

		self.boundaryMap.update()
		self.obstacleMap.update()

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
		

	def computeAlphaBoundary(self, sweep = False):

		if self.hullComputed and self.computedSweep == sweep:
			" already computed "
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
				self.a_vert = self.computeAlpha2(points, radius = 0.02)
			else:
				self.a_vert = self.computeAlpha2(points, radius = 0.2)
			" cut out the repeat vertex "
			self.a_vert = self.a_vert[:-1]
			
			self.a_vert = self.convertAlphaUniform(self.a_vert)
		
		else:
			self.a_vert = []
			
		self.hullComputed = True
		self.computedSweep = sweep
			
			
	def getAlphaBoundary(self):
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
		
		#print inputStr
		
		if True:
			
			" start the subprocess "
			subProc = Popen(["./alpha2.exe"], stdin=PIPE, stdout=PIPE)
			
			" send input and receive output "
			sout, serr = subProc.communicate(inputStr)
	
			#print numPoints
			#print sout
			
			" convert string output to typed data "
			sArr = sout.split(" ")
			
			#print sArr[0]
			numVert = int(sArr[0])
			
			vertices = []
			for i in range(numVert+1):
				vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
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
		
			
		return vertices
							