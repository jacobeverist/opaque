import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from LocalOccMap import *
from LocalBoundaryMap import *
from LocalObstacleMap import *
from LocalVoronoiMap import *
from LocalNavRoadMap import *
from LocalFrontierMap import *

import pose

#import Image


class LocalNode:

	def __init__(self, probe, contacts, nodeID, rootNode, inSim = True):

		self.nodeID = nodeID
		self.probe = probe
		self.contacts = contacts
		self.rootNode = rootNode
		#self.rootPose = [0.0,0.0,0.0]
		
		self._refnodes = []
		self._refent = []
		
		"""
		initial information for the local node
		
		1. estimated pose
		2. profile of the anchor points 
		"""
		
		if inSim:
			self.setEstPose(self.contacts.getClosestPose(self.rootNode))			
			self.setGndPose(self.probe.getActualJointPose(self.rootNode))

		self.poseProfile = pose.PoseProfile(self.contacts,self.rootNode, inSim)

		# MAPS
		self.occMap = LocalOccMap(self)
		self.boundaryMap = LocalBoundaryMap(self)
		self.obstacleMap = LocalObstacleMap(self)
		self.voronoiMap = LocalVoronoiMap(self)
		self.navMap = LocalNavRoadMap(self)
		self.frontierMap = LocalFrontierMap(self)
		
		self.saveCount = 0
	
		" by default, the center points are the estimated positions of joints: "
		cJoints = [2, 6, 10, 14, 18, 22, 26, 30, 34, 38]

		if inSim:
			self.centerPoints = []
			for j in cJoints:
				self.centerPoints.append(self.contacts.getClosestPose(j))
			
			self.localCenterPoints = []
	
			for pnt in self.centerPoints:
				newPoint = self.convertGlobalToLocal(pnt)
				self.localCenterPoints.append(newPoint)

	def setEstPose(self, newPose):

		self.estPose = copy(newPose)
		self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)
		self.vecAng = acos(self.estPose[0]/self.dist)
		if asin(self.estPose[1]/self.dist) < 0:
			self.vecAng = -self.vecAng
		
		self.backR = array([[cos(self.vecAng), sin(self.vecAng)],[-sin(self.vecAng),cos(self.vecAng)]])
		self.foreR = array([[cos(self.vecAng), -sin(self.vecAng)],[sin(self.vecAng),cos(self.vecAng)]])
		
		self.R = array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])

	def setGndPose(self, newPose):
		
		self.gndPose = copy(newPose)
		self.gndDist = sqrt(self.gndPose[0]**2 + self.gndPose[1]**2)
		self.gndVecAng = acos(self.gndPose[0]/self.gndDist)
		if asin(self.gndPose[1]/self.gndDist) < 0:
			self.gndVecAng = -self.gndVecAng
		
		self.gndBackR = array([[cos(self.gndVecAng), sin(self.gndVecAng)],[-sin(self.gndVecAng),cos(self.gndVecAng)]])
		self.gndForeR = array([[cos(self.gndVecAng), -sin(self.gndVecAng)],[sin(self.gndVecAng),cos(self.gndVecAng)]])
		
		self.gndR = array([[cos(self.gndPose[2]), sin(self.gndPose[2])],[-sin(self.gndPose[2]),cos(self.gndPose[2])]])
		
	def convertLocalOffsetToGlobal(self, offset):

	   globalEst = [0.0,0.0,0.0]

	   finalVec = array([[offset[0]], [offset[1]]])
	   transVec = dot(transpose(self.R), finalVec)
	   resVec = dot(self.backR, transVec)
	   resVec[0, 0] += self.dist
	   tempVec = dot(self.foreR, resVec)
	   globalEst[0] = tempVec[0, 0]
	   globalEst[1] = tempVec[1, 0]
	   globalEst[2] = normalizeAngle(self.estPose[2] + offset[2])

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

	   localPose = [localVec[0,0], localVec[1,0], normalizeAngle(pose[2] - self.estPose[2])]

	   return localPose

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
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1
		
		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0		

		joints = range(-1,self.rootNode)
		joints.reverse()

		for i in joints:
			totalAngle = totalAngle + self.probe.getServo(i+1)
			totalAngle = normalizeAngle(totalAngle)
			
			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

		return [xTotal, zTotal, totalAngle]

	def getAftTip(self):
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1
		
		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0		

		joints = range(self.rootNode, 39)

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)

			if i < 38:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)
						
		return [xTotal, zTotal, totalAngle]

	def getJointPose(self, jointI):
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1
		
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
				totalAngle = normalizeAngle(totalAngle)
		
		elif jointI < self.rootNode:
			
			joints = range(jointI,self.rootNode)
			joints.reverse()
	
			for i in joints:
				totalAngle = totalAngle + self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)
				
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

	def getVoronoiMap(self):
		return self.voronoiMap

	def getFrontierMap(self):
		return self.frontierMap
		
	def getNavMap(self):
		return self.navMap
	
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
		self.poseProfile.saveToFile("prof%04u.txt" % self.nodeID)
		
		f = open("estpose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.estPose))
		f.write("\n")
		f.close()

		f = open("gndpose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.gndPose))
		f.write("\n")
		f.close()

		f = open("cntpose%04u.txt" % self.nodeID, 'w')
		f.write(repr(self.localCenterPoints))
		f.write("\n")
		f.close()
	
	def readFromFile(self, nodeID):
		
		self.nodeID = nodeID
		self.poseProfile.readFromFile("prof%04u.txt" % self.nodeID)
		
		" occupancy map "
		self.occMap.readFromFile()
	
		" obstacle map "
		self.obstacleMap.readFromFile()
		
		f = open("estpose%04u.txt" % self.nodeID, 'r')
		estPose = eval(f.read())
		f.close()

		self.setEstPose(estPose)

		f = open("gndpose%04u.txt" % self.nodeID, 'r')
		gndPose = eval(f.read())
		f.close()

		self.setGndPose(gndPose)

		self.poseProfile.setRootPose(self.estPose)
		
		" read in the data "
		f = open("cntpose%04u.txt" % self.nodeID, 'r')
		centerPoints = eval(f.read())
		f.close()
		
		#print "reading in center points", "cntpose%04u.txt" % self.nodeID
		#print centerPoints
		
		self.setCenterPoints(centerPoints)

		self.synch()
		
	def getCenterPoints(self):
		return copy(self.centerPoints)
	
	def setCenterPoints(self, centerPoints):
		self.centerPoints = copy(centerPoints)
		
		"""
		for i in range(len(self.centerPoints)):
			
			## Create the visual representation of active reference nodes
			name = "cnt_node" + str(i) + "_" + str(self.nodeID)
			entity = self.probe._mgr.createEntity(name, "Cube.mesh")
			node = self.probe._mgr.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)

			size = ogre.Vector3(0.05,0.05,0.05)
			node.setScale(size)

			pos = self.centerPoints[i]
			position = ogre.Vector3(pos[0],0.1,pos[1])
			node.setPosition(position)
			# negative the angle because the bodies are somehow flipped upside down
			node.setOrientation(ogre.Quaternion(0.0,ogre.Vector3().UNIT_Y))

			entity.setCastShadows(False)
			entity.setMaterialName("Green")

			#entity.setVisible(False)
			self._refnodes.append(node)
			self._refent.append(entity)
		"""
		
		self.localCenterPoints = []

		for pnt in self.centerPoints:
			newPoint = self.convertGlobalToLocal(pnt)
			self.localCenterPoints.append(newPoint)
			
		

	def getPoseProfile(self):
		return self.poseProfile
		
	def getEstPose(self):
		return copy(self.estPose)

	def obstCallBack(self, direction):
		#self.currNode.obstCallBack(direction)		
		self.obstacleMap.updateContact(direction)
		
	def update(self, isForward = True):
		
		" update the local occupancy "
		#if isForward:
		#	newRoot = self.contacts.getClosestPose(34)
		#else:
		#	newRoot = self.contacts.getClosestPose(5)

		"""
		newRoot = self.contacts.getClosestPose(19)
		
		" transform newRoot to local coordinates"
		relVec = [newRoot[0],newRoot[1]]
		relDist = sqrt(relVec[0]**2 + relVec[1]**2)
		
		tempVec = array([[relVec[0]],[relVec[1]]])
		resVec = dot(self.backR,tempVec)
		resVec[0,0] += -self.dist
		transVec = dot(self.foreR, resVec)
		
		" now, apply rotation correction "
		finalVec = dot(self.R, transVec)
		
		finalAng = normalizeAngle( normalizeAngle(newRoot[2]) - normalizeAngle(self.estPose[2]) )
		
		relRoot = [finalVec[0,0],finalVec[1,0],finalAng]
		
		#if isForward:
		#	self.occMap.setRootPose(relRoot, 34)
		#else:
		#	self.occMap.setRootPose(relRoot, 5)
		
		self.occMap.setRootPose(relRoot, 19)
		"""	
		
		self.occMap.update(isForward)
		#self.boundaryMap.update()
		#self.obstacleMap.update()
		
		#self.voronoiMap.update(self.boundaryMap.getBoundaryPoints(), self.occMap.getImage())
		#self.obstacleMap.update(self.boundaryMap)
		#self.frontierMap.update(self.boundaryMap, self.obstacleMap)
		#self.navRoadMap.update(self.voronoiMap.getGraph())
	
	def synch(self):
		self.occMap.buildMap()
		
		self.computeAlphaBoundary()

		self.boundaryMap.update()
		self.obstacleMap.update()
		self.voronoiMap.update()
		self.frontierMap.update()
		self.navMap.update()
		
	def saveMap(self):
			
		" FIXME: uncomment these to build the other maps "
		
		" synchronize maps first "
		self.synch()

		# save a copy of each version of the map
		#self.saveCount += 1
		self.occMap.saveMap()
		self.boundaryMap.saveMap()
		self.obstacleMap.saveMap()
		self.voronoiMap.saveMap()
		self.frontierMap.saveMap()
		
		
		self.computeGndAlphaBoundary()
		val = repr(self.a_vert)
		f = open("gnd_alpha_bound_%04u.txt" % self.nodeID, 'w')
		f.write(val)
		f.close()

		val = repr(self.gndOccPoints)
		f = open("gnd_occ_points_%04u.txt" % self.nodeID, 'w')
		f.write(val)
		f.close()

				
		self.computeAlphaBoundary()
		val = repr(self.a_vert)
		f = open("alpha_bound_%04u.txt" % self.nodeID, 'w')
		f.write(val)
		f.close()

		val = repr(self.estOccPoints)
		f = open("est_occ_points_%04u.txt" % self.nodeID, 'w')
		f.write(val)
		f.close()

	def computeAlphaBoundary(self):

		" 1. pick out the points "
		numPixel = self.occMap.numPixel
		mapImage = self.occMap.getMap()
		image = mapImage.load()
		
		#print numPixel
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					pnt = self.occMap.gridToReal([j,k])
					points.append(pnt)
		
		self.estOccPoints = points
		#print len(points)
		
		if len(points) > 0:
			#print points
		
			self.a_vert = self.computeAlpha2(points)
	
			" cut out the repeat vertex "
			self.a_vert = self.a_vert[:-1]
			
			self.convertAlphaUniform()
		
		else:
			self.a_vert = []
			
	def getAlphaBoundary(self):
		return self.a_vert

	def computeGndAlphaBoundary(self):

		" 1. pick out the points "
		numPixel = self.occMap.numPixel
		mapImage = self.occMap.getGndMap()
		image = mapImage.load()
		
		#print numPixel
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					pnt = self.occMap.gridToReal([j,k])
					points.append(pnt)
		
		self.gndOccPoints = points
		print len(points)
		
		if len(points) > 0:
			#print points
		
			self.a_vert = self.computeAlpha2(points)
	
			" cut out the repeat vertex "
			self.a_vert = self.a_vert[:-1]
			
			self.convertAlphaUniform()
		
		else:
			self.a_vert = []
				
	def convertAlphaUniform(self):
		
		" make the vertices uniformly distributed "
		
		new_vert = []
		
		max_spacing = 0.04
		#max_spacing = 0.1
		
		for i in range(len(self.a_vert)):
			p0 = self.a_vert[i]
			p1 = self.a_vert[(i+1) % len(self.a_vert)]
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
					
		self.a_vert = new_vert
	
	def computeAlpha2(self, points):
		
		numPoints = len(points)
		inputStr = str(numPoints) + " "

		" radius 0.2 "
		#inputStr += str(1.0) + " "
		inputStr += str(0.8) + " "
		
		for p in points:
			p2 = copy(p)
			p2[0] += gauss(0.0,0.0001)
			p2[1] += gauss(0.0,0.0001)

			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		
		inputStr += "\n"
		
		#print inputStr
		
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
		
		#print vertices
		return vertices
							