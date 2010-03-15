import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from RefNode import *
from RefEdge import *
from ValueStability import *
from PoseProfile import *

import math
import cProfile
import pylab
from copy import *

import ogre.renderer.OGRE as ogre

class AverageContacts:
	
	def __init__(self, probe):

		self.timeInc = 1
		self.count = 0
		
		self.prof = cProfile.Profile()

		# for drawing into the engine the points
		self.markers = []

		# details of probe
		self.probe = probe
		self.numJoints = self.probe.getNumJoints()

		# data output of angles over time
		#self.angle_output = open("frame/angle_output.txt", "w")

		# contact parameter arrays
		num = self.numJoints

		self.errorTime = [0 for i in range(num)]
		self.maskTime = [0 for i in range(num)]
		self.avgErr = [0.0 for i in range(num)]
		self.mask = [False for i in range(num)]
		self.maxErrorReached = [False for i in range(num)]
		self.recvData = [False for i in range(num)]

		# reference node ids
		self.originNode = 0
		self.lastNode = 0

		# active reference nodes
		self.numActiveRef = 0

		self.activeRef = [False for i in range(num)]
		self.activeRefPtr = [0 for i in range(num)]

		# proposed reference points
		self.numRef = 0
		self.refPoints = []
		
		# ref point constraints
		self.numEdge = 0
		self.refEdges = []
		
		# tunable constants
		self.errorThresh = 0.15
		self.maskThresh = 0.98
		self.timeThresh = 50.0

		# variance threshold for determining joint stability
		self.varThresh = 0.001

		# forgetting factor of error average
		self.K = 0.05

		# number of data constructs we produce
		self.graphCount = 0
		self.poseCount = 0

		# variance estimator for determining stability
		self.jointVar = []
		for i in range(self.numJoints):
			self.jointVar.append(ValueStability(self.varThresh))
	
		self.L_saddlePoses = []
		self.R_saddlePoses = []
		self.inflectionPoses = []

		self.poses = []
		self.offsets = []
		self.destinationReached = False
		self.tipPoint = [0,0,0]

		#self.resetPose()
		
		self._allRefNodes = []
		self._allRefEnt = []
		
		"""
		for i in range(0,self.numJoints):
			## Create the visual reprsentation of active reference nodes
			name = "act_node" + str(i)
			entity = self.probe._mgr.createEntity(name, "Cube.mesh")
			node = self.probe._mgr.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)

			size = ogre.Vector3(0.05,0.05,0.05)
			node.setScale(size)

			pos = self.probe.getActualJointPose(i)
			position = ogre.Vector3(pos[0],0.1,pos[1])
			node.setPosition(position)
			# negative the angle because the bodies are somehow flipped upside down
			node.setOrientation(ogre.Quaternion(-pos[2],ogre.Vector3().UNIT_Y))

			entity.setCastShadows(False)

			entity.setMaterialName("Red")

			#entity.setVisible(False)
			self._allRefNodes.append(node)
			self._allRefEnt.append(entity)
		"""
		
	def setTimerAliasing(self, timeInc):
		self.timeInc = timeInc

	def isStep(self):
		if self.count % self.timeInc == 0:
			return True
		
		return False

	def getAverageSegPose(self, targetSegment):
	
		if targetSegment == 0:
			pose = self.getAveragePose(0)

			pose[2] = pose[2] + self.probe.getServo(0)
			pose[2] = self.normalizeAngle(pose[2])
			pose[0] = pose[0] - self.probe.segLength*cos(pose[2])
			pose[1] = pose[1] - self.probe.segLength*sin(pose[2])

			return pose
		
		pose = self.getAveragePose(targetSegment-1)

	
		# adding constant offset so that we are on the joint
		pose[0] = pose[0] + (self.probe.segLength/2)*cos(pose[2])
		pose[1] = pose[1] + (self.probe.segLength/2)*sin(pose[2])

		return pose

	def getAveragePose(self, targetJoint):
		
		# find the first active reference from which to determine the new node's position in space
		refExists = False
		refNumber = 0
		isNewRoot = False
		refCount = 0


		# search for the largest consecutive section of stable reference nodes
		nodes = range(0,self.numJoints-1)

		prevActive = False
		currStart = 0
		currSize = 0
		maxConsec = 0
		maxSize = 0
		
		" select the largest contiguous set of reference points to compute the average pose "
		for i in nodes:
			if self.activeRef[i]:
				currStart = i
				if prevActive:
					currSize += 1
						
					if currSize > maxSize:
						maxSize = currSize
						maxConsec = currStart
				else:
					prevActive = True
					currSize = 1
			else:
				currSize = 0
				currStart = 0
				prevActive = False
		
		if maxSize > 0:
			refExists = True
			refNumber = maxConsec

		# record this as new reference point
		pose = [0.,0.,0.]

		nodes = range(maxConsec-maxSize+1,maxConsec+1)

		refCount = maxSize
		
		if refExists:
			sumX = 0.0
			sumY = 0.0
			sumPx = 0.0
			sumPy = 0.0
			
			for j in nodes:
				
				if j != targetJoint and self.activeRef[j]:
					#print "ref ", j
					# with an existing reference point, we determine the new node's position with respect to it
					refX, refZ, refP = self.activeRefPtr[j].getRefPose()
			
					initPose = [refX,refZ,refP]
					jointID = self.activeRefPtr[j].getJointID()
					pose = self.probe.getJointWRTJointPose(initPose, jointID, targetJoint) 
					sumX += pose[0]
					sumY += pose[1]
					sumPx += cos(pose[2])
					sumPy += sin(pose[2])
					
					# the existing reference node saved for debugging purposes 
					self.lastNode = self.activeRefPtr[refNumber]
			
			sumX /= refCount
			sumY /= refCount
			mag = sqrt(sumPx**2 + sumPy**2)
			sumPx /= mag
			sumPy /= mag

			vecAng = acos(sumPx)
			if asin(sumPy) < 0:
				vecAng = -vecAng
				
			pose = [sumX, sumY, vecAng]
			
		# this is the first node, so reference from 0,0
		elif self.numRef == 0:
			
			if forcePose != []:
				pose = forcePose
				print "taking forced pose of", forcePose
			else:
				# initialize the first reference node with its position in global coordinates
				pose = self.probe.getActualJointPose(targetJoint)

		else:
			# no active node exists and not the first,
			# so lets use the latest dormant one
			recentNode = self.refPoints[self.numRef-1]	

			refX, refZ, refP = recentNode.getRefPose()
			initPose = [refX, refZ, refP]

			jointID = recentNode.getJointID()
			pose = self.probe.getJointWRTJointPose(initPose, jointID, targetJoint) 

		return pose

	
	def getClosestPose(self, targetJoint):
		# find the first active reference from which to determine the new node's position in space
		refExists = False
		refNumber = 0
		isNewRoot = False

		if targetJoint == -1:
			midJoint = 0
		elif targetJoint == self.numJoints:
			midJoint = self.numJoints-1
		else:
			midJoint = targetJoint
			
			if self.activeRef[midJoint]:
				refX, refZ, refP = self.activeRefPtr[midJoint].getRefPose()
				return [refX, refZ, refP]

		# FIXME: find the *closest* recent node, not the *first*
		nodes = range(0,targetJoint-1)
		nodes.reverse()

		for j in nodes:

			if self.activeRef[j]:
				refExists = True
				refNumber = j
				break

		if not refExists:
			nodes = range(0,self.numJoints-1)
			nodes.reverse()
			for j in nodes:
				if self.activeRef[j]:
					refExists = True
					refNumber = j
					break

		# record this as new reference point
		pose = [0.,0.,0.]

		if refExists:

			# with an existing reference point, we determine the new node's position with respect to it
			refX, refZ, refP = self.activeRefPtr[refNumber].getRefPose()

			initPose = [refX,refZ,refP]
			jointID = self.activeRefPtr[refNumber].getJointID()
			pose = self.probe.getJointWRTJointPose(initPose, jointID, targetJoint) 
			return pose
		
		else:
			# no active node exists and not the first,
			# so lets use the latest dormant one
			recentNode = self.refPoints[self.numRef-1]	

			refX, refZ, refP = recentNode.getRefPose()
			initPose = [refX, refZ, refP]

			jointID = recentNode.getJointID()
			pose = self.probe.getJointWRTJointPose(initPose, jointID, targetJoint) 
			return pose
		
	def resetPose(self, estPose = []):
		
		# contact parameter arrays
		num = self.numJoints

		# reference node ids
		self.originNode = 0
		self.lastNode = 0

		# active reference nodes
		self.numActiveRef = 0

		self.activeRef = [False for i in range(num)]
		self.activeRefPtr = [0 for i in range(num)]

		# deference the child nodes now
		for child in self._allRefNodes:
			self.probe._mgr.destroySceneNode(child)

		self._allRefNodes = []
	
		for child in self._allRefEnt:
			self.probe._mgr.destroyEntity(child)
	
		self._allRefEnt = []

		# proposed reference points
		self.numRef = 0
		self.refPoints = []
		
		# ref point constraints
		self.numEdge = 0
		self.refEdges = []

		
		# number of data constructs we produce
		self.graphCount = 0
		self.poseCount = 0

		self.L_saddlePoses = []
		self.R_saddlePoses = []
		self.inflectionPoses = []

		self.poses = []
		self.offsets = []
		self.destinationReached = False
		self.tipPoint = [0,0,0]

		startJoint = 2
		endJoint = self.numJoints-2

		if estPose != []:
			newNode = self.createNewNodeAverage(19, forcePose = estPose)
			self.createEdges(newNode)


		for i in range(startJoint, endJoint):
			#if not self.activeRef[i] and not self.maxErrorReached[i] and (self.avgErr[i] > self.errorThresh) and (self.mask[i] > self.maskThresh) and (self.maskTime[i] > self.timeThresh):
			if not self.activeRef[i]:

				# create a new reference node
				newNode = self.createNewNode(i)

				# create position constraints between all the current active reference nodes
				self.createEdges(newNode)


	def setMask(self, mask):
		self.initMask = copy(mask)
		#print "setting mask:", self.initMask
		
	def step(self):
		self.count += 1
		if not self.isStep():
			return 
	
		"""
		for i in range(self.numJoints):
			pos = self.probe.getActualJointPose(i)
			position = ogre.Vector3(pos[0],0.1,pos[1])
			self._refnodes[i].setPosition(position)

			# negative the angle because the bodies are somehow flipped upside down
			self._refnodes[i].setOrientation(ogre.Quaternion(-pos[2],ogre.Vector3().UNIT_Y))

			self._refent[i].setVisible(self.activeRef[i])
			#if i <= 27:
			#	self._refent[i].setVisible(True)
			#else:
			#	self._refent[i].setVisible(False)


			mPtr = self._refnodes[i].getMaterial()
			mPtr.setAmbient(1.0,0.0,0.0)
			mPtr.setDiffuse(1.0,0.0,0.0, 1.0)
		"""
		
		if True:

			# amplitude max
			for i in range(self.numJoints):
				self.mask[i] = self.initMask[i]

			# update error timers
			for i in range(self.numJoints):
				if self.probe.getError(i) > self.errorThresh: 
					self.errorTime[i] += 1	
				else:
					self.errorTime[i] = 0

			# update mask timers
			for i in range(self.numJoints):
				if self.mask[i] >= 1.0:
					self.maskTime[i] += 1	
				else:
					self.maskTime[i] = 0

			# compute the position of each joint with a moving average
			for i in range(self.numJoints):
				self.avgErr[i] = (1.0 - self.K) * (self.avgErr[i]) + self.K * self.probe.getError(i)

			# joint variance computation for stability checking
			for i in range(self.numJoints):
				self.jointVar[i].addData(self.probe.getServo(i))

			# include all joints for node creation
			startJoint = 0
			endJoint = self.numJoints
			#startJoint = 2
			#endJoint = self.numJoints-2

			# 1. compute the stability of a length between contact points
			# 2. deactive the whole length
			# 3. hypothesize the slip mode
			# 4. reactivate the nodes when things become stable again

			# saddle points that serve as boundaries between closed sections
			# used for determining the stability of points within the sections
			countPoints = [0,4,8,12,16,20,24,28,32,36,38]
			numCountPoints = 10

			for i in range(startJoint,endJoint):

				if self.activeRef[i]:

					self.activeRefPtr[i].computeStabilityError()

					if self.activeRefPtr[i].isMaxErrorReached():

						self.activeRef[i] = False
						self.numActiveRef -= 1
						self.maxErrorReached[i] = True
						#print "caseA", i

						# reset variance calculation
						self.jointVar[i].reset()

				if self.mask[i] < 1.0 and self.maxErrorReached[i]:
					self.maxErrorReached[i] = False

				# add new node from joint if it crosses threshold
				# only check error threshold on joints that are taut, maxed Gaussian
				if not self.activeRef[i] and not self.maxErrorReached[i] and (self.avgErr[i] > self.errorThresh) and (self.mask[i] > self.maskThresh) and (self.maskTime[i] > self.timeThresh):

					# create a new reference node
					newNode = self.createNewNode(i)

					# create position constraints between all the current active reference nodes
					self.createEdges(newNode)

				elif self.activeRef[i] and self.mask[i] < 1.0 :
					self.activeRef[i] = False
					self.numActiveRef -= 1
					#print "caseB", i, self.mask[i]
					# reset variance calculation
					self.jointVar[i].reset()

			# add angle pose data to the active reference node
			for i in range(self.numJoints):

				if i == 2 or i == 35:

					#print self.activeRef
					#print self.activeRefPtr
					if self.activeRef[i]:

						#if self.count % 2 == 0 :
						if True :
							angleData = []
							for j in range(self.numJoints):
								angleData.append(self.probe.getServo(j))

							#print "referencing ", i , self.activeRefPtr[i]

							self.activeRefPtr[i].addAngleData(angleData)

						self.recvData[i] = True

					elif self.recvData[i]:

						self.recvData[i] = False

						# write this data to the map
						#mapping.addOpenSpace(self.activeRefPtr[i])

			for i in range(startJoint,endJoint):
				if not self.activeRef[i]:

					isStable = True

					low = i - 4
					high = i + 4

					if low < 0:
						low = 0
					if high >= self.numJoints:
						high = self.numJoints

					for j in range(low,high):
						if not self.jointVar[j].isStable():
							isStable = False

					# the neighborhood has stabled, now lets reactivate the reference
					if isStable and self.mask[i] > self.maskThresh and self.maskTime[i] > self.timeThresh:

						# create a new reference node
						newNode = self.createNewNode(i)

						# create position constraints between all the current active reference nodes
						self.createEdges(newNode)

			# compute the active node local configuration error
			#self.computeActiveError()
			
			#self.prof.runcall(self.computeActiveError)
			#self.prof.dump_stats("profile_info2")	

			#for i in range(self.numJoints):
			#	self.angle_output.write(str(self.probe.getServo(i)) + " ")

			#for i in range(self.numJoints):
			#	self.angle_output.write(str(self.activeRef[i]) + " ")

			#self.angle_output.write("\n")

		# save the gray map every 5 frames
		if self.count % 100 == 0:
			pass
			#mapping.saveMap()


	def setErrorThresh(self, thresh):
		self.errorThresh = thresh

	def setTimeThresh(self, thresh):
		self.timeThresh = thresh

	def draw(self):
		pass

	def createEdges(self, newNode):

		jointID = newNode.getJointID()

		# create an edge between every active node
		for j in range(self.numJoints):
			if self.activeRef[j] and jointID != j:

				# 1. compute the displacement and rotation between them
				initPose = [0.,0.,0.]
				targetPose = self.probe.getJointWRTJointPose(initPose, j, jointID) 

				# convert the relative position to foward/slide/rotate coordinates for TORO
				# a) rotate stays the same since our init position is 0/0/0
				# b) forward is the distance in the current orientation
				# c) slide is zero since we are oriented correctly

				forward = targetPose[0]
				slide = targetPose[1]
				rotate = targetPose[2]

				# 2. create the new edge
				newEdge = RefEdge(self.activeRefPtr[j], newNode, forward, slide, rotate)

				nodeDistance = abs(jointID-j)

				if nodeDistance >= 8:
					newEdge.setVariance(0.02)
				else:
					newEdge.setVariance(0.08)

				# 3. save this edge so we can access later
				self.refEdges.append(newEdge)
				self.numEdge += 1

	def createNewNode(self, newJointID):

		# find the first active reference from which to determine the new node's position in space
		refExists = False
		refNumber = 0
		isNewRoot = False

		# FIXME: find the *closest* recent node, not the *first*
		
		nodes = range(0,newJointID-1)
		nodes.reverse()

		for j in nodes:

			if self.activeRef[j]:
				refExists = True
				refNumber = j
				break

		if not refExists:
			nodes = range(0,self.numJoints-1)
			nodes.reverse()
			for j in nodes:
				if self.activeRef[j]:
					refExists = True
					refNumber = j
					break
		
		# set the current node as active
		self.activeRef[newJointID] = True
		self.numActiveRef += 1

		# record this as new reference point
		pose = [0.,0.,0.]

		if refExists:

			# with an existing reference point, we determine the new node's position with respect to it
			refX, refZ, refP = self.activeRefPtr[refNumber].getRefPose()

			initPose = [refX,refZ,refP]
			jointID = self.activeRefPtr[refNumber].getJointID()
			pose = self.probe.getJointWRTJointPose(initPose, jointID, newJointID) 

			# the existing reference node saved for debugging purposes 
			self.lastNode = self.activeRefPtr[refNumber]

		# this is the first node, so reference from 0,0
		elif self.numRef == 0:

			# initialize the first reference node with its position in global coordinates
			origin = self.probe.getActualJointPose(newJointID)
			origin[0] += 2.0
			#print "initial pose at joint", newJointID

			#origin[2] = self.normalizeAngle(origin[2]+math.pi)
			originNode = RefNode(-1,newJointID,origin[0],origin[1],origin[2], self.probe)
			self.lastNode = originNode

			refX, refZ, refP = originNode.getRefPose()
			initPose = [refX, refZ, refP]

			jointID = originNode.getJointID()

			#print "jointID =", jointID, "newJointID =", newJointID

			pose = self.probe.getJointWRTJointPose(initPose, jointID, newJointID) 
			#print initPose,pose

		else:
			# no active node exists and not the first,
			# so lets use the latest dormant one
			recentNode = self.refPoints[self.numRef-1]	

			refX, refZ, refP = recentNode.getRefPose()
			initPose = [refX, refZ, refP]

			jointID = recentNode.getJointID()
			pose = self.probe.getJointWRTJointPose(initPose, jointID, newJointID) 
			self.lastNode = recentNode

			isNewRoot = True

		# create the new node
		nodeID = self.numRef
		newNode = RefNode(nodeID, newJointID, pose[0], pose[1], pose[2], self.probe)	
		newNode.setNewRoot(isNewRoot)
		newNode.setGuideNode(self.lastNode.getNodeID())
		gndPose = self.probe.getActualJointPose(newJointID)
		newNode.setGroundTruthPose(gndPose[0],gndPose[1],gndPose[2])

		self.activeRefPtr[newJointID] = newNode	

		" draw the estimated reference nodes "
		if False:
			i = self.numRef
			## Create the visual reprsentation of active reference nodes
			name = "est_node" + str(i)
			entity = self.probe._mgr.createEntity(name, "Cube.mesh")
			node = self.probe._mgr.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)

			size = ogre.Vector3(0.05,0.05,0.05)
			node.setScale(size)

			position = ogre.Vector3(pose[0],0.1,pose[1])
			node.setPosition(position)
			# negative the angle because the bodies are somehow flipped upside down
			node.setOrientation(ogre.Quaternion(-pose[2],ogre.Vector3().UNIT_Y))

			entity.setCastShadows(False)
			entity.setMaterialName("Green")

			#entity.setVisible(False)
			self._allRefNodes.append(node)
			self._allRefEnt.append(entity)

		# record it to the array so we can access it later
		self.refPoints.append(newNode)	
		self.numRef += 1
		
		return newNode


	def createNewNodeAverage(self, newJointID, forcePose = []):

		# find the first active reference from which to determine the new node's position in space
		refExists = False
		refNumber = 0
		isNewRoot = False
		refCount = 0


		# search for the largest consecutive section of stable reference nodes
		nodes = range(0,self.numJoints-1)

		prevActive = False
		currStart = 0
		currSize = 0
		maxConsec = 0
		maxSize = 0
		
		" select the largest contiguous set of reference points to compute the average pose "
		for i in nodes:
			if self.activeRef[i]:
				currStart = i
				if prevActive:
					currSize += 1
						
					if currSize > maxSize:
						maxSize = currSize
						maxConsec = currStart
				else:
					prevActive = True
					currSize = 1
			else:
				currSize = 0
				currStart = 0
				prevActive = False
		
		if maxSize > 0:
			refExists = True
			refNumber = maxConsec
		
		#print maxSize, maxConsec
		
		# FIXME: find the *closest* recent node, not the *first*
		"""
		nodes = range(0,self.numJoints-1)
		for j in nodes:

			if self.activeRef[j]:
				refExists = True
				refNumber = j
				refCount += 1
		"""

		# set the current node as active
		self.activeRef[newJointID] = True
		self.numActiveRef += 1

		# record this as new reference point
		pose = [0.,0.,0.]

		nodes = range(maxConsec-maxSize+1,maxConsec+1)
		print maxConsec-maxSize+1, maxConsec
		refCount = maxSize
		
		if refExists:
			sumX = 0.0
			sumY = 0.0
			sumPx = 0.0
			sumPy = 0.0
			
			for j in nodes:
				
				if j != newJointID and self.activeRef[j]:
					#print "ref ", j
					# with an existing reference point, we determine the new node's position with respect to it
					refX, refZ, refP = self.activeRefPtr[j].getRefPose()
			
					initPose = [refX,refZ,refP]
					jointID = self.activeRefPtr[j].getJointID()
					pose = self.probe.getJointWRTJointPose(initPose, jointID, newJointID) 
					sumX += pose[0]
					sumY += pose[1]
					sumPx += cos(pose[2])
					sumPy += sin(pose[2])
					
					# the existing reference node saved for debugging purposes 
					self.lastNode = self.activeRefPtr[refNumber]
			
			sumX /= refCount
			sumY /= refCount
			mag = sqrt(sumPx**2 + sumPy**2)
			sumPx /= mag
			sumPy /= mag

			vecAng = acos(sumPx)
			if asin(sumPy) < 0:
				vecAng = -vecAng
				
			pose = [sumX, sumY, vecAng]
			
		# this is the first node, so reference from 0,0
		elif self.numRef == 0:

			if forcePose != []:
				origin = forcePose
				print "taking forced pose of", forcePose
			else:
				
				# initialize the first reference node with its position in global coordinates
				origin = self.probe.getActualJointPose(newJointID)
				#print "initial pose at joint", newJointID

			#origin[2] = self.normalizeAngle(origin[2]+math.pi)
			originNode = RefNode(-1,newJointID,origin[0],origin[1],origin[2], self.probe)
			self.lastNode = originNode

			refX, refZ, refP = originNode.getRefPose()
			initPose = [refX, refZ, refP]

			jointID = originNode.getJointID()

			#print "jointID =", jointID, "newJointID =", newJointID

			pose = self.probe.getJointWRTJointPose(initPose, jointID, newJointID) 
			#print initPose,pose

		else:
			# no active node exists and not the first,
			# so lets use the latest dormant one
			recentNode = self.refPoints[self.numRef-1]	

			refX, refZ, refP = recentNode.getRefPose()
			initPose = [refX, refZ, refP]

			jointID = recentNode.getJointID()
			pose = self.probe.getJointWRTJointPose(initPose, jointID, newJointID) 
			self.lastNode = recentNode

			isNewRoot = True

		# create the new node
		nodeID = self.numRef
		newNode = RefNode(nodeID, newJointID, pose[0], pose[1], pose[2], self.probe)	
		newNode.setNewRoot(isNewRoot)
		newNode.setGuideNode(self.lastNode.getNodeID())
		gndPose = self.probe.getActualJointPose(newJointID)
		newNode.setGroundTruthPose(gndPose[0],gndPose[1],gndPose[2])

		self.activeRefPtr[newJointID] = newNode	

		" draw the estimated reference nodes "
		if False:
			i = self.numRef
			## Create the visual reprsentation of active reference nodes
			name = "est_node" + str(i)
			entity = self.probe._mgr.createEntity(name, "Cube.mesh")
			node = self.probe._mgr.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)

			size = ogre.Vector3(0.05,0.05,0.05)
			node.setScale(size)

			position = ogre.Vector3(pose[0],0.1,pose[1])
			node.setPosition(position)
			# negative the angle because the bodies are somehow flipped upside down
			node.setOrientation(ogre.Quaternion(-pose[2],ogre.Vector3().UNIT_Y))

			entity.setCastShadows(False)
			entity.setMaterialName("Green")

			#entity.setVisible(False)
			self._allRefNodes.append(node)
			self._allRefEnt.append(entity)

		# record it to the array so we can access it later
		self.refPoints.append(newNode)	
		self.numRef += 1
		
		return newNode


	def printGraph(self):

		self.graphCount += 1

		fileName1 = "frame/refnodes_info_%04d.txt" % self.graphCount
		fileName2 = "frame/snake_pipe_%04d.graph" % self.graphCount
		fileName3 = "frame/snake_pipe_gnd_%04d.graph" % self.graphCount
		fileName4 = "frame/newRoots_%04d.txt" % self.graphCount

		# output max errors for all references nodes:

		refNodes = open(fileName1, 'w')
		for i in range(0,self.numRef):
			nodeID = self.refPoints[i].getNodeID()
			jointID = self.refPoints[i].getJointID()
			maxError = self.refPoints[i].getMaxStabilityError()
			avgError = self.refPoints[i].getAvgStabilityError()
			refNodes.write(str(nodeID) + " ")
			refNodes.write(str(jointID) + " ")
			refNodes.write(str(maxError) + " ")
			refNodes.write(str(avgError) + "\n")

		refNodes.close()
		
		toroGraph = open(fileName2, 'w')

		# print the nodes
		for i in range(self.numRef):
			nodeID = self.refPoints[i].getNodeID()
			jointID = self.refPoints[i].getJointID()
			x, z, p = self.refPoints[i].getRefPose()

			toroGraph.write("VERTEX " + str(nodeID) + " ")
			toroGraph.write(str(x) + " " + str(z) + " " + str(p) + "\n")

		# print the edges
		for i in range(self.numEdge):
			oID = self.refEdges[i].origin.getNodeID()
			tID = self.refEdges[i].target.getNodeID()
			x, z, p = self.refEdges[i].getRefConstraint()
			var = self.refEdges[i].variance
		
			toroGraph.write("EDGE ")
			toroGraph.write(str(oID) + " ")
			toroGraph.write(str(tID) + " " + str(x) + " " + str(z) + " " + str(p) + " ")
			toroGraph.write(str(var) + " " + str(var) + " " + str(var) + " " + str(var) + " " + str(var) + " " + str(var) + "\n")

		toroGraph.close()

		# now print the ground truth map
		gndGraph = open(fileName3, 'w')

		# print the vertices only, no edges for GND truth
		for i in range(self.numRef):

			nodeID = self.refPoints[i].getNodeID()
			jointID = self.refPoints[i].getJointID()
			x, z, p = self.refPoints[i].getGroundTruthPose()

			gndGraph.write("VERTEX " + str(nodeID) + " ")
			gndGraph.write(str(x) + " " + str(z) + " " + str(p) + "\n")

		gndGraph.close()

		# record the root nodes or orphan nodes, the breaks in the graph
		rootNodes = open(fileName4, 'w')
		for i in range(self.numRef):
			if self.refPoints[i].isNewRoot():
				nodeID = self.refPoints[i].getNodeID()
				rootNodes.write(str(nodeID) + " " + str(self.refPoints[i].getGuideNode()) + "\n")

		rootNodes.close()


	def computeActiveError(self):

		totalErrorX = 0.0
		totalErrorZ = 0.0
		totalErrorP = 0.0

		for i in range(self.numJoints):
			if self.activeRef[i]:
				activeNode = self.activeRefPtr[i]
				activeNode.resetPoseError()
				activeNode.updateTime()

		for i in range(self.numJoints):
			if self.activeRef[i]:

				activeNode = self.activeRefPtr[i]
				aID = activeNode.getNodeID()

				# compute the actual offset
				numEdges = activeNode.getNumEdges()

				# cycle through each edge
				# 1. check if this node is origin or target
				# 2. determine if edge is active
				# 3. determine if we haven't already found this edge

				for j in range(numEdges):
					origin = activeNode.getEdge(j).origin
					target = activeNode.getEdge(j).target

					oID = origin.getNodeID()
					tID = target.getNodeID()
					other = 0
					otherID = 0

					# origin is this node
					if aID == oID:
						otherID = tID	
						other = target

					# target is this node
					else:
						otherID = oID	
						other = origin

					# determine if other joint is active, the current active node, and a joint we haven't come to yet
					otherJoint = other.getJointID()
					otherActive = False

					if self.activeRef[otherJoint] and self.activeRefPtr[otherJoint].getNodeID() == otherID and otherJoint > i:
						otherActive = True

					# compute error of previous offset to current
					if otherActive:

						oX, oZ, oP = activeNode.getEdge(j).getRefConstraint()

						initPose = [0.0,0.0,0.0]
						targetPose = [0.,0.,0.]

						targetPose = self.probe.getJointWRTJointPose(initPose, origin.getJointID(), target.getJointID()) 

						errorX = targetPose[0] - oX
						errorZ = targetPose[1] - oZ
						errorP = targetPose[2] - oP

						totalErrorX += abs(errorX)
						totalErrorZ += abs(errorZ)
						totalErrorP += abs(errorP)

						origin.addPoseError(errorP)	
						target.addPoseError(errorP)
						
		return totalErrorP


	def isJointActive(self, jointID):
		if jointID >= self.numJoints:
			return False

		return self.activeRef[jointID]

	def getActiveJointPose(self, pose, jointID):

		if jointID >= self.numJoints or not self.activeRef[jointID]:
			raise

		return copy(self.activeRefPtr[jointID].getRefPose())


	def getNumActiveRef(self):
		return self.numActiveRef


	def getPoseProfile(self):
		pose = PoseProfile(self, 19, inSim = True)	
		return pose

	def normalizeAngle(self, angle):
		while angle > math.pi:
			angle = angle - 2*math.pi

		while angle <= -math.pi:
			angle = angle + 2*math.pi

		return angle


	def orderInflectionPoints(self,points):
		# order these points starting from the back of the back point of the first pose at points[8]

		unordered = deepcopy(points)
		orderedPoints = []

		origin = unordered.pop(8)
		orderedPoints.append(origin)

		# 1. create a new list of ordered points
		# 2. starting from the origin, find the nearest point and add it next order
		# 3. set nearest point to new origin
		# 4. goto 2

		while len(unordered) > 0:

			# find nearest point to origin
			minIndex = 0
			minDist = 10e100

			for i in range(len(unordered)):
				dist = math.sqrt((unordered[i][0]-origin[0])**2 + (unordered[i][1]-origin[1])**2)
				if dist < minDist:
					minDist = dist
					minIndex = i

			origin = unordered.pop(minIndex)
			orderedPoints.append(origin)

		return orderedPoints










