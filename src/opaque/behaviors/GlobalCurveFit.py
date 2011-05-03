
from functions import Pose
from math import sin, cos, asin, acos, sqrt, fabs, pi
from Behavior import Behavior
from numpy import arange
from copy import copy

globalNodeCount = 0

class GlobalCurveFit(Behavior):

	def __init__(self, robotParam, contacts, curve = 0, localNode = 0):
		global globalNodeCount

		Behavior.__init__(self, robotParam)


		#print "setting currNode =", localNode
		self.currNode = localNode
		self.contacts = contacts
		self.numSegs = robotParam['numSegs']
		self.numJoints = robotParam['numJoints']
		self.segLength = robotParam['segLength']
		self.setCurve(curve)
		self.startNode = 19
		self.endNode = 0
		
		# time divider
		#self.setTimerAliasing(5)
		#self.currNode = localNode
		# nodes for drawing purposes
		self.nodeCount = globalNodeCount
		globalNodeCount += 1
		self.childNodes = []
		self.childEntities = []
		self.childSegNodes = []
		self.childSegEntities = []
		
		self.parentNode = 0

		self._allRefNodes = []
		self._allRefEnt = []

	def __del__(self):

		self.clearDraw()
		#self.probe._mgr.destroySceneNode(self.parentNode)
				
	def setCurve(self, curve):
		self.curve = curve
		self.computePathDirection()
		
		print "globalCurveFit Behavior direction =", self.direction
	
	def clearDraw(self):
		

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []
		
		for child in self.childSegNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childSegNodes = []
	
		for child in self.childSegEntities:
			self.probe._mgr.destroyEntity(child)

		self.childSegEntities = []
		
		# remove all children
		if self.parentNode != 0:
			self.parentNode.removeAllChildren()


	def draw(self):
		
		if self.parentNode == 0:
			self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("globalCurveRoot" + str(self.nodeCount))

		self.clearDraw()

		if self.currNode != 0:
			#onSegPose = Pose(self.contacts.getClosestPose(self.currNode.rootNode))
			onSegPose = Pose(self.contacts.getAveragePose(self.currNode.rootNode))
			onSegActPose = Pose( self.probe.getActualJointPose(self.currNode.rootNode))
	
		# draw the points in the simulation
		#for i in range(len(self.realPath)):
		for i in arange(0.0,1.0, 0.01):
			
			pnt = self.curve.getU(i)

			" convert points from estimated points to actual points in view "
			if self.currNode != 0:
				newPnt = copy(pnt)
				
				localPnt = onSegPose.convertGlobalToLocal(pnt)
				pnt = onSegActPose.convertLocalToGlobal(localPnt)
				
				#localPnt = pose.convertGlobalToLocal(pnt)
				#pnt = pose2.convertLocalToGlobal(localPnt)
			
			
			# take samples from the fitted curve
			#pnt = copy(self.realPath[i])

			childNode = self.parentNode.createChildSceneNode("globalCurvePoint" + str(self.nodeCount) + "_" + str(i))
			self.childNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("globalCurveEnt" + str(self.nodeCount) + "_" +str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)

			#position = ogre.Vector3(pnt[0],0.0,pnt[1])
			#childNode.setPosition(position)
	
			#size = ogre.Vector3(0.03,0.03,0.03)
			#childNode.setScale(size)
	
			childNode.attachObject(currEntity)
	
	def setBoundaries(self, startNode, endNode):
		self.startNode = startNode
		self.endNode = endNode
		
		print "setting boundaries:", self.startNode, self.endNode

	# perform a continuous transition to a continuous fit
	def step(self, probeState):
		Behavior.step(self, probeState)
		
		if not self.isStep():
			return False
		
		#self.clearDraw()

		self.resetJoints()
		
		stateJoints = probeState['joints']
		cmdJoints = probeState['cmdJoints']
		
		# 2 pieces of information
		# 1. the direction of the curve
		# 2. the direction of kinematics dependent on which side of the snake is anchored
		
		# Three cases:
		# 1. front half active (0-2x)
		# 2. middle section active (1x-3x)
		# 3. back half active (2x-38)

		
		#print self.startNode, self.endNode
		
		if self.startNode == -1 or self.endNode == -1:
			return False
		
		if self.startNode < 0 or self.startNode >= (self.numJoints) or self.endNode < 0 or self.endNode >= (self.numJoints):
			print "Nodes out of bounds for fitBody:", self.startNode, self.endNode
			raise

		#self.direction = False

		backwards = not self.direction

		self.closePoints = []
		finalAngles = []
		
		# terminator point and orientation
		tPoint = self.curve.getU(1.0)
		tVec = self.curve.getUVector(0.9)
		tAngle = acos(tVec[0])
		if asin(tVec[1]) < 0:
			tAngle = -tAngle
		#print "tPoint:", tPoint, "tVec:", tVec, "tAngle:", tAngle
		
		#tVec = [cos(math.pi/2), sin(math.pi/2)]
		#tAngle = acos(tVec[0])
		#if asin(tVec[1]) < 0:
		#	tAngle = -tAngle
		#print "tPoint:", tPoint, "tVec:", tVec, "tAngle:", tAngle
		
		# now we need to determine which segments are active
		# start from the inner segment and work to the 0th segment
		if self.startNode == 0:
			#segments = range(0,self.endNode+1)
			segments = range(0,self.endNode+1)
			segments.reverse()
			segmentOrder = -1
			
		# start from inner segment and work out to the 39th segment
		elif self.endNode == self.numJoints-1:
			segments = range(self.startNode+1,self.endNode+2)
			segmentOrder = 1
		
		# use the direction of the curve to guide our kinematic order
		else:
			# something in the middle
			# includes one extra segment depending on fitting direction
			
			if not self.direction:
				segments = range(self.startNode,self.endNode+1)
				segments.reverse()
				segmentOrder = -1
			else:
				segments = range(self.startNode+1,self.endNode+2)
				segmentOrder = 1
		
		#print "global curve fit segments:", segments
		
		#origin = self.probe.getActualSegPose(segments[0])
		origin = self.contacts.getAverageSegPose(segments[0])
		#origin = self.contacts.getSegPose(segments[0])

		#origin = self.probe.getActualSegPose(segments[0])
		xTotal = origin[0]
		yTotal = origin[1]
		totalAngle = origin[2]

		if segmentOrder == -1:
			# adding constant offset so that we are on the joint
			xTotal = xTotal + (self.segLength/2)*cos(totalAngle)
			yTotal = yTotal + (self.segLength/2)*sin(totalAngle)
			totalAngle += -stateJoints[segments[0]]
		else:
			xTotal = xTotal - (self.segLength/2)*cos(totalAngle)
			yTotal = yTotal - (self.segLength/2)*sin(totalAngle)
			totalAngle += stateJoints[segments[0]-1]
			
			
		totalVec = [cos(totalAngle),sin(totalAngle)]

		"""
		if self.currNode != 0:
			onSegPose = Pose(self.contacts.getAveragePose(self.currNode.rootNode))
			onSegActPose = Pose( self.probe.getActualJointPose(self.currNode.rootNode))
		"""
		
		#self.draw()

		# compute the change in joint angle for every consecutive joint
		for i in segments:

			# find the closest point on the roadmap curve for the joint's location
			# - returns the distance, the point that is closest, the orientation of the curve
			dist, linePoint, orientVec = self.curve.findClosestPoint([xTotal, yTotal])

			# check if we've gone beyond the boundaries of the curve
			# by measuring the orthogonal distance of the current joint to the vector at U=1.0
			relPoint = [xTotal-tPoint[0], yTotal-tPoint[1]]
			rotPoint = [relPoint[0]*cos(tAngle) + relPoint[1]*sin(tAngle), -relPoint[0]*sin(tAngle) + relPoint[1]*cos(tAngle)]
			orthoDist = fabs(rotPoint[1])

			# if orthogonal point is in the positive direction of the Uvector
			# and less distance than the closest point already computed
			if rotPoint[0] > 0 and orthoDist < dist:

				# compute the closest point on the vector
				closePoint = [rotPoint[0], 0.0]
				closePoint = [closePoint[0]*cos(tAngle) - closePoint[1]*sin(tAngle), closePoint[0]*sin(tAngle) + closePoint[1]*cos(tAngle)]
				closePoint = [closePoint[0] + tPoint[0], closePoint[1] + tPoint[1]]

				# orientation vector is now the vector at U=1.0
				orientVec = copy(tVec) 

				# location vector is now computed to the orthogonal point
				locVec = [closePoint[0]-xTotal, closePoint[1]-yTotal]
				locMag = sqrt(locVec[0]**2 + locVec[1]**2)
				locVec = [locVec[0]/locMag, locVec[1]/locMag]

				dist = orthoDist

			# now treat everything regularly
			else:

				# closest point is point found on interior of roadmap curve

				# location vector is computed to point on interior of roadmap
				locVec = [linePoint[0]-xTotal, linePoint[1]-yTotal]
				locMag = sqrt(locVec[0]**2 + locVec[1]**2)
				locVec = [locVec[0]/locMag, locVec[1]/locMag]

			# negate the orientation and location vector since the backwards kinematics
			# require us to align the snake segments to 180 degrees from the desired orientation
			# negating the orientation vectors allows us to use the same algorithm for forwards and backwards
			
			if not self.direction:
				orientVec[0] *= -1
				orientVec[1] *= -1
				locVec[0] *= -1
				locVec[1] *= -1
			
			
				
			if i == 39:
				ang1 = acos(orientVec[0])
				if asin(orientVec[1]) < 0:
					ang1 = -ang1
				ang2 = acos(locVec[0])
				if asin(locVec[1]) < 0:
					ang2 = -ang2
				
				#newQuat = ogre.Quaternion(-ogre.Radian(ang1), ogre.Vector3().UNIT_Y)
				#self.vecNode1.setOrientation(newQuat)
				#newQuat = ogre.Quaternion(-ogre.Radian(ang2), ogre.Vector3().UNIT_Y)
				#self.vecNode2.setOrientation(newQuat)
	
				#position = ogre.Vector3(xTotal + (self.segLength/2)*cos(ang1),0.0,yTotal + (self.segLength/2)*sin(ang1))
				#self.vecNode1.setPosition(position)
				#position = ogre.Vector3(xTotal + (self.segLength/2)*cos(ang2),0.0,yTotal + (self.segLength/2)*sin(ang2))
				#self.vecNode2.setPosition(position)

			if i == 39:
				ang3 = acos(totalVec[0])
				if asin(totalVec[1]) < 0:
					ang3 = -ang3
				#newQuat = ogre.Quaternion(-ogre.Radian(ang3), ogre.Vector3().UNIT_Y)
				#self.vecNode3.setOrientation(newQuat)
				#position = ogre.Vector3(xTotal + (self.segLength/2)*cos(ang3),0.0,yTotal + (self.segLength/2)*sin(ang3))
				#self.vecNode3.setPosition(position)
			
				
			# weighted direction we want to point to given orientation and location vectors
			desVec = [(orientVec[0] + dist*locVec[0])/(dist+1.0), (orientVec[1] + dist*locVec[1])/(dist+1.0)]
			desMag = sqrt(desVec[0]**2 + desVec[1]**2)
			desVec = [desVec[0]/desMag, desVec[1]/desMag]

			#print self.direction, i, "orientVec:", orientVec, "locVec:", locVec, "desVec:", desVec

			# Compute the required rotation from current orientation to desired orientation (totalVec, desVec)

			# 1. rotate totalVec to [1,0]
			rotAngle = acos(totalVec[0])
			if asin(totalVec[1]) < 0:
				rotAngle = -rotAngle

			# 2. rotate desVec by the same rotation as totalVec
			desRotVec = [desVec[0]*cos(rotAngle) + desVec[1]*sin(rotAngle), -desVec[0]*sin(rotAngle) + desVec[1]*cos(rotAngle)]
			desMag = sqrt(desRotVec[0]**2 + desRotVec[1]**2)
			desRotVec[0] /= desMag
			desRotVec[1] /= desMag

			# 3. compute angle of desRotVec and this is our desired joint angle
			desAngle = acos(desRotVec[0])
			if asin(desRotVec[1]) < 0:
				desAngle = -desAngle

			k = 0.1
			if segmentOrder == -1:
				angle = cmdJoints[i]
			else:
				angle = cmdJoints[i-1]
				
			#angle = stateJoints[i]
			if segmentOrder == -1:
				finalAngle = k*desAngle + (1-k)*angle
			else:
				finalAngle = -k*desAngle + (1-k)*angle
				#finalAngle = k*desAngle + (1-k)*angle

			#print "setting", i, "to", finalAngle*180.0/pi
			#self.probe.setServo(i, finalAngle*180.0/pi)
			if segmentOrder == -1:
				self.setJoint(i, finalAngle * 180.0 / pi)
			else:
				self.setJoint(i-1, finalAngle * 180.0 / pi)
				
			# FIXME: this works, by I still don't understand why
			if segmentOrder == 1:
				finalAngle = -finalAngle
				
			finalAngles.append(finalAngle)

			totalVec = [cos(finalAngle), sin(finalAngle)]
			totalVec = [totalVec[0]*cos(rotAngle) - totalVec[1]*sin(rotAngle), totalVec[0]*sin(rotAngle) + totalVec[1]*cos(rotAngle)]

			
			self.closePoints.append(copy([xTotal,yTotal]))
			if segmentOrder == -1:
				totalAngle += finalAngle
				xTotal = xTotal - self.segLength*cos(totalAngle)
				yTotal = yTotal - self.segLength*sin(totalAngle)
			
			else:
				totalAngle += finalAngle
				xTotal = xTotal + self.segLength*cos(totalAngle)
				yTotal = yTotal + self.segLength*sin(totalAngle)


			"""
			#print "creating child node " + "globalSegPoint" + str(self.nodeCount) + "_" + str(i)
			childNode = self.parentNode.createChildSceneNode("globalSegPoint" + str(self.nodeCount) + "_" + str(i))
			self.childSegNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("globalSegEnt" + str(self.nodeCount) + "_" +str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Red")
			self.childSegEntities.append(currEntity)

			pnt = [xTotal,yTotal]
			if self.currNode != 0:
				localPnt = onSegPose.convertGlobalToLocal(pnt)
				pnt = onSegActPose.convertLocalToGlobal(localPnt)
	
				#localPnt = pose.convertGlobalToLocal(pnt)
				#pnt = pose2.convertLocalToGlobal(localPnt)


			if segmentOrder == -1:
				pnt[0] = pnt[0] + 0.5*self.segLength*cos(totalAngle)
				pnt[1] = pnt[1] + 0.5*self.segLength*sin(totalAngle)
			
			else:
				pnt[0] = pnt[0] - 0.5*self.segLength*cos(totalAngle)
				pnt[1] = pnt[1] - 0.5*self.segLength*sin(totalAngle)

			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			#position = ogre.Vector3(xTotal,0.0,yTotal)
			
			R = ogre.Quaternion(-ogre.Radian(totalAngle), ogre.Vector3().UNIT_Y)
			
			childNode.setPosition(position)
			childNode.setOrientation(R)
			
			size = ogre.Vector3(ARMLENGTH, 0.01, ARMWIDTH)
			childNode.setScale(size)
	
			childNode.attachObject(currEntity)
			"""

		return False

	def getPathDirection(self):
		print "returning GlobalCurveFit direction:", self.direction
		return self.direction
	
	def computePathDirection(self):
		joints = range(self.numJoints)

		# segment angle vector
		vecSum1 = [0.0,0.0]

		# path vector
		vecSum2 = [0.0,0.0]

		for j in joints:

			#origin = self.contacts.getClosestPose(j)
			origin = self.contacts.getAveragePose(j)
			#origin = self.probe.getActualJointPose(j)
			poseVec = [cos(origin[2]),sin(origin[2])] 
			vecSum1[0] += poseVec[0]
			vecSum1[1] += poseVec[1]

			dist, linePoint, orientVec = self.curve.findClosestPoint([origin[0], origin[1]])
			#print dist, linePoint, orientVec
			vecSum2[0] += orientVec[0]
			vecSum2[1] += orientVec[1]

		# renormalize the snake orientation vector
		mag1 = sqrt(vecSum1[0]**2 + vecSum1[1]**2)
		vecSum1 = [vecSum1[0]/mag1,vecSum1[1]/mag1]

		#print "snake vector:", vecSum1
		#print "path vector:", vecSum2
		# compute angle of snake orientation average
		ang1 = acos(vecSum1[0])
		if asin(vecSum1[1]) < 0:
			ang1 *= -1

		# rotate path vector by the same angle
		vecSum2 = [vecSum2[0]*cos(ang1) + vecSum2[1]*sin(ang1), -vecSum2[0]*sin(ang1) + vecSum2[1]*cos(ang1)]

		# renormalize the path orientation vector
		mag2 = sqrt(vecSum2[0]**2 + vecSum2[1]**2)
		vecSum2 = [vecSum2[0]/mag2,vecSum2[1]/mag2]
		
		# compute the angle of the path orientation average wrt snake orientation average
		ang2 = acos(vecSum2[0])
		if asin(vecSum2[1]) < 0:
			ang2 *= -1

		#print "snake orientation:", ang1
		#print "PATH DIRECTION:", ang2

		# forward, path along ascending joints
		if fabs(ang2) < pi/2:
			self.direction = True
		# backward, path along descending joints
		else:
			self.direction = False

		#self.direction = False

		#print "globalCurveFit direction =", self.direction
