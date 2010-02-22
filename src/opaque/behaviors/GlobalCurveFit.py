import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

import math
from common import *
from Behavior import *
from numpy import arange


globalNodeCount = 0

class GlobalCurveFit(Behavior):

	def __init__(self, probe, contacts, curve = 0):
		global globalNodeCount

		Behavior.__init__(self, probe)

		self.contacts = contacts
		self.numSegs = probe.numSegs
		self.setCurve(curve)
		self.startNode = 19
		self.endNode = 0
		
		# time divider
		self.setTimerAliasing(5)

		# nodes for drawing purposes
		self.nodeCount = globalNodeCount
		globalNodeCount += 1
		self.childNodes = []
		self.childEntities = []
		self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("globalCurveRoot" + str(self.nodeCount))

		

		
		self.vecNode1 = self.parentNode.createChildSceneNode("vecNode1")
		self.vecNode2 = self.parentNode.createChildSceneNode("vecNode2")
		self.vecNode3 = self.parentNode.createChildSceneNode("vecNode3")
		self.vecEnt1 = self.probe._mgr.createEntity("vecEnt1", "Cube.mesh")
		self.vecEnt2 = self.probe._mgr.createEntity("vecEnt2", "Cube.mesh")
		self.vecEnt3 = self.probe._mgr.createEntity("vecEnt3", "Cube.mesh")
		self.vecEnt1.setCastShadows(False)
		self.vecEnt2.setCastShadows(False)
		self.vecEnt3.setCastShadows(False)
		self.vecEnt1.setMaterialName("Red")
		self.vecEnt2.setMaterialName("Green")
		self.vecEnt3.setMaterialName("Blue")

		self.childNodes = [self.vecNode1, self.vecNode2, self.vecNode3]
		self.childEntities = [self.vecEnt1, self.vecEnt2, self.vecEnt3]
				
		position = ogre.Vector3(0.0,0.0,0.0)
		self.vecNode1.setPosition(position)
		self.vecNode2.setPosition(position)
		self.vecNode3.setPosition(position)
		
		newQuat = ogre.Quaternion(-ogre.Radian(0.0), ogre.Vector3().UNIT_Y)
		self.vecNode1.setOrientation(newQuat)
		newQuat = ogre.Quaternion(-ogre.Radian(math.pi/2), ogre.Vector3().UNIT_Y)
		self.vecNode2.setOrientation(newQuat)
		
		size = ogre.Vector3(0.20,0.05,0.05)
		self.vecNode1.setScale(size)
		self.vecNode2.setScale(size)
		self.vecNode3.setScale(size)
		
		self.vecNode1.attachObject(self.vecEnt1)
		self.vecNode2.attachObject(self.vecEnt2)
		self.vecNode3.attachObject(self.vecEnt3)
		
	def __del__(self):

		# remove all children
		self.parentNode.removeAllChildren()

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []		
		
	def setCurve(self, curve):
		self.curve = curve
		self.computePathDirection()
	
	def clearDraw(self):
		
		# remove all children
		self.parentNode.removeAllChildren()

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []

	def draw(self):
		
		return

		# remove all children
		self.parentNode.removeAllChildren()

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []
		
		# draw the points in the simulation
		#for i in range(len(self.realPath)):
		for i in arange(0.0,1.0, 0.01):
			
			pnt = self.curve.getU(i)
	
			# take samples from the fitted curve
			#pnt = copy(self.realPath[i])
	
			childNode = self.parentNode.createChildSceneNode("globalCurvePoint" + str(self.nodeCount) + "_" + str(i))
			self.childNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("globalCurveEnt" + str(self.nodeCount) + "_" +str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)
	
			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			childNode.setPosition(position)
	
			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)
	
			childNode.attachObject(currEntity)
	
	def setBoundaries(self, startNode, endNode):
		self.startNode = startNode
		self.endNode = endNode
		
		print "setting boundaries:", self.startNode, self.endNode

	# perform a continuous transition to a continuous fit
	def step(self):
		Behavior.step(self)
		
		if not self.isStep():
			return False
		
		self.clearDraw()

		self.resetJoints()
		
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
		
		if self.startNode < 0 or self.startNode >= (self.numSegs-1) or self.endNode < 0 or self.endNode >= (self.numSegs-1):
			print "Nodes out of bounds for fitBody:", self.startNode, self.endNode
			raise

		backwards = not self.direction

		self.closePoints = []
		finalAngles = []
		
		# terminator point and orientation
		tPoint = self.curve.getU(1.0)
		tVec = self.curve.getUVector(1.0)
		tAngle = acos(tVec[0])
		if asin(tVec[1]) < 0:
			tAngle = -tAngle
			
		
		# now we need to determine which segments are active
		# start from the inner segment and work to the 0th segment
		if self.startNode == 0:
			#segments = range(0,self.endNode+1)
			segments = range(0,self.endNode+1)
			segments.reverse()
			segmentOrder = -1
			
		# start from inner segment and work out to the 39th segment
		elif self.endNode == self.numSegs-2:
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
		
		print "global curve fit segments:", segments
		
		#origin = self.probe.getActualSegPose(segments[0])
		origin = self.contacts.getAverageSegPose(segments[0])

		#origin = self.probe.getActualSegPose(segments[0])
		xTotal = origin[0]
		yTotal = origin[1]
		totalAngle = origin[2]

		if segmentOrder == -1:
			# adding constant offset so that we are on the joint
			xTotal = xTotal + (self.probe.segLength/2)*cos(totalAngle)
			yTotal = yTotal + (self.probe.segLength/2)*sin(totalAngle)
			totalAngle += -self.probe.getServo(segments[0])
		else:
			xTotal = xTotal - (self.probe.segLength/2)*cos(totalAngle)
			yTotal = yTotal - (self.probe.segLength/2)*sin(totalAngle)
			totalAngle += self.probe.getServo(segments[0]-1)
			
			
		totalVec = [cos(totalAngle),sin(totalAngle)]

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
				
				newQuat = ogre.Quaternion(-ogre.Radian(ang1), ogre.Vector3().UNIT_Y)
				self.vecNode1.setOrientation(newQuat)
				newQuat = ogre.Quaternion(-ogre.Radian(ang2), ogre.Vector3().UNIT_Y)
				self.vecNode2.setOrientation(newQuat)
	
				position = ogre.Vector3(xTotal + (self.probe.segLength/2)*cos(ang1),0.0,yTotal + (self.probe.segLength/2)*sin(ang1))
				self.vecNode1.setPosition(position)
				position = ogre.Vector3(xTotal + (self.probe.segLength/2)*cos(ang2),0.0,yTotal + (self.probe.segLength/2)*sin(ang2))
				self.vecNode2.setPosition(position)
			if i == 39:
				ang3 = acos(totalVec[0])
				if asin(totalVec[1]) < 0:
					ang3 = -ang3
				newQuat = ogre.Quaternion(-ogre.Radian(ang3), ogre.Vector3().UNIT_Y)
				self.vecNode3.setOrientation(newQuat)
				position = ogre.Vector3(xTotal + (self.probe.segLength/2)*cos(ang3),0.0,yTotal + (self.probe.segLength/2)*sin(ang3))
				self.vecNode3.setPosition(position)
			
				
			# weighted direction we want to point to given orientation and location vectors
			desVec = [(orientVec[0] + dist*locVec[0])/(dist+1.0), (orientVec[1] + dist*locVec[1])/(dist+1.0)]
			desMag = sqrt(desVec[0]**2 + desVec[1]**2)
			desVec = [desVec[0]/desMag, desVec[1]/desMag]
		

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
				angle = self.probe.getServoCmd(i)
			else:
				angle = self.probe.getServoCmd(i-1)
				
			#angle = self.probe.getServo(i)
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
				xTotal = xTotal - self.probe.segLength*cos(totalAngle)
				yTotal = yTotal - self.probe.segLength*sin(totalAngle)
			
			else:
				totalAngle += finalAngle
				xTotal = xTotal + self.probe.segLength*cos(totalAngle)
				yTotal = yTotal + self.probe.segLength*sin(totalAngle)
				
			childNode = self.parentNode.createChildSceneNode("globalCurvePoint" + str(i))
			self.childNodes.append(childNode)
	
			currEntity = self.probe._mgr.createEntity("globalCurveEnt" +str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Red")
			self.childEntities.append(currEntity)
	
			position = ogre.Vector3(xTotal,0.0,yTotal)
			childNode.setPosition(position)
	
			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)
	
			childNode.attachObject(currEntity)
				
		"""
		if not self.direction:
			temp = self.startNode
			self.startNode = self.endNode
			self.endNode = temp
			
		origin = self.probe.getActualJointPose(self.startNode)
		xTotal = origin[0]
		yTotal = origin[1]
		totalAngle = origin[2]
		totalVec = [cos(totalAngle),sin(totalAngle)]


		if not self.direction:

			# joints we're visiting in reverse order
			if self.startNode  >= 35:
				joints = range(self.endNode, 34+1)
				joints.reverse()

				tEnd = self.endNode
				if tEnd < 35:
					tEnd = 35

				zeroJoints = range(tEnd,self.startNode+1)
				for i in zeroJoints:
					#print "setting", i, "to zero"
					#self.probe.setServo(i, 0.0)
					#self.joints[i] = 0.0
					self.setJoint(i, 0.0)

			else:
				joints = range(self.endNode, self.startNode+1)
				joints.reverse()
		"""

		return False

	# perform a continuous transition to a continuous fit
	def old_step(self):
		Behavior.step(self)

		# FIXME:  if the tail joint is controlled and is also the start node,
		# we need to reverse the control because normal control will cause the 
		# tail to curl into a ball

		# cases:
		# 1. startNode == endNode => do nothing
		# 2. startNode or endNode out of bounds => raise exception
		# 3. startNode > endNode => backwards fitting
		# 4. startNode < endNode => forwards fitting

		#if startNode == endNode:
		#	return
		
		if not self.isStep():
			return False
		
		self.resetJoints()

		if self.startNode < 0 or self.startNode >= (self.numSegs-1) or self.endNode < 0 or self.endNode >= (self.numSegs-1):
			print "Nodes out of bounds for fitBody:", self.startNode, self.endNode
			raise
		
		if self.startNode == self.endNode:
			return False
		

		if self.startNode > self.endNode:
			backwards = True
		else:
			backwards = False

		origin = self.contacts.getAveragePose(self.startNode)
		#origin = self.probe.getActualJointPose(self.startNode)
		xTotal = origin[0]
		yTotal = origin[1]
		totalAngle = origin[2]
		totalVec = [cos(totalAngle),sin(totalAngle)]

		self.closePoints = []
		finalAngles = []

		# 1.  terminator point and orientation
		tPoint = self.curve.getU(1.0)
		tVec = self.curve.getUVector(1.0)
		tAngle = acos(tVec[0])
		if asin(tVec[1]) < 0:
			tAngle = -tAngle

		if backwards:

			# joints we're visiting in reverse order
			if self.startNode  >= 35:
				joints = range(self.endNode, 34+1)
				joints.reverse()

				tEnd = self.endNode
				if tEnd < 35:
					tEnd = 35

				zeroJoints = range(tEnd,self.startNode+1)
				for i in zeroJoints:
					#print "setting", i, "to zero"
					#self.probe.setServo(i, 0.0)
					#self.joints[i] = 0.0
					self.setJoint(i, 0.0)

			else:
				joints = range(self.endNode, self.startNode+1)
				joints.reverse()

			# compute the change in joint angle for every consecutive joint
			for i in joints:

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
				orientVec[0] *= -1
				orientVec[1] *= -1
				locVec[0] *= -1
				locVec[1] *= -1

				# weighted direction we want to point to given orientation and location vectors
				desVec = [(orientVec[0] + dist*locVec[0])/(dist+1.0), (orientVec[1] + dist*locVec[1])/(dist+1.0)]
				desMag = sqrt(desVec[0]**2 + desVec[1]**2)
				desVec = [desVec[0]/desMag, desVec[1]/desMag]

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
				angle = self.probe.getServoCmd(i)
				#angle = self.probe.getServo(i)
				finalAngle = k*desAngle + (1-k)*angle

				#print "setting", i, "to", finalAngle*180.0/pi
				#self.probe.setServo(i, finalAngle*180.0/pi)
				self.setJoint(i, finalAngle * 180.0 / pi)
				finalAngles.append(finalAngle)

				totalVec = [cos(finalAngle), sin(finalAngle)]
				totalVec = [totalVec[0]*cos(rotAngle) - totalVec[1]*sin(rotAngle), totalVec[0]*sin(rotAngle) + totalVec[1]*cos(rotAngle)]

				self.closePoints.append(copy([xTotal,yTotal]))
				totalAngle += finalAngle
				xTotal = xTotal - self.probe.segLength*cos(totalAngle)
				yTotal = yTotal - self.probe.segLength*sin(totalAngle)

		# forward fitting
		# NOTE:  in this case we are modifying the current frame instead of modifying the next frame
		else:

			# self.startNode < self.endNode

			"""
			if self.startNode  <= 4:
				joints = range(5, self.endNode+1)

				tEnd = self.endNode
				if tEnd > 4:
					tEnd = 4

				zeroJoints = range(self.startNode,tEnd+1)
				for i in zeroJoints:
					#self.probe.setServo(i, 0.0)
					self.setJoint(i,0.0)

			else:
				joints = range(self.startNode, self.endNode+1)
			"""
				
			joints = range(self.startNode, self.endNode+1)

			for i in joints:

				dist, linePoint, orientVec = self.curve.findClosestPoint([xTotal, yTotal])

				# check if we've gone beyond the boundaries of the curve
				relPoint = [xTotal-tPoint[0], yTotal-tPoint[1]]
				rotPoint = [relPoint[0]*cos(tAngle) + relPoint[1]*sin(tAngle), -relPoint[0]*sin(tAngle) + relPoint[1]*cos(tAngle)]
				orthoDist = fabs(rotPoint[1])
				if rotPoint[0] > 0 and orthoDist < dist:

					closePoint = [rotPoint[0], 0.0]
					closePoint = [closePoint[0]*cos(tAngle) - closePoint[1]*sin(tAngle), closePoint[0]*sin(tAngle) + closePoint[1]*cos(tAngle)]
					closePoint = [closePoint[0] + tPoint[0], closePoint[1] + tPoint[1]]


					orientVec = copy(tVec) 
					locVec = [closePoint[0]-xTotal, closePoint[1]-yTotal]
					locMag = sqrt(locVec[0]**2 + locVec[1]**2)
					locVec = [locVec[0]/locMag, locVec[1]/locMag]

					dist = orthoDist
				else:

					# compute the direction from current joint to line point
					locVec = [linePoint[0]-xTotal, linePoint[1]-yTotal]
					locMag = sqrt(locVec[0]**2 + locVec[1]**2)
					locVec = [locVec[0]/locMag, locVec[1]/locMag]

				# weighted direction we want to point to
				desVec = [(orientVec[0] + dist*locVec[0])/(dist+1.0), (orientVec[1] + dist*locVec[1])/(dist+1.0)]
				desMag = sqrt(desVec[0]**2 + desVec[1]**2)
				desVec = [desVec[0]/desMag, desVec[1]/desMag]

				# Compute the rotation to desVec from totalVec and this is our desired joint angle

				# 1. rotate totalVec to [1,0]
				rotAngle = acos(totalVec[0])
				if asin(totalVec[1]) < 0:
					rotAngle = -rotAngle

				# 2. rotate desVec by the same rotation as totalVec
				desRotVec = [desVec[0]*cos(rotAngle) + desVec[1]*sin(rotAngle), -desVec[0]*sin(rotAngle) + desVec[1]*cos(rotAngle)]
				desMag = sqrt(desRotVec[0]**2 + desRotVec[1]**2)
				desRotVec[0] /= desMag
				desRotVec[1] /= desMag

				# 3. compute angle of desRotVec and this is our desired offset from the current joint angle
				offsetAngle = acos(desRotVec[0])
				if asin(desRotVec[1]) < 0:
					offsetAngle = -offsetAngle

				# offsetAngle must be negative because we are going forwards
				desAngle = -offsetAngle + self.probe.getServoCmd(i)
				desAngle = normalizeAngle(desAngle)

				k = 0.1
				angle = self.probe.getServoCmd(i)
				#angle = self.probe.getServo(i)
				finalAngle = k*desAngle + (1-k)*angle
				diffAngle = finalAngle-angle

				#self.probe.setServo(i, finalAngle*180.0/pi)
				self.setJoint(i, finalAngle * 180.0 / pi)
				finalAngles.append(finalAngle)

				if i < self.endNode:
					rotAngle -= diffAngle
					nextAngle = self.probe.getServoCmd(i+1)
					totalVec = [cos(-nextAngle), sin(-nextAngle)]
					totalVec = [totalVec[0]*cos(rotAngle) - totalVec[1]*sin(rotAngle), totalVec[0]*sin(rotAngle) + totalVec[1]*cos(rotAngle)]

					self.closePoints.append(copy([xTotal,yTotal]))
					totalAngle -= diffAngle
					xTotal = xTotal + self.probe.segLength*cos(totalAngle)
					yTotal = yTotal + self.probe.segLength*sin(totalAngle)
					totalAngle -= nextAngle


			#print finalAngles

		return False
	
	def getPathDirection(self):
		return self.direction
	
	def computePathDirection(self):
		joints = range(NUM_SEGS-1)

		# segment angle vector
		vecSum1 = [0.0,0.0]

		# path vector
		vecSum2 = [0.0,0.0]

		for j in joints:

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

		print "globalCurveFit direction =", self.direction
