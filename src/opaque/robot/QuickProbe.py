
#import os
#import sys
#dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
#if not dir in sys.path:
#	sys.path.append(dir)

#from common import *
from math import *
from copy import *
from transform import Transform


class QuickProbe:

	def __init__(self, numSegs, segLength, segWidth, maxTorque):
		
		self.timer = 0.000

		self.transform = Transform()

		# control program for the snake probe
		self.control = 0
		self.wallEnv = 0

		self.numSegs = numSegs
		self.segLength = segLength
		self.segWidth = segWidth
		self.maxTorque = maxTorque

		self.robotParam = {}
		self.robotParam['numJoints'] = self.numSegs-1
		self.robotParam['numSegs'] = self.numSegs
		self.robotParam['segLength'] = self.segLength
		self.robotParam['segWidth'] = self.segWidth
		self.robotParam['maxTorque'] = self.maxTorque

		self.isStateChanged = True
		self.isJointChanged = True
		self.probeState = {}

		self.cmdJoints = [0.0 for i in range(self.getNumJoints())]
		self.joints = [0.0 for i in range(self.getNumJoints())]
		self.torques = [self.maxTorque for i in range(self.getNumJoints())]
		
	def __del(self):
		del self.control

	def setAnchor(self, isAnchor):
		pass

	def setWalls(self, wallEnv):
		self.wallEnv = wallEnv
	
	def createWall(self, points):
		self.wallEnv.createWall(points)
		
	def getWalls(self):
		if self.wallEnv == 0:
			return []
		else:
			return self.wallEnv.getWalls()
	
	def getProbeState(self):
		if self.isStateChanged:
			self.isStateChanged = False
	
			self.probeState = {}
			self.probeState['joints'] = copy(self.joints)
			self.probeState['cmdJoints'] = copy(self.cmdJoints)
			self.probeState['torques'] = copy(self.torques)
			
			errors = []
			for i in range(self.getNumJoints()):
				errors.append(self.getError(i))
	
			self.probeState['errors'] = errors
		
		return self.probeState	
		
	def getNumJoints(self):
		return self.numSegs - 1

	def addControl(self, controlProgram):
		self.control = controlProgram

	def frameStarted(self):

		self.timer += 0.001

		self.isStateChanged = True
		self.isJointChanged = True
		self.probeState = {}

		if self.control:
			self.control.frameStarted()

		for i in range(self.getNumJoints()):
			self.joints[i] = self.cmdJoints[i]
			self.transform.setJoints(self.joints)

	def normalizeAngle(self,angle):
		# this function converts the angle to its equivalent
	# in the range [-pi,pi] 

		while angle>pi:
			angle=angle-2*pi

		while angle<=-pi:
			angle=angle+2*pi

		return angle

	def getError(self, index):
		return abs(self.getServoCmd(index) - self.getServo(index))

	def getServo(self, index):
		if index < self.numSegs :
			return self.joints[index]
		raise

	def getServoCmd(self, index):
		if index < self.numSegs:
			return self.cmdJoints[index]
		raise

	def setJointTorque(self, i, torque):
		self.torques[i] = torque

	def getJointTorque(self, i):
		return self.torques[i]

	def setServo(self, i, angle):
		self.cmdJoints[i] = angle * pi / 180.0

	# gets pose of joint affixed to the i'th segment, should it be the i+1'th segment instead?
	def getActualJointPose(self, i):
		return [0.1, 0.0, 0.0]

	def getJointPose(self, originPose, originJoint, targetJoint):
		#return self.transform.getCJointPose(originPose,originJoint, targetJoint)
		return self.transform.getJointFromJoint(originPose,originJoint, targetJoint)
		
		
	def computeTransform(self, originJoint, targetJoint):

		" We exploit symmetry here.  If the normal calls are not NxN and are 1xN instead, remove symmetric computation "

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			return

		# forward kinematics
		if targetJoint > originJoint:

			if self.jointTransforms[originJoint][targetJoint] == None:

				" combine [originJoint, targetJoint-1] with [targetJoint-1, targetJoint] "
				if targetJoint-originJoint == 1:
					" compute it right now "

					totalAngle = self.joints[targetJoint]
					xTotal = self.segLength
					yTotal = 0.0

					offset = [xTotal, yTotal, totalAngle]
					offsetT = [-xTotal*cos(totalAngle), -xTotal*sin(totalAngle), totalAngle]

					" store result back into the dynamic programming table "
					self.jointTransforms[originJoint][targetJoint] = offset
					self.jointTransforms[targetJoint][originJoint] = offsetT
					#print originJoint,targetJoint,": ", offset
				
				else:
					
					" take two components of the dynamic programming table and combine them together for the solution "
					" store the result back in the table so it doesn't need to be recomputed "
					self.computeTransform(originJoint, targetJoint-1)
					offset1 = self.jointTransforms[originJoint][targetJoint-1]
					
					if self.jointTransforms[targetJoint-1][targetJoint] == None:
						
						totalAngle = self.joints[targetJoint]
						xTotal = self.segLength
						yTotal = 0.0
	
						offset2 = [xTotal, yTotal, totalAngle]
						offsetT = [-xTotal*cos(totalAngle), -xTotal*sin(totalAngle), totalAngle]
	
						" store result back into the dynamic programming table "
						self.jointTransforms[targetJoint-1][targetJoint] = offset2
						self.jointTransforms[targetJoint][targetJoint-1] = offsetT
						#offset2 = self.jointTransforms[targetJoint-1][targetJoint]
						#print targetJoint-1,targetJoint,": ", offset2
					
					else:
						#print "case2"
						offset2 = self.jointTransforms[targetJoint-1][targetJoint]
					
					
					x1 = offset1[0]
					y1 = offset1[1]
					theta_1 = offset1[2]
					
					" rotate x2,y2 by theta_1 "
					x2 = cos(theta_1) * offset2[0] + sin(theta_1) * offset2[1]
					y2 = -sin(theta_1) * offset2[0] + cos(theta_1) * offset2[1]
					theta_2 = offset2[2]
					
					" add x1+x2`, y1+y2`, theta_1 + theta_2 "
					xTotal = x1+x2
					yTotal = y1+y2

					angle = theta_1+theta_2
					while angle>pi:
						angle=angle-2*pi
					while angle<=-pi:
						angle=angle+2*pi
						
					finalOffset = [xTotal, yTotal, angle]
					offsetT = [-xTotal*cos(angle) + yTotal*sin(angle), -xTotal*sin(angle) - yTotal*cos(angle), angle]
				
					" store result back into the dynamic programming table "
					self.jointTransforms[originJoint][targetJoint] = finalOffset
					self.jointTransforms[targetJoint][originJoint] = offsetT
					#print originJoint,targetJoint,": ", finalOffset, "from", x1, x2, y1, y2, theta_1, theta_2
	
		# backward kinematics
		else:

			if self.jointTransforms[originJoint][targetJoint] == None:
				" originJoint > targetJoint "

				" combine [originJoint, targetJoint+1] with [targetJoint+1, targetJoint] "
				if originJoint-targetJoint == 1:
					" compute it right now "

					totalAngle = self.joints[originJoint]
					xTotal = -self.segLength*cos(totalAngle)
					yTotal = -self.segLength*sin(totalAngle)

					offset = [xTotal, yTotal, totalAngle]
					offsetT = [-xTotal*cos(totalAngle) - yTotal*sin(totalAngle), xTotal*sin(totalAngle) - yTotal*cos(totalAngle), totalAngle]

					" store result back into the dynamic programming table "
					self.jointTransforms[originJoint][targetJoint] = offset
					self.jointTransforms[targetJoint][originJoint] = offsetT
					#print originJoint,targetJoint,": ", offset
				
				else:
					
					" take two components of the dynamic programming table and combine them together for the solution "
					" store the result back in the table so it doesn't need to be recomputed "
					self.computeTransform(originJoint, targetJoint+1)
					offset1 = self.jointTransforms[originJoint][targetJoint+1]
					
					if self.jointTransforms[targetJoint+1][targetJoint] == None:
						
						totalAngle = self.joints[targetJoint+1]
						xTotal = -self.segLength*cos(totalAngle)
						yTotal = -self.segLength*sin(totalAngle)
	
						offset2 = [xTotal, yTotal, totalAngle]
						offsetT = [-xTotal*cos(totalAngle) - yTotal*sin(totalAngle), xTotal*sin(totalAngle) - yTotal*cos(totalAngle), totalAngle]
	
						" store result back into the dynamic programming table "
						self.jointTransforms[targetJoint+1][targetJoint] = offset2
						self.jointTransforms[targetJoint][targetJoint+1] = offsetT
						#offset2 = self.jointTransforms[targetJoint-1][targetJoint]
						#print targetJoint+1,targetJoint,": ", offset2
					
					else:
						#print "case2"
						offset2 = self.jointTransforms[targetJoint+1][targetJoint]
					
					
					x1 = offset1[0]
					y1 = offset1[1]
					theta_1 = offset1[2]
					
					" rotate x2,y2 by theta_1 "
					x2 = cos(theta_1) * offset2[0] - sin(theta_1) * offset2[1]
					y2 = sin(theta_1) * offset2[0] + cos(theta_1) * offset2[1]
					theta_2 = offset2[2]
					
					" add x1+x2`, y1+y2`, theta_1 + theta_2 "
					xTotal = x1+x2
					yTotal = y1+y2

					angle = theta_1+theta_2
					while angle>pi:
						angle=angle-2*pi
					while angle<=-pi:
						angle=angle+2*pi
						
					finalOffset = [xTotal, yTotal, angle]
					offsetT = [-xTotal*cos(angle) - yTotal*sin(angle), xTotal*sin(angle) - yTotal*cos(angle), angle]
					
					" store result back into the dynamic programming table "
					self.jointTransforms[originJoint][targetJoint] = finalOffset
					self.jointTransforms[targetJoint][originJoint] = offsetT
					#print originJoint,targetJoint,": ", finalOffset, "from", x1, x2, y1, y2, theta_1, theta_2
	

	def getJointFromJoint(self, originPose, originJoint, targetJoint):

		if self.isJointChanged:
			self.isJointChanged = False
			self.jointTransforms = [[None for j in range(0,39)] for i in range(0,39)]

		targetPose = [0.0,0.0,0.0]

		joints = self.probeState['joints']
		self.joints = joints

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			targetPose = copy(originPose)
			return targetPose

		if self.jointTransforms[originJoint][targetJoint] == None:
			self.computeTransform(originJoint,targetJoint)

			" force us to use forward kinematics which takes advantage of more zeros "
			
			" FIXME:  using this snippet will double the computation time, not sure why "
			#if originJoint > targetJoint:
			#	self.computeTransform(targetJoint,originJoint)
			#else:
			#	self.computeTransform(originJoint,targetJoint)
				
		" get transform from dynamic transform table "
		offset = self.jointTransforms[originJoint][targetJoint]
			
		# forward kinematics
		if targetJoint > originJoint:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2]

			" rotate base offset by current angle "
			xOff = cos(totalAngle) * offset[0] - sin(totalAngle) * offset[1]
			yOff = sin(totalAngle) * offset[0] + cos(totalAngle) * offset[1]			

			xTotal = xTotal + xOff
			zTotal = zTotal + yOff
			if targetJoint != self.numSegs-1:
				totalAngle = totalAngle - offset[2]

			angle = totalAngle
			while angle>pi:
				angle=angle-2*pi	
			while angle<=-pi:
				angle=angle+2*pi
			
			targetPose = [xTotal, zTotal, angle] 

			return targetPose

		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2]
			
			" rotate base offset by current angle "
			xOff = cos(totalAngle) * offset[0] - sin(totalAngle) * offset[1]
			yOff = sin(totalAngle) * offset[0] + cos(totalAngle) * offset[1]			

			xTotal = xTotal + xOff
			zTotal = zTotal + yOff
			if targetJoint != self.numSegs-1:
				totalAngle = totalAngle + offset[2]

			angle = totalAngle
			while angle>pi:
				angle=angle-2*pi	
			while angle<=-pi:
				angle=angle+2*pi
			
			targetPose = [xTotal, zTotal, angle] 
			#print "transforming", originPose, "by", offset, "resulting in", targetPose

			return targetPose
						
	def getJointWRTJointPose(self, originPose, originJoint, targetJoint):

		targetPose = [0.0,0.0,0.0]

		if originJoint >= self.numSegs-1 or targetJoint >= self.numSegs or originJoint < 0 or targetJoint < -1:
			print "ERROR: getJointWRTJointPose joint out of range!" 
			print "Received joint", originJoint , "and" , targetJoint , "with pose"
			print originPose[0], originPose[1], originPose[2]
			raise

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			targetPose = copy(originPose)
			return targetPose

		# forward kinematics
		if targetJoint > originJoint:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2]

			for i in range(originJoint+1, targetJoint+1):
				xTotal = xTotal + self.segLength*cos(totalAngle)
				zTotal = zTotal + self.segLength*sin(totalAngle)
				if i != self.numSegs-1:
					totalAngle = totalAngle - self.getServo(i)

			totalAngle = self.normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2] 

			ind = range(targetJoint+1, originJoint + 1) # (28, 11) 
			ind.reverse()

			for i in ind:
				totalAngle = totalAngle + self.getServo(i)
				xTotal = xTotal - self.segLength*cos(totalAngle)
				zTotal = zTotal - self.segLength*sin(totalAngle)

			totalAngle = self.normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		return targetPose
