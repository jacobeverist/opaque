
#import os
#import sys
#dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
#if not dir in sys.path:
#	sys.path.append(dir)

#from common import *
from math import *
from copy import *

class QuickProbe:

	def __init__(self, numSegs, segLength, segWidth, maxTorque):
		
		self.timer = 0.000

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
		probeState = {}
		probeState['joints'] = copy(self.joints)
		probeState['cmdJoints'] = copy(self.cmdJoints)
		probeState['torques'] = copy(self.torques)
		
		errors = []
		for i in range(self.getNumJoints()):
			errors.append(self.getError(i))

		probeState['errors'] = errors
		
		return probeState	
		
	def getNumJoints(self):
		return self.numSegs - 1

	def addControl(self, controlProgram):
		self.control = controlProgram

	def frameStarted(self):

		self.timer += 0.001

		if self.control:
			self.control.frameStarted()

		for i in range(self.getNumJoints()):
			self.joints[i] = self.cmdJoints[i]

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
