import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from pose import ValueStability
from Transition import Transition

class Curl(Behavior):
	
	def __init__(self, probe):
		Behavior.__init__(self, probe)	

		# angle value for curling joints
		self.val = 0
		self.rail = 25.0
		self.up = True
		self.topJoint = 0

		self.poseWindow = []

		self.stabilityX = ValueStability(1e-7, 10)
		self.stabilityY = ValueStability(1e-7, 10)
		
		self.lastPose = [0.0,0.0,0.0]
		
		self.poseIndex = 2
		
		self.returnedSuccess = False
		self.slipError = 0.0
		
		self.direction = True

		self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = False
	
	def reset(self):
		
		self.val = 0
		self.rail = 25.0
		self.up = True
		self.topJoint = 0

		self.poseWindow = []

		self.stabilityX.reset()
		self.stabilityY.reset()
		
		self.lastPose = [0.0,0.0,0.0]
		
		self.poseIndex = 2
		
		self.returnedSuccess = False
		self.slipError = 0.0
		
		self.direction = True

		self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.timeInc = 1
		self.gaitTimer = 0
		
		self.resetJoints()
		#self._joints = [None for i in range(self.probe.numSegs-1)]
			
	def setDirection(self, val):
		self.direction = val
		#self.stabilityX.reset()
		#self.stabilityY.reset()
		self.resetJoints()
		
		if self.direction:
			for i in range(0, self.topJoint):
				# joint value
				self.setJoint(i, self.val)
		else:
			#print "setting", range(self.topJoint, self.probe.numSegs-1)
			for i in range(self.topJoint, self.probe.numSegs-1):
				# joint value
				self.setJoint(i, self.val)		
				
	def getSlipError(self):
		return self.slipError
	
	def setTopJoint(self, node):
		self.topJoint = node

	def step(self):
		Behavior.step(self)

		if not self.isStep():
			return False
		
		self.resetJoints()

		if self.isTransitioning:

			# first steps
			isDone = self.transition.step()

			#print "transitioning"
			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])

			if isDone:
				self.isTransitioning = False

			return False

		if not self.isTransitioning and self.hasTransitioned:

			if self.direction:
				pose = self.probe.getActualJointPose(-1)
			else:
				pose = self.probe.getActualJointPose(39)
				
			self.stabilityX.addData(pose[0])
			self.stabilityY.addData(pose[1])
	
			varX = self.stabilityX.getVar()
			varY = self.stabilityY.getVar()
			#print varX, varY
		
			if self.stabilityX.isStable() and self.stabilityY.isStable():
				
				diff = [pose[0]-self.lastPose[0], pose[1]-self.lastPose[1]]
				
				if self.val == 30 or self.val == -30:
					if not self.returnedSuccess:
						
						self.slipError = sqrt(diff[0]**2 + diff[1]**2)
						#print "pose error =", self.slipError
						
						
						self.returnedSuccess = True
						return True
					else:
						self.returnedSuccess = False
						
				self.stabilityX.reset()
				self.stabilityY.reset()
				self.lastPose = pose
				self.resetJoints()
				
				if self.val == 0:
					self.stabilityX.setThresh(1e-7)
					self.stabilityY.setThresh(1e-7)
				else:
					self.stabilityX.setThresh(1e-5)
					self.stabilityY.setThresh(1e-5)
				
			else:
				return False
		
			# 1. increment only after the end point becomes stable
			# 2. measure distance between new end point and old end point, and this determines our contact
			
			seq = [-30,-25,0,25,30]
			
			if self.up:
				if self.poseIndex == 0:
					self.poseIndex += 2
				else:
					self.poseIndex += 1
					
				self.val = seq[self.poseIndex]
				#self.val += 25
			else:
				if self.poseIndex == 4:
					self.poseIndex -= 2
				else:
					self.poseIndex -= 1
				self.val = seq[self.poseIndex]
				#self.val += -25
				
			if self.poseIndex >= 4:
				self.up = False
				
			if self.poseIndex <= 0:
				self.up = True
			
			#if self.val >= self.rail:
			#	self.val = self.rail
			#	self.up = False
			
			#if self.val <= -self.rail:
			#	self.val = -self.rail
			#	self.up = True
	
			if self.direction:
				for i in range(0, self.topJoint):
					# joint value
					self.setJoint(i, self.val)
			else:
				#print "setting", range(self.topJoint, self.probe.numSegs-1)
				for i in range(self.topJoint, self.probe.numSegs-1):
					# joint value
					self.setJoint(i, self.val)

			self.hasTransitioned = False

		if not self.hasTransitioned:
	
			cmdJoints = self.getJoints()

			initState = []
			for i in range(self.probe.numSegs-1):
				if cmdJoints[i] != None:
					initState.append(180.0*self.probe.getServo(i)/pi)
				else:
					initState.append(None)

			#print self.hasTransitioned, self.isTransitioning, "joints:", joints2
			targetState = self.getJoints()
			#print "init", initState
			#print "target", targetState

			#for i in range(len(targetState)):
			#	if targetState[i] == None:
			#		targetState[i] = 180.0*self.probe.getServo(i)/pi

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])
			transTime = int(4*errSum)
			
			print "poke transition time =", transTime

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			#self.transition.resetTime(200)		
			self.transition.resetTime(transTime)		

			# first steps
			self.transition.step()

			resJoints = self.transition.getJoints()
			self.resetJoints()
			self.mergeJoints([resJoints])

			self.hasTransitioned = True
			self.isTransitioning = True				

		#if self.val == 0 and self.up:
		#if self.val == self.rail or self.val == -self.rail or self.val == 0:
		#if self.val == self.rail or self.val == -self.rail:
		#if self.poseIndex == 0 or self.poseIndex == 4:
		#	return True
		#else:
		return False
	
	
	
	