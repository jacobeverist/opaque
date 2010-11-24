
from Behavior import *
from copy import copy

class Transition(Behavior):
	
	def __init__(self, robotParam):
		Behavior.__init__(self, robotParam)
	
		self.resetTime()
		self.initJoints = [None for i in range(self.numJoints)]
		self.targetJoints = [None for i in range(self.numJoints)]
	
	def setInit(self, joints):
		self.initJoints = copy(joints)
	
	def setTarget(self, joints):
		self.targetJoints = copy(joints)
	
	def resetTime(self, newTime=500):
		self.transitionTime = 0
		
		#if newTime == 0:
		if newTime <= 1:
			#self.timeLength = 1
			self.timeLength = 2
		else:
			self.timeLength = newTime
	
	def step(self, probeState):
		Behavior.step(self, probeState)
			
		if not self.isStep():
			return False
		
		self.resetJoints()

		self.transitionTime += 1
		
		if self.transitionTime > self.timeLength:
			self.transitionTime = self.timeLength
		
		for i in range(self.numJoints):
			
			if self.initJoints[i] == None:
				self.setJoint(i, None)
				
			else:
				diff = self.targetJoints[i] - self.initJoints[i]
				
				val = 1.0 - float(self.timeLength - self.transitionTime) / self.timeLength
				
				newAngle = self.initJoints[i] + diff * val
				
				self.setJoint(i, newAngle)


		if self.transitionTime >= self.timeLength:
			return True
		else:
			return False
