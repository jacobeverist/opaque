from Behavior import *
from Transition import Transition

from math import pi, cos
from copy import copy

class AnchorTransition(Behavior):
	
	def __init__(self, robotParam):
		Behavior.__init__(self, robotParam)

		# angle increment per segment
		self.per_seg = 2*pi / (self.robotParam['numSegs'] / 5.0)
		self.transition = Transition(self.robotParam)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.isDone = False
		self.count = 0
		
	def step(self, probeState):
		Behavior.step(self, probeState)
		
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step(probeState)

			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])
		
			if isDone:
				self.isDone = True
				self.isTransitioning = False
			
			return False
					
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.robotParam['numJoints']):

			# joint value
			val = 70.0*cos(i*self.per_seg)
			self.setJoint(i, val)
			
		if not self.hasTransitioned:
			
			joints = probeState['joints']

			initState = copy(joints)
			for i in range(self.robotParam['numJoints']):
				initState[i] = initState[i] * 180.0 / pi

			targetState = self.getJoints()

			for i in range(len(targetState)):
				if targetState[i] == None:
					targetState[i] = 180.0*joints[i]/pi

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			#self.transition.resetTime(200)		
			self.transition.resetTime(10)		
			
			# first steps
			self.transition.step(probeState)
			
			resJoints = self.transition.getJoints()
			self.resetJoints()
			self.mergeJoints([resJoints])

			self.hasTransitioned = True
			self.isTransitioning = True			

		if self.isDone:
			self.count += 1

		if self.count > 300:
			self.count = 0
			self.isDone = False
			self.hasTransitioned = False
			self.isTransitioning = False
			return True

		return False