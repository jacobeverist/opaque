from Behavior import *
from Transition import Transition
from math import exp, fabs, pi
from copy import copy

class HoldSlideTransition(Behavior):
	
	def __init__(self, robotParam, direction = True):
		Behavior.__init__(self, robotParam)

		self.transition = Transition(robotParam)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.isDone = False
		self.count = 4

		self.direction = direction

		self.targetState = []
		for i in range(self.numJoints):
			self.targetState.append(None)

		self.positions = []
		for i in range(self.numJoints):
			self.positions.append(0.0)

		self.spliceJoint = 19

		self.mask = [0.0 for i in range(self.numJoints)]

	def getMask(self):
		return self.mask
	
	def setSpliceJoint(self, joint):
		self.spliceJoint = joint
	

	def reset(self, probeState, joints = [], direction = True):
		self.targetState = []
		self.positions = []

		self.direction = direction

		cmdJoints = probeState['cmdJoints']

		self.positions = []

		if len(joints) > 0:
			for i in range(len(joints)):				
				self.targetState.append(joints[i])
				if joints[i] != None:
					self.positions.append(180.0 / pi * cmdJoints[i])
				else:
					self.positions.append(None)
							
		else:
			for i in range(self.numJoints):
				self.targetState.append(180.0 / pi * cmdJoints[i])
				self.positions.append(180.0 / pi * cmdJoints[i])
	

		self.mask = [0.0 for i in range(self.numJoints)]
		
		self.isTransitioning = False
		self.hasTransitioned = False
		self.isDone = False
				
	def step(self, probeState):
		Behavior.step(self, probeState)

		#print self.isTransitioning, self.hasTransitioned
		#print self.isTransitioning, self.hasTransitioned, self.isDone
		
		stateJoints = probeState['joints']
		
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step(probeState)

			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])
		
			if isDone:
				self.isDone = True
				self.isTransitioning = False
			
			#print "transitions =", self.hasTransitioned, self.isTransitioning
			#print "caseA"
			return False
		
		if not self.isStep():
			#print "caseB"
			return False
		
		self.resetJoints()

		for i in range(self.numJoints):
			self.setJoint(i, self.positions[i])

		if not self.hasTransitioned:

			initState = []
			for i in range(self.numJoints):
				if self.positions[i] != None:
					initState.append(180.0*stateJoints[i]/pi)
				else:
					initState.append(None)

			targetState = copy(self.targetState)

			for i in range(len(targetState)):
				if initState[i] != None:
					if targetState[i] != None:
						diff = targetState[i] - initState[i]
						
						if self.direction:
							diff *= 1 / (1 + exp(i-2*self.count))
						else:
							diff *= 1 / (1 + exp((self.numJoints-1-i)-2*self.count))
							
						targetState[i] = initState[i] + diff
					else:
						targetState[i] = initState[i]

			self.mask = [0.0 for i in range(self.numJoints)]
			if self.direction:
				for i in range(self.numJoints):
					if 1. - 1. / (1. + exp(i-2*self.count)) >= 0.99:
						self.mask[i] = 1.0
			else:
				for i in range(self.numJoints):
					if 1. - 1. / (1. + exp((self.numJoints-1-i)-2*self.count)) >= 0.99:
						self.mask[i] = 1.0



			self.positions = copy(targetState)

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])

			transTime = int(2*errSum)

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			self.transition.resetTime(transTime)		
			
			# first steps
			self.transition.step(probeState)
			
			resJoints = self.transition.getJoints()
			self.resetJoints()
			self.mergeJoints([resJoints])

			self.hasTransitioned = True
			self.isTransitioning = True			

		if self.isDone:
			self.count += 1
			self.isDone = False
			self.hasTransitioned = False
			self.isTransitioning = False

		" delay period to allow body to stabilize "
		if self.count > 10:
			self.count = 4
			self.isDone = False
			self.hasTransitioned = False
			self.isTransitioning = False

			return True
		
		return False
	
	