from Behavior import *
from Transition import Transition

from math import pi, fabs

class HoldTransition(Behavior):
	
	def __init__(self, robotParam):
		Behavior.__init__(self, robotParam)

		self.transition = Transition(robotParam)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.isDone = False
		self.count = 0

		self.positions = [0.0 for i in range(self.numJoints)]

	def reset(self, probeState, joints = []):
		
		stateJoints = probeState['joints']
		
		self.positions = []
		
		if len(joints) > 0:
			for i in range(len(joints)):			
				self.positions.append(joints[i])			
		else:
			for i in range(self.numJoints):
				self.positions.append(180.0 / pi * stateJoints[i])		

		self.count = 0
		self.isTransitioning = False
		self.hasTransitioned = False
		self.isDone = False

		print "resetting Hold Transition:", self.positions

	def step(self, probeState):
		Behavior.step(self, probeState)
		
		#print self.isTransitioning, self.hasTransitioned, self.isDone
		
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step(probeState)

			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])
		
			if isDone:
				self.isDone = True
				self.isTransitioning = False

			#print "HoldTransition transitions =", self.hasTransitioned, self.isTransitioning
			
			return False
		
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.numJoints):
			self.setJoint(i, self.positions[i])

		if not self.hasTransitioned:

			joints = probeState['joints']

			initState = []
			for i in range(self.numJoints):
				if self.positions[i] != None:
					initState.append(180.0*joints[i]/pi)
				else:
					initState.append(None)

			#print "new transition "

			#print initState

			targetState = self.getJoints()

			#print targetState

			#for i in range(len(targetState)):
			#	if targetState[i] == None:
			#		targetState[i] = 180.0*stateJoints[i]/pi

			#print targetState

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])

			#transTime = int(2*errSum)
			transTime = int(errSum)
			
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

		if self.count > 300:
			self.count = 0
			self.isDone = False
			self.hasTransitioned = False
			self.isTransitioning = False
			return True

		return False