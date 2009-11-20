import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Transition import Transition

class AnchorTransition(Behavior):
	
	def __init__(self, probe):
		Behavior.__init__(self, probe)

		# angle increment per segment
		self.per_seg = 2*pi / (NUM_SEGS / 5.0)
		self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.isDone = False
		self.count = 0
		
	def step(self):
		Behavior.step(self)
		
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step()

			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])
		
			if isDone:
				self.isDone = True
				self.isTransitioning = False
			
			return False
					
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.probe.numSegs-1):

			# joint value
			val = 70.0*cos(i*self.per_seg)
			self.setJoint(i, val)
			#self.probe.setServo(i, val)
			
		if not self.hasTransitioned:

			initState = []
			for i in range(self.probe.numSegs-1):
				initState.append(180.0*self.probe.getServo(i)/pi)

			targetState = self.getJoints()

			for i in range(len(targetState)):
				if targetState[i] == None:
					targetState[i] = 180.0*self.probe.getServo(i)/pi

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			self.transition.resetTime(200)		
			
			# first steps
			self.transition.step()
			
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