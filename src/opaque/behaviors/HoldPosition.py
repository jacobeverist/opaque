from Behavior import *
from Transition import Transition

from math import pi

#from common import *

class HoldPosition(Behavior):
	
	def __init__(self, robotParam):
		Behavior.__init__(self, robotParam)

		self.positions = [0.0 for i in range(self.numJoints)]
	
	def reset(self, probeState):
		
		joints = probeState['joints']
		
		self.positions = []
		for i in range(self.numJoints):
			self.positions.append(180.0 / pi * joints[i])		
		
	def step(self, probeState):
		Behavior.step(self, probeState)
			
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.numJoints):
			self.setJoint(i, self.positions[i])

		# never finished
		return False