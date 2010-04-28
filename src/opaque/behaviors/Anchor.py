from Behavior import Behavior
from math import pi, cos

class Anchor(Behavior):
	
	def __init__(self, robotParam):
		Behavior.__init__(self, robotParam)

		# angle increment per segment
		self.per_seg = 2*pi / (self.robotParam['numSegs'] / 5.0)
		
	def step(self, probeState):
		Behavior.step(self, probeState)
			
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.numJoints):

			# joint value
			val = 70.0*cos(i*self.per_seg)
			self.setJoint(i, val)

		# never finished
		return False