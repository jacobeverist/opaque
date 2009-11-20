import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *

class Anchor(Behavior):
	
	def __init__(self, probe):
		Behavior.__init__(self, probe)

		# angle increment per segment
		self.per_seg = 2*pi / (NUM_SEGS / 5.0)
		
	def step(self):
		Behavior.step(self)
			
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.probe.numSegs-1):

			# joint value
			val = 70.0*cos(i*self.per_seg)
			self.setJoint(i, val)
			#self.probe.setServo(i, val)


		# never finished
		return False