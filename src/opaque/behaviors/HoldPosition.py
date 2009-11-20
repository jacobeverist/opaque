import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *

class HoldPosition(Behavior):
	
	def __init__(self, probe):
		Behavior.__init__(self, probe)

		self.positions = []
		for i in range(self.probe.numSegs-1):
			self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
		#print "reseting hold position:", self.positions
	
	def reset(self):
		self.positions = []
		#for i in range(self.probe.numSegs-1):
		#	self.positions.append(180.0 / pi * self.probe.getServoCmd(i))		
		for i in range(self.probe.numSegs-1):
			self.positions.append(180.0 / pi * self.probe.getServo(i))		
		#print "reseting hold position:", self.positions
		
	def step(self):
		Behavior.step(self)
			
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.probe.numSegs-1):
			self.setJoint(i, self.positions[i])

		# never finished
		return False