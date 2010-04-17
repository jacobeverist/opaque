import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from HoldSlideTransition import HoldSlideTransition
from math import exp

class FrontExtend(Behavior):
	
	def __init__(self, probe, direction = True):
		Behavior.__init__(self, probe)

		self.transition = HoldSlideTransition(self.probe)

		self.direction = direction

		self.positions = []
		for i in range(self.probe.numSegs-1):
			self.positions.append(180.0 / pi * self.probe.getServoCmd(i))

		resultJoints = [None for i in range(40)]

		if self.direction:
			self.spliceJoint = 7
			for i in range(0,25):
				resultJoints[i] = 0.0
		else:
			self.spliceJoint = 31
			for i in range(15,self.probe.numSegs-1):
				resultJoints[i] = 0.0

		self.holdSlideT.reset(resultJoints, self.direction)
		
	def setDirection(self, isForward):

		self.direction = isForward

		if self.direction:
			self.spliceJoint = 7
		else:
			self.spliceJoint = 31

		self.holdSlideT.reset(resultJoints, self.direction)
		
	def setSpliceJoint(self, joint):
		self.spliceJoint = joint
	
	def getSpliceJoint(self, joint):
		return self.spliceJoint

	def reset(self, joints = [], direction = True):
		self.targetState = []
		self.positions = []

		self.direction = direction

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
			self.setJoint(i, self.positions[i])


		if self.isDone:
			self.isDone = False
			self.isTransitioning = False

			return True

		return False