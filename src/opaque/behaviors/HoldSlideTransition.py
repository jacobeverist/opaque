import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Transition import Transition
from math import exp

class HoldSlideTransition(Behavior):
	
	def __init__(self, probe, direction = True):
		Behavior.__init__(self, probe)

		self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.isDone = False
		#self.count = 0
		self.count = 4

		self.direction = direction

		self.targetState = []
		for i in range(self.probe.numSegs-1):
			self.targetState.append(None)

		self.positions = []
		for i in range(self.probe.numSegs-1):
			self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
		#print "reseting hold transition:", self.positions

		self.spliceJoint = 19

	def setSpliceJoint(self, joint):
		self.spliceJoint = joint
	

	def reset(self, joints = [], direction = True):
		self.targetState = []
		self.positions = []

		self.direction = direction
		
		#print "resetting HoldSlideTransition to:", joints
		
		if len(joints) > 0:
			for i in range(len(joints)):				
				self.targetState.append(joints[i])
				if joints[i] != None:
					self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
				else:
					self.positions.append(None)
					
		else:
			for i in range(self.probe.numSegs-1):
				self.targetState.append(180.0 / pi * self.probe.getServoCmd(i))
				self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
	
				
		#print "resetting:"
		#print self.positions
		#print joints
			
		self.isTransitioning = False
		self.hasTransitioned = False
		self.isDone = False
				
		#print "reseting hold slide transition:", self.positions

	def step(self):
		Behavior.step(self)

		#print self.isTransitioning, self.hasTransitioned
		#print self.isTransitioning, self.hasTransitioned, self.isDone
		
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step()

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

		for i in range(self.probe.numSegs-1):
			self.setJoint(i, self.positions[i])

		if not self.hasTransitioned:

			#print "setting next transition"
			#print "self.positions = ", self.positions
			initState = []
			for i in range(self.probe.numSegs-1):
				#print self.positions[i]
				if self.positions[i] != None:
					initState.append(180.0*self.probe.getServo(i)/pi)
				else:
					initState.append(None)

			#print "initState =", initState
			#print "self.targetState =", self.targetState
			#print initState

			targetState = copy(self.targetState)

			for i in range(len(targetState)):
				if initState[i] != None:
					if targetState[i] != None:
						diff = targetState[i] - initState[i]
						
						if self.direction:
							diff *= 1 / (1 + exp(i-2*self.count))
						else:
							diff *= 1 / (1 + exp((self.probe.numSegs-2-i)-2*self.count))
							
						targetState[i] = initState[i] + diff
					else:
						targetState[i] = initState[i]

			self.positions = copy(targetState)
			#print "self.positions = ", self.positions

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])

			#print "errSum = ", errSum
			#transTime = int(errSum/2.0)
			transTime = int(2*errSum)

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			self.transition.resetTime(transTime)		
			#self.transition.resetTime(200)		
			#self.transition.resetTime(1000)		
			
			# first steps
			self.transition.step()
			
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

		#print self.count

		" delay period to allow body to stabilize "
		if self.count > 10:
			self.count = 4
			self.isDone = False
			self.hasTransitioned = False
			self.isTransitioning = False
			
			#print "caseC"

			return True
		
		#print "transitions =", self.hasTransitioned, self.isTransitioning
		#print "caseD"
		return False