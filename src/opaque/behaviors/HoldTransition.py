import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Transition import Transition

class HoldTransition(Behavior):
	
	def __init__(self, probe):
		Behavior.__init__(self, probe)

		self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.isDone = False
		self.count = 0

		self.positions = []
		for i in range(self.probe.numSegs-1):
			self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
		#print "reseting hold transition:", self.positions

	def reset(self, joints = []):
		self.positions = []
		
		#print "resetting HoldTransition to:", joints
		
		if len(joints) > 0:
				
			for i in range(len(joints)):
			
				self.positions.append(joints[i])
		else:
			for i in range(self.probe.numSegs-1):
			
				self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
	
				
		#print "resetting:"
		#print self.positions
		#print joints
			
		self.isTransitioning = False
		self.hasTransitioned = False
		self.isDone = False
				
		#print "reseting hold transition:", self.positions

	def step(self):
		Behavior.step(self)
		
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

			#print "HoldTransition transitions =", self.hasTransitioned, self.isTransitioning
			
			return False
		
		if not self.isStep():
			return False
		
		self.resetJoints()

		for i in range(self.probe.numSegs-1):
			self.setJoint(i, self.positions[i])

		if not self.hasTransitioned:

			initState = []
			for i in range(self.probe.numSegs-1):
				if self.positions[i] != None:
					initState.append(180.0*self.probe.getServo(i)/pi)
				else:
					initState.append(None)

			#print "new transition "

			#print initState

			targetState = self.getJoints()

			#print targetState

			#for i in range(len(targetState)):
			#	if targetState[i] == None:
			#		targetState[i] = 180.0*self.probe.getServo(i)/pi

			#print targetState

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

		if self.count > 300:
			self.count = 0
			self.isDone = False
			self.hasTransitioned = False
			self.isTransitioning = False
			return True

		return False