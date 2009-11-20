import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Transition import Transition

class Spring(Behavior):
	
	def __init__(self, probe):
		Behavior.__init__(self, probe)

		# angle value for curling joints
		self.val = 0
		self.rail = 25.0
		self.up = True
		self.topJoint = 19

		self.poseWindow = []
		
		self.lastPose = [0.0,0.0,0.0]
		
		self.poseIndex = 2
		
		self.returnedSuccess = False
		
		self.direction = True

		self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = True
	
		self.flip = True
	
	def setDirection(self, val):

		self.direction = val

		print "spring direction =", self.direction 
		self.resetJoints()
		
		#if self.direction:
			#for i in range(0, self.topJoint):
				# joint value
				#self.setJoint(i, self.val)
		#else:
			#print "setting", range(self.topJoint, self.probe.numSegs-1)
			#for i in range(self.topJoint, self.probe.numSegs-1):
				# joint value
				#self.setJoint(i, self.val)		
				
	def setTopJoint(self, node):
		self.topJoint = node

	def step(self):
		Behavior.step(self)

		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step()

			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])

			if isDone:
				self.isTransitioning = False
				return True

			return False

		if not self.isStep():
			return False

		if not self.isTransitioning and self.hasTransitioned:

			#print "setting up transition"
			" 1. setup target configuration and begin the transition "

			if self.direction:
				joints = range(0, self.topJoint)
			else:
				joints = range(self.topJoint, self.probe.numSegs-1)
				
			#print joints
			
			for j in joints:
				if self.flip:
					if j % 2 == 0:
						self.setJoint(j, 145.0)
					else:
						self.setJoint(j, -145.0)
				else:
					if j % 2 == 0:
						self.setJoint(j, 0.0)
					else:
						self.setJoint(j, 0.0)
					
			self.hasTransitioned = False

		if not self.hasTransitioned:
	
			cmdJoints = self.getJoints()

			initState = []
			for i in range(self.probe.numSegs-1):
				if cmdJoints[i] != None:
					initState.append(180.0*self.probe.getServo(i)/pi)
				else:
					initState.append(None)

			targetState = self.getJoints()

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])
			transTime = int(errSum)
			
			#print "poke transition time =", transTime
			#print initState
			#print targetState

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			#self.transition.resetTime(200)		
			self.transition.resetTime(transTime)		

			# first steps
			self.transition.step()

			resJoints = self.transition.getJoints()
			self.resetJoints()
			self.mergeJoints([resJoints])

			self.hasTransitioned = True
			self.isTransitioning = True				

		return False
	
	
	
	