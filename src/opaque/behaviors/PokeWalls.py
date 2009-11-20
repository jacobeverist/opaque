import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Curl import *
from Anchor import *
from Transition import Transition

class PokeWalls(Behavior):

	def __init__(self, probe, direction = False, obstContact = 0):
		Behavior.__init__(self, probe)

		self.dir = direction
		self.obstContact = obstContact

		#self.MIN = 6
		self.MAX = 10
		self.MIN = 6
		#self.MAX = 16

		self.topJoint = self.MIN

		self.anchor = Anchor(self.probe)

		self.curl = Curl(self.probe)
		#self.curl.setTimerAliasing(10)
		self.curl.setTopJoint(self.topJoint)

		self.doneCount = 0

		self.direction = True

		self.reverse = False

		self.state = 0

		#self.transition = Transition(self.probe)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.f = open("pokeJoints.txt",'w')
		
		self.isInit = False

	def hasInitialized(self):
		return self.isInit
	
	def setDirection(self, val):
		#print "setting direction", val
		self.direction = val

		if self.direction:
			self.curl.setTopJoint(self.topJoint)
		else:
			self.curl.setTopJoint(-self.topJoint + 1 + self.probe.numSegs-2)

		self.curl.setDirection(self.direction)

	def step(self):
		Behavior.step(self)
		
		#print "Start:", self.isTransitioning, self.hasTransitioned

		"""
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step()

			#print "transitioning"
			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])

			if isDone:
				self.isTransitioning = False
				
				self.hasTransitioned = False

			return False
		"""
		if not self.isStep():
			return False

		self.resetJoints()

		finalDone = False

		self.anchor.step()
		
		isDone = False

		if self.state == 0:
	
			isDone = self.curl.step()
	
			if isDone:
				self.doneCount += 1
				self.state = 1
				self.isInit = True

	
		elif self.state == 1:
	
			if self.doneCount >= 2:
				if not self.reverse == True:
					self.topJoint += 2
					if self.direction:
						self.curl.setTopJoint(self.topJoint)
						self.curl.step()
					else:
						#print "setting topJoint", -self.topJoint + 1 + self.probe.numSegs-2
						self.curl.setTopJoint(-self.topJoint + 1 + self.probe.numSegs-2)
						self.curl.step()
					self.doneCount = 0
				else:
					self.topJoint -= 2
					if self.direction:
						self.curl.setTopJoint(self.topJoint)
						self.curl.step()
					else:
						#print "setting topJoint", -self.topJoint + 1 + self.probe.numSegs-2
						self.curl.setTopJoint(-self.topJoint + 1 + self.probe.numSegs-2)
						self.curl.step()
					self.doneCount = 0
	
			self.state = 0
	
			if self.curl.getSlipError() < 0.07:
				#print "capturedObstacle: ", self.curl.getSlipError()
				" temporarily disabling the obstacle detection callback for speed "
				self.obstContact(self.direction)
			else:
				pass
				#print "discardedObstacle: ", self.curl.getSlipError()
	
			if self.topJoint > self.MAX:
				#finalDone = True
				self.reverse = True
	
			if self.topJoint < self.MIN and self.reverse:
				finalDone = True
				self.reverse = False
				self.topJoint = self.MIN
	
			#print self.state, self.doneCount, self.waitTimer

		#joints1 = self.anchor.getJoints()
		joints2 = self.curl.getJoints()

		# joint value
		#self.mergeJoints([joints2, joints1])
		self.mergeJoints([joints2])

		#print self.state, self.doneCount, isDone, joints2
		#self.f.write(repr(self.state) + " " + repr(self.doneCount) + " " + repr(isDone) + " "+ repr(joints2))
		#self.f.write("\n")
		
		"""
		if not self.hasTransitioned:

			initState = []
			for i in range(self.probe.numSegs-1):
				if joints2[i] != None:
					initState.append(180.0*self.probe.getServo(i)/pi)
				else:
					initState.append(None)

			#print self.hasTransitioned, self.isTransitioning, "joints:", joints2
			targetState = self.getJoints()
			#print "init", initState
			#print "target", targetState

			#for i in range(len(targetState)):
			#	if targetState[i] == None:
			#		targetState[i] = 180.0*self.probe.getServo(i)/pi

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])
			transTime = int(errSum)
			
			print "poke transition time =", transTime

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
		"""
		# this would be just after we exited state 0
		#elif self.state == 1:
		#	self.hasTransitioned = False

		if finalDone:
			self.hasTransitioned = False
			self.curl.reset()
			self.isInit = False
			#print "done:", self.hasTransitioned, self.isTransitioning, self.topJoint, self.direction

		return finalDone

