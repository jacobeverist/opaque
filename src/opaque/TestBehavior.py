from SnakeControl import SnakeControl
from copy import *
from math import *

from behaviors.HoldPosition import HoldPosition

class TestBehavior(SnakeControl):

	"""
	Unit tests
	
	1. never returns
	2. reset position commands the probe to position
	3. None values in the joints do not affect the probe commands
	4. Changing the joint real values will be compensated by behavior
	
	Setup the probe
	Add the controller
	Drive the joints
	Step and get results
	
	"""

	def __init__(self, probe, drawThings):
		SnakeControl.__init__(self)
	
		self.drawThings = drawThings
		self.probe = probe
		
		self.robotParam = self.probe.robotParam
		
		self.setTimerAliasing(1)
		
		" behaviors "
		self.holdP = HoldPosition(self.robotParam)

		self.nextState = False
		

		self.stateA = 0
		self.prevTime = 0
		
		self.isAnchored = False
		self.renderTransition = False		
		self.isInitialized = False
		self.globalTimer = 1


	def frameStarted(self):
		
		self.probe.setAnchor(self.isAnchored)

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		" get the probe state "
		probeState = self.probe.getProbeState()

		" initialize "
		if self.stateA == 0:
		
			" create and initialize behavior "
			self.holdP = HoldPosition(self.robotParam)
			self.holdP.reset(probeState)
			
			for i in range(self.robotParam['numJoints']):
				modVal = i % 4
				
				if modVal == 0:
					self.probe.setServo(i, 0.0)
				elif modVal == 1:
					self.probe.setServo(i, 90.0)
				elif modVal == 2:
					self.probe.setServo(i, 0.0)
				elif modVal == 3:
					self.probe.setServo(i, -90.0)
					
			self.stateA = 1
			
			print "step 0 complete"
			
		# Anchor to the walls
		elif self.stateA == 1:
			
			isDone = self.holdP.step(probeState)
			joints1 = self.holdP.getJoints()
			
			self.mergeJoints([joints1])
			
			if self.nextState:

				if isDone:
					print "ERROR: isDone is true in step 2"
				
				currJoints = probeState['joints']

				for i in range(self.robotParam['numJoints']):
					modVal = i % 4
					
					if modVal == 0:
						if modVal == 0 and currJoints[i] != 0.0:
							print "ERROR: joint", i, "is", currJoints[i], "not 0.0"
						elif modVal == 1 and currJoints[i] != 0.0:
							print "ERROR: joint", i, "is", currJoints[i], "not 90.0"
						elif modVal == 2 and currJoints[i] != 0.0:
							print "ERROR: joint", i, "is", currJoints[i], "not 0.0"
						elif modVal == 3 and currJoints[i] != 0.0:
							print "ERROR: joint", i, "is", currJoints[i], "not -90.0"

				print "step 2 complete"
				self.testSuccess()
				
			else:
				self.nextState = True
			
				if isDone:
					print "ERROR: isDone is true in step 1"
					
				currJoints = probeState['joints']

				for i in range(self.robotParam['numJoints']):
					modVal = i % 4
					
					if modVal == 0:
						if modVal == 0 and currJoints[i] != 0.0:
							print "ERROR: joint", i, "is", currJoints[i], "not 0.0"
						elif modVal == 1 and currJoints[i] != 90.0:
							print "ERROR: joint", i, "is", currJoints[i], "not 90.0"
						elif modVal == 2 and currJoints[i] != 0.0:
							print "ERROR: joint", i, "is", currJoints[i], "not 0.0"
						elif modVal == 3 and currJoints[i] != -90.0:
							print "ERROR: joint", i, "is", currJoints[i], "not -90.0"

				print "step 1 complete"
				
			
			if isDone:
				self.stateA = 2
				print "going to state ", self.stateA

				self.holdP.reset(probeState)

				self.isAnchored = False
				
