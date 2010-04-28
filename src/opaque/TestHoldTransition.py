from SnakeControl import SnakeControl
from copy import *
from math import *

from behaviors.HoldTransition import HoldTransition

class TestHoldTransition(SnakeControl):

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
		
		self.nextState = False
		
		self.stateA = 0		
		self.globalTimer = 1

	def frameStarted(self):

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		" get the probe state "
		probeState = self.probe.getProbeState()

		" initialize probe "
		if self.stateA == 0:	

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

		elif self.stateA == 1:
			
			newJoints = [0.0 for i in range(39)]
			
			" create and initialize behavior "
			self.holdT = HoldTransition(self.robotParam)
			self.holdT.reset(probeState, newJoints)
								
			self.stateA = 2
			
			print "step 1 complete"
			
		# Anchor to the walls
		elif self.stateA == 2:
			
			isDone = self.holdT.step(probeState)
			joints1 = self.holdT.getJoints()
			
			self.mergeJoints([joints1])
			
			if isDone:
				
				currJoints = probeState['joints']

				for i in range(self.robotParam['numJoints']):					
					if currJoints[i] != 0.0:
						print "ERROR: joint", i, "is", currJoints[i], "not 0.0"

				print "step 2 complete"
				raise
			
								

				
