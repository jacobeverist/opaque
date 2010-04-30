from SnakeControl import SnakeControl
from copy import *
from math import *

from behaviors.FrontExtend import FrontExtend
from pose.AverageContacts import AverageContacts

class TestFrontExtend(SnakeControl):

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
			
			self.contacts = AverageContacts(self.probe)
			self.contacts.setMask( [1.0 for i in range(39)] )	
			self.contacts.step()
	
			if self.globalTimer > 1000:
				self.stateA = 1			
				print "step 0 complete"


		elif self.stateA == 1:
		
			self.contacts.setMask( [1.0 for i in range(39)] )	
			self.contacts.step()

			if self.contacts.isStable():
				self.stateA = 2
				print "isStable()"

		elif self.stateA == 2:
			
			newJoints = [0.0 for i in range(39)]

			" create and initialize behavior "
			self.frontExtend = FrontExtend(self.robotParam, self.contacts)
			self.frontExtend.reset(probeState, newJoints)
								
			self.stateA = 3
			
			print "step 1 complete"
			
		# Anchor to the walls
		elif self.stateA == 3:
			
			isDone = self.frontExtend.step(probeState)
			joints1 = self.frontExtend.getJoints()
			
			self.mergeJoints([joints1])

			val = self.frontExtend.getMask()			
			self.contacts.setMask(val)			
			self.contacts.step()
			
			if isDone:
				
				currJoints = probeState['joints']

				for i in range(0,17):
					
					if fabs(currJoints[i]) >= 0.0001:
						print "ERROR: joint", i, "is", currJoints[i], "not 0.0"

				#print currJoints
				#for i in range(self.robotParam['numJoints']):					
				#	if currJoints[i] != 0.0:
				#		print "ERROR: joint", i, "is", currJoints[i], "not 0.0"

				print "step 2 complete"
				raise
			
								

				
