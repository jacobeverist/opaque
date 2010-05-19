from SnakeControl import SnakeControl
from copy import copy
from math import fabs

from behaviors.AdaptiveStep import AdaptiveStep
from behaviors.PokeWalls import PokeWalls
from behaviors.AnchorTransition import AnchorTransition
from pose.AverageContacts import AverageContacts
from maps.MapGraph import MapGraph

class TestAdaptiveStep(SnakeControl):

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

		self.contacts = AverageContacts(self.probe)
		self.contacts.setMask( [1.0 for i in range(39)] )	
		self.contacts.step()
					
		self.mapGraph = MapGraph(self.probe, self.contacts)
		self.mapGraph.loadFile("testData/behaviorTest", 1)
		
		self.mapGraph.correctPoses2()
		self.mapGraph.synch()
		self.mapGraph.saveMap()
		
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
				self.probe.setServo(i, 0.0)
					
			self.contacts = AverageContacts(self.probe)
			self.contacts.setMask( [1.0 for i in range(39)] )	
			self.contacts.step()
				
			self.stateA = 1			
			print "step 0 complete"

		elif self.stateA == 1:

			self.contacts.setMask( [1.0 for i in range(39)] )	
			self.contacts.step()

			if self.contacts.isStable():
				self.stateA = 2
				print "isStable()"

		elif self.stateA == 2:
			
			" create and initialize behavior "
			self.anchorT = AnchorTransition(self.robotParam)
			
			self.stateA = 3
			
			print "step 1 complete"
			
		# Anchor to the walls
		elif self.stateA == 3:
			
			isDone = self.anchorT.step(probeState)
			joints1 = self.anchorT.getJoints()
			
			self.mergeJoints([joints1])
			
			if self.globalTimer % 10 == 0:
				rootPose = self.contacts.getAveragePose(19)	
				poses = []
				for i in range(self.robotParam['numJoints']):
					poses.append(self.probe.getJointWRTJointPose(rootPose, 19, i))					
				self.drawThings.plotRobotConfiguration(poses)
	
			if isDone:
				self.stateA = 4
				print "step 3 complete"
	
			
		elif self.stateA == 4:

			" create and initialize behavior "
			self.adaptiveStep = AdaptiveStep(self.robotParam, probeState, self.contacts, self.mapGraph, direction = True)
			
			self.stateA = 5
			
			print "step 4 complete"
			
		# Anchor to the walls
		elif self.stateA == 5:
			
			isDone = self.adaptiveStep.step(probeState)
			joints1 = self.adaptiveStep.getJoints()

			self.mergeJoints([joints1])

			if self.globalTimer % 200 == 0:
				rootPose = self.contacts.getClosestPose(19)	
				poses = []
				for i in range(self.robotParam['numJoints']):
					poses.append(self.probe.getJointWRTJointPose(rootPose, 19, i))					
				self.drawThings.plotRobotConfiguration(poses)

					
			if isDone:
												
				currJoints = probeState['joints']

				for i in range(0,17):
					
					if fabs(currJoints[i]) >= 0.0001:
						print "ERROR: joint", i, "is", currJoints[i], "not 0.0"

				print "step 5 complete"
				raise
			
								

				
