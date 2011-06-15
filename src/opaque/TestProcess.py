from SnakeControl import SnakeControl
from copy import *
from math import *
import time

from behaviors.AnchorTransition import AnchorTransition

from pose.AverageContacts import AverageContacts
from maps.MapGraph import MapGraph
from maps.VisualGraph import VisualGraph

import numpy
import sys

class TestProcess(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe, drawThings):
		SnakeControl.__init__(self)
	
		self.drawThings = drawThings
		self.probe = probe
		
		robotParam = {}
		robotParam['numJoints'] = self.probe.numSegs-1
		robotParam['numSegs'] = self.probe.numSegs
		robotParam['segLength'] = self.probe.segLength
		robotParam['segWidth'] = self.probe.segWidth
		robotParam['maxTorque'] = self.probe.maxTorque
		
		self.numJoints = robotParam['numJoints']
		
		self.robotParam = self.probe.robotParam
		
		self.setTimerAliasing(1)
		
		self.direction = True
		
		" get the probe state "
		probeState = self.probe.getProbeState()
		torques = probeState['torques']

		self.torques = [torques[i] for i in range(self.numJoints)]
			
		self.localState = 0
		self.globalState = 0
		self.prevTime = 0
		
		self.stepDist = 0.14
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
		self.targetReached = False
		
		self.isAnchored = False

		self.postureCounter = 0

		self.isCapture = False
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

		if self.globalState == 0 and self.globalTimer > 2:
			
			" anchor the robot initially "

			isDone = self.doInitAnchor()
			#isDone = self.doHoldPosition()
			#isDone = False
			if isDone:
				pass
				self.globalState = 1
				
				#self.probe.savePose()
				
		elif self.globalState == 1:
			
			" create the motion model object "
			self.contacts = AverageContacts(self.probe)
			
			self.globalState = 2
			
		elif self.globalState == 2:

			" get stable reference nodes on the snake body "
			self.contacts.setMask( [1.0 for i in range(self.numJoints)] )	
			self.contacts.step(probeState)
			
			if self.contacts.isStable():
				self.globalState = 3
				self.lastPose = self.contacts.getAveragePose(0)

		elif self.globalState == 3:
			
			" create the mapping object "
			#self.mapGraph = MapGraph(self.probe, self.contacts, isStable = True)
			self.mapGraph = VisualGraph(self.probe, self.contacts)
			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()

			#self.mapGraph.testHypotheses("testData/sensorLocalize2", 175, "cornerHypotheses_2011_06_05.txt")
			#self.mapGraph.loadSeries("testData/sensorLocalize2", 175)
			self.mapGraph.loadSeries("testData/sensorLocalize2", 161)
			
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 175)
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 50)
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 50)

			#self.mapGraph.sensorTest("testData/sensorLocalize2", 50)
			#self.mapGraph.sensorTest("testData/sensorLocalize2", 176)

			exit()
			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			self.mapGraph.synch()
			self.mapGraph.saveMap()
			self.mapGraph.poseGraph.drawConstraints(0)
			#self.mapGraph.saveLocalMap()
	
			#exit()
	
			#self.mapGraph.numNodes
	
			self.mapGraph.poseGraph.correctPoses3()
			self.mapGraph.synch()
			self.mapGraph.saveMap()
			self.mapGraph.poseGraph.drawConstraints(1)
			#self.mapGraph.renderConstraints()

			exit()
			
			
			self.restState = deepcopy(probeState)
			
			#self.globalState = 6
			self.globalState = 4
			#self.globalState = 9

			#self.restState = deepcopy(probeState)
			#self.mapGraph.newNode(self.stepDist, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.globalState = 6

			print "Start Time:", time.clock()

	def doInitAnchor(self):
		
		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			
			print "START:  doInitAnchor()"
			
			self.behavior = AnchorTransition(self.robotParam)
			
			self.localState = 1

		" immediately begin running behavior "
		if self.localState == 1:
			
			isDone = self.behavior.step(probeState)
			joints1 = self.behavior.getJoints()
			
			self.mergeJoints([joints1])
			
			if isDone:
				self.localState = 0
				self.behavior = 0
				
				self.isAnchored = True
				
				print "FINISH:  doInitAnchor()"
				return True
		
		return False
