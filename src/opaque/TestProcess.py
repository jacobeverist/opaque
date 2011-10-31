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

			# cross - junctions, test 6
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
			#wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,3.6], [-3.6,3.6], [-3.6,0.2], [2.0,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]			

			# 45 degree junction, right, test 5
			#WLEN2 = 5.0
			#WLEN3 = 5.5
			#wall1 = [[-14.0, -0.2], [-3.0, -0.2]]
			#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			#w2 = wall2[0]			
			#wall3 = [[w2[0] + 0.4*cos(pi/6) - WLEN3*cos(pi/3), w2[1] - 0.4*sin(pi/6) - WLEN3*sin(pi/3)], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
			#walls = [wall1, wall2, wall3]

			# L - junction, right, test 4
			#wall1 = [[-14.0, -0.2], [-3.6,-0.2], [-3.6,0.2]]
			#wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,3.6], [-3.6,3.6], [-3.6,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]
			
			# L - junction, left, test 3
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [-3.6,0.2]]
			#wall2 = [[-14.0, 0.2], [-3.6,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]

			# T - junctions, test 2
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
			#wall2 = [[-14.0, 0.2], [-1.0, 0.2], [-1.0,3.6], [-0.6,3.6], [-0.6,0.2], [2.0,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]

			# Y - junction, test 1 			

			
			WLEN = 3.0
			WLEN2 = 5.0
			wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
			wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			w1 = wall1[2]
			w2 = wall2[0]
			
			wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
			lp = wall3[0]
			rp = wall3[2]
			
			wall6 = [lp, [lp[0] + WLEN*cos(pi/6), lp[1] + WLEN*sin(pi/6)]]
			wall6.append([wall6[1][0] + 0.4*cos(pi/3), wall6[1][1] - 0.4*sin(pi/3)])
			wall6.append([wall6[2][0] - WLEN*cos(pi/6), wall6[2][1] - WLEN*sin(pi/6)])
			wall6.append([wall6[3][0] + WLEN*cos(pi/3), wall6[3][1] - WLEN*sin(pi/3)])
			wall6.append([wall6[4][0] - 0.4*cos(pi/6), wall6[4][1] - 0.4*sin(pi/6)])
			wall6.append(w1)
			wall6.reverse()
			walls = [wall1, wall2, wall3, wall6]
			


			
			for wall in walls:
				for i in range(len(wall)):
					p = copy(wall[i])
					p[0] += 6.0
					wall[i] = p
			self.mapGraph.loadWalls(walls)

			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()

			#self.mapGraph.testHypotheses("testData/sensorLocalize2", 175, "cornerHypotheses_2011_06_05.txt")
			#self.mapGraph.loadSeries("testData/sensorLocalize2", 175)
			#self.mapGraph.loadSeries("testData/sensorLocalize2", 161)
			#self.mapGraph.loadSeries("testData/sensorLocalize2", 28)
			#self.mapGraph.loadSeries("testData/sensorLocalize2", 10)
			#self.mapGraph.loadSeries("testData/sensorLocalize2", 53)
			#self.mapGraph.loadSeries("testData/junctionTest6", 16)
			#self.mapGraph.loadSeries("testData/junctionTest5", 94)
			#self.mapGraph.loadSeries("testData/junctionTest4", 22)
			#self.mapGraph.loadSeries("testData/junctionTest3", 19)
			#self.mapGraph.loadSeries("testData/junctionTest2", 192)
			#self.mapGraph.loadSeries("testData/junctionTest1", 117)
			#exit()
			
			#self.mapGraph.loadSeries("testData/junctionTest1", 117)

			#self.mapGraph.loadSeries("testData/junctionTest7", 40)
			#self.mapGraph.loadSeries("testData/junctionTest7", 28)
			#self.mapGraph.loadSeries("testData/junctionTest8", 40)
			#self.mapGraph.loadSeries("testData/junctionTest9", 40)
			#self.mapGraph.loadSeries("testData/junctionTest11", 40)
			#self.mapGraph.loadSeries("testData/junctionTest11", 30)
			#self.mapGraph.loadSeries(".", 40)
			#self.mapGraph.loadSeries(".", 3)
			self.mapGraph.loadSeries("testData/junctionTest13", 40)

			#self.mapGraph.loadSeries("testData/junctionTest8", 40)
			#self.mapGraph.loadSeries("testData/junctionTest8", 20)


			#self.mapGraph.loadSeries("testData/junctionTest1", 13)
			#self.mapGraph.loadSeries("testData/junctionTest1", 66)
			#self.mapGraph.loadSeries("tes tData/junctionTest1", 67)

			#self.mapGraph.loadSeries("testData/junctionTest1", 25)
			#self.mapGraph.loadSeries("testData/junctionTest2", 26)
			#self.mapGraph.loadSeries("testData/junctionTest3", 19)
			#self.mapGraph.loadSeries("testData/junctionTest4", 22)
			#self.mapGraph.loadSeries("testData/junctionTest5", 28)
			#self.mapGraph.loadSeries("testData/junctionTest6", 16)
			
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 175)
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 50)
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 50)

			#self.mapGraph.sensorTest("testData/sensorLocalize2", 50)
			#self.mapGraph.sensorTest("testData/sensorLocalize2", 176)

			#self.mapGraph.drawMap()
			exit()
			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			#self.mapGraph.poseGraph.drawConstraints(0)
			#self.mapGraph.saveLocalMap()
	
			exit()
	
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
