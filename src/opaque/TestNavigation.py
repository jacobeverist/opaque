import random
random.seed(0)

from SnakeControl import SnakeControl
from copy import *
from math import *
import time

from behaviors.AnchorTransition import AnchorTransition
from behaviors.HoldTransition import HoldTransition
from behaviors.HoldPosition import HoldPosition
from behaviors.FrontExtend import FrontExtend
from behaviors.PokeWalls import PokeWalls
from behaviors.PathStep import PathStep
from behaviors.AdaptiveStep import AdaptiveStep

from pose.AverageContacts import AverageContacts
from maps.MapGraph import MapGraph

import numpy
import sys

class TestNavigation(SnakeControl):

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
		
		self.travelDir = True
		
		" pose estimation "
		self.contacts = AverageContacts(self.probe)
		self.contacts.setTimerAliasing(5)
		
		" maps "
		#self.mapGraph = MapGraph(self.probe, self.contacts)
		#self.mapGraph.loadFile("testData/correctionTest", 9)
		#self.mapGraph.loadFile("testData/correctionTest", 3)
		
		#self.mapGraph.correctPoses2()
		#self.mapGraph.synch()
		#self.mapGraph.saveMap()
		
				
		" get the probe state "
		probeState = self.probe.getProbeState()
		torques = probeState['torques']

		self.torques = [torques[i] for i in range(self.numJoints)]
		

		#self.targetPoint = [3.0, -1.0]
		self.targetPoint = []
		
		self.exploreRoot = [0.5,0.0]
		self.exploreRoot = [-3.5,0.0]
		self.pathStep = 0
		
		self.localPathDirection = True
		self.localPathState = 0
		self.localState = 0
		self.globalState = 0
		self.prevTime = 0
		
		self.stepDist = 0.14
		
		# error tracking for path-following
		#self.lastDist = 1e100
		self.lastDist1 = 1e100
		self.lastDist2 = 1e100
		self.distCount = 0
		self.targetReached = False
		self.isCollided = False
		
		
		self.isAnchored = False

		self.postureCounter = 0

		self.isCapture = False
		self.renderTransition = False
		
		self.isInitialized = False

		self.globalTimer = 1

		self.termPathID = -1
		
		self.inversion = 1
		
		
	def updateMaps(self):
		pass
	
	def saveMaps(self):
		self.mapGraph.saveMap()

	def grabImage(self):

		#inc = 50
		inc = 200
		#inc = 1000
		#inc = 500

		if self.globalTimer % inc == 0:
			#print "saving scene%06u.png" % (self.globalTimer/inc)
			self.drawThings.saveView("scene%06u.png" % (self.globalTimer/inc))
			
			#poses = self.probe.getPose()		
			#f = open("poses%06u.txt" % (self.globalTimer/inc), 'w')
			#f.write(repr(poses))
			#f.close()

	def grabPosture(self, isForce = False):
		
		if self.isCapture:
			inc = 50
			
			if self.globalTimer % inc == 0 or isForce:
				
				localPosture = []
				
				for j in range(self.probe.numSegs-1):
					localPosture.append(self.probe.getJointPose([0.0,0.0,0.0], self.mapGraph.currNode.rootNode, j))
				
				posture_file = open("posture_snap_%03u_%06u.txt" % (self.mapGraph.numNodes-1, self.postureCounter), 'w')		
				posture_file.write(repr(localPosture))
				posture_file.write('\n')

				f = open("estpose_%03u_%06u.txt" % (self.mapGraph.numNodes-1, self.postureCounter), 'w')
				f.write(repr(self.mapGraph.currNode.estPose))
				f.write("\n")
				f.close()
		
				f = open("gndpose_%03u_%06u.txt" % (self.mapGraph.numNodes-1, self.postureCounter), 'w')
				gndPose = self.probe.getActualJointPose(self.mapGraph.currNode.rootNode)
				f.write(repr(gndPose))
				#f.write(repr(self.mapGraph.currNode.gndPose))
				f.write("\n")
				f.close()

				self.postureCounter += 1

		else:
			self.postureCounter = 0
			
	def grabAngles(self):

		if self.isCapture:
			inc = 10
			
			if self.globalTimer % inc == 0:
				
				angles = []
				for i in range(39):
					angles.append(self.probe.getServo(i))
					
				self.angle_output.write(repr(angles))
				self.angle_output.write('\n')
				
				poses = []
				for i in range(40):
					poses.append(self.probe.getActualJointPose(i-1))

				self.global_output.write(repr(poses))
				self.global_output.write('\n')
		
	def frameStarted(self):
		
		
		#self.grabPosture()
		
		self.grabImage()
		#self.grabAngles()


		if self.globalState > 1 and self.globalTimer % 10 == 0:
			pnts = []
			for i in range(self.numJoints):
				#pnts.append(self.probe.getActualJointPose(i))
				#pnts.append(self.contacts.getAveragePose(i))
				pnts.append(self.contacts.getClosestPose(i))
			self.drawThings.renderPoints(pnts)

		#if (self.globalState == 6 or self.globalState == 8) and self.globalTimer % 50 == 0:
		#	probeState = self.probe.getProbeState()
		# 	print "pokeWalls torques:", probeState['torques']

		#if self.globalState == 5 and self.globalTimer % 100 == 0:
		#	probeState = self.probe.getProbeState()
		#	print "adaptiveStep torques:", probeState['torques']
		
		#if self.globalTimer % 4000 == 0:
		#	if self.stateA <= 0 or self.stateA == 5 or self.renderTransition:
		#		self.mapGraph.drawEstBoundary()
		#		self.renderTransition = False

		self.probe.setAnchor(self.isAnchored)

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		" get the probe state "
		probeState = self.probe.getProbeState()

		#print "state", self.globalState
		#print probeState['joints']
		#if self.globalTimer % 100 == 0:
		#	print "state:", self.globalState

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
				self.lastPose1 = self.contacts.getAveragePose(0)
				self.lastPose2 = self.contacts.getAveragePose(39)

		elif self.globalState == 3:
			
			" create the mapping object "
			self.mapGraph = MapGraph(self.probe, self.contacts, isStable = True)

			#self.mapGraph.loadFile("testData/sweep5", 8)
			#self.mapGraph.loadFile("testData/sweep5", 30)
			#self.mapGraph.loadFile("testData/sweep5", 10)
			#self.mapGraph.correctPoses3()
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			#exit()

			initPose = [-0.83215373754501343, 0.046221695840358734, -2.1473097801208496]
			#initpose = [-0.70651340484619141, 0.02392122894525528, -2.0914869430965584]

			#self.mapGraph.initNode(initPose, self.travelDir)
			#self.mapGraph.forceUpdate()
			#self.mapGraph.relaxCorrect()
			#self.mapGraph.synch()
			#self.mapGraph.saveLocalMap()
			#self.mapGraph.saveMap()

			#exit()
			
			#exit()
			#self.mapGraph.loadFile("testData/correctionTest", 43)
			#self.mapGraph.loadFile("uncorrectedNX", 6)
			#self.mapGraph.loadFile("testData/mapBuild_21June2010", 66)
			#self.mapGraph.loadFile("testData/poseTest", 3)
			#self.mapGraph.loadFile("testData/poseTest", 5)
			#self.mapGraph.loadFile("testData/poseTest", 8)
			#self.mapGraph.loadFile("testData/poseTest", 16)
			#self.mapGraph.loadFile("testData/fixedPoseTest5", 15)
			#self.mapGraph.loadFile("testData/bulletTest2", 17)
			#self.mapGraph.loadFile("testData/constraints", 17)
			#self.mapGraph.loadFile("testData/backtrack1", 7)

			#self.mapGraph.newNode(0.0, self.travelDir)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()

			#self.mapGraph.loadFile("testData/backtrack1", 7)
			#self.mapGraph.loadFile("testData/backtrack2", 47)
			#self.mapGraph.loadFile("testData/backtrack2", 47)
			#self.mapGraph.loadSeries("resultProcess_2012_04_13_3", 14)
			#self.mapGraph.restoreSeries("resultProcess_2012_05_02_2", 4)
			#self.mapGraph.restoreSeries("resultProcess_2012_05_02_2", 4)

			#self.mapGraph.restoreSeries("resultProcess_2012_05_03",10)
			#self.mapGraph.restoreSeries("resultProcess_2012_05_04", 30)
			#self.mapGraph.restoreSeries("resultProcess_2012_05_04", 32)

			
			#self.mapGraph.restoreSeries("resultProcess_2012_05_04", 10)
			#self.mapGraph.insertPose(self.probe.getActualJointPose(19), self.travelDir)

			#self.mapGraph.restoreSeries("result_2012_08_22", 50)
			#self.mapGraph.insertPose(self.probe.getActualJointPose(19), self.travelDir)
			
			
			#exit()
			#E = self.computeCovar()
			#print repr(E)


			
			#self.mapGraph.newNode(0.0, self.travelDir)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			#self.mapGraph.saveLocalMap()
	
			#exit()
	
			#self.mapGraph.correctPoses3()
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			#exit()
			
			
			self.restState = deepcopy(probeState)
			
			#self.globalState = 6
			self.globalState = 4
			#self.globalState = 10

			#self.restState = deepcopy(probeState)
			#self.mapGraph.newNode(self.stepDist, self.travelDir)
			#self.mapGraph.forceUpdate(False)
			#self.globalState = 6
			self.currPose1 = self.contacts.getAverageSegPose(0)

			print "Start Time:", time.clock()


			
		elif self.globalState == 4:
			
			" extend the front of the snake "
			isDone = self.doFrontExtend()
			
			if isDone:
				self.globalState = 5
				#self.globalState = 7

		elif self.globalState == 5:

			" do a concertina gait step "
			isDone = self.doAdaptiveStep()
			
			if isDone:
				self.restState = deepcopy(probeState)

				self.currPose = self.contacts.getAverageSegPose(0)
				deltaDist = sqrt((self.currPose[0]-self.lastPose1[0])**2 + (self.currPose[1]-self.lastPose1[1])**2)
				print "deltaDist =", deltaDist
	
				#print "Stop Time:", time.clock()

				#exit()

				faceDir = True
				self.mapGraph.newNode(faceDir, self.travelDir)
				self.mapGraph.forceUpdate(faceDir)
				
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				
				self.globalState = 6
				#self.globalState = 9
				#self.globalState = 7

				#self.probe.restorePose()
				self.isCapture = True
				#self.grabPosture(True)
				
		elif self.globalState == 6:

			" do a forward poking behavior to sense the environment"
			isDone = self.doPokeWalls(True)
			
			if isDone:
				self.globalState = 7
			
		elif self.globalState == 7:

			" return to the rest position "
			isDone = self.doReturnToRest()
			
			if isDone:
				#self.globalState = 9
				self.globalState = 8

				self.mapGraph.correctPosture()
				#self.mapGraph.localizeCurrentNode()
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)

				#self.mapGraph.synch()
				#self.mapGraph.saveMap()
				#self.mapGraph.saveLocalMap()
				#self.mapGraph.saveState()
				self.mapGraph.drawConstraints()
				
				
				faceDir = False		
				self.mapGraph.newNode(faceDir, self.travelDir)
				self.mapGraph.forceUpdate(faceDir)

		elif self.globalState == 8:
	
			" do a backward poking behavior to sense the environment"
			isDone = self.doPokeWalls(False)
			
			if isDone:
				self.globalState = 9
		
		elif self.globalState == 9:
			
			" return to the rest position "
			isDone = self.doReturnToRest()
			
			if isDone:
				print "Stop Time:", time.clock()

				self.isCapture = False

				#exit()
				
				#self.mapGraph.correctPoses2()
				self.mapGraph.correctPosture()
				#self.mapGraph.localizeCurrentNode()
				self.mapGraph.pairDone()
				self.mapGraph.localizePose()
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose = self.contacts.getAverageSegPose(0)

				#self.mapGraph.relaxCorrect()
				#self.mapGraph.synch()
				#self.mapGraph.saveMap()
				self.mapGraph.saveLocalMap()
				self.mapGraph.saveState()
				self.mapGraph.drawConstraints()				
				

				#if self.mapGraph.numNodes > 24:
				#	exit()
				
				#self.currPose = self.contacts.getAverageSegPose(0)
				#deltaDist = sqrt((self.currPose[0]-self.lastPose[0])**2 + (self.currPose[1]-self.lastPose[1])**2)

				" FORCE to the probe to continue walking forward "
				#deltaDist = 1.0
				#print "deltaDist =", deltaDist
	
				#self.lastPose = self.currPose
				
				
				frontSum = 0.0
				frontProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-2].frontProbeError
				for n in frontProbeError:
					frontSum += n
				foreAvg = frontSum / len(frontProbeError)
				

				backSum = 0.0
				backProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-1].backProbeError
				for n in backProbeError:
					backSum += n
				backAvg = backSum / len(backProbeError)

				print "frameStarted: foreAvg =", foreAvg
				print "frameStarted: backAvg =", backAvg
				print "travelDir =", self.travelDir
								
				#print "foreAvg =", foreAvg
				#if self.mapGraph.numNodes >=6:
				#	foreAvg = 2.0
				#foreAvg = 0.0
				
				if self.travelDir:
					
					#if foreAvg >= 1.4:
					if foreAvg >= 1.1:
						self.globalState = 10
						self.isCollided = True
					else:
						self.globalState = 4

				else:

					#if backAvg >= 1.4:
					if backAvg >= 1.1:
						self.globalState = 10
						self.isCollided = True
					else:
						self.globalState = 4
					

				#self.globalState = 10
					
				#if self.mapGraph.numNodes > 26:
				#	self.travelDir = False	
				#if self.mapGraph.numNodes > 40:
				#	exit()	
					
		elif self.globalState == 10:

			self.contacts.setCautious(True)			
			" select the destination and generate path "
			
			" select the point to go to "
			#frontierPoint = self.mapGraph.selectNextFrontier()
			#frontierPoint = self.mapGraph.selectNextFrontier()
			
			if self.targetPoint == []:
				frontierPoint, self.termPathID = self.mapGraph.selectNextDestination()
			
			else:
				frontierPoint = self.targetPoint
				self.termPathID = -1
				
							
			"""
			self.currPose1 = self.contacts.getAverageSegPose(0)
			self.currPose2 = self.contacts.getAverageSegPose(39)
			dest = self.localWayPoints[0]			
			self.lastDist1 = self.mapGraph.getPathLength(self.currPose1, dest, self.localWayPaths[0])
			self.lastDist2 = self.mapGraph.getPathLength(self.currPose2, dest, self.localWayPaths[0])
			
			" pop the next way point if we are already close to the first one "
			if self.lastDist1 < 0.5:
			
				#frontierPoint = [1.07-0.4, 0.0]
				#649, 358
				#frontierPoint = [6.0, -1.5]
				
				" generate the path "
			"""

			#originPath, goalPath, breakPoint = self.mapGraph.computeHeadPath(self.currPose, frontierPoint, self.exploreRoot)
			#self.wayPoints = [breakPoint, goalPath[-1]]
			#self.wayPaths = [originPath, goalPath]

			frontPose = self.contacts.getAverageSegPose(0)
			backPose = self.contacts.getAverageSegPose(39)

			frontPathPoint = self.mapGraph.getNearestPathPoint(frontPose)
			backPathPoint = self.mapGraph.getNearestPathPoint(backPose)
			frontierPathPoint = self.mapGraph.getNearestPathPoint(frontierPoint)
			
			#goalPath1 = self.mapGraph.computeGraphPath(frontPose, frontierPoint)
			#goalPath2 = self.mapGraph.computeGraphPath(backPose, frontierPoint)
			goalPath1 = self.mapGraph.computeGraphPath(frontPathPoint, frontierPathPoint)
			goalPath2 = self.mapGraph.computeGraphPath(backPathPoint, frontierPathPoint)

			#len1 = self.mapGraph.getPathLength(frontPose, frontierPoint, goalPath1)
			#len2 = self.mapGraph.getPathLength(backPose, frontierPoint, goalPath2)
			len1 = self.mapGraph.getPathLength(frontPathPoint, frontierPathPoint, goalPath1)
			len2 = self.mapGraph.getPathLength(backPathPoint, frontierPathPoint, goalPath2)

			if len1 > len2:
				#self.wayPoints = [frontierPoint]
				self.wayPoints = [frontierPathPoint]
				self.wayPaths = [goalPath1]
			else:
				#self.wayPoints = [frontierPoint]
				self.wayPoints = [frontierPathPoint]
				self.wayPaths = [goalPath2]
				
			
			
			print "self.wayPoints =", self.wayPoints
			print "self.wayPaths =", self.wayPaths



			#goalPath = self.mapGraph.computePath(self.currPose, frontierPoint)
			#self.wayPoints = [goalPath[-1]]
			#self.wayPaths = [goalPath]


			
			self.globalState = 11

		elif self.globalState == 11:

			#if self.mapGraph.numNodes > 40:
			#	exit()
			
			" Instantiate the Path Step behavior and give it the path to follow "
		
			" do the path following behavior "
			isDone = self.doPathFollow(self.wayPoints, self.wayPaths)
			
			if isDone:
				if self.termPathID != -1:
					self.mapGraph.poseGraph.paths.pathTermVisited(self.termPathID)
					self.termPathID = -1
				else:
					" invert our pose so we can hit other junctions "
					self.inversion *= -1

				self.restState = deepcopy(probeState)

				if self.isCollided:
					self.travelDir = not self.travelDir
					self.isCollided = False


				self.currPose = self.contacts.getAverageSegPose(0)
				deltaDist = sqrt((self.currPose[0]-self.lastPose1[0])**2 + (self.currPose[1]-self.lastPose1[1])**2)
				print "deltaDist =", deltaDist
				faceDir = True
				self.mapGraph.newNode(faceDir, self.travelDir)
				self.mapGraph.forceUpdate(faceDir)
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				self.globalState = 6
				self.isCapture = True
				
				self.globalState = 6
				self.contacts.setCautious(False)			

		
		elif self.globalState == 12:
			pass

		for i in range(self.numJoints):
			self.probe.setJointTorque(i, self.torques[i])

		sys.stdout.flush()
		
		#if self.globalTimer % 200 == 0:
		#	print self.globalTimer, "setting torques:", self.torques

	def doHoldPosition(self):
		
		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			
			print "START:  doHoldPosition()"
			
			self.behavior = HoldPosition(self.robotParam)
			self.behavior.reset(probeState)
			
			self.localState = 1

		" immediately begin running behavior "
		if self.localState == 1:
			
		
			
			isDone = self.behavior.step(probeState)
			joints1 = self.behavior.getJoints()
			
			self.mergeJoints([joints1])
			
			"isDone always returns False"
			
			if self.globalTimer > 10:
				self.localState = 0
				self.behavior = 0
				
				self.isAnchored = True
				
				print "FINISH:  doHoldPosition()"
				return True
		
		return False		
		
	def doInitAnchor(self):
		
		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			
			print "START:  doInitAnchor()"
			
			self.behavior = AnchorTransition(self.robotParam, self.inversion)
			
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
		
	def doFrontExtend(self):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			print "START:  doFrontExtend()"

			self.lastPose1 = self.contacts.getAveragePose(0)
			self.lastPose2 = self.contacts.getAveragePose(39)
			
			" extend the front "
			self.behavior = FrontExtend(self.robotParam, self.contacts, self.travelDir)
			self.behavior.reset(probeState, self.travelDir)
			
			" hold the back pose anchor "
			self.holdP = HoldPosition(self.robotParam)
			self.holdP.reset(probeState)
			
			self.localState = 1
			self.isAnchored = False
		
		elif self.localState == 1:
	
			" step the main behavior "
			isDone = self.behavior.step(probeState)

			" get and merge the joints for the two behaviors "
			joints1 = self.holdP.getJoints()
			joints2 = self.behavior.getJoints()
			self.mergeJoints([joints2,joints1])
			
			val = self.behavior.getMask()
			self.contacts.setMask(val)
			
 			self.contacts.step(probeState)
				
			if isDone:
				
				self.topJoint = self.behavior.getTopJoint()
				
				self.localState = 0
				self.behavior = 0
				self.holdP = 0
				
				print "FINISH:  doFrontExtend()"
				return True
			
		return False
		
	def doAdaptiveStep(self):

		" get the probe state "
		probeState = self.probe.getProbeState()
				
		if self.localState == 0:
			print "START:  doAdaptiveStep()"
			
			self.behavior = AdaptiveStep(self.robotParam, probeState, self.contacts, self.mapGraph, self.travelDir, self.inversion)
			self.behavior.reset(probeState)
			self.behavior.setTopJoint(self.topJoint)
			
			self.localState = 1
			self.isAnchored = False

		elif self.localState == 1:

			isDone = self.behavior.step(probeState)

			joints2 = self.behavior.getJoints()
			self.mergeJoints([joints2])
			
			val = self.behavior.getMask()
			self.contacts.setMask(val)

			" set the new torques "
			self.torques = [self.behavior.torques[i] for i in range(self.numJoints)]

			for i in range(len(self.torques)):
				if self.torques[i] == None:
					self.torques[i] = self.probe.maxTorque
			
			self.contacts.step(probeState)			
			
			if isDone:
				self.behavior = 0
				self.localState = 0
				
				print "FINISH:  doAdaptiveStep()"
				return True
			
		return False
		
	def doPokeWalls(self, direction):

		" get the probe state "
		probeState = self.probe.getProbeState()
				
		if self.localState == 0:
			print "START:  doPokeWalls()"

			self.behavior = PokeWalls(self.robotParam, self.contacts, direction, self.mapGraph.obstCallBack)
			
			self.localState = 1
			self.isAnchored = True

		elif self.localState == 1:
			
			isDone = self.behavior.step(probeState)
			
			joints2 = self.behavior.getJoints()
			self.mergeJoints([joints2])
			
			" set the new torques "
			self.torques = [self.behavior.torques[i] for i in range(self.numJoints)]

			for i in range(len(self.torques)):
				if self.torques[i] == None:
					self.torques[i] = self.probe.maxTorque

			#val = self.behavior.getMask()
			#self.contacts.setMask(val)
			
 			#self.contacts.step(probeState)
 			
			if self.behavior.hasInitialized():
				self.mapGraph.keepStablePose()
				
				if self.globalTimer % 20 == 0:
					self.mapGraph.update(direction)
	
			if isDone:

				if direction:
					self.mapGraph.currNode.frontProbeError = self.behavior.curlErrors				
				else:
					self.mapGraph.currNode.backProbeError = self.behavior.curlErrors

				
				self.behavior = 0
				self.localState = 0
			
				print "FINISH:  doPokeWalls()"
				return True
			
		return False
	
	def doReturnToRest(self):
		
		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			
			print "START:  doReturnToRest()"
			
			self.behavior = HoldTransition(self.robotParam)
			self.behavior.reset(self.restState)
			
			self.localState = 1
			self.isAnchored = True
	
			self.torques = [self.probe.maxTorque for i in range(self.numJoints)]

						
						
		" immediately begin running behavior "
		if self.localState == 1:
			
			isDone = self.behavior.step(probeState)
			joints1 = self.behavior.getJoints()
			
			self.mergeJoints([joints1])
			
			if isDone:
				self.localState = 0
				self.behavior = 0

				print "FINISH:  doReturnToRest()"
				return True
								
		return False
	
	def doPathStep(self, wayPaths, direction):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
	
			print "START:  doPathStep()", direction

			" instantiate the behavior "
			self.behavior = PathStep(self.robotParam, probeState, self.contacts, self.mapGraph, self.inversion)

			newPaths = deepcopy(wayPaths)
			if not direction:
				newPaths.reverse()
			self.behavior.setPath(newPaths)
			self.behavior.computeCurve()

			" anchoring turned off "
			self.isAnchored = False

			self.localPathDirection = self.behavior.getDirection()
			
			
			self.localState = 1
			
		elif self.localState == 1:
			
			" step the behavior forward "
			isDone = self.behavior.step(probeState)
			
			" get the joints from the behavior "
			joints = self.behavior.getJoints()
			self.mergeJoints([joints])

			" set the new torques "
			self.torques = [self.behavior.torques[i] for i in range(self.numJoints)]
			
			" set the active nodes mask "
			val = self.behavior.getMask()
			self.contacts.setMask(val)
			
			" update the motion model "
			self.contacts.step(probeState)

			if isDone:
				self.localState = 0

				" anchoring turned on "			
				self.isAnchored = True
	
				self.localState = 0
				self.behavior = 0
	
				print "FINISH:  doPathStep()"
				return True
	
		return False

	
	def doPathFollow(self, wayPoints, wayPaths):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" manage the paths "
		if self.localPathState == 0:
			
			print "START:  doPathFollow()"
			self.targetReached = False
			
			self.localDirection = True
			
			self.localWayPoints = deepcopy(wayPoints)
			self.localWayPaths = deepcopy(wayPaths)



			print len(self.localWayPoints), "wayPoints"
			for i in range(len(self.localWayPoints)):
				print self.localWayPoints[i]

			self.drawThings.drawPath(self.localWayPaths[0])
			self.drawThings.drawPoints(self.localWayPoints)
			
			
			
			self.currPose1 = self.contacts.getAverageSegPose(0)
			self.currPose2 = self.contacts.getAverageSegPose(39)

			" determine distance to the next way point "
			dest = self.localWayPoints[0]			

			frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
			backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)
			
			self.lastDist1 = self.mapGraph.getPathLength(frontPathPoint, dest, self.localWayPaths[0])
			self.lastDist2 = self.mapGraph.getPathLength(backPathPoint, dest, self.localWayPaths[0])

			
			" pop the next way point if we are already close to the first one "
			#if self.lastDist1 < 0.5:
			if self.lastDist1 < 0.3:
				self.localWayPoints = self.localWayPoints[1:]
				self.localWayPaths = self.localWayPaths[1:]

				if len(self.localWayPoints) > 0:
	
					self.drawThings.drawPath(self.localWayPaths[0])
					self.drawThings.drawPoints(self.localWayPoints)
					
					" determine distance to the next way point "
					dest = self.localWayPoints[0]			
					#self.lastDist1 = sqrt((dest[0]-self.currPose1[0])**2 + (dest[1]-self.currPose1[1])**2)
					#self.lastDist2 = sqrt((dest[0]-self.currPose2[0])**2 + (dest[1]-self.currPose2[1])**2)
					self.lastDist1 = self.mapGraph.getPathLength(frontPathPoint, dest, self.localWayPaths[0])
					self.lastDist2 = self.mapGraph.getPathLength(backPathPoint, dest, self.localWayPaths[0])

					print "lastDist1: ", self.lastDist1
					print "lastDist2: ", self.lastDist2
				
				else:
					self.targetReached = True
					self.drawThings.drawPath([])
					self.drawThings.drawPoints([])
					
					#if self.localDirection:
					#	self.travelDir = self.localPathDirection
					#else:
					#	self.travelDir = not self.localPathDirection
					
					return True

			" we are not already there, so turn off the collision flag since we are moving "
			self.isCollided = False				
			
			self.distCount = 0
			self.localPathState = 1
			
		elif self.localPathState == 1:

			" do a path step on this path "
			isDone = self.doPathStep(self.localWayPaths[0], self.localDirection)

			" compute distance to the target destination point "
			if self.globalTimer % 10 == 0:
				self.currPose1 = self.contacts.getAverageSegPose(0)
				self.currPose2 = self.contacts.getAverageSegPose(39)
				frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
				backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)
			

				dest = self.localWayPoints[0]
				#dist1 = sqrt((dest[0]-self.currPose1[0])**2 + (dest[1]-self.currPose1[1])**2)
				#dist2 = sqrt((dest[0]-self.currPose2[0])**2 + (dest[1]-self.currPose2[1])**2)
				dist1 = self.mapGraph.getPathLength(frontPathPoint, dest, self.localWayPaths[0])
				dist2 = self.mapGraph.getPathLength(backPathPoint, dest, self.localWayPaths[0])
				
				" check if we've crossed the destination "
				if dist1 < 0.5 or dist2 < 0.5:
					if not self.targetReached:
						print "target reached at wayPoint:", dest, "with dist1,dist2 =", dist1, dist2
					self.targetReached = True
			
			if isDone:
				
				self.restState = deepcopy(probeState)






				self.currPose1 = self.contacts.getAverageSegPose(0)
				self.currPose2 = self.contacts.getAverageSegPose(39)
				frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
				backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)

				lastFrontPathPoint = self.mapGraph.getNearestPathPoint(self.lastPose1)
				lastBackPathPoint = self.mapGraph.getNearestPathPoint(self.lastPose2)

				deltaDist1 = self.mapGraph.getPathLength(frontPathPoint, lastFrontPathPoint, self.localWayPaths[0])
				deltaDist2 = self.mapGraph.getPathLength(backPathPoint, lastBackPathPoint, self.localWayPaths[0])

				#deltaDist1 = sqrt((self.currPose1[0]-self.lastPose1[0])**2 + (self.currPose1[1]-self.lastPose1[1])**2)
				#deltaDist2 = sqrt((self.currPose2[0]-self.lastPose2[0])**2 + (self.currPose2[1]-self.lastPose2[1])**2)
				print "deltaDist1,deltaDist2 =", deltaDist1, deltaDist2


				self.localPathState = 2

				#self.mapGraph.newNode(True, self.localPathDirection)
				#self.mapGraph.newNode(self.stepDist, self.localDirection)
				#self.mapGraph.forceUpdate(True)

				#self.mapGraph.localizeCurrentNode()
	
				#self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				
				self.currPose1 = copy(self.lastPose1)
				self.currPose2 = copy(self.lastPose2)


				if len(self.localWayPoints) > 0:
					

					dest = self.localWayPoints[0]			
					
					path = self.localWayPaths[0]

					self.drawThings.drawPath(path)
					self.drawThings.drawPoints(self.localWayPoints)

					frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
					backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)
					
					totalDist1 = self.mapGraph.getPathLength(frontPathPoint, dest, path)
					totalDist2 = self.mapGraph.getPathLength(backPathPoint, dest, path)
					
					print "distance from", dest, "to", frontPathPoint, "=", totalDist1
					print "distance from", dest, "to", backPathPoint, "=", totalDist2
					
					" better directional detection when more complicated maps "
					" check if we're going away from our target "
					if totalDist1 > self.lastDist1:
						self.distCount += 1
					else:
						self.distCount = 0
					
					print "lastDist1,lastDist2 =", self.lastDist1, self.lastDist2
					print "distCount =", self.distCount
					self.lastDist1 = totalDist1
					self.lastDist2 = totalDist2

					" if we've reached our target after this step, go to next waypoint "
					if self.targetReached:
						
						" if there is another way point, set the path and begin "
						self.localWayPoints = self.localWayPoints[1:]
						self.localWayPaths = self.localWayPaths[1:]
						
						if len(self.localWayPoints) > 0: 
							#self.behavior.setPath(self.localWayPaths[0])
							#self.behavior.computeCurve()
							self.targetReached = False
						else:
							" all points traveled, end behavior "
							self.localPathState = 0
							
							print "FINISH:  doPathFollow()"
							self.drawThings.drawPath([])
							self.drawThings.drawPoints([])

							if self.localDirection:
								self.travelDir = self.localPathDirection
							else:
								self.travelDir = not self.localPathDirection
							return True
							
						#elif self.distCount >= 2:
					elif self.distCount >= 1:
						print "changing direction from", self.localDirection, "to", not self.localDirection
						" if we're going the wrong way, reverse our direction "
						self.localDirection = not self.localDirection
						self.lastDist1 = 1e100
						self.lastDist2 = 1e100
						self.distCount = 0				

				faceDir = True
				self.mapGraph.newNode(faceDir, self.localPathDirection)
				self.mapGraph.forceUpdate(faceDir)
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())

		elif self.localPathState == 2:

			" do a forward poking behavior to sense the environment"
			isDone = self.doPokeWalls(True)
			
			if isDone:
				self.localPathState = 3
			
		elif self.localPathState == 3:

			" return to the rest position "
			isDone = self.doReturnToRest()
			
			if isDone:
				
				self.mapGraph.correctPosture()
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				self.mapGraph.drawConstraints()				
				faceDir = False
				self.mapGraph.newNode(faceDir, self.localPathDirection)
				self.mapGraph.forceUpdate(faceDir)
				
				self.localPathState = 4			
				
		elif self.localPathState == 4:
	
			" do a backward poking behavior to sense the environment"
			isDone = self.doPokeWalls(False)
			
			if isDone:
				self.localPathState = 5
		
		elif self.localPathState == 5:
			
			" return to the rest position "
			isDone = self.doReturnToRest()

			if isDone:



				self.mapGraph.correctPosture()
				self.mapGraph.pairDone()
				self.mapGraph.localizePose()
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)

				self.mapGraph.saveLocalMap()
				self.mapGraph.saveState()
				self.mapGraph.drawConstraints()				


				self.localPathState = 1
				
				
				frontSum = 0.0
				frontProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-2].frontProbeError
				for n in frontProbeError:
					frontSum += n
				foreAvg = frontSum / len(frontProbeError)

				backSum = 0.0
				backProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-1].backProbeError
				for n in backProbeError:
					backSum += n
				backAvg = backSum / len(backProbeError)
							
				print "PathFollow: foreAvg =", foreAvg
				print "PathFollow: backAvg =", backAvg
															
		return False				

		


	def doPathFollow2(self, wayPoints, wayPaths):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" manage the paths "
		if self.localPathState == 0:
			
			print "START:  doPathFollow()"
			self.targetReached = False
			
			self.localDirection = True
			
			self.localWayPoints = deepcopy(wayPoints)
			self.localWayPaths = deepcopy(wayPaths)

			print len(self.localWayPoints), "wayPoints"
			for i in range(len(self.localWayPoints)):
				print self.localWayPoints[i]
			
			
			
			self.currPose1 = self.contacts.getAverageSegPose(0)
			self.currPose2 = self.contacts.getAverageSegPose(39)

			" determine distance to the next way point "
			dest = self.localWayPoints[0]			
			self.lastDist1 = sqrt((dest[0]-self.currPose1[0])**2 + (dest[1]-self.currPose1[1])**2)
			self.lastDist2 = sqrt((dest[0]-self.currPose2[0])**2 + (dest[1]-self.currPose2[1])**2)
			
			" pop the next way point if we are already close to the first one "
			if self.lastDist1 < 0.5 or self.lastDist2 < 0.5:
				self.localWayPoints = self.localWayPoints[1:]
				self.localWayPaths = self.localWayPaths[1:]

				self.drawThings.drawPath(self.localWayPaths[0])
				self.drawThings.drawPoints(self.localWayPoints)
				
				" determine distance to the next way point "
				dest = self.localWayPoints[0]			
				#self.lastDist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
				self.lastDist1 = sqrt((dest[0]-self.currPose1[0])**2 + (dest[1]-self.currPose1[1])**2)
				self.lastDist2 = sqrt((dest[0]-self.currPose2[0])**2 + (dest[1]-self.currPose2[1])**2)
				print "lastDist1, self.lastDist2: ", self.lastDist1, self.lastDist2
				
			self.distCount = 0
			self.localPathState = 1
			
		elif self.localPathState == 1:

			" do a path step on this path "
			isDone = self.doPathStep(self.localWayPaths[0], self.localDirection)

			" compute distance to the target destination point "
			if self.globalTimer % 10 == 0:
				self.currPose1 = self.contacts.getAverageSegPose(0)
				self.currPose2 = self.contacts.getAverageSegPose(39)

			dest = self.localWayPoints[0]
			dist1 = sqrt((dest[0]-self.currPose1[0])**2 + (dest[1]-self.currPose1[1])**2)
			dist2 = sqrt((dest[0]-self.currPose2[0])**2 + (dest[1]-self.currPose2[1])**2)
			
			" check if we've crossed the destination "
			if dist1 < 0.1 or dist2 < 0.1:
				if not self.targetReached:
					print "target reached at wayPoint:", dest
				self.targetReached = True
			
			if isDone:

				#exit()

				#print "Stop Time:", time.clock()
				#exit()
				
				self.restState = deepcopy(probeState)


				self.currPose1 = self.contacts.getAverageSegPose(0)
				self.currPose2 = self.contacts.getAverageSegPose(39)
				deltaDist1 = sqrt((self.currPose1[0]-self.lastPose1[0])**2 + (self.currPose1[1]-self.lastPose1[1])**2)
				deltaDist2 = sqrt((self.currPose2[0]-self.lastPose2[0])**2 + (self.currPose2[1]-self.lastPose2[1])**2)
				print "deltaDist1, deltaDist2 =", deltaDist1, deltaDist2


				self.localPathState = 2

				self.mapGraph.newNode(True, self.localPathDirection)
				#self.mapGraph.newNode(self.stepDist, self.localDirection)
				self.mapGraph.forceUpdate(True)

				#self.mapGraph.localizeCurrentNode()
	
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose = self.contacts.getAverageSegPose(0)
	
						
				if len(self.localWayPoints) > 0:
					

					dest = self.localWayPoints[0]			
					#dist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
					
					indexA = 0
					distA = 1e100
					indexB = 0
					distB = 1e100
					indexC = 0
					distC = 1e100
					path = self.localWayPaths[0]

					self.drawThings.drawPath(path)
					self.drawThings.drawPoints(self.localWayPoints)

					for i in range(len(path)):
						dist = sqrt((dest[0]-path[i][0])**2 + (dest[1]-path[i][1])**2)
						if dist < distA:
							indexA = i
							distA = dist

						dist = sqrt((self.currPose1[0]-path[i][0])**2 + (self.currPose1[1]-path[i][1])**2)
						if dist < distB:
							indexB = i
							distB = dist

						dist = sqrt((self.currPose2[0]-path[i][0])**2 + (self.currPose2[1]-path[i][1])**2)
						if dist < distC:
							indexC = i
							distC = dist
					
					
					pathLen = 0.0
					if indexB >= indexA:
						for i in range(indexA, indexB):
							pathLen += sqrt((path[i][0]-path[i+1][0])**2 + (path[i][1]-path[i+1][1])**2)
					else:
						for i in range(indexB, indexA):
							pathLen += sqrt((path[i][0]-path[i+1][0])**2 + (path[i][1]-path[i+1][1])**2)
					
						
					
					totalDist = distA + distB + pathLen
					
					print "distance from", dest, "to", self.currPose1, "=", totalDist
					
					" FIXME: need better directional detection when more complicated maps "
					" check if we're going away from our target "
					if totalDist > self.lastDist1 and totalDist > self.lastDist2:
						self.distCount += 1
					else:
						self.distCount = 0
					
					print "lastDist1, lastDist2 =", self.lastDist1, self.lastDist2
					print "distCount =", self.distCount
					self.lastDist = totalDist

					" if we've reached our target after this step, go to next waypoint "
					if self.targetReached:
						
						" if there is another way point, set the path and begin "
						self.localWayPoints = self.localWayPoints[1:]
						self.localWayPaths = self.localWayPaths[1:]
						
						if len(self.localWayPoints) > 0: 
							#self.behavior.setPath(self.localWayPaths[0])
							#self.behavior.computeCurve()
							self.targetReached = False
						else:
							" all points traveled, end behavior "
							self.localPathState = 0
							
							print "FINISH:  doPathFollow()"
							self.drawThings.drawPath([])
							self.drawThings.drawPoints([])

							return True
							
						#elif self.distCount >= 2:
					elif self.distCount >= 1:
						print "changing direction from", self.localDirection, "to", not self.localDirection
						" if we're going the wrong way, reverse our direction "
						self.localDirection = not self.localDirection
						self.lastDist = 1e100
						self.distCount = 0				

		elif self.localPathState == 2:

			" do a forward poking behavior to sense the environment"
			isDone = self.doPokeWalls(True)
			
			if isDone:
				self.localPathState = 3
			
		elif self.localPathState == 3:

			" return to the rest position "
			isDone = self.doReturnToRest()
			
			if isDone:
				#self.mapGraph.update()


				#self.mapGraph.correctPosture()
				#self.mapGraph.localizeCurrentNode()
				#self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				#self.lastPose = self.contacts.getAverageSegPose(0)

				#self.mapGraph.synch()
				#self.mapGraph.saveMap()
				#self.mapGraph.saveLocalMap()
				#self.mapGraph.drawConstraints()
				
				#faceDir = False		
				#self.mapGraph.newNode(faceDir, self.travelDir)
				#self.mapGraph.forceUpdate(faceDir)

				
				self.mapGraph.correctPosture()
				#self.mapGraph.localizeCurrentNode()
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose = self.contacts.getAverageSegPose(0)

				#self.mapGraph.synch()
				#self.mapGraph.saveMap()
				self.mapGraph.saveLocalMap()
				self.mapGraph.saveState()
				self.mapGraph.drawConstraints()				
				
				faceDir = False
				self.mapGraph.newNode(faceDir, self.localPathDirection)
				self.mapGraph.forceUpdate(faceDir)

				#self.mapGraph.newNode(self.stepDist, False, self.travelDir)
				#self.mapGraph.currNode.resetPosture()		
				self.localPathState = 4			
				
		elif self.localPathState == 4:
	
			" do a backward poking behavior to sense the environment"
			isDone = self.doPokeWalls(False)
			
			if isDone:
				self.localPathState = 5
		
		elif self.localPathState == 5:
			
			" return to the rest position "
			isDone = self.doReturnToRest()

			if isDone:

				#self.mapGraph.update()
				self.mapGraph.correctPosture()
				#self.mapGraph.localizeCurrentNode()
				
				self.mapGraph.pairDone()
				self.mapGraph.localizePose()
				#self.mapGraph.relaxCorrect()
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose = self.contacts.getAverageSegPose(0)

				#self.mapGraph.synch()
				#self.mapGraph.saveMap()
				self.mapGraph.saveLocalMap()
				self.mapGraph.saveState()
				self.mapGraph.drawConstraints()				



				self.localPathState = 1
				
				
				frontSum = 0.0
				frontProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-2].frontProbeError
				for n in frontProbeError:
					frontSum += n
				foreAvg = frontSum / len(frontProbeError)
							
				print "PathFollow: foreAvg =", foreAvg
															
		return False				