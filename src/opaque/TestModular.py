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

class TestModular(SnakeControl):

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
		
		
		self.exploreRoot = [0.5,0.0]
		self.exploreRoot = [-3.5,0.0]
		self.pathStep = 0
		
		self.localPathState = 0
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
		
	def updateMaps(self):
		pass
	
	def saveMaps(self):
		self.mapGraph.saveMap()

	def grabImage(self):

		inc = 50
		#inc = 200
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
	
	def computeCovar(self):
		N = len(self.mapGraph.naives)
		
		for i in range(N):
			x = self.mapGraph.ests[i][1][0]
			y = self.mapGraph.ests[i][1][1]
			p = self.mapGraph.ests[i][1][2]

			xGnd = self.mapGraph.gnds[i][1][0]
			yGnd = self.mapGraph.gnds[i][1][1]
			pGnd = self.mapGraph.gnds[i][1][2]

		avgDist= 0.0
		for i in range(N):
			#print self.mapGraph.gnds[i][0], self.mapGraph.ests[i][0], self.mapGraph.naives[i][0]
			
			avgDist += self.mapGraph.gnds[i][0]
			
		avgDist /= N

		#print "compute mean of gnd offset"
		avgX = 0.0
		avgY = 0.0
		avgPx = 0.0
		avgPy = 0.0
		for i in range(N):
			avgX += self.mapGraph.gnds[i][1][0]
			avgY += self.mapGraph.gnds[i][1][1]
			
			P = self.mapGraph.gnds[i][1][2]
			avgPx += cos(P)
			avgPy += sin(P)
			
		avgX /= N
		avgY /= N
		avgPx /= N
		avgPy /= N
		avgP = atan2(avgPy, avgPx)
		
		#print "average x,y:", avgX, avgY
		#print "average dist:", avgDist
		#print "average angle:", avgP

		stdX = 0.0
		stdY = 0.0
		stdP = 0.0

		Exx = 0.0
		Exy = 0.0
		Exp = 0.0
		Eyp = 0.0
		Eyy = 0.0
		Epp = 0.0
		
		for i in range(N):
			
			x = self.mapGraph.ests[i][1][0]
			y = self.mapGraph.ests[i][1][1]
			p = self.mapGraph.ests[i][1][2]
			
			Exx += (x - avgX) * (x - avgX)
			Exy += (x - avgX) * (y - avgY)
			Exp += (x - avgX) * (p - avgP)
			Eyp += (y - avgY) * (p - avgP)
			Eyy += (y - avgY) * (y - avgY)
			Epp += (p - avgP) * (p - avgP)
			
		#print "pre:", Exx, Eyy, Epp, Exy, Exp, Eyp
		Exx /= N
		Exy /= N
		Exp /= N
		Eyp /= N
		Eyy /= N
		Epp /= N
		
		E = numpy.matrix([[Exx, Exy, Exp],
					[Exy, Eyy, Eyp],
					[Exp, Eyp, Epp]])
		#print "covar:", Exx, Eyy, Epp, Exy, Exp, Eyp	

		return E
	
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
		#	print "torques:", probeState['torques']

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
				self.lastPose = self.contacts.getAveragePose(0)

		elif self.globalState == 3:
			
			" create the mapping object "
			self.mapGraph = MapGraph(self.probe, self.contacts, stabilizePose = True)
			
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
			self.mapGraph.loadFile("testData/constraints", 17)

			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			
			exit()
			#E = self.computeCovar()
			#print repr(E)


			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			#self.mapGraph.saveLocalMap()
	
			#exit()
	
			#self.mapGraph.correctPoses3()
			#exit()
			
			
			self.restState = deepcopy(probeState)
			
			#self.globalState = 6
			self.globalState = 4
			#self.globalState = 9

			#self.restState = deepcopy(probeState)
			#self.mapGraph.newNode(self.stepDist, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.globalState = 6

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

				self.mapGraph.newNode(self.stepDist, self.direction)
				self.mapGraph.forceUpdate(False)
				
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				
				self.lastPose = self.contacts.getAverageSegPose(0)
				
				#self.globalState = 6
				self.globalState = 9

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
				self.globalState = 8
				#if self.mapGraph.currNode.nodeID == 4:
				self.mapGraph.currNode.resetPosture()		
			
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
				self.mapGraph.synch()
				#self.mapGraph.saveMap()
				self.mapGraph.saveLocalMap()
				
		
				self.currPose = self.contacts.getAverageSegPose(0)
				deltaDist = sqrt((self.currPose[0]-self.lastPose[0])**2 + (self.currPose[1]-self.lastPose[1])**2)

				" FORCE to the probe to continue walking forward "
				#deltaDist = 1.0
				print "deltaDist =", deltaDist
	
				self.lastPose = self.currPose
				
				deltaDist = 1.0

				#deltaDist = 0.0

				#if deltaDist > 0.4:
				if deltaDist > 0.2:
					self.globalState = 4
				else:
					self.globalState = 10
					
				#self.globalState = 10
		
		elif self.globalState == 10:
			
			" select the destination and generate path "
			
			" select the point to go to "
			#frontierPoint = self.mapGraph.selectNextFrontier()
			
			frontierPoint = [1.07-0.4, 0.0]
			#649, 358
			#frontierPoint = [6.0, -1.5]
			
			" generate the path "
			originPath, goalPath, breakPoint = self.mapGraph.computeHeadPath(self.currPose, frontierPoint, self.exploreRoot)
			self.wayPoints = [breakPoint, goalPath[-1]]
			self.wayPaths = [originPath, goalPath]
			
			self.globalState = 11

		elif self.globalState == 11:
			
			" Instantiate the Path Step behavior and give it the path to follow "
		
			" do the path following behavior "
			isDone = self.doPathFollow(self.wayPoints, self.wayPaths)
			
			if isDone:
				self.restState = deepcopy(probeState)

				" make sure that at least one blind adaptive step is performed "
				self.lastPose = [-100.0, -100.0]
				
				self.globalState = 6
		
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
		
	def doFrontExtend(self):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			print "START:  doFrontExtend()"

			self.lastPose = self.contacts.getAveragePose(0)
			
			" extend the front "
			self.behavior = FrontExtend(self.robotParam, self.contacts, self.direction)
			self.behavior.reset(probeState, True)
			
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
			
			self.behavior = AdaptiveStep(self.robotParam, probeState, self.contacts, self.mapGraph, self.direction)
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
			
 			self.contacts.step(probeState)
 			
			if self.behavior.hasInitialized():
				self.mapGraph.keepStablePose()
				
				if self.globalTimer % 20 == 0:
					self.mapGraph.update(direction)
	
			if isDone:
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
			self.behavior = PathStep(self.robotParam, probeState, self.contacts, self.mapGraph, direction)

			newPaths = deepcopy(wayPaths)
			if not direction:
				newPaths.reverse()
			self.behavior.setPath(newPaths)
			self.behavior.computeCurve()

			" anchoring turned off "
			self.isAnchored = False
			
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
			
			self.currPose = self.contacts.getAverageSegPose(0)

			" determine distance to the next way point "
			dest = self.localWayPoints[0]			
			self.lastDist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
			
			" pop the next way point if we are already close to the first one "
			if self.lastDist < 0.5:
				self.localWayPoints = self.localWayPoints[1:]
				self.localWayPaths = self.localWayPaths[1:]
				
				" determine distance to the next way point "
				dest = self.localWayPoints[0]			
				self.lastDist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
				
			self.distCount = 0
			self.localPathState = 1
			
		elif self.localPathState == 1:

			" do a path step on this path "
			isDone = self.doPathStep(self.localWayPaths[0], self.localDirection)

			" compute distance to the target destination point "
			if self.globalTimer % 10 == 0:
				self.currPose = self.contacts.getAverageSegPose(0)

			dest = self.localWayPoints[0]
			dist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
			
			" check if we've crossed the destination "
			if dist < 0.1:
				if not self.targetReached:
					print "target reached!"
				self.targetReached = True
			
			if isDone:

				#exit()

				self.restState = deepcopy(probeState)

				self.localPathState = 2

				self.mapGraph.newNode(self.stepDist, self.localDirection)
				self.mapGraph.forceUpdate(False)
					
				if len(self.localWayPoints) > 0:
					

					dest = self.localWayPoints[0]			
					#dist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
					
					indexA = 0
					distA = 1e100
					indexB = 0
					distB = 1e100
					path = self.localWayPaths[0]

					for i in range(len(path)):
						dist = sqrt((dest[0]-path[i][0])**2 + (dest[1]-path[i][1])**2)
						if dist < distA:
							indexA = i
							distA = dist

						dist = sqrt((self.currPose[0]-path[i][0])**2 + (self.currPose[1]-path[i][1])**2)
						if dist < distB:
							indexB = i
							distB = dist
					
					
					pathLen = 0.0
					if indexB >= indexA:
						for i in range(indexA, indexB):
							pathLen += sqrt((path[i][0]-path[i+1][0])**2 + (path[i][1]-path[i+1][1])**2)
					else:
						for i in range(indexB, indexA):
							pathLen += sqrt((path[i][0]-path[i+1][0])**2 + (path[i][1]-path[i+1][1])**2)
					
						
					
					totalDist = distA + distB + pathLen
					
					print "distance from", dest, "to", self.currPose, "=", totalDist
					
					" FIXME: need better directional detection when more complicated maps "
					" check if we're going away from our target "
					if totalDist > self.lastDist:
						self.distCount += 1
					else:
						self.distCount = 0
					
					print "lastDist =", self.lastDist
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

				self.mapGraph.correctPoses2()
				self.mapGraph.synch()
				self.mapGraph.saveLocalMap()
				#self.mapGraph.saveMap()

				self.localPathState = 1
											
		return False				

		
