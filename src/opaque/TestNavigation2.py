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
from maps.MapUser import MapUser

import numpy
import sys


class TestNavigation2(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe, drawThings, args):
	
		SnakeControl.__init__(self)
		
		self.drawThings = drawThings
		self.probe = probe
		
		robotParam = self.probe.robotParam

		self.args = args

		self.maxNumPoses = self.args.maxNumPoses


		
		self.numJoints = robotParam['numJoints']
		
		self.robotParam = self.probe.robotParam
		
		self.setTimerAliasing(1)
		
		self.travelDir = True
		
		" pose estimation "
		self.contacts = AverageContacts(self.probe)
		self.contacts.setTimerAliasing(5)
		
		" get the probe state "
		probeState = self.probe.getProbeState()
		torques = probeState['torques']

		self.torques = [torques[i] for i in range(self.numJoints)]
		

		self.frontDivDist = 0.0
		self.backDivDist = 0.0

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
		self.lastDist1 = 1e100
		self.distCount = 0
		self.targetReached = False
		self.isCollided = False
		self.collisionCount = 0		
		
		self.isAnchored = False

		self.postureCounter = 0

		self.isCapture = False
		self.renderTransition = False
		
		self.isInitialized = False

		self.globalTimer = 1

		self.termPathID = -1
		self.termPoint = None
		
		self.inversion = 1
		
		
	def updateMaps(self):
		pass
	
	def grabImage(self):

		inc = 200

		if self.globalTimer % inc == 0:

			self.drawThings.render()

			self.drawThings.saveView("scene%06u.png" % (self.globalTimer/inc))
			
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


		if self.globalState > 2 and self.globalTimer % 10 == 0:
			pnts = []
			for i in range(self.numJoints):
				pnts.append(self.contacts.getClosestPose(i))
			self.drawThings.renderPoints(pnts)

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
			if isDone:
				pass
				self.globalState = 1
				
		elif self.globalState == 1:
			
			" create the motion model object "
			self.contacts = AverageContacts(self.probe)
			
			self.globalState = 2
			
		elif self.globalState == 2:

			" get stable reference nodes on the snake body "
			self.contacts.setMask( [1.0 for i in range(self.numJoints)] )	
			self.contacts.step(probeState)

			
			if self.contacts.isStable():
				self.contacts.resetPose(self.probe.getActualJointPose(19))
				self.globalState = 3
				self.lastPose1 = self.contacts.getAveragePose(0)
				self.lastPose2 = self.contacts.getAveragePose(39)

		elif self.globalState == 3:
			
			""" create the mapping object """
			self.mapGraph = MapUser(self.probe, self.contacts, isStable = True, args = self.args)

			#self.mapGraph.restorePickle("eval_success_cross_2015_01_10", 199)
			print "actual joint pose:", self.probe.getActualJointPose(19)

			#self.mapGraph.insertPose(self.probe.getActualJointPose(19), self.travelDir)
			

			self.restState = deepcopy(probeState)
			
			#self.globalState = 6
			self.globalState = 4
			#self.globalState = 10

			self.restState = deepcopy(probeState)

			faceDir = True
			self.mapGraph.newNode(faceDir, self.travelDir)
			self.mapGraph.forceUpdate(faceDir)

			""" capture the front and back tip points of head and tail segments """
			self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
			self.lastPose1 = self.contacts.getAverageSegPose(0)
			self.lastPose2 = self.contacts.getAverageSegPose(39)

			""" enable screen grabs """
			self.isCapture = True

			""" current position of head segment """
			self.currPose1 = self.contacts.getAverageSegPose(0)

			print "Start Time:", time.clock()
			self.mapGraph.drawNavigation([],[])


			
		elif self.globalState == 4:
			
			""" extend the front of the snake """
			isDone = self.doFrontExtend()
			
			if isDone:
				self.globalState = 5

		elif self.globalState == 5:

			""" do a concertina gait step """
			isDone = self.doAdaptiveStep()
			
			if isDone:
				self.restState = deepcopy(probeState)


				""" current position of head segment """
				self.currPose = self.contacts.getAverageSegPose(0)

				""" distance between last head seg position and current head seg position """
				deltaDist = sqrt((self.currPose[0]-self.lastPose1[0])**2 + (self.currPose[1]-self.lastPose1[1])**2)
				print "deltaDist =", deltaDist
	
				""" create forward facing node for front sweep, give current travel direction """
				faceDir = True
				self.mapGraph.newNode(faceDir, self.travelDir)

				""" update with current posture """
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

				self.mapGraph.drawNavigation([],[])
				
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

				#self.mapGraph.synch()
				#self.mapGraph.saveMap()
				#self.mapGraph.saveLocalMap()
				#self.mapGraph.saveState()
				self.mapGraph.drawConstraints()
				
				
				faceDir = False		
				self.mapGraph.newNode(faceDir, self.travelDir)
				self.mapGraph.forceUpdate(faceDir)

				#self.contacts.resetPose(self.mapGraph.getMaxPose())
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)

				self.mapGraph.drawNavigation([],[])

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

				self.mapGraph.correctPosture()
				self.mapGraph.pairDone()


				""" check for collisions so we know when to stop """
				frontSum = 0.0
				print self.mapGraph.numNodes

				frontProbeError = self.mapGraph.foreNode.frontProbeError
				for n in frontProbeError:
					frontSum += n
				foreAvg = frontSum / len(frontProbeError)
				

				backSum = 0.0
				backProbeError = self.mapGraph.backNode.backProbeError
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
						#if self.collisionCount >= 2:
						if self.collisionCount >= 1:
							self.globalState = 10
							self.isCollided = True
							self.collisionCount = 0		
						else:
							if self.collisionCount == 0:
								self.inversion *= -1
							self.globalState = 4
							self.collisionCount += 1		
					else:
						self.globalState = 4
						self.collisionCount = 0		

				else:

					#if backAvg >= 1.4:
					if backAvg >= 1.1:
						#if self.collisionCount >= 2:
						if self.collisionCount >= 1:
							self.globalState = 10
							self.isCollided = True
							self.collisionCount = 0		
						else:
							if self.collisionCount == 0:
								self.inversion *= -1
							self.globalState = 4
							self.collisionCount += 1		
							
					else:
						self.globalState = 4
						self.collisionCount = 0		
					
				print "collisionCount =", self.collisionCount, self.isCollided, self.globalState, foreAvg, backAvg, self.inversion, self.travelDir
				self.mapGraph.localizePose()
				self.contacts.resetPose(self.mapGraph.getMaxPose())
				self.lastPose = self.contacts.getAverageSegPose(0)

				self.mapGraph.saveLocalMap()
				self.mapGraph.saveState()
				self.mapGraph.drawConstraints()				
				
				self.mapGraph.drawNavigation([],[])
				#self.globalState = 10

				#if self.mapGraph.numNodes > 6:
				#	self.travelDir = False	
				#if self.mapGraph.numNodes > 6:
				#	self.globalState = 10

				if self.mapGraph.numNodes > self.maxNumPoses:
					exit()	
					
		elif self.globalState == 10:

			self.contacts.setCautious(True)			
			" select the destination and generate path "
			
			" select the point to go to "
			if self.targetPoint == []:
				self.termPoint, self.termPathID = self.mapGraph.selectNextDestination()
				#self.wayPoints, self.wayPaths = self.mapGraph.recomputePath(self.termPathID, self.termPoint)
				self.wayPoint, self.wayPath = self.mapGraph.computePathSplice(self.termPathID, self.termPoint)

				" determine if we are already at the destination, otherwise recompute "
				dest = self.wayPoint
				wayPath = self.wayPath		

				destReached, goalDist = self.mapGraph.isDestReached(dest, wayPath)
				if destReached:

					self.mapGraph.pathTermVisited(self.termPathID)
					self.termPathID = -1


			else:
				self.termPoint = self.targetPoint
				self.termPathID = -1


			self.termPoint, self.termPathID = self.mapGraph.selectNextDestination()
			#self.wayPoints, self.wayPaths = self.mapGraph.recomputePath(self.termPathID, self.termPoint)
			self.wayPoint, self.wayPath = self.mapGraph.computePathSplice(self.termPathID, self.termPoint)


							
			#print "self.wayPoints =", self.wayPoints
			#print "self.wayPaths =", self.wayPaths
			print "self.wayPoint =", self.wayPoint
			print "self.wayPath =", self.wayPath

			self.mapGraph.drawNavigation(self.wayPath,self.wayPoint)

			self.globalState = 11

		elif self.globalState == 11:

			if self.mapGraph.numNodes > self.maxNumPoses:
				exit()
			
			" Instantiate the Path Step behavior and give it the path to follow "
		
			" do the path following behavior "
			#isDone = self.doPathFollow(self.wayPoints, self.wayPaths)
			isDone = self.doPathFollow(self.wayPoint, self.wayPath)
			
			if isDone:
				if self.termPathID != -1:
					self.mapGraph.pathTermVisited(self.termPathID)
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
				#self.contacts.resetPose(self.mapGraph.getMaxPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				self.globalState = 6
				self.isCapture = True
				
				self.globalState = 6
				self.contacts.setCautious(False)			

				self.mapGraph.drawNavigation([],[])

		
		elif self.globalState == 12:
			pass

		for i in range(self.numJoints):
			self.probe.setJointTorque(i, self.torques[i])

		#print self.probe.getProbeState()['torques']
		#print self.torques	
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

			""" send data into the posture map, filtered by the stablePose process """
			if self.behavior.hasInitialized():
				self.mapGraph.keepStablePose()
				
				if self.globalTimer % 20 == 0:
					self.mapGraph.update(direction)
	
			if isDone:

				""" curl errors used for detecting collisions """
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
	
	def doPathStep(self, wayPath, direction):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
	
			print "START:  doPathStep()", direction

			" instantiate the behavior "
			self.behavior = PathStep(self.robotParam, probeState, self.contacts, self.mapGraph, self.inversion)

			newPath = deepcopy(wayPath)
			if not direction:
				newPath.reverse()
			self.behavior.setPath(newPath)
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

	
	def doPathFollow(self, wayPoint, wayPath):

		" get the probe state "
		probeState = self.probe.getProbeState()

		" manage the paths "
		if self.localPathState == 0:
			
			print "START:  doPathFollow()"
			self.targetReached = False
			
			self.localDirection = True

			targetDirection = self.mapGraph.getDirection(wayPath)

			
			self.localWayPoint = deepcopy(wayPoint)
			self.localWayPath = deepcopy(wayPath)



			print "wayPoint", self.localWayPoint

			self.drawThings.drawPath(self.localWayPath)
			self.drawThings.drawPoints(self.localWayPoint)
			
			
			self.currPose1 = self.contacts.getAverageSegPose(0)
			self.currPose2 = self.contacts.getAverageSegPose(39)

			" determine distance to the next way point "
			dest = self.localWayPoint
			wayPath = self.localWayPath		

			""" get destination reached status """
			destReached, goalDist = self.mapGraph.isDestReached(dest, wayPath)

			" - get splice of GPAC, is it the same as wayPaths? "
			" - if so, get closest points on wayPath "
			" - have either of the tips passed the waypoint? "

			frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
			backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)

			dist1 = self.mapGraph.getPathLength(frontPathPoint, dest, wayPath)
			dist2 = self.mapGraph.getPathLength(backPathPoint, dest, wayPath)

			self.lastFrontDivDist = 0.0
			self.lastBackDivDist = 0.0

			#self.frontDivDist = self.mapGraph.getDistanceToPath(self.currPose1, wayPath)
			#self.backDivDist = self.mapGraph.getDistanceToPath(self.currPose2, wayPath)
			self.frontDivDist, self.backDivDist, self.contigFrac = self.mapGraph.isFollowing(wayPath)

			self.backtrackCount = 0
			self.isBacktrack = False

			if targetDirection:
				if self.frontDivDist > 0.2 and self.contigFrac < 0.5:
					self.localDirection = False
					self.isBacktrack = True
			else:
				if self.backDivDist > 0.2 and self.contigFrac < 0.5:
					self.localDirection = False
					self.isBacktrack = True

			#if dist1 < dist2:
			#	self.localDirection = False

			#	if self.frontDivDist > 0.2:
			#		self.localDirection = True
			#		self.isBacktrack = True
			
			#else:
			#	self.localDirection = True

			#	if self.backDivDist > 0.2:
			#		self.localDirection = False
			#		self.isBacktrack = True

			print "path divergence distance:", self.frontDivDist, self.backDivDist, self.lastFrontDivDist, self.lastBackDivDist, dist1, dist2, self.isBacktrack, self.localDirection, targetDirection, self.contigFrac

			self.lastDist1 = goalDist
			
			" pop the next way point if we are already close to the first one "
			if destReached:
				self.localWayPoint = None
				self.localWayPath = None

				self.targetReached = True
				self.drawThings.drawPath([])
				self.drawThings.drawPoints([])
					
				return True

			" we are not already there, so turn off the collision flag since we are moving "
			self.isCollided = False				
			
			self.distCount = 0
			self.localPathState = 1

			
		elif self.localPathState == 1:

			" do a path step on this path "
			isDone = self.doPathStep(self.localWayPath, self.localDirection)

			" compute distance to the target destination point "
			if self.globalTimer % 10 == 0:
				self.currPose1 = self.contacts.getAverageSegPose(0)
				self.currPose2 = self.contacts.getAverageSegPose(39)
				frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
				backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)

				dest = self.localWayPoint
				wayPath = self.localWayPath		
				
				dist1 = self.mapGraph.getPathLength(frontPathPoint, dest, wayPath)
				dist2 = self.mapGraph.getPathLength(backPathPoint, dest, wayPath)
				
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

				deltaDist1 = self.mapGraph.getPathLength(frontPathPoint, lastFrontPathPoint, self.localWayPath)
				deltaDist2 = self.mapGraph.getPathLength(backPathPoint, lastBackPathPoint, self.localWayPath)
				print "deltaDist1,deltaDist2 =", deltaDist1, deltaDist2

				self.localPathState = 2

				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				
				self.currPose1 = copy(self.lastPose1)
				self.currPose2 = copy(self.lastPose2)

				self.lastFrontDivDist = self.frontDivDist
				self.lastBackDivDist = self.backDivDist

				faceDir = True
				self.mapGraph.newNode(faceDir, self.localPathDirection)
				self.mapGraph.forceUpdate(faceDir)
				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())

				self.mapGraph.drawNavigation(self.localWayPath,self.localWayPoint)
				
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				
				self.currPose1 = copy(self.lastPose1)
				self.currPose2 = copy(self.lastPose2)

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
				self.mapGraph.drawConstraints()				
				faceDir = False
				self.mapGraph.newNode(faceDir, self.localPathDirection)
				self.mapGraph.forceUpdate(faceDir)

				self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)
				
				self.mapGraph.drawNavigation(self.localWayPath,self.localWayPoint)

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


				self.localPathState = 1
				
				
				frontSum = 0.0
				#frontProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-2].frontProbeError
				frontProbeError = self.mapGraph.foreNode.frontProbeError
				for n in frontProbeError:
					frontSum += n
				foreAvg = frontSum / len(frontProbeError)

				backSum = 0.0
				#backProbeError = self.mapGraph.nodeHash[self.mapGraph.numNodes-1].backProbeError
				backProbeError = self.mapGraph.backNode.backProbeError
				for n in backProbeError:
					backSum += n
				backAvg = backSum / len(backProbeError)
							
				print "PathFollow: foreAvg =", foreAvg
				print "PathFollow: backAvg =", backAvg
															
				self.mapGraph.localizePose()
				self.contacts.resetPose(self.mapGraph.getMaxPose())
				#self.contacts.resetPose(self.mapGraph.currNode.getEstPose())
				self.lastPose1 = self.contacts.getAverageSegPose(0)
				self.lastPose2 = self.contacts.getAverageSegPose(39)

				self.mapGraph.saveLocalMap()
				self.mapGraph.saveState()
				self.mapGraph.drawConstraints()				




				""" map state hypothesis no longer active, terminate path following """
				#if self.mapGraph.mapAlgorithm.activeHypID == None:
				#	self.mapGraph.drawNavigation([],[])
				#	return True

				try:
					""" check where we are, what we need to do """
					print "respliced", self.localWayPoint
					#self.localWayPoint, self.localWayPath = self.mapGraph.recomputePath(self.termPathID, self.termPoint)
					self.localWayPoint, self.localWayPath = self.mapGraph.computePathSplice(self.termPathID, self.termPoint)

					print "respliced after", self.localWayPoint

							
					if self.localWayPoint != None:


						print "case foo A"

						dest = self.localWayPoint
						
						path = self.localWayPath

						self.drawThings.drawPath(path)
						self.drawThings.drawPoints([self.localWayPoint,])

						print "case foo B"
						destReached, goalDist = self.mapGraph.isDestReached(dest,path, foreAvg, backAvg)

						print "case foo C"

						self.frontDivDist, self.backDivDist, self.contigFrac = self.mapGraph.isFollowing(path)
						print "case foo D"

						print "path divergence distance:", self.frontDivDist, self.backDivDist, self.lastFrontDivDist, self.lastBackDivDist, self.isBacktrack, self.localDirection, self.contigFrac


						frontPathPoint = self.mapGraph.getNearestPathPoint(self.currPose1)
						backPathPoint = self.mapGraph.getNearestPathPoint(self.currPose2)
						
						print "distance from", dest, "=", goalDist
						
						" better directional detection when more complicated maps "
						" check if we're going away from our target "
						if goalDist > self.lastDist1:
							self.distCount += 1
						else:
							self.distCount = 0
						
						print "lastDist1 =", self.lastDist1, 
						print "distCount =", self.distCount
						self.lastDist1 = goalDist


						if self.isBacktrack:
							self.backtrackCount += 1
							#if self.backtrackCount > 1:
							if self.frontDivDist <= 0.2 and self.backDivDist <= 0.2 or self.contigFrac > 0.5:
								print "stopping backtracking:", self.frontDivDist, self.backDivDist
								self.isBacktrack = False
								self.localDirection = not self.localDirection
								self.backtrackCount

						elif (self.frontDivDist > 0.2 or self.backDivDist > 0.2) and self.contigFrac < 0.5:
							print "backtracking:", self.frontDivDist, self.backDivDist
							self.isBacktrack = True
							self.localDirection = not self.localDirection
							
							" if we've reached our target after this step, go to next waypoint "

						elif destReached:
							
							" if there is another way point, set the path and begin "
							self.localWayPoint = None
							self.localWayPath = None
							
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
								

					self.contacts.resetPose(self.mapGraph.currNode.getEstPose())

					self.mapGraph.drawNavigation(self.localWayPath,self.localWayPoint)
					
					self.lastPose1 = self.contacts.getAverageSegPose(0)
					self.lastPose2 = self.contacts.getAverageSegPose(39)
					
					self.currPose1 = copy(self.lastPose1)
					self.currPose2 = copy(self.lastPose2)



				except:
					""" we lost our destination from the map """
					self.mapGraph.drawNavigation(self.localWayPath, self.localWayPoint)
					self.targetReached = True
					self.termPathID = -1
					#raise
					return True

				self.mapGraph.drawNavigation(self.localWayPath, self.localWayPoint)

				if self.mapGraph.numNodes > self.maxNumPoses:
					exit()	

		return False				

		
	
