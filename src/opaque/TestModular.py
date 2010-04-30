from SnakeControl import SnakeControl
from copy import *
from math import *

#from behaviors import AnchorTransition, HoldPosition, HoldTransition, AdaptiveStep, FrontExtend, PokeWalls, PathStep
from behaviors import AnchorTransition, HoldPosition, HoldTransition, FrontExtend, PokeWalls
from pose import AverageContacts		

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
		
		self.robotParam = self.probe.robotParam
		
		self.setTimerAliasing(1)
		
		direction = True
		
		" pose estimation "
		self.contacts = AverageContacts(self.probe)
		self.contacts.setTimerAliasing(5)
		
		" maps "
		#self.mapGraph = maps.MapGraph(self.probe, self.contacts)
		#self.mapGraph.loadFile(9)
		#self.mapGraph.loadFile(1)
		
		#self.mapGraph.correctPoses2()
		#self.mapGraph.synch()
		#self.mapGraph.saveMap()
		
		" behaviors "
		self.anchorT = AnchorTransition(robotParam)
		self.holdP = HoldPosition(robotParam)
		self.holdT = HoldTransition(robotParam)
		#self.adaptiveStep = AdaptiveStep(robotParam, self.contacts, self.mapGraph, direction)
		self.frontExtend = FrontExtend(robotParam, self.contacts, direction)
		self.pokeWalls = PokeWalls(robotParam, direction, self.mapGraph.obstCallBack)
		
		self.exploreRoot = [0.5,0.0]
		self.exploreRoot = [-3.5,0.0]
		self.pathStep = 0
		
		self.stateA = -2
		self.prevTime = 0
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
		self.targetReached = False
		
		self.isAnchored = False

		self.isCapture = False
		self.renderTransition = False
		
		self.isInitialized = False

		self.globalTimer = 1
		
	def updateMaps(self):
		pass
	
	def saveMaps(self):
		self.mapGraph.saveMap()

	def grabImage(self):

		inc = 250

		if self.globalTimer % inc == 0:
			self.drawThings.saveView("scene%06u.png" % (self.globalTimer/inc))

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
		
		#self.grabImage()
		#self.grabAngles()

		if self.globalTimer % 4000 == 0:
			if self.stateA <= 0 or self.stateA == 5 or self.renderTransition:
				self.mapGraph.drawEstBoundary()
				self.renderTransition = False

		self.probe.setAnchor(self.isAnchored)

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		" get the probe state "
		probeState = self.probe.getProbeState()

		# behavior steps and return values

		# Anchor to the walls

		if self.stateA == -2:
			
			isDone = self.anchorT.step(probeState)
			joints1 = self.anchorT.getJoints()
			
			self.mergeJoints([joints1])
			
			self.contacts.setMask( [1.0 for i in range(39)] )
			self.contacts.step(probeState)
			
			if isDone:
				self.stateA = -1
				#self.stateA = 1
				print "going to state ", self.stateA

				self.holdP.reset(probeState)
				self.holdT.reset(probeState)

				self.mapGraph.newNode()
				self.mapGraph.forceUpdate(False)
				
				#self.mapGraph.newNode()
				#self.mapGraph.forceUpdate(False)
				#self.mapGraph.saveLocalMap()
				
				self.mapGraph.synch()
				self.mapGraph.saveMap()
				self.mapGraph.saveLocalMap()
				
				self.isAnchored = False

				self.holdP.reset(probeState)
				self.holdT.reset(probeState)
				
				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

				self.lastPose = self.contacts.getAveragePose(0)
				
				#wall7 = [[-4.8+6.0,-0.2],[-4.5+6.0,0.2]]
				#self.probe.createWall(wall7)

		elif self.stateA == -1:

			# make a blind concertina step
			isDone = self.frontExtend.step(probeState)

			joints1 = self.holdP.getJoints()
			joints2 = self.frontExtend.getJoints()
			self.mergeJoints([joints2,joints1])
			
			val = self.frontExtend.getMask()
			
			self.contacts.setMask(val)
 			self.contacts.step(probeState)
			
			if isDone:
				
				self.adaptiveStep.reset(probeState)
				self.adaptiveStep.setTopJoint(self.frontExtend.getTopJoint())
				self.holdP.reset(probeState)
				self.holdT.reset(probeState)
				self.frontExtend.reset(probeState)
				
				self.stateA = 0
				print "going to state 0"

				
		elif self.stateA == 0:
			
			# make a blind concertina step
			isDone = self.adaptiveStep.step(probeState)

			#joints1 = self.holdP.getJoints()
			joints2 = self.adaptiveStep.getJoints()
			#self.mergeJoints([joints2,joints1])
			self.mergeJoints([joints2])
			
			val = self.adaptiveStep.getMask()

			
			self.contacts.setMask(val)
			self.contacts.step(probeState)
			
			if isDone:
				
				self.stateA = 1
				#self.stateA = 4
				#self.stateA = 3

				self.holdP.reset(probeState)
				self.holdT.reset(probeState)
				print "going to state 1"

				" set anchoring "
				self.isAnchored = True

				""" create a new pose node since we are in stable anchor position"""
				self.mapGraph.newNode()

				self.isCapture = True

				#self.adaptiveStep.computeCenterPoints()
				centerPoints = self.adaptiveStep.getCenterPoints()
				print "centerPoints =", centerPoints
				self.mapGraph.setCenterPoints(centerPoints)
				
				self.mapGraph.forceUpdate(False)
				
		elif self.stateA == 1:
			
			self.pokeWalls.setDirection(False)
			isDone = self.pokeWalls.step(probeState)

			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])

			if self.pokeWalls.hasInitialized():
				self.mapGraph.keepStablePose()
				
				if self.globalTimer % 30 == 0:
					self.mapGraph.update(False)
	
			if isDone:
				self.isCapture = False
				self.stateA = 2
				self.prevTime = self.globalTimer
				print "going to state 2"

		elif self.stateA == 2:
			
			isDone = self.holdT.step(probeState)
			joints1 = self.holdT.getJoints()

			self.mergeJoints([joints1])

			if isDone:
				self.isCapture = True

				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

				self.stateA = 3
				print "going to state 3"

		elif self.stateA == 3:
			
			self.pokeWalls.setDirection(True)
			isDone = self.pokeWalls.step(probeState)
			
			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])
			
			if self.pokeWalls.hasInitialized():
				self.mapGraph.keepStablePose()

				if self.globalTimer % 30 == 0:
					self.mapGraph.update(True)
	
			if isDone:
				self.isCapture = False				

				#self.stateA = 6
				self.prevTime = self.globalTimer
				#print "going to state 6"

				" skip the anchoring and just go to locomotion "
				if True:
					
	
					#if self.mapGraph.numNodes >= 18:					
					#	exit()

					self.stateA = 4
					print "going to state 4"
	
					" anchoring turned off "
					self.isAnchored = False

		elif self.stateA == 4:
			
			joints1 = self.holdP.getJoints()
			self.mergeJoints([joints1])
			if self.globalTimer - self.prevTime > 100:

				print "check1"
				#self.mapGraph.saveLocalMap()
				self.mapGraph.correctPoses2()
				print "check2"
				print "self.contacts.activeRef =", self.contacts.activeRef

				self.mapGraph.synch()
				print "check3"
				print "self.contacts.activeRef =", self.contacts.activeRef
				self.mapGraph.saveLocalMap()
				print "check4"
				print "self.contacts.activeRef =", self.contacts.activeRef
				self.mapGraph.saveMap()
				print "check5"
				print "self.contacts.activeRef =", self.contacts.activeRef

				self.renderTransition = True


				currPose = self.contacts.getAveragePose(0)
				deltaDist = sqrt((currPose[0]-self.lastPose[0])**2 + (currPose[1]-self.lastPose[1])**2)

				print "deltaDist =", deltaDist

				if deltaDist > 0.4:

					#if self.mapGraph.isFrontier():
					#if False:
				

					self.stateA = -1
					print "going to state -1"

					self.lastPose = self.contacts.getAveragePose(0)
					print "self.contacts.activeRef =", self.contacts.activeRef

					# inhibit frontier points in front of me

					" anchoring turned off "
					self.isAnchored = False
				else:
					self.stateA = 5
					self.targetReached = False

					self.lastPose = self.contacts.getAveragePose(0)
					#self.mapGraph.saveLocalMap()

					" anchoring turned off "
					self.isAnchored = False
					
					print "check6"
					frontierPoint = self.mapGraph.selectNextFrontier()
					#currPose = self.probe.getActualSegPose(0)
					print "check7"

					currPose = self.contacts.getAverageSegPose(0)
					print "check8"
					print "self.contacts.activeRef =", self.contacts.activeRef
					
					originPath, goalPath, breakPoint = self.mapGraph.computeHeadPath(currPose, frontierPoint, self.exploreRoot)
					self.wayPoints = [breakPoint, goalPath[-1]]
					self.wayPaths = [originPath, goalPath]
					
					#print "originPath:", originPath
					#print "goalPath:", goalPath
					
					#self.wayPaths[0].reverse()

					print "check9"
					#self.pathStep = PathStep(self.probe, self.contacts, self.mapGraph, False)
					print "check10"
					self.pathStep.setPath(self.wayPaths[0])
					print "check11"
					self.pathStep.computeCurve()
					
					" check if we're already there "
					if len(self.wayPoints) > 0:
						
						dest = self.wayPoints[0]
						pose = self.contacts.getAverageSegPose(0)
						
						dist = sqrt((dest[0]-pose[0])**2 + (dest[1]-pose[1])**2)
						print "distance from", dest, "to", pose, "=", dist
						
						if dist > self.lastDist:
							self.distCount += 1
						else:
							self.distCount = 0
							
						self.lastDist = dist
							
					
						if dist < 0.5:
							self.wayPoints = self.wayPoints[1:]
							self.wayPaths = self.wayPaths[1:]
							if len(self.wayPoints) > 0: 
								#curve = VoronoiFit(self.wayPaths[0])
								#self.pathStep.setCurve(curve)
								self.pathStep.setPath(self.wayPaths[0])
								self.pathStep.computeCurve()
							else:
								self.stateA = 1
								self.mapGraph.newNode()
								self.mapGraph.forceUpdate(False)
								self.isAnchored = True
								
								self.holdP.reset(probeState)
								self.holdT.reset(probeState)

						elif self.distCount >= 2:
							self.pathStep.reverseDirection()
							self.lastDist = 1e100
							self.distCount = 0
					
					print "going to state", self.stateA
					
		elif self.stateA == 5:
	
			isDone = self.pathStep.step(probeState)
			joints = self.pathStep.getJoints()
			self.mergeJoints([joints])

			val = self.pathStep.getMask()
			
			self.contacts.setMask(val)
			self.contacts.step(probeState)
			
			dest = self.wayPoints[0]
			pose = self.contacts.getAverageSegPose(0)
			dist = sqrt((dest[0]-pose[0])**2 + (dest[1]-pose[1])**2)
			
			if dist < 0.1:

				if not self.targetReached:
					print "target reached!"
				self.targetReached = True
						
			if isDone:
				print "done!"
				self.stateA = 6


				self.mapGraph.newNode()
				self.mapGraph.forceUpdate(False)
		
				if len(self.wayPoints) > 0:
					dest = self.wayPoints[0]
					#pose = self.probe.getActualSegPose(0)
					pose = self.contacts.getAverageSegPose(0)
			
					dist = sqrt((dest[0]-pose[0])**2 + (dest[1]-pose[1])**2)
					print "distance from", dest, "to", pose, "=", dist
					
					if dist > self.lastDist:
						self.distCount += 1
					else:
						self.distCount = 0
						
					self.lastDist = dist
					
				
					#if dist < 0.5:
					if self.targetReached:
						self.wayPoints = self.wayPoints[1:]
						self.wayPaths = self.wayPaths[1:]
						if len(self.wayPoints) > 0: 
							#curve = VoronoiFit(self.wayPaths[0])
							#self.pathStep.setCurve(curve)
							self.pathStep.setPath(self.wayPaths[0])
							self.pathStep.computeCurve()
						else:
							self.stateA = 1
							
					elif self.distCount >= 2:
						self.pathStep.reverseDirection()
						self.lastDist = 1e100
						self.distCount = 0
				
				self.holdP.reset(probeState)
				self.holdT.reset(probeState)

				self.isAnchored = True
				print "going to state", self.stateA
				
		elif self.stateA == 6:
			
			self.pokeWalls.setDirection(False)
			isDone = self.pokeWalls.step(probeState)

			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])

			if self.pokeWalls.hasInitialized():
				self.mapGraph.keepStablePose()
				
				if self.globalTimer % 30 == 0:
					self.mapGraph.update(False)
	
			if isDone:
				self.isCapture = False
				self.stateA = 7
				self.prevTime = self.globalTimer
				print "going to state 7"

		elif self.stateA == 7:
			
			isDone = self.holdT.step(probeState)
			joints1 = self.holdT.getJoints()
			self.mergeJoints([joints1])

			if isDone:
				self.isCapture = True

				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

				self.holdT.reset(probeState)
				
				self.stateA = 8
				print "going to state 8"

		elif self.stateA == 8:
			
			self.pokeWalls.setDirection(True)
			isDone = self.pokeWalls.step(probeState)

			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])
				
			if self.pokeWalls.hasInitialized():
				self.mapGraph.keepStablePose()

				if self.globalTimer % 30 == 0:
					self.mapGraph.update(True)
	
			if isDone:
				self.isCapture = False				

				self.prevTime = self.globalTimer

				self.stateA = 9
				print "going to state 9"

		elif self.stateA == 9:
			
			isDone = self.holdT.step(probeState)
			joints1 = self.holdT.getJoints()
			self.mergeJoints([joints1])

			#joints1 = self.holdP.getJoints()
			#self.mergeJoints([joints1])
			#if self.globalTimer - self.prevTime > 100:
				
			if isDone:
				self.isCapture = True

				self.mapGraph.correctPoses2()
				self.renderTransition = True
				
				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

				self.mapGraph.synch()
				self.mapGraph.saveLocalMap()
				self.mapGraph.saveMap()

				" anchoring turned off "
				self.isAnchored = False

				self.stateA = 5
				print "going to state 5"
