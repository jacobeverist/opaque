from SnakeControl import SnakeControl
from copy import *
from math import *

from behaviors.AnchorTransition import AnchorTransition
from behaviors.HoldTransition import HoldTransition
from behaviors.HoldPosition import HoldPosition
from behaviors.FrontExtend import FrontExtend
from behaviors.PokeWalls import PokeWalls
from behaviors.PathStep import PathStep
from behaviors.AdaptiveStep import AdaptiveStep

from pose.AverageContacts import AverageContacts
from maps.MapGraph import MapGraph

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

		inc = 1000

		if self.globalTimer % inc == 0:
			pass
			self.drawThings.saveView("scene%06u.png" % (self.globalTimer/inc))
			
			poses = self.probe.getPose()		
			f = open("poses%06u.txt" % (self.globalTimer/inc), 'w')
			f.write(repr(poses))
			f.close()

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
		
		self.grabImage()
		#self.grabAngles()

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


		#print probeState['joints']

		if self.globalState == 0 and self.globalTimer > 2:
			
			" anchor the robot initially "

			isDone = self.doInitAnchor()
			#isDone = self.doHoldPosition()
			
			if isDone:
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
				self.globalState = 3
				self.lastPose = self.contacts.getAveragePose(0)

		elif self.globalState == 3:
			
			" create the mapping object "
			self.mapGraph = MapGraph(self.probe, self.contacts)
			#self.mapGraph.loadFile("testData/correctionTest", 6)

			self.mapGraph.newNode()
			self.mapGraph.forceUpdate(False)

			self.mapGraph.synch()
			self.mapGraph.saveMap()
			self.mapGraph.saveLocalMap()
			
			self.restState = deepcopy(probeState)
			
			#self.globalState = 6
			#self.globalState = 4
			self.globalState = 9
			
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

				self.mapGraph.newNode()
				self.mapGraph.forceUpdate(False)
				
				self.globalState = 6
				
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
			
		elif self.globalState == 8:
	
			" do a backward poking behavior to sense the environment"
			isDone = self.doPokeWalls(False)
			
			if isDone:
				self.globalState = 9
		
		elif self.globalState == 9:
			
			" return to the rest position "
			isDone = self.doReturnToRest()
			
			if isDone:
				self.mapGraph.correctPoses2()
				self.mapGraph.synch()
				self.mapGraph.saveMap()
				self.mapGraph.saveLocalMap()
				
		
				self.currPose = self.contacts.getAveragePose(0)
				deltaDist = sqrt((self.currPose[0]-self.lastPose[0])**2 + (self.currPose[1]-self.lastPose[1])**2)

				print "deltaDist =", deltaDist
	
				self.lastPose = self.currPose
				
				#deltaDist = 0.0

				if deltaDist > 0.4:
					self.globalState = 4
				else:
					self.globalState = 10
		
		elif self.globalState == 10:
			
			" select the destination and generate path "
			
			" select the point to go to "
			frontierPoint = self.mapGraph.selectNextFrontier()
			
			#frontierPoint = [1.07-0.4, 0.0]
			
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
			
			#val = self.behavior.getMask()
			#self.contacts.setMask(val)
			
 			self.contacts.step(probeState)
 			
			if self.behavior.hasInitialized():
				self.mapGraph.keepStablePose()
				
				if self.globalTimer % 200 == 0:
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
				
			self.currPose = self.contacts.getAverageSegPose(0)
			self.distCount = 0
			self.localPathState = 1
			
		elif self.localPathState == 1:

			" do a path step on this path "
			isDone = self.doPathStep(self.localWayPaths[0], self.localDirection)

			" compute distance to the target destination point "
			if self.globalTimer % 100 == 0:
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

				self.mapGraph.newNode()
				self.mapGraph.forceUpdate(False)
		
				if len(self.localWayPoints) > 0:
					
					dest = self.localWayPoints[0]			
					dist = sqrt((dest[0]-self.currPose[0])**2 + (dest[1]-self.currPose[1])**2)
					print "distance from", dest, "to", self.currPose, "=", dist
					
					" FIXME: need better directional detection when more complicated maps "
					" check if we're going away from our target "
					if dist > self.lastDist:
						self.distCount += 1
					else:
						self.distCount = 0
					
					print "lastDist =", self.lastDist
					print "distCount =", self.distCount
					self.lastDist = dist

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
							
					elif self.distCount >= 2:
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
				self.mapGraph.saveMap()

				self.localPathState = 1
											
		return False				

		
