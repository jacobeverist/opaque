import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from SnakeControl import *
from copy import *
from math import *
from numpy import arange

import Image

import opaque.behaviors as behave
import opaque.maps as maps
import opaque.pose as pose

class TestMapping(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe):
		SnakeControl.__init__(self)
		
		self.probe = probe
		
		self.setTimerAliasing(1)
		
		direction = True
		
		
		" pose estimation "
		#self.contacts = pose.ContactReferences(self.probe)
		self.contacts = pose.AverageContacts(self.probe)
		self.contacts.setTimerAliasing(5)
		
		" maps "
		self.mapGraph = maps.MapGraph(self.probe, self.contacts)
		self.mapGraph.loadFile(9)
		
		#self.mapGraph.correctPoses2()
		#self.mapGraph.synch()
		#self.mapGraph.saveMap()
		
		" behaviors "
		self.anchorT = behave.AnchorTransition(self.probe)
		self.holdP = behave.HoldPosition(self.probe)
		self.holdT = behave.HoldTransition(self.probe)
		self.adaptiveStep = behave.FrontAnchorTest(self.probe, self.contacts, self.mapGraph, direction)
		
		#self.sweep = behave.SpaceSweep(self.probe, direction)
		self.pokeWalls = behave.PokeWalls(self.probe, direction, self.mapGraph.obstCallBack)
		
		self.exploreRoot = [0.5,0.0]
		self.exploreRoot = [-3.5,0.0]
		self.pathStep = 0
		
		self.stateA = -1
		self.prevTime = 0
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
		self.targetReached = False
		
		self.isAnchored = False
		self.anchorR = 0
		self.anchorPos = 0
		size = ogre.Vector3(self.probe.segLength, 0.01, self.probe.segWidth)
		self.anchorMass = OgreOde.BoxMass(STD_WEIGHT*10000.0,size)
		self.normalMass = OgreOde.BoxMass(STD_WEIGHT,size)

		#self.angle_output = open("angle_output_%04u.txt" % 0, "w")
		#self.global_output = open("global_pose_output_%04u.txt" % 0, "w")
		self.isCapture = False
		
		self.isInitialized = False
		
	def updateMaps(self):
		pass
	
	def saveMaps(self):
		self.mapGraph.saveMap()

	def grabImage(self):

		inc = 250

		if self.globalTimer % inc == 0:
			self.renderWindow.writeContentsToFile("scene%06u.png" % (self.globalTimer/inc))

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
		

	def adjustCamera(self):

		pose = self.probe.getActualJointPose(19)
		xAvg = pose[0]
		yAvg = pose[1]
		
		prevPose = self.camera.getPosition()
		xPrev = prevPose[0] - 3 
		yPrev = prevPose[2]
		zPrev = prevPose[1]
		
		newPose = [xPrev*0.99 + 0.01*xAvg, yPrev*0.99 + 0.01*yAvg, zPrev*0.99 + 0.01*10]

		self.camera.setPosition(newPose[0]+3,newPose[2],newPose[1])
		#self.camera.lookAt(newPose[0],0.5,newPose[1])

		oriQuat = ogre.Quaternion(ogre.Radian(-math.pi/2.0), ogre.Vector3().UNIT_X)
		#oriQuat = ogre.Quaternion(ogre.Radian(-math.pi/2.0), ogre.Vector3().UNIT_X) * ogre.Quaternion(ogre.Radian(-math.pi/2.0), ogre.Vector3().UNIT_Z)
		self.camera.setOrientation(oriQuat)
		
	def frameStarted(self):
		
		self.adjustCamera()
		#self.grabImage()
		#self.grabAngles()

		if self.globalTimer % 4000 == 0:
			self.mapGraph.drawEstBoundary()
		
		if self.isAnchored:
			body = self.probe._bodies[20]
			body.setMass(self.anchorMass)
			
		else:
			body = self.probe._bodies[20]
			body.setMass(self.normalMass)
						
		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		# behavior steps and return values

		# Anchor to the walls

		if self.stateA == -1:
			
			isDone = self.anchorT.step()
			joints1 = self.anchorT.getJoints()
			
			self.mergeJoints([joints1])
			
			self.contacts.setMask( [1.0 for i in range(39)] )
			self.contacts.step()
			
			if isDone:
				#self.stateA = 0
				print "going to state ", self.stateA
				self.stateA = 4

				self.holdP.reset()
				self.holdT.reset()

				self.mapGraph.newNode()
				self.mapGraph.forceUpdate(False)
				
				#self.mapGraph.newNode()
				#self.mapGraph.forceUpdate(False)
				#self.mapGraph.saveLocalMap()
				
				self.mapGraph.synch()
				self.mapGraph.saveMap()
						
				#self.stateA = 0
				#self.stateA = 1
				self.isAnchored = False

				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

		elif self.stateA == 0:
			
			# make a blind concertina step
			isDone = self.adaptiveStep.step()

			joints1 = self.holdP.getJoints()
			joints2 = self.adaptiveStep.getJoints()
			self.mergeJoints([joints2,joints1])
			
			val = self.adaptiveStep.getMask()
			#print "setting", val
			
			self.contacts.setMask(val)
			self.contacts.step()
			
			if isDone:
				
				self.stateA = 1
				#self.stateA = 4

				self.holdP.reset()
				self.holdT.reset()
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
			isDone = self.pokeWalls.step()

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
			
			isDone = self.holdT.step()
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
			isDone = self.pokeWalls.step()
			
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
				self.mapGraph.synch()
				print "check3"
				self.mapGraph.saveLocalMap()
				print "check4"
				self.mapGraph.saveMap()
				print "check5"
				
				#if self.mapGraph.isFrontier():
				if False:
					self.stateA = 0
					print "going to state 0"


					# inhibit frontier points in front of me

					" anchoring turned off "
					self.isAnchored = False
				else:
					self.stateA = 5
					self.targetReached = False

					#self.mapGraph.saveLocalMap()

					" anchoring turned off "
					self.isAnchored = False
					
					print "check6"
					frontierPoint = self.mapGraph.selectNextFrontier()
					#currPose = self.probe.getActualSegPose(0)
					print "check7"

					currPose = self.contacts.getAverageSegPose(0)
					print "check8"
					
					originPath, goalPath, breakPoint = self.mapGraph.computeHeadPath(currPose, frontierPoint, self.exploreRoot)
					self.wayPoints = [breakPoint, goalPath[-1]]
					self.wayPaths = [originPath, goalPath]
					
					#print "originPath:", originPath
					#print "goalPath:", goalPath
					
					#self.wayPaths[0].reverse()

					print "check9"
					self.pathStep = behave.PathStep(self.probe, self.contacts, self.mapGraph, False)
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
								self.holdP.reset()
								self.holdT.reset()

						elif self.distCount >= 2:
							self.pathStep.reverseDirection()
							self.lastDist = 1e100
							self.distCount = 0
					
					print "going to state", self.stateA
					
		elif self.stateA == 5:
			
			isDone = self.pathStep.step()
			joints = self.pathStep.getJoints()
			self.mergeJoints([joints])

			val = self.pathStep.getMask()
			
			self.contacts.setMask(val)
			self.contacts.step()

			
			dest = self.wayPoints[0]
			pose = self.contacts.getAverageSegPose(0)
			dist = sqrt((dest[0]-pose[0])**2 + (dest[1]-pose[1])**2)
			
			if dist < 0.1:
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
				
				self.holdP.reset()
				self.holdT.reset()

				self.isAnchored = True
				print "going to state", self.stateA
				
		elif self.stateA == 6:
			
			self.pokeWalls.setDirection(False)
			isDone = self.pokeWalls.step()

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
			
			isDone = self.holdT.step()
			joints1 = self.holdT.getJoints()

			self.mergeJoints([joints1])

			if isDone:
				self.isCapture = True

				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

				self.stateA = 8
				print "going to state 8"

		elif self.stateA == 8:
			
			self.pokeWalls.setDirection(True)
			isDone = self.pokeWalls.step()

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
			
			isDone = self.holdT.step()
			joints1 = self.holdT.getJoints()

			self.mergeJoints([joints1])

			if isDone:
				self.isCapture = True

				self.mapGraph.correctPoses2()

				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)

				self.mapGraph.synch()
				self.mapGraph.saveLocalMap()
				self.mapGraph.saveMap()

				" anchoring turned off "
				self.isAnchored = False

				self.stateA = 5
				print "going to state 5"
