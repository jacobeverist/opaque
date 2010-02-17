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
		self.mapGraph.loadFile(11)


		" behaviors "
		self.anchorT = behave.AnchorTransition(self.probe)
		self.holdP = behave.HoldPosition(self.probe)
		self.holdT = behave.HoldTransition(self.probe)
		self.adaptiveStep = behave.FrontAnchorTest(self.probe, self.contacts, self.mapGraph, direction)

		self.sweep = behave.SpaceSweep(self.probe, direction)
		self.pokeWalls = behave.PokeWalls(self.probe, direction, self.mapGraph.obstCallBack)


		self.stateA = -1
		self.prevTime = 0
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
		
		self.isAnchored = False
		self.anchorR = 0
		self.anchorPos = 0
		size = ogre.Vector3(self.probe.segLength, 0.01, self.probe.segWidth)
		self.anchorMass = OgreOde.BoxMass(STD_WEIGHT*10000.0,size)
		self.normalMass = OgreOde.BoxMass(STD_WEIGHT,size)

		#self.angle_output = open("angle_output_%04u.txt" % 0, "w")
		#self.global_output = open("global_pose_output_%04u.txt" % 0, "w")
		self.isCapture = False
		
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

		pose = self.probe.getActualJointPose(1)
		xAvg = pose[0]
		yAvg = pose[1]
		
		prevPose = self.camera.getPosition()
		xPrev = prevPose[0] + 4
		yPrev = prevPose[2]
		
		newPose = [xPrev*0.99 + 0.01*xAvg, yPrev*0.99 + 0.01*yAvg]

		self.camera.setPosition(newPose[0]-4,7,newPose[1])
		self.camera.lookAt(newPose[0]-2.5,0.5,newPose[1])

	def frameStarted(self):
		
		self.adjustCamera()
		#self.grabImage()
		#self.grabAngles()

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
				self.stateA = 0
				print "going to state 0"

				self.holdP.reset()
				self.holdT.reset()
	
				self.mapGraph.newNode()
				self.mapGraph.forceUpdate(False)
				self.mapGraph.saveLocalMap()
				
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

				self.holdP.reset()
				self.holdT.reset()
				print "going to state 1"

				" set anchoring "
				self.isAnchored = True

				""" create a new pose node since we are in stable anchor position"""
				self.mapGraph.newNode()

				self.isCapture = True

				centerPoints = self.adaptiveStep.getCenterPoints()
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
					
					self.mapGraph.saveLocalMap()
	
					if self.mapGraph.numNodes >= 18:					
						exit()

					self.stateA = 0
					print "going to state 0"
	
					" anchoring turned off "
					self.isAnchored = False

		elif self.stateA == 4:
			
			joints1 = self.holdP.getJoints()
			self.mergeJoints([joints1])
			if self.globalTimer - self.prevTime > 100:
			
				#if self.frontierMap.isFrontier():
				if True:
					self.stateA = 0
					print "going to state 0"
					self.mapGraph.saveLocalMap()

					" anchoring turned off "
					self.isAnchored = False
					#exit()

				else:
					self.stateA = 5
					
					frontierPoint = self.frontierMap.selectNextFrontier()
					currPose = self.probe.getActualSegPose(0)
				
					originPath, goalPath, breakPoint = self.navRoadMap.computeHeadPath(currPose, frontierPoint, self.exploreRoot)
					self.wayPoints = [breakPoint, goalPath[-1]]
					self.wayPaths = [originPath, goalPath]
					
					if self.pathStep == 0:
						self.pathStep = behave.PathConcertinaGait(self.probe, False, self.wayPaths[0])
					else:
						self.pathStep.setPath(self.wayPaths[0])
					
					
		elif self.stateA == 5:
					
			isDone = self.pathStep.step()
			joints = self.pathStep.getJoints()
			self.mergeJoints([joints])
			
			if isDone:
				print "done!"
				
				if len(self.wayPoints) > 0:
					dest = self.wayPoints[0]
					pose = self.probe.getActualSegPose(0)
					
					dist = sqrt((dest[0]-pose[0])**2 + (dest[1]-pose[1])**2)
					print "dist = ", dist
					
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
						else:
							self.stateA = 1
							
					elif self.distCount >= 2:
						self.pathStep.reverseDirection()
						self.lastDist = 1e100
						self.distCount = 0
												
