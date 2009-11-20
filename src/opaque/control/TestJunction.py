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

class TestJunction(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe):
		SnakeControl.__init__(self)
			
		self.probe = probe

		self.setTimerAliasing(1)

		direction = True
		
		" pose estimation "
		self.contacts = pose.ContactReferences(self.probe)
		self.contacts.setTimerAliasing(5)

		" maps "
		self.mapGraph = maps.MapGraph(self.probe, self.contacts)

		" behaviors "
		self.anchorT = behave.AnchorTransition(self.probe)
		self.holdP = behave.HoldPosition(self.probe)
		self.holdT = behave.HoldTransition(self.probe)
		self.adaptiveStep = behave.AdaptiveConcertinaGait(self.probe, self.contacts, self.mapGraph, direction)

		self.sweep = behave.SpaceSweep(self.probe, direction)
		self.pokeWalls = behave.PokeWalls(self.probe, direction, self.mapGraph.obstCallBack)

		self.guidedPoke = behave.GuidedPoke(self.probe, self.mapGraph)

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
	
	def updateMaps(self):
		pass
	
	def saveMaps(self):
		self.mapGraph.saveMap()

	def grabImage(self):

		inc = 50

		if self.globalTimer % inc == 0:
			self.renderWindow.writeContentsToFile("scene%06u.png" % (self.globalTimer/inc))
			#self.saveMaps()

	def adjustCamera(self):
		xAvg = 0
		yAvg = 0
		for i in range(NUM_SEGS-1):
			pose = self.probe.getActualJointPose(0)
			xAvg += pose[0]
			yAvg += pose[1]

		xAvg /= NUM_SEGS-1
		yAvg /= NUM_SEGS-1
		
		prevPose = self.camera.getPosition()
		xPrev = prevPose[0] + 4
		yPrev = prevPose[2]
		
		newPose = [xPrev*0.99 + 0.01*xAvg, yPrev*0.99 + 0.01*yAvg]

		self.camera.setPosition(newPose[0]-4,7,newPose[1])
		self.camera.lookAt(newPose[0]-2,0.5,newPose[1])

	def frameStarted(self):
		
		self.adjustCamera()
		#self.grabImage()

		if self.isAnchored:
			for i in range(16,23):
				body = self.probe._bodies[i]
				body.setMass(self.anchorMass)			
		else:
			for i in range(16,23):
				body = self.probe._bodies[i]
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
				self.guidedPoke.setDirection(True)
				#self.isAnchored = True

				#self.mapGraph.loadFile()
				#self.mapGraph.newNode()
				
				#self.mapGraph.synch()
				#self.mapGraph.saveLocalMap()
				
				#self.mapGraph.saveMap()
				#self.mapGraph.correctPoses2()
				#self.mapGraph.saveMap()
				
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
				#self.mapGraph.saveMap()
			
				" set anchoring "
				self.isAnchored = True
				#body = self.probe._bodies[20]
				#self.anchorR = body.getOrientation()
				#self.anchorPos = body.getPosition()
				
				self.mapGraph.saveLocalMap()
				
				" quit after 15 iterations "
				if self.mapGraph.numNodes >= 15:
					exit()
				
				""" create a new pose node since we are in stable anchor position"""
				self.mapGraph.newNode()

				self.adaptiveStep.computeCenterPoints()
				centerPoints = self.adaptiveStep.getCenterPoints()
				self.mapGraph.setCenterPoints(centerPoints)
				
				self.adaptiveStep.reset()
				

		elif self.stateA == 1:

			self.sweep.setDirection(True)
			isDone = self.sweep.step()
			
			joints1 = self.holdP.getJoints()
			joints2 = self.sweep.getJoints()

			#print "sweep:", joints2
			
			self.mergeJoints([joints1, joints2])

			#self.contacts.setMask( [0.0 for i in range(18)] + [1.0 for i in range(18,39)] )
			#self.contacts.step()

			if self.globalTimer % 30 == 0:
				self.mapGraph.update(True)

			if isDone:
				self.stateA = 2
				print "going to state 2"

		elif self.stateA == 2:
			
			self.pokeWalls.setDirection(True)
			isDone = self.pokeWalls.step()

			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])

			#self.contacts.setMask( [0.0 for i in range(18)] + [1.0 for i in range(18,39)] )
			#self.contacts.step()
			
			if self.globalTimer % 30 == 0:
				self.mapGraph.update(True)
	
			if isDone:
				#self.saveMaps()
				
				self.stateA = 3
				self.prevTime = self.globalTimer
				print "going to state 3"

		elif self.stateA == 3:
			
			isDone = self.holdT.step()
			joints1 = self.holdT.getJoints()

			self.mergeJoints([joints1])

			#self.contacts.setMask( [1.0 for i in range(39)] )
			#self.contacts.step()
			
			if isDone:
				""" do an update forward and backward since we are reversing our sensing direction """
				self.mapGraph.update(True)
				self.stateA = 4
				print "going to state 4"

		elif self.stateA == 4:
			isDone = self.guidedPoke.step()
	
			joints1 = self.anchorT.getJoints()
			joints2 = self.guidedPoke.getJoints()
			self.mergeJoints([joints2, joints1])
			
			if self.globalTimer % 30 == 0:
				self.mapGraph.update(True)
				
			if isDone:
				self.isAnchored = False
				self.stateA = 0
	

