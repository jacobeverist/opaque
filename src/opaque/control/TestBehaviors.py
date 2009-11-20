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

class TestBehaviors(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe):
		SnakeControl.__init__(self)
			
		self.probe = probe

		self.setTimerAliasing(1)

		direction = False
		self.contacts = pose.ContactReferences(self.probe)
		self.contacts.setTimerAliasing(5)

		# maps
		self.mapGraph = maps.MapGraph(self.probe, self.contacts)

		# behaviors should be
		self.anchor = behave.Anchor(self.probe)
		self.anchorT = behave.AnchorTransition(self.probe)
		self.sweep = behave.SpaceSweep(self.probe, direction)
		self.blindStep = behave.BlindConcertinaGait(self.probe, direction)
		self.transition = behave.Transition(self.probe)
		self.pokeWalls = behave.PokeWalls(self.probe, self.contacts, direction, self.mapGraph.obstCallBack)

		self.stateA = -1
		self.prevTime = 0
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
	
	def updateMaps(self):
		pass
	
	def saveMaps(self):
		self.mapGraph.saveMap()

	def grabImage(self):

		if self.globalTimer % 500 == 0:
			self.renderWindow.writeContentsToFile("scene%06u.png" % (self.globalTimer/500))
			self.saveMaps()

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
		#self.camera.setPosition(xAvg-4,7,yAvg)
		#self.camera.lookAt(xAvg-2,0.5,yAvg)

	def frameStarted(self):
		
		self.adjustCamera()
		self.grabImage()

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		# behavior steps and return values

		# Anchor to the walls
		isStable = self.anchor.step()
		joints1 = self.anchor.getJoints()
		if self.stateA == -1:
			
			isDone = self.anchorT.step()
			joints1 = self.anchorT.getJoints()
			
			self.mergeJoints([joints1])

			self.contacts.setMask( [1.0 for i in range(39)] )
			self.contacts.step()
			
			if isDone:
			#if self.globalTimer > 500:
				self.stateA = 0
				self.mapGraph.newNode()
				#self.mapGraph.loadFile()
				#exit()
	
		elif self.stateA == 0:
			# make a blind concertina step
			isDone = self.blindStep.step()

			joints1 = self.anchor.getJoints()
			joints2 = self.blindStep.getJoints()
			self.mergeJoints([joints2,joints1])
			
			val = self.blindStep.getMask()
			#print "setting", val
			
			self.contacts.setMask(val)
			self.contacts.step()
			
			if isDone:
				self.stateA = 1
				
				self.mapGraph.saveMap()
				
				" quit after 15 iterations "
				if self.mapGraph.numNodes >= 15:
					exit()
				
				""" create a new pose node since we are in stable anchor position"""
				self.mapGraph.newNode()
								
		elif self.stateA == 1:

			self.sweep.setDirection(True)
			isDone = self.sweep.step()
			
			joints1 = self.anchor.getJoints()
			joints2 = self.sweep.getJoints()
			
			
			self.mergeJoints([joints2,joints1])

			self.contacts.setMask( [0.0 for i in range(18)] + [1.0 for i in range(18,39)] )
			#self.contacts.setMask( [1.0 for i in range(0,39)] )
			self.contacts.step()

			if self.globalTimer % 30 == 0:
				self.mapGraph.update(True)

			if isDone:
				self.stateA = 2

		elif self.stateA == 2:
			
			self.pokeWalls.setDirection(True)
			isDone = self.pokeWalls.step()

			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])

			self.contacts.setMask( [0.0 for i in range(18)] + [1.0 for i in range(18,39)] )
			self.contacts.step()
			
			if self.globalTimer % 30 == 0:
				self.mapGraph.update(True)
	
			if isDone:
				#self.saveMaps()
				
				self.stateA = 3
				self.prevTime = self.globalTimer

		elif self.stateA == 3:
			
			isDone = self.anchorT.step()
			joints1 = self.anchorT.getJoints()
			
			self.mergeJoints([joints1])

			self.contacts.setMask( [1.0 for i in range(39)] )
			self.contacts.step()
			
			if isDone:
				""" do an update forward and backward since we are reversing our sensing direction """
				self.mapGraph.update(True)
				self.mapGraph.update(False)
				self.stateA = 4

		elif self.stateA == 4:
			
			self.sweep.setDirection(False)
			isDone = self.sweep.step()
			
			joints1 = self.anchor.getJoints()
			joints2 = self.sweep.getJoints()
			self.mergeJoints([joints2,joints1])

			self.contacts.setMask( [1.0 for i in range(22)] + [0.0 for i in range(22,39)] )
			self.contacts.step()

			if self.globalTimer % 30 == 0:
				self.mapGraph.update(False)

			if isDone:
				self.stateA = 5

		elif self.stateA == 5:
			
			self.pokeWalls.setDirection(False)
			isDone = self.pokeWalls.step()

			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2])

			self.contacts.setMask( [1.0 for i in range(22)] + [0.0 for i in range(22,39)] )
			self.contacts.step()
			
			if self.globalTimer % 30 == 0:
				self.mapGraph.update(False)
	
			if isDone:
				#self.saveMaps()
				
				self.stateA = 6
				self.prevTime = self.globalTimer


		elif self.stateA == 6:
			
			joints1 = self.anchor.getJoints()
			self.mergeJoints([joints1])
			if self.globalTimer - self.prevTime > 100:
			
				#if self.frontierMap.isFrontier():
				if True:
					self.stateA = 0
				else:
					self.stateA = 7
					
					frontierPoint = self.frontierMap.selectNextFrontier()
					currPose = self.probe.getActualSegPose(0)
				
					originPath, goalPath, breakPoint = self.navRoadMap.computeHeadPath(currPose, frontierPoint, self.exploreRoot)
					self.wayPoints = [breakPoint, goalPath[-1]]
					self.wayPaths = [originPath, goalPath]
					
					if self.pathStep == 0:
						self.pathStep = behave.PathConcertinaGait(self.probe, False, self.wayPaths[0])
					else:
						self.pathStep.setPath(self.wayPaths[0])
					
					
		elif self.stateA == 7:
					
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
