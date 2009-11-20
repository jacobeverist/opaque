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

class BuildMap(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe):
		SnakeControl.__init__(self)
		
		self.probe = probe

		# Components of control
		# 1. pose estimation
		# 2. maps
		# 3. behaviors

		direction = False
		joints = range(0, 20)

		# behaviors should be
		self.anchor = behave.Anchor(self.probe)
		self.sweep = behave.SpaceSweep(self.probe, direction)
		self.blindStep = behave.BlindConcertinaGait(self.probe, direction)
		
		self.setTimerAliasing(1)

		self.occMap = maps.OccupancyMap(probe, 19)
		
		# optionally initialize with an existing occupancy map
		#initImage = Image.open("startMap.png")
		#self.occMap.loadImage(initImage)

		self.boundaryMap = maps.FreeSpaceBoundaryMap(self.occMap)
		self.voronoiMap = maps.VoronoiMap(self.boundaryMap.getBoundaryPoints(), self.occMap.getImage())
		self.navRoadMap = maps.NavRoadMap(self.probe, self.voronoiMap.getGraph())
		self.obstacleMap = maps.ObstacleMap(self.probe, self.boundaryMap)
		self.frontierMap = maps.FrontierMap(self.probe, self.boundaryMap, self.obstacleMap)
		self.obstacleMap.setFrontierMap(self.frontierMap)
		
		self.pokeWalls = behave.PokeWalls(self.probe, direction, self.obstacleMap)
	
		self.exploreRoot = [-8.0,0.0]

		self.pathStep = 0
		self.stateA = 1
		self.topJoint =4
		self.prevTime = 0
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
		
	def updateMaps(self):
		self.boundaryMap.update(self.occMap)
		self.voronoiMap.update(self.boundaryMap.getBoundaryPoints(), self.occMap.getImage())
		self.obstacleMap.update(self.boundaryMap)
		self.frontierMap.update(self.boundaryMap, self.obstacleMap)
		self.navRoadMap.update(self.voronoiMap.getGraph())

	def saveMaps(self):
		self.occMap.saveMap()
		self.boundaryMap.saveMap()
		self.voronoiMap.saveMap()
		self.obstacleMap.saveMap()
		self.frontierMap.saveMap()
	
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
		
		if self.stateA == 0:
			# make a blind concertina step
			isDone = self.blindStep.step()

			joints1 = self.anchor.getJoints()
			joints2 = self.blindStep.getJoints()
			self.mergeJoints([joints2,joints1])
				
			if isDone:
				self.stateA = 1
				
		elif self.stateA == 1:
			isDone = self.sweep.step()
			
			joints1 = self.anchor.getJoints()
			joints2 = self.sweep.getJoints()
			self.mergeJoints([joints2,joints1])
			
			if self.globalTimer % 30 == 0:
				self.occMap.update()

			if isDone:
				self.stateA = 2
				self.updateMaps()
			
				
		elif self.stateA == 2:
			
			isDone = self.pokeWalls.step()

			joints1 = self.anchor.getJoints()
			joints2 = self.pokeWalls.getJoints()
			self.mergeJoints([joints2, joints1])
			
			if self.globalTimer % 30 == 0:
				self.occMap.update()
	
			if isDone:
				self.updateMaps()
				
				self.stateA = 3
				self.prevTime = self.globalTimer


		elif self.stateA == 3:
			
			joints1 = self.anchor.getJoints()
			self.mergeJoints([joints1])
			if self.globalTimer - self.prevTime > 100:
			
				if self.frontierMap.isFrontier():
					self.stateA = 0
				else:
					self.stateA = 4
					
					frontierPoint = self.frontierMap.selectNextFrontier()
					currPose = self.probe.getActualSegPose(0)
				
					originPath, goalPath, breakPoint = self.navRoadMap.computeHeadPath(currPose, frontierPoint, self.exploreRoot)
					self.wayPoints = [breakPoint, goalPath[-1]]
					self.wayPaths = [originPath, goalPath]
					
					if self.pathStep == 0:
						self.pathStep = behave.PathConcertinaGait(self.probe, False, self.wayPaths[0])
					else:
						self.pathStep.setPath(self.wayPaths[0])
					
					
		elif self.stateA == 4:
					
			isDone = self.pathStep.step()
			joints = self.pathStep.getJoints()
			self.mergeJoints([joints])
			
			if isDone:
				
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
							self.pathStep.setPath(self.wayPaths[0])
						else:
							self.stateA = 1
							
					elif self.distCount >= 2:
						self.pathStep.reverseDirection()
						self.lastDist = 1e100
						self.distCount = 0

		