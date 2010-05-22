from Behavior import *
from HoldSlideTransition import HoldSlideTransition
from HoldPosition import HoldPosition
from math import sqrt
from copy import copy

"""
FrontExtend behavior

The purpose of this behavior is to extend the front end of the snake as far as it can reasonably go.
Special attention is taken to detect collisions or dead ends and stop the movement when they are detected.
A short retraction is then performed and the current pose is passed to the Front Anchoring behavior.

A variant of this behavior is made that also performs curve fitting for path following purposes.
Its name will be FrontExtendGlobal.

"""

class FrontExtend(Behavior):
	
	def __init__(self, robotParam, contacts, direction = True):
		Behavior.__init__(self, robotParam)

		self.contacts = contacts
		
		self.holdP = HoldPosition(robotParam)
		self.holdSlideT = HoldSlideTransition(robotParam)

		#self.setDirection(direction)
				
		self.poses = []
	
	def setDirection(self, probeState, isForward):

		self.direction = isForward

		resultJoints = [None for i in range(self.numJoints)]

		if self.direction:
			self.spliceJoint = 7

			for i in range(0,25):
				resultJoints[i] = 0.0

			self.topJoint = self.spliceJoint
			
		else:
			self.spliceJoint = 31

			for i in range(15,self.numJoints):
				resultJoints[i] = 0.0

			self.topJoint = self.spliceJoint
				
		self.holdSlideT.reset(probeState, resultJoints, self.direction)

		self.mask = copy(self.holdSlideT.getMask())

		self.isTransitioning = True
		self.isDone = False
		
	def getTopJoint(self):
		return self.topJoint

	def getSpliceJoint(self, joint):
		return self.spliceJoint

	def getMask(self):
		return self.mask

	def reset(self, probeState, direction = True):
		self.holdSlideT = HoldSlideTransition(self.robotParam)
		
		self.holdP.reset(probeState)
		
		self.setDirection(probeState, direction)
				
		self.poses = []
		
	def step(self, probeState):
		Behavior.step(self, probeState)
		
		if self.isTransitioning:
			
			if self.direction:
				tipPose = self.contacts.getAveragePose(0)
			else:
				tipPose = self.contacts.getAveragePose(self.numJoints-1)
			
			self.poses.append(tipPose)
			
			points = []
			if len(self.poses) >= 20:
				for i in range(0,20):
					points.append([self.poses[len(self.poses)-i-1][0],self.poses[len(self.poses)-i-1][1]])
			
				#spl = SplineFit(points, kp=2)
				
				#uvec = spl.getUVector(1.0)
				#print uvec
			
			avgDist = 1e100
			
			if len(self.poses) > 100:
				
				distSum = 0.0
				for i in range(100):
					dist = sqrt((tipPose[0]-self.poses[-2-i][0])**2 + (tipPose[1]-self.poses[-2-i][1])**2)
					distSum += dist
					
				avgDist = distSum / 100.0
			
			self.resetJoints()
			
			# first steps
			isDone = self.holdSlideT.step(probeState)
			
			resJoints = self.holdSlideT.getJoints()
			self.mergeJoints([resJoints])
			
			self.mask = copy(self.holdSlideT.getMask())
			
			if self.direction:
				self.topJoint = self.spliceJoint
				joints = range(self.spliceJoint,self.numJoints)
				for i in joints:
					if self.mask[i] == 0.0:
						self.topJoint = i
					else:
						break
			else:
				self.topJoint = self.spliceJoint
				joints = range(0, self.spliceJoint+1)
				joints.reverse()
				for i in joints:
					if self.mask[i] == 0.0:
						self.topJoint = i
					else:
						break

			if avgDist < 0.001:
				print "detected obstruction!", avgDist, "with topJoint =", self.topJoint
				isDone = True

			if isDone:
				self.isDone = True
				self.isTransitioning = False
			
			return False
		
		if not self.isStep():
			return False
		
		if self.isDone:
			self.isDone = False
			self.isTransitioning = False

			return True

		return False
	
