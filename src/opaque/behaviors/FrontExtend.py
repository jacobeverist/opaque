import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from HoldSlideTransition import HoldSlideTransition
from math import exp
from math import sqrt

"""
FrontExtend behavior

The purpose of this behavior is to extend the front end of the snake as far as it can reasonably go.
Special attention is taken to detect collisions or dead ends and stop the movement when they are detected.
A short retraction is then performed and the current pose is passed to the Front Anchoring behavior.

A variant of this behavior is made that also performs curve fitting for path following purposes.
Its name will be FrontExtendGlobal.

"""

class FrontExtend(Behavior):
	
	def __init__(self, probe, contacts, direction = True):
		Behavior.__init__(self, probe)

		self.contacts = contacts
		
		self.holdSlideT = HoldSlideTransition(self.probe)

		self.setDirection(direction)
				
		self.poses = []
		
		self.maxDist = 0.0
		self.minDist = 10e100
		
	def setDirection(self, isForward):

		self.direction = isForward

		resultJoints = [None for i in range(self.probe.numSegs-1)]

		if self.direction:
			self.spliceJoint = 7

			for i in range(0,25):
				resultJoints[i] = 0.0

			self.topJoint = self.spliceJoint
			
		else:
			self.spliceJoint = 31

			for i in range(15,self.probe.numSegs-1):
				resultJoints[i] = 0.0

			self.topJoint = self.spliceJoint
				
		self.holdSlideT.reset(resultJoints, self.direction)

		self.mask = copy(self.holdSlideT.getMask())

		self.isTransitioning = True
		self.isDone = False
		
	def getTopJoint(self):
		return self.topJoint

	def getSpliceJoint(self, joint):
		return self.spliceJoint

	def getMask(self):
		return self.mask

	def reset(self, joints = [], direction = True):

		self.setDirection(direction)
		
	def step(self):
		Behavior.step(self)
		
		if self.isTransitioning:
			
			if self.direction:
				tipPose = self.contacts.getAveragePose(0)
			else:
				tipPose = self.contacts.getAveragePose(self.probe.numSegs-2)
			
			self.poses.append(tipPose)
			
			points = []
			if len(self.poses) >= 20:
				for i in range(0,20):
					points.append([self.poses[len(self.poses)-i-1][0],self.poses[len(self.poses)-i-1][1]])
			
				spl = SplineFit(points, kp=2)
				
				uvec = spl.getUVector(1.0)
				#print uvec
			
			avgDist = 1e100
			
			if len(self.poses) > 100:
				
				distSum = 0.0
				for i in range(100):
					dist = sqrt((tipPose[0]-self.poses[-2-i][0])**2 + (tipPose[1]-self.poses[-2-i][1])**2)
					distSum += dist
					
				avgDist = distSum / 100.0
				
				#print "avgDist =", avgDist	
				
				
				
				#dist = sqrt((tipPose[0]-self.poses[-2][0])**2 + (tipPose[1]-self.poses[-2][1])**2)
				
				#if dist > self.maxDist:
				#	self.maxDist = dist
				
				#if dist < self.minDist:
				#	self.minDist = dist
			
			self.resetJoints()
			
			# first steps
			isDone = self.holdSlideT.step()
			
			resJoints = self.holdSlideT.getJoints()
			self.mergeJoints([resJoints])
			
			#self.mask = [0.0 for i in range(self.probe.numSegs-1)]
			
			self.mask = copy(self.holdSlideT.getMask())
			
			if self.direction:
				self.topJoint = self.spliceJoint
				joints = range(self.spliceJoint,self.probe.numSegs-1)
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
		
		#self.resetJoints()

		#for i in range(self.probe.numSegs-1):
		#	self.setJoint(i, self.positions[i])

		if self.isDone:
			self.isDone = False
			self.isTransitioning = False

			print len(self.poses), self.maxDist, self.minDist

			#pylab.clf()
			#xP = []
			#yP = []
			#for p in self.poses:
			#	xP.append(p[0])
			#	yP.append(p[1])
			#pylab.plot(xP,yP)
			#pylab.scatter(xP,yP)
			#pylab.show()

			return True

		return False
	
