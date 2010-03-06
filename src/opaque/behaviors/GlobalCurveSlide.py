import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
from Transition import Transition
from GlobalCurveFit import GlobalCurveFit
from math import exp

class GlobalCurveSlide(Behavior):
	
	def __init__(self, probe, contacts, curve = 0):
		
		Behavior.__init__(self, probe)

		self.contacts = contacts

		self.setCurve(curve)

		self.isDone = False
		self.count = 4
		self.moveCount = 0

		self.direction = True

		self.targetState = []
		for i in range(self.probe.numSegs-1):
			self.targetState.append(None)

		self.positions = []
		for i in range(self.probe.numSegs-1):
			self.positions.append(180.0 / pi * self.probe.getServoCmd(i))

	def setCurve(self, curve):
		self.curve = curve

		self.globalCurveFit = GlobalCurveFit(self.probe, self.contacts, self.curve)
		self.direction = not self.getPathDirection()
	
	def getPathDirection(self):
		return self.globalCurveFit.getPathDirection()

	def draw(self):
		self.globalCurveFit.draw()

	def reset(self, joints = []):

		print "globalCurveSlide reset:", joints
		
		self.targetState = []
		self.positions = []
		
		#print "resetting HoldSlideTransition to:", joints
		
		if len(joints) > 0:
			for i in range(len(joints)):				
				self.targetState.append(joints[i])
				if joints[i] != None:
					self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
				else:
					self.positions.append(None)
					
		else:
			for i in range(self.probe.numSegs-1):
				self.targetState.append(180.0 / pi * self.probe.getServoCmd(i))
				self.positions.append(180.0 / pi * self.probe.getServoCmd(i))
	
		self.isDone = False
				
	def step(self):
		Behavior.step(self)

		self.resetJoints()


		if self.moveCount == 0:
			if self.direction:
				startNode = 0
				endNode = self.count
				self.globalCurveFit.setBoundaries(startNode, endNode)
			else:
				startNode = endNode - self.count
				endNode = self.probe.numSegs-2	
				self.globalCurveFit.setBoundaries(startNode, endNode)
			
		self.globalCurveFit.step()
		targetState = self.globalCurveFit.getJoints()
		
		print "globalCurveFit:", targetState
		
		self.positions = copy(targetState)

		self.resetJoints()
		self.mergeJoints([targetState])

		self.moveCount += 1
		print "moveCount =", self.moveCount
		
		
		if self.moveCount > 5 :

			self.moveCount = 0
	
			self.count += 1
	
			print "count =", self.count
	
			" delay period to allow body to stabilize "
			if self.count > 24:
				self.count = 4
				
				print "GlobalCurveSlide = DONE"
				return True
			
		
		return False
	
	