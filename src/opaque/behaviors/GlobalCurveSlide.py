
from Behavior import Behavior
from Transition import Transition
from GlobalCurveFit import GlobalCurveFit
from math import exp, pi
from copy import copy

class GlobalCurveSlide(Behavior):
	
	def __init__(self, robotParam, contacts, curve = 0, localNode = 0):
		
		Behavior.__init__(self, robotParam)

		self.contacts = contacts
		self.currNode = localNode
		self.direction = True
		self.numJoints = robotParam['numJoints']

		self.setCurve(curve)

		self.isDone = False
		self.count = 4
		self.moveCount = 0


		self.targetState = []
		for i in range(self.numJoints):
			self.targetState.append(None)

		self.positions = []
		for i in range(self.numJoints):
			self.positions.append(None)
			#self.positions.append(180.0 / pi * cmdJoints[i])

	def setCurve(self, curve):
		self.curve = curve

		#print "setting localNode to GlobalCurveFit:", self.currNode
		self.globalCurveFit = GlobalCurveFit(self.robotParam, self.contacts, self.curve, localNode = self.currNode)
		#self.direction = not self.getPathDirection()
		
		self.direction = not self.globalCurveFit.getPathDirection()
	
		#print "globalCurveFit direction =", self.direction
	
	def getPathDirection(self):
		#print "globalCurveSlide returning", self.direction
		return self.direction

	def clearDraw(self):
		self.globalCurveFit.clearDraw()

	def draw(self):
		self.globalCurveFit.draw()

	def reset(self, probeState, joints = []):

		self.probeState = probeState
		
		cmdJoints = self.probeState['cmdJoints']
		#print "globalCurveSlide reset:", joints
		
		self.targetState = []
		self.positions = []
		
		#print "resetting HoldSlideTransition to:", joints
		
		if len(joints) > 0:
			for i in range(len(joints)):				
				self.targetState.append(joints[i])
				if joints[i] != None:
					self.positions.append(180.0 / pi * cmdJoints[i])
				else:
					self.positions.append(None)
					
		else:
			for i in range(self.numJoints):
				self.targetState.append(180.0 / pi * cmdJoints[i])
				self.positions.append(180.0 / pi * cmdJoints[i])
	
		self.isDone = False
				
	def step(self, probeState):
		Behavior.step(self, probeState)

		self.resetJoints()


		if self.moveCount == 0:
			if self.direction:
				startNode = 0
				endNode = self.count
				self.globalCurveFit.setBoundaries(startNode, endNode)
			else:
				endNode = self.numJoints-1	
				startNode = endNode - self.count
				self.globalCurveFit.setBoundaries(startNode, endNode)
		
		self.globalCurveFit.step(self.probeState)
		targetState = self.globalCurveFit.getJoints()
		
		#print "globalCurveFit:", targetState
		
		self.positions = copy(targetState)

		self.resetJoints()
		self.mergeJoints([targetState])




		self.moveCount += 1
		#print "moveCount =", self.moveCount
		
		#self.draw()
		
		if self.moveCount > 5 :

			self.moveCount = 0
			#return True
	
			self.count += 1
	
			#print "count =", self.count
	
			" delay period to allow body to stabilize "
			if self.count > 24:
				self.count = 4
				
				#print "GlobalCurveSlide = DONE"
				return True
			
		
		return False
	
	