
from math import sqrt, pi, fabs, cos, sin
from Behavior import Behavior
from ValueStability import ValueStability
from Transition import Transition

class Curl(Behavior):
	
	def __init__(self, robotParam, contacts, direction):
		Behavior.__init__(self, robotParam)	

		# angle value for curling joints
		self.val = 0
		self.rail = 25.0
		self.up = True
		self.topJoint = 0

		self.poseWindow = []
		
		self.contacts = contacts

		#self.stabilityX = ValueStability(1e-7, 10)
		#self.stabilityY = ValueStability(1e-7, 10)
		#self.stabilityX = ValueStability(1e-5, 10)
		#self.stabilityY = ValueStability(1e-5, 10)
		self.stabilityX = ValueStability(1e-4, 10)
		self.stabilityY = ValueStability(1e-4, 10)
		
		self.lastPose = [0.0,0.0,0.0]
		
		self.poseIndex = 2
		
		self.returnedSuccess = False
		self.slipError = 0.0
		
		self.direction = direction

		self.transition = Transition(robotParam)
		self.isTransitioning = False
		self.hasTransitioned = False
	
	def reset(self):
		
		self.val = 0
		self.rail = 25.0
		self.up = True
		self.topJoint = 0

		self.poseWindow = []

		self.stabilityX.reset()
		self.stabilityY.reset()
		
		self.lastPose = [0.0,0.0,0.0]
		
		self.poseIndex = 2
		
		self.returnedSuccess = False
		self.slipError = 0.0
		
		self.direction = True

		self.transition = Transition(self.robotParam)
		self.isTransitioning = False
		self.hasTransitioned = False

		self.timeInc = 1
		self.gaitTimer = 0
		
		self.resetJoints()
			
	def setDirection(self, val):
		self.direction = val

		self.resetJoints()
		
		if self.direction:
			for i in range(0, self.topJoint):
				# joint value
				self.setJoint(i, self.val)
		else:
			for i in range(self.topJoint, self.numJoints):
				# joint value
				self.setJoint(i, self.val)		
				
	def getSlipError(self):
		return self.slipError
	
	def setTopJoint(self, node):
		self.topJoint = node

	def normalizeAngle(self,angle):
		# this function converts the angle to its equivalent
		# in the range [-pi,pi] 

		while angle>pi:
			angle=angle-2*pi

		while angle<=-pi:
			angle=angle+2*pi

		return angle

	def step(self, probeState):
		Behavior.step(self, probeState)

		if not self.isStep():
			return False
		
		self.resetJoints()

		if self.isTransitioning:

			# first steps
			isDone = self.transition.step(probeState)

			#print "transitioning"
			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])

			if isDone:
				self.isTransitioning = False

			return False

		if not self.isTransitioning and self.hasTransitioned:

			stateJoints = probeState['joints']
			segLength = self.robotParam['segLength']

			if self.direction:
				#pose = self.contacts.getAveragePose(-1)
				pose = self.contacts.getAveragePose(0)

				pose[2] = pose[2] + stateJoints[0]
				pose[0] = pose[0] - segLength*cos(pose[2])
				pose[1] = pose[1] - segLength*sin(pose[2])

				pose[2] = self.normalizeAngle(pose[2])

			else:
				#pose = self.contacts.getAveragePose(39)
				pose = self.contacts.getAveragePose(38)

				pose[0] = pose[0] + segLength*cos(pose[2])
				pose[1] = pose[1] + segLength*sin(pose[2])

				
			self.stabilityX.addData(pose[0])
			self.stabilityY.addData(pose[1])
			
			varX = self.stabilityX.getVar()
			varY = self.stabilityY.getVar()
					
			print "variances =", varX, varY
						
			if self.stabilityX.isStable() and self.stabilityY.isStable():
				
				diff = [pose[0]-self.lastPose[0], pose[1]-self.lastPose[1]]
				
				if self.val == 30 or self.val == -30:
					if not self.returnedSuccess:
						
						self.slipError = sqrt(diff[0]**2 + diff[1]**2)
						
						self.returnedSuccess = True
						return True
					else:
						self.returnedSuccess = False
						
				self.stabilityX.reset()
				self.stabilityY.reset()
				self.lastPose = pose
				self.resetJoints()
				
				if self.val == 0:
					#self.stabilityX.setThresh(1e-7)
					#self.stabilityY.setThresh(1e-7)
					self.stabilityX.setThresh(1e-4)
					self.stabilityY.setThresh(1e-4)
				else:
					self.stabilityX.setThresh(1e-4)
					self.stabilityY.setThresh(1e-4)
				
			else:
				return False
		
			# 1. increment only after the end point becomes stable
			# 2. measure distance between new end point and old end point, and this determines our contact
			
			seq = [-30,-25,0,25,30]
			
			if self.up:
				if self.poseIndex == 0:
					self.poseIndex += 2
				else:
					self.poseIndex += 1
					
				self.val = seq[self.poseIndex]
				#self.val += 25
			else:
				if self.poseIndex == 4:
					self.poseIndex -= 2
				else:
					self.poseIndex -= 1
				self.val = seq[self.poseIndex]
				#self.val += -25
				
			if self.poseIndex >= 4:
				self.up = False
				
			if self.poseIndex <= 0:
				self.up = True
			
			if self.direction:
				for i in range(0, self.topJoint):
					# joint value
					self.setJoint(i, self.val)
			else:
				for i in range(self.topJoint, self.numJoints):
					# joint value
					self.setJoint(i, self.val)

			self.hasTransitioned = False

		if not self.hasTransitioned:
	
			cmdJoints = self.getJoints()
			stateJoints = probeState['joints']
			
			initState = []
			for i in range(self.numJoints):
				if cmdJoints[i] != None:
					initState.append(180.0*stateJoints[i]/pi)
				else:
					initState.append(None)

			targetState = self.getJoints()

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])

			#transTime = int(4*errSum)
			transTime = int(errSum)
			
			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			self.transition.resetTime(transTime)		

			# first steps
			self.transition.step(probeState)

			resJoints = self.transition.getJoints()
			self.resetJoints()
			self.mergeJoints([resJoints])

			self.hasTransitioned = True
			self.isTransitioning = True				

		return False
	
	
	
	