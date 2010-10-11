from Behavior import Behavior
from Curl import Curl
from Anchor import Anchor

class PokeWalls(Behavior):

	def __init__(self, robotParam, contacts, direction = False, obstContact = 0):
		Behavior.__init__(self, robotParam)

		self.direction = direction
		self.obstContact = obstContact


		self.contacts = contacts
		
		#self.MIN = 6
		self.MAX = 10
		self.MIN = 6
		#self.MAX = 16

		self.topJoint = self.MIN

		self.anchor = Anchor(robotParam)

		self.curl = Curl(robotParam, self.contacts, self.direction)
		
		self.updateCurl()
		#self.curl.setTimerAliasing(10)
		self.curl.setTopJoint(self.topJoint)

		self.doneCount = 0

		self.direction = direction

		self.reverse = False

		self.state = 0

		self.isTransitioning = False
		self.hasTransitioned = False

		self.isInit = False

	def hasInitialized(self):
		return self.isInit
	
	def updateCurl(self):
		#print "setting direction", val
		#self.direction = val

		if self.direction:
			self.curl.setTopJoint(self.topJoint)
		else:
			self.curl.setTopJoint(-self.topJoint + 1 + self.numJoints-1)

		self.curl.setDirection(self.direction)

	def step(self, probeState):
		Behavior.step(self, probeState)
		
		self.updateCurl()
		
		if not self.isStep():
			return False

		self.resetJoints()

		finalDone = False

		self.anchor.step(probeState)
		
		isDone = False

		if self.state == 0:
	
			isDone = self.curl.step(probeState)
				
			if isDone:
				self.doneCount += 1
				self.state = 1
				self.isInit = True
	
		elif self.state == 1:
			
			if self.doneCount >= 2:
				if not self.reverse == True:
					self.topJoint += 2
					if self.direction:
						self.curl.setTopJoint(self.topJoint)
						self.curl.step(probeState)
					else:
						self.curl.setTopJoint(-self.topJoint + 1 + self.numJoints-1)
						self.curl.step(probeState)
					self.doneCount = 0
				else:
					self.topJoint -= 2
					if self.direction:
						self.curl.setTopJoint(self.topJoint)
						self.curl.step(probeState)
					else:
						self.curl.setTopJoint(-self.topJoint + 1 + self.numJoints-1)
						self.curl.step(probeState)
					self.doneCount = 0
	
			self.state = 0
	
			if self.curl.getSlipError() < 0.07:
				" temporarily disabling the obstacle detection callback for speed "
				self.obstContact(self.direction)
			else:
				pass
	
			if self.topJoint > self.MAX:
				self.reverse = True
	
			if self.topJoint < self.MIN and self.reverse:
				finalDone = True
				self.reverse = False
				self.topJoint = self.MIN
	

		joints2 = self.curl.getJoints()

		# joint value
		self.mergeJoints([joints2])

		if finalDone:
			self.hasTransitioned = False
			self.curl.reset()
			self.isInit = False

		return finalDone

