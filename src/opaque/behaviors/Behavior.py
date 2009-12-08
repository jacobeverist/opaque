
class Behavior:
	
	def __init__(self, probe, contacts = 0):
		self.probe = probe
		self.contacts = contacts
		self.timeInc = 1
		self.gaitTimer = 0
		
		self._joints = [None for i in range(self.probe.numSegs-1)]
	
	def getJoints(self):
		return self._joints
	
	def resetJoints(self):
		self._joints = [None for i in range(self.probe.numSegs-1)]
	
	def setJoint(self, i, val):
		self._joints[i] = val
		
	def mergeJoints(self, jointList):
		
		result = [None for i in range(self.probe.numSegs-1)]
		
		# merge the joints prioritizing the first to last
		for joints in jointList:
			for i in range(len(joints)):
				if result[i] == None and joints[i] != None:
					result[i] = joints[i]
					
		for i in range(len(result)):
			if result[i] != None:
				self.setJoint(i, result[i])
		
	def setTimerAliasing(self, timeInc):
		self.timeInc = timeInc

	def step(self):
		self.gaitTimer += 1
	
	def isStep(self):
		if self.gaitTimer % self.timeInc == 0:
			return True
		
		return False

	def resetTimer(self):
		self.gaitTimer = 0