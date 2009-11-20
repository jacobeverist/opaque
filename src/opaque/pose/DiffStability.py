from math import *

# number of servo position samples for stability determinatin
SAMPLES = 20

class DiffStability:

	def __init__(self, thresh, sampleCount):
		self.thresh = thresh
		self.sampleCount = sampleCount

		self.currAngle = 0.0
		self.nomAngle = 0.0
		self.isInit = False
		self.count = 0
		
		"""
		1. set nominal angle
		2. for each iteration, check if current angle violates difference threshold
		3. otherwise set as stable
		4. reset
		"""	
		
	def getDiff(self):
		return self.currAngle - self.nomAngle
		
	def setThresh(self, val):
		self.thresh = val
		
	def isStable(self):
		
		if not self.isInit:
			return False
		
		if fabs(self.currAngle - self.nomAngle) > self.thresh:
			return False
		
		return True

	def initialize(self, newValue):
		self.currAngle = newValue
		self.nomAngle = newValue
		self.isInit = True
		self.count = self.sampleCount+1

	
	def reset(self):
		self.currAngle = 0.0
		self.nomAngle = 0.0
		self.isInit = False
		self.count = 0

	def addData(self, newValue):
		
		self.count += 1
		
		if self.count < self.sampleCount:
			return False
		
		if not self.isInit:
			self.nomAngle = newValue
			self.isInit = True
					
		self.currAngle = newValue	
		return self.isStable()
	
