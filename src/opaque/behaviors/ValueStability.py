
# number of servo position samples for stability determinatin
SAMPLES = 20

class ValueStability:

	def __init__(self, thresh, sample_size = SAMPLES):
		self.sampleMean = 0.0
		self.sampleVar = 0.0
		self.sampleCount = 0
		self.sampleIndex = 0
		self.varThresh = thresh
		
		self.sampleSize = sample_size
		
		self.samples=[]
		for i in range(self.sampleSize):
			self.samples.append(0.0)
		
			
	def setThresh(self, val):
		self.varThresh = val
		
	def isStable(self):
		if self.sampleCount < self.sampleSize:
			return False

		if self.sampleVar < self.varThresh:
			return True
		
		return False

	def getMean(self):
		return self.sampleMean

	def getVar(self):
		return self.sampleVar
	
	def getSampleCount(self):
		return self.sampleCount

	def reset(self):
		self.sampleCount = 0
		self.sampleIndex = 0
		self.sampleMean = 0.0
		self.sampleVar = 0.0

	def addData(self, newValue):

		# **** COMPUTE SAMPLE VARIANCE AND MEAN  ****
		#	
		# NOTE: Samples are taken at every time step and are used to determine
		# the stability of the probe position.  Low variance indicates that the
		# position is stable and it is safe to compute a feature point.
		#

		# retrieve data sample
		self.samples[self.sampleIndex] = newValue

		# increment the circular index
		self.sampleIndex += 1
		self.sampleIndex = self.sampleIndex % self.sampleSize

		if self.sampleCount < self.sampleSize:
			self.sampleCount += 1

		# compute the mean
		self.sampleSum = 0.0
		for i in range(self.sampleCount):
			self.sampleSum += self.samples[i]
		self.sampleMean = self.sampleSum / self.sampleCount 

		# compute the variance
		self.sampleSum = 0
		for i in range(self.sampleCount):
			self.sampleSum += (self.samples[i] - self.sampleMean)*(self.samples[i] - self.sampleMean)

		self.sampleVar = self.sampleSum/self.sampleCount

		max = -1e200
		min = 1e200
		for i in range(0,self.sampleCount):
			if self.samples[i] > max:
				max = self.samples[i]
			if self.samples[i] < min:
				min = self.samples[i]

		# should have at least 10 values before we compute the variance
		if self.sampleCount < self.sampleSize :
			# make extremely high for invalid variances
			self.sampleVar = 1e200;

		if self.sampleVar < self.varThresh:
			return True

		return False


