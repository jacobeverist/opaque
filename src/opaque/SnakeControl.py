from ControlError import ControlError

class SnakeControl():

	def __init__(self):
		self.timeInc = 1
		self.globalTimer = 0

	def testSuccess(self):
		raise ControlError("test success")
	
	def testFail(self):
		raise ControlError("test fail")

	def setProbe(self, probe):
		self.probe = probe
		
	def setCamera(self, camera):
		self.camera = camera
	
	def setRenderWindow(self, renderWindow):
		self.renderWindow = renderWindow
	
	def setTimerAliasing(self, timeInc):
		self.timeInc = timeInc
		
	def mergeJoints(self, jointList):
		
		result = [None for i in range(self.probe.numSegs-1)]
		
		# merge the joints prioritizing the first to last
		for joints in jointList:
			for i in range(len(joints)):
				if result[i] == None and joints[i] != None:
					result[i] = joints[i]
					
		for i in range(len(result)):
			if result[i] != None:
				self.probe.setServo(i, result[i])
				
	# run control program
	def frameStarted(self):
		self.globalTimer += 1
		
	def isStep(self):
		if self.globalTimer % self.timeInc == 0:
			return True
		
		return False