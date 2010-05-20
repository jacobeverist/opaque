
from functions import *
from DiffStability import DiffStability

"""
1. reset the root node to (0,0,0)
2. flag the nodes that are currently anchored
3. detect discontinuities and choose largest stable section to be assumed fixed wrt environment
4. recompute the root node position
5. cut out data from unstable sections of time
"""


class StablePose:
	
	def __init__(self, probe):
		self.valStabs = [DiffStability(0.05, 1) for i in range(39)]
		self.probe = probe
		self.node = 0
		self.direction = True
		self.refPoses = [[0.0,0.0,0.0] for i in range(39)]
		self.isStables = [False for i in range(39)]
		self.maxWidth = 0
	
		#self.sectFile = open("sectFile.txt", 'w')
		self.rootFile = open("rootFile.txt", 'w')

		self.lowIndex = 0
		self.highIndex = 0
		self.prevStables = [False for i in range(39)]
		self.count = 0
		
	def isStable(self):
		
		self.count += 1
		
		#print self.maxWidth, self.highIndex, self.lowIndex
		
		if self.maxWidth <= 6:
			"""
			
			if self.count >= 210:
				for i in range(len(self.valStabs)):
					self.valStabs[i].initialize(self.probe.getServo(i))
			"""
			
			#print "maxWidth =", self.maxWidth
			return False

		#if self.direction and self.highIndex > 19:
		if self.direction and self.highIndex >= 38:
			self.count = 200
			return True
		
		#if not self.direction and self.highIndex < 19:
		if not self.direction and self.lowIndex <= 0:
			self.count = 200
			return True

		"""
		if self.count >= 210:
			for i in range(len(self.valStabs)):
				self.valStabs[i].initialize(self.probe.getServo(i))
		"""
		
#		print "default Fail, maxWidth =", self.maxWidth

		return False
	
	def setDirection(self, isForward):
		self.direction = True
	
	def update(self):
		
		for i in range(39):
			val = self.probe.getServo(i)
			self.valStabs[i].addData(val)

			#if self.valStabs[i].sampleCount >= self.valStabs[i].sampleSize:
			if True:
				if self.valStabs[i].isStable():
					self.isStables[i] = True
					self.prevStables[i] = True
				else:
					if self.prevStables[i]:
						self.prevStables[i] = False
						self.valStabs[i].reset()
						
					self.isStables[i] = False
					
		
		" now find the largest number of consecutive stable joints "
		self.lowIndex = 0
		self.highIndex = 0
		maxWidth = 0
		
		currLow = 0
		currHigh = 0
		currWidth = 0
		isContiguous = False
		
		for i in range(39):
			
			if self.isStables[i]:
				if isContiguous:
					currHigh = i
					currWidth = currHigh - currLow
					
					if currWidth > maxWidth:
						maxWidth = currWidth
						self.highIndex = currHigh
						self.lowIndex = currLow
						
				else:
					currLow = i
					currHigh = currLow
					isContiguous = True
					
			else:
				isContiguous = False

		#self.sectFile.write(str(self.lowIndex) + " " + str(self.highIndex) + " " + str(maxWidth) + " " + repr(self.isStables))
		#self.sectFile.write("\n")
		
		#print maxWidth > 0, self.node.rootNode < self.lowIndex , self.node.rootNode > self.highIndex
		
		self.maxWidth = self.highIndex - self.lowIndex
		
		" 1. check if root node is within contiguous stable region or not "
		if False:
		#if maxWidth > 0 and (self.node.rootNode < self.lowIndex or self.node.rootNode > self.highIndex):
		
			" 2. if not, compute the root node's new position with respect to a node in stable region "
			diff = self.highIndex - self.lowIndex
			diff /= 2
			
			refNode = self.lowIndex + diff
			
			"""
			if self.lowIndex < 19 and self.highIndex < 29:
				refNode = 29
			elif self.lowIndex > 9 and self.highIndex > 19:
				refNode = 9
			else:
				refNode = self.lowIndex + diff
			"""
			
			resultPose = self.probe.getJointWRTJointPose(self.refPoses[refNode], refNode, self.node.rootNode)
			
			self.refPoses[self.node.rootNode] = resultPose
			
			self.node.setEstPose(resultPose)
		
		"""
		if self.direction:
			refNode = 29
			resultPose = self.probe.getJointWRTJointPose(self.refPoses[refNode], refNode, self.node.rootNode)
			self.refPoses[self.node.rootNode] = resultPose
			self.node.setEstPose(resultPose)

		else:
			refNode = 9
			resultPose = self.probe.getJointWRTJointPose(self.refPoses[refNode], refNode, self.node.rootNode)
			self.refPoses[self.node.rootNode] = resultPose
			self.node.setEstPose(resultPose)
		"""
		
		" 3. compute all the node positions with respect to the root node "

		" compute reference poses "
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1

		xTotal = self.refPoses[self.node.rootNode][0]
		zTotal = self.refPoses[self.node.rootNode][1]
		totalAngle = self.refPoses[self.node.rootNode][2]

		self.rootFile.write(str(xTotal) + " " + str(zTotal) + " " + str(totalAngle) + "\n")
	
		joints = range(-1,self.node.rootNode)
		joints.reverse()

		for i in joints:

			totalAngle = totalAngle + self.probe.getServo(i+1)
			totalAngle = normalizeAngle(totalAngle)

			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

			self.refPoses[i] = [xTotal,zTotal,totalAngle]

		joints = range(self.node.rootNode+1, self.probe.numSegs-1)

		xTotal = self.refPoses[self.node.rootNode][0]
		zTotal = self.refPoses[self.node.rootNode][1]
		totalAngle = self.refPoses[self.node.rootNode][2]

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			totalAngle = totalAngle - self.probe.getServo(i)
			totalAngle = normalizeAngle(totalAngle)
			
			self.refPoses[i] = [xTotal,zTotal,totalAngle]
		
	
	def reset(self):
		for stab in self.valStabs:
			stab.reset()

		for i in range(39):
			self.isStables[i] = False
		
		self.refPoses = [[0.0,0.0,0.0] for i in range(39)]
		self.isStables = [False for i in range(39)]

		self.lowIndex = 0
		self.highIndex = 0

	def setNode(self, node):
		self.node = node
		
		#self.rootFile.write("\n" + str(node.nodeID) + "\n")
		#self.sectFile.write("\n" + str(node.nodeID) + "\n")
			
		" compute reference poses "
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		self.rootFile.write(str(xTotal) + " " + str(zTotal) + " " + str(totalAngle) + "\n")
	
		joints = range(-1,self.node.rootNode)
		joints.reverse()

		for i in joints:

			totalAngle = totalAngle + self.probe.getServo(i+1)
			totalAngle = normalizeAngle(totalAngle)

			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

			self.refPoses[i] = [xTotal,zTotal,totalAngle]

		joints = range(self.node.rootNode+1, self.probe.numSegs-1)

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			totalAngle = totalAngle - self.probe.getServo(i)
			totalAngle = normalizeAngle(totalAngle)
			
			self.refPoses[i] = [xTotal,zTotal,totalAngle]
					
		
