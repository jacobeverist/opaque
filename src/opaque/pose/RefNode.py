
class RefNode:

	def __init__(self, nid, jid, x, z, p, probe):

		#self.statNode = probe.statNode

		#self.currEntity = probe._mgr.createEntity("ent"+str(nid), "Cube.mesh")
		#self.currEntity.setCastShadows(False)
		#self.currEntity.setMaterialName("Green")

		#position = ogre.Vector3(x,1.0,z)
		#size = ogre.Vector3(0.05,0.05,0.05)
		#self.statNode.addEntity(self.currEntity, position, scale = size)
		#self.statNode.build()

		self.poseError = 0.

		self.edges = []

		self.refX = x
		self.refZ = z
		self.refP = p

		self.gndX = x
		self.gndZ = z
		self.gndP = p

		self.nodeID = nid
		self.jointID = jid

		self.probe = probe

		self.numJoints = self.probe.getNumJoints()

		self.sweep = []

		self.nom_ang = []
		for i in range(self.numJoints):
			self.nom_ang.append(probe.getServo(i))

		self.maxError = 0.
		self.active_time = 0
		self.errSum = 0.
		self.currError = 0.
		self.maxErrorReached = False
		self.newRoot = False


	def getRefPose(self):
		return self.refX, self.refZ, self.refP

	def getNodeID(self):
		return self.nodeID

	def getJointID(self):
		return self.jointID

	def setGroundTruthPose(self, x, z, p):
		self.gndX = x
		self.gndZ = z
		self.gndP = p

	def getGroundTruthPose(self):
		return self.gndX, self.gndZ, self.gndP

	def resetPoseError(self):
		self.poseError = 0.0

	def getPoseError(self):
		return self.poseError

	def addPoseError(self, error):
		self.poseError += error

	def getNumEdges(self):
		return len(self.edges)

	def addEdge(self, newEdge):

		if newEdge.origin != self and newEdge.target != self:
			print "invalid edge attempted to add to node"
			raise

		for i in range(len(self.edges)):
			if self.edges[i] == newEdge:
				print "tried to add duplicate edge!"
				return 

		self.edges.append(newEdge)


	def getEdge(self, i):
		return self.edges[i]

	def computeStabilityError(self):
	
		# saddle points: 4,8,12,16,20,24,28,32,36
		self.interval = self.jointID % 4

		# find the saddle point joints to which to check stability
		if self.interval == 0:
			self.low_index = self.jointID
			self.high_index = self.jointID + 4

		else:

			self.low_index = self.jointID - self.interval
			self.high_index = self.low_index + 4
		
		# set the max and min limits if they go over
		if self.low_index < 0:
			self.low_index = 0

		if self.high_index >= self.numJoints:
			self.high_index = self.numJoints - 1

		error = 0.0
		for i in range(self.low_index,self.high_index + 1):
			error += abs(self.nom_ang[i] - self.probe.getServo(i))

		if error > self.maxError:
			self.maxError = error

		self.currError = error
		self.errSum += error

		if self.currError > 0.1:
			self.maxErrorReached = True

		return error

	def getStabilityError(self):
		return self.currError

	def getMaxStabilityError(self):
		return self.maxError

	def getAvgStabilityError(self):
		return self.errSum/self.active_time

	def isMaxErrorReached(self):
		return self.maxErrorReached

	def updateTime(self):
		self.active_time += 1

	def setNewRoot(self, isRoot):
		self.newRoot = isRoot

	def isNewRoot(self):
		return self.newRoot

	def addAngleData(self, angleSet):
		self.sweep.append(angleSet)

	def getAngleData(self):
		return self.sweep

	def getNumJoints(self):
		return self.numJoints

