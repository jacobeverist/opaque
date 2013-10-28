





class MapState:
	
	def __init__(self, probe, contacts):
		
		self.probe = probe
		self.contacts = contacts
		
		#self.walls = self.probe.getWalls()
		
		self.nodeHash = {}
		self.nodePoses = {}
		
		self.pathClasses = {}
		self.pathClasses[0] = {"parentID" : None, "branchNodeID" : None, "localJunctionPose" : None, 
							"sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : None}

	def restoreState(self):
		pass

	@logFunction
	def saveState(self):
		pass



	def addBranch(self):
		pass

	def computeConsistency(self):
		pass

	def updatePose(self, nodeID, pose):
		pass

