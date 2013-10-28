from MapState import MapState

class BayesMapper:

	def __init__(self, probe, contacts):
		
		self.probe = probe
		self.contacts = contacts
		
		self.walls = self.probe.getWalls()
		
		self.nodeHash = {}
		self.numNodes = 0


		""" initialize to a single map hypothesis """
		self.particleIDs = 0
		self.mapHyps = {}
		self.mapHyps[self.particleIDs] = MapState(probe,contacts)
		self.particleIDs += 1


	@logFunction
	def loadNewNode(self, newNode):

		self.currNode = newNode
		nodeID = self.numNodes
		self.nodeHash[nodeID] = self.currNode
		self.numNodes += 1
		
		""" for each current map hypothesis, integrate the new node """
		currHyps = self.mapHyps
		for mid, mapHyp in currHyps.iteritems():
			self.integrateNode(newNode, nodeID, mapHyp)


	@logFunction
	def integrateNode(self, newNode, nodeID, mapHyp):

		" DIRECTION OF TRAVEL FROM PREVIOUS POSE PAIR "
		direction = newNode.travelDir

		" ensure the medial axes are computed before this check "
		computeHullAxis(nodeID, newNode, tailCutOff = False)


		if nodeID > 0:
			
			" ESTIMATE TRAVEL WITH MEDIAL OVERLAP CONSTRAINT OF EVEN NUMBER POSE "
			if self.numNodes >= 4:

				" Move node along path "
				self.movePath(nodeID, direction)
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
		
		
		" CHECK FOR A BRANCHING EVENT "
		
		if self.numNodes >= 2 and self.numNodes % 2 == 0:

			" DETECT BRANCHING EVENTS FOR THE 2 NODES OF LAST STOP "
			" AFTER ODD-NUMBER OVERLAP OF CURRENT STOP HAS IRONED OUT ERRORS "
			nodeID1 = self.numNodes-2
			nodeID2 = self.numNodes-1
			
			" if these nodes are already path-classified, return"
			isContained1 = False
			isContained2 = False
			
			pathIDs = self.paths.getPathIDs()
			for k in pathIDs:
				if self.paths.getNodes(k).count(nodeID1) > 0:
					isContained1 = True
				if self.paths.getNodes(k).count(nodeID2) > 0:
					isContained2 = True
					
			
			if isContained1 or isContained2:
				return

			" COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
			self.paths.generatePaths()

				
			self.drawConstraints(self.statePlotCount)
			self.statePlotCount += 1
			self.drawPathAndHull()


			self.addToPaths(nodeID1, nodeID2)

			" IF THIS IS THE FIRST NODES IN THE FIRST PATH, JUST ADD THEM AS DEFAULT, DO NOTHING "
			if len(self.paths.paths[0]) > 0:
				
				self.mergePaths()

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				paths = {}
				pathIDs = self.paths.getPathIDs()



				orderedPathIDs1 = self.paths.getOrderedOverlappingPaths(nodeID1)
				print nodeID1, "recompute orderedPathIDs1:", orderedPathIDs1

			
				orderedPathIDs2 = self.paths.getOrderedOverlappingPaths(nodeID2)
				print nodeID2, "recompute orderedPathIDs2:", orderedPathIDs2


				for k in pathIDs:
					path = self.paths.paths[k]
					if len(path) > 0:
						paths[k] = path 
				self.trimmedPaths = self.paths.trimPaths(paths)

				
				self.drawTrimmedPaths(self.trimmedPaths)
				print "trimmed paths:", len(self.trimmedPaths)



				"1)  guess pose in front if there is a departure with local and path junction points (not first node) "
				"2)  guess pose in back if there is a departure with local and path junction points (not first node) "
				"3)  select choice with the lowest cost "
				"4)  guess pose with no departure point with guessed control points len(orderedPathIDs) == 1"
				
				" nodeID1:	is it in a junction or a single path? "
				
				
				hull1, medial1 = computeHullAxis(nodeID1, self.nodeHash[nodeID1], tailCutOff = False)
				hull2, medial2 = computeHullAxis(nodeID2, self.nodeHash[nodeID2], tailCutOff = False)
				estPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()
				estPose2 = self.nodeHash[nodeID2].getGlobalGPACPose()
				self.currSplicePath = self.selectSplice(nodeID1, nodeID2, medial1, medial2, estPose1, estPose2, orderedPathIDs1, orderedPathIDs2)
				
				self.paths.generatePaths()
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				self.drawTrimmedPaths(self.trimmedPaths)

				try:
					self.consistentFit(nodeID1,  self.nodeHash[nodeID1].getGlobalGPACPose(), numGuesses = 11)
				except IndexError:
					print "failed to consistentFit node", nodeID1
					pass
					
				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()

				try:
					self.consistentFit(nodeID2,  self.nodeHash[nodeID2].getGlobalGPACPose(), numGuesses = 11)
				except IndexError:
					print "failed to consistentFit node", nodeID2
					pass

				self.drawConstraints(self.statePlotCount)
				self.statePlotCount += 1
				self.drawPathAndHull()
				
			self.paths.generatePaths()
			self.drawPathAndHull()


			if nodeID1 >= 2:
				newPath = deepcopy(self.paths.paths[0])
				p0 = newPath[0]
				pN = newPath[-1]
				
				rootPose = [-2.0,0.0]
				
				dist0 = sqrt((rootPose[0]-p0[0])*(rootPose[0]-p0[0]) + (rootPose[1]-p0[1])*(rootPose[1]-p0[1]))
				distN = sqrt((rootPose[0]-pN[0])*(rootPose[0]-pN[0]) + (rootPose[1]-pN[1])*(rootPose[1]-pN[1]))
		
				if dist0 > distN:
					newPath.reverse()
							
				" new path "
				newSpline = SplineFit(newPath, smooth = 0.1)
	
	
				posNew = self.nodeHash[nodeID1].getGlobalGPACPose()
				posOld = self.nodeHash[nodeID1-2].getGlobalGPACPose()
				minDist, minU, closestPoint = newSpline.findClosestPoint(posNew)
				arcDistNew = newSpline.dist(0.0, minU)
	
				minDist, minU, closestPoint = newSpline.findClosestPoint(posOld)
				arcDistOld = newSpline.dist(0.0, minU)
	
				print "arcDistNew, arcDistOld, diff =", arcDistNew, arcDistOld, arcDistNew-arcDistOld
		
				
				

		for k in range(self.numNodes):
			self.nodePoses[k] = self.nodeHash[k].getGlobalGPACPose()




	@logFunction
	def saveState(self):
		pass

	def insertPose(self):
		pass

	def restoreState(self):
		pass



