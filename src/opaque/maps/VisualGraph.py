import PoseGraph
from LocalNode import LocalNode
from Pose import Pose
import pylab
import gen_icp

class VisualGraph:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts		
		self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)	

		self.edgeAgeHash = {}

		self.localNodes = []

		self.sensorHypotheses = []

	def sensorTest(self, dirName, num_poses):
		
		for i in range(1, num_poses):
			
			self.sensorHypotheses = self.poseGraph.sensorHypotheses
			
			self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)
			self.loadFile(dirName, i)
			self.poseGraph.sensorHypotheses = self.sensorHypotheses

			
			self.poseGraph.makeHypotheses()
			self.drawConstraints(i)

	def visualize(self, dirName, num_poses):
		
		for i in range(1, num_poses):
			self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)
			self.loadFile(dirName, i)
			
			print "AGE:"
			for k, v in self.poseGraph.edgeHash.items():
				
				tup1 = k
				tup2 = (k[1], k[0])
				
				try:
					self.edgeAgeHash[tup1] += 1
					self.edgeAgeHash[tup2] += 1
				except:
					self.edgeAgeHash[tup1] = 0
					self.edgeAgeHash[tup2] = 0
			
				print tup1, "=", self.edgeAgeHash[tup1]
			
			self.drawConstraints(i)

	def drawConstraints(self, id = []):
		
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())

		pylab.clf()
		for i in range(self.numNodes):

			hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False)
			hull.append(hull[0])

			node1 = self.nodeHash[i]
			currPose = currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	

			hull_trans = []
			for p in hull:
				hull_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	

		for k, v in self.poseGraph.edgeHash.items():	
			node1 = self.nodeHash[k[0]]
			node2 = self.nodeHash[k[1]]
			pose1 = node1.getGlobalGPACPose()
			pose2 = node2.getGlobalGPACPose()

			xP = [pose1[0], pose2[0]]
			yP = [pose1[1], pose2[1]]
			
			try:
				age = self.edgeAgeHash[k]
			except:
				age = 0

			fval = float(age) / 20.0
			if fval > 0.4:
				fval = 0.4
			color = (fval, fval, fval)
			pylab.plot(xP, yP, color=color, linewidth=1)


		xP = []
		yP = []
		for pose in poses:
			xP.append(pose[0])
			yP.append(pose[1])		
		pylab.scatter(xP,yP, color='k', linewidth=1)


			
		pylab.xlim(-5,10)
		#pylab.xlim(-8,12)
		#pylab.ylim(-10,10)
		pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % i)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)

	def localizeCurrentNode(self):
		self.poseGraph.localizeCurrentNode()

	def __getattr__(self, name):

		if name == "numNodes":
			return self.poseGraph.numNodes

		if name == "nodeHash":
			return self.poseGraph.nodeHash

		if name == "edgeHash":
			return self.poseGraph.edgeHash

		if name == "overlap_constraints":
			return self.poseGraph.overlap_constraints
		if name == "motion_constraints":
			return self.poseGraph.motion_constraints
		if name == "sensor_constraints":
			return self.poseGraph.sensor_constraints
		if name == "inplace_constraints":
			return self.poseGraph.inplace_constraints
		if name == "gnd_constraints":
			return self.poseGraph.gnd_constraints
		if name == "merged_constraints":
			return self.poseGraph.merged_constraints

	def loadFile(self, dirName, num_poses):
		
		numNodes = 0	
		PIXELSIZE = 0.05
		nodeHash = {}
		
		for i in range(0,num_poses):
			print "loading node", i			
			
			" don't reload the node if we've already done it before "
			if i < len(self.localNodes):
				currNode = self.localNodes[i]
			else:	
				currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
				currNode.readFromFile(dirName, i)
				self.localNodes.append(currNode)

			nodeHash[i] = currNode
			numNodes += 1

		" add the nodes to the graph "
		self.poseGraph.loadNodes(nodeHash)
		
		" add the constraint edges to the graph "
		self.poseGraph.loadConstraints(dirName, num_poses-1)

		#self.drawConstraints()

		" save the maps "		
		#self.synch()
		#self.saveMap()

	def addFile(self, dirName, num_poses):
		
		numNodes = 0	
		PIXELSIZE = 0.05
		nodeHash = {}
		
		for i in range(0,num_poses):
			print "loading node", i			
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile(dirName, i)
			nodeHash[i] = currNode

			numNodes += 1

		" add the nodes to the graph "
		self.poseGraph.loadNodes(nodeHash)
		
		" add the constraint edges to the graph "
		self.poseGraph.loadConstraints(dirName, num_poses-1)

		#self.drawConstraints()

		" save the maps "		
		#self.synch()
		#self.saveMap()


	def newInPlaceNode(self, direction):
		
		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, direction = direction)
		
		self.poseGraph.addInPlaceNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
	
	def newNode(self, stepDist, direction):

		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, direction = True)		
		
		self.poseGraph.addNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

	def initNode(self, estPose, direction):

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, direction = True)

		self.poseGraph.addInitNode(self.currNode, estPose)
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)	

	def correctPosture(self):
		self.currNode.updateCorrectPosture()

	def saveLocalMap(self):
		" save the local map "
		if self.currNode != 0:
			self.currNode.saveMap()
			
	def saveConstraints(self):
		self.poseGraph.saveConstraints()

