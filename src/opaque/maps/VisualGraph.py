import PoseGraph
from LocalNode import LocalNode
from Pose import Pose
import pylab
import gen_icp
from numpy import matrix
from random import random

class VisualGraph:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts		
		self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)	

		self.edgeAgeHash = {}

		self.localNodes = []

		self.sensorHypotheses = []
		
		self.colors = []
		for i in range(100):
			self.colors.append((random(),random(),random()))

	def loadSeries(self, dirName, num_poses):
		
		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", self.numNodes			
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)

			if i > 0 and i % 2 == 0:
				
				" since estimated poses are modified from the original motion estimation, we need to restore them "

				f = open(dirName + "/motion_constraints_%04u.txt" % i, 'r')
				motion_constraints = eval(f.read().rstrip())
				f.close()				
				transform = motion_constraints[-1][2]
				offset = [transform[0,0], transform[1,0], transform[2,0]]
				
				estPose1 = self.poseGraph.nodeHash[i-1].getEstPose()
				profile1 = Pose(estPose1)
				estPose2 = profile1.convertLocalOffsetToGlobal(offset)
				
				currNode.readFromFile(dirName, i, forcedPose = estPose2)
			else:
				currNode.readFromFile(dirName, i)

			self.localNodes.append(currNode)
			
			self.poseGraph.loadNewNode(currNode)
	
			#if i % 10 == 0:
			#	self.poseGraph.makeCornerBinConsistent()
	
			self.drawConstraints(i)

		" merge the constraint "
		for k, v in self.poseGraph.edgePriorityHash.items():
			
			if len(v) > 0:
		
				id1 = k[0]
				id2 = k[1]
				
				" take the highest priority constraints only "
				maxPriority = -1
				for const in v:
					thisPriority = const[2]
					if thisPriority > maxPriority:
						maxPriority = thisPriority
		
				if maxPriority == -1:
					raise
				
				for const in v:
					if const[2] == maxPriority:
						transform = const[0]
						covE = const[1]
						#print id1, id2, maxPriority
 
 	def testHypotheses(self, dirName, num_poses, hypFile):
		self.sensorHypotheses = self.poseGraph.sensorHypotheses
		
		self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)
		self.loadFile(dirName, num_poses-1)
		self.poseGraph.sensorHypotheses = self.sensorHypotheses
		
		fn = open(hypFile,'r')
		allCornerHypotheses = eval(fn.read())
		self.poseGraph.allCornerHypotheses = allCornerHypotheses
		fn.close() 
		
		self.poseGraph.processCornerHypotheses()
		self.drawConstraints(1)
		
	def instantSensorTest(self, dirName, num_poses):
		
		self.sensorHypotheses = self.poseGraph.sensorHypotheses
		
		self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)
		self.loadFile(dirName, num_poses-1)
		self.poseGraph.sensorHypotheses = self.sensorHypotheses
		
		self.printCorners()
		
		self.drawConstraints(0)
		self.poseGraph.makeHypotheses()
		self.drawConstraints(1)
		
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
	
		cornerTrans = []

		pylab.clf()
		for i in range(self.numNodes):

			hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False)
			hull.append(hull[0])

			node1 = self.nodeHash[i]
			currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))	

			hull_trans = []
			for p in hull:
				hull_trans.append(gen_icp.dispOffset(p, currPose))	
			
			" extract corner points for each node "
			cornerCandidates = node1.cornerCandidates
			cornerPoints = []
			for cand in cornerCandidates:
				cornerPoints.append(cand[0])
			
			for p in cornerPoints:
				cornerTrans.append(gen_icp.dispOffset(p, currPose))	
				
			
			
			#xP = []
			#yP = []
			#for p in posture_trans:
			#	xP.append(p[0])
			#	yP.append(p[1])
			#pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	



		#for k, v in self.poseGraph.edgeHash.items():	
		#	node1 = self.nodeHash[k[0]]
		#	node2 = self.nodeHash[k[1]]
		#	pose1 = node1.getGlobalGPACPose()
		#	pose2 = node2.getGlobalGPACPose()

		#	xP = [pose1[0], pose2[0]]
		#	yP = [pose1[1], pose2[1]]
			
		#	try:
		#		age = self.edgeAgeHash[k]
		#	except:
		#		age = 0

		#	fval = float(age) / 20.0
		#	if fval > 0.4:
		#		fval = 0.4
		#	color = (fval, fval, fval)
		#	pylab.plot(xP, yP, color=color, linewidth=1)


		#xP = []
		#yP = []
		#for pose in poses:
		#	xP.append(pose[0])
		#	yP.append(pose[1])		
		#pylab.scatter(xP,yP, color='k', linewidth=1)

		#colors = ['b','g', 'r', 'c', 'm', 'y', 'k']

		for i in range(len(self.poseGraph.cornerBins)):
			bin = self.poseGraph.cornerBins[i]
			clusterTrans = []
			for item in bin:
				nodeID, cornerID = item
				currPose = self.nodeHash[nodeID].getGlobalGPACPose()
				cornerP = self.nodeHash[nodeID].cornerCandidates[cornerID][0]
				clusterTrans.append(gen_icp.dispOffset(cornerP, currPose))

			if len(clusterTrans) > 0:
				xP = []
				yP = []
				for p in clusterTrans:
					xP.append(p[0])
					yP.append(p[1])
				pylab.scatter(xP,yP, color=self.colors[i], linewidth=1, zorder=10)	

		#if len(cornerTrans) > 0:
		#	xP = []
		#	yP = []
		#	for p in cornerTrans:
		#		xP.append(p[0])
		#		yP.append(p[1])
		#	pylab.scatter(xP,yP, color='k', linewidth=1, zorder=10)	




			
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
		
	def printCorners(self):
		print "CORNERS"
		for i in range(len(self.localNodes)):
			cornerCandidates = self.localNodes[i].cornerCandidates
			for j in range(len(cornerCandidates)):
				cand = cornerCandidates[j]
				pnt2, cornerAngle, ori = cand
				print i, j, pnt2, cornerAngle, ori
	

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

