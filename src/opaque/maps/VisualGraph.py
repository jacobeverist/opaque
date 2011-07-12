import PoseGraph
from LocalNode import LocalNode
from Pose import Pose
import pylab
import gen_icp
from numpy import matrix
from random import random
import functions
from copy import copy
from OccupancyMap import OccupancyMap


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

		MAPSIZE = 20.0
		self.occMap = OccupancyMap(self.probe, self, MAPSIZE)

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
			#self.poseGraph.loadNodeGroundPast(currNode)
	
			#if i % 10 == 0:
			#	self.poseGraph.makeCornerBinConsistent()
	
			self.drawConstraints(i)
			
			#self.drawShapeConstraints(i)

		self.drawMap()
		self.drawConstraints(num_poses)

		
		self.poseGraph.performShapeConstraints()
		self.drawMap()
		self.drawConstraints(num_poses+1)

		for i in range(0, num_poses):
			self.poseGraph.addCornerConstraints(i)			
			self.poseGraph.mergePriorityConstraints()
		
		self.drawMap()
		self.drawConstraints(num_poses+2)
		
		self.poseGraph.resetGraphToGround()
		self.drawMap()
		self.drawConstraints(num_poses+3)

		#self.drawPoses()

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


	def drawPoses(self):
		
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())
	
		for i in range(self.numNodes):

			pylab.clf()

			hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False)
			hull.append(hull[0])

			node1 = self.nodeHash[i]
			currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture2 = node1.getGPACPosture()
			
			posture1_trans = []
			for p in posture1:
				posture1_trans.append(gen_icp.dispOffset(p, currPose))	
			posture2_trans = []
			for p in posture2:
				posture2_trans.append(gen_icp.dispOffset(p, currPose))	

			hull_trans = []
			for p in hull:
				hull_trans.append(gen_icp.dispOffset(p, currPose))	
			
			xP = []
			yP = []
			#for p in posture_trans:
			for p in posture1:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

			xP = []
			yP = []
			for p in posture2:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(1.0,0.0,0.0))	

			xP = []
			yP = []
			#for p in hull_trans:
			for p in hull:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	

			count1 = 0
			count2 = 0

			for p in posture1:
				if functions.point_inside_polygon(p[0],p[1],hull):
					count1 += 1
	
			for p in posture2:
				if functions.point_inside_polygon(p[0],p[1],hull):
					count2 += 1
					
			#pylab.xlim(-5,10)
			#pylab.ylim(-8,8)
			pylab.title("count1 = %d, count2 = %d" % (count1,count2))
			pylab.xlim(-3,3)
			pylab.ylim(-2,2)
			pylab.savefig("plotPose%04u.png" % i)

	def drawShapeConstraints(self, id):

		numEdges = 0

		for edge in self.poseGraph.badHypotheses2:

			pylab.clf()

			nodeID1 = edge[0]
			nodeID2 = edge[1]
			status = edge[4]

			if edge[0] < edge[1]:
				nodeID1 = edge[0]
				nodeID2 = edge[1]
				transform = edge[2]
				offset = [transform[0,0], transform[1,0], transform[2,0]]

			else:
				nodeID1 = edge[1]
				nodeID2 = edge[0]
				transform = edge[2]
				offset = [transform[0,0], transform[1,0], transform[2,0]]
				prof1 = Pose(offset)
				offset = prof1.convertGlobalPoseToLocal([0.0,0.0,0.0])
		

			hull1 = PoseGraph.computeBareHull(self.nodeHash[nodeID1], sweep = False)
			hull1.append(hull1[0])

			node1 = self.nodeHash[nodeID1]
			currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture1_trans = []
			for p in posture1:
				posture1_trans.append(gen_icp.dispOffset(p, currPose))	

			occ1 = node1.getOccPoints()
			occ1_trans = []
			for p in occ1:
				occ1_trans.append(gen_icp.dispOffset(p, currPose))	

			hull1_trans = []
			for p in hull1:
				hull1_trans.append(gen_icp.dispOffset(p, currPose))	
							
			xP = []
			yP = []
			for p in hull1_trans:
				xP.append(p[0])
				yP.append(p[1])
			#pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	
			pylab.plot(xP,yP, color='r', linewidth=2)	

			xP = []
			yP = []
			for p in occ1_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.scatter(xP,yP, color=(1.0,0.5,0.5), faceted = False)



			estPose2 = currProfile.convertLocalOffsetToGlobal(offset)
			
			hull2 = PoseGraph.computeBareHull(self.nodeHash[nodeID2], sweep = False)
			hull2.append(hull2[0])

			node2 = self.nodeHash[nodeID2]
			#currPose = node2.getGlobalGPACPose()
			currProfile = Pose(estPose2)
			posture2 = node2.getStableGPACPosture()
			posture2_trans = []
			for p in posture2:
				posture2_trans.append(gen_icp.dispOffset(p, estPose2))	

			occ2 = node2.getOccPoints()
			occ2_trans = []
			for p in occ2:
				occ2_trans.append(gen_icp.dispOffset(p, estPose2))	

			hull2_trans = []
			for p in hull2:
				hull2_trans.append(gen_icp.dispOffset(p, estPose2))	
			
			xP = []
			yP = []
			for p in hull2_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(52/256.,87/256.,159/256.), linewidth=2)

			xP = []
			yP = []
			for p in occ2_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.scatter(xP,yP, color=(192./256.,227./256.,256./256.), faceted = False)
	
			self.drawWalls()
			
			pylab.xlim(currPose[0]-4, currPose[0]+4)					
			pylab.ylim(currPose[1]-3, currPose[1]+3)
			
			if status == 1:
				statusStr = "ORI_FAIL"
			elif status == 2:
				statusStr = "COST_FAIL"
			elif status == 3:
				statusStr = "HYP_FAIL"
			else:
				statusStr = "UNKN_FAIL"

			pylab.title("%d -> %d, status = %s" % (nodeID1, nodeID2, statusStr))
			
			pylab.savefig("plotBadShapeConstraint%04u_%04u.png" % (id, numEdges))

			numEdges += 1


		numEdges = 0
		edgeHash = self.poseGraph.getPriorityEdges(PoseGraph.SHAPE_PRIORITY)

		
		for k, v in edgeHash.items():
			
			print len(v), "edges printing"
			for i in range(len(v)):

				pylab.clf()
				
				if k[0] < k[1]:
					nodeID1 = k[0]
					nodeID2 = k[1]
					transform = v[i][0]
					offset = [transform[0,0], transform[1,0], transform[2,0]]

				else:
					nodeID1 = k[1]
					nodeID2 = k[0]
					transform = v[i][0]
					offset = [transform[0,0], transform[1,0], transform[2,0]]
					prof1 = Pose(offset)
					offset = prof1.convertGlobalPoseToLocal([0.0,0.0,0.0])
	
				hull1 = PoseGraph.computeBareHull(self.nodeHash[nodeID1], sweep = False)
				hull1.append(hull1[0])
	
				node1 = self.nodeHash[nodeID1]
				currPose = node1.getGlobalGPACPose()
				currProfile = Pose(currPose)
				posture1 = node1.getStableGPACPosture()
				posture1_trans = []
				for p in posture1:
					posture1_trans.append(gen_icp.dispOffset(p, currPose))	

				occ1 = node1.getOccPoints()
				occ1_trans = []
				for p in occ1:
					occ1_trans.append(gen_icp.dispOffset(p, currPose))	
	
				hull1_trans = []
				for p in hull1:
					hull1_trans.append(gen_icp.dispOffset(p, currPose))	
								
				xP = []
				yP = []
				for p in hull1_trans:
					xP.append(p[0])
					yP.append(p[1])
				#pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	
				pylab.plot(xP,yP, color='r', linewidth=2)	

				xP = []
				yP = []
				for p in occ1_trans:
					xP.append(p[0])
					yP.append(p[1])
				pylab.scatter(xP,yP, color=(1.0,0.5,0.5), faceted = False)



				estPose2 = currProfile.convertLocalOffsetToGlobal(offset)
				
				hull2 = PoseGraph.computeBareHull(self.nodeHash[nodeID2], sweep = False)
				hull2.append(hull2[0])
	
				node2 = self.nodeHash[nodeID2]
				#currPose = node2.getGlobalGPACPose()
				currProfile = Pose(estPose2)
				posture2 = node2.getStableGPACPosture()
				posture2_trans = []
				for p in posture2:
					posture2_trans.append(gen_icp.dispOffset(p, estPose2))	

				occ2 = node2.getOccPoints()
				occ2_trans = []
				for p in occ2:
					occ2_trans.append(gen_icp.dispOffset(p, estPose2))	
	
				hull2_trans = []
				for p in hull2:
					hull2_trans.append(gen_icp.dispOffset(p, estPose2))	
				
				xP = []
				yP = []
				for p in hull2_trans:
					xP.append(p[0])
					yP.append(p[1])
				pylab.plot(xP,yP, color=(52/256.,87/256.,159/256.), linewidth=2)

				xP = []
				yP = []
				for p in occ2_trans:
					xP.append(p[0])
					yP.append(p[1])
				pylab.scatter(xP,yP, color=(192./256.,227./256.,256./256.), faceted = False)
		
				self.drawWalls()
				
				pylab.xlim(currPose[0]-4, currPose[0]+4)					
				pylab.ylim(currPose[1]-3, currPose[1]+3)

				pylab.title("%d -> %d" % (nodeID1, nodeID2))
				
				#pylab.xlim(-5,10)
				#pylab.ylim(-8,8)
				pylab.savefig("plotShapeConstraint%04u_%04u.png" % (id, numEdges))
	
				numEdges += 1
	
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
				
			
			
			if i == (self.numNodes - 1):
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


		self.drawWalls()


			
		pylab.xlim(-5,10)
		#pylab.xlim(-8,12)
		#pylab.ylim(-10,10)
		pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % i)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)

	def drawWalls(self):
		for wall in self.walls:
			xP = []
			yP = []
			for i in range(len(wall)):
				p = wall[i]
				xP.append(p[0])
				yP.append(p[1])
	
			pylab.plot(xP,yP, linewidth=2, color = 'g')

	def drawMap(self):

		self.occMap.update()
		self.occMap.saveMap()

	def loadWalls(self, walls):
		self.walls = walls

	def getWalls(self):
		return self.walls

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
		

