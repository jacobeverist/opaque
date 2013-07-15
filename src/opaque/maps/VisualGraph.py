import PoseGraph
from LocalNode import LocalNode, getLongestPath
from Pose import Pose
import pylab
import gen_icp
from numpy import matrix
from random import random, gauss
import functions
from copy import copy
from OccupancyMap import OccupancyMap
from SplineFit import SplineFit
import math
import graph
from subprocess import Popen, PIPE
from medialaxis import computeMedialAxis
import Image
import sys
import alphamod

class VisualGraph:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts		
		self.poseGraph = PoseGraph.PoseGraph(self.probe, self.contacts)	

		self.edgeAgeHash = {}

		self.localNodes = []

		self.sensorHypotheses = []
		
		self.colors = []
		self.colors.append([215, 48, 39])
		self.colors.append([252, 141, 89])
		self.colors.append([254, 224, 144])
		self.colors.append([224, 243, 248])
		self.colors.append([145, 191, 219])
		self.colors.append([69, 117, 180])


		for color in self.colors:
			color[0] = float(color[0])/256.0
			color[1] = float(color[1])/256.0
			color[2] = float(color[2])/256.0

		for i in range(100):
			self.colors.append((random(),random(),random()))

		MAPSIZE = 20.0
		self.occMap = OccupancyMap(self.probe, self, MAPSIZE)

	def restoreSeries(self, dirName, num_poses):

		self.poseGraph.restoreState(dirName, num_poses)

		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", i			
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile2(dirName, i)
			self.localNodes.append(currNode)
			currNode.setGPACPose(self.poseGraph.nodePoses[i])
			self.poseGraph.nodeHash[i] = currNode
			#self.drawConstraints()
			#self.poseGraph.mergePriorityConstraints()

			#self.poseGraph.restoreNode(dirName, currNode)		
			#for k in range(len(self.nodePoses)):
			#self.nodeHash[k].setGlobalGPACPose(self.nodePoses[k])

		self.poseGraph.mergePriorityConstraints()

		#self.drawConstraints()
		

		#return
		#self.mapGraph.insertPose(self.probe.getActualJointPose(19), self.travelDir)


		foreNode = LocalNode(self.probe, self.contacts, num_poses, 19, PIXELSIZE)
		foreNode.readFromFile2(dirName, num_poses)
		backNode = LocalNode(self.probe, self.contacts, num_poses+1, 19, PIXELSIZE)
		backNode.readFromFile2(dirName, num_poses+1)
		self.poseGraph.insertPose2(foreNode, backNode, initLocation = foreNode.getEstPose())
		
		
		for i in range(num_poses+2, num_poses+70):
		#for i in range(num_poses+2, num_poses+4):
		#for i in range(num_poses+2, num_poses+26):
			#for i in range(num_poses+2, num_poses+8):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile2(dirName, i)
			self.localNodes.append(currNode)
			self.poseGraph.loadNewNode(currNode)
			self.poseGraph.mergePriorityConstraints()
			#self.drawConstraints(i)
			self.poseGraph.saveState()
		

		self.drawMap()
		#self.drawConstraints(num_poses)

		return

		for i in range(num_poses, num_poses+0):

			print "loading node", self.numNodes			
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)

			currNode.readFromFile2(dirName, i)

			self.localNodes.append(currNode)
			
			self.poseGraph.loadNewNode(currNode)
			self.poseGraph.mergePriorityConstraints()
			self.drawConstraints(i)
			self.poseGraph.saveState()

		self.poseGraph.paths.generatePaths()
		trimPaths = self.poseGraph.paths.trimPaths(self.poseGraph.paths.paths)  		
		self.poseGraph.drawTrimmedPaths(trimPaths)		

		splices, terminals, junctions = self.poseGraph.paths.getAllSplices()

		self.poseGraph.mergePaths()
		#self.drawConstraints()
		
		#self.poseGraph.paths.comparePaths()


		print "consistency:"
		print self.poseGraph.paths.consistency

		xP1 = []
		yP1 = []
		xP2 = []
		yP2 = []
		
		for k, v in terminals.iteritems():
			pnt = v[1]
			xP1.append(pnt[0])
			yP1.append(pnt[1])


		for k, v in junctions.iteritems():
			pnt = v[1]
			xP2.append(pnt[0])
			yP2.append(pnt[1])
		
		spliceCount = 0
		for k, v in splices.iteritems():
			
			for sPath in v:
				
				splicedPath = sPath['path']   
				pylab.clf()
				
				self.drawWalls()
				
				print "spliced path has", len(splicedPath), "points", sPath['orderedPathIDs']
				xP = []
				yP = []
				for p in splicedPath:
					xP.append(p[0])
					yP.append(p[1])
				
				pylab.plot(xP,yP)
				
				pylab.scatter(xP1,yP1, color='k')
				pylab.scatter(xP2,yP2, color='r')
				
				
				pylab.xlim(-5,10)
				pylab.ylim(-8,8)			
				pylab.title("Spliced Paths %d " % spliceCount)
				pylab.savefig("splicedPath_%04u.png" % spliceCount)
				spliceCount += 1
		print "len(splices) =", len(splices)

		self.poseGraph.drawPathAndHull()

		
		return


	def loadSeries(self, dirName, num_poses):
		
	
		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)

			#if i > 0 and i % 2 == 0:
			if False:	
				" since estimated poses are modified from the original motion estimation, we need to restore them "

				f = open(dirName + "/motion_constraints_%04u.txt" % (i+1), 'r')
				#motion_constraints = eval(f.read().rstrip())
				motion_constraints = eval(f.read().replace('\r\n','\n'))
				f.close()				
				transform = motion_constraints[-1][2]
				offset = [transform[0,0], transform[1,0], transform[2,0]]
				
				estPose1 = self.localNodes[i-1].getEstPose()
				profile1 = Pose(estPose1)
				estPose2 = profile1.convertLocalOffsetToGlobal(offset)
				
				currNode.readFromFile(dirName, i, forcedPose = estPose2)
			else:
				currNode.readFromFile2(dirName, i)

			self.localNodes.append(currNode)
			

			self.poseGraph.loadNewNode(currNode)
			#self.poseGraph.resetGraphToGround()						
			#self.poseGraph.loadNodeGroundPast(currNode)
	
			#if i % 10 == 0:
			#	self.poseGraph.makeCornerBinConsistent()

			#self.poseGraph.addCornerConstraints(i)			
			self.poseGraph.mergePriorityConstraints()

	
			#self.drawMotion(i)
			#self.drawConstraints(i)
			#self.drawTopology2(i)
			#self.drawMedialPath(i)
			
			#self.drawShapeConstraints(i)

			self.poseGraph.saveState()
		

		self.drawMap()
		#self.drawConstraints(num_poses)
		#self.drawPoses()
		
		#self.drawAllConstraints()
		
		#exit()

		#self.poseGraph.performOverlapConstraints()
		#self.drawMap()
		#self.drawConstraints(num_poses)
		#self.drawTopology2(num_poses)
		#self.drawMedialPath(num_poses)
		"""
		self.poseGraph.performShapeConstraints()
		self.drawMap()
		self.drawConstraints(num_poses+1)
		"""		
		
		#for i in range(0, num_poses):
		#	self.poseGraph.addCornerConstraints(i)			
		#	self.poseGraph.mergePriorityConstraints()
		
		#self.drawMap()
		#self.drawConstraints(num_poses)
		#self.drawTopology2(num_poses+1)
		#self.drawMedialPath(num_poses+1)
		
		#self.drawAllConstraints()

		#self.drawTopology()
		#self.poseGraph.resetGraphToGround()
		#self.drawMap()
		#self.drawConstraints(num_poses+2)
		#self.drawMedialPath(num_poses+2)
		#
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

			if self.nodeHash[i].isBowtie:			
				hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False, static = True)
				hull.append(hull[0])
				medial = self.nodeHash[i].getStaticMedialAxis()
			else:
				hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False)
				hull.append(hull[0])
				medial = self.nodeHash[i].getMedialAxis(sweep = False)

			#hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False, static = True)
			#hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False)
			#hull = PoseGraph.computeBareHull(self.nodeHash[6], sweep = False)
			#hull.append(hull[0])

			#medial = gen_icp.getMedialAxis(hull)
			#medial = self.nodeHash[i].getMedialAxis()
			#medial = self.nodeHash[i].getStaticMedialAxis()

			
			edge1 = medial[0:2]
			edge2 = medial[-2:]
			
			frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
			backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
			frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
			backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
			
			frontVec[0] /= frontMag
			frontVec[1] /= frontMag
			backVec[0] /= backMag
			backVec[1] /= backMag
			
			newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
			newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

			edge1 = [newP1, edge1[1]]
			edge2 = [edge2[0], newP2]

			
			interPoints = []
			for j in range(len(hull)-1):
				hullEdge = [hull[j],hull[j+1]]
				isIntersect1, point1 = functions.Intersect(edge1, hullEdge)
				if isIntersect1:
					#print "found", point1, "while intersecting"
					#print edge1
					#print hullEdge
					#print
					interPoints.append(point1)
					break

			for j in range(len(hull)-1):
				hullEdge = [hull[j],hull[j+1]]
				isIntersect2, point2 = functions.Intersect(edge2, hullEdge)
				if isIntersect2:
					#print "found", point2, "while intersecting"
					#print edge2
					#print hullEdge
					#print
					interPoints.append(point2)
					break
				
			medial = medial[1:-2]
			if isIntersect1:
				medial.insert(0, point1)
			if isIntersect2:
				medial.append(point2)
			
			
			medialSpline = SplineFit(medial, smooth=0.1)
			medialPoints = medialSpline.getUniformSamples()

			node1 = self.nodeHash[i]
			#node1 = self.nodeHash[6]
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
			
			medial_trans = []
			for p in medial:
				medial_trans.append(gen_icp.dispOffset(p, currPose))	
			
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

			xP = []
			yP = []
			for p in medialPoints:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(0.5,0.5,1.))	

			xP = []
			yP = []
			for p in medial:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(0.,0.,1.))	

			for edge in node1.lineEdges:
				xP = [edge[0][0],edge[1][0]]
				yP = [edge[0][1],edge[1][1]]
				pylab.plot(xP,yP, color=(1.,0.,1.))	
				
			xP = []
			yP = []
			currX = 0.0
			totalWidth = 0.0
			totalThresh = 0
			for samp in node1.medialWidth:
				#currU = samp[2]
				width = 1.0 + samp[3] + samp[4]
				totalWidth += samp[3] + samp[4]
				if samp[3] + samp[4] > 0.5:
					totalThresh += 1
				xP.append(currX - 3.0)
				yP.append(width)
				currX += 0.04
				#self.medialWidth.append([linePoint[0],linePoint[1],currU, distR, distL])
			pylab.plot(xP,yP, color=(0.0,0.0,1.0))


			#xP = []
			#yP = []
			#for p in interPoints:
			#	xP.append(p[0])
			#	yP.append(p[1])
			#pylab.scatter(xP,yP, color=(0.,0.,1.))	

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
			#pylab.title("count1 = %d, count2 = %d, totalWid = %f, totalThresh = %d" % (count1,count2, totalWidth, totalThresh))
			pylab.title("totalWid = %f, totalThresh = %d" % (totalWidth, totalThresh))
			pylab.xlim(-3,3)
			pylab.ylim(-2,2)
			pylab.savefig("plotPose%04u.png" % i)
			#pylab.show()

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
			
			#print len(v), "edges printing"
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


	def drawTopology2(self, id = []):
		
		realPath = self.poseGraph.getTopology()
		pylab.clf()
		self.drawWalls()

		xP = []
		yP = []
		for p in realPath:
			xP.append(p[0])
			yP.append(p[1])

		pylab.plot(xP,yP, color=(0.0,0.0,1.0))
				

		pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotTopology%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotTopology%04u.png" % id)

	def drawTopology(self):

		random.seed(0)
		
		def convertAlphaUniform(a_vert, max_spacing = 0.04):
			
			" make the vertices uniformly distributed "
			
			new_vert = []
		
			for i in range(len(a_vert)):
				p0 = a_vert[i]
				p1 = a_vert[(i+1) % len(a_vert)]
				dist = math.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
	
				vec = [p1[0]-p0[0], p1[1]-p0[1]]
				vec[0] /= dist
				vec[1] /= dist
				
				new_vert.append(copy(p0))
				
				if dist > max_spacing:
					" cut into pieces max_spacing length or less "
					numCount = int(math.floor(dist / max_spacing))
					
					for j in range(1, numCount+1):
						newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
						new_vert.append(newP)
			
			return new_vert			

		medialPointSoup = []

		for nodeID in range(0, self.poseGraph.numNodes):

			estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()		
	
			if self.nodeHash[nodeID].isBowtie:			
				hull1 = PoseGraph.computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
				#hull1.append(hull1[0])
			else:
				hull1 = PoseGraph.computeBareHull(self.nodeHash[nodeID], sweep = False)
				#hull1.append(hull1[0])

			" set the origin of pose 1 "
			poseOrigin = Pose(estPose1)

			#medialPath1 = self.nodeHash[nodeID].medialPathCut		
			#medialUnif = convertAlphaUniform(medialPath1)

			for p in hull1:
				p1 = poseOrigin.convertLocalToGlobal(p)
				medialPointSoup.append(p1)

			pylab.clf()
			self.drawWalls()

			xP = []
			yP = []
			for p in medialPointSoup:
				xP.append(p[0])
				yP.append(p[1])

			pylab.scatter(xP,yP, color=(1.0,0.5,0.5), faceted = False)

			radius = 0.2

			numPoints = len(medialPointSoup)
			isDone = False
			
			while not isDone:
		
				perturbPoints = []
				
				for p in medialPointSoup:
					p2 = copy(p)
					" add a little bit of noise to avoid degenerate conditions in CGAL "
					p2[0] += random.gauss(0.0,0.000001)
					p2[1] += random.gauss(0.0,0.000001)
		
					perturbPoints.append(p2)
			
				try:			
		
					vertices = alphamod.doAlpha(radius,perturbPoints)
					numVert = len(vertices)
					
					if numVert <= 2:
						print "Failed, hull had only", numVert, "vertices"
						raise
					
					isDone = True
				except:
					print "hull has holes!  retrying..."
					#print sArr	
	
			return vertices


			"""
			numPoints = len(medialPointSoup)
			inputStr = str(numPoints) + " "
	
			" alpha shape circle radius "
			inputStr += str(radius) + " "
			
			#for p in medialPointSoup:
			#	p2 = copy(p)
			#	" add a little bit of noise to avoid degenerate conditions in CGAL "	
			#	inputStr += str(p2[0]) + " " + str(p2[1]) + " "
			#
			#inputStr += "\n"

			isDone = False
			
			while not isDone:
	
				inputStr = str(numPoints) + " "
		
				" alpha shape circle radius "
				inputStr += str(radius) + " "
				
				for p in medialPointSoup:
					p2 = copy(p)
					" add a little bit of noise to avoid degenerate conditions in CGAL "
					p2[0] += random.gauss(0.0,0.000001)
					p2[1] += random.gauss(0.0,0.000001)
		
					inputStr += str(p2[0]) + " " + str(p2[1]) + " "
				
				inputStr += "\n"
				
				try:			
					" start the subprocess "
					if sys.platform == "win32":
						subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
					else:
						subProc = Popen(["./alpha2"], stdin=PIPE, stdout=PIPE)
						
					
					" send input and receive output "
					sout, serr = subProc.communicate(inputStr)
			
					" convert string output to typed data "
					sArr = sout.split(" ")
					
					numVert = int(sArr[0])
					
					sArr = sArr[1:]
					
					
					vertices = []
					for i in range(len(sArr)/2):
						vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
					isDone = True
				except:
					print "hull has holes!  retrying..."
					#print sArr	
			"""
			"""
			vertices = []
			try:			
				" start the subprocess "
				subProc = Popen(["./alpha2.exe"], stdin=PIPE, stdout=PIPE)
				
				print subProc.stdin, subProc.stderr, subProc.stdout
				#print "input:", inputStr
				" send input and receive output "
				sout, serr = subProc.communicate(inputStr)

				" convert string output to typed data "
				sArr = sout.split(" ")
				#print sArr
				
				numVert = int(sArr[0])
				
				sArr = sArr[1:]
				
				vertices = []
				for i in range(len(sArr)/2):
					vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
			except:
				print "hull has holes!"
			"""

			" cut out the repeat vertex "
			vertices = vertices[:-1]
			
			vertices = convertAlphaUniform(vertices)

			vertices.append(vertices[0])
			
			minX = 1e100
			maxX = -1e100
			minY = 1e100
			maxY = -1e100
			for p in vertices:
				if p[0] > maxX:
					maxX = p[0]
				if p[0] < minX:
					minX = p[0]
				if p[1] > maxY:
					maxY = p[1]
				if p[1] < minY:
					minY = p[1]
		
		
			PIXELSIZE = 0.05
			mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
			pixelSize = PIXELSIZE
			numPixel = int(2.0*mapSize / pixelSize + 1.0)
			divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)
			
			def realToGrid(point):
				indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
				indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
				return indexX, indexY
		
			def gridToReal(indices):
				i = indices[0]
				j = indices[1]
				point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
				return point
		
			gridHull = []
			for i in range(len(vertices)):
				p = vertices[i]
				gridHull.append(realToGrid(p))
		
			minX = 1e100
			maxX = -1e100
			minY = 1e100
			maxY = -1e100
			for p in gridHull:
				if p[0] > maxX:
					maxX = p[0]
				if p[0] < minX:
					minX = p[0]
				if p[1] > maxY:
					maxY = p[1]
				if p[1] < minY:
					minY = p[1]

			resultImg = Image.new('L', (numPixel,numPixel))
			resultImg = computeMedialAxis(nodeID, numPixel,numPixel, 0, resultImg, len(gridHull[:-2]), gridHull[:-2])

			imgA = resultImg.load()
		
			points = []
			for i in range(1, numPixel-1):
				for j in range(1, numPixel-1):
					if imgA[i,j] == 0:
						points.append((i,j))
		
			medialGraph = graph.graph()
			for p in points:
				medialGraph.add_node(p, [])
		
			builtGraph = {}
			for i in range(2, numPixel-2):
				for j in range(2, numPixel-2):
					if imgA[i,j] == 0:
						builtGraph[(i,j)] = []
						for k in range(i-1, i+2):
							for l in range(j-1, j+2):
								if imgA[k,l] == 0:
									builtGraph[(i,j)].append((k,l))
									medialGraph.add_edge((i,j), (k,l))
									
			mst = medialGraph.minimal_spanning_tree()
		
			uni_mst = {}
			isVisited = {}
			nodeSum = {}
			for k, v in mst.items():
				uni_mst[k] = []
				isVisited[k] = 0
				nodeSum[k] = 0
		
			for k, v in mst.items():
				if v != None:
					uni_mst[k].append(v)
					uni_mst[v].append(k)
		
			leaves = []
			for k, v in uni_mst.items():
				if len(v) == 1:
					leaves.append(k)
		
		
			" find the longest path from each leaf"
			maxPairDist = 0
			maxPair = None
			maxPath = []
			for leaf in leaves:
		
				isVisited = {}
				nodeSum = {}
				nodePath = {}
				for k, v in uni_mst.items():
					isVisited[k] = 0
					nodeSum[k] = 0
					nodePath[k] = []
		
				getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
		
				maxDist = 0
				maxNode = None
				for k, v in nodeSum.items():
					#print k, v
					if v > maxDist:
						maxNode = k
						maxDist = v
		
				#print leaf, "-", maxNode, maxDist
		
				if maxDist > maxPairDist:
					maxPairDist = maxDist
					maxPair = [leaf, maxNode]
					maxPath = nodePath[maxNode]

			" convert path points to real "
			realPath = []
			for p in maxPath:
				realPath.append(gridToReal(p))			

			#print "maxPath:", maxPath
			xP = []
			yP = []
			for p in vertices:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color=(0.0,0.0,1.0))

			xP = []
			yP = []
			for p in realPath:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP,yP, color=(0.0,0.0,1.0))
					
	
			pylab.xlim(-5,10)
			pylab.ylim(-8,8)
			pylab.savefig("plotTopology%04u.png" % nodeID)

	
	def drawAllConstraints(self):

		overlaps = self.poseGraph.getPriorityEdges(priorityLevel = PoseGraph.OVERLAP_PRIORITY)
		inplaces = self.poseGraph.getPriorityEdges(priorityLevel = PoseGraph.INPLACE_PRIORITY)
		
		allConstraints = self.poseGraph.getPriorityEdges()
		
		#self.edgePriorityHash[(edge[0],edge[1])].append([edge[2], edge[3], priority])
		#self.addPriorityEdge([nodeID-2,nodeID,transform,covE], OVERLAP_PRIORITY)

		const_count = 0

		for k, v in allConstraints.items():

			if len(v) > 0:
				#print k, v
				val = v[0]
				nodeID1 = k[0]
				nodeID2 = k[1]
				transform = val[0]
				covE = val[1]
				priority = val[2]
				
				estPose1 = self.nodeHash[nodeID1].getGlobalGPACPose()		
	
			
				if self.nodeHash[nodeID1].isBowtie:			
					hull1 = PoseGraph.computeBareHull(self.nodeHash[nodeID1], sweep = False, static = True)
					hull1.append(hull1[0])
				else:
					hull1 = PoseGraph.computeBareHull(self.nodeHash[nodeID1], sweep = False)
					hull1.append(hull1[0])
	
				if self.nodeHash[nodeID2].isBowtie:			
					hull2 = PoseGraph.computeBareHull(self.nodeHash[nodeID2], sweep = False, static = True)
					hull2.append(hull2[0])
				else:
					hull2 = PoseGraph.computeBareHull(self.nodeHash[nodeID2], sweep = False)
					hull2.append(hull2[0])
	
	
				offset = [transform[0,0], transform[1,0], transform[2,0]]
	
				#hull2_trans = []
				#for p in hull2:
				#	hull2_trans.append(gen_icp.dispOffset(p, offset))	
	
				" set the origin of pose 1 "
				poseOrigin = Pose(estPose1)
		
				pylab.clf()
				pylab.axes()
				
				xP = []
				yP = []
				for b in hull1:
					p1 = poseOrigin.convertLocalToGlobal(b)
		
					xP.append(p1[0])	
					yP.append(p1[1])
		
				pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
	
	
				xP = []
				yP = []
				for p in hull2:
					p1 = gen_icp.dispOffset(p,offset)
					
					p2 = poseOrigin.convertLocalToGlobal(p1)
					xP.append(p2[0])	
					yP.append(p2[1])
	
				pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))


				" create the ground constraints "
				gndGPAC1Pose = self.nodeHash[nodeID1].getGndGlobalGPACPose()
				currProfile = Pose(gndGPAC1Pose)
				gndGPAC2Pose = self.nodeHash[nodeID2].getGndGlobalGPACPose()
				gndOffset = currProfile.convertGlobalPoseToLocal(gndGPAC2Pose)

				xP = []
				yP = []
				for p in hull2:
					p1 = gen_icp.dispOffset(p,gndOffset)
					
					p2 = poseOrigin.convertLocalToGlobal(p1)
					xP.append(p2[0])	
					yP.append(p2[1])
	
				pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))


				medialPath1 = self.nodeHash[nodeID1].medialPathCut
				medialPath2 = self.nodeHash[nodeID2].medialPathCut
				medial_trans = []
				for p in medialPath2:
					medial_trans.append(gen_icp.dispOffset(p, offset))	

				xP = []
				yP = []
				for p in medialPath1:
					p1 = poseOrigin.convertLocalToGlobal(p)
					xP.append(p1[0])
					yP.append(p1[1])
				pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

				xP = []
				yP = []
				for p in medial_trans:
					p1 = poseOrigin.convertLocalToGlobal(p)
					xP.append(p1[0])
					yP.append(p1[1])
				pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

				xP = []
				yP = []
				for p in medialPath2:
					p1 = gen_icp.dispOffset(p, gndOffset)
					p2 = poseOrigin.convertLocalToGlobal(p1)
					xP.append(p2[0])
					yP.append(p2[1])
				pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))
	
				self.drawWalls()

				typeStr = "UNKNOWN:"
				if priority == PoseGraph.OVERLAP_PRIORITY:
					typeStr = "OVERLAP:"
				if priority == PoseGraph.INPLACE_PRIORITY:
					typeStr = "INPLACE:"
				if priority == PoseGraph.CORNER_PRIORITY:
					typeStr = "CORNER:"
				if priority == PoseGraph.SHAPE_PRIORITY:
					typeStr = "SHAPE:"

				pylab.title(typeStr + "  %u -> %u, covar: %1.5f %1.5f %1.5f" % (nodeID1, nodeID2, covE[0,0],covE[1,1],covE[2,2]))

	
				pylab.xlim(estPose1[0]-4, estPose1[0]+4)					
				pylab.ylim(estPose1[1]-3, estPose1[1]+3)

				if priority == PoseGraph.OVERLAP_PRIORITY:
					pylab.savefig("constraint_%03u_%03u_overlap.png" % (nodeID1,nodeID2))
				elif priority == PoseGraph.INPLACE_PRIORITY:
					pylab.savefig("constraint_%03u_%03u_inplace.png" % (nodeID1,nodeID2))
				elif priority == PoseGraph.CORNER_PRIORITY:
					pylab.savefig("constraint_%03u_%03u_corner.png" % (nodeID1,nodeID2))
				elif priority == PoseGraph.SHAPE_PRIORITY:
					pylab.savefig("constraint_%03u_%03u_shape.png" % (nodeID1,nodeID2))				
				else:
					pylab.savefig("constraint_%03u_%03u_unknown.png" % (nodeID1,nodeID2))
					
				#pylab.savefig("constraints_%04u.png" % const_count)
				pylab.clf()			
	
				const_count += 1

		return
	
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())
	
		cornerTrans = []

		pylab.clf()
		for i in range(self.numNodes):

			if self.nodeHash[i].isBowtie:			
				hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False, static = True)
				hull.append(hull[0])
			else:
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
			
			
			medialPath = self.nodeHash[i].medialPathCut
			medial_trans = []
			for p in medialPath:
				medial_trans.append(gen_icp.dispOffset(p, currPose))	
			
			" extract corner points for each node "
			cornerCandidates = node1.cornerCandidates
			cornerPoints = []
			for cand in cornerCandidates:
				cornerPoints.append(cand[0])
			
			for p in cornerPoints:
				cornerTrans.append(gen_icp.dispOffset(p, currPose))	

			xP = []
			yP = []
			for p in medial_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	


			xP = []
			yP = []
			for p in hull_trans:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(112/256.,147/256.,219/256.))	


		self.drawWalls()
			
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)



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
			
			#print len(v), "edges printing"
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

			#hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False)
			#hull = PoseGraph.computeBareHull(self.nodeHash[i], static= True)
			#hull.append(hull[0])
			if self.nodeHash[i].isBowtie:			
				hull = PoseGraph.computeBareHull(self.nodeHash[i], sweep = False, static = True)
				hull.append(hull[0])
			else:
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
			
			
			medialPath = self.nodeHash[i].medialPathCut
			medial_trans = []
			for p in medialPath:
				medial_trans.append(gen_icp.dispOffset(p, currPose))	
			
			" extract corner points for each node "
			cornerCandidates = node1.cornerCandidates
			cornerPoints = []
			for cand in cornerCandidates:
				cornerPoints.append(cand[0])
			
			for p in cornerPoints:
				cornerTrans.append(gen_icp.dispOffset(p, currPose))	
				
			
			
			#if i == (self.numNodes - 1):
			#	xP = []
			#	yP = []
			#	for p in posture_trans:
			#		xP.append(p[0])
			#		yP.append(p[1])
			#	pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

			if i >= self.numNodes-2:
				xP = []
				yP = []
				for p in medial_trans:
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


			
		#pylab.xlim(-5,10)
		#pylab.xlim(-8,12)
		#pylab.ylim(-10,10)
		#pylab.ylim(-8,8)
		pylab.xlim(-10,10)
		pylab.ylim(-10,10)

		
		if id == []:
			pylab.savefig("plotEstimate%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotEstimate%04u.png" % id)

	
	def drawMotion(self, id = []):
		

		"""
		f = open(dirName + "/motion_constraints_%04u.txt" % (i+1), 'r')
		motion_constraints = eval(f.read().replace('\r\n','\n'))
		f.close()				
		transform = motion_constraints[-1][2]
		offset = [transform[0,0], transform[1,0], transform[2,0]]
		
		estPose1 = self.poseGraph.nodeHash[i-1].getEstPose()
		profile1 = Pose(estPose1)
		estPose2 = profile1.convertLocalOffsetToGlobal(offset)
		
		currNode.readFromFile(dirName, i, forcedPose = estPose2)
		"""
		
		print "drawing", len(self.localNodes), "nodes"
		

		pylab.clf()
		for i in range(len(self.localNodes)):

			node1 = self.localNodes[i]
			posture1 = node1.getStableGPACPosture()
				
			gndPose = node1.getGndGlobalGPACPose()
			posture_gnd = []
			for p in posture1:
				posture_gnd.append(gen_icp.dispOffset(p, gndPose))
							
			xP = []
			yP = []
			for p in posture_gnd:
				xP.append(p[0])
				yP.append(p[1])
			#pylab.plot(xP,yP, color=(.5,0.2,0.2))	
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

		for i in range(len(self.localNodes)):

			node1 = self.localNodes[i]
			currPose = node1.getGlobalGPACPose()
			currProfile = Pose(currPose)
			posture1 = node1.getStableGPACPosture()
			posture_trans = []
			for p in posture1:
				posture_trans.append(gen_icp.dispOffset(p, currPose))
				
			
			xP = []
			yP = []
			for p in posture_trans:
				xP.append(p[0])
				yP.append(p[1])
			#pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	
			pylab.plot(xP,yP, color=(0,0,0))	

		self.drawWalls()


		pylab.xlim(-3,6)
		pylab.ylim(-4,4)
		
		#pylab.xlim(-3,12)
		#pylab.ylim(-4,12)

		#pylab.xlim(-5,10)
		#pylab.xlim(-8,12)
		#pylab.ylim(-10,10)
		#pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotMotion%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotMotion%04u.png" % id)

	def drawMedialPath(self, id = []):
		
		poses = []
		for i in range(self.numNodes):
			poses.append(self.nodeHash[i].getGlobalGPACPose())
	
		medials = []

		overlaps = self.poseGraph.getPriorityEdges(priorityLevel = PoseGraph.OVERLAP_PRIORITY)
		inplaces = self.poseGraph.getPriorityEdges(priorityLevel = PoseGraph.INPLACE_PRIORITY)
		shapes = self.poseGraph.getPriorityEdges(priorityLevel = PoseGraph.SHAPE_PRIORITY)
		pairs = []
		
		for k, v in overlaps.items():
			if len(v) > 0:
				nodeID1 = k[0]
				nodeID2 = k[1]
				pairs.append((nodeID1,nodeID2))

		for k, v in inplaces.items():
			if len(v) > 0:
				nodeID1 = k[0]
				nodeID2 = k[1]
				pairs.append((nodeID1,nodeID2))

		for k, v in shapes.items():
			if len(v) > 0:
				nodeID1 = k[0]
				nodeID2 = k[1]
				pairs.append((nodeID1,nodeID2))

			
		for i in range(self.numNodes):
				
			node1 = self.nodeHash[i]
			currPose = node1.getGlobalGPACPose()
			#currProfile = Pose(currPose)

			medialPath = self.nodeHash[i].medialPathCut
			medial_trans = []
			for p in medialPath:
				medial_trans.append(gen_icp.dispOffset(p, currPose))	
			
			medials.append(medial_trans)		
		
		def getIntersection(medial1, medial2):
			
			for i in range(len(medial1)-1):
				seg1 = [medial1[i],medial1[i+1]]
				for j in range(len(medial2)-1):
					seg2 = [medial2[j],medial2[j+1]]
			
					result, point = functions.Intersect(seg1, seg2)
					
					if result:
						return i, i+1, j, j+1, point
					
			raise


		def performSplices(splices, medials):
			
			medialGraph = graph.graph()
			
			for i in range(len(medials)):
				for j in range(len(medials[i])):
					medialGraph.add_node((i,j),medials[i][j])

				for j in range(len(medials[i])-1):
					medialGraph.add_edge((i,j), (i,j+1))

			" perform splices "
						
			for i in range(len(splices)):
				splice = splices[i]
				#medialSplices.append([n1,n2,n1_i,n1_j,n2_i,n2_j, pnt])
				n1 = splice[0]
				n2 = splice[1]
				n1_i = splice[2]
				n1_j = splice[3]
				n2_i = splice[4]
				n2_j = splice[5]
				splicePoint = splice[6]
				
				#print "splicing:", splice
				
				" remove intersecting edges "
				if medialGraph.has_edge((n1,n1_i),(n1,n1_j)):
					medialGraph.del_edge((n1,n1_i),(n1,n1_j))
				if medialGraph.has_edge((n2,n2_i),(n2,n2_j)):
					medialGraph.del_edge((n2,n2_i),(n2,n2_j))

				" add splice point "
				medialGraph.add_node((-1,i),splicePoint)
				
				" add the splicing edges "
				medialGraph.add_edge((n1,n1_i), (-1,i))
				medialGraph.add_edge((n1,n1_j), (-1,i))
				medialGraph.add_edge((n2,n2_i), (-1,i))
				medialGraph.add_edge((n2,n2_j), (-1,i))

			
			#medialSplices.append([n1,n2,n1_i,n1_j,n2_i,n2_j, pnt])		
			#removeEdge(n1_i,n1_j)
			#removeEdge(n2_i,n2_j)
			#addEdge(n1_i,pnt)
			#addEdge(n1_j,pnt)
			#addEdge(n2_i,pnt)
			#addEdge(n2_j,pnt)									

			def getLongestPath(leaf, medialGraph, node, currSum, currPath, tree, isVisited, nodePath, nodeSum):

				#getLongestPath(medialGraph, leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
			
				if isVisited[node]:
					return
			
				isVisited[node] = 1
				nodePath[node] = copy(currPath) + [node]
			
				if nodeSum[node] < currSum:
					nodeSum[node] = currSum
			
				for childNode in tree[node]:
					pnt1 = medialGraph.get_node_attributes(leaf)
					pnt2 = medialGraph.get_node_attributes(childNode)

					newDist = math.sqrt((pnt1[0]+pnt2[0])**2 + (pnt1[1]+pnt2[1])**2)
					getLongestPath(leaf, medialGraph, childNode, newDist, nodePath[node], tree, isVisited, nodePath, nodeSum)
					#getLongestPath(leaf, medialGraph, childNode, currSum + newDist, nodePath[node], tree, isVisited, nodePath, nodeSum)
			
				isVisited[node] = 0
			
			mst = medialGraph.minimal_spanning_tree()
			#mst = medialGraph
		
			uni_mst = {}
			isVisited = {}
			#nodeSum = {}
			for k, v in mst.items():
				#for k in medialGraph.nodes():
				uni_mst[k] = []
				isVisited[k] = 0
				#nodeSum[k] = 0
		
			for k, v in mst.items():
				#for k in medialGraph.nodes():
				#v = medialGraph.get_node_attributes(k)
				#neighbs = medialGraph.neighbors(k)
				#for v in neighbs:
				if v != None:
					#print k, v
					uni_mst[k].append(v)
					uni_mst[v].append(k)
		
			leaves = []
			for k, v in uni_mst.items():
				if len(v) == 1:
					leaves.append(k)
					
			maxPairDist = 0.0
			maxPair = []
			for i in range(len(leaves)):
				pnt1 = medialGraph.get_node_attributes(leaves[i])
				for j in range(i+1, len(leaves)):
					pnt2 = medialGraph.get_node_attributes(leaves[j])
		
					currDist = math.sqrt((pnt1[0]-pnt2[0])**2 + (pnt1[1]-pnt2[1])**2)
					if currDist > maxPairDist:
						maxPairDist = currDist
						maxPair = [leaves[i],leaves[j]]

			shortestPathSpanTree, shortestDist = medialGraph.shortest_path(maxPair[0])
			
			print "maxPair:", maxPair
			print "pair distance =", maxPairDist
			print "shortestDist =", shortestDist[maxPair[1]]
			#print shortestPathSpanTree
			
			currNode = maxPair[1]
			maxPath = [currNode]
			while currNode != maxPair[0]:
				#print "currNode: ", currNode
				nextNode = shortestPathSpanTree[currNode]
				maxPath.append(nextNode)
				currNode = nextNode

			return maxPath, medialGraph, maxPairDist

		
			" find the longest path from each leaf"
			maxPairDist = 0
			maxPair = None
			maxPath = []
			print "leaves:", leaves
			for leaf in leaves:
		
				isVisited = {}
				nodeSum = {}
				nodePath = {}
				for k, v in uni_mst.items():
					isVisited[k] = 0
					nodeSum[k] = 0
					nodePath[k] = []
		
				getLongestPath(leaf, medialGraph, leaf, 0.0, [], uni_mst, isVisited, nodePath, nodeSum)
		
				maxDist = 0
				maxNode = None
				for k, v in nodeSum.items():
					#print k, v
					if v > maxDist:
						maxNode = k
						maxDist = v
		
				#print leaf, "-", maxNode, maxDist
		
				if maxDist > maxPairDist:
					maxPairDist = maxDist
					maxPair = [leaf, maxNode]
					maxPath = nodePath[maxNode]
			
			return maxPath, medialGraph, maxPairDist
		
		medialSplices = []
		for pair in pairs:
			
			n1 = pair[0]
			n2 = pair[1]
			try:
				n1_i, n1_j, n2_i, n2_j, pnt = getIntersection(medials[n1],medials[n2])
				medialSplices.append([n1,n2,n1_i,n1_j,n2_i,n2_j, pnt])
			except:
				pass
		
			#removeEdge(n1_i,n1_j)
			#removeEdge(n2_i,n2_j)
			#addEdge(n1_i,pnt)
			#addEdge(n1_j,pnt)
			#addEdge(n2_i,pnt)
			#addEdge(n2_j,pnt)
		
		
		
		#print "splices:"
		#for splice in medialSplices:
		#	print splice
		
		maxPath, medialGraph, maxDist = performSplices(medialSplices, medials)
		pntPath = []
		for ind in maxPath:
			pntPath.append(medialGraph.get_node_attributes(ind))
		#print "maxPath:", pntPath


		pylab.clf()
		for i in range(self.numNodes):			
			xP = []
			yP = []
			for p in medials[i]:
				xP.append(p[0])
				yP.append(p[1])
			pylab.plot(xP,yP, color=(255/256.,182/256.,193/256.))	

		xP = []
		yP = []
		for p in pntPath:			
			xP.append(p[0])
			yP.append(p[1])
		pylab.plot(xP,yP, color=(0.,0.,1.))	
		pylab.title("maxDist = %3.2f, maxPair = %s %s" % (maxDist,repr(maxPath[0]), repr(maxPath[-1])))

		self.drawWalls()
			
		pylab.xlim(-5,10)
		pylab.ylim(-8,8)
		if id == []:
			pylab.savefig("plotMedialPath%04u.png" % self.numNodes)
		else:
			pylab.savefig("plotMedialPath%04u.png" % id)


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
		self.poseGraph.walls = walls

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
		#if name == "merged_constraints":
		#	return self.poseGraph.merged_constraints

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

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = direction)
		
		self.poseGraph.addInPlaceNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
	
	def newNode(self, stepDist, direction):

		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = True)		
		
		self.poseGraph.addNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

	def initNode(self, estPose, direction):

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = True)

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
		

