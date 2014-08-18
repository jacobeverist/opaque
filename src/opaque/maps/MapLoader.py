#import PoseGraph
import BayesMapper
from LocalNode import LocalNode, getLongestPath
from Pose import Pose
import pylab
import gen_icp
from numpy import matrix
from random import random, gauss
import functions
from copy import copy
from SplineFit import SplineFit
import math
import graph
from subprocess import Popen, PIPE
from medialaxis import computeMedialAxis
import Image
import sys
import alphamod

import cPickle as pickle

class MapLoader:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts		
		#self.mapAlgorithm = PoseGraph.PoseGraph(self.probe, self.contacts)	
		self.mapAlgorithm = BayesMapper.BayesMapper(self.probe.getWalls())	

		print "mapAlgorithm:", self

		self.localNodes = []

		MAPSIZE = 20.0
		#self.occMap = OccupancyMap(self.probe, self, MAPSIZE)

	def restorePickle(self, dirName, nodeID):
		PIXELSIZE = 0.05
		print self

		with open('map_%04u.obj' % nodeID, 'rb') as inputVal:
			print pickle
			print inputVal
			print self
			print self.mapAlgorithm
			self.mapAlgorithm = pickle.load(inputVal)
			print self.mapAlgorithm
			self.mapAlgorithm.saveState()

			for val in self.mapAlgorithm.mapHyps.values():
				self.mapAlgorithm.drawPathAndHull2(val)

			#self.mapAlgorithm.loadSeries("../results/result_2013_08_24_cross", 238)

			"""
			print "loading node", 12		
			currNode = LocalNode(self.probe, self.contacts, 12, 19, PIXELSIZE)
			currNode.readFromFile2("../results/result_2013_08_24_cross", 12)
			self.localNodes.append(currNode)
			self.mapAlgorithm.loadNewNode(currNode)
			self.mapAlgorithm.saveState()
		
			print "loading node", 13		
			currNode = LocalNode(self.probe, self.contacts, 13, 19, PIXELSIZE)
			currNode.readFromFile2("../results/result_2013_08_24_cross", 13)
			self.localNodes.append(currNode)
			self.mapAlgorithm.loadNewNode(currNode)
			self.mapAlgorithm.saveState()
			"""
		

	def loadSeries(self, dirName, startID, endID):
	
		PIXELSIZE = 0.05
		#initNum = len(self.localNodes)

		for i in range(startID, endID+1):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)

			currNode.readFromFile2(dirName, i)

			self.localNodes.append(currNode)
			
			self.mapAlgorithm.loadNewNode(currNode)
			self.mapAlgorithm.saveState()

			#for val in self.mapAlgorithm.mapHyps.values():
			#	self.mapAlgorithm.drawPathAndHull2(val)
		
 
			with open("map_%04u.obj" % i, 'wb') as output:
				pickle.dump(self.mapAlgorithm, output, pickle.HIGHEST_PROTOCOL)

		#self.drawMap()


	def loadDataSeries(self, dirName, num_poses):
	
		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", i		
			#currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)

			#currNode.readFromFile2(dirName, i)

			#self.localNodes.append(currNode)
			

			self.mapAlgorithm.restoreNode(dirName, i)
			self.mapAlgorithm.saveState()
		

		self.drawMap()


	def restoreSeries(self, dirName, num_poses):

		self.mapAlgorithm.restoreState(dirName, num_poses)

		PIXELSIZE = 0.05
		for i in range(0, num_poses):

			print "loading node", i			
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile2(dirName, i)
			self.localNodes.append(currNode)
			currNode.setGPACPose(self.mapAlgorithm.nodePoses[i])
			self.mapAlgorithm.nodeHash[i] = currNode

		foreNode = LocalNode(self.probe, self.contacts, num_poses, 19, PIXELSIZE)
		foreNode.readFromFile2(dirName, num_poses)
		backNode = LocalNode(self.probe, self.contacts, num_poses+1, 19, PIXELSIZE)
		backNode.readFromFile2(dirName, num_poses+1)
		self.mapAlgorithm.insertPose(foreNode, backNode, initLocation = foreNode.getEstPose())
		
		
		for i in range(num_poses+2, num_poses+70):

			print "loading node", i		
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile2(dirName, i)
			self.localNodes.append(currNode)
			self.mapAlgorithm.loadNewNode(currNode)
			self.mapAlgorithm.saveState()
		

		self.drawMap()

		return


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
		pass

		#self.occMap.update()
		#self.occMap.saveMap()

	def loadWalls(self, walls):
		self.walls = walls
		self.mapAlgorithm.walls = walls

	def getWalls(self):
		return self.walls

	@property
	def numNodes(self):
		return self.mapAlgorithm.numNodes

	@property
	def nodeHash(self):
		return self.mapAlgorithm.nodeHash

	def newNode(self, stepDist, direction):

		if self.numNodes > 0:
			
			if self.currNode.isDirty():
				self.currNode.synch()

			self.currNode.saveToFile()

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = True)		
		
		self.mapAlgorithm.addNode(self.currNode)

		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)

	def initNode(self, estPose, direction):

		self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19, self.pixelSize, stabilizePose = self.isStable, faceDir = True)

		self.mapAlgorithm.addInitNode(self.currNode, estPose)
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)	

	def correctPosture(self):
		self.currNode.updateCorrectPosture()

	def saveLocalMap(self):
		" save the local map "
		if self.currNode != 0:
			self.currNode.saveMap()
			

