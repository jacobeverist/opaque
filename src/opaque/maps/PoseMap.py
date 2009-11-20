import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import graph
import csv
import Image
from LocalOccMap import *

class PoseMap:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts
		
		# pose sizes
		self.pixelSize = PIXELSIZE
		self.mapSize = 20.0
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		self.fileName = "poseMap%04u.png"
		self.saveCount = 0
		
		self.initPose = self.probe.getActualJointPose(19)
		#self.initPose = self.contacts.getClosestPose(19)
		self.poseMaps = []
		
		self.newPose()

		self.poseGraph = graph.graph()

	def newPose(self):
		newMap = LocalOccMap(self.probe, self.contacts, len(self.poseMaps))
		newMap.setRootPose([0.0,0.0,0.0], 19)
		#newMap.fileName = "localOccMap%02u" % len(self.poseMaps) + "_%04u.png"
		if len(self.poseMaps) > 0:
			
			# convert the last one to a file
			
			#self.saveMap()
			poseProfile = self.contacts.captureRef()
			self.poseMaps.append([newMap,self.contacts.getClosestPose(19),poseProfile])
			
		else:
			self.poseMaps.append([newMap,self.initPose, 0])


		print "new pose at", self.poseMaps[-1][1]

	def update(self, isForward = True):
		if isForward:
			newRoot = self.contacts.getClosestPose(34)
			newRoot2 = self.probe.getActualJointPose(34)
		else:
			newRoot = self.contacts.getClosestPose(5)
			newRoot2 = self.probe.getActualJointPose(5)

		origin = self.poseMaps[-1][1]
		relRoot = [newRoot[0]-origin[0], newRoot[1]-origin[1], newRoot[2]]
		relRoot2 = [newRoot2[0]-origin[0], newRoot2[1]-origin[1], newRoot2[2]]

		#relRoot[0] *= -1
		#relRoot[1] *= -1
		
		#relRoot[2] += pi
		#relRoot[2] = normalizeAngle(relRoot[2])
		#relRoot[2] = newRoot[2]

		if isForward:
			#print "actual:", origin, newRoot2, relRoot2
			#print "global:", origin, newRoot
			#print "local:", [0.0,0.0,0.0], relRoot
			self.poseMaps[-1][0].setRootPose(relRoot, 34)
		else:
			#print "actual:", origin, newRoot2, relRoot2
			#print "global:", origin, newRoot
			#print "local:", [0.0,0.0,0.0], relRoot
			self.poseMaps[-1][0].setRootPose(relRoot, 5)
			
		self.poseMaps[-1][0].update(isForward)
	
	def saveToFile(self):

		f = open("poses.txt", 'w')

		for localMap, pose, poseProfile in self.poseMaps:
			localFile = "rect%03u.txt" % localMap.nodeID
			localMap.saveToFile(localFile)

			f.write(str(localMap.nodeID) + " " + str(pose[0]) + " ")
			f.write(str(pose[1]) + " " + str(pose[2]) + " " + localFile + "\n")
		
		
		f.close()
		
	def readFromFile(self, fileName):
		
		fileReader = csv.reader(file(fileName), delimiter=' ')
		self.poseMaps = []

		for row in fileReader:
			nodeID = int(row[0])
			x = float(row[1])
			y = float(row[2])
			t = float(row[3])
			localFile = row[4]

			f2 = open(localFile, 'r')
			localRects = eval(f2.read())
			#print localRects
			#print nodeID, x, y, t, localFile, len(localRects)
			
			newMap = LocalOccMap(self.probe, nodeID)
			newMap.setRootPose([0.0,0.0,0.0], 19)
			newMap.setRectangles(localRects)
			self.poseMaps.append([newMap,[x,y,t]])

	def realToGrid(self, point):
		indexX = int(floor(point[0]*self.divPix)) + self.numPixel/2 + 1
		indexY = int(floor(point[1]*self.divPix)) + self.numPixel/2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0, (j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]
		return point
						
	def saveMap(self):
		self.poseMaps[-1][0].saveMap()
		
		return
			
		#1. for each pose
		# 2. offset rectangles by x/y offset
		# 3. plot those rectangles to map

		self.mapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.image = self.mapImage.load()
		
		allRects = []
		
		for localMap, pose, poseProfile in self.poseMaps:
			rects = localMap.getRectangles()
			
			for snapshot in rects:
				for rect in snapshot:
					for p in rect:
						p[0] += pose[0]
						p[1] += pose[1]
						
			allRects += rects

		for snapshot in allRects:
			for rect in snapshot:
				self.fillOpen(rect)

		self.mapImage.save(self.fileName % self.saveCount)	
		self.saveCount += 1
				
		#print self.poseMaps

	def fillOpen(self, polygon):

		# vertices of the 4-sided convex polygon
		polyX = [polygon[0][0], polygon[1][0], polygon[2][0], polygon[3][0], polygon[0][0]]
		polyY = [polygon[0][1], polygon[1][1], polygon[2][1], polygon[3][1], polygon[0][1]]

		# bounding box of polygon
		maxX = -10000.0
		minX = 10000.0
		maxY = -10000.0
		minY = 10000.0
		for i in range(4):
			if polyX[i] > maxX:
				maxX = polyX[i]
			if polyX[i] < minX:
				minX = polyX[i]
			if polyY[i] > maxY:
				maxY = polyY[i]
			if polyY[i] < minY:
				minY = polyY[i]

		# set boundaries of map grid
		highIndexX, highIndexY = self.realToGrid([maxX, maxY])
		lowIndexX, lowIndexY = self.realToGrid([minX, minY])

		# fill map grid cell if occupied by arm
		indices = []

		for i in range(lowIndexX,highIndexX+1,1):
			for j in range(lowIndexY,highIndexY+1,1):
				point = self.gridToReal([i,j])
				if IsContained(polygon, point):
					# free space overwrites everything
					if i >= 0 and i < self.numPixel and j >= 0 and j < self.numPixel:
						self.image[i,j] = 255


	
	