import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import Image
import pylab
import scipy.optimize
from numpy import *
from subprocess import *
from random import *

class LocalOccMap:

	#def __init__(self, probe, contacts, nodeID):
	def __init__(self, localNode):

		self.localNode = localNode

		self.prevRects = []
		
		self.pixelSize = PIXELSIZE
		self.mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		self.prevJoints = [0.0,0.0,0.0,0.0]
		self.rootNeighbors = [18,19,20,21]

		self.nodeID = self.localNode.getNodeID()
		
		self.probe = self.localNode.getProbe()
		self.contacts = self.localNode.getContacts()
		self.rootNode = self.localNode.getRootNode()
		#self.rootPose = self.localNode.getRootPose()
		self.numJoints = self.probe.numSegs-1
		self.changed = False

		self.gnd_rectangles = []
		#self.gndInitRoot = self.probe.getActualJointPose(self.rootNode)
		self.gndMapImage = Image.new('L', (self.numPixel,self.numPixel),127)
		self.gndImage = self.gndMapImage.load()
		#da = self.gndInitRoot[2]	
		#self.gndT = array([[cos(da), -sin(da), 0.0],[sin(da), cos(da), 0.0],[0.0,0.0,1.0]])

		self.rectangles = []
		self.oldRectangles = []
		self.gndFileName = "gndLocalOccMap%03u" % self.nodeID + "_%04u.png"
		self.fileName = "localOccMap%03u" % self.nodeID + "_%04u.png"
		self.saveCount = 0
		self.mapImage = 0

		#self.costFile = open("costFile%03u.txt" % self.nodeID, 'w')
		self.hullFile = "convexHull%03u_%04u.png"
		self.hullCount = 0
		
		self.T = array([[1,0,0],
			[0,1,0],
			[0,0,1]], float)

	def saveToFile(self, filename):
		
		#1. first line has the offset
		#2. next is the printout of the rectangles
		val = repr(self.oldRectangles + self.rectangles)
		f = open(filename, 'w')
		f.write(val)
		f.close()

	def readFromFile(self):
		
		#self.fileName = "gndLocalOccMap%03u" % self.nodeID + "_%04u.png"
		self.fileName = "localOccMap%03u" % self.nodeID + "_%04u.png"
		
		self.mapImage = Image.open("localOccMap%03u" % self.nodeID + ".png")
		self.image = self.mapImage.load()
		
		self.changed = True
	
		#self.pixelSize = PIXELSIZE
		#self.mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0
		#self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		#self.halfPix = self.pixelSize/2.0
		#self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		#f2 = open(fileName, 'r')
		#val = f2.read()
		#localRects = eval(val)
		#f2.close()
		
		#self.setRootPose([0.0,0.0,0.0], 19)
		#self.setRectangles(localRects)
		#self.buildMap(True)
		#self.saveMap()
			
		#print "read", len(localRects), "rectangles"
		
	def getRectangles(self):
		return deepcopy(self.oldRectangles + self.rectangles)
	
	def setRectangles(self, rectangles):
		self.rectangles = rectangles
		self.oldRectangles = []
		
	"""
	def setRootPose(self, pose, node = 19):
		self.rootPose = pose
		self.rootNode = node
	"""
	
	"""
	Correct poses by centering joint 19 to 0,0
	"""
	def correctPoseBySegment(self, localRects):

		" first, find the transform to put segment 19 in the aligned position "
		
		" desired positions for two points"
		p0 = self.baseRect[1][0]
		p1 = self.baseRect[1][1]
		
		" find the 19th segment rectangle "
		seg19 = []
		for rect in localRects:
			if rect[0] == 19:
				seg19 = rect
				break
		
		if seg19 == []:
			print "pose did not have a rectangle 19"
			raise
		
		" starting positions for two points "
		i0 = seg19[1][0]
		i1 = seg19[1][1]
		
		" create matrices for solving for the parameters "
		b = array([[p0[0]],[p0[1]],[p1[0]],[p1[1]]])
		A = array([[i0[0],-i0[1],1,0],
				[i0[1],i0[0],0,1],
				[i1[0],-i1[1],1,0],
				[i1[1],i1[0],0,1]])
		
		result = linalg.solve(A,b)
		
		A_val = result[0,0]
		B_val = result[1,0]
		dx = result[2,0]
		dy = result[3,0]
		
		" build transform matrix "
		offsetT = array([[A_val, -B_val, dx],[B_val, A_val, dy],[0.0,0.0,1.0]])

		# make the offset version first
		newRects = []
		for seg, rect, isActive in localRects:
			newRect = []
			for p in rect:
				pVec = array([[p[0]],[p[1]],[1.0]])
				result = dot(offsetT,pVec)
				newRect.append([result[0,0], result[1,0]])

			newRects.append([seg,newRect, isActive])	

		return newRects
	
	def computeAlpha2(self, points):
		
		numPoints = len(points)
		inputStr = str(numPoints) + " "

		" radius 0.2 "
		inputStr += str(0.2) + " "
		
		for p in points:
			p2 = copy(p)
			p2[0] += gauss(0.0,0.01)
			p2[1] += gauss(0.0,0.01)

			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		
		inputStr += "\n"
		
		" start the subprocess "
		subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
		
		" send input and receive output "
		sout, serr = subProc.communicate(inputStr)

		#print numPoints
		#print sout
		
		" convert string output to typed data "
		sArr = sout.split(" ")
		
		#print sArr[0]
		numVert = int(sArr[0])
		
		vertices = []
		for i in range(numVert+1):
			vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
		
		#print vertices
		return vertices
	
	def alphaShapeCheck(self, newRects):
		
		" 1. pick out the points "
		points = []
		for j in range(self.numPixel):
			for k in range(self.numPixel):
				if self.image[j,k] == 255:
					pnt = self.gridToReal([j,k])
					points.append(pnt)
		
		" 2. create alpha shape "
		alphaShape = self.computeAlpha2(points)

		" print the new snapshot "
		pylab.clf()
		outside = []
		for seg, rect, isActive in newRects:
			xP = []
			yP = []
			for p in rect:

				pVec = array([[p[0]],[p[1]],[1.0]])
				result = dot(self.T,pVec)
				newP = [result[0,0],result[1,0]]
				
				xP.append(newP[0])
				yP.append(newP[1])
				
				val = point_inside_polygon(newP[0],newP[1],alphaShape[:-1])

				if val == False:
				#if val == 0:
					outside.append(newP)

			pVec = array([[rect[0][0]],[rect[0][1]],[1.0]])
			result = dot(self.T,pVec)
			newP = [result[0,0],result[1,0]]
				
			xP.append(newP[0])
			yP.append(newP[1])
				
			pylab.plot(xP,yP, color = 'b')


		initGuess = [0.0,0.0,0.0]
		offset = scipy.optimize.fmin(self.alphaCost,initGuess,args=(alphaShape, newRects))

		dx = offset[0]
		dy = offset[1]
		da = offset[2]	
		offsetT = array([[cos(da), -sin(da), dx],[sin(da), cos(da), dy],[0.0,0.0,1.0]])
		currT = dot(self.T,offsetT)
		
		# normalize the rotation component
		mag1 = sqrt(currT[0,0]**2 + currT[0,1]**2)
		currT[0,0] /= mag1
		currT[0,1] /= mag1
			
		#mag2 = sqrt(currT[1,0]**2 + currT[1,1]**2)
		currT[1,0] = -currT[0,1]
		currT[1,1] = currT[0,0]
			
		self.T = currT
				

		xP = []
		yP = []
		for p in points:
			xP.append(p[0])
			yP.append(p[1])
		pylab.scatter(xP,yP, faceted=False, linewidth = 0, color = '0.5')

	
		

		for seg, rect, isActive in newRects:
			xP = []
			yP = []
			for p in rect:

				pVec = array([[p[0]],[p[1]],[1.0]])
				result = dot(self.T,pVec)
				newP = [result[0,0],result[1,0]]
				
				xP.append(newP[0])
				yP.append(newP[1])

			pVec = array([[rect[0][0]],[rect[0][1]],[1.0]])
			result = dot(self.T,pVec)
			newP = [result[0,0],result[1,0]]
				
			xP.append(newP[0])
			yP.append(newP[1])
			pylab.plot(xP,yP, color = 'g')


		" mark the outside points of new snapshot "
		xP = []
		yP = []
		for p in outside:
			xP.append(p[0])
			yP.append(p[1])
	
		if len(xP) > 0:
			pylab.scatter(xP,yP,faceted=False, linewidth=0, color='r')

		#pylab.clf()
		xP = []
		yP = []
		print alphaShape[0], alphaShape[-1]
		for p in alphaShape:
			xP.append(p[0])
			yP.append(p[1])
		#xP.append(alphaShape[0][0])
		#yP.append(alphaShape[0][1])
		pylab.plot(xP,yP, color = 'r')


		pylab.xlim(-4,4)
		pylab.ylim(-4,4)
		
		pylab.savefig(self.hullFile % (self.nodeID, self.hullCount))
		self.hullCount += 1
		
		pylab.clf()

		"""
		1. compute the hull fitting cost
		2. twiddle the offset to find optimal fit
		3. adjust root pose
		
		"""



	def alphaCost(self, offset, alphaShape, newRects):
		
		dx = offset[0]
		dy = offset[1]
		da = offset[2]
		offsetT = array([[cos(da), -sin(da), dx],[sin(da), cos(da), dy],[0.0,0.0,1.0]])
		currT = dot(self.T, offsetT)

		sumDist = 0.0
		
		for seg, rect, isActive in newRects:		
			for p in rect:
				pVec = array([[p[0]],[p[1]],[1.0]])
				result = dot(currT,pVec)

				newP = [result[0,0],result[1,0]]

				val = point_inside_polygon(newP[0],newP[1],alphaShape[:-1])
				if not val:
					minDist, closePoint = self.findClosestPoint(newP, alphaShape[:-1])
					sumDist += minDist

		return sumDist

				


	
	def findClosestPoint(self, point, cHull):
		"""
		TWO CANDIDATES
		# 1. check each of the path points
		# 2. check orthogonal distance to each edge if less than both terminating points
		"""
		
		"  point-to-point distances "
		pointDist = []
		for p in cHull:
			dist = sqrt((point[0]-p[0])**2 + (point[1]-p[1])**2)
			pointDist.append(dist)

		" orthogonal distance to each segment "
		orthoDist = []
		closePoints = []

		for i in range(len(cHull)-1):
			p0 = copy(cHull[i])
			p1 = copy(cHull[i+1])

			# rotate each edge to 0 degrees, and compare point distance to the horizontal
			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			mag = sqrt(vec[0]**2 + vec[1]**2)
			if mag == 0:
				orthoDist.append( 1e100 )
				closePoints.append(copy(cHull[i]))
				print "degenerate"
				continue

			eAngle = acos(vec[0]/mag)
			if asin(vec[1]/mag) < 0:
				eAngle = -eAngle

			# rotate the edge
			rotVec = [vec[0] * cos(eAngle) + vec[1] * sin(eAngle), -vec[0] * sin(eAngle) + vec[1] * cos(eAngle)]

			oPoint = [point[0]-p0[0], point[1]-p0[1]]
			oPoint = [oPoint[0] * cos(eAngle) + oPoint[1] * sin(eAngle), -oPoint[0] * sin(eAngle) + oPoint[1] * cos(eAngle)]

			dist = abs(oPoint[1])
			if oPoint[0] < 0.0:
				orthoDist.append( pointDist[i])
				closePoints.append(copy(cHull[i]))

			elif oPoint[0] > mag:
				orthoDist.append(pointDist[i+1])
				closePoints.append(copy(cHull[i+1]))

			else:
				orthoDist.append(dist)
				nearPoint = [oPoint[0],0.0]
				nearPoint = [nearPoint[0] * cos(eAngle) - nearPoint[1] * sin(eAngle), nearPoint[0] * sin(eAngle) + nearPoint[1] * cos(eAngle)]
				closePoints.append(copy(nearPoint))

		minDist = 1e100
		minIndex = 0
		for i in range(len(orthoDist)):
			if orthoDist[i] < minDist:
				minDist = orthoDist[i]
				minIndex = i
				
		return minDist, copy(closePoints[minIndex])
	
	def costFunction(self, offset, rects1, rects2):
		cost = 0.0
		compareCount = 0

		dx = offset[0]
		dy = offset[1]
		da = offset[2]

		offsetT = array([[cos(da), -sin(da), dx],[sin(da), cos(da), dy],[0.0,0.0,1.0]])
		currT = dot(self.T, offsetT)
		
		# make the offset version first
		newRects = []
		for seg, rect, isActive in rects2:
			newRect = []
			for p in rect:
				pVec = array([[p[0]],[p[1]],[1.0]])
				result = dot(currT,pVec)
				newRect.append([result[0,0], result[1,0]])

			newRects.append([seg,newRect, isActive])	

		i1 = rects1.__iter__()
		i2 = newRects.__iter__()	
		
		try:
			r1 = i1.next()
			r2 = i2.next()
			
			while True:
				# only compare rectangles that represent the same segment
				if r1[0] >= 18 and r1[0] < 22 and r1[0] == r2[0] and r1[2] >= 0.9 and r2[2] >= 0.9:
					localCost = self.rectCost(r1[1],r2[1])
					#print "comparing segment", r1[0], "with cost =", localCost

					cost += localCost
					compareCount += 1
					r1 = i1.next()
					r2 = i2.next()
					
				elif r1[0] < r2[0]:
					r1 = i1.next()
					
				elif r1[0] > r2[0]:
					r2 = i2.next()
									
				else:
					r1 = i1.next()
					r2 = i2.next()
		
		except:
			pass
		
		if compareCount > 0:
			return cost/compareCount
		else:
			return 0.0
	
	def rectCost(self, rect1, rect2):
		# ideally, we should be measuring the intersecting area
		# here, we will compute the distance between the end points

		dist = 0.0
		for i in range(4):
			p1 = rect1[i]
			p2 = rect2[i]
			
			dist += sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
		
		cost = 1.0 - 100.0**-dist
		return cost

	def update(self, isForward = True):
		
		"""
		1. add a ground truth map to the LocalOccMap class
		2. requires setting of initial pose
		2. requires setting of ground pose at each iteration
		3. requires plotting of each ground rectangle to map, but skipping all corrections
		
		"""
		
		# compute ground free space first
		gndRects = self.gndComputeFreeSpace()
		gndRects.sort()
		self.gnd_rectangles.append(gndRects)
		
		# compute free space rectangles
		freeRects = self.computeFreeSpace(isForward)
		freeRects.sort()

		jointErr = 0.0
		currJoints = []
		for j in self.rootNeighbors:
			currJoints.append(self.probe.getServo(j))

		for i in range(len(currJoints)):
			jointErr += abs(currJoints[i]-self.prevJoints[i])
		
		cost = self.costFunction([0.0,0.0,0.0], self.prevRects, freeRects)
		
		" FIXME:  turn off local correction "
		cost = 0.01
		
		if cost > 0.05:
			
			" make sure all previous rectangles are written to map"
			self.buildMap()

			self.alphaShapeCheck(freeRects)

		newRects = []
		for seg, rect, isActive in freeRects:	
			newRect = []	
			for p in rect:
				pVec = array([[p[0]],[p[1]],[1.0]])
				result = dot(self.T,pVec)
				newP = [result[0,0],result[1,0]]
				newRect.append(newP)
			newRects.append([seg,newRect, isActive])

		self.prevRects = deepcopy(newRects)		
		self.prevJoints = currJoints


		self.rectangles.append(newRects)

	def gndComputeFreeSpace(self):
							
		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		actualConfig = []
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1
		
		currGnd = self.probe.getActualJointPose(self.rootNode)

		"""
		dx = currGnd[0] - self.gndInitRoot[0]
		dy = currGnd[1] - self.gndInitRoot[1]
		
		pVec = array([[dx],[dy],[1.0]])
		result = dot(self.gndT,pVec)
		newOrigin = [result[0,0],result[1,0], currGnd[2]-self.gndInitRoot[2]]
		
		xTotal = newOrigin[0]
		zTotal = newOrigin[1]
		totalAngle = newOrigin[2]
		"""
		
		globalVec = array([[currGnd[0]], [currGnd[1]]])
		tempVec = dot(self.localNode.gndBackR, globalVec)
		tempVec[0, 0] -= self.localNode.gndDist
		transVec = dot(self.localNode.gndForeR, tempVec)
		finalVec = dot(self.localNode.gndR, transVec)
		
		xTotal = finalVec[0, 0]
		zTotal = finalVec[1, 0]
		totalAngle = normalizeAngle(currGnd[2] - self.localNode.gndPose[2])

		#print xTotal, zTotal, totalAngle
		
		joints = range(-1,self.rootNode)
		joints.reverse()

		for i in joints:

			totalAngle = totalAngle + self.probe.getServo(i+1)
			totalAngle = normalizeAngle(totalAngle)

			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])


		joints = range(self.rootNode, self.probe.numSegs-1)

		#xTotal = newOrigin[0]
		#zTotal = newOrigin[1]
		#totalAngle = newOrigin[2]
		xTotal = finalVec[0, 0]
		zTotal = finalVec[1, 0]
		totalAngle = normalizeAngle(currGnd[2] - self.localNode.gndPose[2])

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])
			
			if i < self.probe.numSegs-2:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)

		actualConfig.sort()
		return actualConfig	
	
	def computeFreeSpace(self, isForward = True):
							
		# 1. for every joint reference, compute the distal desired joint configuration and the actual joint configuration
		# 2. for desired joint configuration, set occupancy to obstacle if not already free space
		# 3. for actual joint configuration, set to free space, even if previously set as obstacle

		actualConfig = []
		
		segLength = self.probe.segLength
		segWidth = self.probe.segWidth
		numJoints = self.probe.numSegs-1

		#xTotal = self.rootPose[0]
		#zTotal = self.rootPose[1]
		#totalAngle = self.rootPose[2]

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		joints = range(-1,self.rootNode)
		joints.reverse()

		for i in joints:

			totalAngle = totalAngle + self.probe.getServo(i+1)
			totalAngle = normalizeAngle(totalAngle)

			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])

		joints = range(self.rootNode, self.probe.numSegs-1)

		#xTotal = self.rootPose[0]
		#zTotal = self.rootPose[1]
		#totalAngle = self.rootPose[2]

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			actualConfig.append([i, [p4,p3,p2,p1], self.contacts.initMask[i]])
			
			if i < self.probe.numSegs-2:
				totalAngle = totalAngle - self.probe.getServo(i+1)
				totalAngle = normalizeAngle(totalAngle)

		actualConfig.sort()
		return actualConfig		
	
	def fillOpen(self, polygon, image):

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
				#if IsContained(newRect, point):
					# free space overwrites everything
					image[i,j] = 255

	def realToGrid(self, point):
		indexX = int(floor(point[0]*self.divPix)) + self.numPixel/2 + 1
		indexY = int(floor(point[1]*self.divPix)) + self.numPixel/2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0, (j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]
		return point

	def buildMap(self, redraw = False):
	
		self.pixelSize = PIXELSIZE
		self.mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		# 1. create the map image
		if self.mapImage == 0 or redraw:
			self.mapImage = Image.new('L', (self.numPixel,self.numPixel),127)
			self.image = self.mapImage.load()

		if redraw:
			# 2. plot all rectangles to the image
			for snapshot in self.oldRectangles + self.rectangles:
				for rect in snapshot:
					self.fillOpen(rect[1], self.image)

		else:
			# 2. plot just the latest rectangles to the image
			for snapshot in self.rectangles:
				for rect in snapshot:
					self.fillOpen(rect[1], self.image)
			
			for snapshot in self.gnd_rectangles:
				#print snapshot
				for rect in snapshot:
					self.fillOpen(rect[1], self.gndMapImage.load())
		
		" so outside maps can trigger changes "
		self.changed = True

		" save the rectangles, and clear the list for new ones "
		#self.oldRectangles += self.rectangles
		self.rectangles = []
		self.gnd_rectangles = []
	
	def saveMap(self, fileName = ""):
		
		self.buildMap()
		
		if fileName != "":
			self.mapImage.save(fileName)
		else:
			self.mapImage.save(self.fileName % self.saveCount)	
			self.gndMapImage.save(self.gndFileName % self.saveCount)
			self.saveCount += 1
		
	def getMap(self):
		return self.mapImage
	
	def getGndMap(self):
		return self.gndMapImage
