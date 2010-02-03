import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *
import traceback

fastLocalNodeCount = 0

class FastLocalCurveFit(Behavior):

	#def __init__(self, probe, startNode = 0, endNode = NUM_SEGS, curve = 0):
	def __init__(self, probe, anterior = True, curve = 0, cutoff = 19):
		global fastLocalNodeCount

		Behavior.__init__(self, probe)

		self.numSegs = self.probe.numSegs
		
		self.startNode = 0
		self.endNode = 39		
		
		self.curve = curve

		" anterior or posterior "
		self.anterior = anterior

		print "setting rootNode to ", cutoff
		self.rootNode = cutoff

		" whether peaks are solid or not "
		self.solid = [False for i in range(30)]

		" the commanded angles for solid joints, None if joint is not solid "
		self.solidJoints = [None for i in range(0,39)]
		
		" the actual angles for solid joints, None if joint is not solid "
		self.solidCompJoints = [None for i in range(0,39)]
		
		" joints classified according to their solid index, does not change"		
		self.solidJointClasses = [[] for i in range(30)]
		
		" location of joint origin of commanded in curve space, does not change "
		self.solidJointPositions = [None for i in range(0,39)]
		
		" classification of joints into peak indices, updated every step "
		self.jointClasses = [[] for i in range(30)]
		self.jointErrorClasses = [[] for i in range(30)]
		
		" last round setting of joint angles and origins "
		self.jointPositions = [None for i in range(0,39)]
		self.jointAngles = [None for i in range(0,39)]

		self.plotCount = 0
		
		self.spliceJoint = 0
		self.spliceAngle = 0.0
		
		"""		
		# There are 2 possibilities to how the local fit is computed
		# if we assume that the curve is always on the extrema, and never the interior of the snake
		# Front half:
		#    - rooted at tip
		# Back Half
		#    - rooted at tip
		"""
	
		# nodes for drawing purposes
		self.nodeCount = fastLocalNodeCount
		fastLocalNodeCount += 1
		self.childNodes = []
		self.childEntities = []
		self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("fastLocalCurveRoot" + str(self.nodeCount))
	
	def setStartNode(self, node):
		self.startNode = node

	def setEndNode(self, node):
		self.endNode = node
	
	def setSide(self, val):
		self.anterior = val

	def setRootNode(self, node):
		self.rootNode = node

	def draw(self):
		
		# remove all children
		self.parentNode.removeAllChildren()

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []
	
		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)
	
		self.childEntities = []

		points = self.curve.getPoints()

		# draw the points in the simulation
		for i in range(len(points)):

			pnt = points[i]

			childNode = self.parentNode.createChildSceneNode("fastLocalCurvePoint" + str(i) + "_" + str(self.nodeCount))
			self.childNodes.append(childNode)

			currEntity = self.probe._mgr.createEntity("fastLocalCurveEnt"+str(i) + "_" + str(self.nodeCount), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)

			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			#position = ogre.Vector3(0.0,0.0,0.0)
			childNode.setPosition(position)

			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)

			childNode.attachObject(currEntity)

	def setCurve(self, curve):
		self.curve = curve

	def setBoundaries(self, startNode, endNode):
		self.startNode = startNode
		self.endNode = endNode

	def getPeakError(self, index):
					
		" last round setting of joint angles and origins "
		#self.jointPositions = [None for i in range(0,39)]
		#self.jointAngles = [None for i in range(0,39)]
		
		cmdPoints = []
		actPoints = []
		
		joints = self.getPeakJoints(index)
		joints.sort()

		if len(joints) == 0:
			return 0.0
		
		minIndex = joints[0]
		
		if minIndex > 0:
			origin = self.jointPositions[minIndex-1]
		else:
			minIndex = 1
			origin = self.jointPositions[minIndex-1]
			
		cmdPoints.append(copy(origin))
		actPoints.append(copy(origin))
			
		for i in joints:
			sampAngle = self.probe.getServo(i)
			totalAngle = origin[2] - sampAngle
			xTotal = origin[0] + self.probe.segLength*cos(totalAngle)
			zTotal = origin[1] + self.probe.segLength*sin(totalAngle)
			pnt = [xTotal, zTotal, totalAngle]
			origin = pnt

			cmdPoints.append(copy(self.jointPositions[i]))
			actPoints.append(copy(origin))
			maxIndex = i

		" add one more joint over the peak "
		sampAngle = self.probe.getServo(maxIndex+1)
		totalAngle = origin[2] - sampAngle
		xTotal = origin[0] + self.probe.segLength*cos(totalAngle)
		zTotal = origin[1] + self.probe.segLength*sin(totalAngle)
		pnt = [xTotal, zTotal, totalAngle]
		origin = pnt

		cmdPoints.append(copy(self.jointPositions[i]))
		actPoints.append(copy(origin))
	
		sum = 0.0
		for i in range(len(cmdPoints)):
			
			pC = cmdPoints[i]
			pA = actPoints[i]
			
			dist = sqrt((pC[0]-pA[0])**2 + (pC[1]-pA[1])**2)
	
			sum += dist
		
		cost = sum / len(cmdPoints)
		
		return cost
		

	def getPeakJoints(self, index):
		return self.jointClasses[index]

	def getPeakErrorJoints(self, index):
		return self.jointErrorClasses[index]

	def isJointSolid(self, joint):
		peakIndex = -1
		
		" We haven't initialized the solid joint classes yet, return false "
		if len(self.solidJointClasses[0]) == 0:
			#print "jointclass 0 is empty"
			return False
		
		for i in range(len(self.solidJointClasses)):
			if self.solidJointClasses[i].count(joint) > 0:
				peakIndex = i
				break
				
		if peakIndex == -1:
			return False
		else:
			
			if self.solid[peakIndex+1]:
				return True
			else:
				if self.anterior:
					maxJoint = max(self.solidJointClasses[peakIndex])
					if maxJoint != joint:
						return True
					else:
						return False
				else:
					minJoint = min(self.solidJointClasses[peakIndex])
					if minJoint != joint:
						return True
					else:
						return False					
			
		
		#return self.solid[peakIndex]
	
	def getSolidIndex(self, joint):
		peakIndex = -1
		
		" We haven't initialized the solid joint classes yet, return false "
		if len(self.solidJointClasses[0]) == 0:
			#print "jointclass 0 is empty"
			return False
		
		for i in range(len(self.solidJointClasses)):
			if self.solidJointClasses[i].count(joint) > 0:
				peakIndex = i
				
		if peakIndex == -1:
			raise
		else:
			return peakIndex	
		
	def getSolidJoint(self, joint):
		
		peakIndex = -1
		for i in range(len(self.solidJointClasses)):
			if self.solidJointClasses[i].count(joint) > 0:
				peakIndex = i
				
		if peakIndex == -1:
			raise
		
		if not self.solid[peakIndex]:
			raise
		
		return self.solidJoints[joint]
	
	def getCenterPoints(self):
		
		points = [[self.curve.ampSpacing/2.0 + i * self.curve.ampSpacing, -self.curve.peakAmp[0]] for i in range(self.curve.numPeaks)]

		" 1. at each inflection point, find the closest segment and it's point on segment"
		" 2. cycle through all 40 segments "
		" 3. transform those points to current estimated positions "
		
		centerPoints = []
		
		for p in points:
			minDist = 1e100
			minPoint = [0.0,0.0]
			jointIndex = 0
			segLen = 0.0
			
			for i in range(len(self.solidJointPositions)-1):
				#for j in [[-ARMLENGTH,0.0]] + self.solidJointPositions:
				p0 = self.solidJointPositions[i]
				p1 = self.solidJointPositions[i+1]
				seg = [p0,p1]
				dist, cp = closestSegPoint(seg, p)

				if dist < minDist:
					minDist = dist
					minPoint = cp
					jointIndex = i
					segLen = sqrt((p0[0]-minPoint[0])**2 + (p0[1]-minPoint[1])**2)

			" return joint number and position on segment "
			centerPoints.append([jointIndex,segLen])
		
		return centerPoints

	def getSolidJointPosition(self, joint):
		return self.solidJointPositions[joint]

	def setSolid(self, index):

		joints = self.getPeakJoints(index)

		if len(joints) == 0:
			return
		
		self.solid[index] = True
		
		#print "setting solid", index, "with joints", joints
		for i in joints:
			self.solidJoints[i] = pi / 180.0 * self.getJoints()[i]
			#self.solidJoints[i] = self.probe.getServoCmd(i)
			
			self.solidCompJoints[i] = self.probe.getServo(i)
			
			if i == self.spliceJoint:
				self.solidCompJoints[i] -= self.spliceAngle
			
			print "setting solid joint", i, self.probe.getServo(i), self.solidCompJoints[i]
			
			self.solidJointClasses[index].append(i)
			self.solidJointPositions[i] = copy(self.jointPositions[i])
			#print self.solidJoints[i], self.solidCompJoints[i]
			#actAngle = self.probe.getServo(i)
			#self.setJoint(i,actAngle*180.0/pi)


		" don't set the previous joint solid if it's not part of this curve "
		if index > 0:
			if self.anterior:
				
				" now set the maximum joint of the previous peak "
				minJoint = min(joints) - 1
				if minJoint >= 0:
					print "altering joint", minJoint, "at index", index
					self.solidJoints[minJoint] = pi / 180.0 * self.getJoints()[minJoint]
					#self.solidJoints[minJoint] = self.probe.getServoCmd(minJoint)
					

					self.solidCompJoints[minJoint] = self.probe.getServo(minJoint)
					if minJoint == self.spliceJoint:
						self.solidCompJoints[minJoint] -= self.spliceAngle

					print "setting solid joint", minJoint, self.probe.getServo(minJoint), self.solidCompJoints[minJoint]
					
					self.solidJointPositions[minJoint] = copy(self.jointPositions[minJoint])
				
			else:
				" now set the maximum joint of the previous peak "
				maxJoint = max(joints) + 1
				if maxJoint <= self.probe.numSegs -2:
					print "altering joint", maxJoint, "at index", index
					self.solidJoints[maxJoint] = pi / 180.0 * self.getJoints()[maxJoint]
					#self.solidJoints[maxJoint] = self.probe.getServoCmd(maxJoint)
					

					self.solidCompJoints[maxJoint] = self.probe.getServo(maxJoint)
					if maxJoint == self.spliceJoint:
						self.solidCompJoints[maxJoint] -= self.spliceAngle

					print "setting solid joint", maxJoint, self.probe.getServo(maxJoint), self.solidCompJoints[maxJoint]

					self.solidJointPositions[maxJoint] = copy(self.jointPositions[maxJoint])
			

		#print self.solid
		
	def resetSolid(self, index):
		self.solid[index] = False
	
	def getSolid(self, index):
		return self.solid[index]

	def classifyJoints(self):

		self.jointClasses = [[] for i in range(30)]
		self.jointErrorClasses = [[] for i in range(30)]
		
		for i in range(len(self.jointPositions)):
		
			pos = self.jointPositions[i]
			
			" don't classify an unfitted joint "
			if pos == None:
				continue
		
			ampIndex = -1
			
			" determine classification index "
			#if pos[0] < 3.0*self.curve.ampSpacing/2.0:
			#	ampIndex = 0
			#else:
			#	ampIndex = int((pos[0] - self.curve.ampSpacing/2.0) / self.curve.ampSpacing)

			if pos[0] >= self.curve.getHeadLength():
				for j in range(len(self.curve.infPoints)):
					infX = self.curve.infPoints[j]
					if pos[0] < infX:
						ampIndex = j
						break
		
			" save this classification "
			#print "appending", i, "to ampIndex=", ampIndex
			if ampIndex >= 0:
				self.jointClasses[ampIndex].append(i)
			
				errIndex = int(pos[0]/self.curve.ampSpacing)	
				#print errIndex, i, pos
				self.jointErrorClasses[errIndex].append(i)
	
		#print self.curve.infPoints
		#print self.jointPositions
		#print self.jointClasses
		
	def forwardFit(self):
		
		# origin is at the tip of segment 0, so we need to set the joint center to -self.segLength
		originPose = [0.0, 0.0, 0.0]
		
		minJoint = 0
		maxJoint = self.rootNode
		
		if minJoint < self.startNode:
			minJoint = self.startNode
		if maxJoint > self.endNode:
			maxJoint = self.endNode

		print "forwardFit:", minJoint, maxJoint

		# indices	
		ind = range(minJoint,maxJoint+1)
		
		maxNode = 0
		maxU = 0.0

		" positions and angles of the joints computed for this step "
		self.jointPositions = [None for i in range(self.probe.numSegs-1)]
		self.jointAngles = [None for i in range(self.probe.numSegs-1)]
		
		radius = self.probe.segLength
		breakFromFit = False
		
		for currJoint in ind:
			
			crossOvers = []
			
			self.jointPositions[currJoint]= originPose
			
			if self.isJointSolid(currJoint):
				"this is a solid joint"
				solidIndex = self.getSolidIndex(currJoint)
			else:
				"not a solid joint"
				solidIndex = -1

			
			" 2 peaks in a row must be solid before we re-use old data "
			" this is to allow segments spanning two neighboring peaks to adjust freely "
			if solidIndex >= 0 and self.solid[solidIndex]:

				" set to actual angle"
				self.setJoint(currJoint,self.solidCompJoints[currJoint] * 180.0 / pi)
				
				" theoretical angle "
				sampAngle = self.solidJoints[currJoint]
				
				totalAngle = originPose[2] - sampAngle
				xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
				zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				originPose = pnt
	
				self.jointAngles[currJoint] = sampAngle
				
				#print "solid:", currJoint, self.solidCompJoints[currJoint]
				#print "solid joint", currJoint, "has index", solidIndex
				
			else:
				self.foo = 0
				
				samples = []
				localMin = 1e10
				localAng = 0
				
				for angle in arange(-150,150,1.0):	
					" looking for cross over events, distance passes threshold of radius "
					sampAngle = angle*pi/180.0

					totalAngle = originPose[2] - sampAngle
					xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
					zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]

					if pnt[0] >= originPose[0]:
						#min, u, curvePoint = self.curve.findClosestPoint(pnt, xNotLessThan=originPose[0])
						
						try:
							curvePoint = self.curve.computeX(pnt, xNotLessThan=originPose[0])	
						
							dist = sqrt((curvePoint[0] - originPose[0])**2 + (curvePoint[1] - originPose[1])**2)
							
							pntDist = sqrt((curvePoint[0] - pnt[0])**2 + (curvePoint[1] - pnt[1])**2)
							if pntDist < localMin:
								localMin = pntDist
								localAng = sampAngle
							
							samples.append([sampAngle, pnt, dist, curvePoint])

							if len(samples) > 1:
								" check if this is a true cross over "
								dist1 = samples[-1][2]
								dist2 = samples[-2][2]
								edge = [samples[-1][3],samples[-2][3]]
								pnt1 = samples[-1][1]
								pnt2 = samples[-2][1]
								avgPoint = [(pnt2[0]+pnt1[0])/2.0, (pnt2[1]+pnt1[1])/2.0]
								
								if dist1 > radius and dist2 < radius:					
									if LeftOnEdge(edge,pnt1) and not LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])
									elif not LeftOnEdge(edge,pnt1) and LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])

								elif dist1 < radius and dist2 > radius:
									if LeftOnEdge(edge,pnt1) and not LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])
									elif not LeftOnEdge(edge,pnt1) and LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])
								
						except:
							" x point is beyond curve range "
							" breaking assumption only works if tail is straight "
							#print "exception"
							#traceback.print_exc()
							breakFromFit = True
							break
							
				if breakFromFit:
					break

				if len(crossOvers) == 0:
					crossOvers.append([1.0*pi/180.0 + localAng, -1.0*pi/180.0 + localAng])
				
				" if no crossover event, take the minimal distance and search the neighborhood "

				samples2 = []		
				candidates = []		
				for events in crossOvers:

					angle1 = events[0] * 180.0 / pi
					angle2 = events[1] * 180.0 / pi
					if angle1 > angle2:
						tempAngle = angle1
						angle1 = angle2
						angle2 = tempAngle
						
					minDist = 1e100
					candidate = []
					
					for angle in arange(angle1,angle2,0.1):
						sampAngle = angle*pi/180.0
						totalAngle = originPose[2] - sampAngle
						xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)								
						pnt = [xTotal, zTotal, totalAngle]
						#print pnt[0], originPose[0]
						
						if pnt[0] > originPose[0]:
							min, u, curvePoint = self.curve.findClosestPoint(pnt, xNotLessThan=originPose[0])
							samples2.append([sampAngle, pnt, min, u, curvePoint])
							
							#print min, u, curvePoint, pnt
							#print minDist, min
							if minDist > min:
								minDist = min
								candidate = [sampAngle, pnt, min, u, curvePoint]
					
					candidates.append(candidate)

				# 4. find closest of all these points and select it for the servo position

				" select the point with lowest x "
				min_i = -1
				minX = 1e100
				for i in range(len(candidates)):
					pnt = candidates[i][1]

					if minX > pnt[0]:
						minX = pnt[0]
						min_i = i

				" new servo position and segment end set to new origin "
				newServo = candidates[min_i][0]  * 180.0/pi
				newOrigin = candidates[min_i][1]

				" determine if we've reached the end of the curve "
				if candidates[min_i][3] >= maxU:
					maxNode = currJoint
					maxU = candidates[min_i][3]

					if maxU == 1.0:
						maxNode = currJoint - 1 

				originPose = newOrigin

				# 5. kinematically update the frame of reference for the next joint and go to 1 for joint i+1 
				if maxU != 1.0:
					self.setJoint(currJoint, newServo)
					self.jointAngles[currJoint] = newServo * pi / 180.0
				else:
					break
	
		" position of last joint iterated through "
		self.jointPositions[currJoint]= originPose
		
		" now classify the joints to their assigned peaks "
		self.classifyJoints()

	" this code has been modified and the changes are untested "
	def backwardFit(self):
					
		# origin is at the tip of segment 39
		originPose = [0.0, 0.0, pi]

		# indices
		minJoint = 0
		maxJoint = self.rootNode
					
		if minJoint < self.startNode:
			minJoint = self.startNode
		if maxJoint > self.endNode:
			maxJoint = self.endNode
			
		print "backwardFit:", minJoint, maxJoint

		ind = range(minJoint,maxJoint+1)
		ind.reverse()

		minNode = self.probe.numSegs-2
		maxU = 0.0

		" positions and angles of the joints computed for this step "
		self.jointPositions = [None for i in range(self.probe.numSegs-1)]
		self.jointAngles = [None for i in range(self.probe.numSegs-1)]
		
		radius = self.probe.segLength
		breakFromFit = False

		for currJoint in ind:

			crossOvers = []

			self.jointPositions[currJoint]= originPose
			
			if self.isJointSolid(currJoint):
				"this is a solid joint"
				solidIndex = self.getSolidIndex(currJoint)
			else:
				"not a solid joint"
				solidIndex = -1
				
			
			" 2 peaks in a row must be solid before we re-use old data "
			" this is to allow segments spanning two neighoring peaks to adjust freely "
			if solidIndex >= 0 and self.solid[solidIndex]:
			#if False:
				" set to actual angle"
				self.setJoint(currJoint,self.solidCompJoints[currJoint] * 180.0 / pi)
				
				" theoretical angle "
				sampAngle = self.solidJoints[currJoint]
				
				totalAngle = originPose[2] + sampAngle
				xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
				zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
				pnt = [xTotal, zTotal, totalAngle]
				originPose = pnt
	
				self.jointAngles[currJoint] = sampAngle
	
			else:
				self.foo = 0
				
				samples = []
				localMin = 1e10
				localAng = 0
				
				for angle in arange(-150,150,1.0):	
					" looking for cross over events, distance passes threshold of radius "
					sampAngle = angle*pi/180.0

					totalAngle = originPose[2] + sampAngle
					xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
					zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
					pnt = [xTotal, zTotal, totalAngle]

					if pnt[0] >= originPose[0]:

						try:
							curvePoint = self.curve.computeX(pnt, xNotLessThan=originPose[0])	
						
							dist = sqrt((curvePoint[0] - originPose[0])**2 + (curvePoint[1] - originPose[1])**2)
							
							pntDist = sqrt((curvePoint[0] - pnt[0])**2 + (curvePoint[1] - pnt[1])**2)
							if pntDist < localMin:
								localMin = pntDist
								localAng = sampAngle
							
							samples.append([sampAngle, pnt, dist, curvePoint])
						
							if len(samples) > 1:
								" check if this is a true cross over "
								dist1 = samples[-1][2]
								dist2 = samples[-2][2]
								edge = [samples[-1][3],samples[-2][3]]
								pnt1 = samples[-1][1]
								pnt2 = samples[-2][1]
								avgPoint = [(pnt2[0]+pnt1[0])/2.0, (pnt2[1]+pnt1[1])/2.0]
								
								if dist1 > radius and dist2 < radius:					
									if LeftOnEdge(edge,pnt1) and not LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])
									elif not LeftOnEdge(edge,pnt1) and LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])

								elif dist1 < radius and dist2 > radius:
									if LeftOnEdge(edge,pnt1) and not LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])
									elif not LeftOnEdge(edge,pnt1) and LeftOnEdge(edge,pnt2):
										crossOvers.append([samples[-1][0],samples[-2][0],avgPoint])

						except:
							" x point is beyond curve range "
							" breaking assumption only works if tail is straight "
							breakFromFit = True
							break

				if breakFromFit:
					break

				#print len(crossOvers), "cross over events on joint", currJoint
				if len(crossOvers) == 0:
					crossOvers.append([1.0*pi/180.0 + localAng, -1.0*pi/180.0 + localAng])
				
				" if no crossover event, take the minimal distance and search the neighborhood "
				
				samples2 = []		
				candidates = []		
				for events in crossOvers:

					angle1 = events[0] * 180.0 / pi
					angle2 = events[1] * 180.0 / pi
					if angle1 > angle2:
						tempAngle = angle1
						angle1 = angle2
						angle2 = tempAngle
						
					minDist = 1e100
					candidate = []
					
					for angle in arange(angle1,angle2,0.1):
						sampAngle = angle*pi/180.0
						totalAngle = originPose[2] + sampAngle
						xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
	
						pnt = [xTotal, zTotal, totalAngle]
						
						#print pnt[0], originPose[0]
						if pnt[0] > originPose[0]:
							min, u, curvePoint = self.curve.findClosestPoint(pnt, xNotLessThan=originPose[0])
							samples2.append([sampAngle, pnt, min, u, curvePoint])

							if minDist > min:
								minDist = min
								candidate = [sampAngle, pnt, min, u, curvePoint]
				
					candidates.append(candidate)
					
				# 4. find closest of all these points and select it for the servo position
				mins2 = [currJoint]
				
				" select the point with lowest x "
				min_i = -1
				minX = 1e100
				for i in range(len(candidates)):
					pnt = candidates[i][1]
					
					if minX > pnt[0]:
						minX = pnt[0]
						min_i = i
				
				" new servo position and segment end set to new origin "
				newServo = candidates[min_i][0]  * 180.0/pi
				newOrigin = candidates[min_i][1]

				if candidates[min_i][3] >= maxU:
					minNode = currJoint
					maxU = candidates[min_i][3]

					if maxU == 1.0:
						maxNode = currJoint + 1 

				originPose = newOrigin

				# 5. kinematically update the frame of reference for the next joint and go to 1 for joint i+1 
				if maxU != 1.0:
					self.setJoint(currJoint, newServo)
					self.jointAngles[currJoint] = newServo * pi / 180.0
				else:
					break		
		
		" position of last joint iterated through "
		self.jointPositions[currJoint]= originPose
		
		" now classify the joints to their assigned peaks "
		self.classifyJoints()
				
	def step(self):
		Behavior.step(self)


		if not self.isStep():
			return False

		if self.startNode == -1 or self.endNode == -1:
			return False
		
		self.resetJoints()
		
		" load the angles of the solid joints "
		"""
		for i in range(len(self.solid)):
			if self.solid[i]:
				joints = self.getPeakJoints(i)
				for j in joints:
					self.setJoint(j,self.solidCompJoints[j]*180.0/pi)
		"""
		
		#print self.curve.control
		
	
		if self.anterior:
			self.forwardFit()

		else:
			self.backwardFit()


		self.plotCount += 1
		
		return False	
		
