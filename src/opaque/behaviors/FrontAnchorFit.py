
from Behavior import *
import traceback
import math

fastLocalNodeCount = 0

class FrontAnchorFit(Behavior):

	def __init__(self, probe, startNode = 0, endNode = 38):
		global fastLocalNodeCount

		Behavior.__init__(self, probe)

		self.numSegs = self.probe.numSegs
		
		self.startNode = startNode
		self.endNode = endNode	
		self.centerNode = endNode - math.floor(endNode-startNode)/2
		
		self.curve = 0
	
	def setCurve(self, curve):
		self.curve = curve

	def setBoundaries(self, startNode, endNode):
		self.startNode = startNode
		self.endNode = endNode
		self.centerNode = endNode - math.floor(endNode-startNode)/2

	def getPeakJoints(self):
		if self.endNode > self.startNode:
			return range(self.startNode, self.endNode+1)
		else:
			return range(self.endNode, self.startNode+1)
			

	def getPeakErrorJoints(self, index):
		return self.jointErrorClasses[index]
		
	def forwardFit(self):
		
		# origin is at the tip of segment 0, so we need to set the joint center to -self.segLength
		originPose = [0.0, 0.0, 0.0]
		
		minJoint = 0
		maxJoint = self.rootNode
		
		if minJoint < self.startNode:
			minJoint = self.startNode
		if maxJoint > self.endNode:
			maxJoint = self.endNode

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
			" this is to allow segments spanning two neighoring peaks to adjust freely "
			if solidIndex >= 0 and self.solid[solidIndex]:
			#if False:
				
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


	def getPeakError(self):
					
		cmdPoints = []
		actPoints = []
		
		joints = self.getPeakJoints()
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
		
