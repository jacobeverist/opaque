import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *

localNodeCount = 0

class LocalCurveFit(Behavior):

	#def __init__(self, probe, startNode = 0, endNode = NUM_SEGS, curve = 0):
	def __init__(self, probe, anterior = True, backwards = False, curve = 0, cutoff = 19):
		global localNodeCount

		Behavior.__init__(self, probe)

		self.numSegs = self.probe.numSegs
		
		self.startNode = 0
		self.endNode = 39		
		
		self.curve = curve

		# anterior or posterior
		self.anterior = anterior

		# forwards or backwards
		self.backwards = backwards

		self.rootNode = cutoff

		# There are 4 possibilities to how the local fit is computed
		# if we assume that the curve is always on the extrema, and never the interior of the snake
		# Front half:
		#    - rooted at tip
		#    - rooted in interior
		# Back Half
		#    - rooted at tip
		#    - rooted in interior

		# nodes for drawing purposes
		self.nodeCount = localNodeCount
		localNodeCount += 1
		self.childNodes = []
		self.childEntities = []
		self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("localCurveRoot" + str(self.nodeCount))
		
	def setSide(self, val):
		self.anterior = val
	def setDir(self, val):
		self.backwards = val
		
	def getDir(self):
		return self.backwards
	
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

			childNode = self.parentNode.createChildSceneNode("localCurvePoint" + str(i) + "_" + str(self.nodeCount))
			self.childNodes.append(childNode)

			currEntity = self.probe._mgr.createEntity("localCurveEnt"+str(i) + "_" + str(self.nodeCount), "Cube.mesh")
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
		

	def step(self):
		Behavior.step(self)


		if not self.isStep():
			return False

		'''
		if self.startNode < 0 or self.startNode >= (self.numSegs-1) or self.endNode < 0 or self.endNode >= (self.numSegs-1):
			print "Nodes out of bounds for fitBody:", self.startNode, self.endNode
			raise

		if self.startNode == self.endNode:
			return False

		if self.startNode > self.endNode:
			backwards = True
		else:
			backwards = False
		'''
	
		if self.startNode == -1 or self.endNode == -1:
			return False

		
		self.resetJoints()

		if self.anterior:

			# segment origin starts
			originPose = [0.0, 0.0, 0.0]

			if self.backwards:

				originPose = [0.0, 0.0, pi]

				# indices
				#ind = range(0, self.rootNode+1) # (0,19+1)
				minJoint = 0
				maxJoint = self.rootNode
				
				if minJoint < self.startNode:
					minJoint = self.startNode
				if maxJoint > self.endNode:
					maxJoint = self.endNode
				
				#ind = range(0, self.rootNode+1) # (0,19+1) 
				ind = range(minJoint,maxJoint+1)
				ind.reverse()

				minNode = self.rootNode
				maxU = 0.0

				# 1. create points of the segment endpoint for the servo set from 90 to -90
				# select a point at 10 degree increments from -90 to 90
				for currJoint in ind:
					samples = []
					for angle in arange(-90,90,10):

						sampAngle = angle*pi/180.0

						totalAngle = originPose[2] + sampAngle
						xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
						pnt = [xTotal, zTotal, totalAngle]

						# 2. for each point, measure the distance to the closest point on the spline curve
						min, u, curvePoint = self.curve.findClosestPoint(pnt)

						# save to a list
						samples.append([sampAngle, pnt, min, u, curvePoint])

					# 3. select the shortest distance point and its neighbors, and create finer points between their boundaries
					min_i = -1
					minDist = 1e100
					for i in range(len(samples)):
						currMin = samples[i][2]

						if minDist > currMin:
							minDist = currMin
							min_i = i

					# select the neighbors
					upper = min_i + 1
					lower = min_i - 1
					if lower < 0:
						lower = 0
					if upper >= len(samples):
						upper = len(samples)-1

					angle1 = samples[lower][0]
					angle2 = samples[upper][0]

					# finer sampling and search
					samples2 = []
					for angle in arange(180.0/pi*angle1,180.0/pi*angle2,0.5):
						sampAngle = angle*pi/180.0
						totalAngle = originPose[2] + sampAngle
						xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
						pnt = [xTotal, zTotal, totalAngle]
						min, u, curvePoint = self.curve.findClosestPoint(pnt)
						samples2.append([sampAngle, pnt, min, u, curvePoint])

					# 4. find closest of all these points and select it for the servo position

					# find new minimum
					min_i_2 = -1
					minDist_2 = 1e100
					for i in range(len(samples2)):
						currMin = samples2[i][2]

						if minDist_2 > currMin:
							minDist_2 = currMin
							min_i_2 = i

					newServo = samples2[min_i_2][0]  * 180.0/pi
					newOrigin = samples2[min_i_2][1]

					#print "U = ", samples[min_i][3]

					if samples2[min_i_2][3] >= maxU:
						minNode = currJoint
						maxU = samples[min_i][3]

						#if maxU == 1.0:
						#	minNode = currJoint +1 

					originPose = newOrigin


					# 5. kinematically update the frame of reference for the next joint and go to 1 for joint i+1 
					self.setJoint(currJoint, newServo)

				# TODO: weaken the joints that were not set
				#print "minNode = ", minNode, "maxU = ", maxU
				#ind = range(0, minNode+1) 
				ind = range(minJoint, minNode+1) 
				for i in ind: 
					#self.setJoint(i, 0.0)
					self.setJoint(i, 0.0)

			else:

				# origin is at the tip of segment 0, so we need to set the joint center to -self.segLength
				#originPose = [self.probe.segLength, 0.0, 0.0]
				originPose = [-self.probe.segLength, 0.0, 0.0]

				#xTotal = pose[0] - self.probe.segLength*cos(pose[2])
				#zTotal = pose[1] - self.probe.segLength*sin(pose[2])
				#pose[0] = xTotal
				#pose[1] = zTotal

				
				minJoint = 0
				maxJoint = self.rootNode
				
				if minJoint < self.startNode:
					minJoint = self.startNode
				if maxJoint > self.endNode:
					maxJoint = self.endNode
					 
				ind = range(minJoint,maxJoint+1)
				
				# indices
				#ind = range(0, self.rootNode+1) # (0,19+1) 
				#ind = range(1, 2) # (0,19+1) 
				maxNode = 0
				maxU = 0.0

				# 1. create points of the segment endpoint for the servo set from 90 to -90
				# select a point at 10 degree increments from -90 to 90
				for currJoint in ind:
					samples = []
					for angle in arange(-90,90,10):

						sampAngle = angle*pi/180.0

						totalAngle = originPose[2] - sampAngle
						xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)
						#totalAngle = totalAngle - sampAngle
						pnt = [xTotal, zTotal, totalAngle]

						# 2. for each point, measure the distance to the closest point on the spline curve
						min, u, curvePoint = self.curve.findClosestPoint(pnt)

						# save to a list
						samples.append([sampAngle, pnt, min, u, curvePoint])

					#for samp in samples:
					#	print samp

					# 3. select the shortest distance point and its neighbors, and create finer points between their boundaries
					min_i = -1
					minDist = 1e100
					for i in range(len(samples)):
						currMin = samples[i][2]

						if minDist > currMin:
							minDist = currMin
							min_i = i
							#print minDist, min_i, maxU, currJoint, samples[i][0]

					# make sure there's still length of curve to be fitted,
					# otherwise break from the fitting process
					#print samples[min_i][3] , maxU
					if samples[min_i][3] >= maxU:
						maxNode = currJoint
						maxU = samples[min_i][3]
						#print maxNode, maxU
					else:
						#print "breaking"
						break

					# select the neighbors
					upper = min_i + 1
					lower = min_i - 1
					if lower < 0:
						lower = 0
					if upper >= len(samples):
						upper = len(samples)-1

					angle1 = samples[lower][0]
					angle2 = samples[upper][0]

					# finer sampling and search
					samples2 = []
					for angle in arange(180.0/pi*angle1,180.0/pi*angle2,0.5):
						sampAngle = angle*pi/180.0

						totalAngle = originPose[2] - sampAngle
						xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)
						#totalAngle = originPose[2] - sampAngle

						pnt = [xTotal, zTotal, totalAngle]
						min, u, curvePoint = self.curve.findClosestPoint(pnt)
						samples2.append([sampAngle, pnt, min, u, curvePoint])

					# 4. find closest of all these points and select it for the servo position

					# find new minimum
					min_i_2 = -1
					minDist_2 = 1e100
					for i in range(len(samples2)):
						currMin = samples2[i][2]

						if minDist_2 > currMin:
							minDist_2 = currMin
							min_i_2 = i

					newServo = samples2[min_i_2][0]  * 180.0/pi
					newOrigin = samples2[min_i_2][1]

					originPose = newOrigin

					# 5. kinematically update the frame of reference for the next joint and go to 1 for joint i+1 
					self.setJoint(currJoint, newServo)

				# TODO: weaken the joints that were not set
				#ind = range(maxNode, self.rootNode+1) # (0,19+1) 
				ind = range(maxNode, maxJoint+1) # (0,19+1) 
				for i in ind: 
					self.setJoint(i, 0.0)

			return False	
		else:

			# segment origin starts
			originPose = [0.0, 0.0, 0.0]

			if self.backwards:

				#originPose = [0.0, 0.0, pi]
				originPose = [-self.probe.segLength, 0.0, pi]

				# indices
				minJoint = self.rootNode+1
				maxJoint = self.probe.numSegs-2
				
				if minJoint < self.startNode:
					minJoint = self.startNode
				if maxJoint > self.endNode:
					maxJoint = self.endNode
					 
				ind = range(minJoint,maxJoint+1)
				#ind = range(self.rootNode+1, self.probe.numSegs-1) # (20,39+1) 
				ind.reverse()

				minNode = self.probe.numSegs-2
				maxU = 0.0

				# 1. create points of the segment endpoint for the servo set from 90 to -90
				# select a point at 10 degree increments from -90 to 90
				for currJoint in ind:
					samples = []
					for angle in arange(-90,90,10):

						sampAngle = angle*pi/180.0

						totalAngle = originPose[2] + sampAngle
						xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
						pnt = [xTotal, zTotal, totalAngle]

						# 2. for each point, measure the distance to the closest point on the spline curve
						min, u, curvePoint = self.curve.findClosestPoint(pnt)

						# save to a list
						samples.append([sampAngle, pnt, min, u, curvePoint])

					# 3. select the shortest distance point and its neighbors, and create finer points between their boundaries
					min_i = -1
					minDist = 1e100
					for i in range(len(samples)):
						currMin = samples[i][2]

						if minDist > currMin:
							minDist = currMin
							min_i = i

					# select the neighbors
					upper = min_i + 1
					lower = min_i - 1
					if lower < 0:
						lower = 0
					if upper >= len(samples):
						upper = len(samples)-1

					angle1 = samples[lower][0]
					angle2 = samples[upper][0]

					# finer sampling and search
					samples2 = []
					for angle in arange(180.0/pi*angle1,180.0/pi*angle2,0.5):
						sampAngle = angle*pi/180.0
						totalAngle = originPose[2] + sampAngle
						xTotal = originPose[0] - self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] - self.probe.segLength*sin(totalAngle)
						pnt = [xTotal, zTotal, totalAngle]
						min, u, curvePoint = self.curve.findClosestPoint(pnt)
						samples2.append([sampAngle, pnt, min, u, curvePoint])

					# 4. find closest of all these points and select it for the servo position

					# find new minimum
					min_i_2 = -1
					minDist_2 = 1e100
					for i in range(len(samples2)):
						currMin = samples2[i][2]

						if minDist_2 > currMin:
							minDist_2 = currMin
							min_i_2 = i

					newServo = samples2[min_i_2][0]  * 180.0/pi
					newOrigin = samples2[min_i_2][1]

					#print "U = ", samples[min_i][3]

					if samples2[min_i_2][3] >= maxU:
						minNode = currJoint
						maxU = samples[min_i][3]

						#if maxU == 1.0:
						#	minNode = currJoint +1 

					originPose = newOrigin


					# 5. kinematically update the frame of reference for the next joint and go to 1 for joint i+1 
					self.setJoint(currJoint, newServo)

				# TODO: weaken the joints that were not set
				#print "minNode = ", minNode, "maxU = ", maxU
				#ind = range(self.rootNode+1, minNode+1) 
				ind = range(minJoint, minNode+1) 
				for i in ind: 
					self.setJoint(i, 0.0)

			else:

				# origin is at the tip of segment 0, so we need to set the joint center to -self.segLength
				originPose = [0.0, 0.0, 0.0]

				#xTotal = pose[0] - self.probe.segLength*cos(pose[2])
				#zTotal = pose[1] - self.probe.segLength*sin(pose[2])
				#pose[0] = xTotal
				#pose[1] = zTotal

				# indices
				minJoint = self.rootNode+1
				maxJoint = self.probe.numSegs-2
				
				if minJoint < self.startNode:
					minJoint = self.startNode
				if maxJoint > self.endNode:
					maxJoint = self.endNode
					 
				ind = range(minJoint,maxJoint+1)
				#ind = range(self.rootNode+1, self.probe.numSegs-1) # (20,39+1) 
				maxNode = 0
				maxU = 0.0

				# 1. create points of the segment endpoint for the servo set from 90 to -90
				# select a point at 10 degree increments from -90 to 90
				for currJoint in ind:
					samples = []
					for angle in arange(-90,90,10):

						sampAngle = angle*pi/180.0

						totalAngle = originPose[2] - sampAngle
						xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)
						#totalAngle = totalAngle - sampAngle
						pnt = [xTotal, zTotal, totalAngle]

						# 2. for each point, measure the distance to the closest point on the spline curve
						min, u, curvePoint = self.curve.findClosestPoint(pnt)

						# save to a list
						samples.append([sampAngle, pnt, min, u, curvePoint])

					#for samp in samples:
					#	print samp

					# 3. select the shortest distance point and its neighbors, and create finer points between their boundaries
					min_i = -1
					minDist = 1e100
					for i in range(len(samples)):
						currMin = samples[i][2]

						if minDist > currMin:
							minDist = currMin
							min_i = i
							#print minDist, min_i, maxU, currJoint, samples[i][0]

					# make sure there's still length of curve to be fitted,
					# otherwise break from the fitting process
					#print samples[min_i][3] , maxU
					if samples[min_i][3] >= maxU:
						maxNode = currJoint
						maxU = samples[min_i][3]
						#print maxNode, maxU
					else:
						#print "breaking"
						break

					# select the neighbors
					upper = min_i + 1
					lower = min_i - 1
					if lower < 0:
						lower = 0
					if upper >= len(samples):
						upper = len(samples)-1

					angle1 = samples[lower][0]
					angle2 = samples[upper][0]

					# finer sampling and search
					samples2 = []
					for angle in arange(180.0/pi*angle1,180.0/pi*angle2,0.5):
						sampAngle = angle*pi/180.0

						totalAngle = originPose[2] - sampAngle
						xTotal = originPose[0] + self.probe.segLength*cos(totalAngle)
						zTotal = originPose[1] + self.probe.segLength*sin(totalAngle)
						#totalAngle = originPose[2] - sampAngle

						pnt = [xTotal, zTotal, totalAngle]
						min, u, curvePoint = self.curve.findClosestPoint(pnt)
						samples2.append([sampAngle, pnt, min, u, curvePoint])

					# 4. find closest of all these points and select it for the servo position

					# find new minimum
					min_i_2 = -1
					minDist_2 = 1e100
					for i in range(len(samples2)):
						currMin = samples2[i][2]

						if minDist_2 > currMin:
							minDist_2 = currMin
							min_i_2 = i

					newServo = samples2[min_i_2][0]  * 180.0/pi
					newOrigin = samples2[min_i_2][1]

					originPose = newOrigin

					# 5. kinematically update the frame of reference for the next joint and go to 1 for joint i+1 
					self.setJoint(currJoint, newServo)

				# TODO: weaken the joints that were not set
				#ind = range(maxNode, self.numSegs-1) # (0,19+1) 
				ind = range(maxNode, maxJoint+1) # (0,19+1) 
				for i in ind: 
					self.setJoint(i, 0.0)

			return False				

			#print "posterior not implemented yet"



		return False