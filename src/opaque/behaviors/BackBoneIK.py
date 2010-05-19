import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import ogre.renderer.OGRE as ogre
from numpy import arange
from math import *
from copy import *

class BackBoneIK:

	# 1. center the start of the curve on the base frame segment
	# 2. orient the tangent of the start of the curve to the orientation of the base frame segment
	# 3. create a single cosine curve that increases in amplitude until impacting walls

	# DRAW component
	# FIT component

	def __init__(self, probe, rootNode = 19, backwards = True):

		self.probe = probe
		self.backwards = backwards
		self.rootNode = rootNode
		self.initDraw = False
		self.numSegs = self.probe.numSegs

		self.childNodes = []
		self.childEntities = []

		# parent node for drawing purposes
		self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("curveRoot")

		self.curve = 0

	#def newCurve(self, points, angle = 0.0):
	#	pass

	#def fitBody(self):
	#	pass

	def draw(self):
		pass

	def computeError(self):

		# Two Approaches
		# 1. compute the global distance between the prescribed positions and the actual
		# 2. compute the local angle difference

		errors = []
		for i in range(0,self.rootNode+1):
			err = self.probe.getError(i)
			errors.append(err)

		#print errors
		# change the colors on the snake according to error

		# change the color of segments to reflect error
		maxError = pi/2
		for i, err in enumerate(errors):
			val = err/maxError
			self.probe.segMaterials[i].setAmbient(val,0.0,1-val)
			self.probe._ent[i].getSubEntity(0).setMaterialName("custom" + str(i))


	def updateOrigin(self, node, isGnd = False):

		# the tip of segment 0
		if not isGnd and node == -1:
			pose = self.probe.getActualJointPose(0)

			totalAngle = pose[2] + self.probe.getServo(0)
			xTotal = pose[0] - self.segLength*cos(totalAngle)
			zTotal = pose[1] - self.segLength*sin(totalAngle)
			pose[0] = xTotal
			pose[1] = zTotal
			pose[2] = totalAngle

			self.parentNode.setPosition(pose[0], 0.2, pose[1])
			newQuat = ogre.Quaternion(-pose[2], ogre.Vector3().UNIT_Y)
			self.parentNode.setOrientation(newQuat)

		elif not isGnd:
			pose = self.probe.getActualJointPose(node)
			self.parentNode.setPosition(pose[0], 0.2, pose[1])
			newQuat = ogre.Quaternion(-pose[2], ogre.Vector3().UNIT_Y)
			self.parentNode.setOrientation(newQuat)

		elif isGnd:

			pose = [0.0,0.0,0.0]

			#totalAngle = pose[2] + self.probe.getServo(0)
			#xTotal = pose[0] - self.segLength*cos(totalAngle)
			#zTotal = pose[1] - self.segLength*sin(totalAngle)
			#pose[0] = xTotal
			#pose[1] = zTotal
			#pose[2] = totalAngle

			self.parentNode.setPosition(pose[0], 0.2, pose[1])
			newQuat = ogre.Quaternion(-pose[2], ogre.Vector3().UNIT_Y)
			self.parentNode.setOrientation(newQuat)


	def getJointWRTJoint(self, originPose, originSeg, targetSeg):

		targetPose = [0.0,0.0,0.0]

		if originSeg >= self.numSegs-1 or targetSeg >= self.numSegs-1 or originSeg < -1 or targetSeg < 0:
			print "ERROR: getJointWRTJointPose joint out of range!" 
			print "Received joint", originSeg , "and" , targetSeg , "with pose"
			print originPose[0], originPose[1], originPose[2]
			raise

		# origin and target are the same, so return origin pose
		if originSeg == targetSeg :
			targetPose = copy(originPose)
			return targetPose

		# forward kinematics
		if targetSeg > originSeg:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2]

			# 
			for i in range(originSeg+1, targetSeg+1):
				xTotal = xTotal + self.segLength*cos(totalAngle)
				zTotal = zTotal + self.segLength*sin(totalAngle)
				totalAngle = totalAngle - self.probe.getServoCmd(i)

			totalAngle = normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2] 

			ind = range(targetSeg+1, originSeg + 1) # (28, 11) 
			ind.reverse()

			for i in ind:
				totalAngle = totalAngle + self.probe.getServoCmd(i)
				xTotal = xTotal - self.segLength*cos(totalAngle)
				zTotal = zTotal - self.segLength*sin(totalAngle)

			totalAngle = normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		return targetPose

	def fitBody(self):
		
		# segment origin starts
		originPose = [0.0, 0.0, 0.0]

		self.servoSettings = [0.0 for i in range(self.rootNode+1)]

		if self.backwards:

			originPose = [0.0, 0.0, pi]

			# indices
			#ind = range(0, self.rootNode+1) # (0,19+1) 
			ind = range(0, self.rootNode+1) # (0,19+1) 
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
					xTotal = originPose[0] - self.segLength*cos(totalAngle)
					zTotal = originPose[1] - self.segLength*sin(totalAngle)
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
					xTotal = originPose[0] - self.segLength*cos(totalAngle)
					zTotal = originPose[1] - self.segLength*sin(totalAngle)
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
				#self.probe.setServo(currJoint, newServo)
				self.servoSettings[currJoint] = newServo

			# TODO: weaken the joints that were not set
			#print "minNode = ", minNode, "maxU = ", maxU
			ind = range(0, minNode+1) 
			for i in ind: 
				self.servoSettings[i] = 0.0

			for i in range(0, len(self.servoSettings)):
				self.probe.setServo(i, self.servoSettings[i])

			#print self.servoSettings

		else:

			# origin is at the tip of segment 0, so we need to set the joint center to -self.segLength
			originPose = [self.segLength, 0.0, 0.0]

			#xTotal = pose[0] - self.segLength*cos(pose[2])
			#zTotal = pose[1] - self.segLength*sin(pose[2])
			#pose[0] = xTotal
			#pose[1] = zTotal

			# indices
			ind = range(0, self.rootNode+1) # (0,19+1) 
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
					xTotal = originPose[0] + self.segLength*cos(totalAngle)
					zTotal = originPose[1] + self.segLength*sin(totalAngle)
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
					print "breaking"
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
					xTotal = originPose[0] + self.segLength*cos(totalAngle)
					zTotal = originPose[1] + self.segLength*sin(totalAngle)
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
				self.probe.setServo(currJoint, newServo)
				self.servoSettings[currJoint] = newServo
			
			# TODO: weaken the joints that were not set
			ind = range(maxNode, self.rootNode+1) # (0,19+1) 
			for i in ind: 
				self.probe.setServo(i, 0.0)
				self.servoSettings[i] = 0.0

		return self.servoSettings


