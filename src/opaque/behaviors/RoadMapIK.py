import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import ogre.renderer.OGRE as ogre
from math import *
from copy import *
import pylab

class RoadMapIK:

	def __init__(self, probe, roadMap):

		self.nodeCount = 0

		# 1. load the shadow map
		# 2. compute the voronoi graph
		# 3. plot the vertices to the 3D world
		self.probe = probe
		self.numSegs = self.probe.numSegs

		self.childNodes = []
		self.childEntities = []

		self.nodeCount += 1

		# parent node for drawing purposes
		# FIXME:  we need a destructor to kill these nodes
		self.parentNode = self.probe._mgr.getRootSceneNode().createChildSceneNode("rootCurveRoot" + str(self.nodeCount))

		self.curve = 0
		self.rootNode = 19

		self.resetRoadMap(roadMap)

	def resetRoadMap(self, roadMap):

		self.rm = roadMap

		self.realPath = []
		self.draw()
		self.segLength = 0.15
		self.closePoints = []
		self.curve = 0

	def gotoLocation(self, start, goal):
		self.realPath = self.rm.computePath(copy(start),copy(goal))
		self.curve = VoronoiFit(self.realPath) 
		self.draw()

	def setPath(self, path):
		self.realPath = path
		#print "setting path:", self.realPath
		self.curve = VoronoiFit(self.realPath) 
		self.draw()

	def printGraph(self):

		self.curve.draw()

		xP = []
		yP = []
		for p in self.closePoints:
			xP.append(p[0])
			yP.append(p[1])
		pylab.scatter(xP,yP)

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


		# draw the points in the simulation
		for i in range(len(self.realPath)):

			# take samples from the fitted curve
			pnt = copy(self.realPath[i])

			childNode = self.parentNode.createChildSceneNode("roadCurvePoint" + str(self.nodeCount) + "_" + str(i))
			self.childNodes.append(childNode)

			currEntity = self.probe._mgr.createEntity("roadCurveEnt" + str(self.nodeCount) + "_" +str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)

			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			#position = ogre.Vector3(0.0,0.0,0.0)
			childNode.setPosition(position)

			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)

			childNode.attachObject(currEntity)

	def getPathDirection(self):
		joints = range(NUM_SEGS-1)

		# segment angle vector
		vecSum1 = [0.0,0.0]

		# path vector
		vecSum2 = [0.0,0.0]

		for j in joints:

			origin = self.probe.getActualJointPose(j)
			poseVec = [cos(origin[2]),sin(origin[2])] 
			vecSum1[0] += poseVec[0]
			vecSum1[1] += poseVec[1]

			dist, linePoint, orientVec = self.curve.findClosestPoint([origin[0], origin[1]])
			#print linePoint
			vecSum2[0] += orientVec[0]
			vecSum2[1] += orientVec[1]

		# renormalize the snake orientation vector
		mag1 = sqrt(vecSum1[0]**2 + vecSum1[1]**2)
		vecSum1 = [vecSum1[0]/mag1,vecSum1[1]/mag1]

		#print "snake vector:", vecSum1
		#print "path vector:", vecSum2
		# compute angle of snake orientation average
		ang1 = acos(vecSum1[0])
		if asin(vecSum1[1]) < 0:
			ang1 *= -1

		# rotate path vector by the same angle
		vecSum2 = [vecSum2[0]*cos(ang1) + vecSum2[1]*sin(ang1), -vecSum2[0]*sin(ang1) + vecSum2[1]*cos(ang1)]

		# renormalize the path orientation vector
		mag2 = sqrt(vecSum2[0]**2 + vecSum2[1]**2)
		vecSum2 = [vecSum2[0]/mag2,vecSum2[1]/mag2]

		# compute the angle of the path orientation average wrt snake orientation average
		ang2 = acos(vecSum2[0])
		if asin(vecSum2[1]) < 0:
			ang2 *= -1

		#print "snake orientation:", ang1
		#print "PATH DIRECTION:", ang2

		# forward, path along ascending joints
		if fabs(ang2) < pi/2:
			return True
		# backward, path along descending joints
		else:
			return False

	# perform a continuous transition to a continuous fit
	def fitBody(self, startNode = 19, endNode = 0 ):

		# FIXME:  if the tail joint is controlled and is also the start node,
		# we need to reverse the control because normal control will cause the 
		# tail to curl into a ball

		# cases:
		# 1. startNode == endNode => do nothing
		# 2. startNode or endNode out of bounds => raise exception
		# 3. startNode > endNode => backwards fitting
		# 4. startNode < endNode => forwards fitting

		#if startNode == endNode:
		#	return

		if startNode < 0 or startNode >= (self.numSegs-1) or endNode < 0 or endNode >= (self.numSegs-1):
			print "Nodes out of bounds for fitBody:", startNode, endNode
			raise

		if startNode > endNode:
			backwards = True
		else:
			backwards = False

		#print startNode, endNode

		origin = self.probe.getActualJointPose(startNode)
		xTotal = origin[0]
		yTotal = origin[1]
		totalAngle = origin[2]
		totalVec = [cos(totalAngle),sin(totalAngle)]

		self.closePoints = []
		finalAngles = []

		# 1.  terminator point and orientation
		tPoint = self.curve.getU(1.0)
		tVec = self.curve.getUVector(1.0)
		tAngle = acos(tVec[0])
		if asin(tVec[1]) < 0:
			tAngle = -tAngle

		if backwards:

			# joints we're visiting in reverse order
			if startNode  >= 35:
				#joints = range(5, endNode+1)
				joints = range(endNode, 34+1)
				joints.reverse()

				tEnd = endNode
				if tEnd < 35:
					tEnd = 35

				zeroJoints = range(tEnd,startNode+1)
				for i in zeroJoints:
					self.probe.setServo(i, 0.0)

			else:
				joints = range(endNode, startNode+1)
				joints.reverse()

			#joints = range(endNode, startNode+1)
			#joints.reverse()

			# compute the change in joint angle for every consecutive joint
			for i in joints:

				# find the closest point on the roadmap curve for the joint's location
				# - returns the distance, the point that is closest, the orientation of the curve
				dist, linePoint, orientVec = self.curve.findClosestPoint([xTotal, yTotal])

				# check if we've gone beyond the boundaries of the curve
				# by measuring the orthogonal distance of the current joint to the vector at U=1.0
				relPoint = [xTotal-tPoint[0], yTotal-tPoint[1]]
				rotPoint = [relPoint[0]*cos(tAngle) + relPoint[1]*sin(tAngle), -relPoint[0]*sin(tAngle) + relPoint[1]*cos(tAngle)]
				orthoDist = fabs(rotPoint[1])

				# if orthogonal point is in the positive direction of the Uvector
				# and less distance than the closest point already computed
				if rotPoint[0] > 0 and orthoDist < dist:
					#print "distance smaller: ", orthoDist, dist

					# compute the closest point on the vector
					closePoint = [rotPoint[0], 0.0]
					closePoint = [closePoint[0]*cos(tAngle) - closePoint[1]*sin(tAngle), closePoint[0]*sin(tAngle) + closePoint[1]*cos(tAngle)]
					closePoint = [closePoint[0] + tPoint[0], closePoint[1] + tPoint[1]]

					#self.closePoints.append(copy(closePoint))

					# orientation vector is now the vector at U=1.0
					orientVec = copy(tVec) 

					# location vector is now computed to the orthogonal point
					locVec = [closePoint[0]-xTotal, closePoint[1]-yTotal]
					locMag = sqrt(locVec[0]**2 + locVec[1]**2)
					locVec = [locVec[0]/locMag, locVec[1]/locMag]

					dist = orthoDist

				# now treat everything regularly
				else:

					# closest point is point found on interior of roadmap curve
					#self.closePoints.append(copy(linePoint))

					# location vector is computed to point on interior of roadmap
					locVec = [linePoint[0]-xTotal, linePoint[1]-yTotal]
					locMag = sqrt(locVec[0]**2 + locVec[1]**2)
					locVec = [locVec[0]/locMag, locVec[1]/locMag]

				# negate the orientation and location vector since the backwards kinematics
				# require us to align the snake segments to 180 degrees from the desired orientation
				# negating the orientation vectors allows us to use the same algorithm for forwards and backwards
				orientVec[0] *= -1
				orientVec[1] *= -1
				locVec[0] *= -1
				locVec[1] *= -1

				# weighted direction we want to point to given orientation and location vectors
				desVec = [(orientVec[0] + dist*locVec[0])/(dist+1.0), (orientVec[1] + dist*locVec[1])/(dist+1.0)]
				desMag = sqrt(desVec[0]**2 + desVec[1]**2)
				desVec = [desVec[0]/desMag, desVec[1]/desMag]

				# Compute the required rotation from current orientation to desired orientation (totalVec, desVec)

				# 1. rotate totalVec to [1,0]
				rotAngle = acos(totalVec[0])
				if asin(totalVec[1]) < 0:
					rotAngle = -rotAngle

				# 2. rotate desVec by the same rotation as totalVec
				desRotVec = [desVec[0]*cos(rotAngle) + desVec[1]*sin(rotAngle), -desVec[0]*sin(rotAngle) + desVec[1]*cos(rotAngle)]
				desMag = sqrt(desRotVec[0]**2 + desRotVec[1]**2)
				desRotVec[0] /= desMag
				desRotVec[1] /= desMag

				# 3. compute angle of desRotVec and this is our desired joint angle
				desAngle = acos(desRotVec[0])
				if asin(desRotVec[1]) < 0:
					desAngle = -desAngle

				k = 0.1
				angle = self.probe.getServoCmd(i)
				#angle = self.probe.getServo(i)
				finalAngle = k*desAngle + (1-k)*angle

				# FIXME: remove this
				#finalAngle = desAngle

				self.probe.setServo(i, finalAngle*180.0/pi)
				finalAngles.append(finalAngle)

				totalVec = [cos(finalAngle), sin(finalAngle)]
				totalVec = [totalVec[0]*cos(rotAngle) - totalVec[1]*sin(rotAngle), totalVec[0]*sin(rotAngle) + totalVec[1]*cos(rotAngle)]

				#totalVec = copy(desVec)
				#totalAngle += finalAngle
				#xTotal = xTotal + self.segLength*cos(totalAngle)
				#yTotal = yTotal + self.segLength*sin(totalAngle)

				self.closePoints.append(copy([xTotal,yTotal]))
				totalAngle += finalAngle
				xTotal = xTotal - self.probe.segLength*cos(totalAngle)
				yTotal = yTotal - self.probe.segLength*sin(totalAngle)

			#print finalAngles

		# forward fitting
		# NOTE:  in this case we are modifying the current frame instead of modifying the next frame
		else:

			#joints = range(0,self.rootNode+1)
			if startNode  <= 4:
				joints = range(5, endNode+1)

				tEnd = endNode
				if tEnd > 4:
					tEnd = 4

				zeroJoints = range(startNode,tEnd+1)
				for i in zeroJoints:
					self.probe.setServo(i, 0.0)

			else:
				joints = range(startNode, endNode+1)

			for i in joints:

				dist, linePoint, orientVec = self.curve.findClosestPoint([xTotal, yTotal])

				# check if we've gone beyond the boundaries of the curve
				relPoint = [xTotal-tPoint[0], yTotal-tPoint[1]]
				rotPoint = [relPoint[0]*cos(tAngle) + relPoint[1]*sin(tAngle), -relPoint[0]*sin(tAngle) + relPoint[1]*cos(tAngle)]
				orthoDist = fabs(rotPoint[1])
				if rotPoint[0] > 0 and orthoDist < dist:
					#print "distance smaller: ", orthoDist, dist

					closePoint = [rotPoint[0], 0.0]
					closePoint = [closePoint[0]*cos(tAngle) - closePoint[1]*sin(tAngle), closePoint[0]*sin(tAngle) + closePoint[1]*cos(tAngle)]
					closePoint = [closePoint[0] + tPoint[0], closePoint[1] + tPoint[1]]

					#self.closePoints.append(copy(closePoint))

					orientVec = copy(tVec) 
					locVec = [closePoint[0]-xTotal, closePoint[1]-yTotal]
					locMag = sqrt(locVec[0]**2 + locVec[1]**2)
					locVec = [locVec[0]/locMag, locVec[1]/locMag]

					dist = orthoDist
				else:

					#self.closePoints.append(copy(linePoint))

					# compute the direction from current joint to line point
					locVec = [linePoint[0]-xTotal, linePoint[1]-yTotal]
					locMag = sqrt(locVec[0]**2 + locVec[1]**2)
					locVec = [locVec[0]/locMag, locVec[1]/locMag]

				#print "joint =", i 
				#print "xTotal =", xTotal
				#print "yTotal =", yTotal
				#print "totalAngle =", totalAngle
				#print "totalVec =", totalVec
				#print "dist =", dist
				#print "orthoDist =", orthoDist
				#print "locVec =", locVec
				#print "orientVec =", orientVec

				#orientVec[0] *= -1
				#orientVec[1] *= -1

				# weighted direction we want to point to
				desVec = [(orientVec[0] + dist*locVec[0])/(dist+1.0), (orientVec[1] + dist*locVec[1])/(dist+1.0)]
				desMag = sqrt(desVec[0]**2 + desVec[1]**2)
				desVec = [desVec[0]/desMag, desVec[1]/desMag]

				#print "desVec =", desVec

				# Compute the rotation to desVec from totalVec and this is our desired joint angle

				# 1. rotate totalVec to [1,0]
				rotAngle = acos(totalVec[0])
				if asin(totalVec[1]) < 0:
					rotAngle = -rotAngle

				#print "rotAngle =", rotAngle

				# 2. rotate desVec by the same rotation as totalVec
				desRotVec = [desVec[0]*cos(rotAngle) + desVec[1]*sin(rotAngle), -desVec[0]*sin(rotAngle) + desVec[1]*cos(rotAngle)]
				desMag = sqrt(desRotVec[0]**2 + desRotVec[1]**2)
				desRotVec[0] /= desMag
				desRotVec[1] /= desMag

				#print "desRotVec =", desRotVec

				# 3. compute angle of desRotVec and this is our desired offset from the current joint angle
				offsetAngle = acos(desRotVec[0])
				if asin(desRotVec[1]) < 0:
					offsetAngle = -offsetAngle

				#print "offsetAngle =", offsetAngle

				#print "cmdAngle =", self.probe.getServoCmd(i)

				# offsetAngle must be negative because we are going forwards
				desAngle = -offsetAngle + self.probe.getServoCmd(i)
				desAngle = normalizeAngle(desAngle)
				#print "desAngle =", desAngle

				# negate desAngle since we are now going forwards
				#desAngle *= -1.0

				k = 0.1
				angle = self.probe.getServoCmd(i)
				#angle = self.probe.getServo(i)
				finalAngle = k*desAngle + (1-k)*angle
				diffAngle = finalAngle-angle

				#print "finalAngle =", finalAngle

				# FIXME: remove this
				#finalAngle = desAngle

				self.probe.setServo(i, finalAngle*180.0/pi)
				finalAngles.append(finalAngle)

				# negate because we are going forward
				#finalAngle *= -1

				#totalVec = [cos(finalAngle), sin(finalAngle)]

				if i < endNode:
					rotAngle -= diffAngle
					nextAngle = self.probe.getServoCmd(i+1)
					totalVec = [cos(-nextAngle), sin(-nextAngle)]
					totalVec = [totalVec[0]*cos(rotAngle) - totalVec[1]*sin(rotAngle), totalVec[0]*sin(rotAngle) + totalVec[1]*cos(rotAngle)]

					self.closePoints.append(copy([xTotal,yTotal]))
					totalAngle -= diffAngle
					#totalAngle -= finalAngle
					#print "desAngle =", desAngle, "angle =", angle, "finalAngle =", finalAngle, "totalAngle =", totalAngle
					#print "totalVec =", totalVec
					xTotal = xTotal + self.probe.segLength*cos(totalAngle)
					yTotal = yTotal + self.probe.segLength*sin(totalAngle)
					totalAngle -= nextAngle
					#print "xTotal, yTotal =", xTotal, yTotal


			#print finalAngles


	def updateOrigin(self, node, isGnd = False):
		# do nothing
		pass

