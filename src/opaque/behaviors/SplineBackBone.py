import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from BackBoneIK import *

import ogre.renderer.OGRE as ogre
from math import *


class SplineBackBone(BackBoneIK):

	def __init__(self, probe, points, rootNode = 19, backwards = True):
		BackBoneIK.__init__(self, probe, rootNode, backwards)

		self.points = points
		self.newCurve(points)

	# 3. rotation of curve around startJoint
	def rotateCurve(self, angle):

		# rotate all the points and rebuild the curve
		newPoints = []
		for point in self.points:
			xRot = point[0]*cos(angle) + point[1]*sin(angle)
			zRot = -point[0]*sin(angle) + point[1]*cos(angle)
			newPoints.append([xRot,zRot])

		# FIXME: DISABLED
		#self.points = newPoints
		#self.curve = SplineFit(self.points)
	
	def newCurve(self, points, angle = 0.0):

		#self.computeError()

		self.points = points

		# rotate all the points and rebuild the curve
		#newPoints = []
		#for point in self.points:
		#	xRot = point[0]*cos(angle) + point[1]*sin(angle)
		#	zRot = -point[0]*sin(angle) + point[1]*cos(angle)
		#	newPoints.append([xRot,zRot])

		# 1. constructor
		# 2. getU
		# 3. getUVector
		# 4. findClosestPoint

		self.curve = SplineFit(self.points)

		origin = self.curve.getU(0.0)

		# rest the origin to 0,0
		newPoints = []
		for point in self.points:
			xRot = point[0] - origin[0]
			zRot = point[1] - origin[1] 
			newPoints.append([xRot,zRot])
		self.points = newPoints
		self.curve = SplineFit(self.points)

		# rotate to [1.0,0.0]
		tanVec = self.curve.getUVector(0.0)
		tanVec = [newPoints[1][0]-newPoints[0][0], newPoints[1][1]-newPoints[0][1]]

		#print "tanVec: ", tanVec
		# normalize
		mag = sqrt(tanVec[0]**2 + tanVec[1]**2)
		#print "mag =", mag
		tanVec = [tanVec[0]/mag, tanVec[1]/mag]
		#print "tanVec: ", tanVec

		x = tanVec[0]
		y = tanVec[1]
		if y != 0.0:
			B = 1 / (x*x/y + y)
			A = B * x / y
		else:
			if x > 0:
				A = 1.0
				B = 0.0
			else:
				A = -1.0
				B = 0.0

		newPoints = []
		for point in self.points:
			xRot = point[0]*A + point[1]*B
			zRot = -point[0]*B + point[1]*A
			newPoints.append([xRot,zRot])

		self.points = newPoints
		#print "new points3: ", self.points
		self.curve = SplineFit(self.points)
	
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
		for i in range(len(self.points)):

			# take samples from the fitted curve
			pnt = self.curve.getU(i/float(len(self.points)))

			childNode = self.parentNode.createChildSceneNode("curvePoint" + str(i))
			self.childNodes.append(childNode)

			currEntity = self.probe._mgr.createEntity("curveEnt"+str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)

			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			#position = ogre.Vector3(0.0,0.0,0.0)
			childNode.setPosition(position)

			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)

			childNode.attachObject(currEntity)


