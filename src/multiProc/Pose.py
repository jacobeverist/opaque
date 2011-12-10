
from functions import *
from copy import copy
from numpy import array, dot, transpose
from math import *

class Pose:
	
	def __init__(self, pose = [0.0,0.0,0.0]):
		self.setEstPose(pose)

	def setEstPose(self, newPose):

	   self.estPose = copy(newPose)
	   self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)

	   if self.dist != 0:
		   self.vecAng = acos(self.estPose[0]/self.dist)
		   if asin(self.estPose[1]/self.dist) < 0:
				   self.vecAng = -self.vecAng
	   else:
			self.vecAng = 0.0

	   self.backR = array([[cos(self.vecAng), sin(self.vecAng)],[-sin(self.vecAng),cos(self.vecAng)]])
	   self.foreR = array([[cos(self.vecAng), -sin(self.vecAng)],[sin(self.vecAng),cos(self.vecAng)]])

	   self.R = array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])

	def doInverse(self, offset):
		
		invOffset = [0.0,0.0,0.0]

		offsetR = array([[cos(offset[2]), sin(offset[2])],[-sin(offset[2]),cos(offset[2])]])
		finalVec = array([[offset[0]], [offset[1]]])
		transVec = dot(offsetR, finalVec)
		
		transVec = -transVec
		#print transVec
		#print -transVec
		
		invOffset[0] = transVec[0, 0]
		invOffset[1] = transVec[1, 0]
		invOffset[2] = -offset[2]
	   	
		#print invOffset
		
		return invOffset
	

	def convertLocalOffsetToGlobal(self, offset):

	   globalEst = [0.0,0.0,0.0]
	   
	   #print "converting offset:", offset

	   finalVec = array([[offset[0]], [offset[1]]])
	   transVec = dot(transpose(self.R), finalVec)
	   resVec = dot(self.backR, transVec)
	   resVec[0, 0] += self.dist
	   tempVec = dot(self.foreR, resVec)
	   globalEst[0] = tempVec[0, 0]
	   globalEst[1] = tempVec[1, 0]
	   globalEst[2] = normalizeAngle(self.estPose[2] + offset[2])

	   return globalEst

	def convertGlobalPoseToLocal(self, pose):

	   " transform pnt to local coordinates"
	   globalVec = array([[pose[0]],[pose[1]]])

	   " perform translation correction "
	   tempVec = dot(self.backR, globalVec)
	   tempVec[0,0] -= self.dist
	   transVec = dot(self.foreR, tempVec)

	   " now, apply rotation correction with respect to origin "
	   localVec = dot(self.R, transVec)

	   localPose = [localVec[0,0], localVec[1,0], normalizeAngle(pose[2] - self.estPose[2])]

	   return localPose

	def convertLocalToGlobal(self, pnt):

	   finalVec = array([[pnt[0]], [pnt[1]]])
	   transVec = dot(transpose(self.R), finalVec)
	   resVec = dot(self.backR, transVec)
	   resVec[0, 0] += self.dist
	   tempVec = dot(self.foreR, resVec)

	   newPoint = [tempVec[0,0],tempVec[1,0]]

	   return newPoint

	def convertGlobalToLocal(self, pnt):

	   " transform pnt to local coordinates"
	   globalVec = array([[pnt[0]],[pnt[1]]])

	   " perform translation correction "
	   tempVec = dot(self.backR, globalVec)
	   tempVec[0,0] -= self.dist
	   transVec = dot(self.foreR, tempVec)

	   " now, apply rotation correction with respect to origin "
	   localVec = dot(self.R, transVec)

	   newPoint = [localVec[0,0], localVec[1,0]]
	   return newPoint
							