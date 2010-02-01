import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from math import *
from copy import *

class FakeProbe:

	def __init__(self, numSegs, segLength, segWidth):
		self.numSegs = numSegs
		self.segLength = segLength
		self.segWidth = segWidth
		
		self.currIndex = 0
		self.angles = [[] for i in range(39)]
		
		self.walls = []

	def getNumJoints(self):
		return self.numSegs - 1

	def getActualJointPose(self, i):
		return [2.0,2.0,0.0]
					
	def loadFile(self, num):
		angleFile = open("angle_output_%04u.txt" % num, 'r')

		self.angles = [[] for i in range(39)]
	
		for line in angleFile:
			angleSnapshot = eval(line)
			for i in range(len(angleSnapshot)):
				self.angles[i].append(angleSnapshot[i])
		
		self.currIndex = 0
		for i in range(len(self.angles)):
			print len(self.angles[i])
		#print len(self.angles[0])
		
	def loadEmpty(self):

		self.angles = [[] for i in range(39)]
	
		for j in range(5000):
			angleSnapshot = [0.0 for i in range(39)]
			for i in range(len(angleSnapshot)):
				self.angles[i].append(angleSnapshot[i])
		
		self.currIndex = 0
			
	def step(self):
		self.currIndex += 1
		
		#print self.currIndex
		if self.currIndex >= len(self.angles[0]):
			raise

	def getServo(self, index):
		#print "index =", index
		if index < self.numSegs :
			return self.angles[index][self.currIndex]

		raise
	
	def getError(self, index):
		return 0.0
	
	def addWalls(self, walls):
		self.walls = copy(walls)
		
	def getWalls(self):
		return self.walls

	def getJointWRTJointPose(self, originPose, originJoint, targetJoint):

		targetPose = [0.0,0.0,0.0]

		if originJoint >= self.numSegs-1 or targetJoint >= self.numSegs or originJoint < 0 or targetJoint < -1:
			print "ERROR: getJointWRTJointPose joint out of range!" 
			print "Received joint", originJoint , "and" , targetJoint , "with pose"
			print originPose[0], originPose[1], originPose[2]
			raise

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			targetPose = copy(originPose)
			return targetPose

		# forward kinematics
		if targetJoint > originJoint:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2]

			# 
			for i in range(originJoint+1, targetJoint+1):
				xTotal = xTotal + self.segLength*cos(totalAngle)
				zTotal = zTotal + self.segLength*sin(totalAngle)
				if i != self.numSegs-1:
					totalAngle = totalAngle - self.getServo(i)

			totalAngle = normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2] 

			ind = range(targetJoint+1, originJoint + 1) # (28, 11) 
			ind.reverse()

			for i in ind:
				totalAngle = totalAngle + self.getServo(i)
				xTotal = xTotal - self.segLength*cos(totalAngle)
				zTotal = zTotal - self.segLength*sin(totalAngle)

			totalAngle = normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		return targetPose
