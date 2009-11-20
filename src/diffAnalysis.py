#!/usr/bin/python


import os
import csv
import numpy
import sys
from math import *
import Image
from matplotlib.pyplot import *
from matplotlib.collections import LineCollection
import pylab

"""
class DiffStability:

	def __init__(self, thresh):
		self.thresh = thresh

		self.currAngle = 0.0
		self.nomAngle = 0.0
		self.isInit = False
		self.count = 0
		
		
		#1. set nominal angle
		#2. for each iteration, check if current angle violates difference threshold
		#3. otherwise set as stable
		#4. reset
			
		
	def getDiff(self):
		return self.currAngle - self.nomAngle
		
	def setThresh(self, val):
		self.thresh = val
		
	def isStable(self):
		
		if not self.isInit:
			return False
		
		if fabs(self.currAngle - self.nomAngle) > self.thresh:
			return False
		
		return True
	
	def reset(self):
		self.currAngle = 0.0
		self.nomAngle = 0.0
		self.isInit = False
		self.count = 0

	def addData(self, newValue):
		
		self.count += 1
		
		if self.count < 50:
			return False
		
		if not self.isInit:
			self.nomAngle = newValue
			self.isInit = True
					
		self.currAngle = newValue	
		return self.isStable()
"""	
	
class DiffStability:

	def __init__(self, thresh, sampleCount):
		self.thresh = thresh
		self.sampleCount = sampleCount

		self.currAngle = 0.0
		self.nomAngle = 0.0
		self.isInit = False
		self.count = 0
		
		"""
		1. set nominal angle
		2. for each iteration, check if current angle violates difference threshold
		3. otherwise set as stable
		4. reset
		"""	
		
	def getDiff(self):
		return self.currAngle - self.nomAngle
		
	def setThresh(self, val):
		self.thresh = val
		
	def isStable(self):
		
		if not self.isInit:
			return False
		
		if fabs(self.currAngle - self.nomAngle) > self.thresh:
			return False
		
		return True

	def initialize(self, newValue):
		self.currAngle = newValue
		self.nomAngle = newValue
		self.isInit = True
		self.count = self.sampleCount+1

	
	def reset(self):
		self.currAngle = 0.0
		self.nomAngle = 0.0
		self.isInit = False
		self.count = 0

	def addData(self, newValue):
		
		self.count += 1
		
		if self.count < self.sampleCount:
			return False
		
		if not self.isInit:
			self.nomAngle = newValue
			self.isInit = True
					
		self.currAngle = newValue	
		return self.isStable()

for k in range(1,15):
#for k in range(2,3):
	pylab.clf()

	angleFile = open("angle_output_%04u.txt" % k, 'r')
	driftFile = open("drift_output_%04u.txt" % k, 'w')

	angles = [[] for i in range(39)]

	for line in angleFile:
		angleSnapshot = eval(line)
		for i in range(len(angleSnapshot)):
			angles[i].append(angleSnapshot[i])

		#angles.append(angleSnapshot)
	
	diffStabs = [DiffStability(0.01, 1) for i in range(39)]
	isStables = [[] for i in range(39)]
	diffYPs = [[] for i in range(39)]
	threshYPs = [[] for i in range(39)]
	numSamples = len(angles[0])
	stableIndices = []

	currStable = [False for i in range(39)]
	prevStable = [False for i in range(39)]
	nomAngles = [0.0 for i in range(39)]
	maxAngles = [0.0 for i in range(39)]
	maxAngleHistory = [[] for i in range(39)]
	stableCount = [0 for i in range(39)]
	entries= []

	xP = []
	for i in range(numSamples):

		if i % 1 == 0:
			xP.append(i/100.0)

			for j in range(39):
				angle = angles[j][i]
				diffStabs[j].addData(angle)

				if diffStabs[j].isStable():
					isStables[j].append(True)
					threshYPs[j].append(1.0)

					prevStable[j] = currStable[j]
					currStable[j] = True
					
				else:
					isStables[j].append(False)
					threshYPs[j].append(0.0)

					prevStable[j] = currStable[j]
					currStable[j] = False
					
				diffYPs[j].append(diffStabs[j].getDiff())

				if not prevStable[j] and currStable[j]:
					nomAngles[j] = angle
					maxAngles[j] = angle
					maxAngleHistory[j] = [diffStabs[j].nomAngle]
					#maxAngleHistory[j].append(maxAngles[j])
					stableCount[j] = 1
					
				if prevStable[j] and not currStable[j]:
					if stableCount[j] > 1:
						#entries.append([stableCount[j], fabs(nomAngles[j] - maxAngles[j]), j, nomAngles[j]])
						entries.append([stableCount[j], fabs(nomAngles[j] - maxAngleHistory[j][-2]), j, nomAngles[j]])
					diffStabs[j].reset()
				
				if prevStable[j] and currStable[j]:
					stableCount[j] += 1
					diff = fabs(nomAngles[j] - angle)
					maxDiff = fabs(nomAngles[j] - maxAngles[j])
					
					if diff > maxDiff:
						maxAngles[j] = angle

					maxAngleHistory[j].append(maxAngles[j])
					
			" now find the largest number of consecutive stable joints "
			lowIndex = 0
			highIndex = 0
			maxWidth = 0
			
			currLow = 0
			currHigh = 0
			currWidth = 0
			isContiguous = False
			
			for j in range(39):
				
				if isStables[j][-1]:
					if isContiguous:
						currHigh = j
						currWidth = currHigh - currLow
						
						if currWidth > maxWidth:
							maxWidth = currWidth
							highIndex = currHigh
							lowIndex = currLow
							
					else:
						currLow = j
						currHigh = currLow
						isContiguous = True
						
				else:
					isContiguous = False
					
			stableIndices.append([lowIndex,highIndex])

	" sort the entries "
	def cmp(list1, list2):
		#print list2[1] >list1[1]
		return int(1e6 * (list2[1] - list1[1]))
		if list2[1] < list1[1]:
			return True
		return False
		
	
	entries.sort(cmp)
	for entry in entries:
		driftFile.write(repr(entry))
		driftFile.write("\n")
		
	for i in range(len(stableIndices)):
		#ind in stableIndices:
		low = stableIndices[i][0]
		high = stableIndices[i][1]

		xP = [i/100.0, i/100.0]
		yP = [low, high]

		pylab.plot(xP,yP, color='0.5')

	xBlue = []
	yBlue = []
	xRed = []
	yRed = []

	for i in range(len(isStables[0])):

		for j in range(len(isStables)):
			if isStables[j][i]:
				xBlue.append(i/100.0)
				yBlue.append(j)
			else:
				xRed.append(i/100.0)
				yRed.append(j)

	maxX = len(isStables[0])/100.0

	pylab.title("Stable Sections, Local Pose Joints %d" % k)
	#pylab.scatter(xBlue,yBlue, color='b', linewidth=1)
	pylab.scatter(xRed,yRed, color='r', linewidth=1)
	pylab.xlim(0,maxX)
	pylab.ylim(0,39)

	#pylab.show()
	pylab.savefig("stableSection_%04u.png" % k)

	"""
	yP = []
	varYP = []
	threshYP = []
	for i in range(len(angles[14])):
		angle = angles[14][i]
		xP.append(i/100.0)
		yP.append(angle)
		valStab.addData(angle)

		if valStab.sampleCount >= valStab.sampleSize:
			varYP.append(valStab.getVar())

		else:
			varYP.append(0.0)

		if valStab.isStable():
			threshYP.append(1.0)
		else:
			threshYP.append(0.0)


	pylab.plot(xP,yP)
	pylab.plot(xP,varYP)
	pylab.plot(xP,threshYP)

	#pylab.show()

	xlabel('time (s)')

	#savefig("jointsHistory_%04u.png" % j)
	show()
	"""
