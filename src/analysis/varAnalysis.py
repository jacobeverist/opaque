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

class ValueStability:

	def __init__(self, thresh, sample_size = 20):
		self.sampleMean = 0.0
		self.sampleVar = 0.0
		self.sampleCount = 0
		self.sampleIndex = 0
		self.varThresh = thresh
		
		self.sampleSize = sample_size
		
		self.samples=[]
		for i in range(self.sampleSize):
			self.samples.append(0.0)
		
			
	def setThresh(self, val):
		self.varThresh = val
		
	def isStable(self):
		if self.sampleCount < self.sampleSize:
			return False

		if self.sampleVar < self.varThresh:
			return True
		
		return False

	def getMean(self):
		return self.sampleMean

	def getVar(self):
		return self.sampleVar
	
	def getSampleCount(self):
		return self.sampleCount

	def reset(self):
		self.sampleCount = 0
		self.sampleIndex = 0
		self.sampleMean = 0.0
		self.sampleVar = 0.0

	def addData(self, newValue):

		# **** COMPUTE SAMPLE VARIANCE AND MEAN  ****
		#	
		# NOTE: Samples are taken at every time step and are used to determine
		# the stability of the probe position.  Low variance indicates that the
		# position is stable and it is safe to compute a feature point.
		#

		# retrieve data sample
		self.samples[self.sampleIndex] = newValue

		# increment the circular index
		self.sampleIndex += 1
		self.sampleIndex = self.sampleIndex % self.sampleSize

		if self.sampleCount < self.sampleSize:
			self.sampleCount += 1

		# compute the mean
		self.sampleSum = 0.0
		for i in range(self.sampleCount):
			self.sampleSum += self.samples[i]
		self.sampleMean = self.sampleSum / self.sampleCount 

		# compute the variance
		self.sampleSum = 0
		for i in range(self.sampleCount):
			self.sampleSum += (self.samples[i] - self.sampleMean)*(self.samples[i] - self.sampleMean)

		self.sampleVar = self.sampleSum/self.sampleCount

		max = -1e200
		min = 1e200
		for i in range(0,self.sampleCount):
			if self.samples[i] > max:
				max = self.samples[i]
			if self.samples[i] < min:
				min = self.samples[i]

		# should have at least 10 values before we compute the variance
		if self.sampleCount < self.sampleSize :
			# make extremely high for invalid variances
			self.sampleVar = 1e200;

		if self.sampleVar < self.varThresh:
			return True

		return False



for k in range(1,15):
#for k in range(1,2):
	pylab.clf()

	angleFile = open("angle_output_%04u.txt" % k, 'r')
	driftFile = open("drift_output_%04u.txt" % k, 'w')

	angles = [[] for i in range(39)]

	for line in angleFile:
		angleSnapshot = eval(line)
		for i in range(len(angleSnapshot)):
			angles[i].append(angleSnapshot[i])

		#angles.append(angleSnapshot)
	
	valStabs = [ValueStability(1e-3, 50) for i in range(39)]
	isStables = [[] for i in range(39)]
	varYPs = [[] for i in range(39)]
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
				valStabs[j].addData(angle)

				if valStabs[j].sampleCount >= valStabs[j].sampleSize:
					if valStabs[j].isStable():
						isStables[j].append(True)
						threshYPs[j].append(1.0)

						prevStable[j] = currStable[j]
						currStable[j] = True
						
					else:
						isStables[j].append(False)
						threshYPs[j].append(0.0)

						prevStable[j] = currStable[j]
						currStable[j] = False

				else:
					isStables[j].append(False)
					threshYPs[j].append(0.0)

					prevStable[j] = currStable[j]
					currStable[j] = False

			
				varYPs[j].append(valStabs[j].getVar())

				if not prevStable[j] and currStable[j]:
					nomAngles[j] = angle
					maxAngles[j] = angle
					maxAngleHistory[j].append(maxAngles[j])
					stableCount[j] = 1
					
				if prevStable[j] and not currStable[j]:
					
					if stableCount[j] > 50:
						#entries.append([stableCount[j], fabs(nomAngles[j] - maxAngles[j]), j, nomAngles[j]])
						entries.append([stableCount[j], fabs(nomAngles[j] - maxAngleHistory[j][-50]), j, nomAngles[j]])
					
					
				
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
