#!/usr/bin/python


import os
import csv
#import pylab
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

for j in range(1,3):
#for j in range(1,15):
#for j in range(10,11):

	pylab.clf()

	angleFile = open("angle_output_%04u.txt" % j, 'r')

	angles = [[] for i in range(39)]
	valStabs = [ValueStability(1e-2, 50) for i in range(39)]
	varYPs = [[] for i in range(39)]
	threshYPs = [[] for i in range(39)]

	for line in angleFile:
		angleSnapshot = eval(line)
		for i in range(len(angleSnapshot)):
			angles[i].append(angleSnapshot[i])
			valStabs[i].addData(angleSnapshot[i])

			if valStabs[i].sampleCount >= valStabs[i].sampleSize:
				varYPs[i].append(valStabs[i].getVar())

			else:
				varYPs[i].append(0.0)

			if valStabs[i].isStable():
				threshYPs[i].append(0.0)
			else:
				threshYPs[i].append(1.0)

		#angles.append(angleSnapshot)

	data = numpy.array(angles)
	data2 = numpy.array(varYPs)
	data3 = numpy.array(threshYPs)

	#xP = []
	#yP = []
	#for i in range(len(angles)):
	#	angle = angles[i][0]
	#	xP.append(i/100.0)
	#	yP.append(angle)

	#pylab.plot(xP,yP)

	#pylab.show()

	numSamples, numRows = len(angles[0]),len(angles)

	data = data.transpose()
	data2 = data2.transpose()
	data3 = data3.transpose()


	t = 0.01 * numpy.arange(numSamples, dtype=float)
	ticklocs = []
	ax = subplot(111)
	xlim(0,t[-1])
	xticks(numpy.arange(0,100,10))

	dmin = data.min()
	dmax = data.max()

	dr = (dmax - dmin)*0.7 # Crowd them a bit.
	y0 = dmin
	y1 = (numRows-1) * dr + dmax
	ylim(y0, y1)

	segs = []
	for i in range(numRows):
		segs.append(numpy.hstack((t[:,numpy.newaxis], data[:,i,numpy.newaxis])))
		ticklocs.append(i*dr)

	offsets = numpy.zeros((numRows,2), dtype=float)
	offsets[:,1] = ticklocs

	lines = LineCollection(segs, offsets=offsets,
												 transOffset=None,
												 )

	ax.add_collection(lines)

	segs2 = []
	for i in range(numRows):
		segs2.append(numpy.hstack((t[:,numpy.newaxis], data2[:,i,numpy.newaxis])))
	lines2 = LineCollection(segs2, offsets=offsets,
												 transOffset=None,
												 )
	#ax.add_collection(lines2)

	# set the yticks to use axes coords on the y axis
	ax.set_yticks(ticklocs)
	labels = [str(i) for i in range(39)]
	ax.set_yticklabels(labels)

	xlabel('time (s)')
	title("Local Pose Joints %d" % j)

	savefig("jointsHistory_%04u.png" % j)


	ax.cla()

	segs3 = []
	for i in range(numRows):
		segs3.append(numpy.hstack((t[:,numpy.newaxis], data3[:,i,numpy.newaxis])))

	lines3 = LineCollection(segs3, offsets=offsets,
												 transOffset=None,
												 )
	ax.add_collection(lines3)


	# set the yticks to use axes coords on the y axis
	ax.set_yticks(ticklocs)
	labels = [str(i) for i in range(39)]
	ax.set_yticklabels(labels)

	xlabel('time (s)')
	title("Local Pose Stability %d" % j)

	savefig("stabHistory_%04u.png" % j)

	#show()


