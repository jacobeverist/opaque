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

def normalizeAngle(angle):

	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle 


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

for j in range(1,15):
#for j in range(1,15):
#for j in range(10,11):

	print j

	pylab.clf()

	angleFile = open("global_pose_output_%04u.txt" % j, 'r')

	angles = [[] for i in range(40)]

	initNom = False
	nomAng = [0.0 for i in range(40)]


	for line in angleFile:

		angleSnapshot = eval(line)

		if not initNom and len(angles[0]) >= 50:

			#print len(angleSnapshot)

			for i in range(len(angleSnapshot)):
				#print angleSnapshot[i]
				nomAng[i] = normalizeAngle(angleSnapshot[i][2])

			initNom = True

		if initNom:
			for i in range(len(angleSnapshot)):
				angles[i].append(normalizeAngle(normalizeAngle(angleSnapshot[i][2])-nomAng[i]))
		else:
			for i in range(len(angleSnapshot)):
				angles[i].append(0.0)


	data = numpy.array(angles)

	numSamples, numRows = len(angles[0]),len(angles)

	#print "numSamples =", numSamples

	data = data.transpose()

	t = 0.01 * numpy.arange(numSamples, dtype=float)
	ticklocs = []
	ax = subplot(111)
	#xlim(0,t[-1])
	xlim(0,t[1000])
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

	# set the yticks to use axes coords on the y axis
	ax.set_yticks(ticklocs)
	labels = [str(i) for i in range(40)]
	ax.set_yticklabels(labels)

	xlabel('time (s)')
	title("Global Pose Joints %d" % j)

	savefig("globalJointsHistory_%04u.png" % j)


	pylab.clf()

	for i in range(40):
		xP = []
		yP = []
		for k in range(numSamples):
			xP=[k*0.01]
			yP=[i]
			val = fabs(angles[i][k])
			if val > 1.0:
				val = 1.0
			dotColor = (1.0-val, 1.0, 1.0)
			pylab.scatter(xP,yP, color=dotColor)

	xlabel('time (s)')
	title("Global Pose Error %d" % j)
	savefig("globalAngleError_%04u.png" % j)


	"""
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
	title("Global Pose Stability %d" % j)

	savefig("globalStabHistory_%04u.png" % j)
	"""

	#show()


