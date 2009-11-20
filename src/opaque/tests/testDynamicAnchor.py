#!/usr/bin/python
import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from numpy import arange
from math import *
from copy import *
import pylab
from random import *

"""
Algorithm for dynamic adaptation anchoring to unknown wall contours with varying widths.

1. Generate a sinusoidal curve starting at tip.
2. Build point-modifiable amplitude mask.
3. Increase amplitude of two anchor points.
4. Detect joint error to see if anchors have made contact.

"""

class AdaptiveCosine:

	def __init__(self, freq, length):

		# cosine x length for one cycle
		self.x_len = length
		self.freq = freq
		self.amp = 1.0
		self.origin = [0.0,0.0]
		self.oVec = [1.0,0.0]

		" amplitude control points, spaced at every crest of the sinusoid "
		self.ampPoints = [0.0 for i in range(40)]
		self.ampPoints[0] = 1.0
		
		self.ampSpacing = 2*pi / self.freq / 2.0

		p1 = self.origin
		p2 = [p1[0] + self.oVec[0]*self.x_len, p1[1] + self.oVec[1]*self.x_len]

		self.control = {'p1':p1, 'p2':p2}

	def setAmpPoint(self, index, amplitude):
		self.ampPoints[index] = amplitude
		
	def getAmpPoint(self, index, amplitude):
		return self.ampPoints[index]

	def draw(self):
		points = self.getPoints()	

		xP = []
		zP = []
		for p in points:
			xP.append(p[0])
			zP.append(p[1])
		pylab.plot(xP,zP)
	
		xP = []
		zP = []		
		for i in range(len(self.ampPoints)):
			xP.append(i*self.ampSpacing)
			zP.append(self.ampPoints[i])
		pylab.plot(xP,zP)
		
		pylab.xlim(-1,5)
		pylab.ylim(-4,2)
		pylab.show()

	def findClosestPoint(self, point):

		# 3 hypotheses: closest point on the trailing edge, the leading edge, or the cosine curve

		# 1. closest point simple if between boundaries of trailing and leading edges
		# assuming oVec = [1,0] and origin = [0.0,0.0]

		p1 = self.control['p1']
		p2 = self.control['p2']

		# Hypothesis 1: closest point on the cosine curve
		# evenly sample the cosine curve to find the closest point
		samples = arange(0.0,self.x_len + 0.01, 0.01)

		min3 = 100000
		curvePoint3 = [0.0,0.0]

		for x_samp in samples:
			z = self.amp * cos(self.freq*x_samp) - self.amp

			dist = sqrt((point[0]-(x_samp))**2 + (point[1]-z)**2)

			if dist < min3:
				min3 = dist
				curvePoint3 = [x_samp,z]

		min = 0.0
		pnt = []

		min = min3
		pnt = curvePoint3

		xP = [pnt[0],point[0]]

		return min, 0.5, copy(pnt)

	def getPoints(self):
		# create a distribution of points along the curve

		points = []

		samples = arange(0.0,self.x_len + 0.01, 0.01)
		for x_samp in samples:
			
			" determine the amplitude control points we'll use "
			ampIndex = int(x_samp / self.ampSpacing)
			
			" point in space between the control points "
			xSet = x_samp % self.ampSpacing
			
			amp0 = self.ampPoints[ampIndex]
			amp1 = self.ampPoints[ampIndex+1]
			
			xVal = xSet / self.ampSpacing
			
			adapAmp = xVal * (amp1-amp0) + amp0

			#z = self.amp * cos(self.freq*x_samp) - self.amp	
			z = adapAmp * cos(self.freq*x_samp) - self.amp
			points.append([x_samp,z])

		return points


if __name__ == "__main__":
	" frequency and length "
	ac = AdaptiveCosine(4*pi, 5.0)
	ac.draw()


