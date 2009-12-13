
from numpy import arange
from math import *
from copy import *
import pylab

class AdaptiveAnchorCurve:

	def __init__(self, freq):

		" the curve is a cosine curve centered at the peak and "
		" straight lines to positive and negative infinity "

		self.initFreq = freq
		self.freq = freq
		self.amp = 1.0
		self.origin = [0.0,0.0]
		self.sampleResolution = 0.01

		self.infLength = 5.0

		self.peakAmp = [0.0 for i in range(40)]

		self.ampSpacing = 5*pi / self.initFreq / 2.0
		self.numPeaks = 1


		dispZ = self.peakAmp[0]			

		pO = self.origin

		p1 = [pO[0] + self.ampSpacing, pO[1] + dispZ]		
		p2 = [p1[0] + self.infLength, p1[1]]		

		self.controlA = {'pO':pO, 'p1':p1, 'p2':p2}

	def getPeakWidth(self):
		return self.ampSpacing
	
	def setPeakAmp(self, amplitude):
		
		self.peakAmp[0] = amplitude
	
		pO = self.controlA['pO']
		
		dispZ = self.peakAmp[0]			

		p1 = [pO[0] + self.ampSpacing, pO[1] + dispZ]		
		p2 = [p1[0] + self.infLength, p1[1]]		

		self.controlA['p1'] = p1
		self.controlA['p2'] = p2
	
	def getPeakAmp(self):
		return self.peakAmp[0]

	def draw(self):
		points = self.getPoints()	

		xP = []
		zP = []
		for p in points:
			xP.append(p[0])
			zP.append(p[1])
		pylab.plot(xP,zP)
	
		pylab.xlim(-1,5)
		pylab.ylim(-4,2)
		pylab.show()

	def isOutOfSegments(self, tipPoint):
		p1 = self.controlA['p1']
		
		#print "comparing tipPoint =", tipPoint, "to inflection point", p1
		if tipPoint[0] < p1[0]:
			return True
		
		return False

	def getPoints(self):

		" create a distribution of points along the curve "

		pO = self.controlA['pO']		
		p1 = self.controlA['p1']
		p2 = self.controlA['p2']
		
		points = []

		" anchor section "
		samples1 = arange(pO[0], p1[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samples1:
			varZ = self.peakAmp[0] * sin(self.initFreq*x_samp)
			points.append([x_samp,varZ])

		" infinite follow "
		samples2 = arange(p1[0], p2[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samples2:
			points.append([x_samp,p2[1]])
		
		return points


	def computeX(self, point, xNotLessThan):
		
		if point[0] < xNotLessThan:
			return 1e100, 0.5, copy(point)
		
		" Hypothesis 1: closest point on the cosine curve "
		" evenly sample the cosine curve to find the closest point "
				
		x_samp = point[0]

		" determine the amplitude control points we'll use "
		if x_samp >= self.ampSpacing:
			p2 = self.controlA['p2']
			return [x_samp,p2[1]]
					
		adapAmp = self.peakAmp[0]
		
		varZ = adapAmp * sin(self.initFreq*x_samp)
		
		
		curvePoint1 = [x_samp,varZ]
		
		return curvePoint1

	def findClosestPoint(self, point, xNotLessThan = 0.0):

		# 3 hypotheses: closest point on the trailing edge, the leading edge, or the cosine curve

		# 1. closest point simple if between boundaries of trailing and leading edges
		# assuming oVec = [1,0] and origin = [0.0,0.0]

		if point[0] < xNotLessThan:
			return 1e100, 0.5, copy(point)		

		pO = self.controlA['pO']		
		p1 = self.controlA['p1']
		p2 = self.controlA['p2']
		
		# Hypothesis 1: closest point on the cosine curve
		# evenly sample the cosine curve to find the closest point
		newSamples = arange(0,self.ampSpacing + self.sampleResolution, self.sampleResolution)

		samples = []
		for i in range(len(newSamples)):
			if newSamples[i] >= xNotLessThan:
				samples.append(newSamples[i])

		min1 = 100000
		curvePoint1 = [0.0,0.0]
		curvePoint2 = [0.0,0.0]

		for x_samp in samples:

			" determine the amplitude control points we'll use "
			if x_samp < self.ampSpacing:
				ampIndex = 0
						
			adapAmp = self.peakAmp[0]
						
			varZ = adapAmp * sin(self.initFreq*x_samp)
			
			dist = sqrt((point[0]-(x_samp))**2 + (point[1]-varZ)**2)

			if dist < min1:
				min1 = dist
				curvePoint1 = [x_samp,varZ]

		# Hypothesis 2: closest point on the trailing edge
		if point[0] <= p1[0]:
			min2 = sqrt((point[0]-p1[0])**2 + (point[1]-p1[1])**2)
			curvePoint2 = p1

		elif point[0] > p1[0] and point[0] < p2[0]:
			min2 = fabs(point[1] - p1[1])
			curvePoint2 = [point[0],p1[1]]

		elif point[0] >= p2[0]:
			min2 = sqrt((point[0]-p2[0])**2 + (point[1]-p2[1])**2)
			curvePoint2 = p2

		else:
			print "points not in expected configuration"
			raise
		
		min = 0.0
		pnt = []

		if min1 <= min2:
			min = min1
			pnt = curvePoint1

		else:
			min = min2
			pnt = curvePoint2


		
		if pnt[0] == pO[0] and pnt[1] == pO[1]:
			return min, 0.0, copy(pnt)
		elif pnt[0] == p2[0] and pnt[1] == p2[1]:
			return min, 1.0, copy(pnt)
		else:
			return min, 0.5, copy(pnt)
		

