
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

		self.ampSpacing = 2*pi / self.initFreq / 2.0
		self.numPeaks = 2


		dispZ = self.peakAmp[0]			

		pO = self.origin

		pA_1 = [pO[0] + self.ampSpacing, pO[1] - 2*dispZ]
		pB_1 = [pO[0] - self.ampSpacing, pO[1] - 2*dispZ]
		
		pA_2 = [pA_1[0] + self.infLength, pA_1[1]]
		pB_2 = [pB_1[0] - self.infLength, pB_1[1]]
		
		self.controlA = {'pO':pO, 'p1':pA_1, 'p2':pA_2}
		self.controlB = {'pO':pO, 'p1':pB_1, 'p2':pB_2}

	def setPeakAmp(self, amplitude):
		
		self.peakAmp[0] = amplitude
	
		pO = self.controlA['pO']
		
		dispZ = self.peakAmp[0]			

		pA_1 = [pO[0] + self.ampSpacing, pO[1] - 2*dispZ]
		pB_1 = [pO[0] - self.ampSpacing, pO[1] - 2*dispZ]
		
		pA_2 = [pA_1[0] + self.infLength, pA_1[1]]
		pB_2 = [pB_1[0] - self.infLength, pB_1[1]]

		self.controlA['p1'] = pA_1
		self.controlA['p2'] = pA_2
		
		self.controlB['p1'] = pB_1
		self.controlB['p2'] = pB_2
	
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


	def getPoints(self):

		" create a distribution of points along the curve "

		pO = self.controlA['pO']		
		pA_1 = self.controlA['p1']
		pA_2 = self.controlA['p2']
		pB_1 = self.controlB['p1']
		pB_2 = self.controlB['p2']
		
		points = []

		" infinite lead "
		samplesB = arange(pB_2[0], pB_1[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samplesB:
			dispZ = -2*self.peakAmp[0]
			points.append([x_samp,dispZ])


		" anchor section "
		samples1 = arange(pB_1[0], pA_1[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samples1:
			adapAmp = self.peakAmp[0]
			varZ = adapAmp * cos(self.initFreq*x_samp) - self.peakAmp[0]

			points.append([x_samp,varZ])

		" infinite follow "
		samplesA = arange(pA_1[0], pA_2[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samplesA:
			dispZ = -2*self.peakAmp[0]
			points.append([x_samp,dispZ])
		
		
		return points


	def computeX(self, point, xNotLessThan):
		
		if point[0] < xNotLessThan:
			return 1e100, 0.5, copy(point)
		
		" Hypothesis 1: closest point on the cosine curve "
		" evenly sample the cosine curve to find the closest point "
				
		x_samp = point[0]

		if x_samp < self.head_len:
			return [x_samp, 0.0]

		" determine the amplitude control points we'll use "
		if (x_samp-self.head_len) < 3.0*self.ampSpacing/2.0:
			ampIndex = 0
		else:
			ampIndex = int(((x_samp-self.head_len) - self.ampSpacing/2.0) / self.ampSpacing)

		for i in range(len(self.infPoints)):
			infX = self.infPoints[i]
			if x_samp < infX:
				ampIndex = i
				break
			
		adapAmp = self.peakAmp[ampIndex]
		
		# FIXME: uncommenting this will cause an error, should have no effect 
		
		#print "ampIndex =", ampIndex

		if ampIndex == 0:
			varZ = adapAmp * cos(self.initFreq*x_samp) - self.peakAmp[0]
		elif ampIndex >= self.numPeaks:
			varZ = -self.peakAmp[0]
		elif ampIndex % 2 == 1:
			xdiff = x_samp - self.infPoints[ampIndex-1]
			peakWidth = self.peakSpacings[ampIndex]
			freq = 2*pi / peakWidth / 2.0
			varZ = adapAmp * sin(freq*xdiff) - self.peakAmp[0]
		else:
			xdiff = x_samp - self.infPoints[ampIndex-1]
			peakWidth = self.peakSpacings[ampIndex]
			freq = 2*pi / peakWidth / 2.0
			varZ = adapAmp * sin(pi + freq*xdiff) - self.peakAmp[0]

		#if ampIndex < 1:
		#print ampIndex, x_samp, varZ
			
		#varZ = adapAmp * cos(self.freq*x_samp) - self.peakAmp[0]

		#if ampIndex < 1:
		#print ampIndex, x_samp, varZ
		#print
		
		dist = sqrt((point[0]-(x_samp))**2 + (point[1]-varZ)**2)
		curvePoint1 = [x_samp,varZ]
		
		" if we go beyond the curve, raise an exception "
		p3 = self.control['p3']
		if curvePoint1[0] > p3[0]:
			raise
		
		return curvePoint1

	def findClosestPoint(self, point, xNotLessThan = 0.0):

		# 3 hypotheses: closest point on the trailing edge, the leading edge, or the cosine curve

		# 1. closest point simple if between boundaries of trailing and leading edges
		# assuming oVec = [1,0] and origin = [0.0,0.0]

		if point[0] < xNotLessThan:
			return 1e100, 0.5, copy(point)		

		p1 = self.control['p1']
		pA = self.control['pA']
		p2 = self.control['p2']
		p3 = self.control['p3']

		#print self.control

		# Hypothesis 1: closest point on the cosine curve
		# evenly sample the cosine curve to find the closest point
		#newSamples = arange(0.0,self.x_len + self.sampleResolution, self.sampleResolution)
		newSamples = arange(self.head_len,self.x_len + self.sampleResolution, self.sampleResolution)

		samples = []
		for i in range(len(newSamples)):
			if newSamples[i] >= xNotLessThan:
				samples.append(newSamples[i])

		min1 = 100000
		curvePoint1 = [0.0,0.0]
		curvePoint2 = [0.0,0.0]
		curvePoint3 = [0.0,0.0]
		termP = []

		for x_samp in samples:

			" determine the amplitude control points we'll use "
			if (x_samp - self.head_len) < 3.0*self.ampSpacing/2.0:
				ampIndex = 0
			else:
				ampIndex = int(((x_samp - self.head_len) - self.ampSpacing/2.0) / self.ampSpacing)

			for i in range(len(self.infPoints)):
				infX = self.infPoints[i]
				if x_samp < infX:
					ampIndex = i
					break
						
			adapAmp = self.peakAmp[ampIndex]
			
			
			if ampIndex == 0:
				varZ = adapAmp * cos(self.initFreq*x_samp) - self.peakAmp[0]
			elif ampIndex >= self.numPeaks:
				varZ = -self.peakAmp[0]
			elif ampIndex % 2 == 1:
				xdiff = x_samp - self.infPoints[ampIndex-1]
				peakWidth = self.peakSpacings[ampIndex]
				freq = 2*pi / peakWidth / 2.0
				varZ = adapAmp * sin(freq*xdiff) - self.peakAmp[0]
			else:
				xdiff = x_samp - self.infPoints[ampIndex-1]
				peakWidth = self.peakSpacings[ampIndex]
				freq = 2*pi / peakWidth / 2.0
				varZ = adapAmp * sin(pi + freq*xdiff) - self.peakAmp[0]
			
			#varZ = adapAmp * cos(self.freq*x_samp) - self.peakAmp[0]

			
			dist = sqrt((point[0]-(x_samp))**2 + (point[1]-varZ)**2)

			if dist < min1:
				min1 = dist
				curvePoint1 = [x_samp,varZ]
			
			if x_samp == samples[-1]:
				termP = [x_samp,varZ]

		# Hypothesis 2: closest point on the trailing edge
		if point[0] <= p2[0]:
			min2 = sqrt((point[0]-p2[0])**2 + (point[1]-p2[1])**2)
			curvePoint2 = p2

		elif point[0] > p2[0] and point[0] < p3[0]:
			min2 = fabs(point[1] - p2[1])
			#min2 = fabs(point[1])
			curvePoint2 = [point[0],p2[1]]

		elif point[0] >= p3[0]:
			min2 = sqrt((point[0]-p3[0])**2 + (point[1]-p3[1])**2)
			curvePoint2 = p3

		else:
			print "points not in expected configuration"
			raise
		

		# Hypothesis 3: closest point on the leading edge
		if point[0] <= p1[0]:
			min3 = sqrt((point[0]-p1[0])**2 + (point[1]-p1[1])**2)
			curvePoint3 = p1

		elif point[0] > p1[0] and point[0] < pA[0]:
			min3 = fabs(point[1] - p1[1])
			curvePoint3 = [point[0],p1[1]]

		elif point[0] >= pA[0]:
			min3 = sqrt((point[0]-pA[0])**2 + (point[1]-pA[1])**2)
			curvePoint3 = pA

		else:
			print "points not in expected configuration"
			raise
		

		min = 0.0
		pnt = []


		if min1 <= min2:

			if min1 <= min3:
				min = min1
				pnt = curvePoint1

			else:
				min = min3
				pnt = curvePoint3

		else:
			if min2 <= min3:
				min = min2
				pnt = curvePoint2
			else:
				min = min3
				pnt = curvePoint3

		#if pnt[0] == termP[0] and pnt[1] == termP[1]:
		#	return min, 1.0, copy(pnt)
		#else:
		#	return min, 0.5, copy(pnt)
		
		#print pnt, p3
		#print "minimum =", min
		if pnt[0] == p1[0] and pnt[1] == p1[1]:
			return min, 0.0, copy(pnt)
		elif pnt[0] == p3[0] and pnt[1] == p3[1]:
			return min, 1.0, copy(pnt)
		else:
			return min, 0.5, copy(pnt)
		

