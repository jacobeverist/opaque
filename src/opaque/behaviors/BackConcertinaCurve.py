
from numpy import arange
from math import *
from copy import *
import pylab

"""
Class for dynamic adaptation anchoring to unknown wall contours with varying widths.

1. Generate a sinusoidal curve starting at tip.
2. Build point-modifiable amplitude mask.
3. Increase amplitude of two anchor points.
4. Detect joint error to see if anchors have made contact.

"""

"""
Stage 2 Adaptive Cosine

1.  Be able to channge the number of peaks.  Start off with 3.
2.  Two sections of curve, the sinusoid and the tail
3.  The tail stays constant and eats up the anchor as it demands more segments

"""

"""
Stage 3

1.  New cosine curve starts at the beginning of each peak division
2.  The first peak includes 1.5*pi of a curve
3.  Subsequent peaks contain only pi
4.  peak 1:  starts at cos(0.0)
    even peaks:  starts at sin(0.0)
    odd peaks:  starts at sin(pi)
    
5.  frequency changes to fit the prescribed width
		freq = 2*pi / peakWidth / 2.0

"""

class BackConcertinaCurve:

	def __init__(self, freq):

		# cosine x length for one cycle
		self.initFreq = freq
		self.amp = 1.0
		self.origin = [0.0,0.0]
		self.oVec = [1.0,0.0]
		self.sampleResolution = 0.01
		self.head_len = 0.0


		self.peakAmp = [0.0 for i in range(40)]
	
		self.ampSpacing = 2*pi / self.initFreq / 2.0
		#self.ampSpacing *= 2.0
		self.numPeaks = 2

		self.peakSpacings = [self.ampSpacing for i in range(40)]

		#self.peakSpacings[0] = 3.0 * self.ampSpacing/2.0
		#self.peakSpacings[3] = self.ampSpacing * 2.0
		#self.peakSpacings[4] = self.ampSpacing * 2.0
		#self.peakSpacings[5] = self.ampSpacing * 2.0

		" peak inflection points, spaced at every cross over the x-axis "
		self.infPoints = []
		for i in range(self.numPeaks):
			totalSum = 0.0
			totalSum += self.head_len

			for j in range(i):
				totalSum += self.peakSpacings[j]

			totalSum += self.peakSpacings[i]
			
			self.infPoints.append(totalSum)	

		p1 = self.origin

		pA = [self.origin[0] + self.head_len, self.origin[1]]
		

		self.x_len = self.infPoints[-1]
		
		#varZ = -self.peakAmp[0]			
		varZ = 0.0			

		p2 = [self.x_len, pA[1] + varZ]

		
		#self.tail_len = 3.0
		self.tail_len = 0.0
		
		followLen = self.tail_len - self.head_len
		if followLen < 0:
			followLen = 0.0

		p3 = [p2[0] + followLen, p2[1]]
		
		self.control = {'p1':p1, 'pA':pA, 'p2':p2, 'p3':p3}

		self.saveCount = 0
		
	def saveAmps(self):
		f = open("amps%04u.txt" % self.saveCount, 'w')
		f.write(repr(self.peakAmp))
		f.close()
		self.saveCount += 1

	def getHeadLength(self):
		return self.head_len

	def setHeadLength(self, length):
		self.head_len = length

		if True:

			" peak inflection points, spaced at every cross over the x-axis "
			self.infPoints = []
			for i in range(self.numPeaks):
				totalSum = 0.0
				totalSum += self.head_len

				for j in range(i):
					totalSum += self.peakSpacings[j]

				totalSum += self.peakSpacings[i]
				
				self.infPoints.append(totalSum)	

			p1 = self.origin
			pA = [self.origin[0] + self.head_len, self.origin[1]]
			self.x_len = self.infPoints[-1]
			#varZ = -self.peakAmp[0]			
			varZ = 0.0			
			p2 = [self.x_len, pA[1] + varZ]
	
			followLen = self.tail_len - self.head_len
			if followLen < 0:
				followLen = 0.0

			p3 = [p2[0] + followLen, p2[1]]

			#print self.control

			self.control['p1'] = p1
			self.control['pA'] = pA
			self.control['p2'] = p2
			self.control['p3'] = p3

			#print self.control

	def setTailLength(self, length):
		
		self.tail_len = length

		if True:
			p1 = self.control['p1']
			pA = self.control['pA']
			p2 = self.control['p2']
			p3 = self.control['p3']
			
			self.x_len = self.infPoints[-1]
		
			#varZ = -self.peakAmp[0]			
			varZ = 0.0		
			p2 = [self.x_len, pA[1] + varZ]

			followLen = self.tail_len - self.head_len
			if followLen < 0:
				followLen = 0.0
			p3 = [p2[0] + followLen, p2[1]]
			
			self.control['p2'] = p2
			self.control['p3'] = p3
			

	def setPeakAmp(self, index, amplitude):

		if index >= self.numPeaks:
			self.numPeaks += 1

			self.infPoints = []
			for i in range(self.numPeaks):
				totalSum = 0.0
				totalSum += self.head_len
				for j in range(i):
					totalSum += self.peakSpacings[j]
	
				totalSum += self.peakSpacings[i]
				
				self.infPoints.append(totalSum)

		self.peakAmp[index] = amplitude
	
		if True:
		#if index == 0 or index >= self.numPeaks:
			p1 = self.control['p1']
			pA = self.control['pA']
			p2 = self.control['p2']
			p3 = self.control['p3']
			

			self.x_len = self.infPoints[-1]

			#varZ = -self.peakAmp[0]			
			varZ = 0.0			

			p2 = [self.x_len, pA[1] + varZ]
			
			
			followLen = self.tail_len - self.head_len
			if followLen < 0:
				followLen = 0.0

			p3 = [p2[0] + followLen, p2[1]]
			
			self.control['p2'] = p2
			self.control['p3'] = p3
		
		#print "infPoints =", self.infPoints
	
	def setPeakWidth(self, index, width):
		
		if index == 3 or index == 4:

		#if False:
			self.peakSpacings[index] = 2.5*self.ampSpacing
		else:
			self.peakSpacings[index] = width

		self.infPoints = []
		for i in range(self.numPeaks):
			totalSum = 0.0
			totalSum += self.head_len
			for j in range(i):
				totalSum += self.peakSpacings[j]

			totalSum += self.peakSpacings[i]
			
			self.infPoints.append(totalSum)

		p1 = self.control['p1']
		pA = self.control['pA']
		p2 = self.control['p2']
		p3 = self.control['p3']

		self.x_len = self.infPoints[-1]
	
		#varZ = -self.peakAmp[0]			
		varZ = 0.0

		p2 = [self.x_len, pA[1] + varZ]
		
		followLen = self.tail_len - self.head_len
		if followLen < 0:
			followLen = 0.0

		p3 = [p2[0] + followLen, p2[1]]
		
		self.control['p2'] = p2
		self.control['p3'] = p3

	def getPeakWidth(self, index):
		return self.peakSpacings[index]
	
	def getPeakAmp(self, index):
		return self.peakAmp[index]

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

	def computeX(self, point, xNotLessThan):
		
		if point[0] < xNotLessThan:
			return 1e100, 0.5, copy(point)
		
		" Hypothesis 1: closest point on the cosine curve "
		" evenly sample the cosine curve to find the closest point "
				
		x_samp = point[0]

		if x_samp < self.head_len:
			return [x_samp, 0.0]

		" determine the amplitude control points we'll use "
		ampIndex = int((x_samp-self.head_len) / self.ampSpacing)

		for i in range(len(self.infPoints)):
			infX = self.infPoints[i]
			if x_samp < infX:
				ampIndex = i
				break
			
		adapAmp = self.peakAmp[ampIndex]
		
		if ampIndex >= self.numPeaks:
			varZ = 0.0
		else:
			varZ = adapAmp * sin(self.initFreq*x_samp)
		
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
			ampIndex = int(((x_samp - self.head_len)) / self.ampSpacing)

			for i in range(len(self.infPoints)):
				infX = self.infPoints[i]
				if x_samp < infX:
					ampIndex = i
					break
						
			adapAmp = self.peakAmp[ampIndex]
			
			
			if ampIndex >= self.numPeaks:
				varZ = 0.0
			else:
				varZ = adapAmp * sin(self.initFreq*x_samp)

			
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

		if pnt[0] == p1[0] and pnt[1] == p1[1]:
			return min, 0.0, copy(pnt)
		elif pnt[0] == p3[0] and pnt[1] == p3[1]:
			return min, 1.0, copy(pnt)
		else:
			return min, 0.5, copy(pnt)
		
	def getPoints(self):

		# create a distribution of points along the curve
		p1 = self.control['p1']
		pA = self.control['pA']
		p2 = self.control['p2']
		p3 = self.control['p3']
		
		points = []

		" head section "
		samples1 = arange(p1[0], pA[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samples1:
			varZ = p1[1]
			points.append([x_samp,varZ])
		

		" curves section "
		samples = arange(self.head_len,self.x_len + self.sampleResolution, self.sampleResolution)
		for x_samp in samples:

			" determine the amplitude control points we'll use "
			ampIndex = int(((x_samp - self.head_len)) / self.ampSpacing)


			for i in range(len(self.infPoints)):
				infX = self.infPoints[i]
				if x_samp < infX:
					ampIndex = i
					break
						
			adapAmp = self.peakAmp[ampIndex]
			
			if ampIndex >= self.numPeaks:
				varZ = 0.0
			else:
				varZ = adapAmp * sin(self.initFreq*x_samp)

			"""		
			4.  peak 1:  starts at cos(0.0)
			    even peaks:  starts at sin(0.0)
			    odd peaks:  starts at sin(pi)
			    
			5.  frequency changes to fit the prescribed width
				freq = 2*pi / peakWidth / 2.0
			"""
			
			points.append([x_samp,varZ])
			
		
		" tail section "
		samples2 = arange(p2[0], p3[0] + self.sampleResolution, self.sampleResolution)
		for x_samp in samples2:
			varZ = p2[1]
			points.append([x_samp,varZ])
		
		
		return points

