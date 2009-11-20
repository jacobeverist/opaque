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

class SineFit:

	def __init__(self, freq, amp, lead_len, trail_len, origin, oVec):

		# control points are:
		# 1. leading edge origin
		# 2. leading edge terminator
		# 3. trailing edge origin
		# 4. trailing edge terminator

		# cosine x length for one cycle
		self.x_len = pi / freq 
		self.freq = freq
		self.amp = amp
		self.origin = origin
		self.oVec = oVec
		self.lead_len = lead_len
		self.trail_len = trail_len

		# normalize oVec
		mag = sqrt(oVec[0]**2 + oVec[1]**2)
		self.oVec = [oVec[0]/mag, oVec[1]/mag]

		# FIXME:  what to do if either lead or trailing lengths are zero?
		p1 = origin
		p2 = [p1[0] + oVec[0]*lead_len, p1[1] + oVec[1]*lead_len]
		p3 = [p2[0] + oVec[0]*self.x_len, p2[1] + oVec[1]*self.x_len]
		p4 = [p3[0] + oVec[0]*trail_len, p3[1] + oVec[1]*trail_len]

		#print p1, p2, p3, p4
		self.control = {'p1':p1, 'p2':p2, 'p3':p3, 'p4':p4}

	def draw(self):
		points = self.getPoints()	

		xP = []
		zP = []
		for p in points:
			xP.append(p[0])
			zP.append(p[1])

		#pylab.plot(xP,zP)
		pylab.xlim(-2,5)
		pylab.ylim(-2,2)
		#pylab.show()

	def findClosestPoint(self, point):

		# 3 hypotheses: closest point on the trailing edge, the leading edge, or the cosine curve

		# 1. closest point simple if between boundaries of trailing and leading edges
		# assuming oVec = [1,0] and origin = [0.0,0.0]
		
		# first rotate point inline with the curve coordinate system
		x = self.oVec[0] # assuming oVec is normalized
		y = self.oVec[1]

		if y != 0.0:
			B = 1 / (x*x/y + y)
			A = B * x / y
		else:
			if x > 0:
				A = 1.0
				B = 0.0
			else:
				A = -1.0
				B = 0.0

		xRot = point[0]*A + point[1]*B
		zRot = -point[0]*B + point[1]*A
		newPoint = [xRot, zRot]
		

		p1 = self.control['p1']
		p2 = self.control['p2']
		p3 = self.control['p3']
		p4 = self.control['p4']

		# Hypothesis 1: closest point on the leading edge
		if newPoint[0] <= p1[0]:
			min1 = sqrt((newPoint[0]-p1[0])**2 + (newPoint[1]-p1[1])**2)
			curvePoint1 = p1

		elif newPoint[0] > p1[0] and newPoint[0] < p2[0]:
			min1 = fabs(newPoint[1])
			curvePoint1 = [newPoint[0],0.0]

		elif newPoint[0] >= p2[0]:
			min1 = sqrt((newPoint[0]-p2[0])**2 + (newPoint[1]-p2[1])**2)
			curvePoint1 = p2
		else:
			print "points not in expected configuration"
			raise

		# Hypothesis 2: closest point on the trailing edge
		if newPoint[0] <= p3[0]:
			min2 = sqrt((newPoint[0]-p3[0])**2 + (newPoint[1]-p3[1])**2)
			curvePoint2 = p3

		elif newPoint[0] > p3[0] and newPoint[0] < p4[0]:
			min2 = fabs(newPoint[1])
			curvePoint2 = [newPoint[0],0.0]

		elif newPoint[0] >= p4[0]:
			min2 = sqrt((newPoint[0]-p4[0])**2 + (newPoint[1]-p4[1])**2)
			curvePoint2 = p4

		else:
			print "points not in expected configuration"
			raise

		# Hypothesis 3: closest point on the cosine curve
		# evenly sample the cosine curve to find the closest point
		samples = arange(0.0,self.x_len + 0.01, 0.01)

		min3 = 100000
		curvePoint3 = [0.0,0.0]

		for x_samp in samples:
			z = self.amp * sin(self.freq*x_samp)

			dist = sqrt((newPoint[0]-(self.lead_len + x_samp))**2 + (newPoint[1]-z)**2)

			if dist < min3:
				min3 = dist
				curvePoint3 = [self.lead_len + x_samp,z]

		min = 0.0
		pnt = []

		#print point, ":", min1, curvePoint1, min2, curvePoint2, min3, curvePoint3

		if min1 <= min2 and min1 <= min3:
			min = min1
			pnt = curvePoint1
			#print "check1"

		elif min2 <= min1 and min2 <= min3:
			min = min2
			pnt = curvePoint2
			#print "check2"

		elif min3 <= min2 and min3 <= min1:
			min = min3
			pnt = curvePoint3
			#print "check3"

		else:
			raise

		#print "min dist =", min, "pnt =", pnt

		xP = [pnt[0],newPoint[0]]
		zP = [pnt[1],newPoint[1]]
		#pylab.scatter(xP,zP)

		# min = minimum distance
		# u = d/c
		# curvePoint = point on the curve
		if pnt == p4:
			return min, 1.0, copy(pnt)
		elif pnt == p1:
			return min, 0.0, copy(pnt)
		else:
			return min, 0.5, copy(pnt)

	def getPoints(self):
		# create a distribution of points along the curve

		points = []

		#p1 = self.control['p1']
		#p2 = self.control['p2']
		#lead_len = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

		samples = arange(0.0,self.lead_len + 0.01, 0.01)
		for samp in samples:
			pnt = [self.origin[0] + samp*self.oVec[0], self.origin[1] + samp*self.oVec[1]]
			points.append(pnt)

		samples = arange(0.0,self.x_len + 0.01, 0.01)
		for x_samp in samples:
			z = self.amp * sin(self.freq*x_samp)
			#z = self.amp * cos(self.freq*x_samp) - self.amp
			points.append([self.lead_len + x_samp,z])

		samples = arange(0.0,self.trail_len + 0.01, 0.01)
		p3 = self.control['p3']
		for samp in samples:
			pnt = [p3[0] + samp*self.oVec[0], p3[1] + samp*self.oVec[1]]
			points.append(pnt)

		return points


