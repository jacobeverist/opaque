import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import Image
from math import *
from copy import *
import pylab


class VoronoiFit:

	def __init__(self, pnts = []):

		self.points = copy(pnts)
		numPoints = len(self.points)
		origin = [0.0,0.0]

		maxDist = 0.0
		minDist = 1e100
		avgDist = 0
		sumDist = 0

		" FIXME: add in points if there are not enough "
		newPoints = self.fillPoints(self.points)
		
		self.points = newPoints
		numPoints = len(self.points)
		

		for i in range(len(self.points)-1):
			p0 = self.points[i]
			p1 = self.points[i+1]
			dist = sqrt((p1[0]-p0[0])**2 + (p1[1]-p0[1])**2)

			sumDist += dist

			if dist> maxDist:
				maxDist = dist

			if dist < minDist:
				minDist = dist

		avgDist = sumDist / (len(self.points)-1)

		#print maxDist, minDist, avgDist

		#self.points.reverse()

		#print "spline fitting:", self.points
		self.curve = SplineFit(self.points, smooth = 0.1)

		self.draw()

	def fillPoints(self, points):
		
		new_points = []
		max_spacing = 0.02
		
		for i in range(len(points)-1):
			p0 = points[i]
			p1 = points[i+1]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			vec[0] /= dist
			vec[1] /= dist
			
			new_points.append(copy(p0))
			
			if dist > max_spacing:
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					new_points.append(newP)
		
		new_points.append(copy(points[-1]))
		
		return new_points
		
	def draw(self):

		xP = []
		yP = []

		for p in self.points:
			xP.append(p[0])
			yP.append(p[1])

		pylab.clf()
		pylab.scatter(xP,yP, color='r')
		self.curve.drawSpline()
		#pylab.scatter(xP,yP)
		
		pylab.xlim(-4.5,4.5)
		pylab.ylim(-4,4)
				
		#pylab.xlim(-1.0,1.0)
		#pylab.ylim(-0.5,0.5)
		
		pylab.savefig("curveFit.png")

	def getU(self,u):
		return self.curve.getU(u)

	def getUVector(self,u, iter=0.01):
		
		vecSum = [0.0,0.0]
		for i in range(5):
			tVec = self.curve.getUVector(0.8-i*0.01)
			#print tVec
			vecSum[0] += tVec[0]
			vecSum[1] += tVec[1]
		
		mag = sqrt(vecSum[0]**2 + vecSum[1]**2)
		vecSum[0] = vecSum[0]/mag
		vecSum[1] = vecSum[1]/mag
			
		return vecSum

	def findClosestPoint(self, point):
		minDist, u, linePoint1 = self.curve.findClosestPoint(point)

		vec = self.curve.getUVector(u)
		#print "u =", u
		#return minDist, copy(linePoint1), vec

		#p0 = self.curve.getU(u)
		#p1 = self.curve.getU(u + 0.01)

		#vec = [p1[0]-p0[0],p1[1]-p0[1]]
		#mag = sqrt(vec[0]**2 + vec[1]**2)
		#vec = [vec[0]/mag, vec[1]/mag]

		return minDist, copy(linePoint1), vec

	def findClosestPoint2(self, point):
		# TWO CANDIDATES
		# 1. check each of the path points
		# 2. check orthogonal distance to each edge if less than both terminating points

		pointDist = []
		for p in self.points:
			dist = sqrt((point[0]-p[0])**2 + (point[1]-p[1])**2)
			pointDist.append(dist)

		orthoDist = []
		tanVecs = []
		closePoints = []

		for i in range(len(self.points)-1):
			p0 = copy(self.points[i])
			p1 = copy(self.points[i+1])

			# rotate each edge to 0 degrees, and compare point distance to the horizontal
			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			mag = sqrt(vec[0]**2 + vec[1]**2)
			if mag == 0:
				orthoDist.append( 1e100 )
				tanVecs.append([1.0, 0.0])
				closePoints.append(copy(self.points[i]))
				print "degenerate"
				continue

			eAngle = acos(vec[0]/mag)
			if asin(vec[1]/mag) < 0:
				eAngle = -eAngle

			# rotate the edge
			rotVec = [vec[0] * cos(eAngle) + vec[1] * sin(eAngle), -vec[0] * sin(eAngle) + vec[1] * cos(eAngle)]

			oPoint = [point[0]-p0[0], point[1]-p0[1]]
			oPoint = [oPoint[0] * cos(eAngle) + oPoint[1] * sin(eAngle), -oPoint[0] * sin(eAngle) + oPoint[1] * cos(eAngle)]

			dist = abs(oPoint[1])
			if oPoint[0] < 0.0:
				orthoDist.append( pointDist[i])
				tanVecs.append([vec[0]/mag, vec[1]/mag])
				closePoints.append(copy(self.points[i]))

			elif oPoint[0] > mag:
				orthoDist.append(pointDist[i+1])
				tanVecs.append([vec[0]/mag, vec[1]/mag])
				closePoints.append(copy(self.points[i+1]))

			else:
				orthoDist.append(dist)

				tanVecs.append([vec[0]/mag, vec[1]/mag])
				nearPoint = [oPoint[0],0.0]
				nearPoint = [nearPoint[0] * cos(eAngle) - nearPoint[1] * sin(eAngle), nearPoint[0] * sin(eAngle) + nearPoint[1] * cos(eAngle)]
				closePoints.append(copy(nearPoint))

		minDist = 1e100
		minIndex = 0
		for i in range(len(orthoDist)):
			if orthoDist[i] < minDist:
				minDist = orthoDist[i]
				minIndex = i
				
		return minDist, copy(closePoints[minIndex]), copy(tanVecs[minIndex])
