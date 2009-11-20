import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from math import *
from scipy.optimize import *
import scipy.interpolate
import pylab

class SplineFit:

	def __init__ (self, points, smooth = 0.1, kp = 5):
		self.pointSet = points
		self.smoothNess = smooth
		self.kp = kp
		
		newP = []
		for p in points:
			if p not in newP:
				newP.append(p)
				
		# unzip the points
		intArray = [[],[]]
#		for p in points:
		for p in newP:
			intArray[0].append(p[0])
			intArray[1].append(p[1])

		#print "smooth = ", self.smoothNess, "input = ", intArray

		dataOutput = open("dataFile.txt", 'w')
		dataOutput.write(repr(intArray))
		dataOutput.close()

		#print len(intArray[0]), len(intArray[1])

		#print "performing spline fit"
		self.tck, self.u = scipy.interpolate.splprep(intArray, s = self.smoothNess, k=self.kp)
		#print "fit complete"

		#print "result =", self.u
		#print self.tck
		#unew = scipy.arange(0, 1.01, 0.01)
		#out = scipy.interpolate.splev(unew,tck)

	def findClosestFromPointSet(self):
		pass

	def getPointOfDist(self, dist):
		# return a point of dist distance along the spline from the origin

		# if spline not long enough!
		totalLen = self.length()
		if dist > totalLen:
			raise

		if dist == totalLen:
			point = self.getU(1.0)
			return point

		# spline distance should not be negative
		if dist < 0:
			raise

		if dist == 0.0:
			point = self.getU(0.0)
			return point

		estU = 0.5
		jumpVal = 0.5

		count = 0

		while True:

			count += 1
			if count > 10:
				return []

			estDist = self.dist(0.0, estU)
			jumpVal /= 2.0

			#print estDist, dist

			if abs(dist-estDist) < 0.01:
				point = self.getU(estU)
				return point

			# increase the estimate
			if estDist < dist:
				estU += jumpVal

			# decrease the estimate
			if estDist > dist:
				estU -= jumpVal

	# compute the distance along the curve between the points parameterized by u1 and u2
	def dist(self, u1, u2, iter = 0.0001):
		totalDist = 0.0

		if u1 >= u2:
			raise


		unew = scipy.arange(u1, u2 + iter, iter)
		out = scipy.interpolate.splev(unew,self.tck)

		points = []
		for i in range(len(out[0])):
			points.append([out[0][i],out[1][i]])

		#print "measuring", u1, "to", u2+iter, "with", len(points), "results"
		#print points

		if len(unew) <= 1:
			return 0.0

		origin = points.pop(0)
		while len(points) > 0:
			pnt = points[0]
			#print pnt, origin
			dist = sqrt( (pnt[0]-origin[0])**2 + (pnt[1]-origin[1])**2 )
			totalDist += dist
			origin = points.pop(0)

		return totalDist

	# compute the tangent vector to u1
	def getUVector(self, u1, iter = 0.01):

		if u1 < 0.0 or u1 > 1.0:
			print "ERROR: u =", u
			raise

		if u1 > 1.0 - iter:
			u1 = u1-iter

		unew1 = [u1]
		newPoint1 = scipy.interpolate.splev(unew1,self.tck)
		unew2 = [u1 + iter]
		newPoint2 = scipy.interpolate.splev(unew2,self.tck)

		# tangent vector
		vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]

		# normalize
		mag = sqrt(vec[0]**2 + vec[1]**2)
		vec = [vec[0]/mag, vec[1]/mag]

		return vec

	# return the length of the curve
	def length(self, iter = 0.0001):
		totalDist = 0

		unew = scipy.arange(0, 1.0 + iter, iter)
		out = scipy.interpolate.splev(unew,self.tck)

		points = []
		for i in range(len(out[0])):
			points.append([out[0][i],out[1][i]])

		origin = points.pop(0)
		while len(points) > 0:
			pnt = points[0]
			dist = sqrt( (pnt[0]-origin[0])**2 + (pnt[1]-origin[1])**2 )
			totalDist += dist
			origin = points.pop(0)

		return totalDist

	" this matching algorithm tries to put the splines directly on top of each other "
	def matchSpline2(self, rearSpline, plot = False, index = 0):

		samples = scipy.arange(0.0,1.0,0.1)
		sample_points = self.getUSet(samples)
		sample_points2 = rearSpline.getUSet(samples)
		distances = []
		
		for i in range(len(sample_points)):
			dist = sqrt((sample_points[i][0]-sample_points2[i][0])**2 + (sample_points[i][1]-sample_points2[i][1])**2)
			distances.append(dist)

		cost = sum(distances)	
		
		if plot:
			pylab.clf()
			self.drawSpline()
			rearSpline.drawSpline()

			for i in range(0,len(sample_points)):
				xP = [sample_points[i][0], sample_points2[i][0]]
				yP = [sample_points[i][1], sample_points2[i][1]]
				pylab.plot(xP,yP,color='0.5')
				pylab.scatter(xP,yP,color='0.5',linewidth='1')

			#pylab.xlim(3,7)
			#pylab.ylim(1,4)
			pylab.xlim(-6,10)
			pylab.ylim(-4,10)
		
			pylab.title("Cost = %f" % cost)
			pylab.savefig("test_%04u.png" % index)	
			
		return cost
	
	def matchSpline(self, rearSpline, plot = False, index = 0):
		# this spline is fore, the new one is rear
		# the assumption here is represented by the function call
		# whether this is forward, or the other is forward

		# parametric boundaries of the spline sections we are comparing
		u0_F = 0.0
		u0_R = 0.0
		u1_F = 1.0
		u1_R = 1.0

		# since this is forward, u = 0.0 for fore, needs to match with some u on rear
		# get the this point for u=0.0
		self.originF = self.getU(0.0)

		# now find the closest point on the rear spline to self.originF
		dist, u0_R, self.originR = rearSpline.findClosestPoint(self.originF)

		# these are now our origins: self.originF, self.originR

		# now find the terminator of fore from rear
		self.termR = rearSpline.getU(1.0)
		dist, u1_F, self.termF = self.findClosestPoint(self.termR)

		# now we sample ten points in the fore spline boundaries
		diff = u1_F - u0_F
		inc = diff/10.0
		#print "u0_F =", u0_F, "u1_F =", u1_F
		#print "u0_R =", u0_R, "u1_R =", u1_R

		# if the end points are the same, lets default to a full length comparison
		if u1_F == u0_F:
			samples = scipy.arange(0.0,1.0,0.1)
		else:
			samples = scipy.arange(u0_F,u1_F,inc)

		sample_points = self.getUSet(samples)

		# for each sample point, compute the min distance to the rear spline
		distances = []
		closePoints = []
		for p in sample_points:
			dist, u_value, closest_p = rearSpline.findClosestPoint(p)
			distances.append(dist)
			closePoints.append(closest_p)

		cost = sum(distances)

		# now we sample ten points in the rear spline boundaries
		diff = u1_R - u0_R
		inc = diff/10.0

		if u1_R == u0_R:
			samples2 = scipy.arange(0.0,1.0,0.1)
		else:
			samples2 = scipy.arange(u0_R,u1_R,inc)

		sample_points2 = rearSpline.getUSet(samples2)

		# for each sample point, compute the min distance to the rear spline
		distances = []
		closePoints2 = []
		for p in sample_points2:
			dist, u_value, closest_p = self.findClosestPoint(p)
			distances.append(dist)
			closePoints2.append(closest_p)

		cost += sum(distances)

		if plot:
			#pylab.clf()
			#self.drawSpline()
			#rearSpline.drawSpline()

			pylab.scatter([self.originF[0], self.originR[0]], [self.originF[1], self.originR[1]])
			pylab.scatter([self.termF[0], self.termR[0]], [self.termF[1], self.termR[1]])

			for i in range(0,len(sample_points)):
				xP = [sample_points[i][0], closePoints[i][0]]
				yP = [sample_points[i][1], closePoints[i][1]]
				pylab.plot(xP,yP,color='0.5')
				pylab.scatter(xP,yP,color='0.5',linewidth='1')

			for i in range(0,len(sample_points2)):
				xP = [sample_points2[i][0], closePoints2[i][0]]
				yP = [sample_points2[i][1], closePoints2[i][1]]
				pylab.plot(xP,yP,color='0.5')
				pylab.scatter(xP,yP,color='0.5',linewidth='1')

			#pylab.xlim(3,7)
			#pylab.ylim(1,4)

			#pylab.title("Cost = %f" % cost)
			#pylab.savefig("pics/fit%04u.png" % index)

		return cost

	def drawSpline(self, clr = '0.5'):
		unew = scipy.arange(0, 1.01, 0.01)
		out = scipy.interpolate.splev(unew,self.tck)


		#pylab.plot(out[0],out[1], color = clr)	

		# flipped to align with Ogre x-z coord system
		#pylab.plot(out[1],out[0], color = clr)	
		pylab.plot(out[0],out[1], color = clr)	

	def getU(self, u):
		if u < 0.0 or u > 1.0:
			print "ERROR: u =", u
			raise

		unew = [u]
		newPoint = scipy.interpolate.splev(unew,self.tck)
		return newPoint

	def getUSet(self, u_set):
		for u in u_set:
			if u < 0.0 or u > 1.0:
				print "ERROR: u =", u
				raise

		newPoints = scipy.interpolate.splev(u_set,self.tck)
		# zip the points together
		zipPoints = []
		for i in range(0,len(newPoints[0])):
			zipPoints.append([newPoints[0][i],newPoints[1][i]])

		return zipPoints

	def findClosestPoint(self, pos):
		# find closest location on spline to pos
		# we do this in three nested searches,
		# sampling 10 uniformly distributed points at each level	

		sample = scipy.arange(0.0,1.01,0.01)
		points = self.getUSet(sample)
		min, min_i = self.findMinDistancePoint(points, pos)

		#print min, min_i

		if min_i == 0:
			# sample above index zero
			high_u = sample[1]
			low_u = sample[0]

		elif min_i == len(sample)-1:
			# sample below index zero 
			high_u = sample[9]
			low_u = sample[8]

		else:
			# sample above and below min index
			high_u = sample[min_i+1]
			low_u = sample[min_i-1]

		diff = high_u - low_u
		#print "high_u =", high_u
		inc = diff/10.0
		sample2 = scipy.arange(low_u,high_u,inc)

		sample3 = []
		for samp in sample2:
			sample3.append(samp)
		if high_u == 1.0:
			sample3 += [1.0]

		points2 = self.getUSet(sample3)
		min, min_i = self.findMinDistancePoint(points2, pos)

		# return distance, u value, and point
		return min, sample2[min_i], points2[min_i]

	def findMinDistancePoint(self, points, pos):

		# compute distances of each point to pos
		min = 10e9
		min_i = 0
		for i in range(0,len(points)):
			p = points[i]
			dist = sqrt((p[0]-pos[0])**2 + (p[1]-pos[1])**2)
			if dist < min:
				min_i = i
				min = dist

		return min, min_i




