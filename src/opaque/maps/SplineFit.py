
from math import sqrt
import scipy.interpolate
import random
from random import gauss
from copy import copy
from math import floor, asin, acos, cos, sin
from time import time
from functions import closestAngle

from scipy.spatial import cKDTree
from numpy import array

class SplineFit:

	def __init__ (self, points, smooth = 0.1, kp = 5):
		
		random.seed(0)		

		if len(points) <= 1:
			print "ERROR, not enough points", points
			raise

		if len(points) <= 5:

			print "SplineFit received points:", points

			max_spacing = 0.08

			newPath3 = []

			" make sure path is geq than 5 points "
			while len(newPath3) <= 5:

				max_spacing /= 2
				
				newPath3 = [copy(points[0])]
									
				for i in range(len(points)-1):
					p0 = points[i]
					p1 = points[(i+1)]
					dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
		
					vec = [p1[0]-p0[0], p1[1]-p0[1]]

					if dist > 0.0:
						vec[0] /= dist
						vec[1] /= dist
					
					#print "dist=", repr(dist), repr(max_spacing)
					
					if dist > max_spacing:
						" cut into pieces max_spacing length or less "
						numCount = int(floor(dist / max_spacing))
						
						for j in range(1, numCount+1):
							newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
							newPath3.append(newP)

					newPath3.append(copy(p1))            

				newP = []
				for p in newPath3:
					if p not in newP:
						newP.append(p)
				newPath3 = newP
					
			points = newPath3
			print "SplineFit recomputed points:", points

		self.pointSet = points
		self.smoothNess = smooth
		self.kp = kp
		
		#print "spline received", len(points), "points"
		
		newP = []
		for p in points:
			if p not in newP:
				newP.append(p)
				
		# unzip the points
		intArray = [[],[]]
#		for p in points:

		" perturb the points to prevent degeneracies "
		for p in newP:
			try:
				intArray[0].append(p[0] + gauss(0.0,0.0001))
			except:
				intArray[0].append(p[0])

			try:
				intArray[1].append(p[1] + gauss(0.0,0.0001))
			except:
				intArray[1].append(p[1])
				
		" performing spline fit "
		try:
			self.tck, self.u = scipy.interpolate.splprep(intArray, s = self.smoothNess, k=self.kp)
		except:
			print "points:", points
			print "intArray:", intArray
			raise


		" precomputed points "
		self.densePoints = []
		self.distPoints = []
		self.denseAngles = []
		self.denseTree = None


	def findU(self, point):

		minDist, uVal, uPoint = self.findClosestPoint(point)
		
		return uVal
		
		#return min, min_i/1000.0, self.densePoints[min_i]
		
		
		samples = scipy.arange(0.0,1.0,0.01)
		sample_points = self.getUVecSet(samples)
		
		minI = -1
		minDist = 1e100
		
		for i in range(len(sample_points)):
			dist = sqrt((sample_points[i][0]-point[0])**2 + (sample_points[i][1]-point[1])**2)
			if dist < minDist:
				minI = i
				minDist = dist
				
		return minI * 0.01

	def getUniformSamples(self, spacing = 0.04, interpAngle = False):
		
		samples = scipy.arange(0.0,1.05,0.05)
		sample_points = self.getUVecSet(samples)
		sample_points = self.makePointsUniform(sample_points, max_spacing = spacing, interpAngle = interpAngle)
		return sample_points
	
	def makePointsUniform(self, points, max_spacing = 0.04, interpAngle = False):
		
		" make the points uniformly distributed "
		
		new_points = []
		
		for i in range(len(points)-1):
			p0 = points[i]
			p1 = points[i+1]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
	
			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			vec[0] /= dist
			vec[1] /= dist
			
			new_points.append(copy(p0))
			
			ang0 = p0[2]
			ang1 = p1[2]
			
			if dist > max_spacing:
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				angStep = (ang1-ang0) / float(numCount+1)
				
				
				for j in range(1, numCount+1):
					if interpAngle:
						newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1], p0[2]+j*angStep]
					else:
						newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1], p0[2]]
						
					new_points.append(newP)
		
		return new_points		

	def getTransformCurve(self, offset = 0.0, compareAngle = None):
		
		#points = self.getUniformSamples(spacing = 0.01)

		samples = scipy.arange(0.0,1.0,0.01)
		sample_points = self.getUVecSet(samples)
		points = self.makePointsUniform(sample_points, max_spacing = 0.01)

		" transform angle to y coordinate, and total distance to the x coordinate "

		totalDist = 0.0
		candAngle = points[0][2]+offset
		if compareAngle != None:
			nextAngle = closestAngle(compareAngle, candAngle)
			print "compareAngle:", compareAngle, candAngle, nextAngle
		else:
			nextAngle = candAngle
		
		points_trans = [[totalDist, nextAngle]]
		for i in range(len(points)-1):
			p0 = points[i]
			p1 = points[i+1]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
			totalDist += dist
			prevAngle = points_trans[i][1]
			candAngle = p1[2] + offset
			nextAngle = closestAngle(prevAngle, candAngle)
			if i == 0:
				print "closestAngle:", offset, prevAngle, candAngle, nextAngle
			points_trans.append([totalDist, nextAngle])
			
		" thin out the points "
		points2 = []
		for i in range(len(points_trans)):
			if i % 2 == 0:
				points2.append(points_trans[i])
			
		return points2



	def getUOfDist(self, startU, dist, distIter = None):

		" search algorithm.   Find a high and a low "
		" subdivide and determine if higher or lower than minimum point "
		" repeat to a desired precision "

		magDist = abs(dist)
		
		returnVal = -1
		#print "dist =", dist, "startU =", startU, "magDist =", magDist
		
		if dist >= 0:
			highU = 1.0
			lowU = startU
			
			if distIter != None:
				termDist = self.dist(startU, highU, iter = distIter)
			else:
				termDist = self.dist(startU, highU)
			#print "termDist =", termDist
			if termDist <= magDist:
				returnVal = highU
				return returnVal

		else:
			highU = startU
			lowU = 0.0
			if distIter != None:
				termDist = self.dist(startU, lowU, iter = distIter)
			else:
				termDist = self.dist(startU, lowU)
			#print "termDist =", termDist

			if termDist <= magDist:
				returnVal = lowU
				return returnVal
		
		if returnVal == -1:
			count = 0
			
			while True:
				" find U between the highest and lowest"
				midU = (highU - lowU)/2.0 + lowU
				if distIter != None:
					midDist = self.dist(startU, midU, iter = distIter)
				else:
					midDist = self.dist(startU, midU)
				
				#print "midU =", midU, "midDist =", midDist
				
				if dist >= 0:
					if midDist > magDist:
						highU = midU
						lowU = lowU
					else:
						highU = highU
						lowU = midU
				else:
					if midDist > magDist:
						highU = highU
						lowU = midU
					else:
						highU = midU
						lowU = lowU
				
				#print "highU =", highU, "lowU =", lowU
				
				" terminate when close enough "
				if abs(magDist-midDist) < 0.01:
					#print "returning u =", midU, "for midU =", midU
					#print
					returnVal = midU
					return returnVal
	
				count += 1
				if count > 30:
					print "Fail getUofDist:"
					print startU, dist, highU, lowU, midU, midDist, magDist, termDist
					#raise

					" brute force search"
					bestDiff = 1e100
					bestIndex = int(startU*1000.)
					initDist = self.distPoints[bestIndex]
			
					print "bestIndex, initDist =", bestIndex, initDist
			
					if dist >= 0:
						setIndices = range(bestIndex,1000)
						for k in setIndices:
							currDist = self.distPoints[k] - initDist
							deltaDiff = abs(magDist-currDist)
							
							if deltaDiff < bestDiff:
								bestDiff = deltaDiff
								bestIndex = k
					else:
						setIndices = range(0,bestIndex+1)
						setIndices.reverse()
						for k in setIndices:
							currDist = initDist - self.distPoints[k]
							deltaDiff = abs(magDist-currDist)
							
							if deltaDiff < bestDiff:
								bestDiff = deltaDiff
								bestIndex = k
			
					uVal = bestIndex * 0.001
					print "currDist, bestIndex, bestDiff, uVal =", currDist, bestIndex, bestDiff, uVal
					returnVal2 = uVal
					
					return uVal

		#print "returnVals:", returnVal, returnVal2
		#if returnVal - returnVal2 > 0.01:
		#	print "WOW"
		
		#return returnVal
	
	

	def getPointOfDist(self, dist):
		# return a point of dist distance along the spline from the origin

		# if spline not long enough!
		
		totalLen = self.dist_u(1.0)
		
		#totalLen = self.length()
		if dist > totalLen:
			#point = self.getU(1.0)
			
			point = self.densePoints[-1]
			
			extraLen = dist - totalLen
			
			angle = point[2]
			
			xd = point[0] + extraLen*cos(angle)
			yd = point[1] + extraLen*sin(angle)
			
			return [xd,yd,angle]
			#raise

		if dist == totalLen:
			point = self.densePoints[-1]
			#point = self.getU(1.0)
			#print "returning 1.0", point
			return point

		# spline distance should not be negative
		if dist < 0:
			point = self.densePoints[0]
			#point = self.getU(0.0)
			
			extraLen = dist
			
			angle = point[2]
			
			xd = point[0] + extraLen*cos(angle)
			yd = point[1] + extraLen*sin(angle)
			
			return [xd,yd,angle]
			#raise

		if dist == 0.0:
			point = self.densePoints[0]
			#point = self.getU(0.0)
			#print "returning 0.0", point
			return point


		uInd = self.dist_search(dist)		
		return copy(self.densePoints[uInd])

		estU = 0.5
		jumpVal = 0.5

		count = 0

		#print "target =", dist

		while True:

			count += 1
			if count > 20:
				#"returning empty"
				return []

			if count < 4:
				estDist = self.dist(0.0, estU, iter = 0.005)
			else:
				estDist = self.dist(0.0, estU, iter = 0.001)
				
			jumpVal /= 2.0
			#print "estU, estDist, jumpVal, count =", estU, estDist, jumpVal, count

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

	def dist_search(self, x):
		if len(self.distPoints) == 0:
			self.precompute()

		lo = 0
		hi = len(self.distPoints)
		
		while lo < hi:
			mid = (lo+hi)//2
			midval = self.distPoints[mid]
			if midval <= x:
				lo = mid+1
			elif midval > x: 
				hi = mid
			else:
				return mid
		
			
		return mid


	def dist_u(self, u):
		if len(self.distPoints) == 0:
			self.precompute()

		if u <= 0.0:
			#print "dist_u(", u, ") = 0"
			return self.distPoints[0]
		
		if u >= 1.0:
			#print "dist_u(", u, ") = -1"
			return self.distPoints[-1]
		
		ind = u * 1000.
		indInt = int(ind)
		#print "dist_u(", u, ") =", indInt
		
		#print "returning", self.distPoints[indInt]
		return self.distPoints[indInt]


	def precompute(self):
		
		iter = 0.001
		unew = scipy.arange(0.0, 1.0 + iter, iter)
		out = scipy.interpolate.splev(unew,self.tck)
		
		" 1001 points "
		self.densePoints = []
		self.denseAngles = []
		kdInput = []
		for i in range(len(out[0])):

			if i == len(out[0])-1:
				newPoint1 = [out[0][i-1],out[1][i-1]]
				newPoint2 = [out[0][i],out[1][i]]
			else:
				newPoint1 = [out[0][i],out[1][i]]
				newPoint2 = [out[0][i+1],out[1][i+1]]
				
			" tangent vector "
			vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]
	
			" normalize "
			mag = sqrt(vec[0]**2 + vec[1]**2)
			vec = [vec[0]/mag, vec[1]/mag]
			
			angle = acos(vec[0])
			if asin(vec[1]) < 0:
				angle = -angle
			
			kdInput.append([out[0][i],out[1][i]])
			self.densePoints.append([out[0][i],out[1][i], angle])
			self.denseAngles.append(angle)
	
		self.denseTree = cKDTree(array(kdInput))

		#print "self.densePoints:", self.densePoints
		
		" 1000 distances "
		totalDist = 0.0
		self.distPoints = [0.0]
		for i in range(len(self.densePoints)-1):
			
			pnt0 = self.densePoints[i]
			pnt1 = self.densePoints[i+1]
			currDist = sqrt( (pnt1[0]-pnt0[0])**2 + (pnt1[1]-pnt0[1])**2 )
			
			totalDist += currDist
			
			self.distPoints.append(totalDist)
			
			
		#print "self.distPoints:", self.distPoints

	# compute the distance along the curve between the points parameterized by u1 and u2
	def dist(self, u1, u2, iter = 0.005):

		if len(self.distPoints) == 0:
			self.precompute()

		
		
		totalDist = 0.0

		if u1 >= u2:
			tempU = u1
			u1 = u2
			u2 = tempU
		
		
		dist1 = self.dist_u(u1)
		dist2 = self.dist_u(u2)

		#print "return dist2-dist1=", dist2, dist1, dist2-dist1
		return dist2-dist1
	

		unew = scipy.arange(u1, u2 + iter, iter)
		out = scipy.interpolate.splev(unew,self.tck)

		points = []
		if len(unew) > 1:
			for i in range(len(out[0])):
				points.append([out[0][i],out[1][i]])
		elif len(unew) == 1:
			points = [out]

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
			#print "totalDist =", totalDist
			origin = points.pop(0)

		return totalDist

	# compute the tangent vector to u1
	def getUVector(self, u1, iter = 0.01):

		#print "getUVector(", u1, ")"

		if u1 < 0.0 or u1 > 1.0:
			print "ERROR: u =", u1
			raise

		iterVal = 1.0-iter
		#print "u1 > 1.0 - iter:", repr(u1), "> ", repr(iterVal)
		if u1 >= iterVal:
			u1 = u1-iter
			#print "u1 changed to", u1
		#else:
		#	print "u1 not changed", u1

		#print "u1=", u1

		unew1 = [u1]
		newPoint1 = scipy.interpolate.splev(unew1,self.tck)
		unew2 = [u1 + iter]
		newPoint2 = scipy.interpolate.splev(unew2,self.tck)
		
		#print unew1, newPoint1, unew2, newPoint2

		# tangent vector
		vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]

		#print "vec=", vec
		# normalize
		mag = sqrt(vec[0]**2 + vec[1]**2)
		vec = [vec[0]/mag, vec[1]/mag]

		#print "returning vec=", vec
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
		
		"""
		if plot:
			pylab.clf()
			self.drawSpline()
			rearSpline.drawSpline()

			for i in range(0,len(sample_points)):
				xP = [sample_points[i][0], sample_points2[i][0]]
				yP = [sample_points[i][1], sample_points2[i][1]]
				pylab.plot(xP,yP,color='0.5')
				pylab.scatter(xP,yP,color='0.5',linewidth='1')

			pylab.xlim(-6,10)
			pylab.ylim(-4,10)
		
			pylab.title("Cost = %f" % cost)
			pylab.savefig("test_%04u.png" % index)	
		"""
			
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

		"""
		if plot:
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
		"""

		return cost

	def drawSpline(self, clr = '0.5'):
		unew = scipy.arange(0, 1.01, 0.01)
		out = scipy.interpolate.splev(unew,self.tck)
		#pylab.plot(out[0],out[1], color = clr)	

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
	
	def getUVecSet(self, u_set):
		
		for u in u_set:
			if u < 0.0 or u > 1.0:
				print "ERROR: u =", u
				raise

		newPoints = scipy.interpolate.splev(u_set,self.tck)
		
		# zip the points together
		zipPoints = []
		
		iter = 0.01
		
		for i in range(len(u_set)):
			u = u_set[i]
			if u > 1.0 - iter:
				u = u-iter
	
			unew1 = [u]
			newPoint1 = scipy.interpolate.splev(unew1,self.tck)
			unew2 = [u + iter]
			newPoint2 = scipy.interpolate.splev(unew2,self.tck)
	
			# tangent vector
			vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]
	
			# normalize
			mag = sqrt(vec[0]**2 + vec[1]**2)
			vec = [vec[0]/mag, vec[1]/mag]
			
			angle = acos(vec[0])
			if asin(vec[1]) < 0:
				angle = -angle
					
			zipPoints.append([newPoints[0][i],newPoints[1][i], angle])
		

		#for i in range(0,len(newPoints[0])):
		#	zipPoints.append([newPoints[0][i],newPoints[1][i]])

		return zipPoints

	def findClosestPoint(self, pos):
		# find closest location on spline to pos
		# we do this in three nested searches,
		# sampling 10 uniformly distributed points at each level	

		#sample = scipy.arange(0.0,1.01,0.01)
		#points = self.getUSet(sample)

		if len(self.distPoints) == 0:
			self.precompute()

		" remove angle if it's there "
		queryPoint = [pos[0],pos[1]]

		minVal, min_i = self.denseTree.query(array(queryPoint))
		
		#minVal, min_i = self.findMinDistancePoint(self.densePoints, pos)

		uVal = min_i/1000.0
		
		if uVal > 1.0:
			uVal = 1.0
		elif uVal < 0.0:
			uVal = 0.0

		return minVal, uVal, self.densePoints[min_i]

	def findMinDistancePoint(self, points, pos):

		# compute distances of each point to pos
		minVal = 10e9
		min_i = 0
		for i in range(0,len(points)):
			p = points[i]
			dist = sqrt((p[0]-pos[0])**2 + (p[1]-pos[1])**2)
			if dist < minVal:
				min_i = i
				minVal = dist

		return minVal, min_i

if __name__ == '__main__':
	
	dirName = "../../testData/poseTest"
	numPoses = 16
	
	localPostures = []
	centerCurves = []
	
	estPoses = []
	gndPoses = []

	points = [[6.393668170242247, 8.942828779300818], [6.433224226989662, 8.948771705734424], [6.472780283737078, 8.954714632168029], [6.5123363404844925, 8.960657558601635], [6.528662149296044, 8.963110358216008]]		

	newSpline = SplineFit(points, smooth = 0.1)
	

	"""
	for i in range(numPoses):

		f = open(dirName + "/posture%04u.txt" % i, 'r')
		localPostures.append(eval(f.read().rstrip()))
		f.close()
	print len(localPostures), numPoses
	print time()
	for i in range(numPoses):    
		" compute a fitted curve of a center point "
		centerCurves.append(SplineFit(localPostures[i], smooth = 0.5, kp = 2))
		vec = centerCurves[i].getUVector(0.5)
		pnt = centerCurves[i].getU(0.5)

	print time()
	"""
        


