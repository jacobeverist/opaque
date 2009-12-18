import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from numpy import *
from scipy.optimize import *
import graph
import csv
from LocalNode import *
import pylab
from matplotlib.patches import Circle

from pose import *
import gen_icp


class Pose:
	
	def __init__(self, pose = [0.0,0.0,0.0]):
		self.setEstPose(pose)

	def setEstPose(self, newPose):

	   self.estPose = copy(newPose)
	   self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)
	   self.vecAng = acos(self.estPose[0]/self.dist)
	   if asin(self.estPose[1]/self.dist) < 0:
			   self.vecAng = -self.vecAng

	   self.backR = array([[cos(self.vecAng), sin(self.vecAng)],[-sin(self.vecAng),cos(self.vecAng)]])
	   self.foreR = array([[cos(self.vecAng), -sin(self.vecAng)],[sin(self.vecAng),cos(self.vecAng)]])

	   self.R = array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])

	def convertLocalOffsetToGlobal(self, offset):

	   globalEst = [0.0,0.0,0.0]

	   finalVec = array([[offset[0]], [offset[1]]])
	   transVec = dot(transpose(self.R), finalVec)
	   resVec = dot(self.backR, transVec)
	   resVec[0, 0] += self.dist
	   tempVec = dot(self.foreR, resVec)
	   globalEst[0] = tempVec[0, 0]
	   globalEst[1] = tempVec[1, 0]
	   globalEst[2] = normalizeAngle(self.estPose[2] + offset[2])

	   return globalEst

	def convertGlobalPoseToLocal(self, pose):

	   " transform pnt to local coordinates"
	   globalVec = array([[pose[0]],[pose[1]]])

	   " perform translation correction "
	   tempVec = dot(self.backR, globalVec)
	   tempVec[0,0] -= self.dist
	   transVec = dot(self.foreR, tempVec)

	   " now, apply rotation correction with respect to origin "
	   localVec = dot(self.R, transVec)

	   localPose = [localVec[0,0], localVec[1,0], normalizeAngle(pose[2] - self.estPose[2])]

	   return localPose

	def convertLocalToGlobal(self, pnt):

	   finalVec = array([[pnt[0]], [pnt[1]]])
	   transVec = dot(transpose(self.R), finalVec)
	   resVec = dot(self.backR, transVec)
	   resVec[0, 0] += self.dist
	   tempVec = dot(self.foreR, resVec)

	   newPoint = [tempVec[0,0],tempVec[1,0]]

	   return newPoint

	def convertGlobalToLocal(self, pnt):

	   " transform pnt to local coordinates"
	   globalVec = array([[pnt[0]],[pnt[1]]])

	   " perform translation correction "
	   tempVec = dot(self.backR, globalVec)
	   tempVec[0,0] -= self.dist
	   transVec = dot(self.foreR, tempVec)

	   " now, apply rotation correction with respect to origin "
	   localVec = dot(self.R, transVec)

	   newPoint = [localVec[0,0], localVec[1,0]]
	   return newPoint



class AlphaMapCost:
	
	def __init__(self):
		self.count = 0
		self.fig = pylab.figure()
		self.ax = self.fig.add_subplot(111)
		
		self.radius = 2.5
		
		self.circle1 = [[0.0, 0.0, 0.0], self.radius]
		self.circle2 = [[0.0, 0.0, 0.0], self.radius]
	
	def getGuess2(self, offset, node1, node2):
		
		est1 = node1.getEstPose()
		est2 = node2.getEstPose()

		" first node stays in its current position "
		node = node1
		cntPoints = deepcopy(node.centerPoints)

		origin1 = node.convertLocalOffsetToGlobal([0.0, 0.0, 0.0])


		xP = []
		yP = []
		cntPoints1 = []
		for j in range(len(cntPoints)):
			pnt = cntPoints[j]

			xP.append(pnt[0])
			yP.append(pnt[1])
			
			cntPoints1.append(copy(pnt))

		spline1 = SplineFit(cntPoints1, kp=3)
		
		samples = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
		totalVec = [0.0, 0.0]
		for i in range(len(samples)):
		
			uvec = spline1.getUVector(samples[i])
			totalVec[0] += uvec[0]
			totalVec[1] += uvec[1]
		
		totalVec[0] /= float(len(samples))
		totalVec[1] /= float(len(samples))
		
		" convert the totalVec into local coordinates "
		globalVec = array([[totalVec[0]], [totalVec[1]]])
		finalVec = dot(node1.R, totalVec)
		#print finalVec
		#print finalVec[0]
		#print finalVec[1]
		
		" modified guess in local coordinates "
		newGuess = copy(offset)
		newGuess[0] -= finalVec[0] * 1.0
		newGuess[1] -= finalVec[1] * 1.0
		
		return newGuess
		
	def guessCost(self, offset, node1, node2):

		est1 = node1.getEstPose()
		est2 = node2.getEstPose()

		" first node stays in its current position "
		node = node1
		cntPoints = deepcopy(node.centerPoints)

		origin1 = node.convertLocalOffsetToGlobal([0.0,0.0,0.0])
		
		xP = []
		yP = []
		cntPoints1 = []
		for j in range(len(cntPoints)):
			pnt = cntPoints[j]
			
			finalVec = array([[pnt[0]], [pnt[1]]])
			transVec = dot(transpose(node.R), finalVec)
			resVec = dot(node.backR, transVec)
			resVec[0, 0] += node.dist
			tempVec = dot(node.foreR, resVec)
			pnt[0] = tempVec[0, 0]
			pnt[1] = tempVec[1, 0]		

			xP.append(pnt[0])
			yP.append(pnt[1])
			
			cntPoints1.append(copy(pnt))

		" second node is offset from is current position "
		newEst2 = copy(est2)

		finalVec = array([[offset[0]], [offset[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0, 0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		newEst2[0] = tempVec[0, 0]
		newEst2[1] = tempVec[1, 0]
		newEst2[2] = normalizeAngle(est1[2] + offset[2])
	
		" second node is offset from is current position "

		"""
		newEst2 = copy(est2)
		origin = copy(est1)
		
		ang = offset[2]				
		newEst2[0] += offset[0]
		newEst2[1] += offset[1]
		
		xdiff = newEst2[0] - origin[0]
		ydiff = newEst2[1] - origin[1]

		xnew = xdiff*cos(ang) - ydiff*sin(ang)
		ynew = xdiff*sin(ang) + ydiff*cos(ang)

		newEst2[0] = origin[0] + xnew
		newEst2[1] = origin[1] + ynew
		newEst2[2] += ang
		"""

		dist = sqrt(newEst2[0] ** 2 + newEst2[1] ** 2)
		vecAng = acos(newEst2[0] / dist)
		if asin(newEst2[1] / dist) < 0:
			vecAng = - vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)], [ - sin(vecAng), cos(vecAng)]])
		foreR = array([[cos(vecAng), - sin(vecAng)], [sin(vecAng), cos(vecAng)]])
		R = array([[cos(newEst2[2]), sin(newEst2[2])], [ - sin(newEst2[2]), cos(newEst2[2])]])
		
		cntPoints = deepcopy(node2.centerPoints)

		origin2 = [0.0, 0.0, 0.0]
		finalVec = array([[origin2[0]], [origin2[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0, 0] += dist
		tempVec = dot(foreR, resVec)
		origin2[0] = tempVec[0, 0]
		origin2[1] = tempVec[1, 0]		

		xP = []
		yP = []
		cntPoints2 = []
		for j in range(len(cntPoints)):
			pnt = cntPoints[j]
			
			finalVec = array([[pnt[0]], [pnt[1]]])
			transVec = dot(transpose(R), finalVec)
			resVec = dot(backR, transVec)
			resVec[0, 0] += dist
			tempVec = dot(foreR, resVec)
			pnt[0] = tempVec[0, 0]
			pnt[1] = tempVec[1, 0]		

			xP.append(pnt[0])
			yP.append(pnt[1])

			cntPoints2.append(copy(pnt))
		
	
		" two splines, transformed from local to global plus offset "
		#print cntPoints1
		spline1 = SplineFit(cntPoints1, kp=3)
		#print cntPoints2
		spline2 = SplineFit(cntPoints2, kp=3)
		
		cost = spline1.matchSpline2(spline2, plot=True, index=self.count)
		
		self.count += 1
				
		return cost

	def cost_func(self, offset, node1, node2):

		self.fig.clf()
		self.ax = self.fig.add_subplot(111)

		#offset = copy(offset)
		#temp = offset[0]
		#offset[0] = offset[2]
		#offset[2] = temp
	
		#radius = 2.8
		#radius = 4
		radius = 2.5
		#radius = 10.0
		
		est1 = node1.getEstPose()
		est2 = node2.getEstPose()
		
		
		" first node stays in its current position "
		node = node1
		vert = deepcopy(node.a_vert)

		origin1 = [0.0, 0.0, 0.0]
		finalVec = array([[origin1[0]], [origin1[1]]])
		transVec = dot(transpose(node.R), finalVec)
		resVec = dot(node.backR, transVec)
		resVec[0, 0] += node.dist
		tempVec = dot(node.foreR, resVec)
		origin1[0] = tempVec[0, 0]
		origin1[1] = tempVec[1, 0]		


		xP = []
		yP = []
		vert1 = []
		for j in range(len(vert)):
			pnt = vert[j]
			
			finalVec = array([[pnt[0]], [pnt[1]]])
			transVec = dot(transpose(node.R), finalVec)
			resVec = dot(node.backR, transVec)
			resVec[0, 0] += node.dist
			tempVec = dot(node.foreR, resVec)
			pnt[0] = tempVec[0, 0]
			pnt[1] = tempVec[1, 0]		

			xP.append(pnt[0])
			yP.append(pnt[1])
			
			vert1.append(copy(pnt))
		

		#xP.append(xP[0])
		#yP.append(yP[0])
		
		pylab.plot(xP, yP, color='0.0')


		" second node is offset from is current position "
		newEst2 = copy(est2)

		finalVec = array([[offset[0]], [offset[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0, 0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		newEst2[0] = tempVec[0, 0]
		newEst2[1] = tempVec[1, 0]
		newEst2[2] = normalizeAngle(est1[2] + offset[2])
		
	
		"""
		ang = offset[2]				
		newEst2[0] += offset[0]
		newEst2[1] += offset[1]

		xdiff = newEst2[0] - origin[0]
		ydiff = newEst2[1] - origin[1]

		xnew = xdiff*cos(ang) - ydiff*sin(ang)
		ynew = xdiff*sin(ang) + ydiff*cos(ang)

		newEst2[0] = origin[0] + xnew
		newEst2[1] = origin[1] + ynew
		newEst2[2] += ang
		"""
		
		dist = sqrt(newEst2[0] ** 2 + newEst2[1] ** 2)
		vecAng = acos(newEst2[0] / dist)
		if asin(newEst2[1] / dist) < 0:
			vecAng = - vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)], [ - sin(vecAng), cos(vecAng)]])
		foreR = array([[cos(vecAng), - sin(vecAng)], [sin(vecAng), cos(vecAng)]])
		R = array([[cos(newEst2[2]), sin(newEst2[2])], [ - sin(newEst2[2]), cos(newEst2[2])]])
		
		vert = deepcopy(node2.a_vert)

		origin2 = [0.0, 0.0, 0.0]
		finalVec = array([[origin2[0]], [origin2[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0, 0] += dist
		tempVec = dot(foreR, resVec)
		origin2[0] = tempVec[0, 0]
		origin2[1] = tempVec[1, 0]		

		xP = []
		yP = []
		vert2 = []
		for j in range(len(vert)):
			pnt = vert[j]
			
			finalVec = array([[pnt[0]], [pnt[1]]])
			transVec = dot(transpose(R), finalVec)
			resVec = dot(backR, transVec)
			resVec[0, 0] += dist
			tempVec = dot(foreR, resVec)
			pnt[0] = tempVec[0, 0]
			pnt[1] = tempVec[1, 0]		

			xP.append(pnt[0])
			yP.append(pnt[1])

			vert2.append(copy(pnt))

		pylab.plot(xP, yP, color='0.0')

		"""
		newEst3 = copy(est2)
		dist = sqrt(newEst3[0]**2 + newEst3[1]**2)
		vecAng = acos(newEst3[0]/dist)
		if asin(newEst3[1]/dist) < 0:
			vecAng = -vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
		foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])
		R = array([[cos(newEst3[2]), sin(newEst3[2])],[-sin(newEst3[2]),cos(newEst3[2])]])
	
		origin3 = [0.0,0.0,0.0]
		finalVec = array([[origin3[0]],[origin3[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0,0] += dist
		tempVec = dot(foreR, resVec)
		origin3[0] = tempVec[0,0]
		origin3[1] = tempVec[1,0]
		"""

		origin3 = self.circle2[0]
		radius3 = self.circle2[1]
		
		#rvert1 = self.pointsInCircle(origin2, radius, vert1)
		rvert1 = self.pointsInCircle(origin3, radius3, vert1)
		fvert1 = self.filterVertices(vert2, rvert1)
		
		xP = []
		yP = []
		for v in fvert1:
			xP.append(v[0])
			yP.append(v[1])

		#if len(xP) > 0:				
		#	pylab.scatter(xP,yP,linewidth=1, color='r')
		
		origin1 = self.circle1[0]
		radius1 = self.circle1[1]

		rvert2 = self.pointsInCircle(origin1, radius1, vert2)
		fvert2 = self.filterVertices(vert1, rvert2)
		
		xP = []
		yP = []
		for v in fvert2:
			xP.append(v[0])
			yP.append(v[1])
		
		#if len(xP) > 0:
		#	pylab.scatter(xP,yP,linewidth=1, color='b')
		
		c1 = Circle([origin1[0], origin1[1]], radius)
		self.ax.add_artist(c1)
		c1.set_clip_box(self.ax.bbox)
		c1.set_alpha(0.2)
		c1.set_facecolor('r')


		#c2 = Circle([origin2[0],origin2[1]], radius)
		c2 = Circle([origin3[0], origin3[1]], radius)
		self.ax.add_artist(c2)
		c2.set_clip_box(self.ax.bbox)
		c2.set_alpha(0.2)
		c2.set_facecolor('b')

		pylab.scatter([origin1[0]], [origin1[1]], linewidth=3)
		pylab.scatter([origin2[0]], [origin2[1]], linewidth=3)


		#xP.append(xP[0])
		#yP.append(yP[0])
		cost = 0.0
		costCount = 0
		
		for p in fvert1:
			
			mp, mi = self.findClosestPoint(p, vert2)
			xP = [p[0], mp[0]]
			yP = [p[1], mp[1]]
			#pylab.plot(xP,yP, color='0.5', linewidth=1)
			
			dist = sqrt((mp[0] - p[0]) ** 2 + (mp[1] - p[1]) ** 2)
			cost += dist
			#cost += 1.0 - 100.0**-dist
			costCount += 1

		for p in fvert2:
			
			mp, mi = self.findClosestPoint(p, vert1)
			xP = [p[0], mp[0]]
			yP = [p[1], mp[1]]
			#pylab.plot(xP,yP, color='0.5', linewidth=1)

			dist = sqrt((mp[0] - p[0]) ** 2 + (mp[1] - p[1]) ** 2)
			#cost += 1.0 - 100.0**-dist
			cost += dist
			costCount += 1
		
		#print "cost =", cost, "costCount =", costCount
		
		#pylab.xlim(-6, 10)
		#pylab.ylim(-4, 10)
		pylab.xlim(-9, 9)
		pylab.ylim(-8, 8)
		
		
		
		" FIXME:  For the time being, [0,0,0] is assumed the location of joint 19 origin "
		#pylab.scatter([origin1[0]], [origin1[1]], linewidth=1, color='0.0')
		#pylab.scatter([origin2[0]], [origin2[1]], linewidth=1, color='0.0')
		
		" closest points to origins "
		#pnt1, i1 = self.findClosestPoint(origin1, vert1)
		#pnt2, i2 = self.findClosestPoint(origin2, vert2)
		#pylab.scatter([pnt1[0]], [pnt1[1]], linewidth=1, color='0.0')
		#pylab.scatter([pnt2[0]], [pnt2[1]], linewidth=1, color='0.0')
		
		#oPnt1 = self.getNthPoint( len(vert1)/2, i1, vert1 )
		#oPnt2 = self.getNthPoint( len(vert2)/2, i2, vert2 )
		#pylab.scatter([oPnt1[0]], [oPnt1[1]], linewidth=1, color='0.0')
		#pylab.scatter([oPnt2[0]], [oPnt2[1]], linewidth=1, color='0.0')
		
		" 1. filter vertices by a circle with origin on joint 19 and radius pre-selected "
		" 2. return set of vertex chains that are within radius, r1 and r2 "
		" 3a. filter vertices of r1 that are contained in b2, resulting in f1"
		" 3b. filter vertices of r2 that are contained in b1, resulting in f2 "
		" 4a. compute the sum of distances of vertices f1 to any vertex in r2 "
		" 4b. compute the sum of distances of vertices f2 to any vertex in r1 "
		
		
		#pylab.show()
		pylab.savefig("test_%04u.png" % self.count)
		self.count += 1
	
		#print vert1
		#print vert2
		#cost /= float(costCount)
		return cost
	
	def setCircles(self, offset, node1, node2):
		
		est1 = node1.getEstPose()
		est2 = node2.getEstPose()

		" second node is offset from is current position "
		newEst2 = copy(est2)

		finalVec = array([[offset[0]], [offset[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0, 0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		newEst2[0] = tempVec[0, 0]
		newEst2[1] = tempVec[1, 0]
		newEst2[2] = normalizeAngle(est1[2] + offset[2])

		"""
		newEst2 = copy(est2)
		origin = copy(est1)
	
		ang = offset[2]				
		newEst2[0] += offset[0]
		newEst2[1] += offset[1]

		xdiff = newEst2[0] - origin[0]
		ydiff = newEst2[1] - origin[1]

		xnew = xdiff*cos(ang) - ydiff*sin(ang)
		ynew = xdiff*sin(ang) + ydiff*cos(ang)

		newEst2[0] = origin[0] + xnew
		newEst2[1] = origin[1] + ynew
		newEst2[2] += ang
		"""
		
		dist = sqrt(newEst2[0] ** 2 + newEst2[1] ** 2)
		vecAng = acos(newEst2[0] / dist)
		if asin(newEst2[1] / dist) < 0:
			vecAng = - vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)], [ - sin(vecAng), cos(vecAng)]])
		foreR = array([[cos(vecAng), - sin(vecAng)], [sin(vecAng), cos(vecAng)]])
		R = array([[cos(newEst2[2]), sin(newEst2[2])], [ - sin(newEst2[2]), cos(newEst2[2])]])
		

		origin2 = [0.0, 0.0, 0.0]
		finalVec = array([[origin2[0]], [origin2[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0, 0] += dist
		tempVec = dot(foreR, resVec)
		origin2[0] = tempVec[0, 0]
		origin2[1] = tempVec[1, 0]		

		origin1 = [0.0, 0.0, 0.0]
		finalVec = array([[origin1[0]], [origin1[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0, 0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		origin1[0] = tempVec[0, 0]
		origin1[1] = tempVec[1, 0]				
		
		self.circle1 = [origin1, self.radius]
		self.circle2 = [origin2, self.radius]
	
	def filterVertices(self, bVert, rVert):
		" body vertices bVert "
		" target containment vertices "
		
		newVert = []
		for v in rVert:
			if not point_inside_polygon(v[0], v[1], bVert):
				newVert.append(copy(v))
	
		return newVert

	def pointsInCircle(self, origin, radius, vertices):
		
		filteredVertices = []
		
		for v in vertices:
			dist = sqrt((v[0] - origin[0]) ** 2 + (v[1] - origin[1]) ** 2)
			
			if dist < radius:
				filteredVertices.append(copy(v))
	
		return filteredVertices
	
	def getNthPoint(self, n, i, vertices):
		" go clockwise or whatever the default orientation is "

		return vertices[(i + n) % len(vertices)]
	
	

	def findClosestPoint(self, point, vertices):
		
		min = 100000
		minPoint = [0.0, 0.0]
		minIndex = - 1

		for i in range(len(vertices)):
			v = vertices[i]
			dist = sqrt((point[0] - v[0]) ** 2 + (point[1] - v[1]) ** 2)

			if dist < min:
				min = dist
				minPoint = copy(v)
				minIndex = i

		return minPoint, minIndex
	
class MapCost:

	def __init__(self):
		self.index = 0
		self.offsets = []

	def getOffsets(self):
		return self.offsets

	def clearHistory(self):
		self.index = 0
		self.offsets = []
		return

	def cost_func(self, offset, pose1, pose2):

		# 1. pose1 is our fixed reference, so compute the splines for each arc
		# 2. get raw points of pose2, and offset them
		# 3. compute splines with resultant points
		# 4. run the matchSpline cost function on the each of the 3 pairs
		# 5. return the sum of the 3

		"""
		if self.index > 1:
			self.offsets.append(offset)
			return 0.0
		"""
		#pose1.printPosePoints(clr='b')
		#pose2.printPosePoints(clr='r')
		#pylab.savefig("tesat.png")

		#self.index += 1

		self.offsets.append(offset)

		L_spline1 = SplineFit(pose1.getLPose(), kp=3)
		R_spline1 = SplineFit(pose1.getRPose(), kp=3)
		I_spline1 = SplineFit(pose1.getIPose(), kp=3)

		# get raw points for pose 2
		L_pose2 = deepcopy(pose2.getLPose())
		R_pose2 = deepcopy(pose2.getRPose())
		I_pose2 = deepcopy(pose2.getIPose())

		# translate all points by (x,y)
		for pose in [L_pose2, R_pose2, I_pose2]:
			for p in pose:
				p[0] += offset[0]
				p[1] += offset[1]

		# rotate all points with respect to p0 of L_pose2 by offset[2]
		origin = L_pose2[0]
		
		origin = pose1.getRootPose()
		
		ang = offset[2]
		for pose in [L_pose2, R_pose2, I_pose2]:
			for p in pose:
				xdiff = p[0] - origin[0]
				ydiff = p[1] - origin[1]

				xnew = xdiff * cos(ang) - ydiff * sin(ang)
				ynew = xdiff * sin(ang) + ydiff * cos(ang)

				p[0] = origin[0] + xnew
				p[1] = origin[1] + ynew

		# compute splines for the resultant offsets of pose2
		L_spline2 = SplineFit(L_pose2, kp=3)
		R_spline2 = SplineFit(R_pose2, kp=3)
		I_spline2 = SplineFit(I_pose2, kp=3)

		# now compute the cost of each	
		#cost1 = L_spline1.matchSpline(L_spline2, True, self.index)
		#pylab.clf()

		cost1 = L_spline1.matchSpline(L_spline2, plot=False, index=self.index)
		cost2 = R_spline1.matchSpline(R_spline2, plot=False, index=self.index)
		cost3 = I_spline1.matchSpline(I_spline2, plot=False, index=self.index)
		#cost1 = L_spline1.matchSpline(L_spline2)
		#cost2 = R_spline1.matchSpline(R_spline2)
		#cost3 = I_spline1.matchSpline(I_spline2)

		cost = cost1 + cost2 + cost3

		#pylab.title("Total Cost = %f" % cost)
		#pylab.savefig("test_0003/pics/fit%04u.png" % self.index)


		pylab.clf()
		pylab.title("Cost = %f" % cost)
		L_spline1.drawSpline('0.5')
		R_spline1.drawSpline('0.5')
		I_spline1.drawSpline('0.5')

		L_spline2.drawSpline('b')
		R_spline2.drawSpline('b')
		I_spline2.drawSpline('b')
		
		pose1.printPosePoints(clr='0.5')
		pose2.printPosePoints(clr='0.0')
		
		pylab.xlim(-10, 6)
		pylab.ylim(-8, 8)
		pylab.savefig("fit%04u.png" % self.index)
		self.index += 1

		return cost


class MapGraph:

	def __init__(self, probe, contacts):

		self.probe = probe
		self.contacts = contacts
		self.stablePose = StablePose(self.probe)
		
		self.initPose = self.probe.getActualJointPose(19)
		self.poseGraph = graph.graph()
		self.numNodes = 0
		self.currNode = 0

		self.pixelSize = PIXELSIZE
		self.mapSize = 20.0
		self.numPixel = int(2.0 * self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize / 2.0
		self.divPix = floor((2.0 * self.mapSize / self.pixelSize) / self.mapSize)
		self.fileName = "mapGraph%04u.png"
		self.saveCount = 0

		#self.newNode()

	def loadFile(self):
		
		self.saveCount = 0
		self.poseGraph = graph.graph()
		self.numNodes = 0
		self.currNode = 0
		
		#for i in range(0,22):
		#for i in range(0, 14):
		for i in range(0,9):
			self.currNode = LocalNode(self.probe, self.contacts, i, 19, inSim = False)
			self.currNode.readFromFile(i)
			
			self.poseGraph.add_node(i, self.currNode)

			if self.numNodes > 0:
				self.poseGraph.add_edge(i-1, i)

			self.numNodes += 1
	
		#self.saveMap()
		#self.correctPoses()
		#self.saveMap()
		
	def setCenterPoints(self, centerPoints):
		self.currNode.setCenterPoints(centerPoints)
	
	def getCurrentNode(self):
		return self.currNode
	
	def newNode(self):
		
		#self.correctPoses()
		
		if self.numNodes > 0:
			self.currNode.saveToFile()

		#self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19)

		if False:
			if self.numNodes % 2 == 0:
				self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 30)
			else:
				self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 10)

		elif False:
			if self.numNodes % 2 == 1:
				self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 30)
			else:
				self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 10)
		
		else:
			self.currNode = LocalNode(self.probe, self.contacts, self.numNodes, 19)
			

		self.poseGraph.add_node(self.numNodes, self.currNode)

		if self.numNodes > 0:
			self.poseGraph.add_edge(self.numNodes - 1, self.numNodes)

		self.numNodes += 1
		
		#self.contacts.resetPose()
		
		self.stablePose.reset()
		self.stablePose.setNode(self.currNode)
	
	def keepStablePose(self):
		self.stablePose.update()
		
	def correctPoses3(self):
			
		# 2. compute the global position of all the poses wrt to first pose
	
		# TUNE ME:  threshold cost difference between iterations to determine if converged
		#costThresh = 0.004
		costThresh = 0.1
	
		# TUNE ME:   minimum match distance before point is discarded from consideration
		minMatchDist = 2.0
	
		# plot the best fit at each iteration of the algorithm?
		plotIteration = True
		#plotIteration = False
	
		# initial guess for x, y, theta parameters
		offset = [0.0,0.0,0.0]
	
		# sample data
		a_data = []
		b_data = []
	
		estPoses = []
		finalOffsets = []

		for i in range(1, self.numNodes-1):
			
			node1 = self.poseGraph.get_node_attributes(i)
			node2 = self.poseGraph.get_node_attributes(i+1)
			
			print "correcting", i , "and", i+1

			node1.computeAlphaBoundary()
			node2.computeAlphaBoundary()
			
			b_data = node1.getAlphaBoundary()
			a_data = node2.getAlphaBoundary()
			
			def decimatePoints(points):
				result = []
			
				for i in range(len(points)):
					if i%2 == 0:
						result.append(points[i])
			
				return result
			
			a_data = decimatePoints(a_data)
			b_data = decimatePoints(b_data)

			estPose1 = node1.getEstPose()
			estPose2 = node2.getEstPose()
			
			# get the offset of node 2 with respect to node 1
			pose1 = Pose(estPose1)
			offset = pose1.convertGlobalPoseToLocal(estPose2)

			# treat the points with the point-to-line constraint
			gen_icp.addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
			gen_icp.addPointToLineCovariance(b_data, high_var=1.0, low_var=0.001)

			gen_icp.addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
			gen_icp.addDistanceFromOriginCovariance(b_data, tan_var=0.1, perp_var=0.01)

			# plot the data without A transformed, plot 997
			gen_icp.draw(a_data, b_data, "rawData.png")

			# transform the points in A by 'offset'
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispOffset(p, offset))

			# plot the data with A transformed, plot 998
			gen_icp.draw(a_trans, b_data, "initialGuess.png") 

			# run generalized ICP (a plot is made for each iteration of the algorithm)
			offset = gen_icp.gen_ICP(offset, a_data, b_data, costThresh, minMatchDist, plotIteration)

			# transform the points of A with parameters determined by algorithm
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispPoint(p, offset))

			print "drawing: ", "finalOutput_%04u.png" % i

			# plot the final result, plot 999
			gen_icp.draw(a_trans, b_data, "finalOutput_%04u.png" % i, fileWrite = True) 

			finalOffsets.append(offset)
			estPoses.append(estPose1)

		"""
		1. get the relative offsets (done)
		2. for each node, recompute the estimated pose of all the subsequent poses
		3. plot them
		"""
		
		for i in range(len(finalOffsets)):
			
			est1 = estPoses[i]
	
			offset = finalOffsets[i]
	
			pose1 = Pose(est1)
	
			newEst2 = pose1.convertLocalOffsetToGlobal(offset)
	
			if i+1 < len(estPoses):
				estPoses[i+1] = newEst2
			else:
				estPoses.append(newEst2)

		pylab.clf()

		for i in range(len(estPoses)):
			est1 = estPoses[i]
			node1 = self.poseGraph.get_node_attributes(i+1)
			node1.setEstPose(est1)

	
		for i in range(1, self.numNodes-1):

			#est1 = estPoses[i-1]
			node1 = self.poseGraph.get_node_attributes(i)
			est1 = node1.getEstPose()
			
			a_data = node1.estOccPoints
			b_data = []
			a_trans = []
			for p in a_data:
				a_trans.append(gen_icp.dispOffset(p, est1))
				
			gen_icp.draw( a_trans, b_data,"foo.png", fileWrite = False) 
	
	
			#gnd1 = gndPoses[i]
			#a_data = eval(val)
			#b_data = []
			#a_trans = []
			#for p in a_data:
			#	a_trans.append(gen_icp.dispOffset(p, gnd1))
			#gen_icp.draw(b_data, a_trans, "foo.png", fileWrite = False) 
			
		def plotEnv():
			WLEN = 3.0
			wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*math.cos(math.pi/3), -0.2 - WLEN*math.sin(math.pi/3)]]
			wall2 = [[-4.0 + WLEN*math.cos(math.pi/3), 0.2 + WLEN*math.sin(math.pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			w1 = wall1[2]
			w2 = wall2[0]
			
			wall3 = [[w1[0] + 0.4*math.cos(math.pi/6), w1[1] + 0.4*math.sin(math.pi/6)], [0.4*math.cos(math.pi/6) - 4, 0.0], [w2[0] + 0.4*math.cos(math.pi/6), w2[1] - 0.4*math.sin(math.pi/6)], w2]
			lp = wall3[0]
			rp = wall3[2]
			
			wall6 = [lp, [lp[0] + WLEN*math.cos(math.pi/6), lp[1] + WLEN*math.sin(math.pi/6)]]
			wall6.append([wall6[1][0] + 0.4*math.cos(math.pi/3), wall6[1][1] - 0.4*math.sin(math.pi/3)])
			wall6.append([wall6[2][0] - WLEN*math.cos(math.pi/6), wall6[2][1] - WLEN*math.sin(math.pi/6)])
			wall6.append([wall6[3][0] + WLEN*math.cos(math.pi/3), wall6[3][1] - WLEN*math.sin(math.pi/3)])
			wall6.append([wall6[4][0] - 0.4*math.cos(math.pi/6), wall6[4][1] - 0.4*math.sin(math.pi/6)])
			wall6.append(w1)
			wall6.reverse()
		
			walls = [wall1, wall2, wall3, wall6]
		
			for wall in walls:
				xP = []
				yP = []
				for i in range(len(wall)):
					p = copy(wall[i])
					p[0] += 6.0
					wall[i] = p
					xP.append(p[0])
					yP.append(p[1])
		
				pylab.plot(xP,yP, linewidth=2, color = 'g')
		
		plotEnv()
	
		pylab.xlim(-3,6)
		pylab.ylim(-4,4)
	
		gen_icp.save("finalMap.png")	
			
	def correctPoses2(self):
		
		if self.numNodes <= 1:
			return 
		
		" compute alpha shapes for each node "
		for i in range(1, self.numNodes):
			localNode = self.poseGraph.get_node_attributes(i)		
			localNode.computeAlphaBoundary()
				
		initGuess = [0.0, 0.0, 0.0]
		offsets = []
		costs = []

		costObject = AlphaMapCost()

		for i in range(1, self.numNodes - 1):

			node1 = self.poseGraph.get_node_attributes(i)
			node2 = self.poseGraph.get_node_attributes(i + 1)

			pose1 = node1.getEstPose()
			pose2 = node2.getEstPose()

			initGuess = node1.getGlobalPoseToLocal(pose2)
						
			guess1 = fmin(costObject.guessCost, initGuess, args=(node1, node2))
			
			guess2 = costObject.getGuess2(guess1, node1, node2)
			costObject.setCircles(guess2, node1, node2)
			
			offset = fmin(costObject.cost_func, guess2, args=(node1, node2))
			
			cost = costObject.cost_func(offset, node1, node2)
			print "initial guess =", initGuess, "final offset =", offset, "cost =", cost

			offsets.append(offset)

		for i in range(1, self.numNodes - 1):

			node1 = self.poseGraph.get_node_attributes(i)
			node2 = self.poseGraph.get_node_attributes(i + 1)

			est1 = node1.getEstPose()
			
			vert1 = deepcopy(node1.a_vert)
			vert2 = deepcopy(node2.a_vert)

			offset = offsets[i - 1]

			newRoot = node1.convertLocalOffsetToGlobal(offset)
	
			pylab.clf()


			xP = []
			yP = []
			for j in range(len(vert1)):
				pnt = vert1[j]
				
				pnt = node1.convertLocalToGlobal(pnt)
	
				xP.append(pnt[0])
				yP.append(pnt[1])

			pylab.plot(xP, yP, linewidth=1, color='0.5')	

			" node 2 "
			dist = sqrt(newRoot[0] ** 2 + newRoot[1] ** 2)
			vecAng = acos(newRoot[0] / dist)
			if asin(newRoot[1] / dist) < 0:
				vecAng = - vecAng
			
			backR = array([[cos(vecAng), sin(vecAng)], [ - sin(vecAng), cos(vecAng)]])
			foreR = array([[cos(vecAng), - sin(vecAng)], [sin(vecAng), cos(vecAng)]])
			R = array([[cos(newRoot[2]), sin(newRoot[2])], [ - sin(newRoot[2]), cos(newRoot[2])]])

			xP = []
			yP = []
			for j in range(len(vert2)):
				pnt = vert2[j]
				
				finalVec = array([[pnt[0]], [pnt[1]]])
				transVec = dot(transpose(R), finalVec)
				resVec = dot(backR, transVec)
				resVec[0, 0] += dist
				tempVec = dot(foreR, resVec)
				pnt[0] = tempVec[0, 0]
				pnt[1] = tempVec[1, 0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])

			pylab.plot(xP, yP, linewidth=1, color='0.5')	
									
			pylab.xlim(-6, 10)
			pylab.ylim(-4, 10)
		
			pylab.savefig("localFit_%04u.png" % i)

		for i in range(1, self.numNodes - 1):

			node1 = self.poseGraph.get_node_attributes(i)
			node2 = self.poseGraph.get_node_attributes(i + 1)

			offset = offsets[i - 1]
			est1 = node1.getEstPose()
			
			newRoot = copy(node2.getEstPose())
				
			finalVec = array([[offset[0]], [offset[1]]])
			transVec = dot(transpose(node1.R), finalVec)
			resVec = dot(node1.backR, transVec)
			resVec[0, 0] += node1.dist
			tempVec = dot(node1.foreR, resVec)
			newRoot[0] = tempVec[0, 0]
			newRoot[1] = tempVec[1, 0]
			newRoot[2] = normalizeAngle(est1[2] + offset[2])
			node2.setEstPose(newRoot)								
			
		pylab.clf()
		
		for i in range(1, self.numNodes):
			print "i =", i
			
			node1 = self.poseGraph.get_node_attributes(i)

			node = node1
			vert = deepcopy(node.a_vert)

			xP = []
			yP = []
			for j in range(len(vert)):
				pnt = vert[j]
				
				finalVec = array([[pnt[0]], [pnt[1]]])
				transVec = dot(transpose(node.R), finalVec)
				resVec = dot(node.backR, transVec)
				resVec[0, 0] += node.dist
				tempVec = dot(node.foreR, resVec)
				pnt[0] = tempVec[0, 0]
				pnt[1] = tempVec[1, 0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])
			
			clr = str(i * 0.2)
			pylab.plot(xP, yP, linewidth=1, color='0.5')	
			#pylab.plot(xP,yP,linewidth=1,color=clr)	
			
			#f = open("correct_pose%04u.txt" % i, 'w')
			#f.write(repr(node1.getEstPose()))
			#f.close()
		
			pylab.xlim(-6, 10)
			pylab.ylim(-4, 10)
		
			pylab.savefig("globalMap_%04u.png" % i)

		
		#f = open("offsets.txt", 'w')
		#f.write(repr(offsets))
		#f.close()
		
	
	def correctPoses(self):
		
		if self.numNodes <= 1:
			return 
		
		initGuess = [0.0, 0.0, 0.0]
		offsets = []
		costs = []

		costObject = MapCost()

		offsetHistories = []

		for i in range(self.numNodes - 1):
			localNode1 = self.poseGraph.get_node_attributes(i)
			localNode2 = self.poseGraph.get_node_attributes(i + 1)

			pose1 = localNode1.getPoseProfile()
			pose2 = localNode2.getPoseProfile()

			# now let's measure the cost of matching two splines
			# for each pair of splines, search for the lowest cost
			offset = fmin(costObject.cost_func, initGuess, args=(pose1, pose2))
	
			#costObject.clearHistory()

			cost = costObject.cost_func(offset, pose1, pose2)

			costs.append(cost)
			offsets.append(offset)

		print "correction results:"
		print costs
		print offsets

		" Now we perform the corrections with the offset results"
		for i in range(self.numNodes - 1):
			currOffset = offsets[i]
			fixedNode = self.poseGraph.get_node_attributes(i)
			fixedPose = fixedNode.getPoseProfile()
			rootPose = fixedPose.getRootPose()
			print "checkA", i
			for j in range(i + 1, self.numNodes):
				print "checkB", j
				localNode = self.poseGraph.get_node_attributes(j)
				pose = localNode.getPoseProfile()
				print "performing offset of", currOffset, "on pose of", pose.getRootPose()
				pose.performOffset(currOffset, rootPose)
				print "resulting in pose of", pose.getRootPose()
		
		for i in range(self.numNodes):

			localNode = self.poseGraph.get_node_attributes(i)
			pose = localNode.getPoseProfile()
			rootPose = pose.getRootPose()
			oldRootPose = copy(localNode.getEstPose())
			localNode.setEstPose(rootPose)
			print "correcting", i, "from", oldRootPose, "to", rootPose
			

	def synch(self):
		self.currNode.synch()

	def forceUpdate(self, isForward=True):

		self.stablePose.setDirection(isForward)
		self.currNode.update(isForward)

	def update(self, isForward=True):

		self.stablePose.setDirection(isForward)

		if self.stablePose.isStable():
			#print "stable"
			self.currNode.update(isForward)
		else:
			#print "not stable"
			pass
		
		#self.currNode.update(isForward)
		
		#self.saveMap()
	
	def obstCallBack(self, direction):
		self.currNode.obstCallBack(direction)
	
	def saveLocalMap(self):
		
		if self.currNode != 0:
			self.currNode.saveMap()
		
		#for i in range(1, self.numNodes):
		#	localNode = self.poseGraph.get_node_attributes(i)
		#	localNode.saveMap()
			
	def saveMap(self):
		
		" ground truth walls of the environment "
		self.gndMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.gndImage = self.gndMapImage.load()
		gndDraw = ImageDraw.Draw(self.gndMapImage)
		
		walls = self.probe.getWalls()
		for wall in walls:
			wallPoints = []
			for p in wall:
				pGrid = self.realToGrid(p)
				wallPoints.append((pGrid[0],pGrid[1]))
			
			gndDraw.line(wallPoints, fill=255)
		
		self.gndMapImage.save("mapGndGraph%04u.png" % self.saveCount)
		
		#print self.poseGraph
			
		self.mapImage = Image.new('L', (self.numPixel, self.numPixel), 127)
		self.image = self.mapImage.load()

		" 1. for each node, copy local image "
		" 2. for each point in image, transform appropriately and plot to main image "

		for i in range(self.numNodes):
			localNode = self.poseGraph.get_node_attributes(i)
			#localNode.saveMap()
			localOccMap = localNode.getOccMap()
			localMap = localOccMap.getMap()
			localImage = localMap.load()
			
			mapSize = localMap.size
		
		
			for j in range(mapSize[0]):
				for k in range(mapSize[1]):
					if localImage[j, k] == 255:

						pnt = localOccMap.gridToReal([j, k])
						
						pnt = localNode.convertLocalToGlobal(pnt)

						indexX, indexY = self.realToGrid(pnt)
						#if indexX >= 0 and indexX < self.numPixel and indexY >= 0 and indexY < self.numPixel:
						#	self.image[indexX,indexY] = 255
						
						" fill in the interstitial spaces to remove aliasing "
						for m in range(indexX - 1, indexX + 2):
							for n in range(indexY - 1, indexY + 2):
								if m >= 0 and m < self.numPixel and n >= 0 and n < self.numPixel:
									self.image[m, n] = 255

							
		self.mapImage.save(self.fileName % self.saveCount)	
		#self.saveCount += 1	
		
		" build global boundary map "
		self.boundMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.boundImage = self.boundMapImage.load()
		
		for i in range(1,self.numPixel-1):
			for j in range(1,self.numPixel-1):

				isBoundary = False

				# in shadow or unexplored, so it could be a boundary
				if self.image[i,j] <= 127:

					# check if any neighbors are free space
					for m in range(-1,2):
						for n in range(-1,2):
							if self.image[i+m,j+n] == 255:
								isBoundary = True
				
				elif self.image[i,j] == 255:
					self.boundImage[i,j] = 127

				# if a neighbor is free space, then this is a boundary point
				if isBoundary:
					self.boundImage[i,j] = 255		

		self.boundMapImage.save("mapBoundGraph%04u.png" % self.saveCount)	

		" alpha shape boundary of global map "
		self.alphaMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.alphaImage = self.alphaMapImage.load()
		alphaDraw = ImageDraw.Draw(self.alphaMapImage)
		
		alpha_vert = []
		#alpha_vert = self.computeAlphaBoundary()
		
		for i in range(len(alpha_vert)):
			p1 = alpha_vert[i]
			p2 = alpha_vert[(i+1) % len(alpha_vert)]

			p1Grid = self.realToGrid(p1)
			p2Grid = self.realToGrid(p2)

			alphaDraw.line([p1Grid,p2Grid], fill=255)
		
		
		self.alphaMapImage.save("mapAlphaGraph%04u.png" % self.saveCount)
		
		" 1. for each node, copy local image "
		" 2. for each point in image, transform appropriately and plot to main image "
		self.obstMapImage = Image.new('L', (self.numPixel, self.numPixel), 0)
		self.obstImage = self.obstMapImage.load()

		for i in range(self.numNodes):
			localNode = self.poseGraph.get_node_attributes(i)
			#localNode.saveMap()
			localObstMap = localNode.getObstacleMap()
			localMap = localObstMap.getMap()
			localImage = localMap.load()
			
			mapSize = localMap.size
		
			for j in range(mapSize[0]):
				for k in range(mapSize[1]):
					if localImage[j, k] == 255:
						pnt = localObstMap.gridToReal([j, k])

						pnt = localNode.convertLocalToGlobal(pnt)

						indexX, indexY = self.realToGrid(pnt)
						#if indexX >= 0 and indexX < self.numPixel and indexY >= 0 and indexY < self.numPixel:
						#	self.image[indexX,indexY] = 255
						
						" fill in the interstitial spaces to remove aliasing "
						#for m in range(indexX-1,indexX+2):
						#	for n in range(indexY-1,indexY+2):
						#		if m >= 0 and m < self.numPixel and n >= 0 and n < self.numPixel:
						#			self.obstImage[m,n] = 255

						self.obstImage[indexX, indexY] = 255
							
		self.obstMapImage.save("mapObstGraph%04u.png" % self.saveCount)	
		
		self.saveCount += 1	

	def fillOpen(self, polygon):

		# vertices of the 4-sided convex polygon
		polyX = [polygon[0][0], polygon[1][0], polygon[2][0], polygon[3][0], polygon[0][0]]
		polyY = [polygon[0][1], polygon[1][1], polygon[2][1], polygon[3][1], polygon[0][1]]

		# bounding box of polygon
		maxX = - 10000.0
		minX = 10000.0
		maxY = - 10000.0
		minY = 10000.0
		for i in range(4):
			if polyX[i] > maxX:
				maxX = polyX[i]
			if polyX[i] < minX:
				minX = polyX[i]
			if polyY[i] > maxY:
				maxY = polyY[i]
			if polyY[i] < minY:
				minY = polyY[i]

		# set boundaries of map grid
		highIndexX, highIndexY = self.realToGrid([maxX, maxY])
		lowIndexX, lowIndexY = self.realToGrid([minX, minY])

		# fill map grid cell if occupied by arm
		indices = []

		for i in range(lowIndexX, highIndexX + 1, 1):
			for j in range(lowIndexY, highIndexY + 1, 1):
				point = self.gridToReal([i, j])
				if IsContained(polygon, point):
					# free space overwrites everything
					if i >= 0 and i < self.numPixel and j >= 0 and j < self.numPixel:
						self.image[i, j] = 255			
			
	def realToGrid(self, point):
		indexX = int(floor(point[0] * self.divPix)) + self.numPixel / 2 + 1
		indexY = int(floor(point[1] * self.divPix)) + self.numPixel / 2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0, (j - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0]
		return point			

	def computeAlphaBoundary(self):
		
		" 1. pick out the points "
		numPixel = self.numPixel
		mapImage = self.mapImage
		image = mapImage.load()
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					pnt = self.gridToReal([j,k])
					points.append(pnt)

		print len(points), "points"
		
		if len(points) > 0:

			a_vert = self.computeAlpha2(points)
	
			" cut out the repeat vertex "
			a_vert = a_vert[:-1]
			
		else:
			a_vert = []
		
		return a_vert
	
	def computeAlpha2(self, points):
		
		numPoints = len(points)
		inputStr = str(numPoints) + " "
		
		" radius 0.2 "
		inputStr += str(2.0) + " "
		
		for p in points:
			p2 = copy(p)
			p2[0] += gauss(0.0,0.0001)
			p2[1] += gauss(0.0,0.0001)

			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		
		inputStr += "\n"
		
		" start the subprocess "
		subProc = Popen(["./alpha2.exe"], stdin=PIPE, stdout=PIPE)
		
		
		" send input and receive output "
		sout, serr = subProc.communicate(inputStr)
	
		" convert string output to typed data "
		sArr = sout.split(" ")

		numVert = int(sArr[0])
		
		print "numVert =", numVert
		print numVert, len(sArr)
		print sArr
		
		vertices = []
		for i in range(numVert+1):
			vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
		
		return vertices
							