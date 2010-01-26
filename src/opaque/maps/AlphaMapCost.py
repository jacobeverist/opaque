import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from numpy import *
from scipy.optimize import *
import pylab
import graph
import csv
from LocalNode import *
from matplotlib.patches import Circle

class AlphaMapCost:
	
	def __init__(self):
		self.count = 0
		self.fig = pylab.figure()
		self.ax = self.fig.add_subplot(111)
		
		self.radius = 2.5
		
		self.circle1 = [[0.0,0.0,0.0], self.radius]
		self.circle2 = [[0.0,0.0,0.0], self.radius]
	
	def getGuess2(self, offset, node1, node2):
		
		est1 = node1.getEstPose()
		est2 = node2.getEstPose()

		" first node stays in its current position "
		node = node1
		cntPoints = deepcopy(node.centerPoints)

		origin1 = [0.0,0.0,0.0]
		finalVec = array([[origin1[0]],[origin1[1]]])
		transVec = dot(transpose(node.R), finalVec)
		resVec = dot(node.backR, transVec)
		resVec[0,0] += node.dist
		tempVec = dot(node.foreR, resVec)
		origin1[0] = tempVec[0,0]
		origin1[1] = tempVec[1,0]		


		xP = []
		yP = []
		cntPoints1 = []
		for j in range(len(cntPoints)):
			pnt = cntPoints[j]
			
			finalVec = array([[pnt[0]],[pnt[1]]])
			transVec = dot(transpose(node.R), finalVec)
			resVec = dot(node.backR, transVec)
			resVec[0,0] += node.dist
			tempVec = dot(node.foreR, resVec)
			pnt[0] = tempVec[0,0]
			pnt[1] = tempVec[1,0]		

			xP.append(pnt[0])
			yP.append(pnt[1])
			
			cntPoints1.append(copy(pnt))

		spline1 = SplineFit(cntPoints1, kp = 3)
		
		samples = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
		totalVec = [0.0,0.0]
		for i in range(len(samples)):
		
			uvec = spline1.getUVector(samples[i])
			totalVec[0] += uvec[0]
			totalVec[1] += uvec[1]
		
		totalVec[0] /= float(len(samples))
		totalVec[1] /= float(len(samples))
		
		" convert the totalVec into local coordinates "
		globalVec = array([[totalVec[0]],[totalVec[1]]])
		finalVec = dot(node1.R, totalVec)
		print finalVec
		print finalVec[0]
		print finalVec[1]
		
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

		origin1 = [0.0,0.0,0.0]
		finalVec = array([[origin1[0]],[origin1[1]]])
		transVec = dot(transpose(node.R), finalVec)
		resVec = dot(node.backR, transVec)
		resVec[0,0] += node.dist
		tempVec = dot(node.foreR, resVec)
		origin1[0] = tempVec[0,0]
		origin1[1] = tempVec[1,0]		


		xP = []
		yP = []
		cntPoints1 = []
		for j in range(len(cntPoints)):
			pnt = cntPoints[j]
			
			finalVec = array([[pnt[0]],[pnt[1]]])
			transVec = dot(transpose(node.R), finalVec)
			resVec = dot(node.backR, transVec)
			resVec[0,0] += node.dist
			tempVec = dot(node.foreR, resVec)
			pnt[0] = tempVec[0,0]
			pnt[1] = tempVec[1,0]		

			xP.append(pnt[0])
			yP.append(pnt[1])
			
			cntPoints1.append(copy(pnt))

		" second node is offset from is current position "
		newEst2 = copy(est2)

		finalVec = array([[offset[0]],[offset[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0,0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		newEst2[0] = tempVec[0,0]
		newEst2[1] = tempVec[1,0]
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

		dist = sqrt(newEst2[0]**2 + newEst2[1]**2)
		vecAng = acos(newEst2[0]/dist)
		if asin(newEst2[1]/dist) < 0:
			vecAng = -vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
		foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])
		R = array([[cos(newEst2[2]), sin(newEst2[2])],[-sin(newEst2[2]),cos(newEst2[2])]])
		
		cntPoints = deepcopy(node2.centerPoints)

		origin2 = [0.0,0.0,0.0]
		finalVec = array([[origin2[0]],[origin2[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0,0] += dist
		tempVec = dot(foreR, resVec)
		origin2[0] = tempVec[0,0]
		origin2[1] = tempVec[1,0]		

		xP = []
		yP = []
		cntPoints2 = []
		for j in range(len(cntPoints)):
			pnt = cntPoints[j]
			
			finalVec = array([[pnt[0]],[pnt[1]]])
			transVec = dot(transpose(R), finalVec)
			resVec = dot(backR, transVec)
			resVec[0,0] += dist
			tempVec = dot(foreR, resVec)
			pnt[0] = tempVec[0,0]
			pnt[1] = tempVec[1,0]		

			xP.append(pnt[0])
			yP.append(pnt[1])

			cntPoints2.append(copy(pnt))
		
	
		" two splines, transformed from local to global plus offset "
		print cntPoints1
		spline1 = SplineFit(cntPoints1, kp = 3)
		print cntPoints2
		spline2 = SplineFit(cntPoints2, kp = 3)
		
		cost = spline1.matchSpline2(spline2, plot=True, index = self.count)
		
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

		origin1 = [0.0,0.0,0.0]
		finalVec = array([[origin1[0]],[origin1[1]]])
		transVec = dot(transpose(node.R), finalVec)
		resVec = dot(node.backR, transVec)
		resVec[0,0] += node.dist
		tempVec = dot(node.foreR, resVec)
		origin1[0] = tempVec[0,0]
		origin1[1] = tempVec[1,0]		


		xP = []
		yP = []
		vert1 = []
		for j in range(len(vert)):
			pnt = vert[j]
			
			finalVec = array([[pnt[0]],[pnt[1]]])
			transVec = dot(transpose(node.R), finalVec)
			resVec = dot(node.backR, transVec)
			resVec[0,0] += node.dist
			tempVec = dot(node.foreR, resVec)
			pnt[0] = tempVec[0,0]
			pnt[1] = tempVec[1,0]		

			xP.append(pnt[0])
			yP.append(pnt[1])
			
			vert1.append(copy(pnt))
		

		#xP.append(xP[0])
		#yP.append(yP[0])
		
		pylab.plot(xP,yP, color='0.0')


		" second node is offset from is current position "
		newEst2 = copy(est2)

		finalVec = array([[offset[0]],[offset[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0,0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		newEst2[0] = tempVec[0,0]
		newEst2[1] = tempVec[1,0]
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
		
		dist = sqrt(newEst2[0]**2 + newEst2[1]**2)
		vecAng = acos(newEst2[0]/dist)
		if asin(newEst2[1]/dist) < 0:
			vecAng = -vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
		foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])
		R = array([[cos(newEst2[2]), sin(newEst2[2])],[-sin(newEst2[2]),cos(newEst2[2])]])
		
		vert = deepcopy(node2.a_vert)

		origin2 = [0.0,0.0,0.0]
		finalVec = array([[origin2[0]],[origin2[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0,0] += dist
		tempVec = dot(foreR, resVec)
		origin2[0] = tempVec[0,0]
		origin2[1] = tempVec[1,0]		

		xP = []
		yP = []
		vert2 = []
		for j in range(len(vert)):
			pnt = vert[j]
			
			finalVec = array([[pnt[0]],[pnt[1]]])
			transVec = dot(transpose(R), finalVec)
			resVec = dot(backR, transVec)
			resVec[0,0] += dist
			tempVec = dot(foreR, resVec)
			pnt[0] = tempVec[0,0]
			pnt[1] = tempVec[1,0]		

			xP.append(pnt[0])
			yP.append(pnt[1])

			vert2.append(copy(pnt))

		pylab.plot(xP,yP, color='0.0')

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
		
		c1 = Circle([origin1[0],origin1[1]], radius)
		self.ax.add_artist(c1)
		c1.set_clip_box(self.ax.bbox)
		c1.set_alpha(0.2)
		c1.set_facecolor('r')


		#c2 = Circle([origin2[0],origin2[1]], radius)
		c2 = Circle([origin3[0],origin3[1]], radius)
		self.ax.add_artist(c2)
		c2.set_clip_box(self.ax.bbox)
		c2.set_alpha(0.2)
		c2.set_facecolor('b')

		pylab.scatter([origin1[0]],[origin1[1]],linewidth=3)
		pylab.scatter([origin2[0]],[origin2[1]], linewidth=3)


		#xP.append(xP[0])
		#yP.append(yP[0])
		cost = 0.0
		costCount = 0
		
		for p in fvert1:
			
			mp, mi = self.findClosestPoint(p, vert2)
			xP = [p[0],mp[0]]
			yP = [p[1],mp[1]]
			#pylab.plot(xP,yP, color='0.5', linewidth=1)
			
			dist = sqrt((mp[0]-p[0])**2 + (mp[1]-p[1])**2)
			cost += dist
			#cost += 1.0 - 100.0**-dist
			costCount += 1

		for p in fvert2:
			
			mp, mi = self.findClosestPoint(p, vert1)
			xP = [p[0],mp[0]]
			yP = [p[1],mp[1]]
			#pylab.plot(xP,yP, color='0.5', linewidth=1)

			dist = sqrt((mp[0]-p[0])**2 + (mp[1]-p[1])**2)
			#cost += 1.0 - 100.0**-dist
			cost += dist
			costCount += 1
		
		#print "cost =", cost, "costCount =", costCount
		
		pylab.xlim(-6,10)
		pylab.ylim(-4,10)
		
		
		
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

		finalVec = array([[offset[0]],[offset[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0,0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		newEst2[0] = tempVec[0,0]
		newEst2[1] = tempVec[1,0]
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
		
		dist = sqrt(newEst2[0]**2 + newEst2[1]**2)
		vecAng = acos(newEst2[0]/dist)
		if asin(newEst2[1]/dist) < 0:
			vecAng = -vecAng
		
		backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
		foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])
		R = array([[cos(newEst2[2]), sin(newEst2[2])],[-sin(newEst2[2]),cos(newEst2[2])]])
		

		origin2 = [0.0,0.0,0.0]
		finalVec = array([[origin2[0]],[origin2[1]]])
		transVec = dot(transpose(R), finalVec)
		resVec = dot(backR, transVec)
		resVec[0,0] += dist
		tempVec = dot(foreR, resVec)
		origin2[0] = tempVec[0,0]
		origin2[1] = tempVec[1,0]		

		origin1 = [0.0,0.0,0.0]
		finalVec = array([[origin1[0]],[origin1[1]]])
		transVec = dot(transpose(node1.R), finalVec)
		resVec = dot(node1.backR, transVec)
		resVec[0,0] += node1.dist
		tempVec = dot(node1.foreR, resVec)
		origin1[0] = tempVec[0,0]
		origin1[1] = tempVec[1,0]				
		
		self.circle1 = [origin1, self.radius]
		self.circle2 = [origin2, self.radius]
	
	def filterVertices(self, bVert, rVert):
		" body vertices bVert "
		" target containment vertices "
		
		newVert = []
		for v in rVert:
			if not point_inside_polygon(v[0],v[1],bVert):
				newVert.append(copy(v))
	
		return newVert

	def pointsInCircle(self, origin, radius, vertices):
		
		filteredVertices = []
		
		for v in vertices:
			dist = sqrt((v[0]-origin[0])**2 + (v[1]-origin[1])**2)
			
			if dist < radius:
				filteredVertices.append(copy(v))
	
		return filteredVertices
	
	def getNthPoint(self, n, i, vertices):
		" go clockwise or whatever the default orientation is "

		return vertices[(i + n) % len(vertices)]
	
	

	def findClosestPoint(self, point, vertices):
		
		min = 100000
		minPoint = [0.0,0.0]
		minIndex = -1

		for i in range(len(vertices)):
			v = vertices[i]
			dist = sqrt((point[0]-v[0])**2 + (point[1]-v[1])**2)

			if dist < min:
				min = dist
				minPoint = copy(v)
				minIndex = i

		return minPoint, minIndex
		

class TestNode:
	
	def __init__(self, nodeID, isCorrect = False):
		
		self.nodeID = nodeID
		
		if not isCorrect:
			" read in the data "
			f = open("estpose%04u.txt" % self.nodeID, 'r')
			self.estPose = eval(f.read())
			f.close()
		else:
			" read in the data "
			f = open("correct_pose%04u.txt" % self.nodeID, 'r')
			self.estPose = eval(f.read())
			f.close()

		" read in the data "
		f = open("cntpose%04u.txt" % self.nodeID, 'r')
		self.centerPoints = eval(f.read())
		f.close()
						
		self.mapImage = Image.open("localOccMap%03u_0000.png" % self.nodeID)
		self.image = self.mapImage.load()
		#self.numPixel = self.mapImage.size[0]

		self.pixelSize = PIXELSIZE
		self.mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)
		
		self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)
		self.vecAng = acos(self.estPose[0]/self.dist)
		if asin(self.estPose[1]/self.dist) < 0:
			self.vecAng = -self.vecAng
		
		self.backR = array([[cos(self.vecAng), sin(self.vecAng)],[-sin(self.vecAng),cos(self.vecAng)]])
		self.foreR = array([[cos(self.vecAng), -sin(self.vecAng)],[sin(self.vecAng),cos(self.vecAng)]])
		
		self.R = array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])

		" 1. pick out the points "
		points = []
		for j in range(self.numPixel):
			for k in range(self.numPixel):
				if self.image[j,k] == 255:
					pnt = self.gridToReal([j,k])
					points.append(pnt)
					
		self.a_vert = self.computeAlpha2(points)
		#self.a_vert.reverse()

		" cut out the repeat vertex "
		self.a_vert = self.a_vert[:-1]
		
		self.convertAlphaUniform()
		
	def convertAlphaUniform(self):
		
		" make the vertices uniformly distributed "
		
		new_vert = []
		
		#max_spacing = 0.04
		max_spacing = 0.1
		
		for i in range(len(self.a_vert)):
			p0 = self.a_vert[i]
			p1 = self.a_vert[(i+1) % len(self.a_vert)]
			dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

			vec = [p1[0]-p0[0], p1[1]-p0[1]]
			vec[0] /= dist
			vec[1] /= dist
			
			new_vert.append(copy(p0))
			
			if dist > max_spacing:
				" cut into pieces max_spacing length or less "
				numCount = int(floor(dist / max_spacing))
				
				for j in range(1, numCount+1):
					newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
					new_vert.append(newP)
					
		self.a_vert = new_vert
	
	def computeAlpha2(self, points):
		
		numPoints = len(points)
		inputStr = str(numPoints) + " "
		
		" radius 0.2 "
		inputStr += str(0.2) + " "
		
		for p in points:
			p2 = copy(p)
			p2[0] += gauss(0.0,0.01)
			p2[1] += gauss(0.0,0.01)

			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		
		inputStr += "\n"
		
		" start the subprocess "
		subProc = Popen(["../../alpha2.exe"], stdin=PIPE, stdout=PIPE)
		
		" send input and receive output "
		sout, serr = subProc.communicate(inputStr)

		#print numPoints
		#print sout
		
		" convert string output to typed data "
		sArr = sout.split(" ")
		
		#print sArr[0]
		numVert = int(sArr[0])
		
		vertices = []
		for i in range(numVert+1):
			vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
		
		#print vertices
		return vertices
	
	def setEstPose(self, estPose):
		self.estPose = copy(estPose)

		self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)
		self.vecAng = acos(self.estPose[0]/self.dist)
		if asin(self.estPose[1]/self.dist) < 0:
			self.vecAng = -self.vecAng
		
		self.backR = array([[cos(self.vecAng), sin(self.vecAng)],[-sin(self.vecAng),cos(self.vecAng)]])
		self.foreR = array([[cos(self.vecAng), -sin(self.vecAng)],[sin(self.vecAng),cos(self.vecAng)]])
		
		self.R = array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])

	
	def getEstPose(self):
		return self.estPose

	def getMap(self):
		return self.mapImage


	def realToGrid(self, point):
		indexX = int(floor(point[0]*self.divPix)) + self.numPixel/2 + 1
		indexY = int(floor(point[1]*self.divPix)) + self.numPixel/2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0, (j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]
		return point
	
class AlphaMap:
	
	def __init__(self):
		self.numNodes = 9
		self.nodes = []
		for i in range(self.numNodes):
			self.nodes.append(TestNode(i+1))
		#self.nodes = [TestNode(4), TestNode(5)]
		#self.numNodes = 2
		
		self.correctPoses()
		
		self.drawCount = 0
	
		#self.draw()
		
	def draw(self):
		
		f = open("offsets.txt", 'r')
		offsets = eval(f.read())
		f.close()

		pylab.clf()

		for i in range(0,self.numNodes-1):

			node1 = self.nodes[i]
			node2 = self.nodes[i+1]
			
			est1 = node1.getEstPose()
			
			vert1 = deepcopy(node1.a_vert)
			vert2 = deepcopy(node2.a_vert)

			offset = offsets[i]
			origin = node1.getEstPose()		

			" second node is offset from is current position "
			newRoot = copy(node2.getEstPose())
	
			finalVec = array([[offset[0]],[offset[1]]])
			transVec = dot(transpose(node1.R), finalVec)
			resVec = dot(node1.backR, transVec)
			resVec[0,0] += node1.dist
			tempVec = dot(node1.foreR, resVec)
			newRoot[0] = tempVec[0,0]
			newRoot[1] = tempVec[1,0]
			newRoot[2] = normalizeAngle(est1[2] + offset[2])
		
			"""
			ang = offset[2]
			
			newRoot = copy(node2.getEstPose())
			
			newRoot[0] += offset[0]
			newRoot[1] += offset[1]

			xdiff = newRoot[0] - origin[0]
			ydiff = newRoot[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			newRoot[0] = origin[0] + xnew
			newRoot[1] = origin[1] + ynew
			newRoot[2] += ang
			"""
			
			
			
			pylab.clf()


			xP = []
			yP = []
			for j in range(len(vert1)):
				pnt = vert1[j]
				
				finalVec = array([[pnt[0]],[pnt[1]]])
				transVec = dot(transpose(node1.R), finalVec)
				resVec = dot(node1.backR, transVec)
				resVec[0,0] += node1.dist
				tempVec = dot(node1.foreR, resVec)
				pnt[0] = tempVec[0,0]
				pnt[1] = tempVec[1,0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])

			pylab.plot(xP,yP,linewidth=1,color='0.5')	

			" node 2 "
			dist = sqrt(newRoot[0]**2 + newRoot[1]**2)
			vecAng = acos(newRoot[0]/dist)
			if asin(newRoot[1]/dist) < 0:
				vecAng = -vecAng
			
			backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
			foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])
			R = array([[cos(newRoot[2]), sin(newRoot[2])],[-sin(newRoot[2]),cos(newRoot[2])]])

			xP = []
			yP = []
			for j in range(len(vert2)):
				pnt = vert2[j]
				
				finalVec = array([[pnt[0]],[pnt[1]]])
				transVec = dot(transpose(R), finalVec)
				resVec = dot(backR, transVec)
				resVec[0,0] += dist
				tempVec = dot(foreR, resVec)
				pnt[0] = tempVec[0,0]
				pnt[1] = tempVec[1,0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])

			pylab.plot(xP,yP,linewidth=1,color='0.5')	
									
			pylab.xlim(-6,10)
			pylab.ylim(-4,10)
		
			pylab.savefig("localFit_%04u.png" % i)

		pylab.clf()
				
		for k in range(self.numNodes):
			node = self.nodes[k]
			vert = deepcopy(node.a_vert)

			xP = []
			yP = []
			for m in range(len(vert)):
				pnt = vert[m]
				
				finalVec = array([[pnt[0]],[pnt[1]]])
				transVec = dot(transpose(node.R), finalVec)
				resVec = dot(node.backR, transVec)
				resVec[0,0] += node.dist
				tempVec = dot(node.foreR, resVec)
				pnt[0] = tempVec[0,0]
				pnt[1] = tempVec[1,0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])
				
			pylab.plot(xP,yP,linewidth=1,color='0.5')	
			
			#f = open("correct_pose%04u.txt" % self.drawCount, 'w')
			#f.write(repr(self.nodes[k].getEstPose()))
			#f.close()
		
		pylab.xlim(-6,10)
		pylab.ylim(-4,10)
	
		pylab.savefig("globalMap_%04u.png" % self.drawCount)

		self.drawCount += 1		
		for i in range(0,self.numNodes-1):

			offset = offsets[i]
			est1 = self.nodes[i].getEstPose()
			node1 = self.nodes[i]
			
			newRoot = copy(self.nodes[i+1].getEstPose())
				
			finalVec = array([[offset[0]],[offset[1]]])
			transVec = dot(transpose(node1.R), finalVec)
			resVec = dot(node1.backR, transVec)
			resVec[0,0] += node1.dist
			tempVec = dot(node1.foreR, resVec)
			newRoot[0] = tempVec[0,0]
			newRoot[1] = tempVec[1,0]
			newRoot[2] = normalizeAngle(est1[2] + offset[2])
			self.nodes[i+1].setEstPose(newRoot)								
			
			"""
			for j in range(i+1, self.numNodes):
				
				newRoot = copy(self.nodes[j].getEstPose())
				
				newRoot[0] += offset[0]
				newRoot[1] += offset[1]
	
				xdiff = newRoot[0] - origin[0]
				ydiff = newRoot[1] - origin[1]
	
				xnew = xdiff*cos(ang) - ydiff*sin(ang)
				ynew = xdiff*sin(ang) + ydiff*cos(ang)
	
				newRoot[0] = origin[0] + xnew
				newRoot[1] = origin[1] + ynew
				newRoot[2] += ang
				self.nodes[j].setEstPose(newRoot)
				"""
				
			
			pylab.clf()
			
			for k in range(self.numNodes):
				node = self.nodes[k]
				vert = deepcopy(node.a_vert)
	
				xP = []
				yP = []
				for m in range(len(vert)):
					pnt = vert[m]
					
					finalVec = array([[pnt[0]],[pnt[1]]])
					transVec = dot(transpose(node.R), finalVec)
					resVec = dot(node.backR, transVec)
					resVec[0,0] += node.dist
					tempVec = dot(node.foreR, resVec)
					pnt[0] = tempVec[0,0]
					pnt[1] = tempVec[1,0]		
		
					xP.append(pnt[0])
					yP.append(pnt[1])
					
				pylab.plot(xP,yP,linewidth=1,color='0.5')	
				
				#f = open("correct_pose%04u.txt" % self.drawCount, 'w')
				#f.write(repr(self.nodes[k].getEstPose()))
				#f.close()
			
			pylab.xlim(-6,10)
			pylab.ylim(-4,10)
		
			pylab.savefig("globalMap_%04u.png" % self.drawCount)

			self.drawCount += 1
	
	def correctPoses(self):
		
		if self.numNodes <= 1:
			return 
		
		initGuess = [0.0, 0.0, 0.0]
		offsets = []
		costs = []

		costObject = AlphaMapCost()

		offsetHistories = []

		for i in range(self.numNodes-1):
			localMap1 = self.nodes[i].getMap()
			localMap2 = self.nodes[i+1].getMap()

			pose1 = self.nodes[i].getEstPose()
			pose2 = self.nodes[i+1].getEstPose()

			node1 = self.nodes[i]

			adiff = normalizeAngle(pose2[2] - pose1[2])
			
			globalVec = array([[pose2[0]],[pose2[1]]])
			tempVec = dot(node1.backR, globalVec)
			tempVec[0,0] -= node1.dist
			transVec = dot(node1.foreR, tempVec)
			finalVec = dot(node1.R, transVec)
			initGuess = [finalVec[0,0], finalVec[1,0], adiff]
			
			guess1 = fmin(costObject.guessCost,initGuess,args=(self.nodes[i], self.nodes[i+1]))
			
			guess2 = costObject.getGuess2(guess1, self.nodes[i], self.nodes[i+1])
			costObject.setCircles(guess2, self.nodes[i], self.nodes[i+1])
			
			offset = fmin(costObject.cost_func,guess2,args=(self.nodes[i], self.nodes[i+1]))
			
			cost = costObject.cost_func(offset, self.nodes[i], self.nodes[i+1])
			print "initial guess =", initGuess, "final offset =", offset, "cost =", cost

			offsets.append(offset)

		for i in range(0,self.numNodes-1):

			node1 = self.nodes[i]
			node2 = self.nodes[i+1]
			
			est1 = node1.getEstPose()
			
			vert1 = deepcopy(node1.a_vert)
			vert2 = deepcopy(node2.a_vert)

			offset = offsets[i]
			origin = self.nodes[i].getEstPose()		

			newRoot = copy(node2.getEstPose())
			finalVec = array([[offset[0]],[offset[1]]])
			transVec = dot(transpose(node1.R), finalVec)
			resVec = dot(node1.backR, transVec)
			resVec[0,0] += node1.dist
			tempVec = dot(node1.foreR, resVec)
			newRoot[0] = tempVec[0,0]
			newRoot[1] = tempVec[1,0]
			newRoot[2] = normalizeAngle(est1[2] + offset[2])

			pylab.clf()


			xP = []
			yP = []
			for j in range(len(vert1)):
				pnt = vert1[j]
				
				finalVec = array([[pnt[0]],[pnt[1]]])
				transVec = dot(transpose(node1.R), finalVec)
				resVec = dot(node1.backR, transVec)
				resVec[0,0] += node1.dist
				tempVec = dot(node1.foreR, resVec)
				pnt[0] = tempVec[0,0]
				pnt[1] = tempVec[1,0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])

			pylab.plot(xP,yP,linewidth=1,color='0.5')	

			" node 2 "
			dist = sqrt(newRoot[0]**2 + newRoot[1]**2)
			vecAng = acos(newRoot[0]/dist)
			if asin(newRoot[1]/dist) < 0:
				vecAng = -vecAng
			
			backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
			foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])
			R = array([[cos(newRoot[2]), sin(newRoot[2])],[-sin(newRoot[2]),cos(newRoot[2])]])

			xP = []
			yP = []
			for j in range(len(vert2)):
				pnt = vert2[j]
				
				finalVec = array([[pnt[0]],[pnt[1]]])
				transVec = dot(transpose(R), finalVec)
				resVec = dot(backR, transVec)
				resVec[0,0] += dist
				tempVec = dot(foreR, resVec)
				pnt[0] = tempVec[0,0]
				pnt[1] = tempVec[1,0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])

			pylab.plot(xP,yP,linewidth=1,color='0.5')	
									
			pylab.xlim(-6,10)
			pylab.ylim(-4,10)
		
			pylab.savefig("localFit_%04u.png" % i)

		for i in range(0,self.numNodes-1):

			offset = offsets[i]
			est1 = self.nodes[i].getEstPose()
			node1 = self.nodes[i]
			
			newRoot = copy(self.nodes[i+1].getEstPose())
				
			finalVec = array([[offset[0]],[offset[1]]])
			transVec = dot(transpose(node1.R), finalVec)
			resVec = dot(node1.backR, transVec)
			resVec[0,0] += node1.dist
			tempVec = dot(node1.foreR, resVec)
			newRoot[0] = tempVec[0,0]
			newRoot[1] = tempVec[1,0]
			newRoot[2] = normalizeAngle(est1[2] + offset[2])
			self.nodes[i+1].setEstPose(newRoot)								
			
		pylab.clf()
		
		for i in range(self.numNodes):
			print "i =", i
			node = self.nodes[i]
			vert = deepcopy(node.a_vert)

			xP = []
			yP = []
			for j in range(len(vert)):
				pnt = vert[j]
				
				finalVec = array([[pnt[0]],[pnt[1]]])
				transVec = dot(transpose(node.R), finalVec)
				resVec = dot(node.backR, transVec)
				resVec[0,0] += node.dist
				tempVec = dot(node.foreR, resVec)
				pnt[0] = tempVec[0,0]
				pnt[1] = tempVec[1,0]		
	
				xP.append(pnt[0])
				yP.append(pnt[1])
			
			clr = str(i * 0.2)
			pylab.plot(xP,yP,linewidth=1,color='0.5')	
			#pylab.plot(xP,yP,linewidth=1,color=clr)	
			
			#f = open("correct_pose%04u.txt" % i, 'w')
			#f.write(repr(self.nodes[i].getEstPose()))
			#f.close()
		
			pylab.xlim(-6,10)
			pylab.ylim(-4,10)
		
			pylab.savefig("globalMap_%04u.png" % i)

		
		#f = open("offsets.txt", 'w')
		#f.write(repr(offsets))
		#f.close()
		

class GlobalMap:
	def __init__(self):
		
		self.numNodes = 13
		self.nodes = []
		for i in range(self.numNodes):
			self.nodes.append(TestNode(i, True))
		#self.nodes = [TestNode(1), TestNode(2)]
		
		self.correctPoses()
		pass