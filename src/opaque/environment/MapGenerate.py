#!/usr/bin/python

from math import *
from random import *
from scipy.optimize import *
import scipy.interpolate
import numpy
import pylab
import getopt, sys, os, shutil, csv
from copy import *
import time
#from common import *
import graph

#from os import *

mapCount = 0

# determines signed area of 3 points (used for solving Point in Polygon problem)
def Area2(Ax,Ay,Bx,By,Cx,Cy):
	return (Bx - Ax) * (Cy - Ay) - (Cx - Ax)*(By - Ay)

# determines if point is left of line segment
def LeftOnEdge(edge, point):
	return (Area2(edge[0][0],edge[0][1],edge[1][0],edge[1][1],point[0],point[1]) >= 0)

# determines if point C is left of line segment AB
def LeftOn(Ax,Ay,Bx,By,Cx,Cy):
	return (Area2(Ax,Ay,Bx,By,Cx,Cy) >= 0)

def segIntersection(seg1, seg2):

	p0 = seg1[0]
	p1 = seg1[1]
	p2 = seg2[0]
	p3 = seg2[1]

	test1 = LeftOn(p0[0],p0[1],p1[0],p1[1],p2[0],p2[1])
	test2 = LeftOn(p0[0],p0[1],p1[0],p1[1],p3[0],p3[1])

	if test1 == test2:
		return False

	test3 = LeftOn(p2[0],p2[1],p3[0],p3[1],p0[0],p0[1])
	test4 = LeftOn(p2[0],p2[1],p3[0],p3[1],p1[0],p1[1])

	if test3 == test4:
		return False

	return True


def Intersect(seg1, seg2):

	denom = (seg2[1][1] - seg2[0][1])*(seg1[1][0] - seg1[0][0]) - (seg2[1][0] - seg2[0][0])*(seg1[1][1] - seg1[0][1])

	nume_a = (seg2[1][0] - seg2[0][0])*(seg1[0][1] - seg2[0][1]) - (seg2[1][1] - seg2[0][1])*(seg1[0][0] - seg2[0][0])

	nume_b = (seg1[1][0] - seg1[0][0])*(seg1[0][1] - seg2[0][1]) - (seg1[1][1] - seg1[0][1])*(seg1[0][0] - seg2[0][0])

	intersection = [0.0,0.0]

	if denom == 0.0:
		if nume_a == 0.0 and nume_b == 0.0:
			return False, intersection
		return False, intersection

	ua = nume_a / denom
	ub = nume_b / denom

	if ua >= 0.0 and ua <= 1.0 and ub >= 0.0 and ub <= 1.0:
		# Get the intersection point.
		intersection[0] = seg1[0][0] + ua*(seg1[1][0] - seg1[0][0])
		intersection[1] = seg1[0][1] + ua*(seg1[1][1] - seg1[0][1])
		return True, intersection

	return False, intersection


def closestSegPoint(seg, tPoint):

	try:
		# distance to each of the end points
		p0 = seg[0]
		p1 = seg[1]

		if tPoint[0] == p1[0] and tPoint[1] == p1[1]:
			return 0.0, deepcopy(p1)

		if tPoint[0] == p0[0] and tPoint[1] == p0[1]:
			return 0.0, deepcopy(p0)


		# vec1 from p0 to p1
		vec1 = [p1[0]-p0[0], p1[1]-p0[1]]

		# vec2 from p1 to p0
		vec2 = [p0[0]-p1[0], p0[1]-p1[1]]

		# vector from seg point to tPoint
		vecA = [tPoint[0]-p0[0], tPoint[1]-p0[1]]
		vecB = [tPoint[0]-p1[0], tPoint[1]-p1[1]]

		dist0 = sqrt((p0[0]-tPoint[0])**2 + (p0[1]-tPoint[1])**2)
		dist1 = sqrt((p1[0]-tPoint[0])**2 + (p1[1]-tPoint[1])**2)

		dotP1 = vecA[0] * vec1[0] + vecA[1] * vec1[1]
		magA = sqrt(vecA[0]*vecA[0] + vecA[1]*vecA[1])
		mag1 = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1])
		val = dotP1/(magA*mag1)
		if val >= 1.0:
			val = 1.0
			
		if val <= -1.0:
			val = -1.0
		a0 = acos(val)

		dotP2 = vecB[0] * vec2[0] + vecB[1] * vec2[1]
		magB = sqrt(vecB[0]*vecB[0] + vecB[1]*vecB[1])
		mag2 = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1])
		
		val = dotP2/(magB*mag2)
		if val >= 1.0:
			val = 1.0
			
		if val <= -1.0:
			val = -1.0
		a1 = acos(val)
		#a1 = acos(dotP2/(magB*mag2))
		
		# obtuse angle, so p0 is closest
		if a0 > pi/2:
			return dist0, deepcopy(p0)
			
		# obtuse angle, so p1 is closest
		if a1 > pi/2:
			return dist1, deepcopy(p1)

		segLength = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)

		# derived Aug, 4, 2008, Jacob's journal
		# otherwise, the closest point is determined by perpedicular line to segment
		v = dist0 
		w = dist1
		d = segLength

		c = (d**2 - v**2 + w**2)/(2*d)
		a = d - c

		# this is the minimum distance to the segment
		val = v**2 - a**2 
		if val < 0.0:
			val = 0.0
		b = sqrt(val)

		# now let's determine the intersection point

		# normalized vec1
		mag1 = sqrt(vec1[0]**2 + vec1[1]**2)
		n1 = [vec1[0]/mag1, vec1[1]/mag1]

		mag2 = sqrt(vec2[0]**2 + vec2[1]**2)
		n2 = [vec2[0]/mag2, vec2[1]/mag2]

		bVec = [a*n1[0], a*n1[1]]

		resultP = [bVec[0] + p0[0], bVec[1] + p0[1]] 

		#print vec1, n1, bVec, resultP

		return b, resultP

	except:

		# find the point on a segment that is closest to "tPoint"
		print "failed finding closest point of seg", seg, "to tPoint", tPoint

		raise


class MapMaker:

	def __init__(self, entr_len, pipe_width):
	
		self.activeSegments = {}
		self.doneSegments = {}

		self.entrLen = entr_len
		self.pipeWidth = pipe_width
		self.wallSoup = []
		self.segCount = 0

		p1 = [0.0,-self.pipeWidth/2]
		p2 = [0.0,self.pipeWidth/2]

		self.wallSoup.append([p1,p2])

		self.walls = []
		self.walls.append([p1,p2])

		#self.walls = {}
		#self.walls[-1] = [p1,p2]

		currAng = 0.0

		# [wallEnd1, wallEnd2, wallAngle]
		self.addSeg(p1,p2,currAng)

		self.addStraight(0,entr_len)

	def done(self):
		for k, v in self.activeSegments.iteritems():
			self.delSeg(k)


	def getSegs(self):
		segIDs = []
		for k, v in self.activeSegments.iteritems():
			segIDs.append(k)

		return segIDs
			
	def addLeftWall(self, segNum, p):
		self.activeSegments[segNum][3].append(p)

	def addRightWall(self, segNum, p):
		self.activeSegments[segNum][4].append(p)

	def updateSeg(self, segNum, p1, p2, currAng):
		self.activeSegments[segNum][0] = p1
		self.activeSegments[segNum][1] = p2
		self.activeSegments[segNum][2] = currAng

	def addSeg(self, p1, p2, currAng):

		segNum = self.segCount
		self.activeSegments[segNum] = [p1, p2, currAng, [p1], [p2]]
		self.segCount += 1
	
		return segNum

	def delSeg(self, segNum):
		self.doneSegments[segNum] = self.activeSegments[segNum]
		del self.activeSegments[segNum]

		self.walls.append(self.doneSegments[segNum][3])
		self.walls.append(self.doneSegments[segNum][4])

	def addCap(self, segNum):

		p1, p2, currAng, wallL, wallR = self.activeSegments[segNum]

		wall = [p1, p2]
		self.addLeftWall(segNum, p2)
		self.wallSoup.append(wall)

		self.delSeg(segNum)

	def addSymJunction(self, segNum, angle):

		turnAngleA = angle
		turnAngleB = -angle

		p1, p2, currAng, wallL, wallR = self.activeSegments[segNum]

		#print "addSymJunction(", p1, p2, currAng, angle

		" center line point between p1 and p2 "
		pC = [(p2[0]+p1[0])/2.0, (p2[1]+p1[1])/2.0]
		angC = currAng


		wallA1 = []
		wallA2 = []

		wallB1 = []
		wallB2 = []

		" SIDE A "

		" turnAngleA "
		currAngA = currAng
	
		" angle of right triangle at the junction "
		pivotAngle = (pi - turnAngleA)/2.0

		" hypotenuse and x side "
		hyp = self.pipeWidth / sin(pivotAngle)
		x = hyp * cos(pivotAngle)
		
		" compensating distance of one side "
		newP1 = [p1[0] + x*cos(currAngA), p1[1] + x*sin(currAngA)]

		currAngA += turnAngleA

		" side continued in junction direction "
		newP1 = [newP1[0] + x*cos(currAngA), newP1[1] + x*sin(currAngA)]

		" two simultaneous equations to find opposing corner point and angle "
		" C is centerline vector and Q is side point and vector "
		pQ = newP1
		angQ = currAngA
		Q = (pC[1] + tan(angC) * (pQ[0]-pC[0]) - pQ[1]) / ( sin(angQ) - cos(angQ) * tan(angC) )

		newPX1 = [newP1[0] + Q*cos(currAngA), newP1[1] + Q*sin(currAngA)]

		wallA1.append(newPX1)
		wallA2.append(p2)

		if angle <= pi/3:

			newPA2 = [newPX1[0] + self.pipeWidth*cos(currAngA+pi/2), newPX1[1] + self.pipeWidth*sin(currAngA+pi/2)]
			wallA2.append(newPA2)

			self.addRightWall(segNum, newPA2)
			self.wallSoup.append(wallA2)
		
			segA = self.addSeg(newPX1,newPA2, currAngA)

		else:
			newPA1 = [p2[0] + self.pipeWidth*cos(currAngA-pi/2), p2[1] + self.pipeWidth*sin(currAngA-pi/2)]
			wallA1.append(newPA1)

			self.addLeftWall(segNum, newPA1)
			self.wallSoup.append(wallA1)
		
			segA = self.addSeg(newPA1,p2, currAngA)

		"  SIDE B "

		" turnAngleB "
		currAngB = currAng

		pivotAngle = (pi + turnAngleB)/2.0
		
		hyp = self.pipeWidth / sin(pivotAngle)
		x = hyp * cos(pivotAngle)
		
		newP2 = [p2[0] + x*cos(currAngB), p2[1] + x*sin(currAngB)]

		currAngB += turnAngleB

		newP2 = [newP2[0] + x*cos(currAngB), newP2[1] + x*sin(currAngB)]

		pQ = newP2
		angQ = currAngB
		Q = (pC[1] + tan(angC) * (pQ[0]-pC[0]) - pQ[1]) / ( sin(angQ) - cos(angQ) * tan(angC) )

		newPX2 = [newP1[0] + Q*cos(currAngA), newP1[1] + Q*sin(currAngA)]

		wallB1.append(p1)
		wallB2.append(newPX2)

		#newPB1 = [newPX2[0] + self.pipeWidth*cos(currAngB-pi/2), newPX2[1] + self.pipeWidth*sin(currAngB-pi/2)]
		#wallB1.append(newPB1)

		#self.wallSoup.append(wallB1)
		#segB = self.addSeg(newPB1,newPX2, currAngB)

		if angle <= pi/3:

			newPB1 = [newPX2[0] + self.pipeWidth*cos(currAngB-pi/2), newPX2[1] + self.pipeWidth*sin(currAngB-pi/2)]
			wallB1.append(newPB1)

			self.addLeftWall(segNum, newPB1)
			self.wallSoup.append(wallB1)

			segB = self.addSeg(newPB1,newPX2, currAngB)

		else:
			newPB2 = [p1[0] + self.pipeWidth*cos(currAngB+pi/2), p1[1] + self.pipeWidth*sin(currAngB+pi/2)]
			wallB2.append(newPB2)

			self.addRightWall(segNum, newPB2)
			self.wallSoup.append(wallB2)

			segB = self.addSeg(p1,newPB2, currAngB)



		self.delSeg(segNum)

		#self.activeSegments[segNum] = [wall1[-1], wall2[-1], currAng]

		#if len(wall1) > 1:
		#	self.wallSoup.append(wall1)
		#if len(wall2) > 1:
		#	self.wallSoup.append(wall2)

		return segA, segB

	def addTurn(self, segNum, turnAngle, ):

		p1, p2, currAng, wallL, wallR = self.activeSegments[segNum]

		#print "addTurn(", p1, p2, currAng, turnAngle

		wall1 = [p1]
		wall2 = [p2]

		if turnAngle > 0.0:	
			pivotAngle = (pi - turnAngle)/2.0

			hyp = self.pipeWidth / sin(pivotAngle)
			x = hyp * cos(pivotAngle)
			
			newP1 = [wall1[-1][0] + x*cos(currAng), wall1[-1][1] + x*sin(currAng)]
			self.addLeftWall(segNum, newP1)
			wall1.append(newP1)


			currAng += turnAngle

			newP1 = [wall1[-1][0] + x*cos(currAng), wall1[-1][1] + x*sin(currAng)]
			self.addLeftWall(segNum, newP1)
			wall1.append(newP1)

		else:
			pivotAngle = (pi + turnAngle)/2.0
			
			hyp = self.pipeWidth / sin(pivotAngle)
			x = hyp * cos(pivotAngle)
			
			newP2 = [wall2[-1][0] + x*cos(currAng), wall2[-1][1] + x*sin(currAng)]
			self.addRightWall(segNum, newP2)
			wall2.append(newP2)

			currAng += turnAngle

			newP2 = [wall2[-1][0] + x*cos(currAng), wall2[-1][1] + x*sin(currAng)]
			self.addRightWall(segNum, newP2)
			wall2.append(newP2)

		#self.activeSegments[segNum] = [wall1[-1], wall2[-1], currAng]
		self.updateSeg(segNum, wall1[-1], wall2[-1], currAng)

		if len(wall1) > 1:
			self.wallSoup.append(wall1)
		if len(wall2) > 1:
			self.wallSoup.append(wall2)


	def addStraight(self, segNum, pipeLen):
		
		p1, p2, currAng, wallL, wallR = self.activeSegments[segNum]

		#print "addStraight(", segNum, p1, p2, currAng, pipeLen

		wall1 = [p1]
		wall2 = [p2]


		newP1 = [p1[0] + pipeLen*cos(currAng), p1[1] + pipeLen*sin(currAng)]
		newP2 = [p2[0] + pipeLen*cos(currAng), p2[1] + pipeLen*sin(currAng)]

		self.addLeftWall(segNum, newP1)
		self.addRightWall(segNum, newP2)
		wall1.append(newP1)
		wall2.append(newP2)

		#self.activeSegments[segNum] = [newP1, newP2, currAng]
		self.updateSeg(segNum, newP1, newP2, currAng)

		self.wallSoup.append(wall1)
		self.wallSoup.append(wall2)

	def draw(self):

		print "walls =", repr(self.walls)

		pylab.clf()
		

		for wall in self.walls:
			xP = []
			yP = []
			for p in wall:
				xP.append(p[0])
				yP.append(p[1])	

			pylab.plot(xP,yP, color='k')
		
		pylab.xlim(-2,18)
		pylab.ylim(-10,10)
		pylab.axis('equal')
		pylab.show()


class CrossJunction:
	
	def __init__(self, pipe_len, entr_len, pipe_width, turn_angle):
		global mapCount

		self.pipeWidth = pipe_width
		self.pipeLen = pipe_len
		self.origin = [0.0,0.0]
		self.entrLen = entr_len
		self.turnAngle = turn_angle

		"""
		WLEN = 3.0
		WLEN2 = 5.0
		wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
		wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
		w1 = wall1[2]
		w2 = wall2[0]
		
		wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
		lp = wall3[0]
		rp = wall3[2]
		
		wall6 = [lp, [lp[0] + WLEN*cos(pi/6), lp[1] + WLEN*sin(pi/6)]]
		wall6.append([wall6[1][0] + 0.4*cos(pi/3), wall6[1][1] - 0.4*sin(pi/3)])
		wall6.append([wall6[2][0] - WLEN*cos(pi/6), wall6[2][1] - WLEN*sin(pi/6)])
		wall6.append([wall6[3][0] + WLEN*cos(pi/3), wall6[3][1] - WLEN*sin(pi/3)])
		wall6.append([wall6[4][0] - 0.4*cos(pi/6), wall6[4][1] - 0.4*sin(pi/6)])
		wall6.append(w1)
		wall6.reverse()
		
		wall7 = [[-4.6+6.0,-0.2],[-4.6+6.0,0.2]]
		"""
		
		WLEN = self.pipeLen
		wall1 = [[-14.0, -self.pipeWidth/2.0], [-4.0, -self.pipeWidth/2.0]]
		angle = 0.0

		cornerPoint = wall1[-1]
		
		angle += pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])


		pPoint = wall1[-1]
		
		angle += self.turnAngle
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])
		
		angle += -pi/2
		dist = (pPoint[0]+self.pipeWidth  - wall1[-1][0]) / cos(angle)
		wall1.append([wall1[-1][0] + dist*cos(angle), wall1[-1][1] - dist*sin(angle)])

		angle += -self.turnAngle
		wall1.append([wall1[-1][0], cornerPoint[1]])
		
		cornerPoint = wall1[-1]

		angle += pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		pPoint = wall1[-1]
		
		angle += self.turnAngle
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])

		angle += -pi/2
		dist = -(pPoint[1]+self.pipeWidth  - wall1[-1][1]) / sin(angle)
		wall1.append([wall1[-1][0] + dist*cos(angle), wall1[-1][1] - dist*sin(angle)])

		angle += -self.turnAngle
		wall1.append([cornerPoint[0], wall1[-1][1]])

		cornerPoint = wall1[-1]

		angle += pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		pPoint = wall1[-1]
		
		angle += self.turnAngle
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])

		angle += -pi/2
		dist = (pPoint[0]-self.pipeWidth  - wall1[-1][0]) / cos(angle)
		wall1.append([wall1[-1][0] + dist*cos(angle), wall1[-1][1] - dist*sin(angle)])

		angle += -self.turnAngle
		wall1.append([wall1[-1][0], cornerPoint[1]])

		wall1.append([-14.0, self.pipeWidth/2.0])
		
		"""
		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])
		"""
		
		#angle += -pi/2
		#wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		"""
		angle += pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += pi/2
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		wall1.append([-14.0, self.pipeWidth/2.0])
		"""
		
		self.walls = [wall1]
		
		
		for wall in self.walls:
			for i in range(len(wall)):
				p = copy(wall[i])
				p[0] += 6.0
				wall[i] = p			
		
		
		self.writeToFile(".")
		self.draw()

		mapCount += 1
		
	def writeToFile(self, testDir):

		global mapCount
		
		fileName = "completeFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		fd.write(repr(self.walls))
		fd.close()
		
		return fileName

	def draw(self):
		global mapCount
		
		pylab.clf()
		
		wall = self.walls[0]
				
		for i in range(len(wall)-1):
			n1 = wall[i]
			n2 = wall[i+1]
			
			xP = [n1[0], n2[0]]
			yP = [n1[1], n2[1]]
			pylab.plot(xP,yP, color="b")


class StraightJunction:
	
	def __init__(self, pipe_len, entr_len, pipe_width, turn_angle):
		global mapCount

		self.pipeWidth = pipe_width
		self.pipeLen = pipe_len
		self.origin = [0.0,0.0]
		self.entrLen = entr_len
		self.turnAngle = turn_angle
		
		WLEN = self.pipeLen
		wall1 = [[-14.0, -self.pipeWidth/2.0], [-4.0, -self.pipeWidth/2.0]]
		angle = 0.0
		
		angle += self.turnAngle
		wall1.append([wall1[-1][0] + WLEN*cos(angle), wall1[-1][1] - WLEN*sin(angle)])

		angle += -pi/2
		wall1.append([wall1[-1][0] + self.pipeWidth*cos(angle), wall1[-1][1] - self.pipeWidth*sin(angle)])
		
		angle += -pi/2
		" compute the length to reach x = 0.2 "
		dist = (wall1[-1][1]-self.pipeWidth/2.0) / sin(angle)
		wall1.append([wall1[-1][0] + dist*cos(angle), wall1[-1][1] - dist*sin(angle)])
		
		wall1.append([-14.0, self.pipeWidth/2.0])
		
		self.walls = [wall1]
		
		
		for wall in self.walls:
			for i in range(len(wall)):
				p = copy(wall[i])
				p[0] += 6.0
				wall[i] = p			
		
		
		self.writeToFile(".")
		self.draw()

		mapCount += 1
		
	def writeToFile(self, testDir):

		global mapCount
		
		fileName = "completeFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		fd.write(repr(self.walls))
		fd.close()
		
		return fileName

	def draw(self):
		global mapCount
		
		pylab.clf()
		
		wall = self.walls[0]
				
		for i in range(len(wall)-1):
			n1 = wall[i]
			n2 = wall[i+1]
			
			xP = [n1[0], n2[0]]
			yP = [n1[1], n2[1]]
			pylab.plot(xP,yP, color="b")

class PipeJunctions:
	
	def __init__(self, pipe_len, entr_len, turn_angle):
		global mapCount

		self.rd = Random()

		self.pipeLen = pipe_len
		self.origin = [0.0,0.0]
		self.entrLen = entr_len
		self.turnAngle = turn_angle
		self.nodeCount = 0
		#self.featureRes= 0.01
		self.featureRes= 0.02
		self.pipeMax = 0.45
		self.pipeMin = 0.35
		self.perturbRange = 0.01
		
		self.genCount = 0
		self.branchCount = 0

		self.nodeCount = 0
		self.pathTree = graph.graph()
		self.pathTree.add_node(self.nodeCount, [self.origin[0],self.origin[1],0.0,0.0,self.branchCount])
		self.nodeCount += 1
		
		self.points =[]
		self.nodes = []
		self.completeWalls = []
		
		mapCount += 1
		self.mapID = mapCount
		
		self.generateTree1()
		

		" establish the spline for the main branch "
		intArray = [[],[]]
		print self.nodes
		for n in self.nodes:
			point = self.pathTree.get_node_attributes(n)
			intArray[0].append(point[0])
			intArray[1].append(point[1])					
		self.tck, self.u = scipy.interpolate.splprep(intArray, k=5)
		
		walls = []
		newWalls = self.createWalls(self.nodes)
		walls.append(newWalls)
		
		self.walls = self.adjustWalls(walls)
		#self.mergeWalls()
		self.reorderWallPoints()
		
		self.writeToFile(".")
		#self.draw()
		
	def writeToFile(self, testDir):

		global mapCount

		fileName = "mapFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		
		#print self.walls
		walls = []		
		for wall in self.walls:
			newWall = [] 
			for p in wall[0].points:
				newWall.append(p)
			walls.append(newWall)
			newWall = [] 
			for p in wall[1].points:
				newWall.append(p)
			walls.append(newWall)
		
		fd.write(repr(walls))
		fd.close()
		
		fileName = "completeFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		fd.write(repr(self.completeWalls))
		fd.close()
		
				
		return fileName
	
	def draw(self):
		global mapCount
		
		pylab.clf()
		
		edges = self.pathTree.edges()
		
		for edge in edges:
			n1 = self.pathTree.get_node_attributes(edge[0])
			n2 = self.pathTree.get_node_attributes(edge[1])
			xP = [n1[0], n2[0]]
			yP = [n1[1], n2[1]]
			if n1[4] == 1 or n2[4] == 1:				
				pylab.plot(xP,yP, color="r")
			else:
				pylab.plot(xP,yP, color="b")

		pylab.savefig("mapTop%04u.png" % mapCount)

	def reorderWallPoints(self):
		
		newWall1 = self.walls[0][0].points
		newWall2 = self.walls[0][1].points
		
		#print newWall1
		#print newWall2
		
		newWall1.reverse()
		
		completeWall = newWall2 + newWall1
		
		self.completeWalls.append(completeWall)
		

	def mergeWalls(self):
		" assume that left wall is merged "
	
		" Left Wall of Main Branch "
		R_Main = self.walls[0][0].points
		L_Main = self.walls[0][1].points
		
		B_Right = self.walls[1][0].points
		B_Left = self.walls[1][1].points
		results = []
		
		rightIntersectFound = False
		leftIntersectFound = False
		
		rightP = [0.0,0.0]
		leftP = [0.0,0.0]
	 	mainNodes = [0.0]
		branchNodes = [0,0]

		"  1. find the intersection of walls[2] and walls[3] with walls[1] "
		for i in range(1,len(R_Main)):
			mainEdge = [R_Main[i-1],R_Main[i]]
			for j in range(1, len(B_Right)):
				branchEdge = [B_Right[j-1],B_Right[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, rightP = Intersect(mainEdge, branchEdge)
				 	rightIntersectFound = True
				 	mainNodes = [i-1,i]
	 				branchNodes = [j-1,j]
				 	break
		 
			if rightIntersectFound:
				break
			
	 	results.append([rightIntersectFound, rightP, mainNodes, mainEdge, branchNodes, branchEdge])
	 	mainNodes = [0.0]
		branchNodes = [0,0]

		for i in range(1,len(R_Main)):
			mainEdge = [R_Main[i-1],R_Main[i]]
			for j in range(1, len(B_Left)):
				branchEdge = [B_Left[j-1],B_Left[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, leftP = Intersect(mainEdge, branchEdge)
				 	leftIntersectFound = True
				 	mainNodes = [i-1,i]
	 				branchNodes = [j-1,j]
				 	break
				 
			if leftIntersectFound:
				break		
			

	 	results.append([leftIntersectFound, leftP, mainNodes, mainEdge, branchNodes, branchEdge])

		rightIntersectFound = False
		leftIntersectFound = False
		rightP = [0.0,0.0]
		leftP = [0.0,0.0]
	 	mainNodes = [0.0]
		branchNodes = [0,0]


		for i in range(1,len(L_Main)):
			mainEdge = [L_Main[i-1],L_Main[i]]
			for j in range(1, len(B_Right)):
				branchEdge = [B_Right[j-1],B_Right[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, rightP = Intersect(mainEdge, branchEdge)
				 	mainNodes = [i-1,i]
				 	branchNodes = [j-1,j]
				 	rightIntersectFound = True
				 
			if rightIntersectFound:
				break

	 	results.append([rightIntersectFound, rightP, mainNodes, mainEdge, branchNodes, branchEdge])
	 	mainNodes = [0.0]
		branchNodes = [0,0]

		for i in range(1,len(L_Main)):
			mainEdge = [L_Main[i-1],L_Main[i]]
			for j in range(1, len(B_Left)):
				branchEdge = [B_Left[j-1],B_Left[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, leftP = Intersect(mainEdge, branchEdge)
				 	mainNodes = [i-1,i]
				 	branchNodes = [j-1,j]
				 	leftIntersectFound = True
				 
			if leftIntersectFound:
				break		

	 	results.append([leftIntersectFound, leftP, mainNodes, mainEdge, branchNodes, branchEdge])

	 	#print results[0][0], results[1][0], results[2][0], results[3][0]
		
		#if results[0][0] and results[1][0]:
		if results[0][0] and results[1][0] and results[0][2][0] < results[1][2][0]:
			" Graft branch onto the right side "	
			print "right graft on ", mapCount, results[0][2], results[1][2]

			"main up to intersection"			
			newWall = deepcopy(R_Main[0:results[0][2][0]])
			"intersection point"
			newWall += [copy(results[0][1])]
			"section of branch" 
			newWall += deepcopy(B_Right[results[0][4][1]:])

			newWall2 = deepcopy(R_Main[results[1][2][1]:])
			newWall2 = [copy(results[1][1])] + newWall2
			temp = B_Left[results[1][4][1]:]
			temp.reverse()
			newWall2 = temp + newWall2
			
			completeWall = newWall + newWall2
			completeWall2 = self.walls[0][1].points
			
			completeWall.reverse()
			self.completeWalls.append(completeWall2 + completeWall)
			
			#self.completeWalls.append(completeWall)
			#self.completeWalls.append(completeWall2)
			
		#if results[2][0] and results[3][0]:
		elif results[2][0] and results[3][0] and results[2][2][0] > results[3][2][0]:
			" Graft branch onto the left side "			

			print "left graft on ", mapCount, results[2][2], results[3][2]

		
		
	def adjustWalls(self, walls):
		
		" rotate the points from the origin until the opening goes straight in the positive Y "
		#pos = scipy.interpolate.splev([0.0], self.tck, der=0)
		#dir = scipy.interpolate.splev([0.0], self.tck, der=1)

		originNode = self.pathTree.get_node_attributes(0)
		nextNode = self.pathTree.get_node_attributes(1)

		org = [originNode[0], originNode[1]]
		dir = [nextNode[0]-originNode[0], nextNode[1]-originNode[1]]
		dirMag = sqrt(dir[0]**2 + dir[1]**2)
		dirVec = [dir[0]/dirMag, dir[1]/dirMag]
		
		#org = [pos[0], pos[1]]
		#dirMag = sqrt(dir[0]**2 + dir[1]**2)
		#dirVec = [dir[0]/dirMag, dir[1]/dirMag]

		rotAng = 0.0
		
		# rotate vector to [1,0]
		if dirVec[1] == 0.0:
			if dirVec[0] < 0.0:
				# rotate 180 degrees
				rotAng = pi
			else:
				rotAng = 0

		else:
			b = -1 / (dirVec[1] + dirVec[0]**2 / dirVec[1])
			a = -b * dirVec[0] / dirVec[1]

			rotAng = acos(a)
			if b < 0.0:
				rotAng = -rotAng

		#print "rotating", rotAng, "from a,b=", a, b, "dirVec =", dirVec

		for wall in walls:
	
			for i in range(0,wall[0].getNumPoints()):
				vec1 = [wall[0].points[i][0] - org[0], wall[0].points[i][1] - org[1]]
	
				newVec = [0.0,0.0]
				newVec[0] = vec1[0] * cos(rotAng) - vec1[1] * sin(rotAng)
				newVec[1] = vec1[0] * sin(rotAng) + vec1[1] * cos(rotAng)
	
				wall[0].points[i][0] = org[0] + newVec[0]
				wall[0].points[i][1] = org[1] + newVec[1]
	
			for i in range(0,wall[1].getNumPoints()):
				vec1 = [wall[1].points[i][0] - org[0], wall[1].points[i][1] - org[1]]
	
				newVec = [0.0,0.0]
				newVec[0] = vec1[0] * cos(rotAng) - vec1[1] * sin(rotAng)
				newVec[1] = vec1[0] * sin(rotAng) + vec1[1] * cos(rotAng)
	
				wall[1].points[i][0] = org[0] + newVec[0]
				wall[1].points[i][1] = org[1] + newVec[1]

		# now center the points along the x-axis
		lY = walls[0][1].points[0][1]
		width = (self.pipeMax + self.pipeMin)/2.0
		dLy = -width/2.0 - lY 
		#print lY, width, dLy

		#print "offseting by", dLx

		for wall in walls:
			for i in range(0,len(wall[0].points)):
				wall[0].points[i][1] += dLy

			for i in range(0,len(wall[1].points)):
				wall[1].points[i][1] += dLy

		# now start the entrance at X = 0
		lX = walls[0][0].points[0][0]
		dLx = -lX 
		
		for wall in walls:
			for i in range(0,len(wall[0].points)):
				wall[0].points[i][0] += dLx
	
			for i in range(0,len(wall[1].points)):
				wall[1].points[i][0] += dLx
		
		# now add entrance corridors
		newP1 = deepcopy(walls[0][0].points[0])
		newP2 = deepcopy(walls[0][1].points[0])

		# entrance corridor length
		newP1[0] += -self.entrLen
		newP2[0] += -self.entrLen
		
		walls[0][0].points.insert(0,newP1)
		walls[0][1].points.insert(0,newP2)

		return walls

	def createWalls(self, nodes):

		intArray = [[],[]]
		for n in nodes:
			point = self.pathTree.get_node_attributes(n)
			#intArray.append([point[0],point[1]])
			intArray[0].append(point[0])
			intArray[1].append(point[1])
					
		tck, u = scipy.interpolate.splprep(intArray,k=5)
		# now sample points at regular intervals
		unew = scipy.arange(0, 1.01, self.featureRes)
		pos = scipy.interpolate.splev(unew, tck, der=0)
		dir = scipy.interpolate.splev(unew, tck, der=1)

		rightP = [[],[]]
		leftP = [[],[]]

		for i in range(0,len(dir[0])):
			mag = sqrt(dir[0][i]**2 + dir[1][i]**2)
			dir[0][i] /= mag
			dir[1][i] /= mag
		
			# rotate 90 degrees
			rightP[0].append(dir[0][i]*cos(pi/2) - dir[1][i]*sin(pi/2))
			rightP[1].append(dir[0][i]*sin(pi/2) + dir[1][i]*cos(pi/2))

			# rotate -90 degrees
			leftP[0].append(dir[0][i]*cos(-pi/2) - dir[1][i]*sin(-pi/2))
			leftP[1].append(dir[0][i]*sin(-pi/2) + dir[1][i]*cos(-pi/2))

		points = []

		# the first point has maximum width for the entrance
		newP =[pos[0][0] + (self.pipeMax/2.0)*rightP[0][0], pos[1][0] + (self.pipeMax/2.0)*rightP[1][0]]
		points.append(newP)

		for i in range(1,len(dir[0])):

			#widthR = self.rd.uniform(self.pipeMin,self.pipeMax)
			widthR = self.rd.uniform(0.39,0.41)

			newP =[pos[0][i] + (widthR/2.0)*rightP[0][i], pos[1][i] + (widthR/2.0)*rightP[1][i]]
			points.append(newP)

		wall1 = Wall(self.pipeMin, self.pipeMax, self.perturbRange)
		wall1.setPoints(points)

		points = []

		# the first point has maximum width for the entrance
		newP =[pos[0][0] + (self.pipeMax/2.0)*leftP[0][0], pos[1][0] + (self.pipeMax/2.0)*leftP[1][0]]
		points.append(newP)

		for i in range(1,len(dir[0])):
		
			#widthL = self.rd.uniform(self.pipeMin,self.pipeMax)
			widthL = self.rd.uniform(0.39,0.41)
			newP =[pos[0][i] + (widthL/2.0)*leftP[0][i], pos[1][i] + (widthL/2.0)*leftP[1][i]]
			points.append(newP)

		wall2 = Wall(self.pipeMin, self.pipeMax, self.perturbRange)
		wall2.setPoints(points)


		return wall1, wall2

	def generateTree1(self):
		
		"""
		1. step straight very shortly
		2. follow back to the woods
		3. visible is the light that follows
		4. beckon thither to the stop
		"""
		
		"""
		1. path goes forward 0.2
		2. take a turn for the remainder of pipe_len
		"""
		
		points = []
		nodes = []
		totalLen = 0.0
		
		originNode = self.pathTree.get_node_attributes(0)
		points.append([originNode[0],originNode[1],originNode[2],0.0, self.branchCount])
		nodes.append(0)
		
		for k in range(0,8):
			" first edge goes straight 0.2 "
			edgeLen = 0.2
			relAng = 0.0
			
			" create new point "
			newP = copy(points[-1])
			newP[2] += relAng
			newP[3] = 0.0
			newP[0] += edgeLen * cos(newP[2])
			newP[1] += edgeLen * sin(newP[2])
			
			" set this angle to the previous point "
			points[-1][2] = newP[2]
			points[-1][3] = edgeLen
			
			" add new point to the list "			
			points.append(newP)
	
			" update path distance "
			totalLen += edgeLen
	
			" add to the graph "
			newP = points[-1]
			self.pathTree.add_node(self.nodeCount, newP)
			self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
			nodes.append(self.nodeCount)
			self.nodeCount += 1
		
		
		" now perform the turn "
		edgeLen = 0.2
		relAng = self.turnAngle
		
		" create new point "
		newP = copy(points[-1])
		newP[2] += relAng
		newP[3] = 0.0
		newP[0] += edgeLen * cos(newP[2])
		newP[1] += edgeLen * sin(newP[2])
			
		" set this angle to the previous point "
		points[-1][2] = newP[2]
		points[-1][3] = edgeLen
		
		" add new point to the list "			
		points.append(newP)

		" update path distance "
		totalLen += edgeLen

		" add to the graph "
		newP = points[-1]
		self.pathTree.add_node(self.nodeCount, newP)
		self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
		nodes.append(self.nodeCount)
		self.nodeCount += 1		

		while totalLen < self.pipeLen:
			" now perform the turn "
			edgeLen = 0.2
			relAng = 0.0
			
			" create new point "
			newP = copy(points[-1])
			newP[2] += relAng
			newP[3] = 0.0
			newP[0] += edgeLen * cos(newP[2])
			newP[1] += edgeLen * sin(newP[2])
			
			" set this angle to the previous point "
			points[-1][2] = newP[2]
			points[-1][3] = edgeLen
			
			" add new point to the list "			
			points.append(newP)
	
			" update path distance "
			totalLen += edgeLen
	
			" add to the graph "
			newP = points[-1]
			self.pathTree.add_node(self.nodeCount, newP)
			self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
			nodes.append(self.nodeCount)
			self.nodeCount += 1		


		self.points = points
		self.nodes = nodes

class ForkEnv:
	
	def __init__(self, pipe_len, seg_len_min, seg_len_max, seg_ang_max, pipe_min, pipe_max, feature_len, entr_len, perturb_range):
		global mapCount
		
		self.rd = Random()
		self.pipeLen = pipe_len
		self.segLenMin = seg_len_min
		self.segLenMax = seg_len_max
		self.segAngMax = seg_ang_max
		self.origin = [0.0,0.0]
		self.featureRes = feature_len/pipe_len
		self.entrLen = entr_len
		self.pipeMax = pipe_max
		self.pipeMin = pipe_min
		self.perturbRange = perturb_range
		self.nodeCount = 0
		
		self.genCount = 0
		self.branchCount = 0

		self.nodeCount = 0
		self.pathTree = graph.graph()
		self.pathTree.add_node(self.nodeCount, [self.origin[0],self.origin[1],0.0,0.0,self.branchCount])
		self.nodeCount += 1
		
		self.completeWalls = []
		
		mapCount += 1
		self.mapID = mapCount

		#self.generateTree1()
		self.generateTree2()
		
		# get the pre-ordering from the DFS search
		dfs = self.pathTree.depth_first_search(0)[1]
		
		branches = []
		newBranch = []
		for n in dfs:
			# root
			if len(newBranch) == 0:
				newBranch.append(n)
			else:
				neighbors = self.pathTree.neighbors(n)
				if neighbors.count(newBranch[-1]) > 0:
					newBranch.append(n)
				else:
					# find the stem
					oldBranch = newBranch
					branches.append(oldBranch)
					newBranch = []
					for neigh in neighbors:
						if oldBranch.count(neigh) > 0:
							newBranch.append(neigh)
							break
					newBranch.append(n)
		
		branches.append(newBranch)

		" establish the spline for the main branch "
		intArray = [[],[]]
		for n in branches[0]:
			point = self.pathTree.get_node_attributes(n)
			intArray[0].append(point[0])
			intArray[1].append(point[1])					
		self.tck, self.u = scipy.interpolate.splprep(intArray, k=5)
	
		walls = []
		for branch in branches:
			newWalls = self.createWalls(branch)
			walls.append(newWalls)
		
		self.walls = self.adjustWalls(walls)
		self.mergeWalls()
		
		self.writeToFile(".")
		#self.draw()
		
	def writeToFile(self, testDir):

		global mapCount

		fileName = "mapFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		
		#print self.walls
		walls = []		
		for wall in self.walls:
			newWall = [] 
			for p in wall[0].points:
				newWall.append(p)
			walls.append(newWall)
			newWall = [] 
			for p in wall[1].points:
				newWall.append(p)
			walls.append(newWall)
		
		fd.write(repr(walls))
		fd.close()
		
		fileName = "completeFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		fd.write(repr(self.completeWalls))
		fd.close()
		
				
		return fileName
	
	def draw(self):
		global mapCount
		
		pylab.clf()
		
		edges = self.pathTree.edges()
		
		for edge in edges:
			n1 = self.pathTree.get_node_attributes(edge[0])
			n2 = self.pathTree.get_node_attributes(edge[1])
			xP = [n1[0], n2[0]]
			yP = [n1[1], n2[1]]
			if n1[4] == 1 or n2[4] == 1:				
				pylab.plot(xP,yP, color="r")
			else:
				pylab.plot(xP,yP, color="b")

		pylab.savefig("mapTop%04u.png" % mapCount)


	def mergeWalls(self):
		" assume that left wall is merged "
	
		" Left Wall of Main Branch "
		R_Main = self.walls[0][0].points
		L_Main = self.walls[0][1].points
		
		B_Right = self.walls[1][0].points
		B_Left = self.walls[1][1].points
		results = []
		
		rightIntersectFound = False
		leftIntersectFound = False
		
		rightP = [0.0,0.0]
		leftP = [0.0,0.0]
	 	mainNodes = [0.0]
		branchNodes = [0,0]

		"  1. find the intersection of walls[2] and walls[3] with walls[1] "
		for i in range(1,len(R_Main)):
			mainEdge = [R_Main[i-1],R_Main[i]]
			for j in range(1, len(B_Right)):
				branchEdge = [B_Right[j-1],B_Right[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, rightP = Intersect(mainEdge, branchEdge)
				 	rightIntersectFound = True
				 	mainNodes = [i-1,i]
	 				branchNodes = [j-1,j]
				 	break
		 
			if rightIntersectFound:
				break
			
	 	results.append([rightIntersectFound, rightP, mainNodes, mainEdge, branchNodes, branchEdge])
	 	mainNodes = [0.0]
		branchNodes = [0,0]

		for i in range(1,len(R_Main)):
			mainEdge = [R_Main[i-1],R_Main[i]]
			for j in range(1, len(B_Left)):
				branchEdge = [B_Left[j-1],B_Left[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, leftP = Intersect(mainEdge, branchEdge)
				 	leftIntersectFound = True
				 	mainNodes = [i-1,i]
	 				branchNodes = [j-1,j]
				 	break
				 
			if leftIntersectFound:
				break		
			

	 	results.append([leftIntersectFound, leftP, mainNodes, mainEdge, branchNodes, branchEdge])

		rightIntersectFound = False
		leftIntersectFound = False
		rightP = [0.0,0.0]
		leftP = [0.0,0.0]
	 	mainNodes = [0.0]
		branchNodes = [0,0]


		for i in range(1,len(L_Main)):
			mainEdge = [L_Main[i-1],L_Main[i]]
			for j in range(1, len(B_Right)):
				branchEdge = [B_Right[j-1],B_Right[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, rightP = Intersect(mainEdge, branchEdge)
				 	mainNodes = [i-1,i]
				 	branchNodes = [j-1,j]
				 	rightIntersectFound = True
				 
			if rightIntersectFound:
				break

	 	results.append([rightIntersectFound, rightP, mainNodes, mainEdge, branchNodes, branchEdge])
	 	mainNodes = [0.0]
		branchNodes = [0,0]

		for i in range(1,len(L_Main)):
			mainEdge = [L_Main[i-1],L_Main[i]]
			for j in range(1, len(B_Left)):
				branchEdge = [B_Left[j-1],B_Left[j]]
				if segIntersection(mainEdge, branchEdge):
				 	flag, leftP = Intersect(mainEdge, branchEdge)
				 	mainNodes = [i-1,i]
				 	branchNodes = [j-1,j]
				 	leftIntersectFound = True
				 
			if leftIntersectFound:
				break		

	 	results.append([leftIntersectFound, leftP, mainNodes, mainEdge, branchNodes, branchEdge])

	 	#print results[0][0], results[1][0], results[2][0], results[3][0]
		
		#if results[0][0] and results[1][0]:
		if results[0][0] and results[1][0] and results[0][2][0] < results[1][2][0]:
			" Graft branch onto the right side "	
			print "right graft on ", mapCount, results[0][2], results[1][2]

			"main up to intersection"			
			newWall = deepcopy(R_Main[0:results[0][2][0]])
			"intersection point"
			newWall += [copy(results[0][1])]
			"section of branch" 
			newWall += deepcopy(B_Right[results[0][4][1]:])

			newWall2 = deepcopy(R_Main[results[1][2][1]:])
			newWall2 = [copy(results[1][1])] + newWall2
			temp = B_Left[results[1][4][1]:]
			temp.reverse()
			newWall2 = temp + newWall2
			
			completeWall = newWall + newWall2
			completeWall2 = self.walls[0][1].points
			
			completeWall.reverse()
			self.completeWalls.append(completeWall2 + completeWall)
			
			#self.completeWalls.append(completeWall)
			#self.completeWalls.append(completeWall2)
			
		#if results[2][0] and results[3][0]:
		elif results[2][0] and results[3][0] and results[2][2][0] > results[3][2][0]:
			" Graft branch onto the left side "			

			print "left graft on ", mapCount, results[2][2], results[3][2]

		
		
	def adjustWalls(self, walls):
		
		# rotate the points from the origin until the opening goes straight in the positive Y
		pos = scipy.interpolate.splev([0.0], self.tck, der=0)
		dir = scipy.interpolate.splev([0.0], self.tck, der=1)

		org = [pos[0], pos[1]]
		dirMag = sqrt(dir[0]**2 + dir[1]**2)
		dirVec = [dir[0]/dirMag, dir[1]/dirMag]

		rotAng = 0.0
		
		# rotate vector to [1,0]
		if dirVec[1] == 0.0:
			if dirVec[0] < 0.0:
				# rotate 180 degrees
				rotAng = pi
			else:
				rotAng = 0

		else:
			b = -1 / (dirVec[1] + dirVec[0]**2 / dirVec[1])
			a = -b * dirVec[0] / dirVec[1]

			rotAng = acos(a)
			if b < 0.0:
				rotAng = -rotAng

		#print "rotating", rotAng, "from a,b=", a, b, "dirVec =", dirVec

		for wall in walls:
	
			for i in range(0,wall[0].getNumPoints()):
				vec1 = [wall[0].points[i][0] - org[0], wall[0].points[i][1] - org[1]]
	
				newVec = [0.0,0.0]
				newVec[0] = vec1[0] * cos(rotAng) - vec1[1] * sin(rotAng)
				newVec[1] = vec1[0] * sin(rotAng) + vec1[1] * cos(rotAng)
	
				wall[0].points[i][0] = org[0] + newVec[0]
				wall[0].points[i][1] = org[1] + newVec[1]
	
			for i in range(0,wall[1].getNumPoints()):
				vec1 = [wall[1].points[i][0] - org[0], wall[1].points[i][1] - org[1]]
	
				newVec = [0.0,0.0]
				newVec[0] = vec1[0] * cos(rotAng) - vec1[1] * sin(rotAng)
				newVec[1] = vec1[0] * sin(rotAng) + vec1[1] * cos(rotAng)
	
				wall[1].points[i][0] = org[0] + newVec[0]
				wall[1].points[i][1] = org[1] + newVec[1]

		# now center the points along the x-axis
		lY = walls[0][1].points[0][1]
		width = (self.pipeMax + self.pipeMin)/2.0
		dLy = -width/2.0 - lY 
		#print lY, width, dLy

		#print "offseting by", dLx

		for wall in walls:
			for i in range(0,len(wall[0].points)):
				wall[0].points[i][1] += dLy

			for i in range(0,len(wall[1].points)):
				wall[1].points[i][1] += dLy

		# now start the entrance at X = 0
		lX = walls[0][0].points[0][0]
		dLx = -lX 
		
		for wall in walls:
			for i in range(0,len(wall[0].points)):
				wall[0].points[i][0] += dLx
	
			for i in range(0,len(wall[1].points)):
				wall[1].points[i][0] += dLx
		
		# now add entrance corridors
		newP1 = deepcopy(walls[0][0].points[0])
		newP2 = deepcopy(walls[0][1].points[0])

		# entrance corridor length
		newP1[0] += -self.entrLen
		newP2[0] += -self.entrLen
		
		walls[0][0].points.insert(0,newP1)
		walls[0][1].points.insert(0,newP2)


		return walls

	def createWalls(self, nodes):

		intArray = [[],[]]
		for n in nodes:
			point = self.pathTree.get_node_attributes(n)
			#intArray.append([point[0],point[1]])
			intArray[0].append(point[0])
			intArray[1].append(point[1])
					
		tck, u = scipy.interpolate.splprep(intArray,k=5)
		# now sample points at regular intervals
		unew = scipy.arange(0, 1.01, self.featureRes)
		pos = scipy.interpolate.splev(unew, tck, der=0)
		dir = scipy.interpolate.splev(unew, tck, der=1)

		rightP = [[],[]]
		leftP = [[],[]]

		for i in range(0,len(dir[0])):
			mag = sqrt(dir[0][i]**2 + dir[1][i]**2)
			dir[0][i] /= mag
			dir[1][i] /= mag
		
			# rotate 90 degrees
			rightP[0].append(dir[0][i]*cos(pi/2) - dir[1][i]*sin(pi/2))
			rightP[1].append(dir[0][i]*sin(pi/2) + dir[1][i]*cos(pi/2))

			# rotate -90 degrees
			leftP[0].append(dir[0][i]*cos(-pi/2) - dir[1][i]*sin(-pi/2))
			leftP[1].append(dir[0][i]*sin(-pi/2) + dir[1][i]*cos(-pi/2))

		points = []

		# the first point has maximum width for the entrance
		newP =[pos[0][0] + (self.pipeMax/2.0)*rightP[0][0], pos[1][0] + (self.pipeMax/2.0)*rightP[1][0]]
		points.append(newP)

		for i in range(1,len(dir[0])):

			widthR = self.rd.uniform(self.pipeMin,self.pipeMax)

			newP =[pos[0][i] + (widthR/2.0)*rightP[0][i], pos[1][i] + (widthR/2.0)*rightP[1][i]]
			points.append(newP)

		wall1 = Wall(self.pipeMin, self.pipeMax, self.perturbRange)
		wall1.setPoints(points)

		points = []

		# the first point has maximum width for the entrance
		newP =[pos[0][0] + (self.pipeMax/2.0)*leftP[0][0], pos[1][0] + (self.pipeMax/2.0)*leftP[1][0]]
		points.append(newP)

		for i in range(1,len(dir[0])):
		
			widthL = self.rd.uniform(self.pipeMin,self.pipeMax)
			newP =[pos[0][i] + (widthL/2.0)*leftP[0][i], pos[1][i] + (widthL/2.0)*leftP[1][i]]
			points.append(newP)

		wall2 = Wall(self.pipeMin, self.pipeMax, self.perturbRange)
		wall2.setPoints(points)


		return wall1, wall2

	def buildBranch(self, stem, pipeLen = 6.0):
		
		
		points = []
		totalLen = 0.0

		originNode = self.pathTree.get_node_attributes(stem)
		#points.append([originNode[0],originNode[1],originNode[2],0.0, self.branchCount])
		points.append([originNode[0],originNode[1],originNode[2]+pi/2.0,0.0, self.branchCount])
		#points.append([originNode[0],originNode[1],pi/2,0.0])

		while totalLen < pipeLen:

			# if this is the first edge, start in the y-direction
			if len(points) == 1:
				edgeLen = self.segLenMax
				relAng = 0.0
			else:
				edgeLen = self.rd.uniform(self.segLenMin,self.segLenMax)
				relAng = self.rd.uniform(-self.segAngMax,self.segAngMax)	

			newP = copy(points[-1])
			newP[2] += relAng
			newP[3] = 0.0

			# set this angle to the previous point
			points[-1][2] = newP[2]
			points[-1][3] = edgeLen

			newP[0] += edgeLen * cos(newP[2])
			newP[1] += edgeLen * sin(newP[2])
			
			points.append(newP)
			totalLen += edgeLen
		
			newEdge = [points[-2],points[-1]]
			isIntersect = False

			# check if there's a self-intersection
			for i in range(1,len(points)-2):
				oldEdge = [points[i-1], points[i]]
			
				# if an intersection, reset and start over
				if segIntersection(newEdge, oldEdge):
					isIntersect = True
					totalLen = 0.0
					points = []
					originNode = self.pathTree.get_node_attributes(stem)
					points.append([originNode[0],originNode[1],0.0,0.0, self.branchCount])
					break				

		self.branchCount += 1
		
		return points
		
	def generateTree1(self):
		
		while True:
						
			mainBranch = self.buildBranch(0,pipeLen = 3.0)
			# add to graph
			
			for i in range(1,len(mainBranch)):
				newP = mainBranch[i]
				#print newP
				self.pathTree.add_node(self.nodeCount, newP)
				self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
				self.nodeCount += 1
	
			stemNode = self.nodeCount/2
			secondBranch = self.buildBranch(stemNode,pipeLen=3.0)
			
			if self.validateTree():
				break
			else:
				self.nodeCount = 0
				self.branchCount = 0
				self.pathTree = graph.graph()
				self.pathTree.add_node(self.nodeCount, [self.origin[0],self.origin[1],0.0,0.0,self.branchCount])
				self.nodeCount += 1
				self.genCount += 1
				
	def generateTree2(self):
		
		while True:
						
			mainBranch = self.buildBranch(0,pipeLen = 3.0)
			# add to graph
			
			for i in range(1,len(mainBranch)):
				newP = mainBranch[i]
				#print newP
				self.pathTree.add_node(self.nodeCount, newP)
				self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
				self.nodeCount += 1
	
			stemNode = self.nodeCount/2
			secondBranch = self.buildBranch(stemNode,pipeLen=3.0)
			
			# add to graph and verify
			newP = secondBranch[1]
			self.pathTree.add_node(self.nodeCount, newP)
			self.pathTree.add_edge(stemNode, self.nodeCount)
			self.nodeCount += 1
	
			for i in range(2,len(secondBranch)):
				newP = secondBranch[i]
				self.pathTree.add_node(self.nodeCount, newP)
				self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
				self.nodeCount += 1

			if self.validateTree():
				break
			else:
				self.nodeCount = 0
				self.branchCount = 0
				self.pathTree = graph.graph()
				self.pathTree.add_node(self.nodeCount, [self.origin[0],self.origin[1],0.0,0.0,self.branchCount])
				self.nodeCount += 1
				self.genCount += 1

	def validateTree(self):

		isIntersect = False

		for edge1 in self.pathTree.edges():
			n1 = self.pathTree.get_node_attributes(edge1[0])
			n2 = self.pathTree.get_node_attributes(edge1[1])
			newEdge = [n1,n2]
			
			# check if there's a self-intersection			
			for edge2 in self.pathTree.edges():
				n3 = self.pathTree.get_node_attributes(edge2[0])
				n4 = self.pathTree.get_node_attributes(edge2[1])
				oldEdge = [n3, n4]
			
				# if an intersection, reset and start over
				if edge1[0] != edge2[0] and edge1[0] != edge2[1]:
					if edge1[1] != edge2[0] and edge1[1] != edge2[1]:
						if segIntersection(newEdge, oldEdge):
							isIntersect = True
							break
		
		return not isIntersect

	def generateTree(self):

		"""
		1. add the branches to the stem
		2. randomly select next direction within a range
		3. add a random number of nodes
		
		"""

		self.nodeCount = 0
		self.pathTree = graph.graph()

		self.pathTree.add_node(self.nodeCount, [self.origin[0],self.origin[1],0.0,0.0])
		self.nodeCount += 1

		while self.nodeCount < 20:

			# if this is the first edge, start in the y-direction
			if self.nodeCount == 1:
				edgeLen = self.segLenMax
				relAng = 0.0
			else:
				edgeLen = self.rd.uniform(self.segLenMin,self.segLenMax)
				relAng = self.rd.uniform(-self.segAngMax,self.segAngMax)

			prevNode = self.pathTree.get_node_attributes(self.nodeCount-1)

			newP = copy(prevNode)
			newP[2] += relAng
			newP[3] = 0.0

			# set this angle to the previous point
			prevNode[2] = newP[2]
			prevNode[3] = edgeLen

			newP[0] += edgeLen * cos(newP[2])
			newP[1] += edgeLen * sin(newP[2])
			
			newEdge = [prevNode, newP]
			isIntersect = False

			# check if there's a self-intersection			
			for edge in self.pathTree.edges():
				n1 = self.pathTree.get_node_attributes(edge[0])
				n2 = self.pathTree.get_node_attributes(edge[1])
				oldEdge = [n1, n2]
			
				# if an intersection, reset and start over
				if edge[0] != self.nodeCount-1 and edge[1] != self.nodeCount-1 and segIntersection(newEdge, oldEdge):
					isIntersect = True
					break
			

			if not isIntersect:
				self.pathTree.add_node(self.nodeCount, newP)
				self.pathTree.add_edge(self.nodeCount-1, self.nodeCount)
				self.nodeCount += 1
				
			else:
				# reset the graph
				self.nodeCount = 0
				self.pathTree = graph.graph()
				self.pathTree.add_node(self.nodeCount, [self.origin[0],self.origin[1],0.0,0.0])
				self.nodeCount += 1
		
		# add a branch ( start at node 10
		edgeLen = self.segLenMax
		branchPoint = self.pathTree.get_node_attributes(10)
		newP = copy(branchPoint)
		newP[2] += pi/4
		newP[3] = 0.0

		# set this angle to the previous point
		prevNode[2] = newP[2]
		prevNode[3] = edgeLen

		newP[0] += edgeLen * cos(newP[2])
		newP[1] += edgeLen * sin(newP[2])
		
		newEdge = [prevNode, newP]
		isIntersect = False
		
		

class Corridor:
	
	def __init__(self, pipe_len, seg_len_min, seg_len_max, seg_ang_max, pipe_min, pipe_max, feature_len, entr_len, perturb_range):

		global mapCount
	
		self.numPoints = 0	
		self.totalLen = 0.0
		self.rd = Random()
		self.points = []
		self.lengths = []
		self.pipeLen = pipe_len
		self.entrWidth = 0.4
		self.segLenMin = seg_len_min
		self.segLenMax = seg_len_max
		self.segAngMax = seg_ang_max
		self.origin = [-14.0,0.0]


		#self.featureRes = feature_len
		self.featureRes = feature_len/pipe_len
		#print "featureRes = ", self.featureRes
		self.entrLen = entr_len
		self.pipeMax = pipe_max
		self.pipeMin = pipe_min
		self.perturbRange = perturb_range

		self.genCount = 0

		mapCount += 1

		self.mapID = mapCount

		while True:

			self.generatePath()

			# unzip the points
			intArray = [[],[]]
			for p in self.points:
				intArray[0].append(p[0])
				intArray[1].append(p[1])
		
			self.tck, self.u = scipy.interpolate.splprep(intArray)

			unew = scipy.arange(0, 1.01, self.featureRes)
			#self.pos = scipy.interpolate.splev([0.0], self.tck, der=0)
			#self.dir = scipy.interpolate.splev([0.0], self.tck, der=1)
			self.pos = scipy.interpolate.splev(unew, self.tck, der=0)
			self.dir = scipy.interpolate.splev(unew, self.tck, der=1)
			


			wallCount = 0
			self.wall1, self.wall2 = self.createWalls()

			print "checkA"

			while not self.checkConstraints():
				wallCount += 1

				print "checkB"
				if wallCount > 10:
					#print "wallCount =", wallCount
					break

				self.wall1, self.wall2 = self.createWalls()

			# perturb each point on the wall
			"""
			for i in range(1,self.wall1.getNumPoints()):
				self.wall1.perturbPoint(i, self.wall2)

			# perturb each point on the wall
			for i in range(1,self.wall2.getNumPoints()):
				self.wall2.perturbPoint(i, self.wall1)
			"""

			# rotate the points from the origin until the opening goes straight in the positive Y
			#pos = scipy.interpolate.splev([0.0], self.tck, der=0)
			#dir = scipy.interpolate.splev([0.0], self.tck, der=1)

			org = [self.pos[0][0], self.pos[1][0]]
			dirMag = sqrt(self.dir[0][0]**2 + self.dir[1][0]**2)
			dirVec = [self.dir[0][0]/dirMag, self.dir[1][0]/dirMag]

			rotAng = 0.0
			
			# rotate vector to [1,0]
			if dirVec[1] == 0.0:
				if dirVec[0] < 0.0:
					# rotate 180 degrees
					rotAng = pi
				else:
					rotAng = 0

			else:
				#b = 1 / (dirVec[0] + dirVec[1]**2 / dirVec[0])
				#a = b * dirVec[1] / dirVec[0]
				b = -1 / (dirVec[1] + dirVec[0]**2 / dirVec[1])
				a = -b * dirVec[0] / dirVec[1]

			rotAng = acos(a)
			if b < 0.0:
				rotAng = -rotAng

			#print "rotating", rotAng, "from a,b=", a, b, "dirVec =", dirVec

			for i in range(0,self.wall1.getNumPoints()):
				vec1 = [self.wall1.points[i][0] - org[0], self.wall1.points[i][1] - org[1]]

				newVec = [0.0,0.0]
				newVec[0] = vec1[0] * cos(rotAng) - vec1[1] * sin(rotAng)
				newVec[1] = vec1[0] * sin(rotAng) + vec1[1] * cos(rotAng)

				self.wall1.points[i][0] = org[0] + newVec[0]
				self.wall1.points[i][1] = org[1] + newVec[1]

			for i in range(0,self.wall2.getNumPoints()):
				vec1 = [self.wall2.points[i][0] - org[0], self.wall2.points[i][1] - org[1]]

				newVec = [0.0,0.0]
				newVec[0] = vec1[0] * cos(rotAng) - vec1[1] * sin(rotAng)
				newVec[1] = vec1[0] * sin(rotAng) + vec1[1] * cos(rotAng)

				self.wall2.points[i][0] = org[0] + newVec[0]
				self.wall2.points[i][1] = org[1] + newVec[1]

			# now center the points along the x-axis
			lY = self.wall2.points[0][1]
			width = (self.pipeMax + self.pipeMin)/2.0
			dLy = -width/2.0 - lY 
			print lY, width, dLy

			#print "offseting by", dLx

			for i in range(0,len(self.wall1.points)):
				self.wall1.points[i][1] += dLy

			for i in range(0,len(self.wall2.points)):
				self.wall2.points[i][1] += dLy

			# now start the entrance at X = 0
			lX = self.wall1.points[0][0]
			dLx = -lX 
			for i in range(0,len(self.wall1.points)):
				self.wall1.points[i][0] += dLx

			for i in range(0,len(self.wall2.points)):
				self.wall2.points[i][0] += dLx
		
			# now add entrance corridors
			newP1 = deepcopy(self.wall1.points[0])
			newP2 = deepcopy(self.wall2.points[0])

			# entrance corridor length
			newP1[0] += -self.entrLen
			newP2[0] += -self.entrLen
			
			self.wall1.points.insert(0,newP1)
			self.wall2.points.insert(0,newP2)

			" adjust the Y so the entrance is centered on the X axis "
			yAdj = self.pipeMin/2.0 - newP1[1]

			for i in range(len(self.wall1.points)):
				self.wall1.points[i][1] += yAdj
			for i in range(len(self.wall2.points)):
				self.wall2.points[i][1] += yAdj


			print "entrance points:", self.wall1.points[0], self.wall2.points[0]

			print "checkC"
			
			# now verify the walls again
			if self.checkConstraints():
				print "generated after", self.genCount

				self.wall1.points.reverse()
				return
			

			print "checkD"

			# start over
			self.genCount += 1

			if self.genCount >= 100:
				print "Generation failed after 100 iterations!"
				raise

	def writeToFile(self, testDir):

		global mapCount

		fileName = "mapFile%04d.txt" % mapCount

		# record the corridor points
		fd = open(testDir + "/" + fileName, 'w')
		
		walls = []
		wall1 = []
		wall2 = []

		for p in self.wall1.points:
			np = copy(p)
			np[0] += 2.0
			wall1.append(np)

		for p in self.wall2.points:
			np = copy(p)
			np[0] += 2.0
			wall2.append(np)
			
		
		#walls = [wall1, wall2]
		walls = [wall2 + wall1]
		fd.write(repr(walls))
		fd.close()
		
		"""
		for p in self.wall1.points:
			fd.write(str(p[0]) + " ")
			fd.write(str(p[1]) + " ")
		fd.write("\n")

		for p in self.wall2.points:
			fd.write(str(p[0]) + " ")
			fd.write(str(p[1]) + " ")
		fd.write("\n")

		fd.close()
		"""
		
		return fileName

	def getTermPoint(self):

		pos = self.wall1.points[-1]
		#pos = scipy.interpolate.splev([1.0], self.tck, der=0)
		return pos


	def checkConstraints(self):
		# check pipe width constraints

		# for every point on wall1, check if they're too close to wall2
		for tP in self.wall1.points:
			dist, p = self.wall2.findClosest(tP)
			if dist > self.pipeMax:
				return False

			if dist < self.pipeMin:
				return False

		# for every point on wall2, check if they're too close to wall1
		for tP in self.wall2.points:
			dist, p = self.wall1.findClosest(tP)
			if dist > self.pipeMax:
				return False

			if dist < self.pipeMin:
				return False

		return True

	def createWalls(self):

		# now sample points at regular intervals

		unew = scipy.arange(0, 1.01, self.featureRes)
		#pos = scipy.interpolate.splev(unew, self.tck, der=0)
		#dir = scipy.interpolate.splev(unew, self.tck, der=1)

		rightP = [[],[]]
		leftP = [[],[]]

		for i in range(0,len(self.dir[0])):
			mag = sqrt(self.dir[0][i]**2 + self.dir[1][i]**2)
			self.dir[0][i] /= mag
			self.dir[1][i] /= mag
		
			# rotate 90 degrees
			rightP[0].append(self.dir[0][i]*cos(pi/2) - self.dir[1][i]*sin(pi/2))
			rightP[1].append(self.dir[0][i]*sin(pi/2) + self.dir[1][i]*cos(pi/2))

			# rotate -90 degrees
			leftP[0].append(self.dir[0][i]*cos(-pi/2) - self.dir[1][i]*sin(-pi/2))
			leftP[1].append(self.dir[0][i]*sin(-pi/2) + self.dir[1][i]*cos(-pi/2))


		print "pos:", self.pos[0][0], self.pos[1][0]

		points = []

		# the first point has maximum width for the entrance
		#newP =[self.pos[0][0] + (self.pipeMax/2.0)*rightP[0][0], self.pos[1][0] + (self.pipeMax/2.0)*rightP[1][0]]
		#newP =[self.pos[0][0] + (self.pipeMin/2.0)*rightP[0][0], self.pos[1][0] + (self.pipeMin/2.0)*rightP[1][0]]
		newP =[self.pos[0][0] + (self.entrWidth/2.0)*rightP[0][0], self.pos[1][0] + (self.entrWidth/2.0)*rightP[1][0]]
		points.append(newP)

		for i in range(1,len(self.dir[0])):

			widthR = self.rd.uniform(self.pipeMin,self.pipeMax)

			newP =[self.pos[0][i] + (widthR/2.0)*rightP[0][i], self.pos[1][i] + (widthR/2.0)*rightP[1][i]]
			points.append(newP)

		wall1 = Wall(self.pipeMin, self.pipeMax, self.perturbRange)
		wall1.setPoints(points)

		points = []

		# the first point has maximum width for the entrance
		#newP =[self.pos[0][0] + (self.pipeMax/2.0)*leftP[0][0], self.pos[1][0] + (self.pipeMax/2.0)*leftP[1][0]]
		newP =[self.pos[0][0] + (self.pipeMin/2.0)*leftP[0][0], self.pos[1][0] + (self.pipeMin/2.0)*leftP[1][0]]
		points.append(newP)

		for i in range(1,len(self.dir[0])):
		
			widthL = self.rd.uniform(self.pipeMin,self.pipeMax)
			newP =[self.pos[0][i] + (widthL/2.0)*leftP[0][i], self.pos[1][i] + (widthL/2.0)*leftP[1][i]]
			points.append(newP)

		wall2 = Wall(self.pipeMin, self.pipeMax, self.perturbRange)
		wall2.setPoints(points)

		return wall1, wall2


	def drawCorridor(self, clr = '0.5'):
		self.wall1.printWall()
		self.wall2.printWall()

	def generatePath(self):

		#print "Generating pipe path!"

		self.points = []
		self.totalLen = 0.0
		self.lengths = []

		self.points.append([self.origin[0],self.origin[1],0.0,0.0])

		while self.totalLen < self.pipeLen:

			# if this is the first edge, start in the y-direction
			if len(self.points) == 1:
				edgeLen = self.segLenMax
				#relAng = pi/2
				relAng = 0.0
			else:
				edgeLen = self.rd.uniform(self.segLenMin,self.segLenMax)
				relAng = self.rd.uniform(-self.segAngMax,self.segAngMax)

			newP = copy(self.points[-1])
			newP[2] += relAng
			newP[3] = 0.0

			# set this angle to the previous point
			self.points[-1][2] = newP[2]
			self.points[-1][3] = edgeLen

			newP[0] += edgeLen * cos(newP[2])
			newP[1] += edgeLen * sin(newP[2])
			
			self.points.append(newP)
			self.totalLen += edgeLen
		
			self.lengths.append(self.totalLen)

			newEdge = [self.points[-2],self.points[-1]]

			# check if there's a self-intersection
			for i in range(1,len(self.points)-2):
				oldEdge = [self.points[i-1], self.points[i]]
			
				# if an intersection, reset and start over
				if segIntersection(newEdge, oldEdge):
					self.points = []
					self.points.append([self.origin[0],self.origin[1],0.0,0.0])
					self.totalLen = 0.0
					self.lengths = []
					#print "wall discarded, restarting!"
					break

class Wall:
	
	def __init__(self, pipe_min, pipe_max, perturb_range):
		self.numPoints = 0	
		self.totalLen = 0.0
		self.rd = Random()
		self.points = []
		self.lengths = []
		self.origin = [0.0,0.0]
		self.pipeMax = pipe_max
		self.pipeMin = pipe_min
		self.perturbRange = perturb_range

	def perturbPoint(self, i, refWall):

		# perturb the i'th point and make sure it follows the pipe width constraints to refWall
		original = self.points[i]

		#print "perturbing", i

		while True:

			#self.points[i][0] = original[0] + self.rd.uniform(-PERTURB_RANGE,PERTURB_RANGE)
			#self.points[i][1] = original[1] + self.rd.uniform(-PERTURB_RANGE,PERTURB_RANGE)
			self.points[i][0] = original[0] + self.rd.gauss(0.0, self.perturbRange/2.0)
			self.points[i][1] = original[1] + self.rd.gauss(0.0, self.perturbRange/2.0)

			minDist, resP = refWall.findClosest(self.points[i])	

			if i > 0:
				edge1 = [deepcopy(self.points[i-1]), deepcopy(self.points[i])]

			if i < len(self.points)-1:
				edge2 = [deepcopy(self.points[i]), deepcopy(self.points[i+1])]

			minDist1 = 1e1000
			minDist2 = 1e1000

			for j in range(0,refWall.getNumPoints()):

				if i > 0:
					newDist1, resP = closestSegPoint(edge1,refWall.points[j])
					if newDist1 < minDist1:
						minDist1 = newDist1

				if i < len(self.points)-1:
					newDist2, resP = closestSegPoint(edge2,refWall.points[j])
					if newDist2 < minDist2:
						minDist2 = newDist2

			if i > 0:
				cond1 = False
			else:
				cond1 = True

			if i < len(self.points)-1:
				cond2 = False
			else:
				cond2 = True

			cond3 = False

			if i > 0 and minDist1 < self.pipeMax and minDist1 > self.pipeMin:
				cond1 = True

			if i < len(self.points)-1 and minDist2 < self.pipeMax and minDist2 > self.pipeMin:
				cond2 = True

			if minDist < self.pipeMax and minDist > self.pipeMin: 
				cond3 = True
	
			#print "minDist =", minDist1, minDist2, minDist
				
			if cond1 and cond2 and cond3:
				#print "returned"
				return

	def getNumPoints(self):
		return len(self.points)

	def getLength(self):
		return self.totalLen

	def setPoints(self, points):

		for i in range(1,len(points)):

			newP = [points[i-1][0], points[i-1][1], 0.0, 0.0]
			nextP = points[i]
			edgeLen = sqrt((newP[0]-nextP[0])**2 + (newP[1]-nextP[1])**2)
			dirVec = [nextP[0]-newP[0], nextP[1]-newP[1]]
			magDir = sqrt(dirVec[0]**2 + dirVec[1]**2)
			dirVec[0] /= magDir
			dirVec[1] /= magDir

			ang = acos(dirVec[0])
			if dirVec[1] < 0.0:
				ang = -ang

			self.lengths.append(edgeLen)
			self.points.append([newP[0], newP[1], edgeLen, ang])

	def checkIntersect(self, seg):

		for i in range(1, len(self.points)):
			wSeg = [self.points[i-1], self.points[i]]
			
			result = segIntersection(seg, wSeg)	

			if result:
				return True

		return False

	def findClosest(self, tP):
		
		# minimum distance so far
		minDist = 1e100
		minPoint = [0.0,0.0]

		# check the distance of each segment to the point
		for i in range(1, len(self.points)):
			edge = [deepcopy(self.points[i]), deepcopy(self.points[i-1])]
			dist, resP = closestSegPoint(edge,tP)

			if dist < minDist:
				minPoint = resP
				minDist = dist

		return minDist, minPoint


	def getLenPos(self, dist):
		
		if dist > self.totalLen or dist < 0.0:
			print "value", dist, "outside the length of wall"
			raise

		for i in range(0,len(self.lengths)):
			if dist < self.lengths[i]:
				# this is the correct starting point 
				startPos = self.points[i]
				xPos = startPos[0]
				yPos = startPos[1]

				print "starting position =", startPos

				lastDist = self.lengths[i-1]
				diffDist = dist - lastDist

				print "lastDist =", lastDist
				print "diffDist =", diffDist

				xPos += diffDist * cos(startPos[2])
				yPos += diffDist * sin(startPos[2])

				print "result =", xPos,yPos

				return [xPos, yPos]
				
	def printWall(self):	

		# plot the points
		xP = []
		yP = []

		for p in self.points:
			xP.append(p[0])
			yP.append(p[1])

		#pylab.scatter(xP,yP)
		pylab.plot(xP,yP)


class MapFile:
	def __init__(self, fileName):
		self.walls = []

		self.readFile(fileName)

	def readFile(self, fileName):

		# read in file
		lineReader = csv.reader(file(fileName), delimiter=' ')

		for row in lineReader:
			newWall = []
			for i in range(0,len(row)/2):
				newWall.append([float(row[2*i]),float(row[2*i+1])])
	
			self.walls.append(newWall)

	def getTermPoint(self):
		pos = self.walls[0][-1]
		return pos

	def printMap(self):
		for wall in self.walls:
			xP = []
			yP = []

			for p in wall:
				xP.append(p[0])
				yP.append(p[1])

			pylab.plot(xP, yP, color='g')


