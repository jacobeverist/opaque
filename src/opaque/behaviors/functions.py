
from math import *
from copy import *
from numpy import array, dot, transpose

# this function converts the angle to its equivalent # in the range [-pi,pi]
def normalizeAngle(angle):

	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle 


# determines signed area of 3 points (used for solving Point in Polygon problem)
def Area2(Ax,Ay,Bx,By,Cx,Cy):
	return (Bx - Ax) * (Cy - Ay) - (Cx - Ax)*(By - Ay)

# determines if point is left of line segment
def LeftOnEdge(edge, point):
	return (Area2(edge[0][0],edge[0][1],edge[1][0],edge[1][1],point[0],point[1]) >= 0)

# determines if point C is left of line segment AB
def LeftOn(Ax,Ay,Bx,By,Cx,Cy):
	return (Area2(Ax,Ay,Bx,By,Cx,Cy) >= 0)

# determine if point is located within a bounding box specified by vertices RectVert
def IsContained(RectVert, Point):
	for i in range(4):
		if not (LeftOn(RectVert[i%4][0],RectVert[i%4][1],
					   RectVert[(i+1)%4][0],RectVert[(i+1)%4][1], Point[0], Point[1])):
			return False
	return True

def point_inside_polygon(x,y,poly):

	n = len(poly)
	inside = False
	
	xinters = -1e30
	
	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y
	
	return inside


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



class Pose:
	
	def __init__(self, pose = [0.0,0.0,0.0]):
		self.setEstPose(pose)

	def setEstPose(self, newPose):

		self.estPose = copy(newPose)
		self.dist = sqrt(self.estPose[0]**2 + self.estPose[1]**2)

		if self.dist != 0:
			self.vecAng = acos(self.estPose[0]/self.dist)

		if asin(self.estPose[1]/self.dist) < 0:
			self.vecAng = -self.vecAng
		else:
			self.vecAng = 0.0

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
							