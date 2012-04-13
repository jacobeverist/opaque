from math import *
from copy import *
from func import IsContained

# this function converts the angle to its equivalent # in the range [-pi,pi]
def normalizeAngle(angle):

	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle 

def diffAngle(angle1, angle2):
	
	midPhi = angle1
	
	tempAngle1 = 0.0
	tempAngle2 = normalizeAngle(angle2 - midPhi)

	return tempAngle1-tempAngle2

def findClosestPoint(a_trans, b):

	minDist = 1e100
	minPoint = None
	min_i = 0

	for i in range(len(a_trans)):
		p = a_trans[i]

		dist = sqrt((p[0]-b[0])**2 + (p[1]-b[1])**2)
		if dist < minDist:
			minPoint = copy(p)
			minDist = dist
			min_i = i

	if minPoint != None:
		return minPoint, min_i, minDist
	else:
		raise
	
def makePointsUniform(points, max_spacing = 0.04):
	
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
		
		if dist > max_spacing:
			" cut into pieces max_spacing length or less "
			numCount = int(floor(dist / max_spacing))
			
			for j in range(1, numCount+1):
				newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
				new_points.append(newP)
	
	return new_points		


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
#def IsContained(RectVert, Point):
#	for i in range(4):
#		if not (LeftOn(RectVert[i%4][0],RectVert[i%4][1],
#					   RectVert[(i+1)%4][0],RectVert[(i+1)%4][1], Point[0], Point[1])):
#			return False
#	return True

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
