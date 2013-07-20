
# determines signed area of 3 points (used for solving Point in Polygon problem)
cdef Area2(Ax,Ay,Bx,By,Cx,Cy):
	return (Bx - Ax) * (Cy - Ay) - (Cx - Ax)*(By - Ay)

# determines if point C is left of line segment AB
cdef LeftOn(Ax,Ay,Bx,By,Cx,Cy):
	return (Area2(Ax,Ay,Bx,By,Cx,Cy) >= 0)

# determine if point is located within a bounding box specified by vertices RectVert
def IsContained(RectVert, Point):

	cdef double Ax = 0.0
	cdef double Ay = 0.0
	cdef double Bx = 0.0
	cdef double By = 0.0
	cdef double Cx = 0.0
	cdef double Cy = 0.0
	Cx = Point[0]
	Cy = Point[1]

	for i in range(4):
		Ax = RectVert[i%4][0]
		Ay = RectVert[i%4][1]
		Bx = RectVert[(i+1)%4][0]
		By = RectVert[(i+1)%4][1]

		result = LeftOn(Ax, Ay, Bx, By, Cx, Cy)

		#result = LeftOn(RectVert[i%4][0],RectVert[i%4][1], RectVert[(i+1)%4][0],RectVert[(i+1)%4][1], Point[0], Point[1])
		if not result:
			return False

	return True


