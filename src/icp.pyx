from libc.stdlib cimport malloc, free

cdef extern from "math.h":
	double cos(double)
	double sin(double)
	double sqrt(double)

cdef findClosestPointInA(int n, double *a_trans, double *b, double *c_minDist):
    
	cdef double minDist = 1e100
	cdef double minPoint[2]
	cdef int min_i = 0
	cdef double dist
	cdef double p[2]
    
	for i in range(n):
		p[0] = a_trans[2*i]
		p[1] = a_trans[2*i+1]
    
		dist = sqrt((p[0]-b[0])**2 + (p[1]-b[1])**2)
		if dist < minDist:
			minPoint[0] = p[0]
			minPoint[1] = p[1]
			minDist = dist
			min_i = i
	
	c_minDist[0] = minDist
    
	return min_i

cdef computeEnclosingCircle(int n, double *a_data, double *centerA):
            
	cdef double maxA = 0.0
	cdef double a_p1[2]
	cdef double a_p2[2]
	cdef double p1[2]
	cdef double p2[2]
	cdef double dist
	cdef double radiusA
	#cdef double centerA[2]

	for i in range(n):
		p1[0] = a_data[2*i]
		p1[1] = a_data[2*i+1]
    
		for j in range(i+1, n):
			p2[0] = a_data[2*j]
			p2[1] = a_data[2*j+1]
			dist = sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]))
    
			if dist > maxA:
				maxA = dist
				a_p1[0] = p1[0]
				a_p1[1] = p1[1]
				a_p2[0] = p2[0]
				a_p2[1] = p2[1]

	centerA[0] = (a_p1[0] + a_p2[0])/2.0
	centerA[1] = (a_p1[1] + a_p2[1])/2.0
	radiusA = maxA/2.0

	#return radiusA, [centerA[0],centerA[1]]
	return radiusA

cdef isInCircle(double *p, double radius, double *center):

	cdef double dist
        
	dist = sqrt((p[0] - center[0])*(p[0]-center[0]) + (p[1] - center[1])*(p[1] - center[1]))
	if dist < radius:
		return 1
    
	return 0
 
def matchPairs(points, points_offset, globalPoints, minMatchDist):

	cdef int numPoints1 = len(points_offset)
	cdef double *c_pointsOffset = <double *>malloc(numPoints1*2*sizeof(double))

	cdef double center2[2]
	cdef double radius2

	cdef double p_1[2]

	cdef int i
	for i in range(numPoints1):
		c_pointsOffset[2*i] = points_offset[i][0]
		c_pointsOffset[2*i+1] = points_offset[i][1]
		#c_pointsOffset[3*i+2] = points_offset[i][2]

    
	" get the circles and radii "
	radius2 = computeEnclosingCircle(numPoints1, c_pointsOffset, center2) 


	cdef double c_minDist
	cdef int numPoints2 = len(globalPoints)
	cdef int i_2

	match_pairs = []

	for i in range(numPoints2):
		#p_1 = [globalPoints[i][0],globalPoints[i][1]]
		p_1[0] = globalPoints[i][0]
		p_1[1] = globalPoints[i][1]
    
		if isInCircle(p_1, radius2, center2) == 1:

			" for every point of B, find it's closest neighbor in transformed A "
			#p_2, i_2, minDist = findClosestPointInA(points_offset, [p_1[0],p_1[1]])
			i_2 = findClosestPointInA(numPoints1, c_pointsOffset, p_1, &c_minDist)

			if c_minDist <= minMatchDist:
	    
				" add to the list of match pairs less than 1.0 distance apart "
				" keep A points and covariances untransformed "
				C2 = points[i_2][2]
				C1 = globalPoints[i][2]

				" we store the untransformed point, but the transformed covariance of the A point "
				match_pairs.append([points[i_2],globalPoints[i],C2,C1])
			
	free(c_pointsOffset)
	#print "match_pairs:", match_pairs
	return match_pairs



def shapeCostC(offset, match_pairs):

	cdef double sum = 0.0
	cdef double xd, yd, theta = 0.0
	cdef double ax = 0.0
	cdef double ay = 0.0
	cdef double bx = 0.0
	cdef double by1 = 0.0
	cdef double c11, c12, c21, c22 = 0.0
	cdef double b11, b12, b21, b22 = 0.0

	cdef double offset_c[3]
	cdef double a_c[2]
	cdef double b_c[2]
	cdef double Ca_c[4]
	cdef double Cb_c[4]

	for pair in match_pairs:

		a = pair[0]
		b = pair[1]
		Ca = pair[2]
		Cb = pair[3]

		xd = offset[0]
		yd = offset[1]
		theta = offset[2]
		
		ax = a[0]
		ay = a[1]		
		bx = b[0]
		by1 = b[1]

		c11 = Ca[0][0]
		c12 = Ca[0][1]
		c21 = Ca[1][0]
		c22 = Ca[1][1]
				
		b11 = Cb[0][0]
		b12 = Cb[0][1]
		b21 = Cb[1][0]
		b22 = Cb[1][1]	

		# [xd, yd, theta]
		offset_c[0] = xd
		offset_c[1] = yd
		offset_c[2] = theta

		# [ax,ay]
		a_c[0] = ax
		a_c[1] = ay
		
		# [bx,by1]
		b_c[0] = bx
		b_c[1] = by1
		
		# [c11,c12,c21,c22]
		Ca_c[0] = c11
		Ca_c[1] = c12
		Ca_c[2] = c21
		Ca_c[3] = c22

		# [b11,b12,b21,b22]
		Cb_c[0] = b11
		Cb_c[1] = b12
		Cb_c[2] = b21
		Cb_c[3] = b22
	
		sum += computeMatchErrorC(offset_c, a_c, b_c, Ca_c, Cb_c)
		
	return sum


def computeMatchErrorP(offset, a, b, Ca, Cb):

	cdef double xd, yd, theta = 0.0
	cdef double ax = 0.0
	cdef double ay = 0.0
	cdef double bx = 0.0
	cdef double by1 = 0.0
	cdef double tx, ty, dx, dy = 0.0
	cdef double r11, r12, r21, r22 = 0.0
	cdef double c11, c12, c21, c22 = 0.0
	cdef double b11, b12, b21, b22 = 0.0
	cdef double res11, res12, res21, res22, resDet = 0.0
	cdef double q11, q12, q21, q22 = 0.0
	cdef errVal = 0.0

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]
	
	ax = a[0]
	ay = a[1]
	
	bx = b[0]
	by1 = b[1]

	tx = ax*cos(theta) - ay*sin(theta) + xd
	ty = ax*sin(theta) + ay*cos(theta) + yd
	dx = bx - tx
	dy = by1 - ty

	r11 = cos(theta)
	r12 = -sin(theta)
	r21 = sin(theta)
	r22 = r11

	c11 = Ca[0]
	c12 = Ca[1]
	c21 = Ca[2]
	c22 = Ca[3]
	
	b11 = Cb[0]
	b12 = Cb[1]
	b21 = Cb[2]
	b22 = Cb[3]
	
	res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
	res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
	res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
	res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)
	
	resDet = res22*res11 - res12*res21
	
	q11 = res22/resDet
	q12 = -res12/resDet
	q21 = -res21/resDet
	q22 = res11/resDet

	errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22)

	return errVal

cdef computeMatchErrorC(double *offset, double *a, double *b, double *Ca, double *Cb):

	cdef double xd, yd, theta = 0.0
	cdef double ax = 0.0
	cdef double ay = 0.0
	cdef double bx = 0.0
	cdef double by1 = 0.0
	cdef double tx, ty, dx, dy = 0.0
	cdef double r11, r12, r21, r22 = 0.0
	cdef double c11, c12, c21, c22 = 0.0
	cdef double b11, b12, b21, b22 = 0.0
	cdef double res11, res12, res21, res22, resDet = 0.0
	cdef double q11, q12, q21, q22 = 0.0
	cdef errVal = 0.0

	xd = offset[0]
	yd = offset[1]
	theta = offset[2]
	
	ax = a[0]
	ay = a[1]
	
	bx = b[0]
	by1 = b[1]

	tx = ax*cos(theta) - ay*sin(theta) + xd
	ty = ax*sin(theta) + ay*cos(theta) + yd
	dx = bx - tx
	dy = by1 - ty

	r11 = cos(theta)
	r12 = -sin(theta)
	r21 = sin(theta)
	r22 = r11

	c11 = Ca[0]
	c12 = Ca[1]
	c21 = Ca[2]
	c22 = Ca[3]
	
	b11 = Cb[0]
	b12 = Cb[1]
	b21 = Cb[2]
	b22 = Cb[3]
	
	res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
	res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
	res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
	res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)
	
	resDet = res22*res11 - res12*res21
	
	q11 = res22/resDet
	q12 = -res12/resDet
	q21 = -res21/resDet
	q22 = res11/resDet

	errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22)

	return errVal
	
