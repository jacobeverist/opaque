cdef extern from "math.h":
	double cos(double)
	double sin(double)


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
	cdef double q11, q12, g21, q22 = 0.0
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
	cdef double q11, q12, g21, q22 = 0.0
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
	