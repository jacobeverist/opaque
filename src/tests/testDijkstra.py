from numpy import *
from math import *
import scipy.linalg

def normalizeAngle(angle):

	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle 

def doTransform(T1, T2, E1, E2):

	x1 = T1[0,0]	
	y1 = T1[1,0]
	p1 = T1[2,0]

	x2 = T2[0,0]	
	y2 = T2[1,0]
	p2 = T2[2,0]

	newT = matrix([[x1 + x2*cos(p1) - y2*sin(p1)],
		[y1 + x2*sin(p1) + y2*cos(p1)],
		[normalizeAngle(p1+p2)]],dtype=float)

	J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]],dtype=float)
	J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]],dtype=float)
						
	newE = J1 * E1 * J1.T + J2 * E2 * J2.T
	
	return newT, newE



Th1 = matrix([[1.0],[ 0.0],[ pi/2.0]])
Ch1 = matrix([[ 0.0001 ,0.  ,   0.    ],[ 0.   ,  0.0001, 0.    ],[ 0.   ,  0.  ,   0.02  ]])

Tp1 = matrix([[1.0],[ 0.0],[ pi/2.0]])
Ep1 = matrix([[ 0.1 , 0.  , 0.  ],[ 0. ,  0.01 ,0.  ],[ 0. ,  0. ,  0.02]])

Th2 = matrix([[1.0],[ 0.0],[ pi/2.0]])
Ch2 = matrix([[ 0.0001 ,0.    , 0.    ],[ 0.  ,   0.0001 ,0.    ],[ 0.  ,   0. ,    0.02  ]])

Tp2 = matrix([[1.0],[ 0.0],[ pi/2.0]])
Ep2 = matrix([[ 0.1  ,0.  , 0.  ],[ 0.  , 0.01 ,0.  ],[ 0. ,  0. ,  0.02]])


covE = Ch1
result1 = Th1
result2, cov2 = doTransform(result1, Tp1, covE, Ep1)
result3, cov3 = doTransform(result2, Th2, cov2, Ch2)
result4, cov4 = doTransform(result3, Tp2, cov3, Ep2)


invMat = scipy.linalg.inv(cov4)

print "invMat:", invMat

err = sqrt(transpose(result4) * invMat * result4)
#results.append([err, i, j])



print result4, cov4
print err

