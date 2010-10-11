from pulp import *
from math import *

import pylab

contactCount = 0

def convertFriction(u, norm, point):
	global contactCount
		
	contactCount += 1

	angle = atan(u)

	# pyramid vectors
	vec1 = [cos(angle)*norm[0] - sin(angle)*norm[1], sin(angle)*norm[0] + cos(angle)*norm[1], 0.0]
	vec2 = [cos(angle)*norm[0] + sin(angle)*norm[1], -sin(angle)*norm[0] + cos(angle)*norm[1], 0.0]

	# skew symmetric matrix of point vector
	pMat = [[0, -point[2], point[1]],
					[point[2], 0, -point[0]],
					[-point[1], point[0], 0]]

	# compute the cross product
	prod1 = [pMat[0][0]*vec1[0] + pMat[0][1]*vec1[1] + pMat[0][2]*vec1[2],
						pMat[1][0]*vec1[0] + pMat[1][1]*vec1[1] + pMat[1][2]*vec1[2],
						pMat[2][0]*vec1[0] + pMat[2][1]*vec1[1] + pMat[2][2]*vec1[2]]

	prod2 = [pMat[0][0]*vec2[0] + pMat[0][1]*vec2[1] + pMat[0][2]*vec2[2],
						pMat[1][0]*vec2[0] + pMat[1][1]*vec2[1] + pMat[1][2]*vec2[2],
						pMat[2][0]*vec2[0] + pMat[2][1]*vec2[1] + pMat[2][2]*vec2[2]]

	val1 = LpVariable("f_" + str(contactCount) + "_1", 0, 1e10)
	val2 = LpVariable("f_" + str(contactCount) + "_2", 0, 1e10)
	X = [val1, val2]

	F = [ [vec1[0],  vec2[0]],
				[vec1[1],  vec2[1]],		
				[vec1[2],  vec2[2]],		
				[prod1[0], prod2[0]],
				[prod1[1], prod2[1]],
				[prod1[2], prod2[2]]]

	return F, X

def processContacts(contacts):

	C1 = [[], [], [], [], [], []]
	X = []

	for contact in contacts:
		# friction coefficient, normal vector, contact point
		F1, X1 = convertFriction(contact[0], contact[1], contact[2])

		for j in range(6):
			C1[j] = C1[j] + F1[j]

		X = X + X1

	return C1, X

def checkStable(contacts, forceVec, forcePoint):

	global contactCount 

	# Create the 'prob' variable to contain the problem data
	prob = LpProblem("Stability problem",LpMinimize)

	xc = forcePoint[0]
	yc = forcePoint[1]

	fx = forceVec[0]
	fy = forceVec[1]

	# matrices built from the inputs
	Y = [xc, yc, 0]
	C2 = [[0,0,0],
				[0,0,0],
				[0,0,0],
				[0,0,-fy],
				[0,0,fx],
				[fy,-fx,0]]

	d = [-fx,-fy,0,0,0,0]

	C1, X = processContacts(contacts)

	# minimums
	#for i in range(4):
	#	prob += X[i] >= 0.0

	# maximums
	#for i in range(4):
	#	prob += X[i] <= 10.0

	for i in range(6):
		prob += lpDot(C1[i],X) + lpDot(C2[i],Y) == d[i]

	status = prob.solve()

	contactCount = 0

	#help(prob)

	if status == 1:
		#print "STABLE"
		solutions = []
		for val in X:
			solutions.append(value(val))
		return True, solutions
	else:
		solutions = []
		for val in X:
			solutions.append(value(val))

		#print "NOT STABLE"
		return False, solutions

	#print prob.infeasibilityGap()


# friction coefficient, normal vector, contact point
contacts = []
contacts.append([0.5, [0.0,1.0,0.0], [0.0,0.0,0.0]])
contacts.append([0.5, [0.0,-1.0,0.0], [1.0,1.0,0.0]])
contacts.append([0.5, [0.0,1.0,0.0], [2.0,0.0,0.0]])
#contacts.append([1.0, [0.0,-1.0,0.0], [1.1,1.0,0.0]])

# inputs to the test

#forceOrigins = [[0.0,0.0],[0.5,0.5],[1.0,1.0]]
forceOrigins = [[0.5,0.5]]
forceMags = [1.0, 10.0, 100.0]
angles = [0.0, pi/4, pi/2, 3*pi/4, pi, -3*pi/4, -pi/2, -pi/4]

results = []

for p in forceOrigins:
	for ang in angles:
		for mag in forceMags:
			fx = mag * cos(ang)
			fy = mag * sin(ang)

			result, X = checkStable(contacts, [fx, fy], [p[0], p[1]])

			results.append([result, ang, mag, p[0], p[1], X])

success = 0
fails = 0


for result in results:
	if result[0]:
		success += 1
	else:
		fails += 1
		#print result

print success, "successes"
print fails, "failures"

#print checkStable(contacts, [fx, fy], [xc, yc])

#C1 = [[-sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, sqrt(2)/2],
#			[sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, -sqrt(2)/2],
#			[0,0,0,0],
#			[0,0,0,0],
#			[0,0,0,0],
#			[sqrt(2)/2, sqrt(2)/2, 0, -sqrt(2)]]

