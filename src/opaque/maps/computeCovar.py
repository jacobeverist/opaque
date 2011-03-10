
import numpy
from math import *
from functions import normalizeAngle

def computeCovar(gnds, ests):
	
	N = len(ests)
	
	for i in range(N):
		x = ests[i][0]
		y = ests[i][1]
		p = ests[i][2]

		xGnd = gnds[i][0]
		yGnd = gnds[i][1]
		pGnd = gnds[i][2]

	#print "compute mean of gnd offset"
	avgX = 0.0
	avgY = 0.0
	avgPx = 0.0
	avgPy = 0.0
	for i in range(N):
		avgX += gnds[i][0]
		avgY += gnds[i][1]
		
		P = gnds[i][2]
		avgPx += cos(P)
		avgPy += sin(P)
		
	avgX /= N
	avgY /= N
	avgPx /= N
	avgPy /= N
	avgP = atan2(avgPy, avgPx)
	
	#print "average x,y:", avgX, avgY
	#print "average dist:", avgDist
	#print "average angle:", avgP

	stdX = 0.0
	stdY = 0.0
	stdP = 0.0

	Exx = 0.0
	Exy = 0.0
	Exp = 0.0
	Eyp = 0.0
	Eyy = 0.0
	Epp = 0.0
	
	for i in range(N):
		
		x = ests[i][0]
		y = ests[i][1]
		p = ests[i][2]
		
		Exx += (x - avgX) * (x - avgX)
		Exy += (x - avgX) * (y - avgY)
		Exp += (x - avgX) * (p - avgP)
		Eyp += (y - avgY) * (p - avgP)
		Eyy += (y - avgY) * (y - avgY)
		Epp += (p - avgP) * (p - avgP)
		
	#print "pre:", Exx, Eyy, Epp, Exy, Exp, Eyp
	Exx /= N
	Exy /= N
	Exp /= N
	Eyp /= N
	Eyy /= N
	Epp /= N
	
	E = numpy.matrix([[Exx, Exy, Exp],
				[Exy, Eyy, Eyp],
				[Exp, Eyp, Epp]])
	#print "covar:", Exx, Eyy, Epp, Exy, Exp, Eyp	

	return E


def getMotionConstraint(gndpose1, gndpose2, estpose1, estpose2):	

	xA = gndpose1[0]
	yA = gndpose1[1]
	pA = gndpose1[2]
	
	xB = gndpose2[0]
	yB = gndpose2[1]
	pB = gndpose2[2]

	xT = cos(pA)*(xB-xA) + sin(pA)*(yB-yA)
	yT = -sin(pA)*(xB-xA) + cos(pA)*(yB-yA)
	pT = normalizeAngle(pB - pA)
	
	gndOffset = [xT, yT, pT]
	gndDist = sqrt(xT*xT + yT*yT)			

	xA = estpose1[0]
	yA = estpose1[1]
	pA = estpose1[2]
	
	xB = estpose2[0]
	yB = estpose2[1]
	pB = estpose2[2]

	xT = cos(pA)*(xB-xA) + sin(pA)*(yB-yA)
	yT = -sin(pA)*(xB-xA) + cos(pA)*(yB-yA)
	pT = normalizeAngle(pB - pA)
	
	estOffset = [xT, yT, pT]
	estDist = sqrt(xT*xT + yT*yT)

	return gndOffset, estOffset

if __name__ == '__main__':
	
	#dirName = "../../testData/correctedNX"
	#numPoses = 30
	dirName = "../../testData/poseTest"
	numPoses = 16
	
	estPoses = []
	gndPoses = []
	
	# import est and gnd poses
	for i in range(numPoses):
	
		f = open(dirName + "/estpose%04u.txt" % i, 'r')
		estPose = eval(f.read().rstrip())
		f.close()
	
		estPoses.append(estPose)
	
		f = open(dirName + "/gndpose%04u.txt" % i, 'r')
		gndPose = eval(f.read().rstrip())
		f.close()
	
		gndPoses.append(gndPose)
	
	gnds = []
	ests = []

	for i in range(numPoses-1):
		gndOffset, estOffset = getMotionConstraint(gndPoses[i], gndPoses[i+1], estPoses[i], estPoses[i+1])
		gnds.append(gndOffset)
		ests.append(estOffset)
		
	print computeCovar(gnds, ests)
