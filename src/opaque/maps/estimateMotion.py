from SplineFit import SplineFit
import computeCovar
from math import *

def makeGuess(centerCurve1, centerCurve2, stepDist):

	" align the root vectors "
	vec1 = centerCurve1.getUVector(0.5)
	vec2 = centerCurve2.getUVector(0.5)
	
	angle1 = acos(vec1[0])
	if asin(vec1[1]) < 0:
		angle1 = -angle1

	angle2 = acos(vec2[0])
	if asin(vec2[1]) < 0:
		angle2 = -angle2
		
	angle = angle2 - angle1
	
	" center the root node positions "
	
	xOff = stepDist * cos(angle1)
	yOff = stepDist * sin(angle1)
	
	return [xOff, yOff, -angle]


if __name__ == '__main__':

	dirName = "../../testData/poseTest"
	numPoses = 16
	
	localPostures = []
	centerCurves = []

	estPoses = []
	gndPoses = []
	
	for i in range(numPoses):
			
		f = open(dirName + "/posture%04u.txt" % i, 'r')
		localPostures.append(eval(f.read().rstrip()))
		f.close()
		
		" compute a fitted curve of a center point "
		centerCurves.append(SplineFit(localPostures[i], smooth = 0.5, kp = 2))
		
		
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
	naives = []

	for i in range(numPoses-1):
		gndOffset, estOffset = computeCovar.getMotionConstraint(gndPoses[i], gndPoses[i+1], estPoses[i], estPoses[i+1])
		gnds.append(gndOffset)
		ests.append(estOffset)
		
		naiveOffset = makeGuess(centerCurves[i], centerCurves[i+1], 0.14)

		naives.append(naiveOffset)
		
	#print gnds[0], naives[0]

	print computeCovar.computeCovar(gnds, ests)
	print computeCovar.computeCovar(gnds, gnds)
	print computeCovar.computeCovar(gnds, naives)

	for i in range(numPoses-1):
		print gnds[i]	
	for i in range(numPoses-1):
		print naives[i]