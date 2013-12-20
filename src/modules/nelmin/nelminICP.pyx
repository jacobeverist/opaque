from libc.stdlib cimport malloc, free

cdef extern from "runNelmin.h":

	void doTest(double *matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum)	
	void doCost(double *matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset)	
	void doICPIndAngle(double *h_matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *pose1, double *poses_2, int numPoses, double *resultParam, double *resultSum)

def ICPcost(matchPairs, int numPairs, initGuess, uHigh, uLow, poses_1, poses_2, int numPoses):

	#cdef double c_matchPairs[12*numPairs]
	cdef double *c_matchPairs = <double *>malloc(12*numPairs * sizeof(double))
	cdef double *c_initGuess = <double *>malloc(3 * sizeof(double))
	cdef double *c_poses_2 = <double *>malloc(3*numPoses * sizeof(double))
	
	for i in range(12*numPairs):
		c_matchPairs[i] = matchPairs[i]
	
	for i in range(3*numPoses):
		c_poses_2[i] = poses_2[i]
	
	c_initGuess[0] = initGuess[0]
	c_initGuess[1] = initGuess[1]
	c_initGuess[2] = initGuess[2]
	
	cdef double resultParam[2]
	cdef double resultSum[1]
	cdef double resultOffset[3]
	cdef double u1 = initGuess[0]
	cdef double c_pose1[3]
	cdef double angNom = initGuess[2]
	cdef double angLim = 1e100

	c_pose1[0] = 0.0
	c_pose1[1] = 0.0
	c_pose1[2] = 0.0


	if u1 >= 1.0:
		c_pose1[0] = poses_1[3*(numPoses-1)+0]
		c_pose1[1] = poses_1[3*(numPoses-1)+1]
		c_pose1[2] = poses_1[3*(numPoses-1)+2]
	elif u1 < 0.0:
		c_pose1[0] = poses_1[0+0]
		c_pose1[1] = poses_1[0+1]
		c_pose1[2] = poses_1[0+2]
	else:
		c_pose1[0] = poses_1[3*int(u1*(numPoses-1))+0]
		c_pose1[1] = poses_1[3*int(u1*(numPoses-1))+1]
		c_pose1[2] = poses_1[3*int(u1*(numPoses-1))+2]


	
	doCost(c_matchPairs, numPairs, c_initGuess, angNom, angLim, uHigh, uLow, c_pose1, c_poses_2, numPoses, resultParam, resultSum, resultOffset)	
	
	free(c_matchPairs)
	free(c_initGuess)
	free(c_poses_2)
	return resultSum[0], [resultParam[0],resultParam[1]], [resultOffset[0], resultOffset[1], resultOffset[2]]
	
	

def ICPmin(matchPairs, int numPairs, initGuess, uHigh, uLow, poses_1, poses_2, int numPoses):

	#cdef double c_matchPairs[12*numPairs]
	cdef double *c_matchPairs = <double *>malloc(12*numPairs * sizeof(double))
	cdef double *c_initGuess = <double *>malloc(3 * sizeof(double))
	cdef double *c_poses_2 = <double *>malloc(3*numPoses * sizeof(double))

	for i in range(12*numPairs):
		c_matchPairs[i] = matchPairs[i]

	for i in range(3*numPoses):
		c_poses_2[i] = poses_2[i]
	
	c_initGuess[0] = initGuess[0]
	c_initGuess[1] = initGuess[1]
	c_initGuess[2] = initGuess[2]

	cdef double resultParam[2]
	cdef double resultSum[1]
	cdef double u1 = initGuess[0]
	cdef double c_pose1[3]
	cdef double angNom = initGuess[2]
	cdef double angLim = 1e100
	c_pose1[0] = 0.0
	c_pose1[1] = 0.0
	c_pose1[2] = 0.0


	if u1 >= 1.0:
		c_pose1[0] = poses_1[3*(numPoses-1)+0]
		c_pose1[1] = poses_1[3*(numPoses-1)+1]
		c_pose1[2] = poses_1[3*(numPoses-1)+2]
	elif u1 < 0.0:
		c_pose1[0] = poses_1[0+0]
		c_pose1[1] = poses_1[0+1]
		c_pose1[2] = poses_1[0+2]
	else:
		c_pose1[0] = poses_1[3*int(u1*(numPoses-1))+0]
		c_pose1[1] = poses_1[3*int(u1*(numPoses-1))+1]
		c_pose1[2] = poses_1[3*int(u1*(numPoses-1))+2]


	doTest(c_matchPairs, numPairs, c_initGuess, angNom, angLim, uHigh, uLow, c_pose1, c_poses_2, numPoses, resultParam, resultSum)	

	free(c_matchPairs)
	free(c_initGuess)
	free(c_poses_2)
	return [resultParam[0], resultParam[1]], resultSum[0]



def ICPbyPose(matchPairs, int numPairs, initGuess, angNom, angLim, uHigh, uLow, pose1, poses_2, int numPoses):

	#cdef double c_matchPairs[12*numPairs]
	cdef double *c_matchPairs = <double *>malloc(12*numPairs * sizeof(double))
	cdef double *c_initGuess = <double *>malloc(3 * sizeof(double))
	cdef double *c_poses_2 = <double *>malloc(3*numPoses * sizeof(double))

	for i in range(12*numPairs):
		c_matchPairs[i] = matchPairs[i]

	for i in range(3*numPoses):
		c_poses_2[i] = poses_2[i]
	
	c_initGuess[0] = initGuess[0]
	c_initGuess[1] = initGuess[1]
	c_initGuess[2] = initGuess[2]

	cdef double resultParam[2]
	cdef double resultSum[1]
	cdef double u1 = initGuess[0]
	cdef double c_pose1[3]
	c_pose1[0] = pose1[0]
	c_pose1[1] = pose1[1]
	c_pose1[2] = pose1[2]


	doICPIndAngle(c_matchPairs, numPairs, c_initGuess, angNom, angLim, uHigh, uLow, c_pose1, c_poses_2, numPoses, resultParam, resultSum)	

	free(c_matchPairs)
	free(c_initGuess)
	free(c_poses_2)

	return [resultParam[0], resultParam[1]], resultSum[0]


