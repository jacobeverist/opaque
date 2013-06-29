from libc.stdlib cimport malloc, free

cdef extern from "runNelmin.h":
	void doTest(double *matchPairs, int numPairs, double *initGuess, double uHigh, double uLow, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum)	
	void doCost(double *matchPairs, int numPairs, double *initGuess, double uHigh, double uLow, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset)	

def ICPcost(matchPairs, int numPairs, initGuess, uHigh, uLow, poses_1, poses_2, int numPoses):

	#cdef double c_matchPairs[12*numPairs]
	cdef double *c_matchPairs = <double *>malloc(12*numPairs * sizeof(double))
	cdef double *c_initGuess = <double *>malloc(3 * sizeof(double))
	cdef double *c_poses_1 = <double *>malloc(3*numPoses * sizeof(double))
	cdef double *c_poses_2 = <double *>malloc(3*numPoses * sizeof(double))
	
	for i in range(12*numPairs):
		c_matchPairs[i] = matchPairs[i]
	
	for i in range(3*numPoses):
		c_poses_1[i] = poses_1[i]
		c_poses_2[i] = poses_2[i]
	
	c_initGuess[0] = initGuess[0]
	c_initGuess[1] = initGuess[1]
	c_initGuess[2] = initGuess[2]
	
	cdef double resultParam[2]
	cdef double resultSum[1]
	cdef double resultOffset[3]
	
	
	doCost(c_matchPairs, numPairs, c_initGuess, uHigh, uLow, c_poses_1, c_poses_2, numPoses, resultParam, resultSum, resultOffset)	
	
	free(c_matchPairs)
	free(c_initGuess)
	free(c_poses_1)
	free(c_poses_2)
	return resultSum[0], [resultParam[0],resultParam[1]], [resultOffset[0], resultOffset[1], resultOffset[2]]
	
	

def ICPmin(matchPairs, int numPairs, initGuess, uHigh, uLow, poses_1, poses_2, int numPoses):

	#cdef double c_matchPairs[12*numPairs]
	cdef double *c_matchPairs = <double *>malloc(12*numPairs * sizeof(double))
	cdef double *c_initGuess = <double *>malloc(3 * sizeof(double))
	cdef double *c_poses_1 = <double *>malloc(3*numPoses * sizeof(double))
	cdef double *c_poses_2 = <double *>malloc(3*numPoses * sizeof(double))

	for i in range(12*numPairs):
		c_matchPairs[i] = matchPairs[i]

	for i in range(3*numPoses):
		c_poses_1[i] = poses_1[i]
		c_poses_2[i] = poses_2[i]
	
	c_initGuess[0] = initGuess[0]
	c_initGuess[1] = initGuess[1]
	c_initGuess[2] = initGuess[2]

	cdef double resultParam[2]
	cdef double resultSum[1]
	

	doTest(c_matchPairs, numPairs, c_initGuess, uHigh, uLow, c_poses_1, c_poses_2, numPoses, resultParam, resultSum)	

	free(c_matchPairs)
	free(c_initGuess)
	free(c_poses_1)
	free(c_poses_2)
	return [resultParam[0], resultParam[1]], resultSum[0]


