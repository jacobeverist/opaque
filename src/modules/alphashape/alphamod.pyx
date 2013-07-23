cdef extern from "alpha2.h":
	int runAlg(int, double, double*, int*, double*)

#cdef extern from "skeletonizeCV.h":
	#int runAlg(int, int, char*, char*)
	#int runAlg(int, int, int, int, double*, double*, char*, char*)
	#int runAlg(int, int, int, double*, double*, char*, char*)
	#int runAlg(int, int, int, int, int, double*, double*, char*)

cdef extern from "stdlib.h":
	void free(void* ptr)
	void* malloc(size_t size)
	void* realloc(void* ptr, size_t size)

def doAlpha(aRadius, points):

	cdef int nPoints = len(points)
	cdef double cRadius = <double> aRadius
	cdef double *cPoints = <double *>malloc(2 * nPoints * sizeof(double))

	cdef double *cVerts = <double *>malloc(2 * nPoints * sizeof(double))
	cdef int *cNumVert = <int *>malloc(sizeof(int))

	cNumVert[0] = 0

	for i in range(nPoints):
		cPoints[2*i] = <double> points[i][0]
		cPoints[2*i+1] = <double> points[i][1]

		cVerts[2*i] = 0
		cVerts[2*i+1] = 0

	runAlg(nPoints,cRadius,cPoints,cNumVert,cVerts)


	alphaPoly = []
	numVert = cNumVert[0]

	for i in range(numVert):
		alphaPoly.append([cVerts[2*i],cVerts[2*i+1]])

	free(cPoints)
	free(cVerts)
	free(cNumVert)

	return alphaPoly

"""
def doToro(v_list, e_list):

	cdef int nPoints = len(v_list)
	cdef double *cVertPoses = <double *>malloc(3 * nPoints * sizeof(double))
	cdef int *cVertIDs = <int *>malloc(nPoints * sizeof(int))

	cdef int nEdges = len(e_list)
	cdef int *cVertPairs = <int *>malloc(2 * nEdges * sizeof(int))
	cdef double *cOffsets = <double *>malloc(3 * nEdges * sizeof(double))
	cdef double *cCovars = <double *>malloc(6 * nEdges * sizeof(double))

	#ls >> id1 >> id2 >> p.x() >> p.y() >> p.theta();
	#e_list.append([node1,node2,offset,prec])			

	for i in range(nPoints):
		cVertIDs[i] = <int> v_list[i][0]
		cVertPoses[3*i] = <double> v_list[i][1][0]
		cVertPoses[3*i+1] = <double> v_list[i][1][1]
		cVertPoses[3*i+2] = <double> v_list[i][1][2]

	for i in range(nEdges):
		cVertPairs[2*i+0] = <int> e_list[i][0]
		cVertPairs[2*i+1] = <int> e_list[i][1]
		cOffsets[3*i+0] = <double> e_list[i][2][0]
		cOffsets[3*i+1] = <double> e_list[i][2][1]
		cOffsets[3*i+2] = <double> e_list[i][2][2]
		cCovars[6*i+0] = <double> e_list[i][3][0]
		cCovars[6*i+1] = <double> e_list[i][3][1]
		cCovars[6*i+2] = <double> e_list[i][3][2]
		cCovars[6*i+3] = <double> e_list[i][3][3]
		cCovars[6*i+4] = <double> e_list[i][3][4]
		cCovars[6*i+5] = <double> e_list[i][3][5]

	v_list2 = []
	e_list2 = []

	runAlg(nPoints,nEdges,cVertIDs,cVertPoses,cVertPairs,cOffsets,cCovars)

	#print cVertIDs[0]

	for i in range(nPoints):
		nodeID = cVertIDs[i]
		pose = [cVertPoses[3*i], cVertPoses[3*i+1], cVertPoses[3*i+2]]
		v_list2.append([nodeID, pose])

	for i in range(nEdges):
		nodeID1 = cVertPairs[2*i+0]
		nodeID2 = cVertPairs[2*i+1]
		offset = [cOffsets[3*i], cOffsets[3*i+1], cOffsets[3*i+2]]
		covar = [cCovars[6*i], cCovars[6*i+1], cCovars[6*i+2], cCovars[6*i+3], cCovars[6*i+4], cCovars[6*i+5]]

		e_list2.append([nodeID1,nodeID2,offset,covar])



	free(cVertPoses)
	free(cVertIDs)
	free(cVertPairs)
	free(cOffsets)
	free(cCovars)

	return (v_list2, e_list2)

"""
