cdef extern from "skeletonizeCV.h":
	#int runAlg(int, int, char*, char*)
	#int runAlg(int, int, int, int, double*, double*, char*, char*)
	#int runAlg(int, int, int, double*, double*, char*, char*)
	int runAlg(int, int, int, int, double*, double*, char*)

cdef extern from "stdlib.h":
	void free(void* ptr)
	void* malloc(size_t size)
	void* realloc(void* ptr, size_t size)


#def computeMedialAxis(index, sizeX, sizeY, inputImg, resultImg, numPoints, points):
def computeMedialAxis(index, sizeX, sizeY, resultImg, numPoints, points):
	
	#sizeX, sizeY = inputImg.size 
	#inputStr = inputImg.tostring()

	cdef int nPoints = numPoints
	cdef double *cPolyX = <double *>malloc(nPoints * sizeof(double))
	cdef double *cPolyY = <double *>malloc(nPoints * sizeof(double))
	for i in range(nPoints):
		cPolyX[i] = <double> points[i][0]
		cPolyY[i] = <double> points[i][1]


	cdef int totalSize = sizeX * sizeY
	#cdef char *cInputStr = <char *>malloc(totalSize * sizeof(char))
	cdef char *cResultStr = <char *>malloc(totalSize * sizeof(char))

	cdef char charData
	cdef unsigned short int intData
	

	#imgAccess = inputImg.load()
	#for i in range(sizeX):
	#	for j in range(sizeY):
	#		intData = <unsigned short int> imgAccess[i,j]
	#		charData = <char> intData
	#		cInputStr[i*sizeX + j] = charData


	runAlg(index, sizeX, sizeY, nPoints, cPolyX, cPolyY, cResultStr)
	#runAlg(index, sizeX, sizeY, nPoints, cPolyX, cPolyY, cInputStr, cResultStr)
	

	imgAccess = resultImg.load()
	for i in range(sizeX):
		for j in range(sizeY):
			charData = cResultStr[i*sizeX + j] 
			intData = <unsigned short int> charData
			imgAccess[i,j] = intData

	#free(cInputStr)
	free(cResultStr)


	return resultImg

