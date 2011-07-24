cdef extern from "skeletonizeCV.h":
	int runAlg(int, int, char*, char*)

cdef extern from "stdlib.h":
	void free(void* ptr)
	void* malloc(size_t size)
	void* realloc(void* ptr, size_t size)


def computeMedialAxis(inputImg, resultImg):
	
	sizeX, sizeY = inputImg.size 
	inputStr = inputImg.tostring()


	cdef int totalSize = sizeX * sizeY


	cdef char *cInputStr = <char *>malloc(totalSize * sizeof(char))
	cdef char *cResultStr = <char *>malloc(totalSize * sizeof(char))

	cdef char charData
	cdef unsigned short int intData
	

	imgAccess = inputImg.load()
	for i in range(sizeX):
		for j in range(sizeY):
			intData = <unsigned short int> imgAccess[i,j]
			charData = <char> intData
			cInputStr[i*sizeX + j] = charData

	runAlg(sizeX, sizeY, cInputStr, cResultStr)

	imgAccess = resultImg.load()
	for i in range(sizeX):
		for j in range(sizeY):
			charData = cResultStr[i*sizeX + j] 
			intData = <unsigned short int> charData
			imgAccess[i,j] = intData

	#cdef char cInputStr[totalSize]
	#cdef char cResultStr[totalSize]

	free(cInputStr)
	free(cResultStr)


	return resultImg

