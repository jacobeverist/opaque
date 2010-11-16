cdef extern from "NxSnake.h":
	ctypedef struct c_NxSnake "NxSnake":

		void Step()
		void createWall(int numPoints, double *points)
		void setServo(int i, double angle)
		double getServo(int i)
		double getServoCmd(int i)
		void setJointTorque(int i, double torque)
		double getJointTorque(int i)
		void *getGlobalPosition(int i, double *x, double *y, double *z)
		void getGlobalOrientationQuat(int i, double *x, double *y, double *z, double *w)
		void savePose()
		void restorePose()

		void addWall(int numPoints, double *points)
		void createWalls()
		

	c_NxSnake *new_NxSnake "new NxSnake" (double *quatR, double *pos, int numSegs, double segLength, double segHeight, double segWidth, double friction)
	void del_NxSnake "delete" (c_NxSnake *refNode)
	
cdef class NxSnakeProbe:
	cdef c_NxSnake *thisptr      # hold a C++ instance which we're wrapping
	cdef int numJ
	cdef double js[100]

	#def __cinit__(self, double thresh, int sample_size):
	#	self.thisptr = new_ValueStability(thresh, sample_size)
	
	def __cinit__(self, quatR, pos, int numSegs, double segLength, double segHeight, double segWidth, double friction):
		
		self.numJ = numSegs

		cdef double c_quatR[4]
		c_quatR[0] = quatR[0]
		c_quatR[1] = quatR[1]
		c_quatR[2] = quatR[2]
		c_quatR[3] = quatR[3]

		cdef double c_pos[3]
		c_pos[0] = pos[0]
		c_pos[1] = pos[1]
		c_pos[2] = pos[2]

		for i in range(self.numJ):
			self.js[i] = 0.0

		self.thisptr = new_NxSnake(c_quatR, c_pos, numSegs, segLength, segHeight, segWidth, friction)
		
	def __dealloc__(self):
		del_NxSnake(self.thisptr)

	def frameStarted(self):
		self.thisptr.Step()

	def getServoCmd(self, int i):
		return self.thisptr.getServoCmd(i)

	def getServo(self, int i):
		return self.thisptr.getServo(i)

	def setServo(self, int i, double angle):
		self.thisptr.setServo(i, angle)

	def setJointTorque(self, int i, double torque):
		self.thisptr.setJointTorque(i, torque)

	def getJointTorque(self, int i):
		return self.thisptr.getJointTorque(i)

	def addWall(self, points):

		numPoints = len(points)
		cdef double cPoints[1000]
		for i in range(numPoints):
			cPoints[i*2] = points[i][0]
			cPoints[i*2+1] = points[i][1]

		self.thisptr.addWall(numPoints, cPoints)

	def createWalls(self):
		self.thisptr.createWalls()

	def createWall(self, points):

		numPoints = len(points)
		cdef double cPoints[1000]
		for i in range(numPoints):
			cPoints[i*2] = points[i][0]
			cPoints[i*2+1] = points[i][1]

		self.thisptr.createWall(numPoints, cPoints)
		
	def getGlobalPosition(self, int i):

		cdef double coordx
		cdef double coordy
		cdef double coordz

		self.thisptr.getGlobalPosition(i, &coordx, &coordy, &coordz)

		py_coord = [coordx, coordy, coordz]

		return py_coord

	def getGlobalOrientation(self, int i):

		cdef double x
		cdef double y
		cdef double z
		cdef double w

		self.thisptr.getGlobalOrientationQuat(i, &x, &y, &z, &w)

		py_coord = [x, y, z, w]

		return py_coord

	def savePose(self):
		self.thisptr.savePose()

	def restorePose(self):
		self.thisptr.restorePose()



