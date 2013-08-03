
#from libcpp.string cimport string
#cimport stdlib
cimport libc.stdlib

cdef extern from "ProbeApp.h":

	ctypedef struct c_ProbeApp "ProbeApp":

		void addWall(int numPoints, double *points)
		void createWalls()
		void updatePose(double *positions, double *quaternions)
		void updateCamera(double *pos, double *quat)
		void render()
		void shutdown()
		void saveView(char *name)

	c_ProbeApp *new_ProbeApp "new ProbeApp" (int numSegs, double segLength, double segHeight, double segWidth)
	void del_ProbeApp "delete" (c_ProbeApp *probe)

cdef class ProbeApp:

	cdef c_ProbeApp *thisptr      # hold a C++ instance which we're wrapping
	cdef int numSegments

	def __cinit__(self, int numSegs, double segLength, double segHeight, double segWidth):
		self.thisptr = new_ProbeApp(numSegs, segLength, segHeight, segWidth)
		
		self.numSegments = numSegs

	def __dealloc__(self):
		del_ProbeApp(self.thisptr)

	def addWall(self, wallPoints):
		cdef int numPoints = len(wallPoints)
		cdef double points[1000]

		for k in range(len(wallPoints)):
			points[2*k] = wallPoints[k][0]
			points[2*k+1] = wallPoints[k][1]

		self.thisptr.addWall(numPoints, points)

	def createWalls(self):
		self.thisptr.createWalls()

	def saveView(self, name):

		self.thisptr.saveView(name)

		#cdef bytes byte_string = <bytes> name.encode('UTF-8')
		#cdef char* c_string = byte_string
		#cdef string cpp_string = c_string

		#cpp_string = name.encode('UTF-8')
		#self.thisptr.saveView(c_string)

		#stdlib.free(c_string)

		#cdef char c_name[1000]

	def render(self):
		self.thisptr.render()

	def shutdown(self):
		self.thisptr.shutdown() 

	def updatePose(self, segPositions, segQuaternions):

		cdef double points[900]
		cdef double quats[1200]

		if len(segPositions) != self.numSegments:
			print "ERROR: incorrect number of positions. received", len(segPositions), "was expecting", self.numSegments
		if len(segQuaternions) != self.numSegments:
			print "ERROR: incorrect number of quaternions. received", len(segQuaternions), "was expecting", self.numSegments

		for k in range(len(segPositions)):
			points[3*k] = segPositions[k][0]
			points[3*k+1] = segPositions[k][1]
			points[3*k+2] = segPositions[k][2]

		for k in range(len(segQuaternions)):
			quats[4*k] = segQuaternions[k][0]
			quats[4*k+1] = segQuaternions[k][1]
			quats[4*k+2] = segQuaternions[k][2]
			quats[4*k+3] = segQuaternions[k][3]

		self.thisptr.updatePose(points, quats)

	def updateCamera(self, position, quaternion):

		cdef double pos[3]
		cdef double quat[4]

		pos[0] = position[0]
		pos[1] = position[1]
		pos[2] = position[2]

		quat[0] = quaternion[0]
		quat[1] = quaternion[1]
		quat[2] = quaternion[2]
		quat[3] = quaternion[3]

		self.thisptr.updateCamera(pos, quat)


