cdef extern from "RefNode.h":
	ctypedef struct c_RefNode "RefNode":

		double getRefPoseX()
		double getRefPoseZ()
		double getRefPoseP()
		int getNodeID()
		int getJointID()
		void setGroundTruthPose(double x, double z, double p)
		double getGroundTruthPoseX()
		double getGroundTruthPoseZ()
		double getGroundTruthPoseP()
		double computeStabilityError(double *joints)
		double getStabilityError()
		double getMaxStabilityError()
		double getAvgStabilityError()
		int isMaxErrorReached()
		void updateTime()
		void setNewRoot(int isRoot)
		int isNewRoot()
	
	c_RefNode *new_RefNode "new RefNode" (int nid, int jid, double x, double z, double p, int numJoints, double *joints)
	void del_RefNode "delete" (c_RefNode *refNode)
	
cdef class RefNode:
	cdef c_RefNode *thisptr      # hold a C++ instance which we're wrapping
	cdef int numJ
	cdef double js[100]

	#def __cinit__(self, double thresh, int sample_size):
	#	self.thisptr = new_ValueStability(thresh, sample_size)
	
	def __cinit__(self, int nid, int jid, double x, double z, double p, int numJoints, joints):
		

		self.numJ = numJoints
		
		#cdef double[numJoints] js

		for i in range(self.numJ):
			self.js[i] = joints[i]
			
		self.thisptr = new_RefNode(nid, jid, x, z, p, numJoints, self.js)

		
	def __dealloc__(self):
		del_RefNode(self.thisptr)

	def getRefPoseX(self):
		return self.thisptr.getRefPoseX()
	
	def getRefPoseZ(self):
		return self.thisptr.getRefPoseZ()

	def getRefPoseP(self):
		return self.thisptr.getRefPoseP()
	
	def getNodeID(self):
		return self.thisptr.getNodeID()

	def getJointID(self):
		return self.thisptr.getJointID()

	def setGroundTruthPose(self, x, z, p):
		self.thisptr.setGroundTruthPose(x,z,p)

	def getGroundTruthPoseX(self):
		return self.thisptr.getGroundTruthPoseX()

	def getGroundTruthPoseZ(self):
		return self.thisptr.getGroundTruthPoseZ()

	def getGroundTruthPoseP(self):
		return self.thisptr.getGroundTruthPoseP()
	
	def computeStabilityError(self, joints):

		#cdef double js[100]
		
		#cdef double[numJoints] js
		for i in range(self.numJ):
			self.js[i] = joints[i]
			
		return self.thisptr.computeStabilityError(self.js)

	def getStabilityError(self):
		return self.thisptr.getStabilityError()

	def getMaxStabilityError(self):
		return self.thisptr.getMaxStabilityError()

	def getAvgStabilityError(self):
		return self.thisptr.getAvgStabilityError()

	def isMaxErrorReached(self):
		return self.thisptr.isMaxErrorReached()

	def updateTime(self):
		self.thisptr.updateTime()

	def setNewRoot(self, isRoot):
		self.thisptr.setNewRoot(isRoot)
	
	def isNewRoot(self):
		return self.thisptr.isNewRoot()

	