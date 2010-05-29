cdef extern from "math.h":
	double sin(double)
	double cos(double)

cdef float pi = 3.14159265


cdef class Transform:
	cdef double joints[39]
	cdef double jointTransforms[39][39][3]

	def setJoints(self, joints):
		cdef int i
		for i in range(39):
			self.joints[i] = joints[i]

	cpdef double computeTransform(self, int originJoint, int targetJoint) except *:
		return 0

cdef double normalizeAngle(double angle):
	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle


def getCJointPose(joints, originPose, originJoint, targetJoint):

	cdef int numSegs = 40
	cdef double segLength = 0.15
	cdef double xTotal
	cdef double zTotal
	cdef double totalAngle
	cdef int i

	cdef int oj = originJoint
	cdef int tj = targetJoint
	cdef double cj[39]

	for i in range(39):
		cj[i] = joints[i]

	targetPose = [0.0,0.0,0.0]

	# origin and target are the same, so return origin pose
	if oj == tj :
		targetPose[0] = originPose[0]
		targetPose[1] = originPose[1]
		targetPose[2] = originPose[2]
		return targetPose

	# forward kinematics
	if tj > oj:

		# segment origin starts
		xTotal = originPose[0] 
		zTotal = originPose[1] 
		totalAngle = originPose[2]

		for i in range(oj+1, tj+1):
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
			if i != numSegs-1:
				totalAngle = totalAngle - cj[i]

		totalAngle = normalizeAngle(totalAngle)

		targetPose[0] = xTotal
		targetPose[1] = zTotal
		targetPose[2] = totalAngle

	# backward kinematics
	else:

		# segment origin starts
		xTotal = originPose[0] 
		zTotal = originPose[1] 
		totalAngle = originPose[2] 

		#ind = range(tj+1, oj + 1) # (28, 11) 
		#ind.reverse()

		#for i in ind:
		for i from oj >= i > tj:

			totalAngle = totalAngle + cj[i]
			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

		totalAngle = normalizeAngle(totalAngle)

		targetPose[0] = xTotal
		targetPose[1] = zTotal
		targetPose[2] = totalAngle

	return targetPose

