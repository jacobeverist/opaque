cdef extern from "math.h":
	double sin(double)
	double cos(double)

cdef float pi = 3.14159265

cdef inline double normalizeAngle(double angle):
	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle

cdef inline void copyToFrom(double *vec2, double *vec1):
		vec2[0] = vec1[0]
		vec2[1] = vec1[1]
		vec2[2] = vec1[2]


cdef class Transform:
	cdef double joints[40]
	cdef double jointTransforms[40][40][3]
	cdef int transformComputed[40][40]
	cdef double segLength
	cdef int numSegs

	def setJoints(self, joints):
		cdef int i, j
		
		self.numSegs = 40
		self.segLength = 0.15
		

		for i in range(39):
			self.joints[i+1] = joints[i]

		#print "set joint 36 to:", self.joints[36]
			
		for i in range(40):
			for j in range(40):
				self.transformComputed[i][j] = 0

	
	def getCJointPose(self, originPose, originJoint, targetJoint):
	
		cdef double xTotal
		cdef double zTotal
		cdef double totalAngle
		cdef int i
	
		cdef int oj = originJoint + 1
		cdef int tj = targetJoint + 1
	
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
				xTotal = xTotal + self.segLength*cos(totalAngle)
				zTotal = zTotal + self.segLength*sin(totalAngle)
				#if i != self.numSegs-1:
				if i != self.numSegs:
					totalAngle = totalAngle - self.joints[i]
	
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
	
				totalAngle = totalAngle + self.joints[i]
				xTotal = xTotal - self.segLength*cos(totalAngle)
				zTotal = zTotal - self.segLength*sin(totalAngle)
	
			totalAngle = normalizeAngle(totalAngle)
	
			targetPose[0] = xTotal
			targetPose[1] = zTotal
			targetPose[2] = totalAngle
	
		return targetPose

		
	cdef computeTransform(self, int originJoint, int targetJoint):

		" NOTE:  Memory is corrupted if we attempt to reference joint -1 "

		" We exploit symmetry here.  If the normal calls are not NxN and are 1xN instead, remove symmetric computation "

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			return

		cdef int oj, tj
		#cdef double jointTransforms[40][40][3]
		#cdef int transformComputed[40][40]
		cdef double totalAngle, xTotal, yTotal, angle
		cdef double offset[3], offsetT[3], offset1[3], offset2[3], finalOffset[3]
		cdef double x1, x2, y1, y2, theta_1, theta_2

		oj = originJoint
		tj = targetJoint 

		# forward kinematics
		if tj > oj:
#
			if self.transformComputed[oj][tj] == 0:

				" combine [oj, tj-1] with [tj-1, tj] "
				if tj-oj == 1:
					" compute it right now "

					totalAngle = self.joints[tj]
					xTotal = self.segLength
					yTotal = 0.0

					#for i in range(oj+1, tj+1):
					#	xTotal = xTotal + self.segLength*cos(totalAngle)
					#	zTotal = zTotal + self.segLength*sin(totalAngle)
					#	if i != self.numSegs-1:
					#		totalAngle = totalAngle - self.joints[i]
	
					#totalAngle = normalizeAngle(totalAngle)


					offset[0] = xTotal
					offset[1] = yTotal
					offset[2] = totalAngle
					
					offsetT[0] = -xTotal*cos(totalAngle)
					offsetT[1] = -xTotal*sin(totalAngle)
					offsetT[2] = totalAngle

					" store result back into the dynamic programming table "
					self.transformComputed[oj][tj] = 1					
					copyToFrom(self.jointTransforms[oj][tj], offset)

					self.transformComputed[tj][oj] = 1
					copyToFrom(self.jointTransforms[tj][oj], offsetT)

					#if tj == 36:
					#	print oj,tj,": ", offset[0], offset[1], offset[2]
				
				else:
					
					" take two components of the dynamic programming table and combine them together for the solution "
					" store the result back in the table so it doesn't need to be recomputed "
					self.computeTransform(oj, tj-1)
					copyToFrom(offset1, self.jointTransforms[oj][tj-1])

					if self.transformComputed[tj-1][tj] == 0:
						
						totalAngle = self.joints[tj]
						xTotal = self.segLength
						yTotal = 0.0
	
						offset2[0] = xTotal
						offset2[1] = yTotal
						offset2[2] = totalAngle
						
						offsetT[0] = -xTotal*cos(totalAngle)
						offsetT[1] = -xTotal*sin(totalAngle)
						offsetT[2] = totalAngle
	
						" store result back into the dynamic programming table "
						self.transformComputed[tj-1][tj] = 1
						copyToFrom(self.jointTransforms[tj-1][tj], offset2)

						self.transformComputed[tj][tj-1] = 1
						copyToFrom(self.jointTransforms[tj][tj-1], offsetT)
						
						#offset2 = self.jointTransforms[tj-1][tj]
						#print tj-1,tj,": ", offset2
					
					else:
						#print "case2"
						copyToFrom(offset2, self.jointTransforms[tj-1][tj])
		
					
					x1 = offset1[0]
					y1 = offset1[1]
					theta_1 = offset1[2]
					
					" rotate x2,y2 by theta_1 "
					x2 = cos(theta_1) * offset2[0] + sin(theta_1) * offset2[1]
					y2 = -sin(theta_1) * offset2[0] + cos(theta_1) * offset2[1]
					theta_2 = offset2[2]
					
					" add x1+x2`, y1+y2`, theta_1 + theta_2 "
					xTotal = x1+x2
					yTotal = y1+y2

					angle = normalizeAngle(theta_1+theta_2)
						
					finalOffset[0] = xTotal
					finalOffset[1] = yTotal
					finalOffset[2] = angle
					
					offsetT[0] = -xTotal*cos(angle) + yTotal*sin(angle)
					offsetT[1] = -xTotal*sin(angle) - yTotal*cos(angle)
					offsetT[2] = angle
				
					" store result back into the dynamic programming table "
					copyToFrom(self.jointTransforms[oj][tj], finalOffset)
					self.transformComputed[oj][tj] = 1

					copyToFrom(self.jointTransforms[tj][oj], offsetT)
					self.transformComputed[tj][oj] = 1

					#if tj > 19:
					#	print oj,tj,": ", finalOffset[0], finalOffset[1], finalOffset[2], "from", x1, x2, y1, y2, theta_1, theta_2

					#	result = self.getCJointPose([0.0,0.0,0.0], oj, tj)
					#	print result[0], result[1], result[2]

		# backward kinematics
		else:

			if self.transformComputed[oj][tj] == 0:
				" oj > tj "

				" combine [oj, tj+1] with [tj+1, tj] "
				if oj-tj == 1:
					" compute it right now "

					totalAngle = self.joints[oj]
					xTotal = -self.segLength*cos(totalAngle)
					yTotal = -self.segLength*sin(totalAngle)

					offset[0] = xTotal
					offset[1] = yTotal
					offset[2] = totalAngle
					
					offsetT[0] = -xTotal*cos(totalAngle) - yTotal*sin(totalAngle)
					offsetT[1] = xTotal*sin(totalAngle) - yTotal*cos(totalAngle)
					offsetT[2] = totalAngle

					" store result back into the dynamic programming table "
					copyToFrom(self.jointTransforms[oj][tj], offset)
					self.transformComputed[oj][tj] = 1

					copyToFrom(self.jointTransforms[tj][oj], offsetT)
					self.transformComputed[tj][oj] = 1
					#print oj,tj,": ", offset
				
				else:
					
					" take two components of the dynamic programming table and combine them together for the solution "
					" store the result back in the table so it doesn't need to be recomputed "
					self.computeTransform(oj, tj+1)
					copyToFrom(offset1, self.jointTransforms[oj][tj+1])
					
					if self.transformComputed[tj+1][tj] == 0:
						
						totalAngle = self.joints[tj+1]
						xTotal = -self.segLength*cos(totalAngle)
						yTotal = -self.segLength*sin(totalAngle)
	
						offset2[0] = xTotal
						offset2[1] = yTotal
						offset2[2] = totalAngle
						
						offsetT[0] = -xTotal*cos(totalAngle) - yTotal*sin(totalAngle)
						offsetT[1] = xTotal*sin(totalAngle) - yTotal*cos(totalAngle)
						offsetT[2] = totalAngle
	
						" store result back into the dynamic programming table "
						copyToFrom(self.jointTransforms[tj+1][tj], offset2)
						self.transformComputed[tj+1][tj] = 1
						
						copyToFrom(self.jointTransforms[tj][tj+1], offsetT)
						self.transformComputed[tj][tj+1] = 1
						#offset2 = self.jointTransforms[tj-1][tj]
						#print tj+1,tj,": ", offset2
					
					else:
						#print "case2"
						copyToFrom(offset2, self.jointTransforms[tj+1][tj])
					
					
					x1 = offset1[0]
					y1 = offset1[1]
					theta_1 = offset1[2]
					
					" rotate x2,y2 by theta_1 "
					x2 = cos(theta_1) * offset2[0] - sin(theta_1) * offset2[1]
					y2 = sin(theta_1) * offset2[0] + cos(theta_1) * offset2[1]
					theta_2 = offset2[2]
					
					" add x1+x2`, y1+y2`, theta_1 + theta_2 "
					xTotal = x1+x2
					yTotal = y1+y2

					angle = normalizeAngle(theta_1+theta_2)
						
					finalOffset[0] = xTotal
					finalOffset[1] = yTotal
					finalOffset[2] = angle

					offsetT[0] = -xTotal*cos(angle) - yTotal*sin(angle)
					offsetT[1] = xTotal*sin(angle) - yTotal*cos(angle)
					offsetT[2] = angle
	
						
					" store result back into the dynamic programming table "
					copyToFrom(self.jointTransforms[oj][tj], finalOffset)
					self.transformComputed[oj][tj] = 1

					copyToFrom(self.jointTransforms[tj][oj], offsetT)
					self.transformComputed[tj][oj] = 1
					#print oj,tj,": ", finalOffset[0], finalOffset[1], finalOffset[2], "from", x1, x2, y1, y2, theta_1, theta_2
	
	def getJointFromJoint(self, originPose, originJoint, targetJoint):

		cdef int oj = originJoint + 1
		cdef int tj = targetJoint + 1
		cdef double offset[3] 
		cdef double xTotal, yTotal, totalAngle, xOff, yOff, angle

		result = [0.0, 0.0, 0.0]
		
		# origin and target are the same, so return origin pose
		if oj == tj :
			#return [originPose[0], originPose[1], originPose[2]]
			return originPose

		#return self.getCJointPose(originPose, originJoint, targetJoint)

		if self.transformComputed[oj][tj] == 0:
			self.computeTransform(oj,tj)

			" force us to use forward kinematics which takes advantage of more zeros "
			
			" FIXME:  using this snippet will double the computation time, not sure why "
			#if oj > tj:
			#	self.computeTransform(tj,oj)
			#else:
			#	self.computeTransform(oj,tj)
				
		" get transform from dynamic transform table "
		#offset = self.jointTransforms[oj][tj]
		copyToFrom(offset, self.jointTransforms[oj][tj])	
			
		# forward kinematics
		if tj > oj:

			# segment origin starts
			xTotal = originPose[0] 
			yTotal = originPose[1] 
			totalAngle = originPose[2]

			" rotate base offset by current angle "
			xOff = cos(totalAngle) * offset[0] - sin(totalAngle) * offset[1]
			yOff = sin(totalAngle) * offset[0] + cos(totalAngle) * offset[1]			

			xTotal = xTotal + xOff
			yTotal = yTotal + yOff
			#if tj != self.numSegs-1:
			if tj != self.numSegs:
				totalAngle = totalAngle - offset[2]

			angle = totalAngle
			while angle>pi:
				angle=angle-2*pi	
			while angle<=-pi:
				angle=angle+2*pi

			#result = self.getCJointPose(originPose, oj, tj)
			#print
			#print oj, tj, result[0], result[1], result[2]
			#print oj, tj, xTotal, yTotal, angle


			result[0] = xTotal
			result[1] = yTotal
			result[2] = angle

			return result
			

		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			yTotal = originPose[1] 
			totalAngle = originPose[2]
			
			" rotate base offset by current angle "
			xOff = cos(totalAngle) * offset[0] - sin(totalAngle) * offset[1]
			yOff = sin(totalAngle) * offset[0] + cos(totalAngle) * offset[1]			

			xTotal = xTotal + xOff
			yTotal = yTotal + yOff
			#if tj != self.numSegs-1:
			if tj != self.numSegs:
				totalAngle = totalAngle + offset[2]

			angle = totalAngle
			while angle>pi:
				angle=angle-2*pi	
			while angle<=-pi:
				angle=angle+2*pi

			#result = self.getCJointPose(originPose, oj, tj)
			#print
			#print oj, tj, result[0], result[1], result[2]
			#print oj, tj, xTotal, yTotal, angle
			
			result[0] = xTotal
			result[1] = yTotal
			result[2] = angle
			return result
			#print "transforming", originPose, "by", offset, "resulting in", targetPose

