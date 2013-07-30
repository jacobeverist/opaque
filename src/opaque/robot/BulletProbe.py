

#import ogre.renderer.OGRE as ogre
#import ogre.physics.OgreOde as OgreOde
#from Servo import Servo
from math import *
from copy import *
from transform import Transform
from bulletsnakeprobe import BulletSnakeProbe



STD_WEIGHT = 1.11


class BulletProbe:

	def __del(self):
		del self.control

	def __init__(self, angQuat, pos, numSegs, segLength, segHeight, segWidth, maxTorque, friction):
		mode = -1

		self.transform = Transform()
		
		print "initial quaternion:", angQuat
		print "initial starting position:", pos
		
		self.bullet_snake = BulletSnakeProbe(angQuat, pos, numSegs, segLength, segHeight, segWidth, friction)
		#self.bullet_snake.frameStarted()


		#self.transform.setJoints(self.joints)

		self.timer = 0.000
		#self._world.setCollisionListener(self)

		# control program for the snake probe
		self.control = 0
		self.wallEnv = 0
		
		# reset estimated pose
		estimatedPose = [0.0,0.0,0.0]
		
		goForwardState = False
		goBackwardState = False
		goRightState = False
		goLeftState = False

		stopCount = 0
		stopTimer = False

		self.numSegs = numSegs
		self.segLength = segLength
		self.segHeight = segHeight
		self.segWidth = segWidth
		self.maxTorque = maxTorque
		self.friction = friction

		self.wallCount = 0

		self.robotParam = {}
		self.robotParam['numJoints'] = self.numSegs-1
		self.robotParam['numSegs'] = self.numSegs
		self.robotParam['segLength'] = self.segLength
		self.robotParam['segHeight'] = self.segHeight
		self.robotParam['segWidth'] = self.segWidth
		self.robotParam['maxTorque'] = self.maxTorque
		
		self.isStateChanged = True
		self.isJointChanged = True
		self.probeState = {}

		self.cmdJoints = [self.bullet_snake.getServo(i) for i in range(self.getNumJoints())]
		self.joints = [self.bullet_snake.getServoCmd(i) for i in range(self.getNumJoints())]
		self.torques = [self.maxTorque for i in range(self.getNumJoints())]
		self.errors = [fabs(self.cmdJoints[i] - self.joints[i]) for i in range(self.getNumJoints())]

		for i in range(numSegs-1):
			self.setJointTorque(i, self.torques[i])
		
		self.wallEnv = 0
	
		self.setAnchor(False)

		self.walls = []
		
	def perturbProbe(self, doForce):
		if doForce:
			print "PERTURB"
			self.bullet_snake.perturb()
	
	def setAnchor(self, isAnchor):
		
		pass
	
		"""
		if isAnchor:
			body = self._bodies[20]
			body.setMass(self.anchorMass)		
		else:
			body = self._bodies[20]
			body.setMass(self.normalMass)
		"""
	
	def setWalls(self, walls):
		self.walls = walls
	
		print "creating walls:"
		print walls
		for pnts in walls:
			print "wall: ",pnts
			self.bullet_snake.addWall(pnts)
			#self.bullet_snake.createWall(pnts)

		self.bullet_snake.createWalls()

	def getWalls(self):
		return self.walls
	
	def getNumJoints(self):
		return self.numSegs - 1

	def addControl(self, controlProgram):
		self.control = controlProgram

	def frameStarted(self):

		self.timer += 1./60.
		#self.timer += 0.001

		self.isStateChanged = True
		self.isJointChanged = True
		self.probeState = {}

		self.bullet_snake.frameStarted()

		if self.control:
			self.control.frameStarted()

		" update the color based on the joint error "

		numJoints = self.getNumJoints()
		
		for i in range(numJoints):
			self.joints[i] = self.bullet_snake.getServo(i)
			self.cmdJoints[i] = self.bullet_snake.getServoCmd(i)
			self.errors[i] = abs(self.cmdJoints[i] - self.joints[i])
		
		self.transform.setJoints(self.joints)

	def translatePose(self, transDist):
		pass
		
	def savePose(self):
		self.bullet_snake.savePose()
		
	def restorePose(self):
		
		self.bullet_snake.restorePose()

		"""
		1. layout the snake in a straight line and create the joints 
		2. re-position the segments to restore a pose
		"""

		"""
		for i in range(self.numSegs):
			
			segPose = poses[i]

			print segPose
			
			pos = ogre.Vector3(segPose[0], segPose[1], segPose[2])			
			self._bodies[i].setPosition(pos)

			vec = ogre.Vector3(segPose[4],segPose[5],segPose[6])
			radAngle = ogre.Radian(segPose[3])			
			R = ogre.Quaternion(radAngle, vec)
			self._bodies[i].setOrientation(R)
		"""

	def getPose(self):
		
		poses = []
		
		"""
		for i in range(self.numSegs):
			
			pos = self._bodies[i].getPosition()
			
			oriQuat = self._bodies[i].getOrientation()			
			vec = ogre.Vector3(0.0,0.0,0.0)
			radAngle = ogre.Radian(0.0)
			val = oriQuat.ToAngleAxis(radAngle,vec)
			curAngle = radAngle.valueRadians()
			curAngle = self.normalizeAngle(curAngle)
			
			poses.append([pos[0],pos[1],pos[2],curAngle,vec[0],vec[1],vec[2]])
			
		return poses
		"""
		
	def normalizeAngle(self,angle):
		# this function converts the angle to its equivalent
	# in the range [-pi,pi] 

		while angle>pi:
			angle=angle-2*pi

		while angle<=-pi:
			angle=angle+2*pi

		return angle

	def getError(self, index):
		return abs(self.getServoCmd(index) - self.getServo(index))

	def getProbeState(self):

		if self.isStateChanged:
			self.isStateChanged = False
			
			self.probeState = {}
			
			torques = []
			
			for i in range(self.getNumJoints()):
				torques.append(self.getJointTorque(i))
			
			self.probeState['joints'] = copy(self.joints)
			self.probeState['cmdJoints'] = copy(self.cmdJoints)
			self.probeState['torques'] = torques
			self.probeState['errors'] = copy(self.errors)

		return self.probeState
	
	def getServo(self, index):
		if index < self.numSegs :
			return self.bullet_snake.getServo(index)
		raise

	def getServoCmd(self, index):
		if index < self.numSegs:
			return self.bullet_snake.getServoCmd(index)
		raise

	def getWorldPose(self, i):
		currPos = self.bullet_snake.getGlobalPosition(i)
		#if i == 0:
		#	print "sim position:", currPos
		return currPos

	def getWorldOrientation(self, i):
		quat_vars = self.bullet_snake.getGlobalOrientation(i)
		x = quat_vars[3]
		y = quat_vars[0]
		z = quat_vars[1]
		w = quat_vars[2]
		angQuat = [x,y,z,w]
		return angQuat

	def setJointTorque(self, i, torque):
		#self.bullet_snake.setJointTorque(i, 0.0005)
		self.bullet_snake.setJointTorque(i, torque)
		
	def getJointTorque(self, i):
		return self.bullet_snake.getJointTorque(i)

	def setServo(self, i, angle):
			
		self.bullet_snake.setServo(i, angle*pi/180.0)
		#self._joints[i].go_to(ogre.Math.AngleUnitsToRadians(angle))

	def getActualSegPose(self, i, phi = 0.0):

		currPos = self.bullet_snake.getGlobalPosition(i)
		quat_vars = self.bullet_snake.getGlobalOrientation(i)
		x = quat_vars[3]
		y = quat_vars[0]
		z = quat_vars[1]
		w = quat_vars[2]
		angQuat = [x,y,z,w]

		currVec, currAngle = self.quatToAngleAxis(angQuat)

		if currVec[1] < 0:
			#pose = [currPos[0], currPos[2], curAngle - pi]
			pose = [currPos[0], currPos[2], curAngle]
		else:
			#pose = [currPos[0], currPos[2], -curAngle - pi]
			pose = [currPos[0], currPos[2], -curAngle]
		#pose = [currPos[0], currPos[2], -curAngle]  # angle negated for display reasons

		return pose

	def quatToAngleAxis(self, angQuat):

		x = angQuat[0]
		y = angQuat[1]
		z = angQuat[2]
		w = angQuat[3]

		fSqrLength = x*x+y*y+z*z

		if fSqrLength > 0.0:
			rfAngle = 2.0*acos(w)
			fInvLength = 1 / sqrt(fSqrLength)
			vec = [x*fInvLength, y*fInvLength, z*fInvLength]
		else:
			# angle is 0 (mod 2*pi), so any axis will do
			rfAngle = 0.0
			vec = [1.0,0.0,0.0]

		currAngle = self.normalizeAngle(rfAngle)

		return vec, currAngle

	" gets pose of local axis centered on i'th joint but affixed to the i+1'th segment "
	def getActualJointPose(self, i, phi = 0.0):
		
		if i == 39:
			currPos = self.bullet_snake.getGlobalPosition(i)
			quat_vars = self.bullet_snake.getGlobalOrientation(i)

			x = quat_vars[3]
			y = quat_vars[0]
			z = quat_vars[1]
			w = quat_vars[2]
			angQuat = [x,y,z,w]
			currVec, currAngle = self.quatToAngleAxis(angQuat)
			
			#vec = ogre.Vector3(0.0,0.0,0.0)
			#radAngle = ogre.Radian(0.0)
			#val = angQuat.ToAngleAxis(radAngle,vec)
			#curAngle = radAngle.valueRadians()
			#curAngle = self.normalizeAngle(curAngle)
			
			if currVec[1] < 0:
				#pose = [currPos[0], currPos[2], curAngle - pi]
				pose = [currPos[0], currPos[2], curAngle]
			else:
				#pose = [currPos[0], currPos[2], -curAngle - pi]
				pose = [currPos[0], currPos[2], -curAngle]

			
			pose[0] = pose[0] + (self.segLength/2)*cos(pose[2])
			pose[1] = pose[1] + (self.segLength/2)*sin(pose[2])
			
			return pose

		" we use the i+1'th body as the fixed frame of reference for the joint i "
		" NEGATE any angle going in or coming out of a Quaternion "
		currPos = self.bullet_snake.getGlobalPosition(i+1)
		quat_vars = self.bullet_snake.getGlobalOrientation(i+1)

		x = quat_vars[3]
		y = quat_vars[0]
		z = quat_vars[1]
		w = quat_vars[2]
		angQuat = [x,y,z,w]
		currVec, currAngle = self.quatToAngleAxis(angQuat)

		#vec = ogre.Vector3(0.0,0.0,0.0)
		#radAngle = ogre.Radian(0.0)
		#val = angQuat.ToAngleAxis(radAngle,vec)
		#curAngle = radAngle.valueRadians()
		#curAngle = self.normalizeAngle(curAngle)

		#print vec, curAngle

		# curAngle is negated because the returned angle of the body is upside down wrt to the global frame
		if currVec[1] < 0:
			#pose = [currPos[0], currPos[2], curAngle - pi]
			pose = [currPos[0], currPos[2], currAngle]
		else:
			#pose = [currPos[0], currPos[2], -curAngle - pi]
			pose = [currPos[0], currPos[2], -currAngle]

		# adding constant offset so that we are on the joint
		pose[0] = pose[0] - (self.segLength/2)*cos(pose[2])
		pose[1] = pose[1] - (self.segLength/2)*sin(pose[2])

		return pose

	def getJointPose(self, originPose, originJoint, targetJoint):
		
		#result1 = self.transform.getCJointPose(originPose,originJoint, targetJoint)
		result2 = self.transform.getJointFromJoint(originPose,originJoint, targetJoint)
		return result2


	" DEPRECATED "
	def getJointFromJoint(self, originPose, originJoint, targetJoint):

		if self.isJointChanged:
			self.isJointChanged = False
			self.jointTransforms = [[None for j in range(0,39)] for i in range(0,39)]

		targetPose = [0.0,0.0,0.0]

		probeState = self.getProbeState()
		joints = probeState['joints']

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			targetPose = copy(originPose)
			return targetPose

		# forward kinematics
		if targetJoint > originJoint:

			if self.jointTransforms[originJoint][targetJoint] != None:

				offset = self.jointTransforms[originJoint][targetJoint]

				# segment origin starts
				xTotal = originPose[0] 
				zTotal = originPose[1] 
				totalAngle = originPose[2]

				xTotal = xTotal + offset[0]
				zTotal = zTotal + offset[1]
				if i != self.numSegs-1:
					totalAngle = totalAngle + offset[2]

				totalAngle = self.normalizeAngle(totalAngle)	
				targetPose = [xTotal, zTotal, totalAngle] 

			else:
				
				" combine [originJoint, targetJoint-1] with [targetJoint-1, targetJoint] "
				if abs(targetJoint-originJoint) == 1:
					" compute it right now "
					
				else:
					
					" take two components of the dynamic programming table and combine them together for the solution "
					" store the result back in the table so it doesn't need to be recomputed "
					offset1 = self.getJointFromJoint([0.0,0.0,0.0], originJoint, targetJoint-1)
					offset2 = self.getJointFromJoint([0.0,0.0,0.0], targetJoint-1, targetJoint)
					
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
					totalAngle = self.normalizeAngle(theta_1+theta_2)
					finalOffset = [xTotal, yTotal, totalAngle]
					
					" store result back into the dynamic programming table "
					self.jointTransforms[originJoint][targetJoint] = finalOffset
					
					return targetPose


		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2] 

			ind = range(targetJoint+1, originJoint + 1) # (28, 11) 
			ind.reverse()

			for i in ind:
				totalAngle = totalAngle + joints[i]
				xTotal = xTotal - self.segLength*cos(totalAngle)
				zTotal = zTotal - self.segLength*sin(totalAngle)

			totalAngle = self.normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		return targetPose



	def getJointWRTJointPose(self, originPose, originJoint, targetJoint):

		targetPose = [0.0,0.0,0.0]

		if originJoint >= self.numSegs-1 or targetJoint >= self.numSegs or originJoint < 0 or targetJoint < -1:
			print "ERROR: getJointWRTJointPose joint out of range!" 
			print "Received joint", originJoint , "and" , targetJoint , "with pose"
			print originPose[0], originPose[1], originPose[2]
			raise

		#joints = self.getProbeState()

		# origin and target are the same, so return origin pose
		if originJoint == targetJoint :
			targetPose = copy(originPose)
			return targetPose

		# forward kinematics
		if targetJoint > originJoint:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2]

			

			#joints = range(originJoint, targetJoint)			
			#xTotal = xTotal + segLength*cos(totalAngle)
			#zTotal = zTotal + segLength*sin(totalAngle)
			#totalAngle = totalAngle - self.probe.getServo(i+1)
			#totalAngle = normalizeAngle(totalAngle)

			# 
			for i in range(originJoint+1, targetJoint+1):
				xTotal = xTotal + self.segLength*cos(totalAngle)
				zTotal = zTotal + self.segLength*sin(totalAngle)
				if i != self.numSegs-1:
					totalAngle = totalAngle - self.getServo(i)

			totalAngle = self.normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		# backward kinematics
		else:

			# segment origin starts
			xTotal = originPose[0] 
			zTotal = originPose[1] 
			totalAngle = originPose[2] 

			ind = range(targetJoint+1, originJoint + 1) # (28, 11) 
			ind.reverse()

			for i in ind:
				totalAngle = totalAngle + self.getServo(i)
				xTotal = xTotal - self.segLength*cos(totalAngle)
				zTotal = zTotal - self.segLength*sin(totalAngle)

			totalAngle = self.normalizeAngle(totalAngle)

			targetPose = [xTotal, zTotal, totalAngle] 

		return targetPose
