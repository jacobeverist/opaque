

import ogre.renderer.OGRE as ogre
import ogre.physics.OgreOde as OgreOde
from Servo import Servo
from math import *
from copy import *
from transform import Transform
from nxsnakeprobe import NxSnakeProbe



STD_WEIGHT = 1.11


class PhysXProbe:

	def __del(self):
		del self.control

	def __init__(self, sceneManager, R, pos, numSegs, segLength, segHeight, segWidth, maxTorque, friction):
		mode = -1

		self.transform = Transform()
		
		print "initial quaternion:", R.x, R.y, R.z, R.w
		print "initial starting position:", pos.x, pos.y, pos.z
		
		self.nx_snake = NxSnakeProbe([R.x,R.y,R.z,R.w], [pos.x, pos.y, pos.z], numSegs, segLength, segHeight, segWidth, friction)
		#self.nx_snake.frameStarted()

		#self.transform.setJoints(self.joints)

		self.timer = 0.000
		self._mgr = sceneManager
		#self._world.setCollisionListener(self)

		self.statNode = self._mgr.createStaticGeometry("contacts")
		self.statNode.setOrigin(ogre.Vector3().ZERO)
		
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
		self.robotParam['segWidth'] = self.segWidth
		self.robotParam['maxTorque'] = self.maxTorque
		
		self.isStateChanged = True
		self.isJointChanged = True
		self.probeState = {}

		self.cmdJoints = [self.nx_snake.getServo(i) for i in range(self.getNumJoints())]
		self.joints = [self.nx_snake.getServoCmd(i) for i in range(self.getNumJoints())]
		self.torques = [self.maxTorque for i in range(self.getNumJoints())]
		self.errors = [fabs(self.cmdJoints[i] - self.joints[i]) for i in range(self.getNumJoints())]
		
		# list of bodies
		self._ent = []
		self._nodes = []

		self.setupMyWorld(R, pos)

		self.wallEnv = 0
	
		size = ogre.Vector3(self.segLength, 0.01, self.segWidth)
		self.anchorMass = OgreOde.BoxMass(STD_WEIGHT*10000.0,size)
		self.normalMass = OgreOde.BoxMass(STD_WEIGHT,size)

		self.setAnchor(False)

		self.walls = []
		self.npnts = []
		
	def perturbProbe(self, doForce):
		if doForce:
			print "PERTURB"
			self.nx_snake.perturb()
	
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
		self.npnts = []
	
		print "creating walls:"
		print walls
		for pnts in walls:
			print "wall: ",pnts
			self.createWall(pnts)
			self.nx_snake.addWall(pnts)
			#self.nx_snake.createWall(pnts)
			#self.npnts.append(pnts)

		self.nx_snake.createWalls()

	def createWall(self, points):
				
		self.npnts.append(points)
		
		# create the static geometry
		s = self._mgr.createStaticGeometry("corridor_%04u" % self.wallCount)

		self.wallCount += 1
		
		# region dimensions (FIXME: what does this mean?)
		s.setRegionDimensions((160.0, 100.0, 160.0))

		## Set the region origin so the center is at 0 world
		s.setOrigin(ogre.Vector3().ZERO)

		# start with a ManualObject for the mesh
		# create the wall1 mesh
		mo1_int = ogre.ManualObject("WallObject_Interior" + str(self.wallCount))
		mo1_int.begin("WallMaterial_Int", ogre.RenderOperation.OT_TRIANGLE_LIST)

		mo1_ext = ogre.ManualObject("WallObject_Exterior" + str(self.wallCount))
		mo1_ext.begin("WallMaterial_Ext", ogre.RenderOperation.OT_TRIANGLE_LIST)

		mo1_top = ogre.ManualObject("WallObject_Top" + str(self.wallCount))
		mo1_top.begin("WallMaterial_Top", ogre.RenderOperation.OT_TRIANGLE_LIST)

		topPoints = []
		for i in range(len(points)-1):
			vec = [points[i][0]-points[i+1][0], points[i][1]-points[i+1][1]]
			vecMag = sqrt(vec[0]**2 + vec[1]**2)
			vec[0] /= vecMag
			vec[1] /= vecMag
			
			orthoVec = [vec[0]*cos(pi/2) - vec[1]*sin(pi/2), vec[0]*sin(pi/2) + vec[1]*cos(pi/2)]
			
			topPoints.append(ogre.Vector3(points[i][0], 0.2, points[i][1]))
			topPoints.append(ogre.Vector3(points[i][0] + 0.05*orthoVec[0], 0.2, points[i][1] + 0.05*orthoVec[1]))
			topPoints.append(ogre.Vector3(points[i+1][0], 0.2, points[i+1][1]))
			topPoints.append(ogre.Vector3(points[i+1][0] + 0.05*orthoVec[0], 0.2, points[i+1][1] + 0.05*orthoVec[1]))
		
		for p in topPoints:
			mo1_top.position(p[0], p[1], p[2])

		for i in range(len(points) - 1):
			mo1_top.triangle(4*i+2, 4*i+1, 4*i)
			mo1_top.triangle(4*i+1, 4*i+2, 4*i+3)
			if i < len(points)-2:
				mo1_top.triangle(4*i+2, 4*i+5, 4*i+3)
				
		wall1Points = []
		for i in range(len(points)):
			wall1Points.append(ogre.Vector3(points[i][0], -100.0, points[i][1]))
			wall1Points.append(ogre.Vector3(points[i][0], 1.5, points[i][1]))

		# build the triangle meshes of the wall
		for p in wall1Points:
			#mo1_int.position(p[0],p[1],p[2])
			#mo1_ext.position(p[0],p[1],p[2])
			if p[1] < 0.0:
				mo1_int.position(p[0],0.0,p[2])
				mo1_ext.position(p[0],0.0,p[2])
			else:
				mo1_int.position(p[0],0.2,p[2])
				mo1_ext.position(p[0],0.2,p[2])

		# indices for ODE geometry
		indices1 = []
		indices2 = []

		for i in range(len(points) - 1):
			indices1 += [2*i+2, 2*i+1, 2*i]
			indices1 += [2*i+1, 2*i+2, 2*i+3]
			mo1_int.triangle(2*i+2, 2*i+1, 2*i)
			mo1_int.triangle(2*i+1, 2*i+2, 2*i+3)

			indices2 += [2*i, 2*i+1, 2*i+2]
			indices2 += [2*i+3, 2*i+2, 2*i+1]
			mo1_ext.triangle(2*i, 2*i+1, 2*i+2)
			mo1_ext.triangle(2*i+3, 2*i+2, 2*i+1)


		mo1_int.end()
		mo1_ext.end()
		mo1_top.end()

		# convert to mesh
		mo1_int.convertToMesh("WallMesh_Interior_" + str(self.wallCount))
		mo1_ext.convertToMesh("WallMesh_Exterior_" + str(self.wallCount))
		mo1_top.convertToMesh("WallMesh_Top_" + str(self.wallCount))
		# create Ogre entity
		print "creating entity " + str(self.wallCount)
		entity1_int = self._mgr.createEntity("WallEntity_Interior_" + str(self.wallCount), "WallMesh_Interior_" + str(self.wallCount))
		entity1_ext = self._mgr.createEntity("WallEntity_Exterior_" + str(self.wallCount), "WallMesh_Exterior_" + str(self.wallCount))
		entity1_top = self._mgr.createEntity("WallEntity_Top_" + str(self.wallCount), "WallMesh_Top_" + str(self.wallCount))

		# add ODE geometry to the entity
		entity1_int.setCastShadows(False)
		entity1_ext.setCastShadows(False)
		entity1_top.setCastShadows(False)

		# create the ODE geometry
		#self.trimeshes.append(OgreOde.makeTriangleMeshGeometry(wall1Points, len(wall1Points), indices1, len(indices1), self._world, self._space))
		#self.trimeshes.append(OgreOde.makeTriangleMeshGeometry(wall1Points, len(wall1Points), indices2, len(indices2), self._world, self._space))

		# add the entity to the Ogre static geometry
		s.addEntity(entity1_int, ogre.Vector3(0.0,0.0,0.0))
		s.addEntity(entity1_ext, ogre.Vector3(0.0,0.0,0.0))
		s.addEntity(entity1_top, ogre.Vector3(0.0,0.0,0.0))

		# now build the StaticGeometry so it will be rendered
		s.build()

	def getWalls(self):
		return self.walls
	
	def getNumJoints(self):
		return self.numSegs - 1

	def addControl(self, controlProgram):
		self.control = controlProgram

	def frameStarted(self, evt):

		self.timer += 0.001

		self.isStateChanged = True
		self.isJointChanged = True
		self.probeState = {}

		self.nx_snake.frameStarted()

		if self.control:
			self.control.frameStarted()

		" update the color based on the joint error "

		numJoints = self.getNumJoints()
		
		for i in range(numJoints):
			self.joints[i] = self.nx_snake.getServo(i)
			self.cmdJoints[i] = self.nx_snake.getServoCmd(i)
			self.errors[i] = abs(self.cmdJoints[i] - self.joints[i])
		
		self.transform.setJoints(self.joints)

		"render the positions"
		for i in range(len(self._nodes)):
			currPos = self.getWorldPose(i)
			#print i, currPos
			self._nodes[i].setPosition(currPos[0], currPos[1], currPos[2])
			self._nodes[i].setOrientation(self.getWorldOrientation(i))
			
	def translatePose(self, transDist):

		pass
		"""
		poses = []

		for i in range(self.numSegs):
			
			pos = self._bodies[i].getPosition()
			newPose = ogre.Vector3(pos[0] + transDist, pos[1], pos[2])
			self._bodies[i].setPosition(newPose)			
			
		return poses
		"""
		
	def savePose(self):
		
		self.nx_snake.savePose()
		
	def restorePose(self):
		
		self.nx_snake.restorePose()
		pass
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
		
	def setupMyWorld(self, R, pos):

		self.segMaterials = []

		for i in range(self.numSegs):

			## Create the visual representation (the Ogre entity and scene node)
			name = "segment" + str(len(self._ent))
			entity = self._mgr.createEntity(name, "Cube.mesh")
			node = self._mgr.getRootSceneNode().createChildSceneNode(name)
			node.attachObject(entity)
			entity.setCastShadows(False)

			newMaterial = entity.getSubEntity(0).getMaterial().clone("custom" + str(i))
			val = float(self.numSegs - i)/self.numSegs

			newMaterial.setAmbient(0.0,0.0,1.0)
			#newMaterial.setAmbient(val,0.0,1.0-val)
			#newMaterial.setDiffuse(0.9,0.9,0.0, 0.1)
			#newMaterial.setSpecular(0.5,0.5,0.0, 0.5)

			self.segMaterials.append(newMaterial)
			entity.getSubEntity(0).setMaterialName("custom" + str(i))

			## Pick a size
			size = ogre.Vector3(self.segLength, self.segHeight, self.segWidth)
			node.setScale(size.x,size.y,size.z)
			
			currPos = self.getWorldPose(i)
			node.setPosition(currPos[0], currPos[1], currPos[2])
			node.setOrientation(self.getWorldOrientation(i))
			
			self._ent.append(entity)
			self._nodes.append(node)

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
			return self.nx_snake.getServo(index)
		raise

	def getServoCmd(self, index):
		if index < self.numSegs:
			return self.nx_snake.getServoCmd(index)
		raise

	def getWorldPose(self, i):
		currPos = self.nx_snake.getGlobalPosition(i)
		#if i == 0:
		#	print "sim position:", currPos
		return currPos

	def getWorldOrientation(self, i):
		quat_vars = self.nx_snake.getGlobalOrientation(i)
		angQuat = ogre.Quaternion(quat_vars[3], quat_vars[0], quat_vars[1], quat_vars[2])
		return angQuat

	def setJointTorque(self, i, torque):
		#self.nx_snake.setJointTorque(i, 0.0005)
		self.nx_snake.setJointTorque(i, torque)

	def getJointTorque(self, i):
		return self.nx_snake.getJointTorque(i)

	def setServo(self, i, angle):
			
		self.nx_snake.setServo(i, ogre.Math.AngleUnitsToRadians(angle))
		#self._joints[i].go_to(ogre.Math.AngleUnitsToRadians(angle))

	def getActualSegPose(self, i, phi = 0.0):

		currPos = self.nx_snake.getGlobalPosition(i)
		quat_vars = self.nx_snake.getGlobalOrientation(i)
		angQuat = ogre.Quaternion(quat_vars[3], quat_vars[0], quat_vars[1], quat_vars[2])

		vec = ogre.Vector3(0.0,0.0,0.0)
		radAngle = ogre.Radian(0.0)
		angQuat.ToAngleAxis(radAngle,vec)
		curAngle = radAngle.valueRadians()
		curAngle = self.normalizeAngle(curAngle)

		if vec[1] < 0:
			pose = [currPos[0], currPos[2], curAngle - pi]
		else:
			pose = [currPos[0], currPos[2], -curAngle - pi]
		#pose = [currPos[0], currPos[2], -curAngle]  # angle negated for display reasons

		return pose

	" gets pose of local axis centered on i'th joint but affixed to the i+1'th segment "
	def getActualJointPose(self, i, phi = 0.0):
		
		if i == 39:
			currPos = self.nx_snake.getGlobalPosition(i)
			quat_vars = self.nx_snake.getGlobalOrientation(i)
			angQuat = ogre.Quaternion(quat_vars[3], quat_vars[0], quat_vars[1], quat_vars[2])
			
			vec = ogre.Vector3(0.0,0.0,0.0)
			radAngle = ogre.Radian(0.0)
			val = angQuat.ToAngleAxis(radAngle,vec)
			curAngle = radAngle.valueRadians()
			curAngle = self.normalizeAngle(curAngle)
			
			if vec[1] < 0:
				pose = [currPos[0], currPos[2], curAngle - pi]
			else:
				pose = [currPos[0], currPos[2], -curAngle - pi]

			
			pose[0] = pose[0] + (self.segLength/2)*cos(pose[2])
			pose[1] = pose[1] + (self.segLength/2)*sin(pose[2])
			
			return pose

		" we use the i+1'th body as the fixed frame of reference for the joint i "
		" NEGATE any angle going in or coming out of a Quaternion "
		currPos = self.nx_snake.getGlobalPosition(i+1)
		quat_vars = self.nx_snake.getGlobalOrientation(i+1)
		angQuat = ogre.Quaternion(quat_vars[3], quat_vars[0], quat_vars[1], quat_vars[2])

		vec = ogre.Vector3(0.0,0.0,0.0)
		radAngle = ogre.Radian(0.0)
		val = angQuat.ToAngleAxis(radAngle,vec)
		curAngle = radAngle.valueRadians()
		curAngle = self.normalizeAngle(curAngle)

		#print vec, curAngle

		# curAngle is negated because the returned angle of the body is upside down wrt to the global frame
		if vec[1] < 0:
			pose = [currPos[0], currPos[2], curAngle - pi]
		else:
			pose = [currPos[0], currPos[2], -curAngle - pi]

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
