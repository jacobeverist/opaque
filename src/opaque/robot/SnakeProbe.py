
import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *


import ogre.renderer.OGRE as ogre
import ogre.physics.OgreOde as OgreOde
from Servo import Servo
from math import *
from copy import *

class SnakeProbe:

	def __del(self):
		del self.control

	def __init__(self, world, R, pos, numSegs, segLength, segWidth, maxTorque, friction):
		mode = -1

		self.timer = 0.000
		self._world = world
		self._mgr = self._world.getSceneManager()
		self._space = self._world.getDefaultSpace()
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
		self.segWidth = segWidth
		self.maxTorque = maxTorque
		self.friction = friction

		self.robotParam = {}
		self.robotParam['numJoints'] = self.numSegs-1
		self.robotParam['numSegs'] = self.numSegs
		self.robotParam['segLength'] = self.segLength
		self.robotParam['segWidth'] = self.segWidth
		self.robotParam['maxTorque'] = self.maxTorque
		
		# list of bodies
		self._bodies = []
		self._geoms = []
		self._joints = []
		self._ent = []

		self.setupMyWorld(R, pos)

		self.wallEnv = 0
	
		size = ogre.Vector3(self.probe.segLength, 0.01, self.probe.segWidth)
		self.anchorMass = OgreOde.BoxMass(STD_WEIGHT*10000.0,size)
		self.normalMass = OgreOde.BoxMass(STD_WEIGHT,size)

		self.setAnchor(False)
	
	def setAnchor(self, isAnchor):
		
		if isAnchor:
			body = self._bodies[20]
			body.setMass(self.anchorMass)		
		else:
			body = self._bodies[20]
			body.setMass(self.normalMass)
						
	def setJointsRigid(self, lowJ, highJ):
		bodies = self._bodies[lowJ,highJ+2]
		
		for b in bodies:
			mass = b.getMass()
		
	def setWalls(self, wallEnv):
		self.wallEnv = wallEnv
	
	def createWall(self, points):
		self.wallEnv.createWall(points)
		
	def getWalls(self):
		if self.wallEnv == 0:
			return []
		else:
			return self.wallEnv.getWalls()
		
	def getNumJoints(self):
		return self.numSegs - 1

	def addControl(self, controlProgram):
		self.control = controlProgram

	def frameStarted(self, evt):

		self.timer += 0.001

		if self.control:
			self.control.frameStarted()

		for joint in self._joints:
			joint.frameStarted(evt)

		" update the color based on the joint error "

		for i in range(39):
			err = self._joints[i].error()
			self.segMaterials[i].setAmbient(err,0.0,1.0-err)
		
		#for i in range(16,23):
		#	self.segMaterials[i].setAmbient(0.0,1.0,0.0)

		# lets also created meshes to visualize our joint positions

		#pos = self.getActualJointPose(4)
		#self._dnodes[4].setPosition(ogre.Vector3(pos[0],0.3,pos[1]))
		#self._dnodes[4].setOrientation(ogre.Quaternion(pos[2],ogre.Vector3().UNIT_Y))

		#for i in range(len(self._joints)):
		#	pos = self.getActualJointPose(i)
		#	self._dnodes[i].setPosition(ogre.Vector3(pos[0],0.3,pos[1]))
		#	self._dnodes[i].setOrientation(ogre.Quaternion(pos[2],ogre.Vector3().UNIT_Y))


	def setupMyWorld(self, R, pos):

		self.segMaterials = []

		for i in range(self.numSegs):

			## Create the visual representation (the Ogre entity and scene node)
			name = "segment" + str(len(self._bodies))
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
			size = ogre.Vector3(self.segLength, 0.01, self.segWidth)

			## Set the mass and geometry to match the visual representation
			size *= 1.0

			## Create a body associated with the node we created
			body = OgreOde.Body(self._world) 
			node.attachObject(body)

			mass = OgreOde.BoxMass(STD_WEIGHT,size)
			#mass.setDensity(1.0,size)

			geom = OgreOde.BoxGeometry(size, self._world, self._space)
			node.setScale(size.x,size.y,size.z)
			body.setMass(mass)

			segPos = ogre.Vector3(i*self.segLength, 0.011*(i%2), 0.0)
			actPos = R*segPos

			body.setPosition(actPos + pos)
			body.setOrientation(R)

			geom.setBody(body)
			entity.setUserObject(geom)

			if len(self._bodies) > 0:

				# joint positions
				armJointPos = ogre.Vector3(self.segLength/2.0 + (i-1)*self.segLength, 0.0, 0.0)
				actJointPos = R*armJointPos

				anch = actJointPos + pos
				axis = R*ogre.Vector3(0.0,1.0,0.0)

				b1 = self._bodies[-1]
				b2 = body

				# attach the joint with reference to b2, 
				# i.e. positive angle is positive globally with respect to b2
				joint = Servo(self._world, b1, b2, anch, axis, 0.0)
				self._joints.append(joint)
				joint.setMaxTorque(self.maxTorque)

			self._bodies.append(body)
			self._geoms.append(geom)
			self._ent.append(entity)

		self._dnodes = []

		# lets also created meshes to visualize our joint positions
		#for i in range(len(self._joints)):
		#	name = "test_j" + str(i)
		#	entity = self._mgr.createEntity(name, "Cube.mesh")
		#	entity.setMaterialName("WallMaterial_Int")
		#	node = self._mgr.getRootSceneNode().createChildSceneNode(name)
		#	node.setScale(size.x*2,size.y,size.z)

		#	pos = self.getActualJointPose(i)
		#	node.setPosition(ogre.Vector3(pos[0],0.3,pos[1]))
		#	node.setOrientation(ogre.Quaternion(pos[2],ogre.Vector3().UNIT_Y))

		#	node.attachObject(entity)

		#	self._dnodes.append(node)

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
		probeState = {}
		
		joints = []
		cmdJoints = []
		torques = []
		errors = []
		
		for i in range(self.getNumJoints()):
			joints.append(self.getServo(i))
			cmdJoints.append(self.getSeroCmd(i))
			torques.append(self.getJointTorque())
			errors.append(self.getError(i))
		
		probeState['joints'] = joints
		probeState['cmdJoints'] = cmdJoints
		probeState['torques'] = torques
		probeState['errors'] = errors

		return probeState
	
	def getServo(self, index):
		if index < self.numSegs :
			return self._joints[index].phi

		raise

	def getServoCmd(self, index):
		if index < self.numSegs:
			#return ogre.Math.RadiansToAngleUnits(self._joints[index].goal_phi)
			return self._joints[index].goal_phi

		raise

	def getWorldPose(self, i):
		currPos = self._bodies[i].getPosition()
		return currPos

	def setJointTorque(self, i, torque):
		self._joints[i].setMaxTorque(torque)

	def getJointTorque(self, i):
		return self._joints[i].getMaxTorque()

	def setServo(self, i, angle):
		self._joints[i].go_to(ogre.Math.AngleUnitsToRadians(angle))

	def getActualSegPose(self, i):

		currPos = self._bodies[i].getPosition()
		angQuat = self._bodies[i].getOrientation()

		vec = ogre.Vector3(0.0,0.0,0.0)
		radAngle = ogre.Radian(0.0)
		angQuat.ToAngleAxis(radAngle,vec)
		curAngle = radAngle.valueRadians()
		curAngle = self.normalizeAngle(curAngle)

		if vec[1] < 0:
			pose = [currPos[0], currPos[2], curAngle]
		else:
			pose = [currPos[0], currPos[2], -curAngle]
		#pose = [currPos[0], currPos[2], -curAngle]  # angle negated for display reasons

		return pose

	# gets pose of joint affixed to the i'th segment, should it be the i+1'th segment instead?
	def getActualJointPose(self, i):
		
		if i == 39:
			currPos = self._bodies[i].getPosition()
			angQuat = self._bodies[i].getOrientation()
			
			vec = ogre.Vector3(0.0,0.0,0.0)
			radAngle = ogre.Radian(0.0)
			val = angQuat.ToAngleAxis(radAngle,vec)
			curAngle = radAngle.valueRadians()
			curAngle = self.normalizeAngle(curAngle)
			
			if vec[1] < 0:
				pose = [currPos[0], currPos[2], curAngle]
			else:
				pose = [currPos[0], currPos[2], -curAngle]
			
			pose[0] = pose[0] + (self.segLength/2)*cos(pose[2])
			pose[1] = pose[1] + (self.segLength/2)*sin(pose[2])
			
			return pose

		# we use the i+1'th body as the fixed frame of reference for the joint i
		# NEGATE any angle going in or coming out of a Quaternion
		currPos = self._bodies[i+1].getPosition()
		angQuat = self._bodies[i+1].getOrientation()

		vec = ogre.Vector3(0.0,0.0,0.0)
		radAngle = ogre.Radian(0.0)
		val = angQuat.ToAngleAxis(radAngle,vec)
		curAngle = radAngle.valueRadians()
		curAngle = self.normalizeAngle(curAngle)

		# curAngle is negated because the returned angle of the body is upside down wrt to the global frame
		if vec[1] < 0:
			pose = [currPos[0], currPos[2], curAngle]
		else:
			pose = [currPos[0], currPos[2], -curAngle]

		# adding constant offset so that we are on the joint
		pose[0] = pose[0] - (self.segLength/2)*cos(pose[2])
		pose[1] = pose[1] - (self.segLength/2)*sin(pose[2])

		return pose

	def getJointWRTJointPose(self, originPose, originJoint, targetJoint):

		targetPose = [0.0,0.0,0.0]

		if originJoint >= self.numSegs-1 or targetJoint >= self.numSegs or originJoint < 0 or targetJoint < -1:
			print "ERROR: getJointWRTJointPose joint out of range!" 
			print "Received joint", originJoint , "and" , targetJoint , "with pose"
			print originPose[0], originPose[1], originPose[2]
			raise

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
