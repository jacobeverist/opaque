
import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)
	
from common import *

import ogre.renderer.OGRE as ogre
import ogre.physics.OgreOde as OgreOde


class Servo():

	def __init__(self, world, b1, b2, anch, axis, angle):

		self._world = world
		self.b1 = b1
		self.b2 = b2
		self.anch = anch
		self.axis = axis
		self.angle = angle

		# world step size
		#self.dt = world.getStepSize()

		# whether desired angle is reached
		self._done = True

		# PID error values
		self.errsum = 0.0
		self.last_err = 0.0

		# PID gains
		self.pgain = SERVO_PGAIN
		self.igain = SERVO_IGAIN
		self.dgain = SERVO_DGAIN

		# boundaries of servo
		self.histop = HI_SERVO_STOP
		self.lostop = LO_SERVO_STOP

		# max joint velocity
		self.maxVel = SERVO_VMAX*30

		# max joint torque
		self.fmax = SERVO_FMAX*30

		# tolerance on the angle
		self.tol = SERVO_ERRMIN

		# starting angle
		self.init_phi = self.lim(angle,self.lostop+self.tol, self.histop-self.tol)

		# goal angle
		self.goal_phi = self.init_phi

		# actual angle
		self.phi = self.init_phi

		self.joint = OgreOde.HingeJoint(self._world)
		self.joint.attach(self.b2,self.b1)
		self.joint.setAnchor(self.anch)
		self.joint.setAxis(self.axis)

		self.rigid = False

		self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_MaximumForce, self.fmax) 
		self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_LowStop, self.lostop-self.init_phi) 
		self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_HighStop, self.histop-self.init_phi) 
		self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_Bounceyness, 0.0) 
		self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_CFM, 0.9e-5) 
		#self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_ERP, 0.8) 

		"""
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_MaximumForce), "=", self.fmax
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_LowStop), "=", self.lostop-self.init_phi
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_HighStop), "=", self.histop-self.init_phi 
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_Bounceyness), "=", 0.0
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_CFM), "=", 0.9e-5
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_ERP), "=", 0.8
		print self.joint.getParameter(OgreOde.Joint.Parameter.Parameter_ERP), "=", 0.8
		print
		print int(OgreOde.Joint.Parameter.Parameter_ERP)
		print int(OgreOde.Joint.Parameter.Parameter_CFM)
		"""
		#pos = body.getPosition()
		#R = body.getOrientation()
		#print self.b1.getPosition()
		#print self.b1.getOrientation()
		
	def go_to(self, angle):
		if abs(self.phi-angle) > self.tol :
			self.goal_phi = self.lim(angle, self.lostop, self.histop)
			self._done = False

	def error(self):
		return self.goal_phi - (self.joint.getAngle() + self.init_phi)

	def lim(self, value, low, high):
		if value > high:
			return high

		if value < low:
			return low

		return value

	def setRigid(self):
		self.rigid = True
		
		self.joint.detach()
		self.joint = OgreOde.HingeJoint(self._world)
		self.joint.attach(self.b2,self.b1)
		
		pos = body.getPosition()
		R = body.getOrientation()
		
		#axis = R*ogre.Vector3(0.0,1.0,0.0)
	
		self.joint.setAnchor(self.anch)
		self.joint.setAxis(self.axis)
	
	def setRegular(self):
		self.rigid = False

	def frameStarted(self, evt):
		
		if not self.rigid:
			self.phi = self.joint.getAngle()
	
			#print "STOPS = ", HI_SERVO_STOP, LO_SERVO_STOP
	
			dt = WORLD_STEP
			vel = 0.0
	
			self._done = False
	
			if not self._done:
			#if True:
				err = self.error()
				vel = 0.0
	
				if abs(err) > self.tol:
					self.errsum += dt*err
					pterm = self.pgain*err
					iterm = self.igain*self.errsum
					dterm = self.dgain*(err-self.last_err)/dt
					self.last_err = err
					vel = self.lim(pterm+iterm+dterm, -self.maxVel, self.maxVel)
	
				else:
					self._done=True
					self.errsum=0.0
					self.last_err=0.0
	
			self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_MotorVelocity, vel) 

		else:
			
			pass
		
		
	def setMaxTorque(self, torque):
		self.joint.setParameter(OgreOde.Joint.Parameter.Parameter_MaximumForce, torque) 

	def setMaxVelocity(self, vel):
		self.maxVel = vel

