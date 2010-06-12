# /*
# SnakeApp.cpp
# --------------------------
# The main applicatin that handles switching between the
# different scenes in this demo, as well as any common 
# setup and input handling.
# */

ANY_QUERY_MASK                  = 1<<0
STATIC_GEOMETRY_QUERY_MASK      = 1<<4
WORLD_STEP = 0.001

import ogre.renderer.OGRE as ogre
import ogre.physics.OgreOde as OgreOde
import ogre.io.OIS as OIS
import sf_OIS as sf

from opaque.robot.SnakeProbe import SnakeProbe
from opaque.TestModular import TestModular
from opaque.TestHoldPosition import TestHoldPosition
from opaque.TestHoldTransition import TestHoldTransition
from opaque.TestHoldTransition import TestHoldTransition
from opaque.TestHoldSlideTransition import TestHoldSlideTransition
from opaque.TestAnchorTransition import TestAnchorTransition
from opaque.TestFrontExtend import TestFrontExtend
from opaque.TestPokeWalls import TestPokeWalls
from opaque.TestAnchor import TestAnchor
from opaque.environment.WallSections import WallSections
from DrawThings import DrawThings


import traceback 
import cProfile
import sys
from math import cos, sin, pi
from copy import copy


targetMapFile = " "
probeDisplacement = 0.0

## We need a frame listener class
class SnakeListener (OgreOde.CollisionListener):

	def __init__(self, world, plane, probe):
		OgreOde.CollisionListener.__init__(self)
		self._world = world
		self._world.setCollisionListener(self)
		self._plane = plane
		self.probe = probe

	def collision(self, contact):

		# print "Simple Scenes Collision"
		## Check for collisions between things that are connected and ignore them
		g1 = contact.getFirstGeometry()
		g2 = contact.getSecondGeometry()

		if (g1 and g2):
			b1 = g1.getBody()
			b2 = g2.getBody()

			if (b1 and b2 and OgreOde.Joint.areConnected(b1, b2)):
				return False 

		## Set the contact friction mode
		contact.setFrictionMode(OgreOde.Contact.Flag.Flag_BothFrictionPyramids)

		
		# plane friction is very low
		if g1 == self._plane or g2 == self._plane:
			contact.setCoulombFriction( 0.001 )    ### OgreOde.Utility.Infinity)

		# wall friction is higher
		else:
			contact.setCoulombFriction( self.probe.friction )    ### OgreOde.Utility.Infinity)

		try:
			if self.probe.control.isAnchored:
				for i in range(16,23):
					if g1 == self.probe._geoms[i] and g2 == self._plane:
						contact.setCoulombFriction( 100.0 )					
	
					if g1 == self._plane and g2 == self.probe._geoms[i]:
						contact.setCoulombFriction( 100.0 )					

		except:
			pass

		## Yes, this collision is valid
		return True


class SimpleScenesFrameListener ( sf.FrameListener ):
	def __init__( self, demo, renderWindow, camera, probe ):
		self._demo = demo
		self.probe = probe
		self.rWindow = renderWindow
		self.camera = camera

		sf.FrameListener.__init__(self, renderWindow, camera)

		# slow the window navigation to a reasonable speed
		self.moveSpeed = 6.0

	def __del__(self):
		sf.FrameListener.__del__(self)

	def frameStarted(self, evt):
		result = sf.FrameListener.frameStarted(self, evt)

		#self.rWindow.writeContentsToFile("testfile.png")

		# render window
		self._demo.frameStarted(evt, self.Keyboard, self.Mouse)

		self.probe.frameStarted(evt)
		
		self.adjustCamera()
		
		return result

	def adjustCamera(self):

		pose = self.probe.getActualJointPose(19)
		xAvg = pose[0]
		yAvg = pose[1]
		
		prevPose = self.camera.getPosition()
		xPrev = prevPose[0] - 2
		yPrev = prevPose[2]
		zPrev = prevPose[1]
		
		newPose = [xPrev*0.99 + 0.01*xAvg, yPrev*0.99 + 0.01*yAvg, zPrev*0.99 + 0.01*10]

		self.camera.setPosition(newPose[0]+2,newPose[2],newPose[1])

		oriQuat = ogre.Quaternion(ogre.Radian(-pi/2.0), ogre.Vector3().UNIT_X)
		self.camera.setOrientation(oriQuat)
		
	def frameEnded(self, evt):
		self._demo.frameEnded(evt, self.Keyboard, self.Mouse)
		return sf.FrameListener.frameEnded(self, evt)

# /*
# Create the scene from an ogre point of view
# and create the common OgreOde things we'll need
# */
class SnakeApp(sf.Application):
	def __init__ ( self ):
		sf.Application.__init__(self)
		self._plane = 0
		self._stepper = 0   
		self._world = 0
		self._spot=None
		self._time_elapsed = 0.0

		#self._time_step = WORLD_STEP ## SimpleScenes::WORLD_STEP

		self._looking = self._chasing = False
		self._paused = False
		self.probe = 0

	def __del__ ( self ):
		del self._plane
		del self._stepper
		del self._world
		del self.probe

	def _createScene(self):
		
		global targetMapFile, probeDisplacement

		sceneManager = self.sceneManager
		ogre.MovableObject.setDefaultQueryFlags (ANY_QUERY_MASK)
		self.shadowtype=0

		## Set up shadowing
		sceneManager.setShadowTechnique(ogre.SHADOWTYPE_TEXTURE_MODULATIVE)
		sceneManager.setShadowColour((0.5, 0.5, 0.5))
		sceneManager.setShadowFarDistance(30)

		if self.root.getRenderSystem().getName().startswith ('direct'): 
			sceneManager.setShadowTextureSettings(1024, 2)
		else: 
			sceneManager.setShadowTextureSettings(512, 2)

		## Add some default lighting to the scene
		sceneManager.setAmbientLight( (.8, .8, .8) )
		light = sceneManager.createLight('MainLight')
		light.setPosition (0, 0, 1)
		light.CastShadows=True

		## Give us some sky
		sceneManager.setSkyBox(True,"kk3d/DesertVII", 5000, True) # Examples/SpaceSkyBox",5000,True)

		## Position and orient the camera
		self.camera.setPosition(2,10,0)
		self.camera.setNearClipDistance(0.1)

		## Create the ODE world
		self._world = OgreOde.World(sceneManager)

		self._world.setGravity( (0,-9.80665,0) )
		self._world.setCFM(9e-6) 
		self._world.setERP(0.8)

		## Create something that will step the world, but don't do it automatically
		time_scale = 1.0
		frame_rate = WORLD_STEP

		# max time step
		max_frame_time = WORLD_STEP*2

		stepModeType = OgreOde.StepHandler.QuickStep

		self._stepper = OgreOde.StepHandler(self._world, stepModeType, WORLD_STEP, max_frame_time,  time_scale)

		self.root.getSingleton().setFrameSmoothingPeriod(5.0)

		## Create a default plane to act as the ground
		self._plane = OgreOde.InfinitePlaneGeometry(ogre.Plane(ogre.Vector3(0,1,0),0),self._world, self._world.getDefaultSpace())
		s = sceneManager.createStaticGeometry("StaticFloor")
		s.setRegionDimensions((160.0, 100.0, 160.0))

		## Set the region origin so the center is at 0 world
		s.setOrigin(ogre.Vector3(0,0,0))

		## Use a load of meshes to represent the floor
		i = 0
		for z in range (-80, 80, 20 ):
			for x in range (-80, 80, 20):
				name = "Plane_" + str(i)
				i += 1

				entity = sceneManager.createEntity(name, "plane.mesh")
				entity.setQueryFlags (STATIC_GEOMETRY_QUERY_MASK)
				entity.setUserObject(self._plane)
				entity.setCastShadows(False)
				s.addEntity(entity, ogre.Vector3(x,0,z))

		s.build()

		yRot = ogre.Quaternion(ogre.Degree(180.0), ogre.Vector3().UNIT_Y)
		zRot = ogre.Quaternion(ogre.Degree(0.0), ogre.Vector3().UNIT_Z)

		pos = ogre.Vector3(1.65,0.04,0.0)

		#self.probe = SnakeProbe(self._world,yRot,pos,40,0.15,0.05,30.0*2.5,0.9)
		self.probe = SnakeProbe(self._world,yRot,pos,40,0.15,0.15,30.0*2.5,0.9)

		f = open("anchorPose.txt", 'r')
		str_f = f.read()
		poses = eval(str_f)
		self.probe.restorePose(poses)
		
		self.probe.translatePose(probeDisplacement)
		
		self.drawThings = DrawThings(self.probe.robotParam)
		self.drawThings.setSim()
		self.drawThings.setRenderView(self.renderWindow)

		currControl = TestAnchor(self.probe, self.drawThings)	
		
		currControl.setRenderWindow(self.renderWindow)
		currControl.setCamera(self.camera)
		self.probe.addControl(currControl)
		
		" junction test "
		f = open(targetMapFile, 'r')
		str_f = f.read()
		walls = eval(str_f)

		self.walls = WallSections(self._world, walls)

		" make accessible to our mapping algorithms for comparison purposes "
		self.probe.setWalls(self.walls)


		self.collisionHandler = SnakeListener(self._world, self._plane, self.probe) 

		# add coordinate system markers 
		xPnt = ogre.Vector3(5.0,0.0,0.0)
		zPnt = ogre.Vector3(0.0,0.0,5.0)

		self.xEntity = self.probe._mgr.createEntity("X_ent"+str(i), "Cube.mesh")
		self.xEntity.setCastShadows(False)
		self.xEntity.setMaterialName("Red")
		self.zEntity = self.probe._mgr.createEntity("Z_ent"+str(i), "Cube.mesh")
		self.zEntity.setCastShadows(False)
		self.zEntity.setMaterialName("Green")

		size = ogre.Vector3(0.1,0.1,0.1)
		self.probe.statNode.addEntity(self.xEntity, xPnt, scale = size)
		self.probe.statNode.addEntity(self.zEntity, zPnt, scale = size)

		self.probe.statNode.build()

	def frameStarted (self, evt, Keyboard, Mouse):
		## Set the shadow distance according to how far we are from the plane that receives them
		self.sceneManager.setShadowFarDistance((abs(self.camera.getPosition().y) + 1.0) * 3.0)

	def frameEnded(self, evt, Keyboard, Mouse):
		time = evt.timeSinceLastFrame
		sceneManager = self.sceneManager

		## Step the world and then synchronise the scene nodes with it, 
		## we could get this to do this automatically, but we 
		## can't be sure of what order the framelisteners will fire in
		if (self._stepper.step(WORLD_STEP)):
			self._world.synchronise()

	## we need to register the framelistener
	def _createFrameListener(self):
		## note we pass ourselves as the demo to the framelistener
		self.frameListener = SimpleScenesFrameListener(self, self.renderWindow, self.camera, self.probe)
		self.root.addFrameListener(self.frameListener)


# set output to "test_%04d/output.txt" % i


if __name__ == '__main__':

	global targetMapFile, probeDisplacement
	
	targetMapFile = sys.argv[1]
	probeDisplacement = float(sys.argv[2])

	try:
		#prof = Profile()

		application = SnakeApp()
		cProfile.run('application.go()', 'prof_sim')
		
		#prof.runcall(application.go)
		#prof.dump_stats("profile_info2") 
	except ogre.OgreException, e:
		traceback.print_exc()
		print e

