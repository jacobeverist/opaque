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
#import ogre.physics.OgreOde as OgreOde
import ogre.io.OIS as OIS
import sf_OIS as sf

#from opaque.robot.SuperbotSnake import SuperbotSnake
from opaque.robot.PhysXProbe import PhysXProbe
from opaque.robot.BulletProbe import BulletProbe
from opaque.TestModular import TestModular
from opaque.TestHoldPosition import TestHoldPosition
from opaque.TestHoldTransition import TestHoldTransition
from opaque.TestHoldTransition import TestHoldTransition
from opaque.TestHoldSlideTransition import TestHoldSlideTransition
from opaque.TestAnchorTransition import TestAnchorTransition
from opaque.TestFrontExtend import TestFrontExtend
from opaque.TestPokeWalls import TestPokeWalls
#from opaque.environment.WallSections import WallSections
from DrawThings import DrawThings


import traceback 
import cProfile
import sys
from math import cos, sin, pi
from copy import copy




import opaque.control as controls


targetMapFile = ""
probeDisplacement = 0.0

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

		self.probe.perturbProbe(self.Keyboard.isKeyDown(OIS.KC_O))
		self.probe.frameStarted(evt)
		
		
		self.adjustCamera()

		self.rWindow.update()
		
		return result

	def adjustCamera(self):

		pose = self.probe.getActualJointPose(9)
		#pose = self.probe.getActualJointPose(19)
		xAvg = pose[0]
		yAvg = pose[1]
		
		prevPose = self.camera.getPosition()
		xPrev = prevPose[0] + 1
		yPrev = prevPose[2]
		zPrev = prevPose[1]
		
		newPose = [xPrev*0.99 + 0.01*xAvg, yPrev*0.99 + 0.01*yAvg, 12]

		#self.camera.setPosition(newPose[0]+1,newPose[2],newPose[1])
		#newPose[0] = -6.0
		#newPose[2] = 0.0

		self.camera.setPosition(newPose[0]-1,newPose[2],newPose[1])
		#self.camera.setPosition(-2.0,9.0,0.0)

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
		#self.camera.lookAt(-7,0.5,0)
		self.camera.setNearClipDistance(0.1)

		self.root.getSingleton().setFrameSmoothingPeriod(5.0)

		## Create a default plane to act as the ground
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
				#entity.setUserObject(self._plane)
				entity.setCastShadows(False)
				s.addEntity(entity, ogre.Vector3(x,0,z))

		s.build()

		yRot = ogre.Quaternion(ogre.Degree(180.0), ogre.Vector3().UNIT_Y)
		zRot = ogre.Quaternion(ogre.Degree(90.0), ogre.Vector3().UNIT_Z)

		#pos = ogre.Vector3(1.65,0.04,0.0)
		pos = ogre.Vector3(1.65,0.04,0.0)

		#self.probe = PhysXProbe(self.sceneManager,yRot,pos,40,0.15,0.2,0.15,30.0*2.5,0.9)
		#self.probe = PhysXProbe(self.sceneManager,yRot,pos,40,0.15,0.2,0.15,1000.0,0.9)

		
		#self.probe = PhysXProbe(self.sceneManager,yRot,pos,40,0.15,0.2,0.15,100.0,0.9)
		#self.probe = PhysXProbe(self.sceneManager,yRot,pos,40,0.15,0.1,0.15,1000.0,0.9)
		#self.probe = BulletProbe(self.sceneManager,yRot,pos,40,0.15,0.1,0.15,1000.0,1.0)
		self.probe = BulletProbe(self.sceneManager,yRot,pos,40,0.15,0.1,0.15,1000.0,1.0)
		
		#exit()
		

		#self.probe = SuperbotSnake(self._world,yRot,pos,40,0.15,0.15,30.0*2.5,0.9)
		#self.probe = SuperbotSnake(self._world,zRot,pos,40,0.15,0.15,30.0*2.5,0.9)

		self.probe.translatePose(probeDisplacement)

		#f = open("poses000160.txt", 'r')
		#str_f = f.read()
		#poses = eval(str_f)
		#print poses
		#self.probe.restorePose(poses)
		
		self.renderWindow.setAutoUpdated(False)
		
		self.drawThings = DrawThings(self.probe.robotParam)
		self.drawThings.setSim(self.sceneManager)
		self.drawThings.setRenderView(self.renderWindow)

		currControl = TestModular(self.probe, self.drawThings)	
		#currControl = TestAnchorTransition(self.probe, self.drawThings)
		#currControl = TestHoldPosition(self.probe, self.drawThings)
		#currControl = TestHoldTransition(self.probe, self.drawThings)
		#currControl = TestHoldSlideTransition(self.probe, self.drawThings)
		#currControl = TestFrontExtend(self.probe, self.drawThings)
		#currControl = TestPokeWalls(self.probe, self.drawThings)
		
		currControl.setRenderWindow(self.renderWindow)
		currControl.setCamera(self.camera)
		self.probe.addControl(currControl)
		
		"""
		WLEN = 3.0
		wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
		wall2 = [[-4.0 + WLEN*cos(pi/3), 0.2 + WLEN*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
		wall3 = [[0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), -WLEN*sin(pi/3)], [0.4*cos(pi/6) - 4.0, 0.0], [0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), WLEN*sin(pi/3)]]
		# caps to the corridors
		wall4 = [[-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)],[0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), -WLEN*sin(pi/3)]]
		wall5 = [[0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), WLEN*sin(pi/3)], [-4.0 + WLEN*cos(pi/3), 0.2 + WLEN*sin(pi/3)]]
		walls = envs.WallSections(self._world, [wall1, wall2, wall3, wall4, wall5])
		"""
		
		walls = []
		
		if targetMapFile != "":
			" junction test "
			f = open(targetMapFile, 'r')
			str_f = f.read()
			walls = eval(str_f)

		else:
			"""
			WLEN = 3.0
			wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
			wall2 = [[-4.0 + WLEN*cos(pi/3), 0.2 + WLEN*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			wall3 = [[0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), -WLEN*sin(pi/3)], [0.4*cos(pi/6) - 4.0, 0.0], [0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), WLEN*sin(pi/3)]]
			# caps to the corridors
			wall4 = [[-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)],[0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), -WLEN*sin(pi/3)]]
			wall5 = [[0.4*cos(pi/6) - 4.0 + WLEN*cos(pi/3), WLEN*sin(pi/3)], [-4.0 + WLEN*cos(pi/3), 0.2 + WLEN*sin(pi/3)]]
			walls = envs.WallSections(self._world, [wall1, wall2, wall3, wall4, wall5])
			"""
			
			" junction test "

			# Y - junction, test 1 				
			#WLEN = 3.0
			#WLEN2 = 5.0
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
			#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			#w1 = wall1[2]
			#w2 = wall2[0]

			#wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
			#lp = wall3[0]
			#rp = wall3[2]

			#wall6 = [lp, [lp[0] + WLEN*cos(pi/6), lp[1] + WLEN*sin(pi/6)]]
			#wall6.append([wall6[1][0] + 0.4*cos(pi/3), wall6[1][1] - 0.4*sin(pi/3)])
			#wall6.append([wall6[2][0] - WLEN*cos(pi/6), wall6[2][1] - WLEN*sin(pi/6)])
			#wall6.append([wall6[3][0] + WLEN*cos(pi/3), wall6[3][1] - WLEN*sin(pi/3)])
			#wall6.append([wall6[4][0] - 0.4*cos(pi/6), wall6[4][1] - 0.4*sin(pi/6)])
			#wall6.append(w1)
			#wall6.reverse()
			
			#walls = [wall1, wall2, wall3, wall6]
			
	
			# T - junctions, test 2
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
			#wall2 = [[-14.0, 0.2], [-1.0, 0.2], [-1.0,3.6], [-0.6,3.6], [-0.6,0.2], [2.0,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]

			# L - junction, left, test 3
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [-3.6,0.2]]
			#wall2 = [[-14.0, 0.2], [-3.6,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]


			# L - junction, right, test 4
			wall1 = [[-14.0, -0.2], [-3.6,-0.2], [-3.6,0.2]]
			wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,3.6], [-3.6,3.6], [-3.6,0.2]]
			wall2.reverse()
			walls = [wall1, wall2]

			# 45 degree junction, right, test 5
			#WLEN2 = 5.0
			#WLEN3 = 5.5
			#wall1 = [[-14.0, -0.2], [-3.0, -0.2]]
			#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			#w2 = wall2[0]			
			#wall3 = [[w2[0] + 0.4*cos(pi/6) - WLEN3*cos(pi/3), w2[1] - 0.4*sin(pi/6) - WLEN3*sin(pi/3)], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
			#walls = [wall1, wall2, wall3]
			
			# cross - junctions, test 6
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
			#wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,3.6], [-3.6,3.6], [-3.6,0.2], [2.0,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]
			
			
			
			#walls = [wall1]
			
			
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0, 0.2] ,[-14.0, 0.2]]
			#walls = [wall1]
	
			
			#walls = [wall1, wall2]
			for wall in walls:
				for i in range(len(wall)):
					p = copy(wall[i])
					p[0] += 6.0
					wall[i] = p
		
			" wobbly walls, curvature test "
			#walls = [[[6.89713661206, -5.66561062189, 0.24785753498134197, -1.5109569205313309], [6.8241517119, -5.44133868031, 0.23584889109403384, -1.0825876596236035], [6.90974230088, -5.1869941312, 0.26835964411070223, -1.7218194876345572], [6.74684126752, -4.99084861917, 0.25497021111463414, -0.70413835356226928], [6.7981173201, -4.74709754961, 0.24908596403833969, -1.6045479942862555], [6.6588994767, -4.55607205709, 0.23637332064884961, -0.76741953146669328], [6.51820501888, -4.37260091024, 0.23120681691414346, -0.74301192422482298], [6.51856746182, -4.14991692188, 0.22268428331215911, -1.3988363709960387], [6.42882138554, -3.95923992868, 0.21074172330951529, -0.95729865666946523], [6.4437259694, -3.73288257475, 0.22684752213712564, -1.4629592187303215], [6.36184984128, -3.54437246979, 0.20552313745644871, -0.98745975135285557], [6.23945850024, -3.37695252482, 0.20738630218531304, -0.76595030092700478], [6.09921472391, -3.22348587441, 0.20789499653697213, -0.65680032352478357], [6.12175071734, -2.99356994889, 0.23101775648241607, -1.4949150609245985], [6.07551508786, -2.79567131389, 0.20322795862302798, -1.1676927599007694], [5.80371644378, -2.72873963638, 0.27991847451904373, -0.067862959256361016], [5.69306677792, -2.57717830937, 0.18765442814230615, -0.76658840560181862], [5.71128501393, -2.34176281058, 0.23611937912146727, -1.4744424079680583], [5.54885903977, -2.22721419158, 0.19875508344200657, -0.44064375145067408], [5.38932493391, -2.11750062715, 0.19361869009713004, -0.42884391727482452], [5.31556424911, -1.94302098725, 0.18943015430245297, -0.99724826124458976], [5.26373154799, -1.74502967833, 0.20466359548535537, -1.1411621471966282], [5.07072842265, -1.67427628895, 0.20556324695108996, -0.17779135270718005], [5.0176058058, -1.46915201647, 0.21189143347691725, -1.1437984885315977], [4.87296885201, -1.3538613129, 0.18496430664331243, -0.49938506556454787], [4.71559083835, -1.2549333914, 0.18588860329443913, -0.38759679139575448], [4.53168633274, -1.19260035973, 0.19418103414887583, -0.15320261209418212], [4.45014159099, -0.995355535151, 0.21343632710659885, -1.0051881424526683], [4.27963330182, -0.915295832213, 0.18836834317260706, -0.26539283580628281], [4.10167680739, -0.849219033873, 0.18982796735230026, -0.18194299420329721], [3.92503600282, -0.783615317058, 0.18842988483853026, -0.18201983793850754], [3.77826172908, -0.662271610035, 0.19043892108802737, -0.51724507879628634], [3.65499814439, -0.481642991134, 0.21867923832985683, -0.7983829182172768], [3.4253540408, -0.519587600577, 0.23275783059355615, 0.3373402598620559], [3.27165086034, -0.384892437471, 0.20437087523695768, -0.54599688649071842], [3.09580643386, -0.290966677883, 0.19935724376297634, -0.31699854097398589], [2.88006455621, -0.304022207808, 0.21613654164005794, 0.23402844268071704], [2.67432985484, -0.290203574913, 0.20619825888009988, 0.10652105586515677], [2.48224692693, -0.222738140667, 0.20358643375424129, -0.16418328523496706], [2.30645300363, -0.0667522417068, 0.23502149719817805, -0.55217906259275129], [2.0912130176, -0.0541290021936, 0.21560982760647221, 0.11500738722041563], [1.88500208371, 0.021004290032, 0.21947200472925255, -0.17581503454269459], [1.67530985111, 0.112274243425, 0.22869419932264154, -0.23693844192417901], [1.42864517932, -0.0313331428735, 0.28542344281648224, 0.70081362224649324], [1.21757362317, 0.130450374512, 0.26594192657693183, -0.48038173500926434], [0.976750147272, 0.0594582283283, 0.25106937560040643, 0.46025692809620605], [0.744519524022, 0.145671326421, 0.2477170980329105, -0.18188184728438828], [0.497741860969, 0.059557941025, 0.26137086701895479, 0.50932782225887885], [0.250712493309, 0.0236477085096, 0.24962582656009352, 0.31794465764303342], [0.0, 0.0125, 0.25096020741543767, 0.21802240999427072], [-8.0, 0.0125, 0.25096020741543767, 0.21802240999427072]], [[-8.0, -0.2375, 0.30300632406581263, -0.47034371381816614], [0.0, -0.2375, 0.30300632406581263, -0.47034371381816614], [0.242326830198, -0.419408053121, 0.24064868479792056, 0.31642809118957987], [0.480524658816, -0.385150441756, 0.23587852483141811, 0.29471430097558599], [0.7146749304, -0.356649059437, 0.23769857348480375, -0.19992236744789607], [0.935984748056, -0.443381847043, 0.22158207316365433, -0.0097349078800936274], [1.1538538685, -0.483775677046, 0.24160824063768785, 0.51714928592812504], [1.38134274111, -0.402391681169, 0.22265096046906507, -0.27159142873004249], [1.58229271433, -0.498269512076, 0.20935700691984246, -0.15959801725288372], [1.78013619111, -0.566740783731, 0.22279108595498778, 0.41262698724922686], [1.99659240384, -0.513990656211, 0.20805775460929393, -0.34285250930915723], [2.17751579313, -0.616727000617, 0.20108826926316684, 0.18441569681185893], [2.37859227386, -0.614549632977, 0.1986038796257289, -0.38352896670747022], [2.54716391928, -0.719559638993, 0.18635458922667528, -0.10072975809353976], [2.72655077602, -0.770041205479, 0.18341806769950758, -0.24865129671941094], [2.89385991957, -0.84520661005, 0.19278692198705649, 0.16901218226509346], [3.08664482364, -0.846088681019, 0.19549086598436111, 0.17820396315719511], [3.28213360657, -0.845186220643, 0.18051806245747723, -0.42178477827811689], [3.43159168417, -0.946423831794, 0.18109471774646069, -0.054741191865903392], [3.60798627701, -0.987414616285, 0.20319724284019314, -0.87929992815900249], [3.7085819774, -1.16396382914, 0.17393757391221121, -0.61178611609702394], [3.83157742646, -1.28695325604, 0.18624573124232888, 0.052594125862925591], [4.01646155597, -1.30943282591, 0.17250546724429636, -0.22188619650029992], [4.17565204267, -1.37588976817, 0.18508056447919094, -0.91499511645720533], [4.2614815876, -1.53986569495, 0.165952201524284, -0.63298221549521327], [4.37631699873, -1.65966954106, 0.21604325759752341, 0.14892306276706016], [4.59229454603, -1.66499760033, 0.17271259776252598, -0.80566149430804368], [4.68860704284, -1.80836267412, 0.21216691203684418, -1.2905567235224151], [4.71119220306, -2.01932406656, 0.16380659988228116, -0.63632924421549331], [4.82414646655, -2.13795768965, 0.1656166033769273, -0.59911253845720969], [4.94273266096, -2.25356985946, 0.16666976627162552, -0.68201770397903394], [5.05202830537, -2.37940019405, 0.1743386675534454, -0.5019510463153829], [5.18807701684, -2.48841725285, 0.18182893285748522, -1.146678544651335], [5.23315562003, -2.66456968786, 0.17229195064423009, -0.9538083269255645], [5.30707120059, -2.82020067194, 0.26400395734459625, -0.0079346824559171324], [5.56673757968, -2.86786051948, 0.18154091480465775, -0.97866231228880374], [5.64052176226, -3.03373094526, 0.18322726719937588, -0.86222790563125817], [5.73393552499, -3.19135740363, 0.26023647763310404, -1.7295854934297927], [5.64902280849, -3.43735101101, 0.1842817150556364, -0.89255523294068162], [5.73812378952, -3.59866054477, 0.18718461438862505, -0.97540017614375485], [5.8147592642, -3.76943842356, 0.24864872568310079, -0.33393332302588158], [6.03206633279, -3.89028471737, 0.27050262068053155, -1.807167735285254], [5.92425164292, -4.13837272684, 0.20509837900940195, -0.79419142064819559], [6.04056913131, -4.30729752361, 0.23883674579235817, -0.53669950768529495], [6.22164911436, -4.4630313699, 0.26336914478588253, -1.7813885096950792], [6.12293867426, -4.70720253003, 0.2090119733793209, -1.0866872559201393], [6.18680337838, -4.90621837007, 0.21536452414075913, -1.3136813027167624], [6.2047713195, -5.12083204938, 0.24843781283846345, -0.6663460824900026], [6.3706065987, -5.30581855435, 0.23194345387930004, -0.91392610300670363], [6.47838802736, -5.51119851792, 0.26767268404651573, -1.7668025565347625], [6.38169483395, -5.76079637876, 0.22982649873860456, -1.2190927307834423]]]
			
			print "setting walls:", repr(walls)
	
			#self.walls = WallSections(self._world, walls)
			#self.walls.createWall(wall7)
	
		" make accessible to our mapping algorithms for comparison purposes "
		self.probe.setWalls(walls)

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

	## we need to register the framelistener
	def _createFrameListener(self):
		## note we pass ourselves as the demo to the framelistener
		self.frameListener = SimpleScenesFrameListener(self, self.renderWindow, self.camera, self.probe)
		self.root.addFrameListener(self.frameListener)


# set output to "test_%04d/output.txt" % i


if __name__ == '__main__':

	try:
		#prof = Profile()
		#global targetMapFile, probeDisplacement
		
		#if len(sys.argv) > 1:
		#	targetMapFile = sys.argv[1]
		#if len(sys.argv) > 2:
		#	probeDisplacement = float(sys.argv[2])

		#print "targetMapFile =", targetMapFile

		application = SnakeApp()
		cProfile.run('application.go()', 'prof_sim')
		
		
		#prof.runcall(application.go)
		#prof.dump_stats("profile_info2") 
	except ogre.OgreException, e:
		traceback.print_exc()
		print e

