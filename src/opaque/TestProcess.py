import random
random.seed(0)

from SnakeControl import SnakeControl
from copy import *
from math import *
import time

from behaviors.AnchorTransition import AnchorTransition

from pose.AverageContacts import AverageContacts
from maps.MapGraph import MapGraph
from maps.VisualGraph import VisualGraph

import numpy
import sys

class TestProcess(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe, drawThings):
		SnakeControl.__init__(self)
	
		self.drawThings = drawThings
		self.probe = probe
		
		robotParam = {}
		robotParam['numJoints'] = self.probe.numSegs-1
		robotParam['numSegs'] = self.probe.numSegs
		robotParam['segLength'] = self.probe.segLength
		robotParam['segWidth'] = self.probe.segWidth
		robotParam['maxTorque'] = self.probe.maxTorque
		
		self.numJoints = robotParam['numJoints']
		
		self.robotParam = self.probe.robotParam
		
		self.setTimerAliasing(1)
		
		self.direction = True
		
		" get the probe state "
		probeState = self.probe.getProbeState()
		torques = probeState['torques']

		self.torques = [torques[i] for i in range(self.numJoints)]
			
		self.localState = 0
		self.globalState = 0
		self.prevTime = 0
		
		self.stepDist = 0.14
		
		# error tracking for path-following
		self.lastDist = 1e100
		self.distCount = 0
		self.targetReached = False
		
		self.isAnchored = False

		self.postureCounter = 0

		self.isCapture = False
		self.renderTransition = False
		
		self.isInitialized = False

		self.globalTimer = 1
		
	def frameStarted(self):

		self.probe.setAnchor(self.isAnchored)

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		" get the probe state "
		probeState = self.probe.getProbeState()

		if self.globalState == 0 and self.globalTimer > 2:
			
			" anchor the robot initially "

			isDone = self.doInitAnchor()
			#isDone = self.doHoldPosition()
			#isDone = False
			if isDone:
				pass
				self.globalState = 1
				
				#self.probe.savePose()
				
		elif self.globalState == 1:
			
			" create the motion model object "
			self.contacts = AverageContacts(self.probe)
			
			self.globalState = 2
			
		elif self.globalState == 2:

			" get stable reference nodes on the snake body "
			self.contacts.setMask( [1.0 for i in range(self.numJoints)] )	
			self.contacts.step(probeState)
			
			if self.contacts.isStable():
				self.globalState = 3
				self.lastPose = self.contacts.getAveragePose(0)

		elif self.globalState == 3:
			
			" create the mapping object "
			#self.mapGraph = MapGraph(self.probe, self.contacts, isStable = True)
			self.mapGraph = VisualGraph(self.probe, self.contacts)

			# cross - junctions, test 6
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
			#wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,3.6], [-3.6,3.6], [-3.6,0.2], [2.0,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]			

			# 45 degree junction, right, test 5
			#WLEN2 = 5.0
			#WLEN3 = 5.5
			#wall1 = [[-14.0, -0.2], [-3.0, -0.2]]
			#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			#w2 = wall2[0]			
			#wall3 = [[w2[0] + 0.4*cos(pi/6) - WLEN3*cos(pi/3), w2[1] - 0.4*sin(pi/6) - WLEN3*sin(pi/3)], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
			#walls = [wall1, wall2, wall3]

			# L - junction, right, test 4
			#wall1 = [[-14.0, -0.2], [-3.6,-0.2], [-3.6,0.2]]
			#wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,3.6], [-3.6,3.6], [-3.6,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]
			
			# L - junction, left, test 3
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [-3.6,0.2]]
			#wall2 = [[-14.0, 0.2], [-3.6,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]

			# T - junctions, test 2
			#wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
			#wall2 = [[-14.0, 0.2], [-1.0, 0.2], [-1.0,3.6], [-0.6,3.6], [-0.6,0.2], [2.0,0.2]]
			#wall2.reverse()
			#walls = [wall1, wall2]

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

			"""
			# triple-Y - multi-junction, test 1 				
			WLEN = 2.0
			WWID = 0.4
			wall1 = [[-11.0, -WWID/2.0], [-4.0, -WWID/2.0], [-4.0 + WLEN*cos(-pi/3), -WWID/2.0 + WLEN*sin(-pi/3)]]
			wall2 = [[-11.0, WWID/2.0], [-4.0, WWID/2.0], [-4.0 + WLEN*cos(pi/3), WWID/2.0 + WLEN*sin(pi/3)]]
			#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			w1 = wall1[2]
			w2 = wall2[2]

			wall3 = [[w1[0] + WWID*cos(pi/6), w1[1] + WWID*sin(pi/6)], [WWID*cos(pi/6) - 4, 0.0], [w2[0] + WWID*cos(pi/6), w2[1] - 0.4*sin(pi/6)]]

			wall4 = [w1, [w1[0] + WLEN*cos(-2*pi/3), w1[1] + WLEN*sin(-2*pi/3)]]			
			wall5 = [w2, [w2[0] + WLEN*cos(2*pi/3), w2[1] + WLEN*sin(2*pi/3)]]
			
			wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi/2), wall4[-1][1] + WWID*sin(-2*pi/3 + pi/2)])
			wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi)])
			wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3)])
			wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi - pi/3 + pi/2.), wall4[-1][1] + WWID*sin(-2*pi/3 + pi - pi/3 + pi/2.)])
			wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3 + pi)])

			wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi/2), wall5[-1][1] + WWID*sin(2*pi/3 - pi/2)])
			wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi)])
			wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3)])
			wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi + pi/3 - pi/2.), wall5[-1][1] + WWID*sin(2*pi/3 - pi + pi/3 - pi/2.)])
			wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3 - pi)])
			
			wall5.reverse()
			
			wall6 = [[-11.0, WWID/2.0],[-11.0, -WWID/2.0]]
			
			walls = [wall1, wall2, wall3, wall4, wall5, wall6]
			"""
			# triple-Y - multi-junction, test 1				 
			WLEN = 5.0
			WLEN2 = 1.0
			WWID = 0.4
			wall1 = [[-11.0, -WWID/2.0], [-4.0, -WWID/2.0], [-4.0 + WLEN2*cos(-pi/3), -WWID/2.0 + WLEN2*sin(-pi/3)]]
			wall2 = [[-11.0, WWID/2.0], [-4.0, WWID/2.0], [-4.0 + WLEN2*cos(pi/3), WWID/2.0 + WLEN2*sin(pi/3)]]
			#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
			w1 = wall1[2]
			w2 = wall2[2]
		
			wall3 = [[w1[0] + WWID*cos(pi/6), w1[1] + WWID*sin(pi/6)], [WWID*cos(pi/6) - 4, 0.0], [w2[0] + WWID*cos(pi/6), w2[1] - 0.4*sin(pi/6)]]
		
			wall4 = [w1, [w1[0] + WLEN*cos(-2*pi/3), w1[1] + WLEN*sin(-2*pi/3)]]			
			wall5 = [w2, [w2[0] + WLEN*cos(2*pi/3), w2[1] + WLEN*sin(2*pi/3)]]
			
			wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi/2), wall4[-1][1] + WWID*sin(-2*pi/3 + pi/2)])
			wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi)])
			wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3)])
			wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi - pi/3 + pi/2.), wall4[-1][1] + WWID*sin(-2*pi/3 + pi - pi/3 + pi/2.)])
			wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3 + pi)])
		
			wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi/2), wall5[-1][1] + WWID*sin(2*pi/3 - pi/2)])
			wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi)])
			wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3)])
			wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi + pi/3 - pi/2.), wall5[-1][1] + WWID*sin(2*pi/3 - pi + pi/3 - pi/2.)])
			wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3 - pi)])
			
			wall5.reverse()
			
			wall6 = [[-11.0, WWID/2.0],[-11.0, -WWID/2.0]]			

			walls = [wall1, wall2, wall3, wall4, wall5, wall6]

			
			for wall in walls:
				for i in range(len(wall)):
					p = copy(wall[i])
					p[0] += 6.0
					wall[i] = p
					
			#walls = [[[-8.0, -0.20000000000000129, 0.15829496853428088, 0.84475807934538483], [2.0, -0.20000000000000129, 0.15829496853428088, 0.84475807934538483], [2.1491352023878836, -0.25306588802605989, 0.15153984232658485, 1.0389829067016203], [2.2990266880820016, -0.27535638664713469, 0.16223038286434097, 1.819599366340775], [2.4298271808898519, -0.17938801578707508, 0.15096383840624542, 1.0290373057576421], [2.5789206969510499, -0.20307775987299803, 0.14567876323053131, 1.7501380557210138], [2.7020741431089608, -0.12526046901930754, 0.13976622365758457, 1.3445200902091372], [2.8401014492026047, -0.10328184194586787, 0.13939723507313898, 1.4358326374053769], [2.9751919951821204, -0.068899673122594429, 0.14365996366134806, 1.954266787449926], [3.0785612580725008, 0.030864953116767124, 0.13988332876165999, 1.4548919017936155], [3.213440714690039, 0.067944298418431859, 0.13692602878725316, 1.711764384434147], [3.3319156186123617, 0.1365914724319649, 0.13793849144120834, 1.9418785943888628], [3.4323470820014403, 0.2311459465863076, 0.13898092300153525, 1.6885437053010135], [3.5541853700828661, 0.29801244429439983, 0.13849407055161722, 1.9888055502142468], [3.6504569469441019, 0.39757343598596395, 0.14194136880222177, 1.7447962554061569], [3.7708541926763051, 0.47275225878565785, 0.15493374408444299, 2.4416972933453001], [3.8189599326710155, 0.6200285394135473, 0.14890927512833427, 1.6795497686634473], [3.9501410240855108, 0.69049479711571371, 0.1453985143878502, 1.869949195253829], [4.0628931980370986, 0.78229716775713454, 0.15972862290731449, 2.5559906934897847], [4.0948482502496084, 0.93879671284900634, 0.1442214398650824, 2.0467023927665404], [4.1889338291165537, 1.0481023740291548, 0.14661462073143006, 2.3029699818031935], [4.2532914358160454, 1.1798367487084112, 0.14665349000837216, 2.13303529124541], [4.3390233975524879, 1.2988211042702114, 0.14825684957977955, 2.2736445275740156], [4.4079797737947342, 1.4300655770028146, 0.14990506322579972, 2.329575863898401], [4.4701753105499886, 1.5664592835601803, 0.15527996594555032, 2.0196614069071757], [4.574620076903555, 1.6813642817343257, 0.156935107685024, 2.0532756106950716], [4.6762156232998926, 1.8009760410915878, 0.15469615296504483, 2.3186122102679638], [4.7419382591975037, 1.9410168739684839, 0.15593526800826429, 2.3551930211954328], [4.8029803227932, 2.0845078990564243, 0.16008719909272789, 2.5338042542350734], [4.8384789059654665, 2.2406096668430715, 0.15857818374162336, 2.4639145390823374], [4.8843552997465878, 2.3924068896517632, 0.17746788375825573, 1.8876586816356336], [5.0199704721056069, 2.5068768629708291, 0.16000923400227843, 2.3942332753479185], [5.0768127715076012, 2.6564492800960791, 0.17085379604249054, 2.7264515672074641], [5.0821009876694863, 2.8272212168771391, 0.16406015557608758, 2.1814150123179963], [5.1714593062642873, 2.9648105534754379, 0.16206331582401587, 2.4879008044408324], [5.2146097436916286, 3.1210237350117724, 0.16544111256503299, 2.1824999239069789], [5.304569643633199, 3.2598688938371998, 0.1630539520307118, 2.5244760406376274], [5.3422076318965672, 3.4185193685050059, 0.16217817663148901, 2.3409774106860004], [5.407808521103961, 3.5668375545467512, 0.16494954677525836, 2.6126462373189501], [5.4316036316570102, 3.7300617736293722, 0.16357857877427753, 2.2477083029672378], [5.5114161271278075, 3.8728480369761531, 0.1645893200142288, 2.6302380047106801], [5.5322906147038129, 4.0361082588037176, 0.16020986068805165, 2.4206966907790992], [5.5852215606145972, 4.1873217318366276, 0.17482338178471646, 1.9947933674745948], [5.7059920743664918, 4.3137246504213973, 0.17098177464150199, 2.8020233796303433], [5.6983662249185674, 4.4845362825772757, 0.15761111364541036, 2.4676837449861551], [5.7433938586802, 4.6355785791116411, 0.15965226851791614, 2.200173247773352], [5.8278246010751822, 4.7710787509991661, 0.15755726734594705, 2.2104872380457286], [5.9097638041114759, 4.9056531151971772, 0.15506991401723921, 2.5744590621699524], [5.9379757195676213, 5.0581351340046741, 0.15107649298853043, 2.3817459382264494], [5.9934039781080823, 5.1986762841667709, 0.16059340954605006, 2.7678876364252583], [5.9917210829848973, 5.3592608757489941, 0.16763721429866371, 1.8816494307951745], [6.1204714356669232, 5.466618140533777, 0.14465121162759012, 2.430276512644868], [6.1669519445929968, 5.6035981949771259, 0.1428542585695563, 2.2429408515111424], [6.2372464153189071, 5.7279606755603876, 0.16896686113514811, 2.9466839295171585], [6.2054557070094258, 5.8939099178817838, 0.14346201439344203, 2.0411263472690853], [6.2996506861664834, 6.0021164586980653, 0.13918598692270986, 2.0393809101434091], [6.3912211896713007, 6.106938125559088, 0.13032326485209694, 2.27523006759401], [6.4516534018661496, 6.222402843507675, 0.12705834457978127, 2.2310863982047202], [6.515481970207011, 6.3322651973460591, 0.1277999207292731, 2.5133630094258428], [6.5463622565168524, 6.4562782114043848, 0.12213969061854742, 2.0659074223452243], [6.6242499977560936, 6.5503611736719698, 0.11740276463677371, 2.0676251269360462], [6.6989615804384375, 6.6409237959697851, 0.12434521272756251, 2.6122013174853711], [6.7169539787702472, 6.7639603983713972, 0.11582567660715642, 2.4496540415700303], [6.7520397877239535, 6.8743441533880194, 0.11775120829763262, 1.6748827906359343], [6.8560312136140542, 6.9295811916585137, 0.11513343371564191, 2.5393127638972439], [6.8809427019688698, 7.0419872649246038, 0.10277216767099399, 1.789140193558943], [6.9656172448813916, 7.1002310648551585, 0.10376500822807315, 1.5588311726882214], [7.0622766871445659, 7.1379687006309229, 0.091663170323326687, 1.7003508190056178], [7.142107347634969, 7.1830153702013728, 0.085482104087736835, 1.9553659313035521], [7.203549957888443, 7.2424459663998491, 0.082155084929628586, 1.8278341586278517], [7.2693861860028166, 7.2915891149301046, 0.080013783229957655, 1.8585333988812356], [7.3320071322305962, 7.3413969717431424, 0.10253784276582643, 2.3630069043865563], [7.3714078655685711, 7.4360626540757551, 0.080628767817108038, 1.6359675313394673], [7.4440323889309603, 7.4710866002567773, 0.082183339215728138, 1.2560342105391045], [7.5260177680638529, 7.476787376790119, 0.081476864573314753, 1.125178294340788], [7.6073409306183519, 7.4717851080215585, 0.07783897334067133, 1.1679964625702268], [7.6851664174194907, 7.4703361860445607, 0.087553424269080388, 1.752864785656302], [7.7590542508029037, 7.5173062848218297, 0.08318335214393742, 1.0962045553820732], [7.8418978854519068, 7.5097961402130462, 0.086522095398283677, 0.97011757346009619], [7.9264002519021046, 7.4912105797849549, 0.094254199800337932, 1.3685981754980467], [8.0190979518117285, 7.5082690267092129, 0.10728041962669597, 0.51425785631405352], [8.1030296114088802, 7.4414516266359705, 0.13194302677060882, 1.6258698047665712], [8.2224469563864595, 7.4975627666341484, 0.12488410652951128, 0.3705522877293157], [8.3080046718488418, 7.4065906973571867, 0.11534234725097728, 0.88790460443053187], [8.4182393963347515, 7.3726471768137838, 0.13780208588154508, 1.2063289567397442], [8.5560146969238176, 7.375364053862425, 0.15656112974092473, 0.15188359932276799], [8.6359797887045744, 7.2407647543016544, 0.15070829386475165, 1.0454163787228334], [8.7851883034709743, 7.2195560541804848, 0.15218304346945627, 0.75991590154915956], [8.923726331025895, 7.1565727735353946, 0.16077118862014445, 0.52853109300689594], [9.0509232659495353, 7.0582451841397162, 0.17161447271571684, 0.68039324073546603], [9.2010145780653492, 6.9750338551463233, 0.19995950259128939, 0.067083509792382992], [9.2882182976850469, 6.7950913255402163, 0.20562234282265465, 0.8935933276619612], [9.4850762651058513, 6.7356986615923251, 0.20146977411271474, 0.4919581490188385], [9.6398605355468909, 6.6067339396523188, 0.2256389092973855, 0.087846188449274207], [9.7424575982461157, 6.4057693653751464, 0.23569296138441245, 0.75446001388072781], [9.9564825024588739, 6.3070550684224038, 0.23520256034857429, 0.35551120785269402], [10.115023875166019, 6.1333171119055301, 0.25274008200939213, 0.15076690785947641], [10.243870562539854, 5.9158864788132632, 0.2656251807273578, 0.59847580710674064], [10.464864403400611, 5.7685146397867673, 0.27729246173706895, 0.12527347119470961], [10.600101065133217, 5.5264358206810371, 0.28298731913627206, 0.36294078719279549], [10.792400147149959, 5.3188234755580837, 0.29587407272147648, 0.33357516177729585], [11.189231974606987, 5.657295850615256, 0.30387147160822248, 0.37779717634705201], [10.959983436997485, 5.8429241845575017, 0.29497927106538935, 0.50596897207579672], [10.774606290312104, 6.0499742957101681, 0.27791084009362488, 0.3460417943845559], [10.605806002519024, 6.2565356643759813, 0.26676044718833469, 0.30095184471824882], [10.45413366160626, 6.4647812214265006, 0.25762513664089376, 0.24530569480798251], [10.256530855188974, 6.6125517313122408, 0.24674479264186894, 0.54451010456346571], [10.083507001819168, 6.7686927597948543, 0.23306066723182628, 0.45245869367678865], [9.9443091619747896, 6.9466357070581584, 0.22591974481680896, 0.27964610179286908], [9.7698773287300114, 7.0707552132084199, 0.21408436714556983, 0.56816658547154719], [9.6327746716956106, 7.2243878622610307, 0.20591291708101966, 0.34441905815223717], [9.4881136461882605, 7.3560906000062403, 0.19563339037705663, 0.44806792835718834], [9.3258114613677474, 7.4495939209280166, 0.1873095572065466, 0.66394654754736893], [9.1646521529833684, 7.5275153855656463, 0.17900859568804292, 0.7362464739337129], [9.0482968255806071, 7.6620389650737071, 0.17786274387478002, 0.32892151393086172], [8.9073980831435922, 7.7455240973496089, 0.16377491545244882, 0.65171154378011287], [8.7374047443931726, 7.7542580636282441, 0.17021755898399807, 1.135278740116147], [8.6150721568397017, 7.843975735174709, 0.15170538080518881, 0.55382527399629389], [8.4664097863946211, 7.8536915849991917, 0.14897952249944146, 1.1213495627740482], [8.3608490181120416, 7.9769150254293244, 0.16225625433815322, 0.32416423758976748], [8.2202895466597923, 7.9723642263879082, 0.14063312122985364, 1.218976907682167], [8.0948896338884104, 8.0113195831043704, 0.13131130164603869, 0.88541423391496776], [7.9697321980127516, 8.0461148704829188, 0.12990417921962974, 0.91544747160085904], [7.8433342761824534, 7.9714689757673955, 0.14679388352686637, 1.7200632792179156], [7.720448754883555, 8.0369782877399079, 0.13925631511717501, 0.69684232073956121], [7.6099406836055117, 7.9584424488572161, 0.13557253337826708, 1.8044695092790177], [7.4777614330168838, 8.0199860012457496, 0.14580453740117705, 0.75085544506041579], [7.3821688215376291, 7.9248682312745604, 0.13485302196731178, 1.9695202006826944], [7.2840429384437888, 7.8631546827684531, 0.11591915718362246, 1.7480267082983583], [7.1547158767159686, 7.8716997475900339, 0.12960905457550567, 1.1206344964736286], [7.0538285340520197, 7.8136837337947735, 0.11637918098390702, 1.7084893113225443], [6.9627713175713399, 7.7406867933513528, 0.1167046271032541, 1.8623642889004028], [6.8615590376193225, 7.681749378401916, 0.1171219214921253, 1.7139261055042585], [6.7652883122297531, 7.6131259514799137, 0.11822532423111991, 1.8058882908829823], [6.6693549366080287, 7.5408034098284098, 0.12014059509221801, 1.8325934934672858], [6.5954050618120128, 7.4442465869223611, 0.12162156072036105, 2.1038280809789609], [6.5211022217598575, 7.3474236539560724, 0.122047500539794, 2.1028578046930608], [6.4451115137073387, 7.2503585393828285, 0.12327296612574695, 2.0931914509292415], [6.3726232028995682, 7.1483896614716294, 0.12510878173107337, 2.1394131027157384], [6.3275367994499003, 7.0256261675631269, 0.13078095890705629, 2.4054384762947656], [6.2489597468052747, 6.9236286495597117, 0.12875498779143693, 2.1009897231878951], [6.1397460869206277, 6.8362536713279081, 0.13986425678638431, 1.8613754441480761], [6.1040031725041022, 6.7025056226594391, 0.13844167166575974, 2.4962700311794759], [6.0495969258054512, 6.5778840528920473, 0.13598005490189943, 2.3457772288773469], [5.9507807704313418, 6.4734838873893299, 0.14374987693878108, 1.999481297294802], [5.9183199957302417, 6.3331816998747676, 0.14400835293678993, 2.5300454251893836], [5.8716630504342717, 6.1982349510858512, 0.14278478754072593, 2.4245305050331605], [5.7614119441156717, 6.0894853806365195, 0.15486050341319874, 1.9651538526201442], [5.7670215662975046, 5.9277081958919347, 0.16187441232296743, 2.792069317866384], [5.6623644871497589, 5.8110219444706868, 0.156744331592822, 2.0263028579227842], [5.6313726364468586, 5.661828187823426, 0.15237871187443697, 2.5525922140715811], [5.5707417787644697, 5.5229310649850483, 0.15155365926325157, 2.3458237378723812], [5.5314096960754817, 5.3742072630675963, 0.15383686809550762, 2.4988634364671478], [5.4444454157024662, 5.2419716363579134, 0.15826890734474458, 2.1756764552561365], [5.4134186974871135, 5.0874905915502016, 0.1575660193318876, 2.5592004573326803], [5.3669543963896906, 4.9376615317236041, 0.15686834749239384, 2.4566972527282678], [5.2627499904191559, 4.8075221535812709, 0.16671777339847482, 2.082232134291385], [5.2608494003453661, 4.6399238872913333, 0.1676090424351295, 2.7460685480945739], [5.1711342488108212, 4.5028958858214203, 0.16378486377472493, 2.1777211347574768], [5.1674568975594006, 4.3349150909809016, 0.1680210413831506, 2.7355202135368684], [5.1244412798481367, 4.1804735114480867, 0.16032013240910462, 2.4857691440434371], [5.0504483215162512, 4.0368148535660424, 0.16159445493943861, 2.2817846064036549], [4.9763722560938941, 3.8932613395273368, 0.16153908152926311, 2.2810292012905586], [4.9058350376451116, 3.7487371165999646, 0.1608189982534782, 2.3033540930089895], [4.8439882177689295, 3.6015321376598743, 0.15966945529295623, 2.3596605298294038], [4.8241709699603579, 3.4392601475440694, 0.1634775889437152, 2.6358863143515752], [4.7721585015475956, 3.2894522992582296, 0.15858022663126017, 2.4232348667882939], [4.6889658982808875, 3.1525643310560172, 0.16018528358379738, 2.2113153162261439], [4.6486818157216927, 3.000096783036251, 0.15769958943128434, 2.4990969190985766], [4.5929219823699849, 2.8548161149974201, 0.15561372536257587, 2.3909382078497434], [4.5391521149420182, 2.7099984072180057, 0.15447772373292265, 2.4018910700519158], [4.4857219533546697, 2.5664074682264149, 0.15320946422377646, 2.4011825959308943], [4.4158658868326235, 2.4317121178465344, 0.15173235463766888, 2.2789736800646843], [4.3686109652183163, 2.2885278254653612, 0.15078053323115637, 2.4386348006971352], [4.2597086857124911, 2.1769046397664811, 0.15594692067216784, 1.9843477118240447], [4.1879713079240091, 2.0503696595757139, 0.14545567222987635, 2.2416545120592231], [4.1655080964580424, 1.9001718320375287, 0.15186830896060954, 2.6089510596821532], [4.1031987079927621, 1.7725368538078732, 0.1420322060611674, 2.3032577941607331], [4.0334500277762082, 1.6515204468477103, 0.13967766158372971, 2.2345545450955204], [3.9073807502697004, 1.5683522369059799, 0.15103116855776702, 1.7697775949777257], [3.8231004369682235, 1.4643929597743714, 0.13383087279070488, 2.0761737907504134], [3.8229197786105242, 1.3053998391493289, 0.15899322326291618, 2.7562719508076983], [3.7069330834250458, 1.2296035380046377, 0.13855682129460728, 1.7654375659373738], [3.6714419911223768, 1.0964106300047343, 0.13784037280240966, 2.4969952043565669], [3.5324483414172736, 1.0507052512378525, 0.14631546844598803, 1.5043046117909233], [3.4877838514346777, 0.92832594303743254, 0.13027513861446813, 2.4074618968528525], [3.4112967544803414, 0.83585099812737695, 0.12000788072714148, 2.066353052453779], [3.3056115325353979, 0.77743676928478278, 0.1207542474152651, 1.691540216900451], [3.1898790011312119, 0.7385208139635745, 0.12210024735348954, 1.5109920639639185], [3.1591343269396877, 0.59643978987354473, 0.14536936540280598, 2.5443054539015977], [3.0219168848869771, 0.59885126339674177, 0.13723863015944843, 1.1690395942253959], [2.9438529075229845, 0.52093324902179328, 0.1102959723925617, 1.9710742842337705], [2.8585226926123113, 0.45279440432415213, 0.10919774600893085, 1.8604566207509297], [2.7553599177555359, 0.41987608408339816, 0.10828745967849912, 1.4954901150612716], [2.6442684992150998, 0.41219503377657818, 0.11135664240242722, 1.2556437316381615], [2.5648420256245323, 0.31466569987257076, 0.12577971091870277, 2.073956982667994], [2.4532147352057319, 0.30820409416682942, 0.11181415077949629, 1.2444329029037637], [2.3527989221686045, 0.24465065106826006, 0.11883760195153077, 1.7508738990659334], [2.2314890720909695, 0.29080146229726322, 0.12979205331203175, 0.82308267297065218], [2.1178024675927336, 0.28025716874304929, 0.11417454256048429, 1.279096078266734], [2.0, 0.20000000000000001, 0.1425434477820261, 1.7846673789603902], [-8.0, 0.20000000000000001, 0.1425434477820261, 1.7846673789603902]]]
		
			self.mapGraph.loadWalls(walls)

			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()

			#self.mapGraph.testHypotheses("testData/sensorLocalize2", 175, "cornerHypotheses_2011_06_05.txt")
			#self.mapGraph.loadSeries("resultProcess_2012_04_05", 20)
			#self.mapGraph.loadSeries("resultProcess_2012_04_04", 38)
			#self.mapGraph.loadSeries("resultProcess_2012_04_11", 52)
			#self.mapGraph.loadSeries("resultProcess_2012_04_12_2", 47)
			#self.mapGraph.restoreSeries("resultProcess_2012_05_02_2", 4)
			#self.mapGraph.restoreSeries("resultProcess_2012_05_03", 6)
			#self.mapGraph.loadSeries("resultProcess_2012_05_03", 6)

			#self.mapGraph.loadSeries("resultProcess_2012_05_04", 49)

			
			self.mapGraph.restoreSeries("resultProcess_2012_05_10", 12)


			" parent-child consistency check "
			#self.mapGraph.restoreSeries("resultProcess_2012_05_04", 49)

			" parent-child consistency check "
			#self.mapGraph.restoreSeries("resultProcess_2012_05_04", 24)

			#self.mapGraph.loadSeries("resultProcess_2012_04_13_1", 14)
			#self.mapGraph.loadSeries("resultProcess_2012_04_13_2", 16)
			#self.mapGraph.loadSeries("resultProcess_2012_04_13_3", 14)
			#self.mapGraph.loadSeries("resultProcess_2012_04_07", 38)
			
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 175)
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 50)
			#self.mapGraph.instantSensorTest("testData/sensorLocalize2", 50)

			#self.mapGraph.sensorTest("testData/sensorLocalize2", 50)
			#self.mapGraph.sensorTest("testData/sensorLocalize2", 176)

			#self.mapGraph.drawMap()
			exit()
			
			#self.mapGraph.newNode(0.0, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.mapGraph.synch()
			#self.mapGraph.saveMap()
			#self.mapGraph.poseGraph.drawConstraints(0)
			#self.mapGraph.saveLocalMap()
	
			exit()
	
			#self.mapGraph.numNodes
	
			self.mapGraph.poseGraph.correctPoses3()
			self.mapGraph.synch()
			self.mapGraph.saveMap()
			self.mapGraph.poseGraph.drawConstraints(1)
			#self.mapGraph.renderConstraints()

			exit()
			
			
			self.restState = deepcopy(probeState)
			
			#self.globalState = 6
			self.globalState = 4
			#self.globalState = 9

			#self.restState = deepcopy(probeState)
			#self.mapGraph.newNode(self.stepDist, self.direction)
			#self.mapGraph.forceUpdate(False)
			#self.globalState = 6

			print "Start Time:", time.clock()

	def doInitAnchor(self):
		
		" get the probe state "
		probeState = self.probe.getProbeState()

		" create the first behavior "
		if self.localState == 0:
			
			print "START:  doInitAnchor()"
			
			self.behavior = AnchorTransition(self.robotParam)
			
			self.localState = 1

		" immediately begin running behavior "
		if self.localState == 1:
			
			isDone = self.behavior.step(probeState)
			joints1 = self.behavior.getJoints()
			
			self.mergeJoints([joints1])
			
			if isDone:
				self.localState = 0
				self.behavior = 0
				
				self.isAnchored = True
				
				print "FINISH:  doInitAnchor()"
				return True
		
		return False
