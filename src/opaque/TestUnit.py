import random
random.seed(0)

from SnakeControl import SnakeControl
from copy import *
from math import *
import time

from behaviors.AnchorTransition import AnchorTransition

from pose.AverageContacts import AverageContacts
from maps.PoseGraph import PoseGraph
from maps.LocalNode import LocalNode
from maps.Pose import Pose
from maps.SplineFit import SplineFit

from maps.functions import *
import hashlib

import scipy
import pca_module

from maps import scsgp

print "foo"

def computeBareHull(node1, sweep = False, static = False):
	
	if static:
		node1.computeStaticAlphaBoundary()

		a_data = node1.getAlphaBoundary(static=True)
		a_data = decimatePoints(a_data)

		" convert hull points to GPAC coordinates before adding covariances "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
		
	else:
				
		" Read in data of Alpha-Shapes without their associated covariances "
		node1.computeAlphaBoundary(sweep = sweep)
		a_data = node1.getAlphaBoundary(sweep = sweep)
		a_data = decimatePoints(a_data)
		
		" convert hull points to GPAC coordinates "
		localGPACPose = node1.getLocalGPACPose()
		localGPACProfile = Pose(localGPACPose)
		
		a_data_GPAC = []
		for pnt in a_data:
			a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
	
	return a_data_GPAC

def computeHullAxis(nodeID, node2, tailCutOff = False):

	medial2 = node2.getBestMedialAxis()

	if tailCutOff:
		medial2 = node2.medialTailCuts[0]
	else:
		medial2 = node2.medialLongPaths[0]


	if node2.isBowtie:			
		hull2 = computeBareHull(node2, sweep = False, static = True)
		hull2.append(hull2[0])

	else:
		hull2 = computeBareHull(node2, sweep = False)
		hull2.append(hull2[0])
	
	
	return hull2, medial2

class TestUnit(SnakeControl):

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

			if isDone:
				pass
				self.globalState = 1
								
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
		
			self.testMedialHull()
			self.testSpline()

			self.testToro()

			self.testPCA()
			
			self.testSCSGP()

			exit()

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

	def testMedialHull(self):
			PIXELSIZE = 0.05

			currNode = LocalNode(self.probe, self.contacts, 10, 19, PIXELSIZE)
			currNode.readFromFile2("result_2012_08_27", 10)			
			hull1, axis1 = computeHullAxis(10, currNode)

			currNode = LocalNode(self.probe, self.contacts, 10, 19, PIXELSIZE)
			currNode.readFromFile2("result_2012_08_27", 10)			
			hull2, axis2 = computeHullAxis(10, currNode)

			hull1_md5_nom = 235711552505359876739222943177420098447
			hull2_md5_nom = 153873940548353262496568220315837306136
			axis1_md5_nom = 218829303125138152013279906890643647977
			axis2_md5_nom = 133308086809380981344520323743369451609		

			m = hashlib.md5()
			m.update(repr(hull1))
			hull1_md5 = int(m.digest().encode('hex'),16)
			m = hashlib.md5()
			m.update(repr(hull2))
			hull2_md5 = int(m.digest().encode('hex'),16)

			m = hashlib.md5()
			m.update(repr(axis1))
			axis1_md5 = int(m.digest().encode('hex'),16)
			m = hashlib.md5()
			m.update(repr(axis2))
			axis2_md5 = int(m.digest().encode('hex'),16)

			if hull1_md5_nom == hull1_md5:
				print "hull1 PASS"
			else:
				print "hull1 FAIL"
				
			if hull2_md5_nom == hull2_md5:
				print "hull2 PASS"
			else:
				print "hull2 FAIL"

			if axis1_md5_nom == axis1_md5:
				print "axis1 PASS"
			else:
				print "axis1 FAIL"
				
			if axis2_md5_nom == axis2_md5:
				print "axis2 PASS"
			else:
				print "axis2 FAIL"
			
			return	


	def testSpline(self):
			PIXELSIZE = 0.05

			currNode = LocalNode(self.probe, self.contacts, 10, 19, PIXELSIZE)
			currNode.readFromFile2("result_2012_08_27", 10)			
			hull1, axis1 = computeHullAxis(10, currNode)

			spline1 = SplineFit(axis1, smooth=0.1)
			spline2 = SplineFit(axis1, smooth=0.1)
			
			points1 = spline1.getUniformSamples()
			points2 = spline2.getUniformSamples()

			m = hashlib.md5()
			m.update(repr(points1))
			points1_md5 = int(m.digest().encode('hex'),16)
			m = hashlib.md5()
			m.update(repr(points2))
			points2_md5 = int(m.digest().encode('hex'),16)
								
			points1_md5_nom = 264552187170657706332967220739035361890
			points2_md5_nom = 261841834995931355572820968043753689877			

			if points1_md5_nom == points1_md5:
				print "points1 PASS"
			else:
				print "points1 FAIL"
				
			if points2_md5_nom == points2_md5:
				print "points2 PASS"
			else:
				print "points2 FAIL"

			return
		
	def testToro(self):
		
		num_poses = 10
		
		poseGraph = PoseGraph(self.probe, self.contacts)					
		poseGraph.restoreState("result_2012_08_27", num_poses)

		PIXELSIZE = 0.05
		for i in range(0, num_poses):
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile2("result_2012_08_27", i)		
			poseGraph.nodeHash[i] = currNode		
		
		
		poseGraph.mergePriorityConstraints()

		m = hashlib.md5()
		m.update(repr(poseGraph.edgePriorityHash))
		edges1_md5 = int(m.digest().encode('hex'),16)
		
		poseGraph = PoseGraph(self.probe, self.contacts)					
		poseGraph.restoreState("result_2012_08_27", num_poses)

		PIXELSIZE = 0.05
		for i in range(0, num_poses):
			currNode = LocalNode(self.probe, self.contacts, i, 19, PIXELSIZE)
			currNode.readFromFile2("result_2012_08_27", i)		
			poseGraph.nodeHash[i] = currNode		
		
		
		poseGraph.mergePriorityConstraints()
		
		m = hashlib.md5()
		m.update(repr(poseGraph.edgePriorityHash))
		edges2_md5 = int(m.digest().encode('hex'),16)
		
		edges1_md5_nom = 146731739815347910490180057980616635724
		edges2_md5_nom = 146731739815347910490180057980616635724

		if edges1_md5_nom == edges1_md5:
			print "edges1 PASS"
		else:
			print "edges1 FAIL"
			
		if edges2_md5_nom == edges2_md5:
			print "edges2 PASS"
		else:
			print "edges2 FAIL"


		#print edges1_md5
		#print edges2_md5
		
		#print repr(poseGraph.edgePriorityHash)
	

	def testPCA(self):


		posture = [[0.32402250170707703, -1.9376647472381592, 2.5240442752838135], [0.20172742009162903, -1.8508090972900391, 2.5791583061218262], [0.074833497405052185, -1.7708220481872559, 2.5749032497406006], [-0.051718935370445251, -1.6902958154678345, 1.0533870458602905], [0.022475600242614746, -1.5599303245544434, 0.95854592323303223], [0.10868218541145325, -1.4371768236160278, 1.8182663917541504], [0.071939393877983093, -1.2917464971542358, -3.0718615055084229], [-0.077696070075035095, -1.3021976947784424, 2.8612296581268311], [-0.22183933854103088, -1.2606920003890991, 1.2641770839691162], [-0.17656373977661133, -1.1176880598068237, 0.83170008659362793], [-0.075520709156990051, -1.0068264007568359, 2.0582084655761719], [-0.14577186107635498, -0.87429428100585938, 2.5009474754333496], [-0.26602840423583984, -0.78463733196258545, 2.0668435096740723], [-0.33742135763168335, -0.65271675586700439, 0.47807663679122925], [-0.2042391449213028, -0.58370590209960938, 0.21110701560974121], [-0.057569209486246109, -0.55227452516555786, 1.2009902000427246], [-0.0033539943397045135, -0.41241493821144104, 1.7855242490768433], [-0.035316232591867447, -0.26585978269577026, 1.9033694267272949], [-0.084287673234939575, -0.12407895922660828, 0.97409385442733765], [0.0, 0.0, 0.0], [0.15000000596046448, 0.0, 0.076145566999912262], [0.29956534504890442, 0.011410800740122795, 0.86578977108001709], [0.39677116274833679, 0.12565189599990845, 1.946380615234375], [0.34174874424934387, 0.26519590616226196, 1.7029930353164673], [0.32197695970535278, 0.41388711333274841, 1.3524518013000488], [0.35446903109550476, 0.5603257417678833, 0.19331979751586914], [0.5016748309135437, 0.58914345502853394, -0.0050837248563766479], [0.65167289972305298, 0.58838087320327759, 0.57785046100616455], [0.77731871604919434, 0.67031455039978027, 1.8788957595825195], [0.73183149099349976, 0.8132513165473938, 2.0654768943786621], [0.66061890125274658, 0.94526934623718262, 1.1665383577346802], [0.71961939334869385, 1.0831785202026367, -0.035087704658508301], [0.86952710151672363, 1.0779165029525757, -0.065629377961158752], [1.0192041397094727, 1.0680791139602661, 0.95869624614715576], [1.1053922176361084, 1.1908456087112427, 1.8644802570343018], [1.0619701147079468, 1.3344231843948364, 1.7305004596710205], [1.0381162166595459, 1.4825143814086914, 1.3925573825836182], [1.0647107362747192, 1.6301380395889282, 0.31836950778961182], [1.2071727514266968, 1.6770907640457153, -1.2506951093673706]]

		x_list = []
		y_list = []
	
		for p in posture:
			x_list.append(p[0])
			y_list.append(p[1])
	
		cov_a = scipy.cov(x_list,y_list)
	
		loadings = []
	
		" NOTE:  seems to create opposing colinear vectors if data is colinear, not orthogonal vectors "
	
		scores, loadings, E = pca_module.nipals_mat(cov_a, 2, 0.000001, False)	

		m = hashlib.md5()
		m.update(repr(loadings))
		loadings1_md5 = int(m.digest().encode('hex'),16)

		scores, loadings, E = pca_module.nipals_mat(cov_a, 2, 0.000001, False)	

		m = hashlib.md5()
		m.update(repr(loadings))
		loadings2_md5 = int(m.digest().encode('hex'),16)

		loadings1_md5_nom = 4604468617513472537582329011925291178
		loadings2_md5_nom = 4604468617513472537582329011925291178

		if loadings1_md5_nom == loadings1_md5:
			print "loadings1 PASS"
		else:
			print "loadings1 FAIL"
			
		if loadings2_md5_nom == loadings2_md5:
			print "loadings2 PASS"
		else:
			print "loadings2 FAIL"

	def testSCSGP(self):
		totalHypotheses = 3
		A = matrix(totalHypotheses*totalHypotheses*[0.0], dtype=float)
		A.resize(totalHypotheses, totalHypotheses)

		errPac = 1.0
		results = []
		for i in range(totalHypotheses):
			for j in range(i+1, totalHypotheses):
				results.append([errPac, i, j])
				errPac += 1.0
				
		for result in results:
			i = result[1]
			j = result[2]
			A[i,j] = result[0]
			A[j,i] = result[0]
	

		e, lmbda = scsgp.dominantEigenvectors2(A)

		m = hashlib.md5()
		m.update(repr(e))
		e1_md5 = int(m.digest().encode('hex'),16)

		m = hashlib.md5()
		m.update(repr(lmbda))
		lmbda1_md5 = int(m.digest().encode('hex'),16)

		e1_md5_nom = 74415382024037744883544301949964047172
		lmbda1_md5_nom = 91214331286978604462059783908825041166

		if e1_md5_nom == e1_md5:
			print "e1 PASS"
		else:
			print "e1 FAIL"
			
		if lmbda1_md5_nom == lmbda1_md5:
			print "lamda1 PASS"
		else:
			print "lamda1 FAIL"

