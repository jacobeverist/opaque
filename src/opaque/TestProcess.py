import random
random.seed(0)

from SnakeControl import SnakeControl
from copy import *
from math import *
import time
import math
import os

from behaviors.AnchorTransition import AnchorTransition

from pose.AverageContacts import AverageContacts
from maps.MapLoader import MapLoader

import numpy
import sys

class TestProcess(SnakeControl):

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe, drawThings, args):
		SnakeControl.__init__(self)
	
		self.drawThings = drawThings
		self.probe = probe

		self.args = args

		self.maxNumPoses = self.args.maxNumPoses
		
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
			self.mapGraph = MapLoader(self.probe, self.contacts, args = self.args)

			targetMapFile = "../mapLibrary/cross_junction.py"
			#targetMapFile = "../mapLibrary/Y_junction.py"
			#targetMapFile = "mapFile0002.txt"
			f = open(targetMapFile, 'r')
			str_f = f.read()
			f.close()

			exec(str_f)

			self.mapGraph.loadWalls(walls)

			#self.mapGraph = MapLoader(self.probe, self.contacts)
			#self.mapGraph.loadWalls(walls)
			#self.mapGraph.loadSeries("../results/result_2013_08_16_Y_junction_complete", 290)
			#retvalue = os.system("mkdir data2")
			#retvalue = os.system("mv *.png data2")

			#self.mapGraph = MapLoader(self.probe, self.contacts)
			#self.mapGraph.loadWalls(walls)
			#self.mapGraph.loadSeries("../results/result_2013_08_18_cross", 236)
			#retvalue = os.system("mkdir data4")
			#retvalue = os.system("mv *.png data4")

			#self.mapGraph = MapLoader(self.probe, self.contacts)
			#self.mapGraph.loadWalls(walls)

			#self.mapGraph.loadSeries("../testDir5/results/2014_09_24/", 0, 206)
			#self.mapGraph.restorePickle("../testDir5", 47)
			#self.mapGraph.loadSeries("../testDir5", 48, 52)
			#self.mapGraph.restorePickle("../testDir5/results/2014_09_28", 49)
			#self.mapGraph.loadSeries("../testDir5/result/2014_09_28", 50, 52)
			#self.mapGraph.restorePickle("../testDir5/", 83)
			#self.mapGraph.loadSeries("../testDir5/", 84, 100)
			#self.mapGraph.loadData("../testDir5/results/2014_09_29_cross_junction", 250)
			#self.mapGraph.loadData("../testDir5/", 250)
			#self.mapGraph.loadSeries("../testDir5/", 0, 250)
			#self.mapGraph.restorePickle("../testDir5", 143)
			#self.mapGraph.loadSeries("../testDir5", 144, 190)
			#self.mapGraph.restorePickle("../testDir6/results/2014_09_25_divergence_and_branching", 103)
			#self.mapGraph.loadSeries("../testDir6/results/2014_09_25_divergence_and_branching", 104, 108)

			#self.mapGraph.restorePickle("../testDir6/", 17)
			#self.mapGraph.loadSeries("../testDir5/", 18, 250)
			#self.mapGraph.restorePickle("../testDir4/results_L_right_junction.py/", 9)
			#self.mapGraph.loadSeries("../testDir4/results_L_right_junction.py/", 10, 14)

			#self.mapGraph.restorePickle("../../doc/results/tip_results_2014_11_22_23:00/cross_0.0/testDir", 15)
			#self.mapGraph.loadSeries("../../doc/results/tip_results_2014_11_22_23:00/cross_0.0/testDir", 16, 1000)
			#self.mapGraph.restorePickle("../../doc/results/tip_results_2014_11_23/cross_0.0/testDir", 21)
			#self.mapGraph.loadSeries("../../doc/results/tip_results_2014_11_23/cross_0.0/testDir", 22, 1000)
			#self.mapGraph.restorePickle("../../doc/results/eval_2014_11_29c/cross_0.4/testDir", 17)

			#self.mapGraph.restorePickle("../../doc/results/eval_2014_12_01/cross_0.2/testDir", 47)
			#self.mapGraph.restorePickle("../../doc/results/eval_2014_12_01/cross_0.2/testDir", 19)
			#for pID, mapHyp in self.mapGraph.mapAlgorithm.mapHyps.iteritems():
			#	mapHyp.generatePaths()
			#	mapHyp.drawPoseParticles()

			#self.mapGraph.loadSeries("../../doc/results/eval_2014_12_01/cross_0.2/testDir", 20, 40)
			#self.mapGraph.loadSeries("../../doc/results/eval_2014_11_29c/cross_0.4/testDir", 18, 1000)
			#for pID, mapHyp in self.mapGraph.mapAlgorithm.mapHyps.iteritems():
			#	mapHyp.generatePaths()
			#	mapHyp.drawPoseParticles()
			#self.mapGraph.restorePickle("./notrim_eval_2014_12_06b/", 13)
			#self.mapGraph.loadSeries("../../doc/results/eval_2014_12_01/cross_0.0/testDir", 14, 1000)

			#self.mapGraph.restorePickle('../../doc/results/eval_2014_12_08/results_0.2/testDir/results_T_89_junction.py', 9)
			#for pID, mapHyp in self.mapGraph.mapAlgorithm.mapHyps.iteritems():
			#	mapHyp.computeEval()
			#	mapHyp.generatePaths()
			#	mapHyp.drawPoseParticles()
			 
			#self.mapGraph.loadSeries('../../doc/results/eval_2014_12_08/results_0.2/testDir/results_T_89_junction.py', 10, 20)

			#self.mapGraph.restorePickle('./eval_longrun_2014_12_09', 3)
			#self.mapGraph.loadSeries('./eval_longrun_2014_12_09', 4, 6)
			#self.mapGraph.loadSeries('./eval_longrun_2014_12_09', 4, 1000)
			#self.mapGraph.loadSeries('./eval_longrun_2014_12_09', 0, 1000)
			#self.mapGraph.restorePickle('./eval_maxbranch_2014_12_14', 17)
			#self.mapGraph.loadSeries('./eval_maxbranch_2014_12_14', 18, 1000)
			#self.mapGraph.restorePickle('./eval_nondiverge_2014_12_14', 17)
			#self.mapGraph.loadSeries('./eval_nondiverge_2014_12_14', 18, 1000)
			#self.mapGraph.restorePickle('./eval_longrun_2014_12_09', 3)
			#self.mapGraph.loadSeries('./eval_longrun_2014_12_09', 4, 1000)

			#self.mapGraph.restorePickle('./eval_maxbranch_2014_12_14', 85)
			#self.mapGraph.loadSeries('./eval_maxbranch_2014_12_14', 86, 1000)
			#self.mapGraph.restorePickle('./eval_maxbranch_2014_12_14', 83)
			#self.mapGraph.loadSeries('./eval_maxbranch_2014_12_14', 84, 1000)

			#self.mapGraph.restorePickle('./eval_nondiverge_2014_12_14b', 105)
			#self.mapGraph.loadSeries('./eval_nondiverge_2014_12_14b', 106, 1000)

			#self.mapGraph.loadSeries('../testDir6/eval_longrun_2014_12_09', 0, 1000)

			#self.mapGraph.restorePickle('./eval_control_2014_12_15', 19)
			#self.mapGraph.loadSeries('./eval_control_2014_12_15', 20, 1000)

			#self.mapGraph.restorePickle('./eval_control_2014_12_16', 5)
			#self.mapGraph.loadSeries('./eval_control_2014_12_16', 6, 1000)

			#self.mapGraph.restorePickle('./eval_control_2014_12_16', 3)
			#self.mapGraph.loadSeries('./eval_control_2014_12_16', 4, 1000)
			#self.mapGraph.restorePickle('./eval_merge_2014_12_17', 17)
			#self.mapGraph.loadSeries('./eval_merge_2014_12_17', 18, 1000)


			#self.mapGraph.restorePickle('./eval_dist_2014_12_21', 19)
			#self.mapGraph.loadSeries('./eval_dist_2014_12_21', 20, 2200)
			#self.mapGraph.loadSeries('./eval_dist_2014_12_21', 0, 5400)
			#self.mapGraph.loadSeries('./eval_dist_2014_12_21', 0, 1000)
			#self.mapGraph.restorePickle('./eval_dist_2014_12_21', 147)
			#self.mapGraph.loadSeries('./eval_dist_2014_12_21', 148, 1000)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2014_12_19/result2/results_cross_junction.py', 0, 1000)
			#self.mapGraph.loadData('/home/everist/opaque/doc/results/eval_2014_12_19/result2/results_cross_junction.py', 1000)
			#self.mapGraph.restorePickle('.', 7)
			#self.mapGraph.loadSeries('.', 8,100)

			#self.mapGraph.restorePickle('./eval_dist_2014_12_21', 51)
			#self.mapGraph.loadSeries('./eval_dist_2014_12_21', 52, 5400)


			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/eval_2014_12_24/results_0.0/testDir/results_Y_junction.py', 45)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2014_12_24/results_0.0/testDir/results_Y_junction.py', 46, 10000)

			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/eval_2014_12_24/results_0.2/testDir/results_T_89_junction.py', 77)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2014_12_24/results_0.2/testDir/results_T_89_junction.py', 78, 10000)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2014_12_24/results_0.2/testDir/results_T_89_junction.py', 0, 10000)

			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/T_side_2014_12_29/testDir', 27)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/T_side_2014_12_29/testDir', 28, 30)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/T_side_2014_12_29/testDir', 0, 3000)

			#self.mapGraph.restorePickle('./results_Y_junction.py/', 153)
			#self.mapGraph.loadSeries('./results_Y_junction.py/', 154, 10000)


			""" displace/localize failure modes """
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/results_0.2/testDir', 11)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/results_0.2/testDir', 12, 14)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/results_0.2/testDir/results_T_side_junction.py', 47)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/results_0.2/testDir/results_T_side_junction.py', 48, 50)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/results_0.2/testDir/results_T_side_junction.py', 0, 50000)
			#self.mapGraph.restorePickle('./results_merge_error_2015_01_01', 37)
			#self.mapGraph.loadSeries('./results_merge_error_2015_01_01', 38, 40)
			#self.mapGraph.restorePickle('./results_merge_error_2015_01_01', 17)
			#self.mapGraph.loadSeries('./results_merge_error_2015_01_01', 18, 20)
			#self.mapGraph.restorePickle('./results_merge_error_2015_01_01b', 47)
			#self.mapGraph.loadSeries('./results_merge_error_2015_01_01b', 48, 5000)
			
			#self.mapGraph.restorePickle('./results_merge_error_2015_01_02', 19)
			#self.mapGraph.loadSeries('./results_merge_error_2015_01_02', 20, 5000)
			#self.mapGraph.restorePickle('./eval_merge_error_2015_01_03c', 37)
			#self.mapGraph.loadSeries('./eval_merge_error_2015_01_03c', 38, 5000)

			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/results_0.2/testDir/results_cross_junction.py', 0, 10000)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/results_0.2/testDir/results_Y_junction.py', 0, 10000)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2015_01_04/testDir/results_Y_junction.py', 0, 10000)
			#self.mapGraph.loadSeries('./success_T_side_2015_01_05', 0, 10000)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2015_01_05/testDir/results_cross_junction.py', 0, 10000)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/eval_2015_01_06/testDir/results_cross_junction.py', 89)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2015_01_06/testDir/results_cross_junction.py', 90, 10000)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/merge_grandchild_2015_01_07', 55)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/merge_grandchild_2015_01_07', 56, 10000)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/T_side_0.4_2015_01_08/testDir', 67)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/T_side_0.4_2015_01_08/testDir', 68, 1000)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/T_side_0.4_2015_01_08/testDir', 0, 1000)
			#self.mapGraph.restorePickle('./eval_T_side_0.4_2015_01_08', 33)
			#self.mapGraph.loadSeries('./eval_T_side_0.4_2015_01_08', 34, 36)
			#self.mapGraph.loadSeries('./results_T_89_2015_01_06', 0, 10000)
			#self.mapGraph.restorePickle('../testDir4', 57)
			#self.mapGraph.loadSeries('../testDir4', 0, 6000)
			#self.mapGraph.restorePickle('eval_2015_01_18', 173)
			#self.mapGraph.loadSeries('eval_2015_01_18', 174,176)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/eval_2015_01_13/results_0.0/testDir/results_Y_junction.py', 41)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2015_01_13/results_0.0/testDir/results_Y_junction.py', 42,44)
			#self.mapGraph.restorePickle('eval_2015_01_13', 227)
			#self.mapGraph.loadSeries('eval_2015_01_13', 228,400)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/eval_2015_01_13/results_0.2/testDir/results_T_side_junction.py', 125)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2015_01_13/results_0.2/testDir/results_T_side_junction.py', 126,128)
			#self.mapGraph.restorePickle('/home/everist/opaque/doc/results/eval_2015_01_13/results_0.2/testDir/results_Y_junction.py', 41)
			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2015_01_13/results_0.2/testDir/results_Y_junction.py', 42,43)
			self.mapGraph.restorePickle('./dataTempResults/', 231)
			self.mapGraph.loadSeries('./dataTempResults/', 232,1000)

			#self.mapGraph.restorePickle('./eval_mergecase_2015_01_06/', 41)
			#self.mapGraph.loadSeries('./eval_mergecase_2015_01_06/', 0, 1000)
			#self.mapGraph.loadSeries('./eval_mergecase_2015_01_06/', 0, 1000)

			#self.mapGraph.loadSeries('/home/everist/opaque/doc/results/eval_2014_12_24/results_0.0/testDir/results_Y_junction.py', 46, 10000)

			#self.mapGraph.restorePickle('./eval_motion_2014_12_29/', 65)
			#self.mapGraph.loadSeries('./eval_motion_2014_12_29/', 66, 10000)

			#self.mapGraph.restorePickle('./eval_internal_2014_12_27', 95)
			#self.mapGraph.loadSeries('./eval_internal_2014_12_27', 96, 5400)

			#self.mapGraph.loadSeries('./eval_dist_2014_12_21', 0, 1000)
			exit()

			#self.mapGraph.restorePickle("./results_cross_junction.py/", 17)

			
			for pID, mapHyp in self.mapGraph.mapAlgorithm.mapHyps.iteritems():
				#if 1 in mapHyp.pathClasses.keys():
				#	mapHyp.mergePath(1)
				#	mapHyp.generatePaths()
				#	mapHyp.drawPoseParticles()
				print pID, "before mapHyp.branchDiverges:", mapHyp.branchDiverges
				mapHyp.branchDiverges[1] = False

				print pID, "after mapHyp.branchDiverges:", mapHyp.branchDiverges

				print "caseA"

				while False in mapHyp.branchDiverges.values():

					print pID, "branchDiverges:", mapHyp.branchDiverges

					pathIDs = mapHyp.branchDiverges.keys()
					
					print "caseB"

					
					for pathID, isDiverge in mapHyp.branchDiverges.iteritems():
						print "caseC", pathID
						if not isDiverge:
							print "caseD"
							mapHyp.mergePath(pathID)
							mapHyp.generatePaths()
							mapHyp.drawPoseParticles()
							print "caseE:", mapHyp.branchDiverges
							
							break

		
			print pID, "mapHyp.branchDiverges:", mapHyp.branchDiverges

			#self.mapGraph.loadSeries("./results_cross_junction.py/", 18, 22)


			#self.mapGraph.restorePickle("../../doc/results/tip_results_2014_11_22/cross_0.0/testDir", 23)
			#self.mapGraph.loadSeries("../../doc/results/tip_results_2014_11_22/cross_0.0/testDir", 24, 1000)
			#self.mapGraph.restorePickle("../../doc/results/Y_0.0/testDir", 33)
			#self.mapGraph.loadSeries("../../doc/results/Y_0.0/testDir", 34, 1000)
			#self.mapGraph.restorePickle("../testDir4/cross_0.2_2014_11_21", 27)
			#self.mapGraph.loadSeries("../testDir4/cross_0.2_2014_11_21", 28, 1000)
			#self.mapGraph.loadSeries("../testDir4/cross_0.2_2014_11_21", 0, 1000)
			#self.mapGraph.loadSeries("../testDir4/", 0, 1000)
			#self.mapGraph.restorePickle("../testDir6/", 33)
			#self.mapGraph.loadSeries("../testDir5/", 34, 250)
			#self.mapGraph.loadData("/home/everist/opaque/doc/results/tip_results_2014_11_23/60_left_0.0", 250)
			#self.mapGraph.restorePickle("/home/everist/opaque/doc/results/tip_results_2014_11_23/T_side_0.2/testDir", 29)
			#self.mapGraph.loadSeries("/home/everist/opaque/doc/results/tip_results_2014_11_23/T_side_0.2/testDir", 30, 38)
			#self.mapGraph.restorePickle("/home/everist/opaque/doc/results/tip_results_2014_11_23/T_bottom_0.4/testDir", 29)
			#self.mapGraph.loadSeries("/home/everist/opaque/doc/results/tip_results_2014_11_23/T_bottom_0.4/testDir", 0, 34)
			#self.mapGraph.loadSeries("/home/everist/opaque/doc/results/tip_results_2014_11_23/cross_0.0/testDir", 0, 34)
			#self.mapGraph.restorePickle(".", 25)
			#self.mapGraph.mapAlgorithm
			#for pID, mapHyp in self.mapGraph.mapAlgorithm.mapHyps.iteritems():
			#	mapHyp.generatePaths()

			#self.mapGraph.loadSeries("/home/everist/opaque/doc/results/tip_results_2014_11_23/cross_0.0/testDir", 26, 34)
			

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
