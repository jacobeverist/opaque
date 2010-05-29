from SnakeControl import SnakeControl
from copy import *
from math import *
import time

class TestTransform(SnakeControl):

	"""
	Unit tests
	
	1. never returns
	2. reset position commands the probe to position
	3. None values in the joints do not affect the probe commands
	4. Changing the joint real values will be compensated by behavior
	
	Setup the probe
	Add the controller
	Drive the joints
	Step and get results
	
	"""

	def __init__(self, probe, drawThings):
		SnakeControl.__init__(self)
	
		self.drawThings = drawThings
		self.probe = probe
		
		self.robotParam = self.probe.robotParam
		
		self.setTimerAliasing(1)
		
		self.nextState = False
		
		self.stateA = 0		
		self.globalTimer = 1

	def frameStarted(self):

		# alias the control sequence by some time step
		SnakeControl.frameStarted(self)
		if not self.isStep():
			return

		" get the probe state "
		probeState = self.probe.getProbeState()

		" initialize probe "
		if self.stateA == 0:	

			for i in range(self.robotParam['numJoints']):

				modVal = i % 4
				
				if modVal == 0:
					self.probe.setServo(i, 0.0)
				elif modVal == 1:
					self.probe.setServo(i, 90.0)
				elif modVal == 2:
					self.probe.setServo(i, 0.0)
				elif modVal == 3:
					self.probe.setServo(i, -90.0)

			self.probe.getProbeState()
			self.stateA = 1		
			self.t1 = time.time()
				
			print "step 0 complete"

		elif self.stateA == 1:

			isDone = False
			
			#if self.globalTimer % 10 == 0:
			#	rootPose = [0.0,0.0,0.0]
				
			#	poses = []
			#	for i in range(self.robotParam['numJoints']):
			#		poses.append(self.probe.getJointFromJoint(rootPose, 19, i))					
			#	self.drawThings.plotRobotConfiguration(poses)

			rootPose = [1.0,1.0,pi/4]		
			for j in range(39):		
				for i in range(self.robotParam['numJoints']):
					self.probe.getJointFromJoint(rootPose, j, i)				
					#self.probe.getJointWRTJointPose(rootPose, j, i)					

			#print self.probe.jointTransforms[18][19], self.probe.jointTransforms[19][18]
			#print self.probe.jointTransforms[0][19], self.probe.jointTransforms[19][0]
			#print self.probe.jointTransforms[0][20], self.probe.jointTransforms[20][0]
			#print self.probe.jointTransforms[1][20], self.probe.jointTransforms[20][1]
			#print self.probe.jointTransforms[0][38], self.probe.jointTransforms[38][0]

			#for i in range(self.robotParam['numJoints']):

			"""
			poses = []
			for i in range(19,39):
				poses.append(self.probe.getJointFromJoint(rootPose, 19, i))					
			self.drawThings.plotRobotConfiguration(poses)

			poses = []
			for i in range(19,39):
				poses.append(self.probe.getJointWRTJointPose(rootPose, 19, i))					
			self.drawThings.plotRobotConfiguration(poses)
			"""
			"""
			poses1 = []
			for i in range(self.robotParam['numJoints']):
				poses1.append(self.probe.getJointFromJoint(rootPose, 19, i))					
			self.drawThings.plotRobotConfiguration(poses1)

			poses2 = []
			for i in range(self.robotParam['numJoints']):
				poses2.append(self.probe.getJointWRTJointPose(rootPose, 19, i))					
			self.drawThings.plotRobotConfiguration(poses2)
			
			for i in range(len(poses1)):
				for j in range(3):
					if fabs(poses1[i][j] - poses2[i][j]) > 1e-10:
						print poses1[i][j], "!=", poses2[i][j]
			"""
			if self.globalTimer == 1000:
				#if self.globalTimer == 3:

				t2 = time.time()

				result = (t2-self.t1)*1000.0
				print result
				exit()
							
			if isDone:
				print "step 1 complete"
		
				self.stateA = 2


				self.testSuccess()

			
								

				
