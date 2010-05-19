import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *

from Anchor import Anchor
from LocalCurveFit import LocalCurveFit
from Transition import Transition

class SpaceSweep(Behavior):
	
	def __init__(self, probe, direction = False):
		Behavior.__init__(self, probe)

		self.dir = direction
		
		self.A_timer = 0
		self.curveSign = 1
		lead_len = 0.0
		trail_len = 0.0
		self.freq = 4*pi
		self.amp =0.23
		origin = [0.0,0.0]
		oVec = [1.0,0.0]
		curve2 = CosineFit(self.freq, self.amp, lead_len, trail_len, origin, oVec)

		# behaviors should be
		#self.anchor = Anchor(self.probe)
		self.anchor = Anchor(self.probe)
		self.localFit = LocalCurveFit(self.probe, True, False, curve2, 19)
		self.transition = Transition(self.probe)
		
		self.isTransitioning = False
		self.hasTransitioned = False
		
		#self.localFit = LocalCurveFit(self.probe, True, False, curve2, 10)
		#self.localFit = LocalCurveFit(self.probe, False, True, curve2, 19)

		self.setTimerAliasing(20)
		
		self.isInit = False
	
	def hasInitialized(self):
		return self.isInit
	
	def setDirection(self, isForward):
		if isForward:
			self.localFit.setSide(True)
			self.localFit.setDir(False)
			self.localFit.setRootNode(10)
		else:
			self.localFit.setSide(False)
			self.localFit.setDir(True)
			self.localFit.setRootNode(28)
		
	def step(self):
		Behavior.step(self)
		
		if self.isTransitioning:

			self.resetJoints()

			# first steps
			isDone = self.transition.step()
			
			resJoints = self.transition.getJoints()
			self.mergeJoints([resJoints])
		
			if isDone:
				self.isTransitioning = False
				
				if not self.isInit:
					self.isInit = True
		
			return False
		
		if not self.isStep():
			return False
		
		self.resetJoints()

		# update curve fitting parameters for this iteration
		lead_len = ((50.0 - self.A_timer)/50.0)*1.7
		lead_points = []

		trail_len = (self.A_timer/50.0)*1.7
		trail_points = []

		origin = [0.0,0.0]
		oVec = [1.0,0.0]
		curve2 = CosineFit(self.freq, self.curveSign * self.amp, lead_len, trail_len, origin, oVec)
		self.localFit.setCurve(curve2)
		
		#self.localFit.draw()

		self.A_timer += 1

		self.localFit.step()
		self.anchor.step()
		
		
		joints1 = self.localFit.getJoints()
		#joints2 = self.anchor.getJoints()
		
		#self.mergeJoints([joints1,joints2])
		self.mergeJoints([joints1])

		if not self.hasTransitioned:

			initState = []
			for i in range(self.numJoints):
				initState.append(180.0*stateJoints[i]/pi)

			targetState = self.getJoints()
			for i in range(len(targetState)):
				if targetState[i] == None:
					targetState[i] = 180.0*stateJoints[i]/pi

			errSum = 0.0
			for i in range(len(initState)):
				if initState[i] != None:
					errSum += fabs(initState[i]-targetState[i])
			#transTime = int(errSum*2)
			transTime = int(4*errSum)
			
			print "sweep transition time =", transTime

			# set target and initial poses
			self.transition.setInit(initState)
			self.transition.setTarget(targetState)
			#self.transition.resetTime()		
			self.transition.resetTime(transTime)		
			
			# first steps
			self.transition.step()
			
			resJoints = self.transition.getJoints()
			self.resetJoints()
			self.mergeJoints([resJoints])

			self.hasTransitioned = True
			self.isTransitioning = True				

		if self.curveSign == 1 and self.A_timer > 50:
			self.A_timer = 0
			self.A_state = 0
			self.curveSign = -1
			self.hasTransitioned = False

		elif self.curveSign == -1 and self.A_timer > 50:
			self.A_timer = 0
			self.A_state = 0
			self.curveSign = 1
			self.hasTransitioned = False
			
			self.isInit = False
			return True

		return False