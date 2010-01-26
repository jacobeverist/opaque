import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
    sys.path.append(dir)

from common import *
from Behavior import *

class Transition(Behavior):
    
    def __init__(self, probe):
        Behavior.__init__(self, probe)
    
        self.resetTime()
        self.initJoints = [None for i in range(self.probe.numSegs-1)]
        self.targetJoints = [None for i in range(self.probe.numSegs-1)]
    
    def setInit(self, joints):
        self.initJoints = copy(joints)
    
    def setTarget(self, joints):
        self.targetJoints = copy(joints)
    
    def resetTime(self, newTime = 500):
        self.transitionTime = 0
        
        if newTime == 0:
        	self.timeLength = 1
        else:
			self.timeLength = newTime
    
    def step(self):
        Behavior.step(self)
            
        if not self.isStep():
            return False
        
        self.resetJoints()

        self.transitionTime += 1
        
        if self.transitionTime > self.timeLength:
            self.transitionTime = self.timeLength

        diffs = []
        
        for i in range(self.probe.numSegs-1):
            
			if self.initJoints[i] == None:
				self.setJoint(i, None)
            	
			else:
				diff = self.targetJoints[i] - self.initJoints[i]
				
				#diffs.append(diff)
				
				val = 1.0 - float(self.timeLength - self.transitionTime) / self.timeLength
				
				newAngle = self.initJoints[i] + diff*val
				
				self.setJoint(i, newAngle)


        if self.transitionTime >= self.timeLength:
            return True
        else:
            return False

