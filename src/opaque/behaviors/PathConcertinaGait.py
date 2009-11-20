import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from copy import *
from common import *
from Behavior import *
from BlindConcertinaGait import BlindConcertinaGait
from GlobalCurveFit import GlobalCurveFit

class PathConcertinaGait(Behavior):

	#def __init__(self, probe, direction = False, curve = 0):
	def __init__(self, probe, direction = False, path = []):
		Behavior.__init__(self, probe)

		self.blindConcertina = BlindConcertinaGait(probe, direction)
		
		self.globalCurve = 0
		
		self.isInit = False
		self.startCount = 0
		
		#self.setCurve(curve)
		self.setPath(path)

	def setDirection(self, direction):
		self.direction = direction
		self.blindConcertina.setDirection(direction)
		
	def setPath(self, path):
		self.path = deepcopy(path)
		self.computeCurve()
		
	def reverseDirection(self):
		self.path.reverse()
		self.computeCurve()
		
	def computeCurve(self):
		self.curve = VoronoiFit(self.path)
		
		if self.globalCurve == 0:
			self.globalCurve = GlobalCurveFit(self.probe, self.curve)
		else:
			self.globalCurve.setCurve(self.curve)
		
		self.setDirection(self.globalCurve.getPathDirection())
			
		self.globalCurve.draw()

	def step(self):
		
		#if not self.isInit:
		#	self.startCount += 1
		#	isDone = self.blindConcertina.step()
		#	
		#	if self.startCount >= 5:
		#		self.isInit = True

		# perform concertina gait
		isDone = self.blindConcertina.step()
		#isDone = False
		
		# mask is used to determine the segments that need to be fitted to the curve
		mask = self.blindConcertina.getMask()
		#mask = [1.0 for i in range(39)]
		
		#for i in range(20,39):
		#	mask[i] = 0.0
		
		# compute the bounds for behavior for global curve fitting
		startNode, endNode = self.computeBounds(mask)
		self.globalCurve.setBoundaries(startNode, endNode)
		
		# execute the global curve fitting
		self.globalCurve.step()
		
		# collect joint settings for both behaviors
		joints1 = self.blindConcertina.getJoints()
		joints2 = self.globalCurve.getJoints()
		
		# Set the joints, with joints2 as priority over joints1
		self.mergeJoints([joints2, joints1])
		
		return isDone


	def computeBounds(self, mask):
		
		# find indices of the upper and lower bound of high mask threshold
		thresh = 0.8
		lowIndex = len(mask)-1
		highIndex = 0
		for i in range(len(mask)):
			if i < lowIndex and mask[i] < thresh:
				lowIndex = i

			if i > highIndex and mask[i] < thresh:
				highIndex = i

		startNode = -1
		endNode = -1
		
		if lowIndex != len(mask)-1 or highIndex != 0:
			startNode = lowIndex
			endNode = highIndex

		return startNode, endNode
	
	

	
	
	
	
	
	