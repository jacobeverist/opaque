import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

from BackBoneIK import *
import ogre.renderer.OGRE as ogre
from math import *

class CosineBackBone(BackBoneIK):

	def __init__(self, probe, freq, amp, lead_len, trail_len, rootNode = 19, backwards = True):
		BackBoneIK.__init__(self, probe, rootNode, backwards)
		
		self.newCurve(freq, amp, lead_len, trail_len)

	def newCurve(self, freq, amp, lead_len, trail_len):

		# leading length
		self.lead_len = lead_len

		# trailing length
		self.trail_len = trail_len
	
		# frequency
		self.freq = freq

		# amplitude
		self.amp = amp

		# 1. constructor
		# 2. getU
		# 3. getUVector
		# 4. findClosestPoint

		# 1. set origin to 0,0
		# 2. rotate so origin vector is [1,0]
		self.origin = [0.0,0.0]
		self.oVec = [1.0,0.0]
		self.curve = CosineFit(self.freq, self.amp, self.lead_len, self.trail_len, self.origin, self.oVec)

	def invert(self):
		self.amp *= -1
		self.curve = CosineFit(self.freq, self.amp, self.lead_len, self.trail_len, self.origin, self.oVec)

	def setLead(self, len):
		self.lead_len = len
		self.curve = CosineFit(self.freq, self.amp, self.lead_len, self.trail_len, self.origin, self.oVec)
	
	def setTrail(self, len):
		self.trail_len = len
		self.curve = CosineFit(self.freq, self.amp, self.lead_len, self.trail_len, self.origin, self.oVec)

	def setAmp(self, amp):
		self.amp = amp
		self.curve = CosineFit(self.freq, self.amp, self.lead_len, self.trail_len, self.origin, self.oVec)
		
	def setFreq(self, freq):
		self.freq = freq
		self.curve = CosineFit(self.freq, self.amp, self.lead_len, self.trail_len, self.origin, self.oVec)

	def draw(self):

		self.curve.draw()

		# remove all children
		self.parentNode.removeAllChildren()

		# deference the child nodes now
		for child in self.childNodes:
			self.probe._mgr.destroySceneNode(child)

		self.childNodes = []

		for child in self.childEntities:
			self.probe._mgr.destroyEntity(child)

		self.childEntities = []

		points = self.curve.getPoints()

		# draw the points in the simulation
		for i in range(len(points)):

			pnt = points[i]

			childNode = self.parentNode.createChildSceneNode("curvePoint" + str(i))
			self.childNodes.append(childNode)

			currEntity = self.probe._mgr.createEntity("curveEnt"+str(i), "Cube.mesh")
			currEntity.setCastShadows(False)
			currEntity.setMaterialName("Green")
			self.childEntities.append(currEntity)

			position = ogre.Vector3(pnt[0],0.0,pnt[1])
			#position = ogre.Vector3(0.0,0.0,0.0)
			childNode.setPosition(position)

			size = ogre.Vector3(0.05,0.05,0.05)
			childNode.setScale(size)

			childNode.attachObject(currEntity)

