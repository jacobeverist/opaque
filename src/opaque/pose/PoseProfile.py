#!/usr/bin/python

from math import cos, sin
from copy import copy


class PoseProfile:

	def __init__ (self, contacts, rootNode, inSim = True):

		self.rootNode = rootNode
		self.contacts = contacts
		self.allNodes = []
		
		if inSim:
			for i in range(-1,self.contacts.numJoints+1):
				self.allNodes.append(self.contacts.getAveragePose(i))		
			
			self.rootPose = self.allNodes[self.rootNode+1]

	def saveToFile(self, filename):

		val = repr(self.allNodes)
		
		f = open(filename, 'w')
		f.write(val)
		f.close()
	
	def readFromFile(self, filename):
		f = open(filename, 'r')
		val = eval(f.read())
		f.close()
		
		self.allNodes = val
		
	def printPosePoints(self, clr):

		xP = []
		yP = []
		
		for p in self.allNodes:
			xP.append(p[0])
			yP.append(p[1])

		" plot all points by (x,y) "

		" flipped to align with Ogre x-z coord system "
		#pylab.scatter(xP, yP, color=clr, faceted=False, linewidth=1)

	def performOffset(self, offset, refPoint):

		origin = refPoint
		ang = offset[2]
		
		for p in self.allNodes:
			p[0] += offset[0]
			p[1] += offset[1]

			xdiff = p[0] - origin[0]
			ydiff = p[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			p[0] = origin[0] + xnew
			p[1] = origin[1] + ynew
			p[2] += ang

		self.separatePoints()
	
	def getRootPose(self):
		return copy(self.rootPose)

	def setRootPose(self, newRoot):
		self.rootPose = copy(newRoot)

