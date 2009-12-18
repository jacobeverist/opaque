import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import traceback
import Image
import ImageDraw
from copy import *
from math import *
from numpy import arange
from maps import *
from pose import *
from robot import *
#from AlphaMapCost import AlphaMapCost

class TestMap:
	
	def __init__(self):
		self.probe = FakeProbe(40,0.15,0.05)
		self.probe.loadEmpty()
		
		WLEN = 3.0
		wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
		wall2 = [[-4.0 + WLEN*cos(pi/3), 0.2 + WLEN*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
		w1 = wall1[2]
		w2 = wall2[0]
		
		wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
		lp = wall3[0]
		rp = wall3[2]
		
		wall6 = [lp, [lp[0] + WLEN*cos(pi/6), lp[1] + WLEN*sin(pi/6)]]
		wall6.append([wall6[1][0] + 0.4*cos(pi/3), wall6[1][1] - 0.4*sin(pi/3)])
		wall6.append([wall6[2][0] - WLEN*cos(pi/6), wall6[2][1] - WLEN*sin(pi/6)])
		wall6.append([wall6[3][0] + WLEN*cos(pi/3), wall6[3][1] - WLEN*sin(pi/3)])
		wall6.append([wall6[4][0] - 0.4*cos(pi/6), wall6[4][1] - 0.4*sin(pi/6)])
		wall6.append(w1)
		wall6.reverse()
		
		walls = [wall1, wall2, wall3, wall6]
		for wall in walls:
			for i in range(len(wall)):
				p = copy(wall[i])
				p[0] += 6.0
				wall[i] = p
				
		self.probe.addWalls(walls)
				
		self.contacts = ContactReferences(self.probe)
		self.mapGraph = MapGraph(self.probe, self.contacts)
		self.count = 0
		self.num = 0
		
		self.initializeContacts()
		self.mapGraph.loadFile()
		
		self.mapGraph.correctPoses3()
		
		self.mapGraph.saveMap()
		
		
	def step(self):
		self.count += 1

		#print self.count
		
		self.probe.step()
		
		self.mapGraph.keepStablePose()
		#if self.count % 3 == 0:
		if True:
			#print "update", self.count
			if self.num %2 == 1:
				self.mapGraph.update(True)
			else:
				self.mapGraph.update(False)
				
	def load(self, num):
		self.num = num
		self.probe.loadFile(num)
		print "loading file", num, "with", len(self.probe.angles[0]), "snapshots"
		
		self.count = 0

	def newNode(self):
		self.mapGraph.newNode()

	def saveMap(self):
		self.mapGraph.saveLocalMap()

	def initializeContacts(self):
		
		self.contacts.setMask([1.0 for i in range(39)])
		for i in range(100):
			self.contacts.step()
		
if __name__ =="__main__":

	fileNum = 1
	testMap = TestMap()
	#testMap.load(fileNum)
	#testMap.initializeContacts()
	#testMap.newNode()
	#testMap.newNode()

	
	#while True:
	#	try:
	#		testMap.step()
	#	except:
	#		traceback.print_exc()
	#		fileNum += 1
	#		testMap.saveMap()
	#		testMap.load(fileNum)
	#		testMap.newNode()
	
	