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
		self.contacts = ContactReferences(self.probe)
		self.mapGraph = MapGraph(self.probe, self.contacts)
		self.count = 0
		self.num = 0
		
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
	testMap.load(fileNum)
	testMap.initializeContacts()
	testMap.newNode()
	testMap.newNode()

	"""
	try:
		while True:
			testMap.step()

	except:
		traceback.print_exc()
		fileNum += 1
		testMap.saveMap()
	"""
	
	while True:
		try:
			testMap.step()
		except:
			traceback.print_exc()
			fileNum += 1
			testMap.saveMap()
			testMap.load(fileNum)
			testMap.newNode()
			
	
