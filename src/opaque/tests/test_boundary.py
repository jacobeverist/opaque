#!/usr/bin/python
import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from copy import *
from math import *
import graph
import pylab

import Image
from maps import *


class OccMap:
	
	def __init__(self, fileName):
	
		self.mapImage = Image.open(fileName)
		self.image = self.mapImage.load()

		self.xMin = 10e6
		self.xMax = 0
		self.yMin = 10e6
		self.yMax = 0

		for i in range(self.mapImage.size[0]):
			for j in range(self.mapImage.size[1]):
				if self.image[i,j] == 255:
					if i < self.xMin:
						self.xMin = i
					if i > self.xMax:
						self.xMax = i
					if j < self.yMin:
						self.yMin = j
					if j > self.yMax:
						self.yMax = j	
		
	def getImage(self):
		return self.mapImage
	
		
if __name__ =="__main__":

	occMap = OccMap("occ0008.png")
	boundMap = FreeSpaceBoundaryMap(occMap)
	boundMap.setMapSize(10.0)
	boundMap.update(occMap)
	bGraph = boundMap.getBoundaryGraph()
	
	mst = bGraph.minimal_spanning_tree(boundMap.rootNode)
	print len(mst)
	
	"""	
	vNodes = copy(bGraph.nodes()).tolist()
	while len(vNodes) > 0:
		mst = bGraph.minimal_spanning_tree(vNodes[0])
		print len(mst)
		for key in mst.keys():
			vNodes.remove(key)
	"""
	print mst.edges()
	#vEdges = bGraph.edges()
	vEdges = mst
	for edge in vEdges:
		print edge
		p0 = eval(edge[0])
		p1 = eval(edge[1])
		
		xP = [p0[0],p1[0]]
		yP = [p0[1],p1[1]]
		pylab.plot(xP,yP, color='0.0')
				
	pylab.show()
	
	
	#print len(vEdges)
	
	#print boundMap.rootNode
	#boundMap.saveMap()
	
	
	
	
	
	