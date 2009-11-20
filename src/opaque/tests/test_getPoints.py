#!/usr/bin/python
import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from copy import *
from math import *
from random import *
import graph
import pylab

import Image
from maps import *


class OccMap:
	
	def __init__(self, fileName):
	
		self.mapImage = Image.open(fileName)
		self.image = self.mapImage.load()

		mapSize = 10.0
		self.pixelSize = PIXELSIZE
		self.mapSize = mapSize
		self.numPixel = int(2.0*self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize/2.0
		self.divPix = floor((2.0*self.mapSize/self.pixelSize)/self.mapSize)

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
	
	def realToGrid(self, point):
		indexX = int(floor(point[0]*self.divPix)) + self.numPixel/2 + 1
		indexY = int(floor(point[1]*self.divPix)) + self.numPixel/2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0, (j - self.numPixel/2 - 1)*(self.pixelSize/2.0) + self.pixelSize/2.0]
		return point	
	
	def getPoints(self):
		points = []
		for i in range(self.numPixel):
			for j in range(self.numPixel):
				if self.image[i,j] == 255:
					p = self.gridToReal((i,j))
					#print gauss(0.0,0.1)
					p[0] += gauss(0.0,0.01)
					p[1] += gauss(0.0,0.01)
					points.append(p)
		
		return points
	
if __name__ =="__main__":

	occMap = OccMap("occ0008.png")
	points = occMap.getPoints()
	print len(points)
	f = open("point.txt", 'w')
	
	for p in points:
		text = str(p[0]) + "," + str(p[1]) + ","
		f.write(text)
		#print p[0], ",", p[1], ","

	f.write("\n")
	f.close()
	
	f = open("alpha.txt", 'r')
	pnt_str = f.read()
	
	edges = eval(pnt_str)
	print len(edges)
	
	xP = []
	yP = []
	for e in edges:
		xP = [e[0][0], e[1][0]]
		yP = [e[0][1], e[1][1]]
		pylab.plot(xP,yP, color = '0.5')
	pylab.show()
	
	"""
	xP = []
	yP = []
	for p in points:
		xP.append(p[0])
		yP.append(p[1])
		
	pylab.scatter(xP,yP, color='0.0')
				
	pylab.show()
	"""

	
	
	
	
	