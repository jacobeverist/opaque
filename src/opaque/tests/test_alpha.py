import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import Image
import ImageDraw
from copy import *
from math import *
from numpy import arange
from maps import *

class Alpha:
	
	def __init__(self, filename):
		self.pixelSize = PIXELSIZE
		self.mapSize = ARMLENGTH*NUM_SEGS + 2.0 + 2.0
		self.numPixel = int(2.0 * self.mapSize / self.pixelSize + 1.0)
		self.halfPix = self.pixelSize / 2.0
		self.divPix = floor((2.0 * self.mapSize / self.pixelSize) / self.mapSize)
		self.mapImage = Image.open(filename)
		self.image = self.mapImage.load()
			
	def realToGrid(self, point):
		indexX = int(floor(point[0] * self.divPix)) + self.numPixel / 2 + 1
		indexY = int(floor(point[1] * self.divPix)) + self.numPixel / 2 + 1
		return indexX, indexY

	def gridToReal(self, indices):
		i = indices[0]
		j = indices[1]
		point = [(i - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0, (j - self.numPixel / 2 - 1) * (self.pixelSize / 2.0) + self.pixelSize / 2.0]
		return point			


	def getImagePoints(self):
		
		" 1. pick out the points "
		numPixel = self.numPixel
		mapImage = self.mapImage
		image = mapImage.load()
		
		points = []
		for j in range(numPixel):
			for k in range(numPixel):
				if image[j,k] == 255:
					pnt = self.gridToReal([j,k])
					points.append(pnt)

		return points

	def computeAlpha2(self, points, radius = 0.2, variance = 0.01):
		
		numPoints = len(points)
		inputStr = str(numPoints) + " "
		
		" radius 0.2 "
		inputStr += str(radius) + " "
		
		for p in points:
			p2 = copy(p)
			p2[0] += gauss(0.0,variance)
			p2[1] += gauss(0.0,variance)

			inputStr += str(p2[0]) + " " + str(p2[1]) + " "
		
		inputStr += "\n"
		
		" start the subprocess "
		subProc = Popen(["../../alpha2.exe"], stdin=PIPE, stdout=PIPE)
		
		
		" send input and receive output "
		sout, serr = subProc.communicate(inputStr)
	
		" convert string output to typed data "
		sArr = sout.split(" ")

		numVert = int(sArr[0])
		
		#print numVert, len(sArr)
		
		vertices = []
		for i in range(numVert+1):
			vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
		
		return vertices


if __name__ =="__main__":
	count = 0
	
	samples = [ 2.0 - i*0.01 for i in range(200)]

	for i in range(1,15):
		
		alpha = Alpha("localOccMap%03u_0000.png" % i)

		points = alpha.getImagePoints()

		for r in samples:

			vert = alpha.computeAlpha2(points, radius = r, variance = 0.0001)
			
			xP = []
			yP = []
			for p in vert:
				xP.append(p[0])
				yP.append(p[1])
			
			pylab.clf()
			pylab.plot(xP,yP)	
			pylab.xlim(-9, 9)
			pylab.ylim(-8, 8)
			pylab.title("Map = %d, Radius = %f" % (i,r))
			pylab.savefig("test_%04u.png" % count)
			count += 1
	
	
