import Image
from medialaxis import computeMedialAxis
from SplineFit import SplineFit
import pylab
import random
import graph
from copy import copy
from math import sqrt
#import numpy
import functions
import math

#img1 = Image.new('L', (5,5))
#img2 = Image.new('L', (5,5))
img1 = Image.open("testInput1.png")
img2 = Image.open("testInput2.png")
img3 = Image.open("testInput3.png")


def getLongestPath(node, currSum, currPath, tree, isVisited, nodePath, nodeSum):

	if isVisited[node]:
		return

	isVisited[node] = 1
	nodePath[node] = copy(currPath) + [node]

	if nodeSum[node] < currSum:
		nodeSum[node] = currSum

	for childNode in tree[node]:
		getLongestPath(childNode, currSum + 1, nodePath[node], tree, isVisited, nodePath, nodeSum)

	isVisited[node] = 0

def getMedialAxis(hull):

	PIXELSIZE = INC = 0.05
	mapSize = 0.15*40 + 2.0 + 2.0
	pixelSize = PIXELSIZE
	numPixel = int(2.0*mapSize / pixelSize + 1.0)
	divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)

	def realToGrid(point):
		indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
		indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
		return indexX, indexY

	def gridToReal(indices):
		i = indices[0]
		j = indices[1]
		point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
		return point

	print "hull:", hull

	gridHull = []
	for p in hull:
		gridHull.append(realToGrid(p))

	print "gridHull:", gridHull
	
	minX = 1e100
	maxX = -1e100
	minY = 1e100
	maxY = -1e100
	for p in gridHull:
		if p[0] > maxX:
			maxX = p[0]
		if p[0] < minX:
			minX = p[0]
		if p[1] > maxY:
			maxY = p[1]
		if p[1] < minY:
			minY = p[1]

	print "minMax:", minX, maxX, minY, maxY
			
	xRange = range(minX, maxX + 1)
	yRange = range(minY, maxY + 1)

	print "xRange:", xRange
	print "yRange:", yRange
	interior = []

	for i in xRange:
		for j in yRange:
			if functions.point_inside_polygon(i, j, gridHull):
				interior.append((i,j))


	inputImg = Image.new('L', (numPixel,numPixel), 255)
	imga = inputImg.load()
	
	for p in gridHull:
		imga[p[0],p[1]] = 0
	
	for p in interior:
		imga[p[0],p[1]] = 0
		#i, j = realToGrid(p)
		#imga[i,j] = 0

	resultImg = Image.new('L', inputImg.size)
	resultImg = computeMedialAxis(inputImg, resultImg)
	imgA = resultImg.load()

	points = []
	for i in range(1, inputImg.size[0]-1):
		for j in range(1, inputImg.size[1]-1):
			if imgA[i,j] == 0:
				points.append((i,j))

	medialGraph = graph.graph()
	for p in points:
		medialGraph.add_node(p, [])




	builtGraph = {}
	for i in range(2, img1.size[0]-2):
		for j in range(2, img1.size[1]-2):
			if imgA[i,j] == 0:
				builtGraph[(i,j)] = []
				for k in range(i-1, i+2):
					for l in range(j-1, j+2):
						if imgA[k,l] == 0:
							builtGraph[(i,j)].append((k,l))
							medialGraph.add_edge((i,j), (k,l))
							


	mst = medialGraph.minimal_spanning_tree()

	uni_mst = {}
	isVisited = {}
	nodeSum = {}
	for k, v in mst.items():
		uni_mst[k] = []
		isVisited[k] = 0
		nodeSum[k] = 0

	for k, v in mst.items():
		if v != None:
			uni_mst[k].append(v)
			uni_mst[v].append(k)

	leaves = []
	for k, v in uni_mst.items():
		if len(v) == 1:
			leaves.append(k)


	" find the longest path from each leaf"
	maxPairDist = 0
	maxPair = None
	maxPath = []
	for leaf in leaves:

		isVisited = {}
		nodeSum = {}
		nodePath = {}
		for k, v in uni_mst.items():
			isVisited[k] = 0
			nodeSum[k] = 0
			nodePath[k] = []

		getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)

		maxDist = 0
		maxNode = None
		for k, v in nodeSum.items():
			#print k, v
			if v > maxDist:
				maxNode = k
				maxDist = v

		#print leaf, "-", maxNode, maxDist

		if maxDist > maxPairDist:
			maxPairDist = maxDist
			maxPair = [leaf, maxNode]
			maxPath = nodePath[maxNode]


	frontVec = [0.,0.]
	backVec = [0.,0.]
	indic = range(6)
	indic.reverse()

	for i in indic:
		p1 = maxPath[i+2]
		p2 = maxPath[i]
		vec = [p2[0]-p1[0], p2[1]-p1[1]]
		frontVec[0] += vec[0]
		frontVec[1] += vec[1]

		p1 = maxPath[-i-3]
		p2 = maxPath[-i-1]
		vec = [p2[0]-p1[0], p2[1]-p1[1]]
		backVec[0] += vec[0]
		backVec[1] += vec[1]

	frontMag = sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
	backMag = sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])

	frontVec[0] /= frontMag
	frontVec[1] /= frontMag
	backVec[0] /= backMag
	backVec[1] /= backMag

	newP1 = (maxPath[0][0] + frontVec[0]*50, maxPath[0][1] + frontVec[1]*50)
	newP2 = (maxPath[-1][0] + backVec[0]*50, maxPath[-1][1] + backVec[1]*50)

	maxPath.insert(0,newP1)
	maxPath.append(newP2)


	" convert path points to real "
	realPath = []
	for p in maxPath:
		realPath.append(gridToReal(p))


	return realPath


hull = [(2.0,1.0), (2.0,-1.0), (-1.0,-1.0), (-1.0,1.0)]
maxPath  = getMedialAxis(hull)

for i in range(len(hull)-1):
	xP = [hull[i][0], hull[i+1][0]]
	yP = [hull[i][1], hull[i+1][1]]
	pylab.plot(xP,yP,color='blue')
xP = [hull[-1][0], hull[0][0]]
yP = [hull[-1][1], hull[0][1]]
pylab.plot(xP,yP,color='blue')

for i in range(len(maxPath)-1):
	xP = [maxPath[i][0], maxPath[i+1][0]]
	yP = [maxPath[i][1], maxPath[i+1][1]]
	pylab.plot(xP,yP,color='red')

#spl = SplineFit(points)
#splPoints = spl.getUniformSamples()
#xP = []
#yP = []
#for p in splPoints:
#	xP.append(p[0])
#	yP.append(p[1])
#pylab.plot(xP,yP)

#pylab.xlim(0,400)
#pylab.ylim(0,400)
pylab.xlim(-4,4)
pylab.ylim(-4,4)
pylab.show()

#print points
exit()


#resultImg.save("testOutput1.png")

resultImg = computeMedialAxis(img2, resultImg)
#resultImg.save("testOutput2.png")

resultImg = computeMedialAxis(img3, resultImg)
#resultImg.save("testOutput3.png")
print "done"
