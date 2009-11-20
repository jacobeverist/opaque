#!/usr/bin/python

import Image
import ImageDraw
from copy import *
from math import *
from numpy import arange


# GOALS
# 1. probe boundaries to make them confirmed boundaries
# 2. tip of cosine curve, inverse kinematics
# 3. * determine snake's position on voronoi graph
# 4. given desired location, compute closest point to graph, and a shortest path to follow
# 5. * select location to explore from frontier points
# 6. execute path following of snake
# 7. return from exploration if encounter dead-end or wide open space

# 1. couple the frontiers and boundary points with the voronoi graph

# 1. first, make a voronoi graph data structure
# 2. prune the graph according to our criteria
# 3. find the snake's location on the graph
# 4. by position or by features?

	
if __name__ =="__main__":

	class Probe(object):
		numSegs = 40

	probe = Probe()

	im = Image.open("testImage.png")
	rm = RoadMap(im)
	#sm = ShadowMap(probe, 19, oldShadowMap = im)

	#dist, cPoint, minEdge = self.getClosestPoint([-20.0,10.0])
	realPath = rm.computePath([0,0], [5,5])
	#print realPath

	exit()
	#realPath = [[5,0], [0,7]]
	#self.printGraph(realPath)

	# start point
	start = [0.0,0.0]
	currPoint = start
	path = [currPoint]
	try:
		while True:
			#print currPoint
			currPoint = rm.selectExplorationPoint(currPoint)
			rm.markExplored(currPoint)
			path.append(currPoint)
	except:
		im = rm.boundaryImage
		draw = ImageDraw.Draw(im)
		for i in range(len(path)-1):
			indX1, indY1 = rm.realToGrid(path[i])
			indX2, indY2 = rm.realToGrid(path[i+1])
			draw.line([(indX1,indY1),(indX2,indY2)], fill = 255)
		im.save("exploreMap.png")


	#print p1, p2
	#print dist1, dist2

	im.save("selectTest.png") 




