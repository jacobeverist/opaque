#!/usr/bin/python


import os
import csv
import pylab
import sys
from math import *
import Image

files = os.listdir(".")

#index = files.index("occupancyMap0000.png")

occupancyMaps = []
voronoiMaps = []
obstacleMaps = []
frontierMaps = []
boundaryMaps = []
checkMaps = []
densityMaps = []
scenes = []

for nameF in files:
	if "occupancyMap" == nameF[0:12]:
		occupancyMaps.append(nameF)

	if "densityMap" == nameF[0:10]:
		densityMaps.append(nameF)

	if "frontierMap" == nameF[0:11]:
		frontierMaps.append(nameF)

	if "obstacleMap" == nameF[0:11]:
		obstacleMaps.append(nameF)

	if "checkMap" == nameF[0:8]:
		checkMaps.append(nameF)

	if "voronoiMap" == nameF[0:10]:
		voronoiMaps.append(nameF)

	if "boundaryMap" == nameF[0:11]:
		boundaryMaps.append(nameF)

	if "scene" == nameF[0:5]:
		scenes.append(nameF)

occupancyMaps.sort()
densityMaps.sort()
frontierMaps.sort()
obstacleMaps.sort()
checkMaps.sort()
voronoiMaps.sort()
boundaryMaps.sort()
scenes.sort()

#maps = [occupancyMaps, frontierMaps, obstacleMaps, checkMaps, voronoiMaps, boundaryMaps]
maps = [occupancyMaps, obstacleMaps, voronoiMaps, boundaryMaps]

img = Image.open(occupancyMaps[-1])
size = img.size
img = Image.open(scenes[-1])
size2 = img.size

xWidth = 370-60
yWidth = 630-250

for i in range(len(occupancyMaps)):
	xSize = xWidth * len(maps)/2
	ySize = 2*yWidth + size2[0]
	newImage = Image.new(img.mode, (xSize,ySize))

	for j in range(len(maps)):
		thisImg = Image.open(maps[j][i])
		result = thisImg.crop([60,250,370,630])

		if j % 2 == 0:
			newImage.paste(result,((j/2)*xWidth,0))
		else:
			newImage.paste(result,((j/2)*xWidth,yWidth))

	finalImage = newImage.rotate(90)

	screenImg = Image.open(scenes[i])
	finalImage.paste(screenImg, (2*yWidth,0))

	finalImage.save("composite%04u.png" % i)


exit()


if True:

	str_num = test[5:]
	#print str_num

	try:
		im1 = Image.open(test + "/mapGround.png")
		im2 = Image.open(test + "/mapEstimate.png")
		im3 = Image.open(test + "/mapCorrection.png")

		x1,y1 = im1.size
		x2,y2 = im2.size
		x3,y3 = im3.size

		xNew = x1 + x2 + x3

		newImage.paste(im1,(0,0))
		newImage.paste(im2,(x1,0))
		newImage.paste(im3,(x1+x2,0))
		newImage.save(test + "/map_" + str_num + ".png")
		newImage.save("resultsDir/map_" + str_num + ".png")
	except IOError:
		print "failed " + test


	#os.system("cp " + test + "/poseFitAfter.png resultsDir/map_" + str_num + ".png")
	

