#!/usr/bin/python

import os
import csv
import pylab
import sys
from math import *
import Image

files = os.listdir(".")

mapOccMaps = []
mapObstMaps = []
scenes = []

sceneCount = 0

for nameF in files:
	if "mapOccupancyMap" == nameF[0:15]:
		mapOccMaps.append(nameF)

	if "mapObstGraph" == nameF[0:12]:
		mapObstMaps.append(nameF)

	if "scene" == nameF[0:5]:
		if sceneCount % 10 == 0:
			scenes.append(nameF)		
		sceneCount += 1


mapOccMaps.sort()
mapObstMaps.sort()
scenes.sort()

img = Image.open(mapOccMaps[0])
size1 = img.size
img = Image.open(mapObstMaps[0])
size2 = img.size
img = Image.open(scenes[0])
size3 = img.size

xSize = size3[0]
ySize = size3[1] + size3[0]/2

print xSize, ySize

saveOcc = newImage = Image.new(img.mode, size1)
saveObst = newImage = Image.new(img.mode, size2)

for i in range(len(scenes)):
	newImage = Image.new(img.mode, (xSize,ySize))

	sceneStat = os.stat(scenes[i])
	sceneTime = sceneStat[8]

	if len(mapObstMaps) > 0:
		obstStat = os.stat(mapObstMaps[0])
		obstTime = obstStat[8]
	else:
		obstTime = 1e100

	if len(mapOccMaps) > 0:
		occStat = os.stat(mapOccMaps[0])
		occTime = occStat[8]
	else:
		occTime = 1e100

	sceneImg = Image.open(scenes[i])

	if occTime < sceneTime:
		occImg = Image.open(mapOccMaps[0])
		#draw = ImageDraw.Draw(frontierDensity)
		#draw.ellipse((xMax-5, yMax-5, xMax+5, yMax+5), outline=255)

		occResizeImg = occImg.resize((size3[0]/2,size3[0]/2))
		saveOcc = occImg
		mapOccMaps = mapOccMaps[1:]
	else:
		occResizeImg = saveOcc.resize((size3[0]/2,size3[0]/2))

	if obstTime < sceneTime:
		obstImg = Image.open(mapObstMaps[0])
		obstResizeImg = obstImg.resize((size3[0]/2,size3[0]/2))
		saveObst = obstImg
		mapObstMaps = mapObstMaps[1:]
	else:
		obstResizeImg = saveObst.resize((size3[0]/2,size3[0]/2))

	newImage.paste(sceneImg,(0,0))
	newImage.paste(occResizeImg,(0,size3[1]))
	newImage.paste(obstResizeImg,(size3[0]/2,size3[1]))

	#print  (size3[0]/2,size3[1])

	newImage.save("composite%04u.png" % i)

exit()

for i in range(len(localOccMaps)):
	xSize = xWidth
	ySize = 2*yWidth
	newImage = Image.new(img.mode, (xSize,ySize))

	thisImg = Image.open(localOccMaps[i])
	#result = thisImg.crop([60,250,370,630])
	#newImage.paste(result,(0,0))
	newImage.paste(thisImg,(0,0))

	finalImage = newImage.rotate(90)

	screenImg = Image.open(scenes[i])
	finalImage.paste(screenImg.resize(size), (yWidth,0))

	finalImage.save("composite%04u.png" % i)




