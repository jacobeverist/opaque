#!/usr/bin/python


import os
import csv
import pylab
import sys
from math import *
import Image

files = os.listdir(".")

occupancyMaps = []
localOccMaps = []
scenes = []

for nameF in files:
	if "occupancyMap" == nameF[0:12]:
		occupancyMaps.append(nameF)

	if "localOccMap" == nameF[0:11]:
		localOccMaps.append(nameF)

	if "scene" == nameF[0:5]:
		scenes.append(nameF)

occupancyMaps.sort()
localOccMaps.sort()
scenes.sort()

#img = Image.open(occupancyMaps[-1])
img = Image.open(localOccMaps[-1])
size = img.size
img = Image.open(scenes[-1])
size2 = img.size

#xWidth = 370-60
#yWidth = 630-250
xWidth = size[0]
yWidth = size[1]

print xWidth, yWidth

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


