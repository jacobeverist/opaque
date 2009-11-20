#!/usr/bin/python


import os
import csv
import pylab
import sys
from math import *
import Image

files = os.listdir(".")

localOccMaps = []
gndLocalOccMaps = []

for nameF in files:
	if "localOccMap" == nameF[0:11]:
		localOccMaps.append(nameF)

	if "gndLocalOccMap" == nameF[0:14]:
		gndLocalOccMaps.append(nameF)

localOccMaps.sort()
gndLocalOccMaps.sort()

img = Image.open(localOccMaps[-1])
size = img.size
img = Image.open(gndLocalOccMaps[-1])
size2 = img.size

xWidth = size[0]
yWidth = size[1]

for i in range(len(localOccMaps)):
	xSize = 2*xWidth
	ySize = yWidth
	newImage = Image.new(img.mode, (xSize,ySize))

	thisImg = Image.open(localOccMaps[i])
	newImage.paste(thisImg,(0,0))
	thisImg = Image.open(gndLocalOccMaps[i])
	newImage.paste(thisImg,(size[0],0))

	newImage.save("comparison%04u.png" % i)

	#finalImage = newImage.rotate(90)

	#screenImg = Image.open(scenes[i])
	#finalImage.paste(screenImg.resize(size), (yWidth,0))

	#finalImage.save("composite%04u.png" % i)


