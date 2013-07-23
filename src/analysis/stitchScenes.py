#!/usr/bin/python


import os
import csv
import sys
from math import *
import Image
import ImageDraw
import ImageFont

# use a truetype font
font = ImageFont.truetype("arial.ttf", 48)
font2 = ImageFont.truetype("arial.ttf", 32)

#fileName = "mapBoundaryMap%04u.png"
fileName = "scene%06u.png"

for i in range(1, 105):

	if i >= 64:
		im1 = Image.open(fileName % 63)
	else:
		im1 = Image.open(fileName % i)
	im2 = Image.open("PathStepVideo/" + fileName % i)
	size = (im1.size[0] + im2.size[0], im1.size[1])

	newImage = Image.new('L', size, 0)
	newImage.paste(im1, (0, 0))
	newImage.paste(im2, (im1.size[0], 0))
	newImage.save("pathScene%06u.png" % i)
