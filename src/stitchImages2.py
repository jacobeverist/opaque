#!/usr/bin/python


import os
import csv
import sys
from math import *
import Image
import ImageDraw
import ImageFont

# use a truetype font
font = ImageFont.truetype("arial.ttf", 24)
font2 = ImageFont.truetype("arial.ttf", 16)

numPoses = 3

im1 = Image.open("localOccMap%03u_0000.png" % 0)
im2 = Image.open("stableLocalOccMap%03u_0000.png" % 0)
im3 = Image.open("gndLocalOccMap%03u_0000.png" % 0)
size = (im1.size[0] + im2.size[0] + im3.size[0], im1.size[1]*numPoses)

newImage = Image.new('L', size, 0)

for i in range(0,numPoses):

	im1 = Image.open("localOccMap%03u_0000.png" % i)
	im2 = Image.open("stableLocalOccMap%03u_0000.png" % i)
	im3 = Image.open("gndLocalOccMap%03u_0000.png" % i)

	newImage.paste(im1, (0,i*im1.size[1]))
	newImage.paste(im2, (im1.size[0],i*im1.size[1]))
	newImage.paste(im3, (im1.size[0] + im2.size[0],i*im1.size[1]))

	draw = ImageDraw.Draw(newImage)
	draw.text((size[0]/2-60,0), "Local Map %d" % i, font=font)
	draw.text((130,60), "Uncorrected Local Map", font=font2)
	draw.text((130+im1.size[0],60), "Corrected Local Map", font=font2)
	draw.text((130+im1.size[0]*2,60), "Ground Truth Local Map", font=font2)
	#draw.text((130+im1.size[0],60), "Ground Truth Local Map", font=font2)

	newImage.save("stitchImage_%03u.png" % i)
