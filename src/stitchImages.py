#!/usr/bin/python


import glob
import os
import csv
import sys
from math import *
import Image
import ImageDraw
import ImageFont

# use a truetype font
#font = ImageFont.truetype("arial.ttf", 24)
#font2 = ImageFont.truetype("arial.ttf", 16)


sceneFiles = glob.glob(os.path.join("scene*.png"))
navFiles = glob.glob(os.path.join("navEst*.png"))

xCropMargin = 120 
cropMargin = 80
im1 = Image.open(sceneFiles[0])
im2 = Image.open(navFiles[0])
size = (im1.size[0] + im2.size[0]-2*xCropMargin, max(im1.size[1],im2.size[1]-2*cropMargin))

files = sceneFiles + navFiles

files.sort(key=lambda x: os.path.getmtime(x))

srcImage1 = Image.new('RGB', im1.size, 0)
srcImage2 = Image.new('RGB', im2.size, 0)

indexCount = 0


for inFile in files:

	if inFile[:5] == "scene":
		srcImage1 = Image.open(inFile)
		print "scene:", inFile

	if inFile[:5] == "navEs":
		srcImage2 = Image.open(inFile)
		#srcImage2 = srcImage2.crop((0,0,im2.size[0],im2.size[1]))
		srcImage2 = srcImage2.crop((xCropMargin,cropMargin,im2.size[0]-xCropMargin,im2.size[1]-cropMargin))

 		#left, upper, right, and lower pixel coordinate.

		print "nav:", inFile

	newImage = Image.new('RGB', size, 0)
	newImage.paste(srcImage1, (0,0))
	newImage.paste(srcImage2, (srcImage1.size[0],0))

	newImage.save("stitchImage_%05u.png" % indexCount)
	indexCount += 1

#for inFile in files:
#	print inFile

exit()


#for i in range(1,31):
for i in range(0,10):

	im1 = Image.open("localOccMap%03u_0000.png" % i)
	im2 = Image.open("stableLocalOccMap%03u_0000.png" % i)
	im3 = Image.open("gndLocalOccMap%03u_0000.png" % i)


	size = (im1.size[0] + im2.size[0] + im3.size[0], im1.size[1])
	#size = (im1.size[0] + im3.size[0], im1.size[1])

	newImage = Image.new('L', size, 0)
	newImage.paste(im1, (0,0))
	#newImage.paste(im3, (im1.size[0],0))
	newImage.paste(im2, (im1.size[0],0))
	newImage.paste(im3, (im1.size[0] + im2.size[0],0))

	draw = ImageDraw.Draw(newImage)
	draw.text((size[0]/2-60,0), "Local Map %d" % i, font=font)
	draw.text((130,60), "Uncorrected Local Map", font=font2)
	draw.text((130+im1.size[0],60), "Corrected Local Map", font=font2)
	draw.text((130+im1.size[0]*2,60), "Ground Truth Local Map", font=font2)
	#draw.text((130+im1.size[0],60), "Ground Truth Local Map", font=font2)

	newImage.save("stitchImage_%03u.png" % i)
