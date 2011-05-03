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
fileName = "mapOccupancyMap%04u.png"

im1 = Image.open(fileName % 0)
im2 = Image.open(fileName % 1)
im3 = Image.open(fileName % 2)
size = (im1.size[0] + im2.size[0] + im3.size[0], im1.size[1])

newImage = Image.new('L', size, 0)
newImage.paste(im1, (0, 0))
newImage.paste(im2, (im1.size[0], 0))
newImage.paste(im3, (im1.size[0] + im2.size[0], 0))

draw = ImageDraw.Draw(newImage)
draw.text((800,0), "35 poses, 70 nodes, Global Map with Constraints", font=font)
draw.text((250,120), "Only Overlap and Motion Constraints", font=font2)
draw.text((250+im1.size[0],120), "+ Hypothesized Overlap Constraints" , font=font2)
draw.text((250+im1.size[0]*2,120), "+ Hypothesized Sensor Constraints", font=font2)
newImage.save("stitchImage.png")
