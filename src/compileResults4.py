#!/usr/bin/python


import os
import csv
import pylab
import sys
from math import *

globalCount = 0

dirs = []
files = os.listdir(".")
for nameF in files:
	if "results" == nameF[0:7]:
		dirs.append(nameF)
		
dirs.sort()

#for dirName in dirs:

for i in range(17):
		
	angle = i * pi/16 - pi/2
	
	for j in range(0,4):
		
		dist = 0.2*j + 0.5
		#testDir = "results%02u_%f" % (i,dist)
		
		dirName = "./" + "results%02u_%f" % (i, dist)
		#print "./" + "results%02u_%f" % (i, dist)
		files = os.listdir(dirName)
		
		scenes = []
		for nameF in files:
			if "scene" == nameF[0:5]:
				scenes.append(nameF)	
	
		scenes.sort()
	
		for scene in scenes:
			#print "cp " + dirName + "/" + scene + " scene%06u.png" % globalCount
			os.system("cp " + dirName + "/" + scene + " scene%06u.png" % globalCount)
			globalCount += 1
			
		anchorCount = 0
		f = open(dirName + "/out.txt", 'r')
		for line in f:
			if line[:13] == "anchor errors":
				anchorCount += 1
		
		#print dirName, "failed anchors:", anchorCount-1
		print angle, dist, anchorCount-1
	
	#exit()
	
