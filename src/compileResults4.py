#!/usr/bin/python


import os
import csv
import pylab
import sys
from math import *

globalCount = 0

for i in range(0,16):
	for j in range(0,4):
		dist = 0.2*j + 0.5
		
		dirName = "./" + "results%02u_%f" % (i, dist)
		#print "./" + "results%02u_%f" % (i, dist)
		files = os.listdir(dirName)
		
		scenes = []
		for nameF in files:
			if "scene" == nameF[0:5]:
				scenes.append(nameF)	

		scenes.sort()

		for scene in scenes:
			print "cp " + dirName + "/" + scene + " scene%06u.png" % globalCount
			#os.system("cp " + dirName + "/" + scene + " scene%06u.png" % globalCount)
			globalCount += 1

#exit()

