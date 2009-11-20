#!/usr/bin/python


import os
import csv
#import pylab
import numpy
from numpy import *
import sys
from math import *
import Image
from matplotlib.pyplot import *
from matplotlib.collections import LineCollection
import pylab

def normalizeAngle(angle):

	while angle>pi:
		angle=angle-2*pi

	while angle<=-pi:
		angle=angle+2*pi

	return angle 


"""
1. Plot the global configuration of each segment at each time increment
2. Plot the local configuration with respect to a root node

"""

segLength = 0.15
segWidth = 0.05
rootNode = 19
numSegs = 40


#for j in range(1,16):
for j in range(7,16):
	os.system("rm poseChange*")

	angleFile = open("angle_output_%04u.txt" % j, 'r')
	poseFile = open("global_pose_output_%04u.txt" % j, 'r')

	poses = []
	angles = []

	count = 0
	for line in poseFile:
		if count % 2 == 0:
			poseSnapshot = eval(line)
			poses.append(poseSnapshot)
		count += 1

	count = 0
	for line in angleFile:
		if count % 2 == 0:
			angleSnapshot = eval(line)
			angles.append(angleSnapshot)
		count += 1
	
	#print len(poses), len(angles)

	estPose = poses[0][rootNode+1] 

	#estAngle = normalizeAngle(estPose[2] + pi/2)

	estAngle = -estPose[2]
	dist = sqrt(estPose[0]**2 + estPose[1]**2)
	vecAng = acos(estPose[0]/dist)
	if asin(estPose[1]/dist) < 0:
		vecAng = -vecAng

	backR = array([[cos(vecAng), sin(vecAng)],[-sin(vecAng),cos(vecAng)]])
	foreR = array([[cos(vecAng), -sin(vecAng)],[sin(vecAng),cos(vecAng)]])

	#R = array([[cos(estPose[2]), sin(estPose[2])],[-sin(estPose[2]),cos(estPose[2])]])
	R = array([[cos(estAngle), sin(estAngle)],[-sin(estAngle),cos(estAngle)]])


	for k in range(len(poses)):
		pose = poses[k]
		pylab.clf()

		for i in range(numSegs):

			xTotal = pose[i][0]
			zTotal = pose[i][1]
			totalAngle = pose[i][2]

			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal + segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal + segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal + segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]
			points = [p1,p2,p3,p4]

			for p in points:
				finalVec = array([[p[0]], [p[1]]])
				transVec = dot(transpose(R), finalVec)
				resVec = dot(backR, transVec)
				resVec[0, 0] += dist
				tempVec = dot(foreR, resVec)
				p[0] = tempVec[0, 0]
				p[1] = tempVec[1, 0]	

			xP = []
			yP = []
			for m in range(5):
				xP.append(points[m%4][0])
				yP.append(points[m%4][1])

			#xP = [p4[0], p3[0], p2[0], p1[0], p4[0]]
			#yP = [p4[1], p3[1], p2[1], p1[1], p4[1]]

			if i == (rootNode+1):
				#pylab.plot(xP,yP, color='r')
				pylab.fill(xP,yP, 'r')
			else:
				pylab.plot(xP,yP, color='b')


		xlim(-3,3)
		ylim(-3,3)
		title("Global Pose %d" % k)
		#savefig("globalPose_%04u.png" % k)
		savefig("globalPose.png")


		" now draw the local pose "
		angle = angles[k]
		pylab.clf()

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0
		actualConfig = []
	
		joints = range(-1,rootNode)
		joints.reverse()

		for i in joints:

			totalAngle = totalAngle + angle[i+1]
			totalAngle = normalizeAngle(totalAngle)

			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			xTotal = xTotal - segLength*cos(totalAngle)
			zTotal = zTotal - segLength*sin(totalAngle)

			actualConfig.append([i, [p4,p3,p2,p1]])

		joints = range(rootNode, numSegs-1)

		xTotal = 0.0
		zTotal = 0.0
		totalAngle = 0.0

		for i in joints:
			xTotal = xTotal + segLength*cos(totalAngle)
			zTotal = zTotal + segLength*sin(totalAngle)
		
			p1 = [xTotal - 0.5*segWidth*sin(totalAngle), zTotal + 0.5*segWidth*cos(totalAngle)]
			p2 = [xTotal + 0.5*segWidth*sin(totalAngle), zTotal - 0.5*segWidth*cos(totalAngle)]
			p3 = [xTotal - segLength*cos(totalAngle) + 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) - 0.5*segWidth*cos(totalAngle)]
			p4 = [xTotal - segLength*cos(totalAngle) - 0.5*segWidth*sin(totalAngle), zTotal - segLength*sin(totalAngle) + 0.5*segWidth*cos(totalAngle)]

			actualConfig.append([i, [p4,p3,p2,p1]])
			
			if i < numSegs-2:
				totalAngle = totalAngle - angle[i+1]
				totalAngle = normalizeAngle(totalAngle)

		actualConfig.sort()

		#for rect in actualConfig:
		for i in range(len(actualConfig)):
			rect = actualConfig[i]
			xP = []
			yP = []
			for m in range(5):
				xP.append(rect[1][m%4][0])
				yP.append(rect[1][m%4][1])

			if i == (rootNode+1):
				#pylab.plot(xP,yP, color='r')
				pylab.fill(xP,yP, 'r')
			else:
				pylab.plot(xP,yP, color='b')
				#pylab.fill(xP,yP, 'b')

			xlim(-3,3)
			ylim(-3,3)
			title("Local Pose %d" % k)
			#savefig("localPose_%04u.png" % k)
			savefig("localPose.png")

	
		img1 = Image.open("globalPose.png")
		img2 = Image.open("localPose.png")

		xSize = img1.size[0] + img2.size[0]
		ySize = img1.size[1]

		imgFinal = Image.new('RGBA', (xSize-220,ySize-40))
		#imgFinal = Image.new('RGBA', (xSize,ySize))
		imgFinal.paste(img2, (img1.size[0]-160,-20))
		imgFinal.paste(img1, (-70,-20))
		imgFinal.save("poseChange_%04u.png" % k)

		imgFinal = 0
		img1 = 0
		img2 = 0


	#system("ffmpeg -i poseChange_\%4d.png globalLocalCompare_%02u.avi" % j)
	os.system('ffmpeg -i poseChange_%4d.png '+ "globalLocalCompare_%02u.avi" % j)	
