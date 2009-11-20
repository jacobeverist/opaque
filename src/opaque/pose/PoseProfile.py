#!/usr/bin/python

import math
from math import *
from random import *
from scipy.optimize import *
import scipy.interpolate
import numpy
import pylab
import getopt, sys, os, shutil, csv
from copy import *

import traceback



class PoseProfile:

	def __init__ (self, contacts, rootNode):

		self.rootNode = rootNode
		self.contacts = contacts
		self.allNodes = []
		for i in range(-1,self.contacts.numJoints+1):
			self.allNodes.append(self.contacts.getClosestPose(i))		


		self.separatePoints()

	def separatePoints(self):

		sadPoints = [0, 4, 8, 12, 16, 20, 24, 28, 32, 36]
		infPoints = [2, 6, 10, 14, 18, 22, 26, 30, 34, 38]
		
		self.L_pose = []
		self.R_pose = []
		self.I_pose = []

		for i in range(9):
			j = infPoints[i]
			pose = self.allNodes[j+1]
			self.I_pose.append(copy(pose))

		for i in range(9):
			j = sadPoints[i]
			pose = self.allNodes[j+1]
			
			if i % 2 == 1 :
				self.L_pose.append(copy(pose))
			else:
				self.R_pose.append(copy(pose))
		
		self.rootPose = self.allNodes[self.rootNode+1]

		# point information is respectively for joints:
		# 4,8,12,16,20,24,28,32,36 

	def saveToFile(self, filename):

		#val = repr([self.L_pose, self.R_pose, self.I_pose])
		val = repr(self.allNodes)
		
		f = open(filename, 'w')
		f.write(val)
		f.close()
	
	def readFromFile(self, filename):
		f = open(filename, 'r')
		val = eval(f.read())
		f.close()
		
		self.allNodes = val
		
		self.separatePoints()
		
	def printPosePoints(self, clr):

		xP = []
		yP = []
		
		for p in self.allNodes:
			xP.append(p[0])
			yP.append(p[1])

		# plot all points by (x,y)
		#for pose in [self.L_pose, self.R_pose, self.I_pose]:
		#	for p in pose:
		#		xP.append(p[0])
		#		yP.append(p[1])

		#pylab.scatter(xP, yP, color=clr, faceted=False)

		# flipped to align with Ogre x-z coord system
		pylab.scatter(xP, yP, color=clr, faceted=False, linewidth=1)

	def createConstraints(self, profile2, v_list):
		# create constraints between this and profile2

		nomVar = 0.08

		L_pose2 = profile2.getLPose()
		R_pose2 = profile2.getRPose()
		I_pose2 = profile2.getIPose()

		edgeConstraints = []

		for pose in [self.L_pose, self.R_pose, self.I_pose]:
			for p1 in pose:
				
				xNom = -p1[0]
				yNom = -p1[1]

				# rotation angle is not in the pose, needs to be retrieved from the vertex list
				nodeID1 = p1[2]
				pose1 = v_list[nodeID1][1]
				rNom = -pose1[2]

				for pose2 in [L_pose2, R_pose2, I_pose2]:
					for p2 in pose2:
						tPose = deepcopy(p2)

						nodeID2 = p2[2]
						pose2 = v_list[nodeID2][1]
						tPose[2] = pose2[2]

						tPose[0] += xNom
						tPose[1] += yNom
						tPose[2] += rNom

						# rotate the p2 around the origin by rNom
						xnew = tPose[0]*cos(rNom) - tPose[1]*sin(rNom)
						ynew = tPose[0]*sin(rNom) + tPose[1]*cos(rNom)

						xDiff = xnew
						yDiff = ynew
						#xDiff = tPose[0]
						#yDiff = tPose[1]

						# FIXME: normalize angle?
						rDiff = tPose[2]

						# add the constraint
						constraint = [nodeID1, nodeID2, xDiff, yDiff, rDiff, nomVar, nomVar, nomVar, nomVar, nomVar, nomVar]
						edgeConstraints.append(constraint)

		return edgeConstraints

	def getMaxMin(self):

		xHigh = -10e100
		xLow = 10e100
		yHigh = -10e100
		yLow = 10e100

		for pose in [self.L_pose, self.R_pose, self.I_pose]:
			for p in pose:
				if p[0] > xHigh:
					xHigh = p[0]
				if p[0] < xLow:
					xLow = p[0]
				if p[1] > yHigh:
					yHigh = p[1]
				if p[1] < yLow:
					yLow = p[1]
		
		return xHigh, xLow, yHigh, yLow


	def performOffset(self, offset, refPoint):

		origin = refPoint
		ang = offset[2]
		
		for p in self.allNodes:
			p[0] += offset[0]
			p[1] += offset[1]

			xdiff = p[0] - origin[0]
			ydiff = p[1] - origin[1]

			xnew = xdiff*cos(ang) - ydiff*sin(ang)
			ynew = xdiff*sin(ang) + ydiff*cos(ang)

			p[0] = origin[0] + xnew
			p[1] = origin[1] + ynew
			p[2] += ang

		self.separatePoints()
		
		"""

		# translate all points by (x,y)
		for pose in [self.L_pose, self.R_pose, self.I_pose]:
			for p in pose:
				p[0] += offset[0]
				p[1] += offset[1]

		self.rootPose[0] += offset[0]
		self.rootPose[0] += offset[1]

		# rotate all points with respect to p0 of L_pose by offset[2]
		origin = refPoint
		ang = offset[2]
		for pose in [self.L_pose, self.R_pose, self.I_pose]:
			for p in pose:
				xdiff = p[0] - origin[0]
				ydiff = p[1] - origin[1]

				xnew = xdiff*cos(ang) - ydiff*sin(ang)
				ynew = xdiff*sin(ang) + ydiff*cos(ang)

				p[0] = origin[0] + xnew
				p[1] = origin[1] + ynew

		xdiff = self.rootPose[0] - origin[0]
		ydiff = self.rootPose[1] - origin[1]
		
		xnew = xdiff*cos(ang) - ydiff*sin(ang)
		ynew = xdiff*sin(ang) + ydiff*cos(ang)
		
		self.rootPose[0] = origin[0] + xnew
		self.rootPose[1] = origin[1] + ynew
		self.rootPose[2] += ang
		"""
		
	def getLPose(self):
		return self.L_pose

	def getRPose(self):
		return self.R_pose

	def getIPose(self):
		return self.I_pose
	
	def getRootPose(self):
		return copy(self.rootPose)

	def setRootPose(self, newRoot):
		self.rootPose = copy(newRoot)
