#!/usr/bin/python
from copy import *
from math import *
from operator import itemgetter
import pylab
import matplotlib.pyplot as plt
import re

import os
import sys

relPath = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(1,relPath)
sys.path.insert(1,relPath + "/modules/")
sys.path.insert(1,relPath + "/modules/nelmin")
sys.path.insert(1,relPath + "/modules/bulletprobe")
sys.path.insert(1,relPath + "/modules/medialaxis")
sys.path.insert(1,relPath + "/modules/alphashape")
sys.path.insert(1,relPath + "/modules/ogre")

from opaque.maps.Pose import Pose

landmarks_G = [([0.96806901939298162, 0.048843927068442317], 0.05, 'bloomPoint'), ([0.81306950256992039, 0.26361237137896043], 0.9, 'bendPoint'), ([0.93724884235706052, 0.24624709762295777], 0.9, 'bendPoint'), ([0.90975576242926481, 0.34376059543918064], 0.9, 'bendPoint'), ([0.81009612809822829, 0.39441732399057705], 0.9, 'bendPoint'), ([0.75421546419128904, 0.39066879832905799], 0.9, 'bendPoint'), ([0.78438799459538622, 0.40917698626228294], 0.9, 'bendPoint'), ([0.77032421148408858, 0.33374972957263604], 0.9, 'bendPoint'), ([0.81685493633563711, 0.39136594771127647], 0.9, 'bendPoint'), ([0.65672740205892666, 0.35468083949259266], 0.9, 'bendPoint'), ([0.77358490516371603, 0.3775871428150322], 0.9, 'bendPoint'), ([0.7560999263991306, 0.3340831765232033], 0.9, 'bendPoint'), ([0.71645308353069437, 0.40593586936737613], 0.9, 'bendPoint'), ([1.1649867773380371, 0.1950193739918461], 0.05, 'bloomPoint'), ([0.96086754721010181, 0.20711062755828463], 0.05, 'bloomPoint'), ([1.0654070300437992, 0.75800085613184354], 0.05, 'archPoint'), ([1.0921994772580448, 0.32991657117451795], 0.05, 'archPoint')]

LANDMARK_THRESH = 7.0
#OUTLIER_THRESH = 1.5
OUTLIER_THRESH = 0.7
#OUTLIER_THRESH = 0.4

#f = open('landmark_data2.txt', 'r')

#f = open('/home/everist/opaque/doc/results/tip_results_2014_11_23/cross_0.0/testDir/out.txt', 'r')
f = open('/home/everist/opaque/doc/results/tip_results_2014_11_23/cross_0.4/testDir/out.txt', 'r')
#f = open('/home/everist/opaque/doc/results/tip_results_2014_11_23/T_bottom_0.4/testDir/out.txt', 'r')
#f = open('testDir6/out.txt', 'r')
plotCount = 0

#prevLine = f.readline()
#parsedLine = re.sub(r'^[0-9]+ landmarks: ', r'', dataLine)


resultLine = f.readline()
resultStrings = []
while resultLine != "":

	#print resultLine

	if re.search(r'^[0-9]+ landmarks: ', resultLine):
		resultStrings.append(resultLine)

	if re.search(r'^finalPoses: ', resultLine):
		resultStrings.append(resultLine)

	resultLine = f.readline()

for result in resultStrings:
	print result


for k in range(len(resultStrings)-1):

	firstLine = resultStrings[k]
	secondLine = resultStrings[k+1]

	#firstLine = f.readline()
	#secondLine = f.readline()

	if re.search(r'^finalPoses: ', secondLine):

		parsedLine = re.sub(r'^[0-9]+ landmarks: ', r'', firstLine)
		parsedLine2 = re.sub(r'^finalPoses: ', r'', secondLine)

		landmarksDict = eval(parsedLine)
		controlPoses_G = eval(parsedLine2)

		landmarks_G = []
		for pathID, val in landmarksDict.iteritems():
			controlPose_G = controlPoses_G[pathID]
			currFrame_G = Pose(controlPose_G)
			for entry in val:
				landmarkPoint_L, thresh0, landName = entry
				landmarkPoint_G = (currFrame_G.convertLocalToGlobal(landmarkPoint_L), thresh0, landName)
				landmarks_G.append(landmarkPoint_G)



		hypID = int(firstLine[0])
		#print "hypID =", hypID
		#print parsedLine

		#print "check for outliers:"

		for k in range(len(landmarks_G)):
			landmark_G = landmarks_G[k]
			#print k, landmark_G

		#print "computeEval:", self.hypothesisID

		isOutlier = False
		poseSum = 0.0

		results = []

		for j in range(len(landmarks_G)):

			#print "landmark:", landmarks_G[j]
			p1 = landmarks_G[j][0]
			thresh1 = landmarks_G[j][1]
			p2 = deepcopy(p1)
			thresh2 = deepcopy(landmarks_G[j][1])


			#print "%d %1.2f" % (j, landmarks_G[j][1])

			for k in range(j+1, len(landmarks_G)):

				if j != k:

					p2 = landmarks_G[k][0]
					thresh2 = landmarks_G[k][1]
					dist = sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)


					if dist < LANDMARK_THRESH:
						diffVal = sqrt(dist*dist/(thresh1*thresh1 + thresh2*thresh2))
						poseSum += diffVal
						#print "%d %d %1.2f %1.2f %1.2f %1.2f" % (j, k, dist, diffVal, landmarks_G[j][1], landmarks_G[k][1])
						results.append((j, k, dist, diffVal, landmarks_G[j][1], landmarks_G[k][1], landmarks_G[j][0], landmarks_G[k][0]))

		""" outlier detection protocol

		1) k nearest neighbors between bloom points.  if one greater than some small threshold, fail
		2) k nearest neighbors between all points.  if one greater than some larger threshold, fail


		"""

		allNearestNeighbors = []
		bloomNearestNeighbors = []

		for j in range(len(landmarks_G)):

			p1 = landmarks_G[j][0]
			thresh1 = landmarks_G[j][1]
			landmarkType1 = landmarks_G[j][2]

			thisLandmarksBloomNeighbors = []
			thisLandmarksNeighbors = []

			for k in range(0, len(landmarks_G)):
				if j != k:
					p2 = landmarks_G[k][0]
					thresh2 = landmarks_G[k][1]
					landmarkType2 = landmarks_G[k][2]

					dist = sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)
					diffVal = sqrt(dist*dist/(thresh1*thresh1 + thresh2*thresh2))


					if (landmarkType1 == "archPoint" or landmarkType1 == "bloomPoint") and (landmarkType2 == "archPoint" or landmarkType2 == "bloomPoint"):
						print landmarkType1, landmarkType2
						""" bloom only nearest neighbors """
						thisLandmarksBloomNeighbors.append((j, k, dist, diffVal, landmarks_G[j][1], landmarks_G[k][1], landmarks_G[j][0], landmarks_G[k][0]))

					""" all nearest neighbors """
					thisLandmarksNeighbors.append((j, k, dist, diffVal, landmarks_G[j][1], landmarks_G[k][1], landmarks_G[j][0], landmarks_G[k][0]))

			thisLandmarksBloomNeighbors = sorted(thisLandmarksBloomNeighbors, key=itemgetter(2), reverse=False)
			thisLandmarksNeighbors = sorted(thisLandmarksNeighbors, key=itemgetter(2), reverse=False)
		
			allNearestNeighbors.append(thisLandmarksNeighbors)

			if len(thisLandmarksBloomNeighbors) > 0:
				bloomNearestNeighbors.append(thisLandmarksBloomNeighbors)

		minResults = []

		for j in range(len(landmarks_G)):

			p1 = landmarks_G[j][0]
			thresh1 = landmarks_G[j][1]

			minDist = 1e100
			minDiff = 1e100
			minK = 0
			for k in range(0, len(landmarks_G)):
				if j != k:
					p2 = landmarks_G[k][0]
					thresh2 = landmarks_G[k][1]
					dist = sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

					if thresh1 <= 0.5 and thresh2 <= 0.5:
						diffVal = sqrt(dist*dist/(thresh1*thresh1 + thresh2*thresh2))

					if diffVal < minDiff:
						minDiff = diffVal
						minDist = dist
						minK = k

			minResults.append((j, minK, minDist, minDiff, landmarks_G[j][1], landmarks_G[minK][1], landmarks_G[j][0], landmarks_G[minK][0]))


		results = sorted(results, key=itemgetter(3), reverse=False)
		minResults = sorted(minResults, key=itemgetter(3), reverse=False)




			#if minDist <= LANDMARK_THRESH and minDist > OUTLIER_THRESH:
			#	print "outlier:", minDist, p1, landmarks_G[minK][0]
			#	isOutlier = True
			#else:
			#	print "inlier:", minDist, p1, landmarks_G[minK][0]
		
		if len(minResults) > 0:
			maxDiff = minResults[-1][3]



			#print "isOutlier =", isOutlier
			#print "poseSum =", poseSum


			#print "one-to-one results"
			#for result in results:
			#	print result


			print "max:", maxDiff

			maxDiff = 0.0
			print "minResults:"
			for minResult in minResults:
				print minResult[3]

				if minResult[3] > maxDiff:
					maxDiff = minResult[3]
			print "max:", maxDiff

			isReject = False
			print "bloom neighbors:"
			for result in bloomNearestNeighbors:

				if len(result) > 1:
					print result[0][2], result[1][2]
					if result[0][2] > 0.4 or result[1][2] > 0.4:
						isReject = True
				else:
					print result[0][2]
					if result[0][2] > 0.4:
						isReject = True


				#allNearestNeighbors.append(thisLandmarksNeighbors)
				#bloomNearestNeighbors.append(thisLandmarksBloomNeighbors)



			fig = plt.gcf()
			axes = fig.gca()


			xP = []
			yP = []
			for landmark_G in landmarks_G:

				xP.append(landmark_G[0][0])
				yP.append(landmark_G[0][1])

				circle2=plt.Circle((landmark_G[0][0],landmark_G[0][1]),landmark_G[1],color='b',fill=False)
				axes.add_artist(circle2)

			axes.scatter(xP,yP)

			for result in minResults:
				
				if maxDiff > 0.0:
					normDiff = result[3]/maxDiff
				else:
					normDiff = 1.0

				p1 = result[6]
				p2 = result[7]

				xP = [p1[0], p2[0]]
				yP = [p1[1], p2[1]]

				axes.plot(xP,yP, color='k', alpha=normDiff)



			#axes.set_ylim(-2.0,2.0)
			axes.set_title("hypID = %d, sum = %1.2f, max = %f, REJECT = %d" % (hypID, poseSum, maxDiff, isReject))
			axes.axis("equal")
			#axes.set_xlim(-2.0,4.0)
			fig.savefig("landmarks_%02u_%04u.png" % (hypID, plotCount))

			axes.cla()
			plotCount += 1

















