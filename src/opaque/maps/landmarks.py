#!/usr/bin/python
from Pose import Pose
from shoots import computeGlobalControlPoses
from math import sqrt




def nodeToLocalLandmarks(pathID, nodeLandmarks, pathClasses, exemptNodes = []):

	""" collect landmarks that we can localize against """
	landmarks_L = []

	nodeSet = pathClasses[pathID]["nodeSet"]

	for nodeID in nodeSet:

		if nodeLandmarks[pathID][nodeID] != None:
			landmarkPoint_N, threshN, landName = nodeLandmarks[pathID][nodeID]

			nodePose_L = pathClasses[pathID]["localNodePoses"][nodeID]
			poseFrame_L = Pose(nodePose_L)

			if not nodeID in exemptNodes:
				landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
				landmarkPoint_tuple = (landmarkPoint_L, threshN, landName)
				landmarks_L.append(landmarkPoint_tuple)



	return landmarks_L


def nodeToGlobalLandmarks(controlPoses, pathIDs, parentHash, nodeLandmarks, pathClasses, exemptNodes = []):

	#print "nodeLandmarks:", nodeLandmarks

	controlPoses_G = computeGlobalControlPoses(controlPoses, parentHash)

	""" collect landmarks that we can localize against """
	landmarks_G = []
	for pathID in pathIDs:

		nodeSet = pathClasses[pathID]["nodeSet"]
		currFrame_L = Pose(controlPoses[pathID])

		for nodeID in nodeSet:

			if nodeLandmarks[pathID][nodeID] != None:
				landmarkPoint_N, threshN, landName = nodeLandmarks[pathID][nodeID]


				currFrame_G = Pose(controlPoses_G[pathID])
				nodePose_L = pathClasses[pathID]["localNodePoses"][nodeID]
				poseFrame_L = Pose(nodePose_L)

				#print "landmark_N:", nodeID, landmarkPoint_N, nodePose_L

				if not nodeID in exemptNodes:
					#nodeID != nodeID0 and nodeID != nodeID1:
					landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
					#print "landmark_L:", nodeID, landmarkPoint_L, controlPoses_G[pathID]
					landmarkPoint_G = (currFrame_G.convertLocalToGlobal(landmarkPoint_L), threshN, landName)
					#print "landmark_G:", nodeID, landmarkPoint_G[0]
					landmarks_G.append(landmarkPoint_G)



	return landmarks_G


def getNodeLandmark(nodeID, poseData):

	BEND_THRESH = 0.9
	ARCH_THRESH = 0.3
	BLOOM_THRESH = 0.3

	targetNodeLandmark_N = None
	spatialFeature = poseData.spatialFeatures[nodeID][0]

	if spatialFeature["inflectionPoint"] != None:
		landmarkPoint_N = spatialFeature["inflectionPoint"]
		targetNodeLandmark_N = (landmarkPoint_N, BEND_THRESH, "bendPoint")
	elif spatialFeature["bloomPoint"] != None:
		landmarkPoint_N = spatialFeature["bloomPoint"]
		targetNodeLandmark_N = (landmarkPoint_N, BLOOM_THRESH, "bloomPoint")
	elif spatialFeature["archPoint"] != None:
		landmarkPoint_N = spatialFeature["archPoint"]
		targetNodeLandmark_N = (landmarkPoint_N, ARCH_THRESH, "archPoint")



	return targetNodeLandmark_N


def computeConsistency(landmarks_G):

	LANDMARK_THRESH = 7.0
	poseSum = 0.0


	for i in range(len(landmarks_G)):
		p1 = landmarks_G[i][0]
		thresh1 = landmarks_G[i][1]

		for j in range(i+1, len(landmarks_G)):
			p2 = landmarks_G[j][0]
			thresh2 = landmarks_G[j][1]

			maxThresh = thresh1
			if thresh2 > thresh1:
				maxThresh = thresh2

			dist0 = sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


			if dist0 < LANDMARK_THRESH:
				""" square of the distance """
				#poseSum += sqrt(dist0*dist0/(maxThresh*maxThresh))

				poseSum += sqrt(dist0*dist0/(thresh1*thresh1 + thresh2*thresh2))
				
				#if dist0 > maxThresh:
				#	poseSum += 10.0*dist0*dist0
				#else:
				#	poseSum += dist0*dist0


	return poseSum
