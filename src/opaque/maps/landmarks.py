#!/usr/bin/python
from Pose import Pose
from shoots import computeGlobalControlPoses




def nodeToGlobalLandmarks(controlPoses, pathIDs, parentHash, nodeLandmarks, pathClasses, exemptNodes = []):


	controlPoses_G = computeGlobalControlPoses(controlPoses, parentHash)

	""" collect landmarks that we can localize against """
	landmarks_G = []
	for pathID in pathIDs:

		nodeSet = pathClasses[pathID]["nodeSet"]
		currFrame_L = Pose(controlPoses[pathID])

		for nodeID in nodeSet:

			landmarkPoint_N = nodeLandmarks[pathID][nodeID]

			currFrame_G = Pose(controlPoses_G[pathID])
			nodePose_L = pathClasses[pathID]["localNodePoses"][nodeID]
			poseFrame_L = Pose(nodePose_L)

			if landmarkPoint_N != None and not nodeID in exemptNodes:
				#nodeID != nodeID0 and nodeID != nodeID1:
				landmarkPoint_L = poseFrame_L.convertLocalToGlobal(landmarkPoint_N)
				landmarkPoint_G = currFrame_G.convertLocalToGlobal(landmarkPoint_L)
				landmarks_G.append(landmarkPoint_G)



	return landmarks_G


