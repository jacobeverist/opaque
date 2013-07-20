#include "RefNode.h"
#include <cmath>

RefNode::RefNode(int nid, int jid, double x, double z, double p, int numJoints, double *joints) {

	refX = x;
	refZ = z;
	refP = p;

	gndX = x;
	gndZ = z;
	gndP = p;

	nodeID = nid;
	jointID = jid;
	this->numJoints = numJoints;

	for (int i = 0 ; i < numJoints ; i++ ) {
		nom_ang[i] = joints[i];
	}

	maxError = 0.0;
	active_time = 0;
	errSum = 0.0;
	currError = 0.0;
	maxErrorReached = 0;
	newRoot = 0;

}

RefNode::~RefNode() {

}

double RefNode::getRefPoseX() {
	return refX;
}

double RefNode::getRefPoseZ() {
	return refZ;
}

double RefNode::getRefPoseP() {
	return refP;
}

int RefNode::getNodeID() {
	return nodeID;
}

int RefNode::getJointID() {
	return jointID;
}

void RefNode::setGroundTruthPose(double x, double z, double p) {
	gndX = x;
	gndZ = z;
	gndP = p;
}

double RefNode::getGroundTruthPoseX() {
	return gndX;
}

double RefNode::getGroundTruthPoseZ() {
	return gndZ;
}

double RefNode::getGroundTruthPoseP() {
	return gndP;
}

double RefNode::computeStabilityError(double *joints) {

	// saddle points: 4,8,12,16,20,24,28,32,36
	int interval = jointID % 4;
	int low_index, high_index;

	// find the saddle point joints to which to check stability
	if (interval == 0) {
		low_index = jointID;
		high_index = jointID + 4;
	}
	else {
		low_index = jointID - interval;
		high_index = low_index + 4;
	}
	
	// set the max and min limits if they go over
	if (low_index < 0) 
		low_index = 0;

	if (high_index >= numJoints)
		high_index = numJoints - 1;

	double error = 0.0;
	for (int i = low_index; i <= high_index ; i++ ) 
		error += fabs(nom_ang[i] - joints[i]);

	if (error > maxError)
		maxError = error;

	currError = error;
	errSum += error;

	if (currError > 0.1)
		maxErrorReached = 1; 

	return error;

}

double RefNode::getStabilityError() {
	return currError;
}

double RefNode::getMaxStabilityError() {
	return maxError;
}

double RefNode::getAvgStabilityError() {
	return errSum/active_time;
}

int RefNode::isMaxErrorReached() {
	return maxErrorReached;
}

void RefNode::updateTime() {
	active_time += 1;
}

void RefNode::setNewRoot(int isRoot) {
	newRoot = isRoot;
}

int RefNode::isNewRoot() {
	return newRoot;
}



