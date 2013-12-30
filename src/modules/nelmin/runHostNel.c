// Includes
#include <stdio.h>
#include <math.h>
#include "runNelmin.h"

double objFunc(double *params, double *d_matchPairs, double *d_offset, double *d_sum, int N, double *pose1, double *poses_2, int numPoses, double angNom, double angLim, double uHigh, double uLow, double *resultOffset);
double objFuncIndAngle(double *params, double *d_matchPairs, double *d_offset, double *d_sum, int N, double *pose1, double *poses_2, int numPoses, double angNom, double angLim, double uHigh, double uLow, double *resultOffset);

#include "nelmin.h"

//#include "data.h"

#ifndef M_PI
#define M_PI 3.141592654
#endif

#define threadsPerBlock 256

void icpError(const double* A, const double* offset, double *sumVal, int N);

__inline double normalizeAngle(double angle);
__inline void dotProduct(double *R, double *vec, double *result);
__inline void buildMatrix(double *R, double ang, int isForward);
__inline void transposeMatrix(double *R, double *Rt);


void icpError(const double* A, const double* offset, double *sumVal, int N) {

	int tid = 0;
	double xd, yd, theta = 0.0;
	double ax = 0.0;
	double ay = 0.0;
	double bx = 0.0;
	double by1 = 0.0;
	double tx, ty, dx, dy = 0.0;
	double r11, r12, r21, r22 = 0.0;
	double c11, c12, c21, c22 = 0.0;
	double b11, b12, b21, b22 = 0.0;
	double res11, res12, res21, res22, resDet = 0.0;
	double q11, q12, q21, q22 = 0.0;
	double errVal = 0.0;

	double tempVal = 0.0;

	while (tid < N) {

		xd = offset[0];
		yd = offset[1];
		theta = offset[2];

		ax = A[tid*12];
		ay = A[tid*12 + 1];
		c11 = A[tid*12 + 2];
		c12 = A[tid*12 + 3];
		c21 = A[tid*12 + 4];
		c22 = A[tid*12 + 5];
		bx = A[tid*12 + 6];
		by1 = A[tid*12 + 7];
		b11 = A[tid*12 + 8];
		b12 = A[tid*12 + 9];
		b21 = A[tid*12 + 10];
		b22 = A[tid*12 + 11];
	
		tx = ax*cos(theta) - ay*sin(theta) + xd;
		ty = ax*sin(theta) + ay*cos(theta) + yd;
		dx = bx - tx;
		dy = by1 - ty;
	
		r11 = cos(theta);
		r12 = -sin(theta);
		r21 = sin(theta);
		r22 = r11;
	
		res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12);
		res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22);
		res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12);
		res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22);
		
		resDet = res22*res11 - res12*res21;
		
		q11 = res22/resDet;
		q12 = -res12/resDet;
		q21 = -res21/resDet;
		q22 = res11/resDet;
	
		errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22);
	
		tempVal += errVal;

		tid += 1;
	}

	sumVal[0] = tempVal;
}

__inline double normalizeAngle(double angle) {

	while (angle>M_PI)
		angle=angle-2*M_PI;

	while (angle<=-M_PI)
		angle=angle+2*M_PI;

	return angle;
}

__inline void dotProduct(double *R, double *vec, double *result) {
	result[0] = R[0]*vec[0] + R[1]*vec[1];
	result[1] = R[2]*vec[0] + R[3]*vec[1];
}

__inline void transposeMatrix(double *R, double *Rt) {
	Rt[0] = R[0];
	Rt[3] = R[3];

	Rt[1] = R[2];
	Rt[2] = R[1];

}

__inline void buildMatrix(double *R, double ang, int isForward) {

	if ( isForward ) {
		R[0] = cos(ang); 
		R[1] = -sin(ang); 
		R[2] = sin(ang); 
		R[3] = cos(ang);
	}
	else {
		R[0] = cos(ang); 
		R[1] = sin(ang); 
		R[2] = -sin(ang); 
		R[3] = cos(ang);
	}

}

double objFunc(double *params, double *d_matchPairs, double *d_offset, double *d_sum, int N, double *pose1, double *poses_2, int numPoses, double angNom, double angLim, double uHigh, double uLow, double *resultOffset) {


	double currU = params[0];
	double currAng = params[1];

	// double pose1[3] = {0.0,0.0,0.0};
	double pose2[3] = {0.0,0.0,0.0};
	double point1[2];
	double point2[2];
	double ang1;
	double ang2;
	double ang2_off;
	double overlapPoint1Pose[3];
	double overlapPoint2Pose[3];
	double estPose1[3] = {0.0,0.0,0.0};
	double dist1;
	double vecAng1 = 0.0;
	double backR1[4];
	double foreR1[4];
	double R1[4];
	double estPose2[3];
	double dist2;
	double vecAng2 = 0.0;
	double backR2[4];
	double foreR2[4];
	double R2[4];
	double negCurve2Pose[3] = {0.0,0.0,0.0};
	double offsetR2[4];
	double finalVec2[2] = {0.0,0.0};
	double transVec2[2] = {0.0,0.0};
	double resPose2[3] = {0.0,0.0,0.0};
	double finalVec[2] = {0.0,0.0};
	double R2t[4];
	double tempVec1[2] = {0.0,0.0};
	double resVec[2] = {0.0,0.0};
	double tempVec2[2] = {0.0,0.0};
	double globalVec[2] = {0.0,0.0};
	double localVec[2] = {0.0,0.0};
	double localOffset[3] = {0.0,0.0,0.0};
	double h_offset[3] = {0.0,0.0,0.0};
	double totalSum;

	double tempVal1;

	//printf("%f %f :", currU, currAng);

	if (currU < 0.0) {
		//printf("%f\n",100000.0 * (0.15-currU));
		return 10000000.0 * (0.15-currU);
	}
	if (currU > 1.0 ) {
		//printf("%f\n",100000.0 * (currU-0.85));
		return 10000000.0 * (currU - 0.85);
	}

	if (currU < uLow) {
		//printf("%f\n",100000.0 * (uLow+0.05-currU));
		return 10000000.0 * (uLow+0.05-currU);
	}
	if (currU > uHigh) {
		//printf("%f\n",100000.0 * (currU-uHigh+0.05));
		return 10000000.0 * (currU - uHigh+0.05);
	}

	if (fabs(currAng-angNom) > angLim) {
		return 10000000.0 * (fabs(currAng-angNom)-angLim+0.05);
	}

	//printf("u1 = %lf\n", u1);

	// compute the point and angles from the parameters 
	/* 
	if ( u1 >= 1.0 ) {
		double *ptr = &(poses_1[3*(numPoses-1)]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
		//printf("foo %d %d\n", numPoses, 3*(numPoses-1));
	}
	else if (u1 < 0.0) {
		double *ptr = &(poses_1[0]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_1[3*(int)(u1*(numPoses-1))]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	*/

	if (currU >= 1.0) {
		double *ptr = &(poses_2[3*(numPoses-1)]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}
	else if (currU < 0.0) {
		double *ptr = &(poses_2[0]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_2[3*(int)(currU*(numPoses-1))]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}

	//printf("pose1 = %lf %lf %lf\n", pose1[0], pose1[1], pose1[2]);
	//printf("pose2 = %lf %lf %lf\n", pose2[0], pose2[1], pose2[2]);


	//double point1[2];
	point1[0] = pose1[0];
	point1[1] = pose1[1];

	//double point2[2];
	point2[0] = pose2[0];
	point2[1] = pose2[1];

	//double ang1 = pose1[2];
	//double ang2 = pose2[2];
	//double ang2_off = ang2 + currAng;
	ang1 = pose1[2];
	ang2 = pose2[2];
	ang2_off = ang2 + currAng;

	//offset = computeOffset(point1, point2, ang1, ang2 + currAng)

	//" corner points and orientations "
	//double overlapPoint1Pose[3];
       	overlapPoint1Pose[0] = point1[0];
       	overlapPoint1Pose[1] = point1[1];
       	overlapPoint1Pose[2] = ang1;

	//double overlapPoint2Pose[3];
       	overlapPoint2Pose[0] = point2[0];
       	overlapPoint2Pose[1] = point2[1];
       	overlapPoint2Pose[2] = ang2_off;

	//printf("overlapPoint1Pose = %f %f %f\n", overlapPoint1Pose[0], overlapPoint1Pose[1], overlapPoint1Pose[2]);
	//printf("overlapPoint2Pose = %f %f %f\n", overlapPoint2Pose[0], overlapPoint2Pose[1], overlapPoint2Pose[2]);


	/// CHECK

	//" convert the desired intersection point on curve 1 into global coordinates "
	//poseProfile1 = Pose([0.0,0.0,0.0])
	
	//double estPose1[3] = {0.0,0.0,0.0};
	//double dist1 = sqrt(estPose1[0]*estPose1[0] + estPose1[1]*estPose1[1]);
	//double vecAng1 = 0.0;
	dist1 = sqrt(estPose1[0]*estPose1[0] + estPose1[1]*estPose1[1]);
	vecAng1 = 0.0;

	//printf("estPose1 = %lf %lf %lf\n", estPose1[0], estPose1[1], estPose1[2]);

	//printf("dist1 = %lf\n", dist1);

	if (dist1 > 0.0) {
		vecAng1 = acos(estPose1[0]/dist1);
		if (asin(estPose1[1]/dist1) < 0) {
			vecAng1 = -vecAng1;
		}
	}
	else {
		vecAng1 = 0.0;
	}
	
	//printf("vecAng1 = %lf\n", vecAng1);

	//double backR1[4];
	buildMatrix(backR1, vecAng1, 0);

	//double foreR1[4];
	buildMatrix(foreR1, vecAng1, 1);

	//double R1[4];
	buildMatrix(R1, estPose1[2], 0);

	// now convert this point into a pose, and perform the inverse transform using corner2Pose 

	//desGlobalPose2 = Pose(overlapPoint1Pose)
	//double estPose2[3];
	estPose2[0] = overlapPoint1Pose[0];
	estPose2[1] = overlapPoint1Pose[1];
	estPose2[2] = overlapPoint1Pose[2];

	//double dist2 = sqrt(estPose2[0]*estPose2[0] + estPose2[1]*estPose2[1]);
	//double vecAng2 = 0.0;
	dist2 = sqrt(estPose2[0]*estPose2[0] + estPose2[1]*estPose2[1]);
	
	//printf("estPose2 = %f %f %f\n", estPose2[0], estPose2[1], estPose2[2]);

	if (dist2 > 0.0) {
		vecAng2 = acos(estPose2[0]/dist2);
		if (asin(estPose2[1]/dist2) < 0) {
			vecAng2 = -vecAng2;
		}
	}
	else {
		vecAng2 = 0.0;
	}

	//double backR2[4];
	buildMatrix(backR2, vecAng2, 0);

	//double foreR2[4];
	buildMatrix(foreR2, vecAng2, 1);

	//double R2[4];
	buildMatrix(R2, estPose2[2], 0);

	// CHECK2

	//" perform inverse offset from the destination pose "
	//negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
	//double negCurve2Pose[3] = {0.0,0.0,0.0};

	tempVal1 = overlapPoint2Pose[2];
	//printf("overlapPoint2Pose[2] = %lf\n", overlapPoint2Pose[2]);

	//double offsetR2[4];
	buildMatrix(offsetR2, overlapPoint2Pose[2], 0);

	//printf("offsetR2 = %f %f %f %f\n", cos(tempVal1), sin(tempVal1), -sin(tempVal1), cos(tempVal1));
	//printf("offsetR2 = %f %f %f %f\n", offsetR2[0], offsetR2[1], offsetR2[2], offsetR2[3]);

	//double finalVec2[2] = {0.0,0.0};
	finalVec2[0] = overlapPoint2Pose[0];
	finalVec2[1] = overlapPoint2Pose[1];

	//printf("finalVec2 = %f %f\n", finalVec2[0], finalVec2[1]);

	// transVec = dot(offsetR2, finalVec2)
	//double transVec2[2] = {0.0,0.0};
	//transVec2[0] = offsetR2[0]*finalVec2[0] + offsetR2[1]*finalVec2[1];
	//transVec2[1] = offsetR2[2]*finalVec2[0] + offsetR2[3]*finalVec2[1];
	dotProduct(offsetR2, finalVec2, transVec2);

	negCurve2Pose[0] = -transVec2[0];
	negCurve2Pose[1] = -transVec2[1];
	negCurve2Pose[2] = -overlapPoint2Pose[2];
	

	//printf("negCurve2Pose = %f %f %f\n", negCurve2Pose[0], negCurve2Pose[1], negCurve2Pose[2]);

	// CHECK3
	
	//" relative pose between pose 1 and pose 2 to make corners coincide and same angle "
	//resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
	//double resPose2[3] = {0.0,0.0,0.0};

	// pss
	//double finalVec[2] = {0.0,0.0};
	finalVec[0] = negCurve2Pose[0];
	finalVec[1] = negCurve2Pose[1];

	//finalVec = array([[negCurve2Pose[0]], [negCurve2Pose[1]]])

	//transVec = dot(transpose(self.R), finalVec)
	//double R2t[4];
	transposeMatrix(R2, R2t);

	//double tempVec1[2] = {0.0,0.0};
	dotProduct(R2t, finalVec, tempVec1);

	//double tempVec1[2] = {0.0,0.0};
	//dotProduct(R2, finalVec, tempVec1);

	//resVec = dot(self.backR, transVec)
	//double resVec[2] = {0.0,0.0};
	dotProduct(backR2, tempVec1, resVec);

	//resVec[0, 0] += self.dist
	resVec[0] += dist2;

	//tempVec = dot(self.foreR, resVec)
	//double tempVec2[2] = {0.0,0.0};
	dotProduct(foreR2, resVec, tempVec2);

	resPose2[0] = tempVec2[0];
	resPose2[1] = tempVec2[1];
	resPose2[2] = normalizeAngle(estPose2[2] + negCurve2Pose[2]);

	//printf("resPose2 = %f %f %f\n", resPose2[0], resPose2[1], resPose2[2]);
	

	//localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)

	//def convertGlobalPoseToLocal(self, pose):

	//" transform pnt to local coordinates"
	//double globalVec[2] = {0.0,0.0};
	globalVec[0] = resPose2[0];
	globalVec[1] = resPose2[1];

	//" perform translation correction "
	//tempVec = dot(self.backR, globalVec)
	dotProduct(backR1, globalVec, tempVec1);
	tempVec1[0] -= dist1;

	//transVec = dot(self.foreR, tempVec)
	dotProduct(foreR1, tempVec1, tempVec2);

	//" now, apply rotation correction with respect to origin "
	//localVec = dot(self.R, transVec)
	//double localVec[2] = {0.0,0.0};
	dotProduct(R1, tempVec2, localVec);

	//localPose = [localVec[0,0], localVec[1,0], normalizeAngle(resPose2[2] - self.estPose[2])]
	//return localPose
	//double localOffset[3] = {0.0,0.0,0.0};
	localOffset[0] = localVec[0];
	localOffset[1] = localVec[1];
	localOffset[2] = normalizeAngle(resPose2[2]-estPose1[2]);

	//double h_offset[3] = {0.0,0.0,0.0};
	h_offset[0] = localOffset[0];
	h_offset[1] = localOffset[1];
	h_offset[2] = localOffset[2];


	//printf("offset = %f %f %f\n", h_offset[0], h_offset[1], h_offset[2]);
	
	resultOffset[0] = h_offset[0];
	resultOffset[1] = h_offset[1];
	resultOffset[2] = h_offset[2];

	d_offset[0] = h_offset[0];
	d_offset[1] = h_offset[1];
	d_offset[2] = h_offset[2];


	icpError(d_matchPairs, d_offset, d_sum, N);

	// Copy result from device memory to host memory
	// h_sum contains the result in host memory

	//double totalSum = d_sum[0];
	totalSum = d_sum[0];

	//printf("totalSum = %lf\n");

	return totalSum;

}

/*
double objFunc2(double *params, double *d_matchPairs, double *d_offset, double *d_sum, int N, double *poses_1, double *poses_2, int numPoses, double uHigh, double uLow, double u1, double *resultOffset) {

	// def medialOverlapCostFunc_GPU(params, input_GPU, N, medialSpline1, medialSpline2, uHigh, uLow, u1):

	double currU = params[0];
	double currAng = params[1];

	double pose1[3] = {0.0,0.0,0.0};
	double pose2[3] = {0.0,0.0,0.0};

	//printf("%f %f :", currU, currAng);

	if (currU < 0.0) {
		//printf("%f\n",100000.0 * (0.15-currU));
		return 100000.0 * (0.15-currU);
	}
	if (currU > 1.0 ) {
		//printf("%f\n",100000.0 * (currU-0.85));
		return 100000.0 * (currU - 0.85);
	}

	if (currU < uLow) {
		//printf("%f\n",100000.0 * (uLow+0.05-currU));
		return 100000.0 * (uLow+0.05-currU);
	}
	if (currU > uHigh) {
		//printf("%f\n",100000.0 * (currU-uHigh+0.05));
		return 100000.0 * (currU - uHigh+0.05);
	}

	// compute the point and angles from the parameters 
	if ( u1 >= 1.0 ) {
		double *ptr = &(poses_1[3*numPoses-1]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else if (u1 < 0.0) {
		double *ptr = &(poses_1[0]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_1[3*(int)(u1*numPoses)]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}

	if (currU >= 1.0) {
		double *ptr = &(poses_2[3*numPoses-1]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}
	else if (currU < 0.0) {
		double *ptr = &(poses_2[0]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_2[3*(int)(currU*numPoses)]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}

	//resultOffset[0] = currU;
	//resultOffset[1] = currAng;
	//resultOffset[2] = u1;


	double point1[2];
	point1[0] = pose1[0];
	point1[1] = pose1[1];

	double point2[2];
	point2[0] = pose2[0];
	point2[1] = pose2[1];

	double ang1 = pose1[2];
	double ang2 = pose2[2];
	double ang2_off = ang2 + currAng;

	//offset = computeOffset(point1, point2, ang1, ang2 + currAng)

	//" corner points and orientations "
	double overlapPoint1Pose[3];
       	overlapPoint1Pose[0] = point1[0];
       	overlapPoint1Pose[1] = point1[1];
       	overlapPoint1Pose[2] = ang1;

	double overlapPoint2Pose[3];
       	overlapPoint2Pose[0] = point2[0];
       	overlapPoint2Pose[1] = point2[1];
       	overlapPoint2Pose[2] = ang2_off;

	//resultOffset[0] = overlapPoint2Pose[0];
	//resultOffset[1] = overlapPoint2Pose[1];
	//resultOffset[2] = overlapPoint2Pose[2];

	/// CHECK

	//" convert the desired intersection point on curve 1 into global coordinates "
	//poseProfile1 = Pose([0.0,0.0,0.0])
	
	double estPose1[3] = {0.0,0.0,0.0};
	double dist1 = sqrt(estPose1[0]*estPose1[0] + estPose1[1]*estPose1[1]);
	double vecAng1 = 0.0;

	if (dist1 > 0.0) {
		vecAng1 = acos(estPose1[0]/dist1);
		if (asin(estPose1[1]/dist1) < 0) {
			vecAng1 = -vecAng1;
		}
	}
	else {
		vecAng1 = 0.0;
	}
	
	double backR1[4];
	buildMatrix(backR1, vecAng1, 0);

	double foreR1[4];
	buildMatrix(foreR1, vecAng1, 1);

	double R1[4];
	buildMatrix(R1, estPose1[2], 0);

	// now convert this point into a pose, and perform the inverse transform using corner2Pose 

	//desGlobalPose2 = Pose(overlapPoint1Pose)
	double estPose2[3];
       	estPose2[0] = overlapPoint1Pose[0];
       	estPose2[1] = overlapPoint1Pose[1];
       	estPose2[2] = overlapPoint1Pose[2];
	double dist2 = sqrt(estPose2[0]*estPose2[0] + estPose2[1]*estPose2[1]);
	double vecAng2 = 0.0;

	if (dist2 > 0.0) {
		vecAng2 = acos(estPose2[0]/dist2);
		if (asin(estPose2[1]/dist2) < 0) {
			vecAng2 = -vecAng2;
		}
	}
	else {
		vecAng2 = 0.0;
	}

	double backR2[4];
	buildMatrix(backR2, vecAng2, 0);

	double foreR2[4];
	buildMatrix(foreR2, vecAng2, 1);

	double R2[4];
	buildMatrix(R2, estPose2[2], 0);

	// CHECK2

	//" perform inverse offset from the destination pose "
	//negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
	double negCurve2Pose[3] = {0.0,0.0,0.0};

	double offsetR2[4];
	buildMatrix(offsetR2, overlapPoint2Pose[2], 0);

	double finalVec2[2] = {0.0,0.0};
	finalVec2[0] = overlapPoint2Pose[0];
	finalVec2[1] = overlapPoint2Pose[1];

	// transVec = dot(offsetR2, finalVec2)
	double transVec2[2] = {0.0,0.0};
	//transVec2[0] = offsetR2[0]*finalVec2[0] + offsetR2[1]*finalVec2[1];
	//transVec2[1] = offsetR2[2]*finalVec2[0] + offsetR2[3]*finalVec2[1];
	dotProduct(offsetR2, finalVec2, transVec2);

	negCurve2Pose[0] = -transVec2[0];
	negCurve2Pose[1] = -transVec2[1];
	negCurve2Pose[2] = -overlapPoint2Pose[2];

	//resultOffset[0] = negCurve2Pose[0];
	//resultOffset[1] = negCurve2Pose[1];
	//resultOffset[2] = negCurve2Pose[2];
	
	// CHECK3
	
	//" relative pose between pose 1 and pose 2 to make corners coincide and same angle "
	//resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
	double resPose2[3] = {0.0,0.0,0.0};

	// pss
	double finalVec[2] = {0.0,0.0};
	finalVec[0] = negCurve2Pose[0];
	finalVec[1] = negCurve2Pose[1];

	//finalVec = array([[negCurve2Pose[0]], [negCurve2Pose[1]]])

	//void dotProduct(double *R, double *vec, double *result) {
	//result[0] = R[0]*vec[0] + R[1]*vec[1];
	//result[1] = R[2]*vec[0] + R[3]*vec[1];

	//transVec = dot(transpose(self.R), finalVec)
	double R2t[4];
	transposeMatrix(R2, R2t);

	double tempVec1[2] = {0.0,0.0};
	dotProduct(R2t, finalVec, tempVec1);

	//resVec = dot(self.backR, transVec)
	double resVec[2] = {0.0,0.0};
	dotProduct(backR2, tempVec1, resVec);

	//resVec[0, 0] += self.dist
	resVec[0] += dist2;

	//tempVec = dot(self.foreR, resVec)
	double tempVec2[2] = {0.0,0.0};
	dotProduct(foreR2, resVec, tempVec2);

	resPose2[0] = tempVec2[0];
	resPose2[1] = tempVec2[1];
	resPose2[2] = normalizeAngle(estPose2[2] + negCurve2Pose[2]);

	resultOffset[0] = resPose2[0];
	resultOffset[1] = resPose2[1];
	resultOffset[2] = resPose2[2];

	//localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)

	//def convertGlobalPoseToLocal(self, pose):

	//" transform pnt to local coordinates"
	double globalVec[2] = {0.0,0.0};
	globalVec[0] = resPose2[0];
	globalVec[1] = resPose2[1];

	//" perform translation correction "
	//tempVec = dot(self.backR, globalVec)
	dotProduct(backR1, globalVec, tempVec1);
	tempVec1[0] -= dist1;

	//transVec = dot(self.foreR, tempVec)
	dotProduct(foreR1, tempVec1, tempVec2);

	//" now, apply rotation correction with respect to origin "
	//localVec = dot(self.R, transVec)
	double localVec[2] = {0.0,0.0};
	dotProduct(R1, tempVec2, localVec);

	//localPose = [localVec[0,0], localVec[1,0], normalizeAngle(resPose2[2] - self.estPose[2])]
	//return localPose
	double localOffset[3] = {0.0,0.0,0.0};
	localOffset[0] = localVec[0];
	localOffset[1] = localVec[1];
	localOffset[2] = normalizeAngle(resPose2[2]-estPose1[2]);

	double h_offset[3] = {0.0,0.0,0.0};
	h_offset[0] = localOffset[0];
	h_offset[1] = localOffset[1];
	h_offset[2] = localOffset[2];

	//resultOffset[0] = h_offset[0];
	//resultOffset[1] = h_offset[1];
	//resultOffset[2] = h_offset[2];

	//printf("offset = %f %f %f\n", h_offset[0], h_offset[1], h_offset[2]);
	
	d_offset[0] = h_offset[0];
	d_offset[1] = h_offset[1];
	d_offset[2] = h_offset[2];


	icpError(d_matchPairs, d_offset, d_sum, N);

	// Copy result from device memory to host memory
	// h_sum contains the result in host memory

	double totalSum = d_sum[0];

	return totalSum;

}
*/

void doPointCost(double *h_matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *pose1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset) {

	double* d_matchPairs;
	double* d_offset;
	double* d_sum;
	double* h_sum;
	int nParam = 1;
	int icount, ifault, numres;
	double reqmin = 1.0E-18;
	int konvge = 10;
	int kcount = 5000;
	double start[2] = {0.0,0.0};
	double step[2] = {1.5, 0.3};
	double xmin[2] = {0.0,0.0};
	double ynewlo = 0.0;

	size_t size;
	size_t size2;
	size_t size3;

	int i, j;

	double u1, u2, currU, currAng;
	double offset[3];
	double newSum;

	size = 3 * sizeof(double);
	size2 = numPairs * 12 * sizeof(double);
	size3 = sizeof(double);

	h_sum = (double*)malloc(size3);

	u1 = initGuess[0];
	u2 = initGuess[1];
	currU = u2;
	currAng = initGuess[2];

	start[0] = currU;
	start[1] = currAng;

	d_offset = (double*)malloc(size);
	d_matchPairs = (double*)malloc(size2);
	d_sum = (double*)malloc(size3);

	// Copy vectors from host memory to device memory

	for ( i = 0 ; i < 12*numPairs ; i++ ) {
		d_matchPairs[i] = h_matchPairs[i];
	}

	newSum = objFuncIndAngle(&start[0], d_matchPairs, d_offset, d_sum, numPairs, pose1, poses_2, numPoses, angNom, angLim, uHigh, uLow, offset);

	resultSum[0] = newSum;
	resultParam[0] = currU;
	resultParam[1] = currAng;

	resultOffset[0] = offset[0];
	resultOffset[1] = offset[1];
	resultOffset[2] = offset[2];

	// Free device memory
	if (d_matchPairs)
		free(d_matchPairs);
	if (d_offset)
		free(d_offset);
	if (d_sum)
		free(d_sum);

	// Free host memory
	//if (h_matchPairs)
	//	free(h_matchPairs);
	if (h_sum)
		free(h_sum);


}


void doICPIndAngle(double *h_matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *pose1, double *poses_2, int numPoses, double *resultParam, double *resultSum) {

	double* d_matchPairs;
	double* d_offset;
	double* d_sum;
	//double* h_matchPairs;
	double* h_sum;
	//int nParam = 2;
	int nParam = 1;
	int icount, ifault, numres;
	double reqmin = 1.0E-18;
	int konvge = 10;
	int kcount = 5000;
	double start[2] = {0.0,0.0};
	double step[2] = {1.5, 0.3};
	double xmin[2] = {0.0,0.0};
	double ynewlo = 0.0;

	size_t size;
	size_t size2;
	size_t size3;

	int i, j;

	//double u1, u2, uHigh, uLow, currU, currAng;
	double u1, u2, currU, currAng;

	//double pose1[3] = {0.0,0.0,0.0};

	//int blocksPerGrid = (numPairs + threadsPerBlock - 1) / threadsPerBlock;

	size = 3 * sizeof(double);
	size2 = numPairs * 12 * sizeof(double);
	size3 = sizeof(double);

	//h_matchPairs = (double*)malloc(size2);
	h_sum = (double*)malloc(size3);


	//for ( i = 0 ; i < 12 * numPairs ; i++ ) {
	//	h_matchPairs[i] = matchPairs[i];
	//}


	u1 = initGuess[0];
	u2 = initGuess[1];
	//uHigh = u2 + 0.2;
	//uLow = u2 - 0.2;
	currU = u2;
	currAng = initGuess[2];

	start[0] = currU;
	start[1] = currAng;

	d_offset = (double*)malloc(size);
	d_matchPairs = (double*)malloc(size2);
	d_sum = (double*)malloc(size3);

	// Copy vectors from host memory to device memory

	for ( i = 0 ; i < 12*numPairs ; i++ ) {
		d_matchPairs[i] = h_matchPairs[i];
	}
	/*
	for ( i = 0 ; i < numPairs ; i++ ) {
		for ( j = 0 ; j < 12 ; j ++ ) {
			printf("%lf ", d_matchPairs[12*i+j]);
		}
		printf("\n");
	}
	*/

	/*
	if ( u1 >= 1.0 ) {
		double *ptr = &(poses_1[3*(numPoses-1)]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
		//printf("foo %d %d\n", numPoses, 3*(numPoses-1));
	}
	else if (u1 < 0.0) {
		double *ptr = &(poses_1[0]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_1[3*(int)(u1*(numPoses-1))]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	*/

	/*
	nelmin(nParam, &start[0], &xmin[0], &ynewlo, reqmin, &step[0],
		konvge, kcount, &icount, &numres, &ifault,
		numPoses, numPairs, d_matchPairs, d_offset, d_sum, poses_1, poses_2,
		uHigh, uLow, u1);
	*/

	nelminIndAngle(nParam, &start[0], &xmin[0], &ynewlo, reqmin, &step[0],
		konvge, kcount, &icount, &numres, &ifault,
		numPoses, numPairs, d_matchPairs, d_offset, d_sum, pose1, poses_2,
		angNom, angLim, uHigh, uLow);

	resultSum[0] = ynewlo;
	resultParam[0] = xmin[0];
	resultParam[1] = xmin[1];

	//printf("sum = %f\n", ynewlo);	
	//printf("minX = %f %f\n", xmin[0], xmin[1]);
	//printf("icount = %d\n", icount);

	// Free device memory
	if (d_matchPairs)
		free(d_matchPairs);
	if (d_offset)
		free(d_offset);
	if (d_sum)
		free(d_sum);

	// Free host memory
	//if (h_matchPairs)
	//	free(h_matchPairs);
	if (h_sum)
		free(h_sum);

}


//void doTest(double *h_matchPairs, int numPairs, double *initGuess, double uHigh, double uLow, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum) {
void doTest(double *h_matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *pose1, double *poses_2, int numPoses, double *resultParam, double *resultSum) {

	double* d_matchPairs;
	double* d_offset;
	double* d_sum;
	//double* h_matchPairs;
	double* h_sum;
	int nParam = 2;
	int icount, ifault, numres;
	double reqmin = 1.0E-18;
	int konvge = 10;
	int kcount = 5000;
	double start[2] = {0.0,0.0};
	double step[2] = {1.5, 0.3};
	double xmin[2] = {0.0,0.0};
	double ynewlo = 0.0;

	size_t size;
	size_t size2;
	size_t size3;

	int i, j;

	//double u1, u2, uHigh, uLow, currU, currAng;
	double u1, u2, currU, currAng;

	//double pose1[3] = {0.0,0.0,0.0};

	//int blocksPerGrid = (numPairs + threadsPerBlock - 1) / threadsPerBlock;

	size = 3 * sizeof(double);
	size2 = numPairs * 12 * sizeof(double);
	size3 = sizeof(double);

	//h_matchPairs = (double*)malloc(size2);
	h_sum = (double*)malloc(size3);


	//for ( i = 0 ; i < 12 * numPairs ; i++ ) {
	//	h_matchPairs[i] = matchPairs[i];
	//}


	u1 = initGuess[0];
	u2 = initGuess[1];
	//uHigh = u2 + 0.2;
	//uLow = u2 - 0.2;
	currU = u2;
	currAng = initGuess[2];

	start[0] = currU;
	start[1] = currAng;

	d_offset = (double*)malloc(size);
	d_matchPairs = (double*)malloc(size2);
	d_sum = (double*)malloc(size3);

	// Copy vectors from host memory to device memory

	for ( i = 0 ; i < 12*numPairs ; i++ ) {
		d_matchPairs[i] = h_matchPairs[i];
	}
	/*
	for ( i = 0 ; i < numPairs ; i++ ) {
		for ( j = 0 ; j < 12 ; j ++ ) {
			printf("%lf ", d_matchPairs[12*i+j]);
		}
		printf("\n");
	}
	*/

	/*
	if ( u1 >= 1.0 ) {
		double *ptr = &(poses_1[3*(numPoses-1)]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
		//printf("foo %d %d\n", numPoses, 3*(numPoses-1));
	}
	else if (u1 < 0.0) {
		double *ptr = &(poses_1[0]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_1[3*(int)(u1*(numPoses-1))]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	*/

	/*
	nelmin(nParam, &start[0], &xmin[0], &ynewlo, reqmin, &step[0],
		konvge, kcount, &icount, &numres, &ifault,
		numPoses, numPairs, d_matchPairs, d_offset, d_sum, poses_1, poses_2,
		uHigh, uLow, u1);
	*/

	nelmin(nParam, &start[0], &xmin[0], &ynewlo, reqmin, &step[0],
		konvge, kcount, &icount, &numres, &ifault,
		numPoses, numPairs, d_matchPairs, d_offset, d_sum, pose1, poses_2,
		angNom, angLim, uHigh, uLow);

	resultSum[0] = ynewlo;
	resultParam[0] = xmin[0];
	resultParam[1] = xmin[1];

	//printf("sum = %f\n", ynewlo);	
	//printf("minX = %f %f\n", xmin[0], xmin[1]);
	//printf("icount = %d\n", icount);

	// Free device memory
	if (d_matchPairs)
		free(d_matchPairs);
	if (d_offset)
		free(d_offset);
	if (d_sum)
		free(d_sum);

	// Free host memory
	//if (h_matchPairs)
	//	free(h_matchPairs);
	if (h_sum)
		free(h_sum);

}


//void doCost(double *h_matchPairs, int numPairs, double *initGuess, double uHigh, double uLow, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset) {
void doCost(double *h_matchPairs, int numPairs, double *initGuess, double angNom, double angLim, double uHigh, double uLow, double *pose1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset) {

	double* d_matchPairs;
	double* d_offset;
	double* d_sum;
	//double* h_matchPairs;
	double* h_sum;
	int nParam = 2;
	int icount, ifault, numres;
	double reqmin = 1.0E-18;
	int konvge = 10;
	int kcount = 500;
	double start[2] = {0.0,0.0};
	double step[2] = {0.1, 0.1};
	double xmin[2] = {0.0,0.0};
	double ynewlo = 0.0;

	size_t size;
	size_t size2;
	size_t size3;

	int i;

	//double u1, u2, uHigh, uLow, currU, currAng;
	double u1, u2, currU, currAng;
	double offset[3];
	double newSum;

	//double pose1[3] = {0.0,0.0,0.0};

	size = 3 * sizeof(double);
	size2 = numPairs * 12 * sizeof(double);
	size3 = sizeof(double);

	//h_matchPairs = (double*)malloc(size2);
	h_sum = (double*)malloc(size3);


	//for ( i = 0 ; i < 12 * numPairs ; i++ ) {
	//	h_matchPairs[i] = matchPairs[i];
	//}


	u1 = initGuess[0];
	u2 = initGuess[1];
	//uHigh = u2 + 0.2;
	//uLow = u2 - 0.2;
	currU = u2;
	currAng = initGuess[2];

	start[0] = currU;
	start[1] = currAng;

	d_offset = (double*)malloc(size);
	d_matchPairs = (double*)malloc(size2);
	d_sum = (double*)malloc(size3);

	for ( i = 0 ; i < 12*numPairs ; i++ ) {
		d_matchPairs[i] = h_matchPairs[i];
	}

	/*
	if ( u1 >= 1.0 ) {
		double *ptr = &(poses_1[3*(numPoses-1)]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else if (u1 < 0.0) {
		double *ptr = &(poses_1[0]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_1[3*(int)(u1*(numPoses-1))]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	*/


	// Copy vectors from host memory to device memory

	//nelmin(nParam, &start[0], &xmin[0], &ynewlo, reqmin, &step[0],
	//	konvge, kcount, &icount, &numres, &ifault,
	//	numPoses, numPairs, d_matchPairs, d_offset, d_sum, poses_1, poses_2,
	//	uHigh, uLow, u1);

	//double offset[3];

	//double newSum = objFunc2(&start[0], d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset);
	//newSum = objFunc(&start[0], d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset);
	newSum = objFunc(&start[0], d_matchPairs, d_offset, d_sum, numPairs, pose1, poses_2, numPoses, angNom, angLim, uHigh, uLow, offset);


	resultSum[0] = newSum;
	resultParam[0] = currU;
	resultParam[1] = currAng;

	resultOffset[0] = offset[0];
	resultOffset[1] = offset[1];
	resultOffset[2] = offset[2];

	//printf("sum = %f\n", ynewlo);	
	//printf("minX = %f %f\n", xmin[0], xmin[1]);
	//printf("icount = %d\n", icount);

	// Free device memory
	if (d_matchPairs)
		free(d_matchPairs);
	if (d_offset)
		free(d_offset);
	if (d_sum)
		free(d_sum);

	// Free host memory
	//if (h_matchPairs)
	//	free(h_matchPairs);
	if (h_sum)
		free(h_sum);

}

double objFuncIndAngle(double *params, double *d_matchPairs, double *d_offset, double *d_sum, int N, double *pose1, double *poses_2, int numPoses, double angNom, double angLim, double uHigh, double uLow, double *resultOffset) {


	double currU = params[0];
	//double currAng = params[1];
	double currAng = angNom;

	// double pose1[3] = {0.0,0.0,0.0};
	double pose2[3] = {0.0,0.0,0.0};
	double point1[2];
	double point2[2];
	double ang1;
	double ang2;
	double ang2_off;
	double overlapPoint1Pose[3];
	double overlapPoint2Pose[3];
	double estPose1[3] = {0.0,0.0,0.0};
	double dist1;
	double vecAng1 = 0.0;
	double backR1[4];
	double foreR1[4];
	double R1[4];
	double estPose2[3];
	double dist2;
	double vecAng2 = 0.0;
	double backR2[4];
	double foreR2[4];
	double R2[4];
	double negCurve2Pose[3] = {0.0,0.0,0.0};
	double offsetR2[4];
	double finalVec2[2] = {0.0,0.0};
	double transVec2[2] = {0.0,0.0};
	double resPose2[3] = {0.0,0.0,0.0};
	double finalVec[2] = {0.0,0.0};
	double R2t[4];
	double tempVec1[2] = {0.0,0.0};
	double resVec[2] = {0.0,0.0};
	double tempVec2[2] = {0.0,0.0};
	double globalVec[2] = {0.0,0.0};
	double localVec[2] = {0.0,0.0};
	double localOffset[3] = {0.0,0.0,0.0};
	double h_offset[3] = {0.0,0.0,0.0};
	double totalSum;

	double tempVal1;

	//printf("%f %f :", currU, currAng);

	if (currU < 0.0) {
		//printf("%f\n",100000.0 * (0.15-currU));
		return 10000000.0 * (0.15-currU);
	}
	if (currU > 1.0 ) {
		//printf("%f\n",100000.0 * (currU-0.85));
		return 10000000.0 * (currU - 0.85);
	}

	if (currU < uLow) {
		//printf("%f\n",100000.0 * (uLow+0.05-currU));
		return 10000000.0 * (uLow+0.05-currU);
	}
	if (currU > uHigh) {
		//printf("%f\n",100000.0 * (currU-uHigh+0.05));
		return 10000000.0 * (currU - uHigh+0.05);
	}

	if (fabs(currAng-angNom) > angLim) {
		return 10000000.0 * (fabs(currAng-angNom)-angLim+0.05);
	}

	//printf("u1 = %lf\n", u1);

	// compute the point and angles from the parameters 
	/* 
	if ( u1 >= 1.0 ) {
		double *ptr = &(poses_1[3*(numPoses-1)]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
		//printf("foo %d %d\n", numPoses, 3*(numPoses-1));
	}
	else if (u1 < 0.0) {
		double *ptr = &(poses_1[0]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_1[3*(int)(u1*(numPoses-1))]);
		pose1[0] = ptr[0];
		pose1[1] = ptr[1];
		pose1[2] = ptr[2];
	}
	*/

	if (currU >= 1.0) {
		double *ptr = &(poses_2[3*(numPoses-1)]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}
	else if (currU < 0.0) {
		double *ptr = &(poses_2[0]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}
	else {
		double *ptr = &(poses_2[3*(int)(currU*(numPoses-1))]);
		pose2[0] = ptr[0];
		pose2[1] = ptr[1];
		pose2[2] = ptr[2];
	}

	//printf("pose1 = %lf %lf %lf\n", pose1[0], pose1[1], pose1[2]);
	//printf("pose2 = %lf %lf %lf\n", pose2[0], pose2[1], pose2[2]);


	//double point1[2];
	point1[0] = pose1[0];
	point1[1] = pose1[1];

	//double point2[2];
	point2[0] = pose2[0];
	point2[1] = pose2[1];

	//double ang1 = pose1[2];
	//double ang2 = pose2[2];
	//double ang2_off = ang2 + currAng;
	ang1 = pose1[2];
	ang2 = pose2[2];
	//ang2_off = ang1 + currAng;
	ang2_off = ang1;

	//offset = computeOffset(point1, point2, ang1, ang2 + currAng)

	//" corner points and orientations "
	//double overlapPoint1Pose[3];
	overlapPoint1Pose[0] = point1[0];
	overlapPoint1Pose[1] = point1[1];
	overlapPoint1Pose[2] = ang1;

	//double overlapPoint2Pose[3];
	overlapPoint2Pose[0] = point2[0];
	overlapPoint2Pose[1] = point2[1];
	overlapPoint2Pose[2] = ang2_off;

	//printf("overlapPoint1Pose = %f %f %f\n", overlapPoint1Pose[0], overlapPoint1Pose[1], overlapPoint1Pose[2]);
	//printf("overlapPoint2Pose = %f %f %f\n", overlapPoint2Pose[0], overlapPoint2Pose[1], overlapPoint2Pose[2]);


	/// CHECK

	//" convert the desired intersection point on curve 1 into global coordinates "
	//poseProfile1 = Pose([0.0,0.0,0.0])
	
	//double estPose1[3] = {0.0,0.0,0.0};
	//double dist1 = sqrt(estPose1[0]*estPose1[0] + estPose1[1]*estPose1[1]);
	//double vecAng1 = 0.0;
	dist1 = sqrt(estPose1[0]*estPose1[0] + estPose1[1]*estPose1[1]);
	vecAng1 = 0.0;

	//printf("estPose1 = %lf %lf %lf\n", estPose1[0], estPose1[1], estPose1[2]);

	//printf("dist1 = %lf\n", dist1);

	if (dist1 > 0.0) {
		vecAng1 = acos(estPose1[0]/dist1);
		if (asin(estPose1[1]/dist1) < 0) {
			vecAng1 = -vecAng1;
		}
	}
	else {
		vecAng1 = 0.0;
	}
	
	//printf("vecAng1 = %lf\n", vecAng1);

	//double backR1[4];
	buildMatrix(backR1, vecAng1, 0);

	//double foreR1[4];
	buildMatrix(foreR1, vecAng1, 1);

	//double R1[4];
	buildMatrix(R1, estPose1[2], 0);

	// now convert this point into a pose, and perform the inverse transform using corner2Pose 

	//desGlobalPose2 = Pose(overlapPoint1Pose)
	//double estPose2[3];
	estPose2[0] = overlapPoint1Pose[0];
	estPose2[1] = overlapPoint1Pose[1];
	estPose2[2] = overlapPoint1Pose[2];

	//double dist2 = sqrt(estPose2[0]*estPose2[0] + estPose2[1]*estPose2[1]);
	//double vecAng2 = 0.0;
	dist2 = sqrt(estPose2[0]*estPose2[0] + estPose2[1]*estPose2[1]);
	
	//printf("estPose2 = %f %f %f\n", estPose2[0], estPose2[1], estPose2[2]);

	if (dist2 > 0.0) {
		vecAng2 = acos(estPose2[0]/dist2);
		if (asin(estPose2[1]/dist2) < 0) {
			vecAng2 = -vecAng2;
		}
	}
	else {
		vecAng2 = 0.0;
	}

	//double backR2[4];
	buildMatrix(backR2, vecAng2, 0);

	//double foreR2[4];
	buildMatrix(foreR2, vecAng2, 1);

	//double R2[4];
	buildMatrix(R2, estPose2[2], 0);

	// CHECK2

	//" perform inverse offset from the destination pose "
	//negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
	//double negCurve2Pose[3] = {0.0,0.0,0.0};

	tempVal1 = overlapPoint2Pose[2];
	//printf("overlapPoint2Pose[2] = %lf\n", overlapPoint2Pose[2]);

	//double offsetR2[4];
	buildMatrix(offsetR2, overlapPoint2Pose[2], 0);

	//printf("offsetR2 = %f %f %f %f\n", cos(tempVal1), sin(tempVal1), -sin(tempVal1), cos(tempVal1));
	//printf("offsetR2 = %f %f %f %f\n", offsetR2[0], offsetR2[1], offsetR2[2], offsetR2[3]);

	//double finalVec2[2] = {0.0,0.0};
	finalVec2[0] = overlapPoint2Pose[0];
	finalVec2[1] = overlapPoint2Pose[1];

	//printf("finalVec2 = %f %f\n", finalVec2[0], finalVec2[1]);

	// transVec = dot(offsetR2, finalVec2)
	//double transVec2[2] = {0.0,0.0};
	//transVec2[0] = offsetR2[0]*finalVec2[0] + offsetR2[1]*finalVec2[1];
	//transVec2[1] = offsetR2[2]*finalVec2[0] + offsetR2[3]*finalVec2[1];
	dotProduct(offsetR2, finalVec2, transVec2);

	negCurve2Pose[0] = -transVec2[0];
	negCurve2Pose[1] = -transVec2[1];
	negCurve2Pose[2] = -overlapPoint2Pose[2];
	

	//printf("negCurve2Pose = %f %f %f\n", negCurve2Pose[0], negCurve2Pose[1], negCurve2Pose[2]);

	// CHECK3
	
	//" relative pose between pose 1 and pose 2 to make corners coincide and same angle "
	//resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
	//double resPose2[3] = {0.0,0.0,0.0};

	// pss
	//double finalVec[2] = {0.0,0.0};
	finalVec[0] = negCurve2Pose[0];
	finalVec[1] = negCurve2Pose[1];

	//finalVec = array([[negCurve2Pose[0]], [negCurve2Pose[1]]])

	//transVec = dot(transpose(self.R), finalVec)
	//double R2t[4];
	transposeMatrix(R2, R2t);

	//double tempVec1[2] = {0.0,0.0};
	dotProduct(R2t, finalVec, tempVec1);

	//double tempVec1[2] = {0.0,0.0};
	//dotProduct(R2, finalVec, tempVec1);

	//resVec = dot(self.backR, transVec)
	//double resVec[2] = {0.0,0.0};
	dotProduct(backR2, tempVec1, resVec);

	//resVec[0, 0] += self.dist
	resVec[0] += dist2;

	//tempVec = dot(self.foreR, resVec)
	//double tempVec2[2] = {0.0,0.0};
	dotProduct(foreR2, resVec, tempVec2);

	resPose2[0] = tempVec2[0];
	resPose2[1] = tempVec2[1];
	resPose2[2] = normalizeAngle(estPose2[2] + negCurve2Pose[2]);

	//printf("resPose2 = %f %f %f\n", resPose2[0], resPose2[1], resPose2[2]);
	

	//localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)

	//def convertGlobalPoseToLocal(self, pose):

	//" transform pnt to local coordinates"
	//double globalVec[2] = {0.0,0.0};
	globalVec[0] = resPose2[0];
	globalVec[1] = resPose2[1];

	//" perform translation correction "
	//tempVec = dot(self.backR, globalVec)
	dotProduct(backR1, globalVec, tempVec1);
	tempVec1[0] -= dist1;

	//transVec = dot(self.foreR, tempVec)
	dotProduct(foreR1, tempVec1, tempVec2);

	//" now, apply rotation correction with respect to origin "
	//localVec = dot(self.R, transVec)
	//double localVec[2] = {0.0,0.0};
	dotProduct(R1, tempVec2, localVec);

	//localPose = [localVec[0,0], localVec[1,0], normalizeAngle(resPose2[2] - self.estPose[2])]
	//return localPose
	//double localOffset[3] = {0.0,0.0,0.0};
	localOffset[0] = localVec[0];
	localOffset[1] = localVec[1];
	localOffset[2] = normalizeAngle(resPose2[2]-estPose1[2]);

	//double h_offset[3] = {0.0,0.0,0.0};
	h_offset[0] = localOffset[0];
	h_offset[1] = localOffset[1];
	h_offset[2] = localOffset[2];


	//printf("offset = %f %f %f\n", h_offset[0], h_offset[1], h_offset[2]);
	
	resultOffset[0] = h_offset[0];
	resultOffset[1] = h_offset[1];
	resultOffset[2] = h_offset[2];

	d_offset[0] = h_offset[0];
	d_offset[1] = h_offset[1];
	d_offset[2] = h_offset[2];


	icpError(d_matchPairs, d_offset, d_sum, N);

	// Copy result from device memory to host memory
	// h_sum contains the result in host memory

	//double totalSum = d_sum[0];
	totalSum = d_sum[0];

	//printf("totalSum = %lf\n");

	return totalSum;

}
