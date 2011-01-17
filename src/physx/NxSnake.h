#pragma once

#ifndef WIN32
#define WIN32
#endif

#ifndef NXSNAKE_H
#define NXSNAKE_H

//#define NxReal double;

# if defined(_MSC_VER) && !defined(GLUT_NO_WARNING_DISABLE)
#  pragma warning (disable:4244)  /* Disable bogus VC++ 4.2 conversion warnings. */
#  pragma warning (disable:4305)  /* VC++ 5.0 version of above warning. */
# endif

#include <stdio.h>

#include "NxPhysics.h"
#include "NxCooking.h"
#include "Stream.h"
//#include "NxOgreMatrix.h"

//#include "NxPhysics.h"
//#include "DrawObjects.h"
//#include "UserData.h"
//#include "HUD.h"
//#include <GL/glut.h>



//#include "NxCooking.h"

//#include "DebugRenderer.h"
//#include "UserAllocator.h"
//#include "ErrorStream.h"

//#include "CommonCode.h"

//using namespace NxOgre;

class NxSnake
{
public:
	NxSnake(double *quatR, double *pos, int numSegs, double segLength, double segHeight, double segWidth, double friction);
	~NxSnake(void);

	// parameters
	NxQuat ori;
	NxVec3 tPos;
	int numSegs;
	double segLength;
	double segWidth;
	double segHeight;
	double friction;

	// material indices
	int dIndex;
	int pIndex;

	// Physics SDK globals
	NxPhysicsSDK*     gPhysicsSDK;
	NxScene*          gScene;
	NxVec3            gDefaultGravity;

	// Actor globals
	NxActor* groundPlane;

	int numWalls;
	NxVec3 *wallSet[100];
	int wallSize[100];

	NxActor **segments;
	NxRevoluteJoint **joints;
	NxReal *targetAngle;
	NxVec3 *savePoses;
	NxQuat *saveOrientations;

	bool isPerturb;
	NxVec3 ApplyForceToActor(NxActor* actor, const NxVec3& forceDir, const NxReal forceStrength);
	void perturb();

	void RunPID();
	void Step();
	void ExitCallback();
	NxActor* CreateGroundPlane();
	NxCCDSkeleton* CreateCCDSkeleton(float sizex, float sizey, float sizez);
	NxCCDSkeleton* CreateCCDSkeleton2();
	NxActor* CreateTriangleMeshWall(int size, NxVec3 *points);
	NxActor* CreateSnake();
	NxRevoluteJoint* CreateRevoluteJoint(NxActor* a0, NxActor* a1, NxVec3 globalAnchor, NxVec3 globalAxis);
	void InitNx();
	void ReleaseNx();
	void ResetNx();
	void StartPhysics();
	void GetPhysicsResults();

	//void setServo();
	//void getServo();
	//void getServoCmd();
	//void setJointTorque();
	//void getJointTorque();

	void addWall(int numPoints, double *points);
	void createWalls();

	void createWall(int numPoints, double *points);
	void setServo(int i, double angle);
	double getServo(int i);
	double getServoCmd(int i);
	void setJointTorque(int i, double torque);
	double getJointTorque(int i);
	//void getGlobalPose(int i, NxF32 *d);

	//NxVec3 getGlobalPosition(int i);
	//NxQuat getGlobalOrientationQuat(int i);
	void getGlobalPosition(int i, double *x, double *y, double *z);
	void getGlobalOrientationQuat(int i, double *x, double *y, double *z, double *w);

	void savePose();
	void restorePose();

	//def getActualSegPose(self, i):
	//def getActualJointPose(self, i):


};

#endif
