
#ifdef _WIN32
#include "stdafx.h"
#endif

#include <iostream>
#include "BulletSnake.h"


BulletSnake::BulletSnake(double *quatR, double *pos, int numSegs, double segLength, double segHeight, double segWidth, double friction)
{

	this->numSegs = numSegs;
	this->segLength = segLength;
	this->segWidth = segWidth;
	this->segHeight = segHeight;
	this->friction = friction;

	std::cout << this->numSegs << " " << this->segLength << std::endl;

	numWalls = 0;
	
	tPos = btVector3(pos[0],pos[1],pos[2]);
	ori = btQuaternion(quatR[0], quatR[1], quatR[2], quatR[3]);

	initPhysics();

}

void BulletSnake::perturb() {

}

void BulletSnake::savePose() {

}

void BulletSnake::restorePose() {

}


void BulletSnake::setServo(int i, double angle) {
	targetAngle[i] = angle;
}

double BulletSnake::getServo(int i) {
	return joints[i]->getHingeAngle();
}

double BulletSnake::getServoCmd(int i) {
	return targetAngle[i];
}

void BulletSnake::setJointTorque(int i, double torque){ 
	joints[i]->setMaxMotorImpulse(torque);
}

double BulletSnake::getJointTorque(int i) {
	return joints[i]->getMaxMotorImpulse();
}

void BulletSnake::getGlobalPosition(int i, double *x, double *y, double *z) {

	btVector3 result = segBodies[i]->getCenterOfMassPosition();
	*x = result[0];
	*y = result[1];
	*z = result[2];
}

void BulletSnake::getGlobalOrientationQuat(int i, double *x, double *y, double *z, double *w) {

	btQuaternion result = segBodies[i]->getOrientation();

	*x = result[0];
	*y = result[1];
	*z = result[2];
	*w = result[3];
}


void BulletSnake::initPhysics() {

	// Build the broadphase
	broadphase = new btDbvtBroadphase();

	// Set up the collision configuration and dispatcher
	collisionConfiguration = new btDefaultCollisionConfiguration();
	dispatcher = new btCollisionDispatcher(collisionConfiguration);

	// The actual physics solver
	solver = new btSequentialImpulseConstraintSolver;

	dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher,broadphase,solver,collisionConfiguration);
	dynamicsWorld->setGravity(btVector3(0,-9.8,0));

	createGroundPlane();
	createSnake();

}

BulletSnake::~BulletSnake() {
	// Clean up behind ourselves like good little programmers
	delete dynamicsWorld;
	delete solver;
	delete dispatcher;
	delete collisionConfiguration;
	delete broadphase;
}
void BulletSnake::Step() {
	RunPID();

	dynamicsWorld->stepSimulation(1/60.f,10,1/(4*60.f));
}

void BulletSnake::RunPID() {

	for ( int i = 0 ; i < numSegs-1 ; i++ ) {
		joints[i]->setMotorTarget(targetAngle[i], 0.1);
		//joints[i]->enableAngularMotor(true, fDesiredAngularVel, m_fMuscleStrength);
	}

	/*
	float pgain = 10;
	float dgain = .001f;
	float tolerance = 0.01 * M_PI/180.0;
	float m_fMuscleStrength = 1.0f;

	for ( int i = 0 ; i < numSegs-1 ; i++ ) {

		double init_phi = 0.0;
		double phi = joints[i]->getHingeAngle();
		double err = targetAngle[i] - (phi + init_phi);

		if (fabs(err) > tolerance) {

			//printf("%d %f %f\n", i, phi, targetAngle[i]);
			double vel = pgain*err + dgain*joints[i]->getVelocity();
			revDesc.motor.velTarget = vel;
			//printf("%f %f %f %f\n", pgain*err, dgain*joints[i]->getVelocity(), phi, vel);
			//printf("Angle: %f\n", phi);
			//printf("Velocity = %f\n", vel);

			btScalar fDesiredAngularVel = vel;
			joints[i]->enableAngularMotor(true, fDesiredAngularVel, m_fMuscleStrength);
		}

	}
	*/
}

void BulletSnake::createSnake()
{

	float  m_fMuscleStrength = 10.0;

	// Set the box starting height to 3.5m so box starts off falling onto the ground
	float startHeight = 0.1;

	std::cout << "numSegs = " << numSegs << std::endl;

	for ( int i = 0 ; i < numSegs ; i++ ) {

		segShapes[i] = new btBoxShape(btVector3(segLength/2.0,segHeight/2.0,segWidth/2.0));

		// dimensions: segLength/2.0,segHeight/2.0,segWidth/2.0
		// mass: 1

		// starting position and orientation
		btVector3 segPose(i*segLength, startHeight + 0*(i%2 * 0.2),0);

		// rotate by quaternion
		// give offset and orientation of motion state
		btTransform trans = btTransform(ori);
		segPose = trans*segPose;

		// offset by tPos
		segPose = segPose + tPos;

		std::cout << i << " " << segPose[0] << " " << segPose[1] << " " << segPose[2] <<  std::endl;

		segMotionStates[i] = new btDefaultMotionState(btTransform(ori,segPose));


		// mass and inertia
		btScalar mass = 1;
		btVector3 segInertia(0,0,0);
		segShapes[i]->calculateLocalInertia(mass,segInertia);

		// construction info
		btRigidBody::btRigidBodyConstructionInfo segRigidBodyCI(mass,segMotionStates[i],segShapes[i],segInertia);
		segRigidBodyCI.m_friction = this->friction;
		segRigidBodyCI.m_linearDamping = 0.8;

		// rigid body 
		segBodies[i] = new btRigidBody(segRigidBodyCI);

		dynamicsWorld->addRigidBody(segBodies[i]);
		segBodies[i]->setActivationState(DISABLE_DEACTIVATION);
		//segBodies[i]->setSleepingThresholds (0.0, 0.0);

		std::cout << "sleeping thresholds: " << segBodies[i]->getLinearSleepingThreshold() << " " << segBodies[i]->getAngularSleepingThreshold() << std::endl;
	}

	
	for ( int i = 0 ; i < numSegs-1 ; i++ ) {
		btRigidBody *seg1 = segBodies[i];
		btRigidBody *seg2 = segBodies[i+1];

		btVector3 armJointPos = btVector3(segLength/2.0+i*segLength,startHeight,0);
		btTransform trans = btTransform(ori);
		armJointPos = trans*armJointPos;

		btVector3 globalAnchor = tPos + armJointPos;
		btVector3 globalAxis = btVector3(0,-1,0);

		btVector3 pivotInA = btVector3(segLength/2.0, 0.0, 0.0);
		btVector3 pivotInB = btVector3(-segLength/2.0, 0.0, 0.0);
		btVector3 axisInA = btVector3(0,-1,0);
		btVector3 axisInB = btVector3(0,-1,0);

		joints[i] = new btHingeConstraint(*seg1, *seg2, pivotInA, pivotInB, axisInA, axisInB);
		joints[i]->setLimit(btScalar(-M_PI/2.0), btScalar(M_PI/2.0));
		joints[i]->enableMotor(true);
		joints[i]->setMaxMotorImpulse(m_fMuscleStrength);
		dynamicsWorld->addConstraint(joints[i], true);

		targetAngle[i] = 0.0;

		//btScalar fDesiredAngularVel = 1000000.f * fAngleError/ms;
		//hingeC->enableAngularMotor(true, fDesiredAngularVel, m_fMuscleStrength);
	}
	

}
void BulletSnake::createGroundPlane() {

	groundShape = new btStaticPlaneShape(btVector3(0,1,0),1);
	groundMotionState = new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),btVector3(0,-1,0)));
	btRigidBody::btRigidBodyConstructionInfo groundRigidBodyCI(0,groundMotionState,groundShape,btVector3(0,0,0));
	groundRigidBodyCI.m_friction = 0.0;
	groundRigidBody = new btRigidBody(groundRigidBodyCI);
	dynamicsWorld->addRigidBody(groundRigidBody);
}


void BulletSnake::addWall(int numPoints, double *points) {

	btVector3 *wallPoints = new btVector3[numPoints];
	
	for ( int i = 0 ; i < numPoints ; i++ ) {
		printf("%f %f %f\n", points[2*i], 0.0, points[2*i+1]);
		wallPoints[i] = btVector3(points[i*2], 0.0, points[2*i+1]);
	}

	wallSet[numWalls] = wallPoints;
	wallSize[numWalls] = numPoints;
	numWalls++;
}


void BulletSnake::createWalls() {

	// create rigid static walls
	for ( int k = 0 ; k < numWalls ; k++ ) {

		int size = wallSize[k];
		btVector3 *points = wallSet[k];
		btVector3 *verts = new btVector3[ 2 * size];

		for ( int i = 0 ; i < size ; i++ ) {
			verts[2*i] = btVector3(points[i][0], -1.0, points[i][2]);
			verts[2*i+1] = btVector3(points[i][0], 1.0, points[i][2]);
		}
		
		btTriangleMesh *mTriMesh = new btTriangleMesh();
		for ( int i = 0 ; i < size-1 ; i++ ) {
			mTriMesh->addTriangle(verts[2*i+2],verts[2*i+1],verts[2*i]);
			mTriMesh->addTriangle(verts[2*i+1],verts[2*i+2],verts[2*i+3]);
		}

		delete verts;
		wallShape[k] = new btBvhTriangleMeshShape(mTriMesh,true);

		// Now use mTriMeshShape as your collision shape.
		// Everything else is like a normal rigid body
		wallMotionState[k] = new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),btVector3(0,0,0)));
		btRigidBody::btRigidBodyConstructionInfo wallRigidBodyCI(0,wallMotionState[k],wallShape[k],btVector3(0,0,0));
		wallRigidBodyCI.m_friction = 1.0;
		wallRigidBody[k] = new btRigidBody(wallRigidBodyCI);
		dynamicsWorld->addRigidBody(wallRigidBody[k]);
	}
}
