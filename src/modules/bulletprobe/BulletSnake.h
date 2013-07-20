#include <btBulletDynamicsCommon.h>

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

class BulletSnake
{
public:
	BulletSnake(double *quatR, double *pos, int numSegs, double segLength, double segHeight, double segWidth, double friction);
	~BulletSnake(void);

	void initPhysics();
	void addWall(int numPoints, double *points);
	void createWalls();
	void createGroundPlane();
	void createSnake();
	void RunPID();
	void Step();


	void perturb();
	void savePose();
	void restorePose();

	void setServo(int i, double angle);
	double getServo(int i);
	double getServoCmd(int i);
	void setJointTorque(int i, double torque);
	double getJointTorque(int i);
	void getGlobalPosition(int i, double *x, double *y, double *z);
	void getGlobalOrientationQuat(int i, double *x, double *y, double *z, double *w);


	// parameters
	//NxQuat ori;
	//NxVec3 tPos;
	btVector3 tPos;
	btQuaternion ori;

	int numSegs;
	double segLength;
	double segWidth;
	double segHeight;
	double friction;

	btDiscreteDynamicsWorld* dynamicsWorld;
	
	// Build the broadphase
	btBroadphaseInterface* broadphase;

	// Set up the collision configuration and dispatcher
	btDefaultCollisionConfiguration* collisionConfiguration;
	btCollisionDispatcher* dispatcher;

	// The actual physics solver
	btSequentialImpulseConstraintSolver* solver;
	
	// snake segments
	btRigidBody *segBodies[100];
	btCollisionShape *segShapes[100];
	btDefaultMotionState *segMotionStates[100];

	// snake joints
	btHingeConstraint *joints[100];
	btScalar targetAngle[100];



	// walls
	int numWalls;
	int wallSize[100];
	btVector3 *wallSet[100];
	btCollisionShape *wallShape[100];
	btDefaultMotionState *wallMotionState[100];
	btRigidBody *wallRigidBody[100];
	
	// ground plane
	btCollisionShape* groundShape;
	btDefaultMotionState* groundMotionState;
	btRigidBody* groundRigidBody;

};

