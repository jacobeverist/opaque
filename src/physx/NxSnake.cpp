#ifndef WIN32
#define WIN32
#endif


#include "NxSnake.h"


class MyContactReport : public NxUserContactReport
{

	private:
		int numContacts;
		NxActor *plane;
		int defaultIndex;
		int planeIndex;

	public:
		MyContactReport() {
			numContacts = 0;
		}

		void setIndices(int def, int pl) {
			defaultIndex = def;
			planeIndex = pl;
		}

		void setPlane(NxActor *groundPlane) {
			plane = groundPlane;
		}

		int getNumContacts() {
			int copyNum = numContacts;
			numContacts = 0;
			return copyNum;
		}

	virtual void  onContactNotify(NxContactPair& pair, NxU32 events)
	{

		return;
		/*
		if(!pair.isDeletedActor[0])
			{
			MyCubeObject* Object = (MyCubeObject*)pair.actors[0]->userData;
			if(Object)	Object->events = events;
			}

		if(!pair.isDeletedActor[1])
			{
			MyCubeObject* Object = (MyCubeObject*)pair.actors[1]->userData;
			if(Object)	Object->events = events;
			}
		*/

		//if ( pair.actors[0] != plane && pair.actors[1] != plane ) {
		if ( 1 ) {

			//printf("Non-Plane contact %d %d\n", pair.actors[0], pair.actors[1]);
			//NxVec3 fricForce = pair.sumFrictionForce;
			//double mag = sqrt(fricForce.x * fricForce.x  + fricForce.y * fricForce.y + fricForce.z * fricForce.z);
			
			//planeDesc.materialIndex	= floorMaterial->getMaterialIndex();
			//actors[0];
			//actors[1];

			
			//double staticFric = pair.actors[0]->getScene().getMaterialFromIndex(defaultIndex)->getStaticFriction();
			//double dynFric = pair.actors[0]->getScene().getMaterialFromIndex(defaultIndex)->getDynamicFriction();
			//printf("%d Frictions: %f %f\n", defaultIndex, staticFric, dynFric);
			//double staticFric2 = pair.actors[0]->getScene().getMaterialFromIndex(planeIndex)->getStaticFriction();
			//double dynFric2 = pair.actors[0]->getScene().getMaterialFromIndex(planeIndex)->getDynamicFriction();
			//printf("%d Frictions: %f %f\n\n", planeIndex, staticFric2, dynFric2);
			
			//printf("mag = %f\n", mag);

			//if ( mag > 0.00001 )
			//	printf("mag = %f\n", mag);

		}
		else {
			printf("Plane contact\n\n");
		}

		//printf("Normal Force: %f, %f, %f\n", pair.sumNormalForce.x, pair.sumNormalForce.y, pair.sumNormalForce.z);
		//printf("Friction Force: %f, %f, %f\n\n", pair.sumFrictionForce.x, pair.sumFrictionForce.y, pair.sumFrictionForce.z);
//		if(events & NX_NOTIFY_ON_START_TOUCH)	printf("Start touch\n");
//		if(events & NX_NOTIFY_ON_END_TOUCH)		printf("End touch\n");

		// Iterate through contact points
		NxContactStreamIterator i(pair.stream);
		//user can call getNumPairs() here
		while(i.goNextPair())
			{
			//user can also call getShape() and getNumPatches() here
			while(i.goNextPatch())
				{
				//user can also call getPatchNormal() and getNumPoints() here
				const NxVec3& contactNormal = i.getPatchNormal();
				while(i.goNextPoint())
					{

					numContacts++;

					//user can also call getPoint() and getSeparation() here
					const NxVec3& contactPoint = i.getPoint();

					//printf("%f  - %f\n", i.getPointNormalForce(), pair.sumNormalForce.magnitude());
					//NxVec3 contactForce = (showPointContacts /*&& i.getShapeFlags()&NX_SF_POINT_CONTACT_FORCE*/) ?  contactNormal * i.getPointNormalForce() : pair.sumNormalForce; 
					NxVec3 contactForce = contactNormal * i.getPointNormalForce(); 
					NxVec3 contactArrowForceTip = contactPoint + contactForce * 0.1f;
					NxVec3 contactArrowFrictionTip = contactPoint + pair.sumFrictionForce* 0.1f;
					NxVec3 contactArrowPenetrationTip = contactPoint - contactNormal * i.getSeparation() * 20.0f;

					}
				}
			}
	}

}gMyContactReport;

//#include "Timing.h"

NxSnake::NxSnake(double *quatR, double *pos, int numSegs, double segLength, double segHeight, double segWidth, double friction)
{

	printf("Quat: %f %f %f %f\n", quatR[0], quatR[1], quatR[2], quatR[3]);
	printf("Pos: %f %f %f\n", pos[0], pos[1], pos[2]);

	//this->pos.set(pos);

	numWalls = 0;
	dIndex = 0;
	pIndex = 0;

	tPos = NxVec3(pos[0],pos[1],pos[2]);
	ori = NxQuat();
	ori.setXYZW(quatR[0], quatR[1], quatR[2], quatR[3]);

	printf("NxVec3: %f, %f, %f\n", tPos.x, tPos.y, tPos.z);

	//this->pos.sety(pos[1]);
	//this->pos.setz(pos[2]);

	this->numSegs = numSegs;
	this->segLength = segLength;
	this->segWidth = segWidth;
	this->segHeight = segHeight;
	this->friction = friction;

	// Physics SDK globals
	gPhysicsSDK = NULL;
	gScene = NULL;
	gDefaultGravity = NxVec3(0,-9.8f,0);

	// Actor globals
	groundPlane = NULL;

	segments = new NxActor * [this->numSegs];
	joints = new NxRevoluteJoint * [this->numSegs-1];
	targetAngle = new NxReal[this->numSegs-1];

	savePoses = new NxVec3[this->numSegs];
	saveOrientations = new NxQuat[this->numSegs];

	for ( int i = 0 ; i < this->numSegs-1 ; i++ ) {
		segments[i] = NULL;
		joints[i] = NULL;
		targetAngle[i] = 0.0;
	}
	for ( int i = 0 ; i < this->numSegs ; i++ ) {
		savePoses[i].x = 0.0;
		savePoses[i].y = 0.0;
		savePoses[i].z = 0.0;

		saveOrientations[i].x = 0.0;
		saveOrientations[i].y = 1.0;
		saveOrientations[i].z = 0.0;
		saveOrientations[i].w = 0.0;
	}


    InitNx();

}

NxSnake::~NxSnake(void)
{
	delete segments;
	delete joints;
	delete targetAngle;
	ReleaseNx();
}

void NxSnake::RunPID() {

	float pgain = 10;
	float dgain = .001f;

	for ( int i = 0 ; i < numSegs-1 ; i++ ) {
		NxRevoluteJointDesc revDesc;
		joints[i]->saveToDesc(revDesc);
	
		NxReal init_phi = 0.0;
		NxReal phi = joints[i]->getAngle();

		NxReal err = targetAngle[i] - (phi + init_phi);
		//printf("%d %f %f\n", i, phi, targetAngle[i]);
		NxReal vel = pgain*err + dgain*joints[i]->getVelocity();
		revDesc.motor.velTarget = vel;
		//printf("%f %f %f %f\n", pgain*err, dgain*joints[i]->getVelocity(), phi, vel);
		//printf("Angle: %f\n", phi);
		//printf("Velocity = %f\n", vel);

		joints[i]->loadFromDesc(revDesc);
	}
}


void NxSnake::Step()
{

    if (gScene)
	{
		for ( int i = 0 ; i < 1 ; i++ ) {
			GetPhysicsResults();

			//printf("Contacts: %d\n", gMyContactReport.getNumContacts());
			
			
			RunPID();

			StartPhysics();
		}
	}
}


void NxSnake::ExitCallback()
{
	ReleaseNx();
}


NxActor* NxSnake::CreateGroundPlane()
{
    // Create a plane with default descriptor
    NxPlaneShapeDesc planeDesc;

	
	NxMaterial* floorMaterial = gScene->getMaterialFromIndex(pIndex); 
	//floorMaterial->setRestitution(0.5);
	//floorMaterial->setStaticFriction(0.1);
	//floorMaterial->setDynamicFriction(0.1);
	//floorMaterial->setFrictionCombineMode(NX_CM_MIN);
	planeDesc.materialIndex	= floorMaterial->getMaterialIndex();

    NxActorDesc actorDesc;
    actorDesc.shapes.pushBack(&planeDesc);

	return gScene->createActor(actorDesc);
}

NxCCDSkeleton* NxSnake::CreateCCDSkeleton(float sizex, float sizey, float sizez)
{
	NxU32 triangles[3 * 12] = { 
		0,1,3,
		0,3,2,
		3,7,6,
		3,6,2,
		1,5,7,
		1,7,3,
		4,6,7,
		4,7,5,
		1,0,4,
		5,1,4,
		4,0,2,
		4,2,6
	};

	NxVec3 points[8];

	// Static mesh
	points[0].set( sizex, -sizey, -sizez);
	points[1].set( sizex, -sizey,  sizez);
	points[2].set( sizex,  sizey, -sizez);
	points[3].set( sizex,  sizey,  sizez);

	points[4].set(-sizex, -sizey, -sizez);
	points[5].set(-sizex, -sizey,  sizez);
	points[6].set(-sizex,  sizey, -sizez);
	points[7].set(-sizex,  sizey,  sizez);

	NxSimpleTriangleMesh stm;
	stm.numVertices = 8;
	stm.numTriangles = 6*2;
	stm.pointStrideBytes = sizeof(NxVec3);
	stm.triangleStrideBytes = sizeof(NxU32)*3;

	stm.points = points;
	stm.triangles = triangles;
	stm.flags |= NX_MF_FLIPNORMALS;
	return gPhysicsSDK->createCCDSkeleton(stm);
}

void NxSnake::addWall(int numPoints, double *points) {

	NxVec3 *wallPoints = new NxVec3[numPoints];
	
	for ( int i = 0 ; i < numPoints ; i++ ) {
		printf("%f %f %f\n", points[2*i], 0.0, points[2*i+1]);
		wallPoints[i] = NxVec3(points[i*2], 0.0, points[2*i+1]);
	}

	wallSet[numWalls] = wallPoints;
	wallSize[numWalls] = numPoints;
	numWalls++;
}


void NxSnake::createWalls() {

	NxInitCooking();

	for ( int k = 0 ; k < numWalls ; k++ ) {

		int size = wallSize[k];
		NxVec3 *points = wallSet[k];
		
		NxTriangleMeshDesc* triangleMeshDesc = NULL;
		NxVec3 boxDim(1.0f, 1.0f, 1.0f);

		NxVec3 *verts = new NxVec3[ 2 * size];

		for ( int i = 0 ; i < size ; i++ ) {
			verts[2*i] = NxVec3(points[i].x, -1.0, points[i].z);
			verts[2*i+1] = NxVec3(points[i].x, 1.0, points[i].z);
		}

		NxU32 *indices = new NxU32[ 6 * (size-1)];

		for ( int i = 0 ; i < size-1 ; i++ ) {
			indices[6*i + 0] = 2*i+2;
			indices[6*i + 1] = 2*i+1;
			indices[6*i + 2] = 2*i;
			indices[6*i + 3] = 2*i+1;
			indices[6*i + 4] = 2*i+2;
			indices[6*i + 5] = 2*i+3;
		}

		// Create descriptor for triangle mesh
		if (!triangleMeshDesc)
		{
			triangleMeshDesc	= new NxTriangleMeshDesc();
			assert(triangleMeshDesc);
		}
		triangleMeshDesc->numVertices			= 2*size;
		triangleMeshDesc->pointStrideBytes		= sizeof(NxVec3);
		triangleMeshDesc->points				= verts;
		triangleMeshDesc->numTriangles			= 2*(size-1);
		triangleMeshDesc->flags					= 0;
		triangleMeshDesc->triangles				= indices;
		triangleMeshDesc->triangleStrideBytes	= 3 * sizeof(NxU32);

		// The actor has one shape, a triangle mesh
		MemoryWriteBuffer buf;

		for ( int i = 0 ; i < 2*size ; i++ ) {

			NxVec3 temp = verts[i];
			printf("Mesh Point: %f %f %f\n", temp[0], temp[1], temp[2]);
		}

		bool status = NxCookTriangleMesh(*triangleMeshDesc, buf);
		NxTriangleMesh* pMesh = NULL;
		if (status)
		{
			int tryCount = 1000;

			// keep trying to create mesh until it succeeds
			//while ( pMesh == NULL && tryCount != 0) {
			//while ( pMesh == NULL ) {
			//	pMesh = gPhysicsSDK->createTriangleMesh(MemoryReadBuffer(buf.data));
				//tryCount = tryCount - 1;
			//}

			pMesh = gPhysicsSDK->createTriangleMesh(MemoryReadBuffer(buf.data));

			if ( pMesh == NULL ) {
				printf("createTriangleMesh failed!\n");
				assert(false);
			}

		}
		else
		{
			assert(false);
			pMesh = NULL;
		}
		// Create TriangleMesh above code segment.

		NxTriangleMeshShapeDesc tmsd;
		tmsd.meshData		= pMesh;
		tmsd.userData		= triangleMeshDesc;
		tmsd.localPose.t	= NxVec3(0, boxDim.y, 0);
		tmsd.meshPagingMode = NX_MESH_PAGING_AUTO;
		tmsd.materialIndex	= 0;

		NxActorDesc actorDesc;
		NxBodyDesc  bodyDesc;
		
		assert(tmsd.isValid());
		actorDesc.shapes.pushBack(&tmsd);
		//Dynamic triangle mesh don't be supported anymore. So body = NULL
		actorDesc.body			= NULL;		
		actorDesc.globalPose.t	= NxVec3(0.0f, 0.0f, 0.0f);

		NxActor *actor = NULL;
		if (pMesh)
		{
			// Save mesh in userData for drawing
			pMesh->saveToDesc(*triangleMeshDesc);
			//
			assert(actorDesc.isValid());
			actor = gScene->createActor(actorDesc);
			assert(actor);
		}

		for ( int i = 0 ; i < this->numSegs ; i++ ) {
			gScene->setActorPairFlags(*actor,*segments[i],NX_NOTIFY_ON_TOUCH|NX_NOTIFY_ON_START_TOUCH|NX_NOTIFY_ON_END_TOUCH|NX_NOTIFY_FORCES);
		}
	}

	NxCloseCooking();

}

void NxSnake::createWall(int numPoints, double *points) {

	//NxVec3 *wallPoints = new NxVec3[numPoints];
	NxVec3 *wallPoints = new NxVec3[numPoints];
	
	for ( int i = 0 ; i < numPoints ; i++ ) {
		printf("%f %f %f\n", points[2*i], 0.0, points[2*i+1]);
		wallPoints[i] = NxVec3(points[i*2], 0.0, points[2*i+1]);
	}

	printf("create wall with %d points\n", numPoints);

	NxActor *actor = CreateTriangleMeshWall(numPoints, wallPoints);

	delete wallPoints;
	
	for ( int i = 0 ; i < this->numSegs ; i++ ) {
		//gScene->setActorPairFlags(*actor,*segments[i],NX_NOTIFY_ON_START_TOUCH|NX_NOTIFY_ON_TOUCH|NX_NOTIFY_ON_END_TOUCH);
		gScene->setActorPairFlags(*actor,*segments[i],NX_NOTIFY_ON_TOUCH|NX_NOTIFY_ON_START_TOUCH|NX_NOTIFY_ON_END_TOUCH|NX_NOTIFY_FORCES);
	}

}

NxActor* NxSnake::CreateTriangleMeshWall(int size, NxVec3 *points)
{

	NxTriangleMeshDesc* triangleMeshDesc = NULL;
	NxVec3 boxDim(1.0f, 1.0f, 1.0f);

	NxVec3 *verts = new NxVec3[ 2 * size];

	for ( int i = 0 ; i < size ; i++ ) {
		verts[2*i] = NxVec3(points[i].x, -1.0, points[i].z);
		verts[2*i+1] = NxVec3(points[i].x, 1.0, points[i].z);
	}

	NxU32 *indices = new NxU32[ 6 * (size-1)];

	for ( int i = 0 ; i < size-1 ; i++ ) {
		indices[6*i + 0] = 2*i+2;
		indices[6*i + 1] = 2*i+1;
		indices[6*i + 2] = 2*i;
		indices[6*i + 3] = 2*i+1;
		indices[6*i + 4] = 2*i+2;
		indices[6*i + 5] = 2*i+3;
	}

	// Create descriptor for triangle mesh
    if (!triangleMeshDesc)
	{
		triangleMeshDesc	= new NxTriangleMeshDesc();
		assert(triangleMeshDesc);
	}
	triangleMeshDesc->numVertices			= 2*size;
	triangleMeshDesc->pointStrideBytes		= sizeof(NxVec3);
	triangleMeshDesc->points				= verts;
	triangleMeshDesc->numTriangles			= 2*(size-1);
	triangleMeshDesc->flags					= 0;
	triangleMeshDesc->triangles				= indices;
	triangleMeshDesc->triangleStrideBytes	= 3 * sizeof(NxU32);

	// The actor has one shape, a triangle mesh
	NxInitCooking();
	MemoryWriteBuffer buf;

	for ( int i = 0 ; i < 2*size ; i++ ) {

		NxVec3 temp = verts[i];
		printf("Mesh Point: %f %f %f\n", temp[0], temp[1], temp[2]);
	}

	bool status = NxCookTriangleMesh(*triangleMeshDesc, buf);
	NxTriangleMesh* pMesh = NULL;
	if (status)
	{
		int tryCount = 1000;

		// keep trying to create mesh until it succeeds
		//while ( pMesh == NULL && tryCount != 0) {
		while ( pMesh == NULL ) {
			pMesh = gPhysicsSDK->createTriangleMesh(MemoryReadBuffer(buf.data));
			tryCount = tryCount - 1;
		}

		if ( pMesh == NULL ) {
			printf("createTriangleMesh failed for 1000th time!\n");
			assert(false);
		}

	}
	else
	{
		assert(false);
		pMesh = NULL;
	}
	NxCloseCooking();
	// Create TriangleMesh above code segment.

	NxTriangleMeshShapeDesc tmsd;
	tmsd.meshData		= pMesh;
	tmsd.userData		= triangleMeshDesc;
	tmsd.localPose.t	= NxVec3(0, boxDim.y, 0);
	tmsd.meshPagingMode = NX_MESH_PAGING_AUTO;
	tmsd.materialIndex	= 0;

	NxActorDesc actorDesc;
	NxBodyDesc  bodyDesc;
	
	assert(tmsd.isValid());
	actorDesc.shapes.pushBack(&tmsd);
	//Dynamic triangle mesh don't be supported anymore. So body = NULL
	actorDesc.body			= NULL;		
	actorDesc.globalPose.t	= NxVec3(0.0f, 0.0f, 0.0f);

	if (pMesh)
	{
		// Save mesh in userData for drawing
		pMesh->saveToDesc(*triangleMeshDesc);
		//
		assert(actorDesc.isValid());
		NxActor *actor = gScene->createActor(actorDesc);
		assert(actor);

		return actor;
	}

	return NULL;
}

NxActor* NxSnake::CreateSnake()
{

	// Set the box starting height to 3.5m so box starts off falling onto the ground
	NxReal boxStartHeight = 0.0; 

	// Add a single-shape actor to the scene
	NxActorDesc actorDesc;
	NxBodyDesc bodyDesc;

	NxActor *pActor;

	for ( int i = 0 ; i < numSegs ; i++ ) {

		// The actor has one shape, a box, 1m on a side
		NxBoxShapeDesc boxDesc;
		boxDesc.dimensions.set(segLength/2.0,segHeight/2.0,segWidth/2.0);
		boxDesc.localPose.t = NxVec3(0, 0, 0);
		boxDesc.materialIndex	= dIndex;

		boxDesc.ccdSkeleton = CreateCCDSkeleton(segLength/2.0, 0.1f, segWidth/2.0);
	    if (1)  boxDesc.shapeFlags |= NX_SF_DYNAMIC_DYNAMIC_CCD;  

		actorDesc.shapes.pushBack(&boxDesc);

		actorDesc.body			= &bodyDesc;
		actorDesc.density		= 0.001f;
		
		
		//" control the position of each segment "
		//segPos = ogre.Vector3(i*self.segLength, 0.011*(i%2), 0.0)
		//actPos = R*segPos
		//body.setPosition(actPos + pos)
		//body.setOrientation(R)

		//NxQuat ori;
		//NxVec3 pos;

		//NxQuat q;
		//NxMat33 orient;
		//orient.fromQuat(q);

		printf("adding x offset %f\n", tPos.x);
		NxVec3 segPos = NxVec3(i*segLength,boxStartHeight + 0*(i%2 * 0.2),0);	
		ori.rotate(segPos);

		actorDesc.globalPose.t	= tPos + segPos;	


		assert(actorDesc.isValid());

		pActor = gScene->createActor(actorDesc);	
		assert(pActor);

		segments[i] = pActor;

		//gScene->setActorPairFlags(*segments[i],*groundPlane,NX_NOTIFY_ON_START_TOUCH|NX_NOTIFY_ON_TOUCH|NX_NOTIFY_ON_END_TOUCH);

		//pActor->setAngularDamping(0);
		//pActor->setMaxAngularVelocity(700);
	}

	float per_seg = 2*NxPi / (numSegs / 5.0);
	
	for ( int i = 0 ; i < numSegs-1 ; i++ ) {
		NxActor *seg1 = segments[i];
		NxActor *seg2 = segments[i+1];

		//# joint positions
		//armJointPos = ogre.Vector3(self.segLength/2.0 + (i-1)*self.segLength, 0.0, 0.0)
		//actJointPos = R*armJointPos
		//anch = actJointPos + pos
		//NxVec3 segPos = NxVec3(i*segLength,boxStartHeight + 0*(i%2 * 0.2),0);	
		//ori.rotate(segPos);
		//actorDesc.globalPose.t	= tPos + segPos;	

		NxVec3 armJointPos = NxVec3(segLength/2.0+i*segLength,boxStartHeight,0);
		ori.rotate(armJointPos);


		NxVec3 globalAnchor = tPos + armJointPos;
		NxVec3 globalAxis = NxVec3(0,-1,0);
		NxRevoluteJoint* revJoint = CreateRevoluteJoint(seg1, seg2, globalAnchor, globalAxis);
		joints[i] = revJoint;


		//float val = NxPi/180.0*70.0*cos(i*per_seg);

		//targetAngle[i] = val;
		targetAngle[i] = 0.0;
		//revJoint.setBreakable(False);
	}

	return pActor;
}


NxRevoluteJoint* NxSnake::CreateRevoluteJoint(NxActor* a0, NxActor* a1, NxVec3 globalAnchor, NxVec3 globalAxis)
{
    NxRevoluteJointDesc revDesc;

    revDesc.actor[0] = a0;
    revDesc.actor[1] = a1;
    revDesc.setGlobalAnchor(globalAnchor);
    revDesc.setGlobalAxis(globalAxis);

    //revDesc.jointFlags |= NX_JF_COLLISION_ENABLED;

	revDesc.projectionMode = NX_JPM_POINT_MINDIST;
	revDesc.projectionDistance = 0.01f;
	revDesc.projectionAngle = 0.0872f;	//about 5 degrees in radians.

	// Force at which the joint breaks
	revDesc.maxForce = 1e20;
	revDesc.maxTorque = 1e20;
	
	// Enable and set up joint limits
    revDesc.flags |= NX_RJF_LIMIT_ENABLED;
    revDesc.limit.high.value = 90.0*NxPi/180.0;
    revDesc.limit.high.restitution = 0.9f;
    revDesc.limit.low.value = -90.0*NxPi/180.0;
    revDesc.limit.low.restitution = 0.9f;

	revDesc.flags |= NX_RJF_MOTOR_ENABLED;

	NxMotorDesc motorDesc;
    motorDesc.velTarget = 0;
    motorDesc.maxForce = 0.0005;

    //motorDesc.freeSpin = true;

    revDesc.motor = motorDesc;

    return (NxRevoluteJoint*)gScene->createJoint(revDesc);
}


void NxSnake::InitNx()
{

	// Create the physics SDK
    gPhysicsSDK = NxCreatePhysicsSDK(NX_PHYSICS_SDK_VERSION);
    if (!gPhysicsSDK)  return;

	// Set the physics parameters
	gPhysicsSDK->setParameter(NX_SKIN_WIDTH, 0.01);
	//gPhysicsSDK->setParameter(NX_CONTINUOUS_CD, 1);
	//gPhysicsSDK->setParameter(NX_CCD_EPSILON, 0.001f);
	//gPhysicsSDK->setParameter(NX_SKIN_WIDTH, 0.0005f);

	// Set the debug visualization parameters
	gPhysicsSDK->setParameter(NX_VISUALIZATION_SCALE, 1);
	gPhysicsSDK->setParameter(NX_VISUALIZE_COLLISION_SHAPES, 1);
	gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES, 1);

	gPhysicsSDK->setParameter(NX_BOUNCE_THRESHOLD, -10);


    // Create the scene
    NxSceneDesc sceneDesc;
 	sceneDesc.simType				= NX_SIMULATION_HW;
    sceneDesc.gravity               = gDefaultGravity;
	sceneDesc.userContactReport		= &gMyContactReport; 
	gScene = gPhysicsSDK->createScene(sceneDesc);	
	if(!gScene)
	{ 
		printf("\nHardware failed\n");
		sceneDesc.simType			= NX_SIMULATION_SW; 
		gScene = gPhysicsSDK->createScene(sceneDesc);  
		if(!gScene) return;
	}


	// Create the default material
	NxMaterial* defaultMaterial = gScene->getMaterialFromIndex(0); 
	defaultMaterial->setRestitution(0.01);
	defaultMaterial->setStaticFriction(friction);
	defaultMaterial->setDynamicFriction(0.01);
	//defaultMaterial->setStaticFriction(0.7);
	//defaultMaterial->setDynamicFriction(0.7);


	NxMaterialDesc materialDesc;
	materialDesc.restitution = 0.01;
	materialDesc.dynamicFriction = 0.01;
	materialDesc.staticFriction = 0.01;
	materialDesc.frictionCombineMode = NX_CM_MIN;
	pIndex = gScene->createMaterial(materialDesc)->getMaterialIndex();
	printf("New material = %d\n", pIndex);

	gMyContactReport.setIndices(dIndex,pIndex); 




	// Create the objects in the scene
	groundPlane		= CreateGroundPlane();

	//gMyContactReport.setPlane(groundPlane);

	CreateSnake();
	
	NxVec3 points[4] = {NxVec3(-10.0, 0.0, -0.2),
						NxVec3(10.0, 0.0, -0.2),
						NxVec3(10.0, 0.0, 0.2),
						NxVec3(-10.0, 0.0, 0.2)
						};

	//CreateTriangleMeshWall(4, points);

	// Get the current time
	//getElapsedTime();

	// Start the first frame of the simulation
	if (gScene)  StartPhysics();
}

void NxSnake::ReleaseNx()
{
    if (gScene)
	{
		GetPhysicsResults();  // Make sure to fetchResults() before shutting down
		gPhysicsSDK->releaseScene(*gScene);
	}
	if (gPhysicsSDK)  gPhysicsSDK->release();
}

void NxSnake::ResetNx()
{
	ReleaseNx();
	InitNx();
}

void NxSnake::StartPhysics()
{
	// Update the time step
	//NxReal gDeltaTime = getElapsedTime();

	// Start collision and dynamics for delta time since the last frame
	gScene->simulate(0.1);
	//gScene->simulate(gDeltaTime);

	gScene->flushStream();

}

void NxSnake::GetPhysicsResults()
{
	// Get results from gScene->simulate(gDeltaTime)
	while (!gScene->fetchResults(NX_RIGID_BODY_FINISHED, false));
}

void NxSnake::setServo(int i, double angle) {
	targetAngle[i] = angle;
}

double NxSnake::getServo(int i) {
	return joints[i]->getAngle();
}

double NxSnake::getServoCmd(int i) {
	return targetAngle[i];
}

void NxSnake::setJointTorque(int i, double torque) {

	NxRevoluteJointDesc revDesc;
	joints[i]->saveToDesc(revDesc);
	revDesc.motor.maxForce = torque;
	joints[i]->loadFromDesc(revDesc);

}

double NxSnake::getJointTorque(int i) {
	NxMotorDesc motorDesc;
	joints[i]->getMotor(motorDesc);
	return motorDesc.maxForce;
}

//void NxSnake::getGlobalPose(int i, NxF32 *d) {
//	segments[i]->getGlobalPose().getRowMajor44(d);
//}

void NxSnake::getGlobalPosition(int i, double *x, double *y, double *z) {

	NxVec3 result = segments[i]->getGlobalPosition();
	*x = result.x;
	*y = result.y;
	*z = result.z;
}

void NxSnake::getGlobalOrientationQuat(int i, double *x, double *y, double *z, double *w) {
	NxQuat result = segments[i]->getGlobalOrientationQuat();

	*x = result.x;
	*y = result.y;
	*z = result.z;
	*w = result.w;
}

void NxSnake::savePose() {
	for ( int i = 0 ; i < numSegs ; i++ ) {
		NxVec3 pos = segments[i]->getGlobalPosition();
		NxQuat quat = segments[i]->getGlobalOrientationQuat();

		savePoses[i] = pos;
		saveOrientations[i] = quat;
	}
}

void NxSnake::restorePose() {
	for ( int i = 0 ; i < numSegs ; i++ ) {
		NxVec3 pos = savePoses[i];
		NxQuat quat = saveOrientations[i];

		segments[i]->setGlobalPosition(pos);
		segments[i]->setGlobalOrientationQuat(quat);
	}
}


//NxVec3 NxSnake::getGlobalPosition(int i) {
//	return segments[i]->getGlobalPosition();
//}

//NxQuat NxSnake::getGlobalOrientationQuat(int i) {
//	return segments[i]->getGlobalOrientationQuat();
//}

//Vector3 Position = Convert(m_pActor->getGlobalPosition());
//Quaternion Orientation = Convert(m_pActor->getGlobalOrientationQuat());
//m_pNode->setPosition(Position);
//m_pNode->setOrientation(Orientation);

/*
NxOgre::Matrix44 NxSnake::getGlobalPose(NxActor* actor) {

	Matrix44 matrix;
	actor->getGlobalPose().getRowMajor44(matrix.ptr());
	return matrix;

}
*/

