// ===============================================================================
//						  NVIDIA PHYSX SDK TRAINING PROGRAMS
//						     LESSON 101 : PRIMARY SHAPE
//
//						    Written by QA BJ, 6-2-2008
// ===============================================================================

#include "Lesson101.h"
#include "Timing.h"

#include "NxSnake.h"

// Physics SDK globals
NxPhysicsSDK*     gPhysicsSDK = NULL;
NxScene*          gScene = NULL;
NxVec3            gDefaultGravity(0,-9.8,0);

// User report globals
DebugRenderer     gDebugRenderer;

// HUD globals
HUD hud;

// Display globals
int gMainHandle;
int mx = 0;
int my = 0;

// Camera globals
float	gCameraAspectRatio = 1.0f;
NxVec3	gCameraPos(0,5,-15);
NxVec3	gCameraForward(0,0,1);
NxVec3	gCameraRight(-1,0,0);
const NxReal gCameraSpeed = 10;

// Force globals
NxVec3	gForceVec(0,0,0);
NxReal	gForceStrength	= 20000;
bool	bForceMode		= true;

// Keyboard globals
#define MAX_KEYS 256
bool gKeys[MAX_KEYS];

// Simulation globals
NxReal	gDeltaTime			= 1.0/60.0;
bool	bHardwareScene		= false;
bool	bPause				= false;
bool	bShadows			= true;
bool	bDebugWireframeMode = false;

// Actor globals
NxActor* groundPlane = NULL;

// Focus actor
NxActor* gSelectedActor = NULL;

#define NUMSEGS 40

NxActor *segments[NUMSEGS];
NxRevoluteJoint *joints[NUMSEGS-1];
NxReal targetAngle[NUMSEGS-1];

void PrintControls()
{
	printf("\n Flight Controls:\n ----------------\n w = forward, s = back\n a = strafe left, d = strafe right\n q = up, z = down\n");
    printf("\n Force Controls:\n ---------------\n i = +z, k = -z\n j = +x, l = -x\n u = +y, m = -y\n");
	printf("\n Miscellaneous:\n --------------\n p   = Pause\n b   = Toggle Debug Wireframe Mode\n x   = Toggle Shadows\n r   = Select Actor\n F10 = Reset scene\n");
}

NxVec3 ApplyForceToActor(NxActor* actor, const NxVec3& forceDir, const NxReal forceStrength)
{
	NxVec3 forceVec = forceStrength*forceDir*gDeltaTime;
	actor->addForce(forceVec);
	return forceVec;
}

void ProcessCameraKeys()
{
	NxReal deltaTime;

    if (bPause) 
	{
		deltaTime = 0.0005;
	}
	else 
	{
		deltaTime = gDeltaTime;
	}
		
	// Process camera keys
	for (int i = 0; i < MAX_KEYS; i++)
	{	
		if (!gKeys[i])  { continue; }

		switch (i)
		{
			// Camera controls
			case 'w':{ gCameraPos += gCameraForward*gCameraSpeed*deltaTime; break; }
			case 's':{ gCameraPos -= gCameraForward*gCameraSpeed*deltaTime; break; }
			case 'a':{ gCameraPos -= gCameraRight*gCameraSpeed*deltaTime;	break; }
			case 'd':{ gCameraPos += gCameraRight*gCameraSpeed*deltaTime;	break; }
			case 'z':{ gCameraPos -= NxVec3(0,1,0)*gCameraSpeed*deltaTime;	break; }
			case 'q':{ gCameraPos += NxVec3(0,1,0)*gCameraSpeed*deltaTime;	break; }
		}
	}
}

void SetupCamera()
{
	gCameraAspectRatio = (float)glutGet(GLUT_WINDOW_WIDTH) / (float)glutGet(GLUT_WINDOW_HEIGHT);
	
	// Setup camera
	glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0f, gCameraAspectRatio, 1.0f, 10000.0f);
	gluLookAt(gCameraPos.x,gCameraPos.y,gCameraPos.z,gCameraPos.x + gCameraForward.x, gCameraPos.y + gCameraForward.y, gCameraPos.z + gCameraForward.z, 0.0f, 1.0f, 0.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void RenderActors(bool shadows)
{
    // Render all the actors in the scene
    NxU32 nbActors = gScene->getNbActors();
    NxActor** actors = gScene->getActors();
    while (nbActors--)
    {
        NxActor* actor = *actors++;
        DrawActor(actor, gSelectedActor, false);

        // Handle shadows
        if (shadows)
        {
			DrawActorShadow(actor, false);
        }
    }
}

void DrawForce(NxActor* actor, NxVec3& forceVec, const NxVec3& color)
{
	// Draw only if the force is large enough
	NxReal force = forceVec.magnitude();
	if (force < 0.1)  return;

	forceVec = 3*forceVec/force;

	NxVec3 pos = actor->getCMassGlobalPosition();
	DrawArrow(pos, pos + forceVec, color);
}

bool IsSelectable(NxActor* actor)
{
   NxShape*const* shapes = gSelectedActor->getShapes();
   NxU32 nShapes = gSelectedActor->getNbShapes();
   while (nShapes--)
   {
       if (shapes[nShapes]->getFlag(NX_TRIGGER_ENABLE)) 
       {           
           return false;
       }
   }

   if(!actor->isDynamic())
	   return false;

   if (actor == groundPlane)
       return false;

   return true;
}

void SelectNextActor()
{
   NxU32 nbActors = gScene->getNbActors();
   NxActor** actors = gScene->getActors();
   for(NxU32 i = 0; i < nbActors; i++)
   {
       if (actors[i] == gSelectedActor)
       {
           NxU32 j = 1;
           gSelectedActor = actors[(i+j)%nbActors];
           while (!IsSelectable(gSelectedActor))
           {
               j++;
               gSelectedActor = actors[(i+j)%nbActors];
           }
           break;
       }
   }
}

void ProcessForceKeys()
{
	// Process force keys
	for (int i = 0; i < MAX_KEYS; i++)
	{	
		if (!gKeys[i])  { continue; }

		switch (i)
		{
			// Force controls
			case 'i': { gForceVec = ApplyForceToActor(gSelectedActor,NxVec3(0,0,1),gForceStrength);		break; }
			case 'k': { gForceVec = ApplyForceToActor(gSelectedActor,NxVec3(0,0,-1),gForceStrength);	break; }
			case 'j': { gForceVec = ApplyForceToActor(gSelectedActor,NxVec3(1,0,0),gForceStrength);		break; }
			case 'l': { gForceVec = ApplyForceToActor(gSelectedActor,NxVec3(-1,0,0),gForceStrength);	break; }
			case 'u': { gForceVec = ApplyForceToActor(gSelectedActor,NxVec3(0,1,0),gForceStrength);		break; }
			case 'm': { gForceVec = ApplyForceToActor(gSelectedActor,NxVec3(0,-1,0),gForceStrength);	break; }

			case '+': { for ( int j = 0 ; j < NUMSEGS-1; j++ ) targetAngle[j] += NxPi/180.0; break; }
			case '-': { for ( int j = 0 ; j < NUMSEGS-1; j++ ) targetAngle[j] -= NxPi/180.0; break; }

		    // Return box to (0,5,0)
			case 't': 
			{ 
				if (gSelectedActor) 
				{
					gSelectedActor->setGlobalPosition(NxVec3(0,5,0)); 
					gScene->flushCaches();
				}
				break; 
			}
		}
	}
}

void ProcessInputs()
{

    ProcessForceKeys();

    // Show debug wireframes
	if (bDebugWireframeMode)
	{
		if (gScene)  gDebugRenderer.renderData(*gScene->getDebugRenderable());
	}
}

/*
void commandServo(Phi, step, goalPhi, initPhi, tolerance, PGain, IGain, DGain, lastErr, maxVelocity):	

	double errSum = 0.0
	double vel = 0.0
	double pterm, iterm, dterm, new_last_err
	int isDone = 0
	
	# input arguments
	double pgain = PGain
	double igain = IGain
	double dgain = DGain
	double phi = Phi
	double dt = step
	double goal_phi = goalPhi
	double init_phi = initPhi
	double tol = tolerance
	double last_err = lastErr
	double maxVel = maxVelocity
	
	# compute the servo joint error
	double err = goal_phi - (phi + init_phi)

	if fabs(err) > tol:
		errSum += dt*err
		pterm = pgain*err
		iterm = igain*errSum
		dterm = dgain*(err-last_err)/dt
		new_last_err = err
		vel = lim(pterm+iterm+dterm, -maxVel, maxVel)

	else:
		isDone=True
		errSum=0.0
		new_last_err=0.0
		
	return vel;
}
*/

void RunPID() {

	float pgain = 10;
	float dgain = .001;

	for ( int i = 0 ; i < NUMSEGS-1 ; i++ ) {
		NxRevoluteJointDesc revDesc;
		joints[i]->saveToDesc(revDesc);
	
		NxReal init_phi = 0.0;
		NxReal phi = joints[i]->getAngle();

		NxReal err = targetAngle[i] - (phi + init_phi);

		NxReal vel = pgain*err + dgain*joints[i]->getVelocity();
		revDesc.motor.velTarget = vel;
		//printf("%f %f %f %f\n", pgain*err, dgain*joints[i]->getVelocity(), phi, vel);
		//printf("Angle: %f\n", phi);
		//printf("Velocity = %f\n", vel);

		/*
		if (revDesc.motor.velTarget == -10)
		   revDesc.motor.velTarget = 10;
		else
		   revDesc.motor.velTarget = -10;
		*/

		joints[i]->loadFromDesc(revDesc);
	}
}


void RenderCallback()
{
    // Clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    ProcessCameraKeys();
	SetupCamera();

    if (gScene && !bPause)
	{
		GetPhysicsResults();
        ProcessInputs();

		targetAngle[20] = NxPi * 90.0;
		//ns->setServo(20,NxPi * 90.0);

		RunPID();

		StartPhysics();


	}

    // Display scene
 	RenderActors(bShadows);

	DrawForce(gSelectedActor, gForceVec, NxVec3(1,1,0));
	gForceVec = NxVec3(0,0,0);

	// Render the HUD
	hud.Render();

    glFlush();
    glutSwapBuffers();
}

void ReshapeCallback(int width, int height)
{
    glViewport(0, 0, width, height);
    gCameraAspectRatio = float(width)/float(height);
}

void IdleCallback()
{
    glutPostRedisplay();
}

void KeyboardCallback(unsigned char key, int x, int y)
{
	gKeys[key] = true;

	switch (key)
	{
		case 'r':	{ SelectNextActor(); break; }
		default:	{ break; }
	}
}

void KeyboardUpCallback(unsigned char key, int x, int y)
{
	gKeys[key] = false;

	switch (key)
	{
		case 'p': 
		{ 
			bPause = !bPause; 
			if (bPause)
				hud.SetDisplayString(0, "Paused - Hit \"p\" to Unpause", 0.3f, 0.55f);
			else
				hud.SetDisplayString(0, "", 0.0f, 0.0f);	
			getElapsedTime(); 
			break; 
		}
		case 'x': { bShadows = !bShadows; break; }
		case 'b': { bDebugWireframeMode = !bDebugWireframeMode; break; }
		case 27 : { exit(0); break; }
		default : { break; }
	}
}

void SpecialCallback(int key, int x, int y)
{
	switch (key)
    {
		// Reset PhysX
		case GLUT_KEY_F10: ResetNx(); return; 
	}
}

void MouseCallback(int button, int state, int x, int y)
{
    mx = x;
    my = y;
}

void MotionCallback(int x, int y)
{
    int dx = mx - x;
    int dy = my - y;
    
    gCameraForward.normalize();
    gCameraRight.cross(gCameraForward,NxVec3(0,1,0));

    NxQuat qx(NxPiF32 * dx * 20 / 180.0f, NxVec3(0,1,0));
    qx.rotate(gCameraForward);
    NxQuat qy(NxPiF32 * dy * 20 / 180.0f, gCameraRight);
    qy.rotate(gCameraForward);

    mx = x;
    my = y;
}

void ExitCallback()
{
	ReleaseNx();
}

void InitGlut(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(512, 512);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    gMainHandle = glutCreateWindow("Lesson 101: Primary Shape");
    glutSetWindow(gMainHandle);
    glutDisplayFunc(RenderCallback);
    glutReshapeFunc(ReshapeCallback);
    glutIdleFunc(IdleCallback);
    glutKeyboardFunc(KeyboardCallback);
    glutKeyboardUpFunc(KeyboardUpCallback);
	glutSpecialFunc(SpecialCallback);
    glutMouseFunc(MouseCallback);
    glutMotionFunc(MotionCallback);
	MotionCallback(0,0);
	atexit(ExitCallback);

    // Setup default render states
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_CULL_FACE);

    // Setup lighting
    glEnable(GL_LIGHTING);
    float AmbientColor[]    = { 0.0f, 0.1f, 0.2f, 0.0f };         glLightfv(GL_LIGHT0, GL_AMBIENT, AmbientColor);
    float DiffuseColor[]    = { 0.2f, 0.2f, 0.2f, 0.0f };         glLightfv(GL_LIGHT0, GL_DIFFUSE, DiffuseColor);
    float SpecularColor[]   = { 0.5f, 0.5f, 0.5f, 0.0f };         glLightfv(GL_LIGHT0, GL_SPECULAR, SpecularColor);
    float Position[]        = { 100.0f, 100.0f, -400.0f, 1.0f };  glLightfv(GL_LIGHT0, GL_POSITION, Position);
    glEnable(GL_LIGHT0);
}

NxActor* CreateGroundPlane()
{
    // Create a plane with default descriptor
    NxPlaneShapeDesc planeDesc;
    NxActorDesc actorDesc;
    actorDesc.shapes.pushBack(&planeDesc);
    return gScene->createActor(actorDesc);
}

NxCCDSkeleton* CreateCCDSkeleton(float sizex, float sizey, float sizez)
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

NxActor* CreateTriangleMeshWall(int size, NxVec3 *points)
{

	NxTriangleMeshDesc* triangleMeshDesc = NULL;
	NxVec3 boxDim(1.0f, 1.0f, 1.0f);

	NxVec3 *verts = new NxVec3[ 2 * size];

	for ( int i = 0 ; i < size ; i++ ) {
		verts[2*i] = NxVec3(points[i].x, -1.0, points[i].z);
		verts[2*i+1] = NxVec3(points[i].x, 1.0, points[i].z);
	}

	//indices1 += [2*i+2, 2*i+1, 2*i]
	//indices1 += [2*i+1, 2*i+2, 2*i+3]

	NxU32 *indices = new NxU32[ 6 * (size-1)];

	for ( int i = 0 ; i < size-1 ; i++ ) {
		indices[6*i + 0] = 2*i+2;
		indices[6*i + 1] = 2*i+1;
		indices[6*i + 2] = 2*i;
		indices[6*i + 3] = 2*i+1;
		indices[6*i + 4] = 2*i+2;
		indices[6*i + 5] = 2*i+3;
	}

	/*
	// Supply hull
	 NxVec3 verts[4] = {NxVec3(-10.0, -1.0, 0.0), 
						NxVec3(-10.0, 1.0, 0.0), 
						NxVec3(10.0, -1.0, 0.0), 
						NxVec3(10.0, 1.0, 0.0), 
						};
	 */

	// Triangles is 12*3
	 //NxU32 indices[2*3] = {2*0+2, 2*0+1, 2*0,        
	//						2*0+1, 2*0+2, 2*0+3  
	//						};

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

	bool status = NxCookTriangleMesh(*triangleMeshDesc, buf);
	NxTriangleMesh* pMesh;
	if (status)
	{
		pMesh = gPhysicsSDK->createTriangleMesh(MemoryReadBuffer(buf.data));
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

NxActor* CreateWall()
{

	// Set the box starting height to 3.5m so box starts off falling onto the ground
	NxReal boxStartHeight = 0.0; 

	// Add a single-shape actor to the scene
	NxActorDesc actorDesc;

	// The actor has one shape, a box, 1m on a side
	NxBoxShapeDesc boxDesc;
	boxDesc.dimensions.set(10.0,0.2,0.1);
	boxDesc.localPose.t = NxVec3(0, 0.2, 0);

	actorDesc.shapes.pushBack(&boxDesc);    

	actorDesc.globalPose.t	= NxVec3(0,0,0.3);	
	assert(actorDesc.isValid());
	NxActor *pActor = gScene->createActor(actorDesc);	
	assert(pActor);


	NxBoxShapeDesc boxDesc2;
	boxDesc2.dimensions.set(10.0,0.2,0.1);
	boxDesc2.localPose.t = NxVec3(0, 0.2, 0);

	actorDesc.shapes.pushBack(&boxDesc2);    

	actorDesc.globalPose.t	= NxVec3(0,0,-0.3);	
	assert(actorDesc.isValid());
	NxActor *pActor2 = gScene->createActor(actorDesc);	
	assert(pActor2);

	return pActor;


	/*
	// Create a plane with default descriptor
    NxPlaneShapeDesc wall1Desc;
	wall1Desc.normal = NxVec3(0.0,0.0,1.0);
    NxActorDesc actorDesc1;
    actorDesc1.shapes.pushBack(&wall1Desc);
	actorDesc1.globalPose.t = NxVec3(0,0,-3.5);
    gScene->createActor(actorDesc1);

    NxPlaneShapeDesc wall2Desc;
	wall2Desc.normal = NxVec3(0.0,0.0,-1.0);
    NxActorDesc actorDesc2;
    actorDesc2.shapes.pushBack(&wall2Desc);
	actorDesc2.globalPose.t = NxVec3(0,0,3.5);
    return gScene->createActor(actorDesc2);

	*/

	/*
	
	points = []
	wall1Points = []
	for i in range(len(points)):
		wall1Points.append(ogre.Vector3(points[i][0], -100.0, points[i][1]))
		wall1Points.append(ogre.Vector3(points[i][0], 1.5, points[i][1]))
	for i in range(len(points) - 1):
		indices1 += [2*i+2, 2*i+1, 2*i]
		indices1 += [2*i+1, 2*i+2, 2*i+3]
	# create the ODE geometry
	self.trimeshes.append(OgreOde.makeTriangleMeshGeometry(wall1Points, len(wall1Points), indices1, len(indices1), self._world, self._space))
	*/

}

NxActor* CreateBox()
{
	// Set the box starting height to 3.5m so box starts off falling onto the ground
	NxReal boxStartHeight = 0.5; 

	// Add a single-shape actor to the scene
	NxActorDesc actorDesc;
	NxBodyDesc bodyDesc;

	// The actor has one shape, a box, 1m on a side
	NxBoxShapeDesc boxDesc;
	boxDesc.dimensions.set(0.5,0.5,0.5);
	boxDesc.localPose.t = NxVec3(0, 0, 0);

	boxDesc.ccdSkeleton = CreateCCDSkeleton(0.5f, 0.5f, 0.5f);
    if (1)  boxDesc.shapeFlags |= NX_SF_DYNAMIC_DYNAMIC_CCD;  

	actorDesc.shapes.pushBack(&boxDesc);    

	actorDesc.body			= &bodyDesc;
	actorDesc.density		= 10.0f;
	actorDesc.globalPose.t	= NxVec3(0,boxStartHeight,0);	
	assert(actorDesc.isValid());
	NxActor *pActor = gScene->createActor(actorDesc);	
	assert(pActor);

	// //create actor with no shapes
	//NxShape* const *shape = pActor->getShapes();
	//NxBoxShape *boxShape = shape[0]->isBox();
	//assert(boxShape);
	//pActor->releaseShape(*boxShape);

	return pActor;
}

NxActor* CreateSnake()
{

	// Set the box starting height to 3.5m so box starts off falling onto the ground
	NxReal boxStartHeight = 0.0; 

	// Add a single-shape actor to the scene
	NxActorDesc actorDesc;
	NxBodyDesc bodyDesc;

	//box1 = CreateBox(NxVec3(0,5,0), NxVec3(0.5,2,1), 10);
	//box1->raiseBodyFlag(NX_BF_KINEMATIC);
	//box2 = CreateBox(NxVec3(0,1,0), NxVec3(0.5,2,1), 10);

	NxActor *pActor;

	for ( int i = 0 ; i < NUMSEGS ; i++ ) {

		// The actor has one shape, a box, 1m on a side
		NxBoxShapeDesc boxDesc;
		boxDesc.dimensions.set(0.075,0.1,0.075);
		boxDesc.localPose.t = NxVec3(0, -0.1, 0);

		boxDesc.ccdSkeleton = CreateCCDSkeleton(0.075f, 0.1f, 0.075f);
	    if (1)  boxDesc.shapeFlags |= NX_SF_DYNAMIC_DYNAMIC_CCD;  

		actorDesc.shapes.pushBack(&boxDesc);

		actorDesc.body			= &bodyDesc;
		actorDesc.density		= 1.0f;
		actorDesc.globalPose.t	= NxVec3(i*0.15,boxStartHeight + 0*(i%2 * 0.2),0);	
		assert(actorDesc.isValid());
		pActor = gScene->createActor(actorDesc);	
		assert(pActor);

		segments[i] = pActor;

		//pActor->setAngularDamping(0);
		//pActor->setMaxAngularVelocity(700);
	}

	float per_seg = 2*NxPi / (NUMSEGS / 5.0);
	
	for ( int i = 0 ; i < NUMSEGS-1 ; i++ ) {
		NxActor *seg1 = segments[i];
		NxActor *seg2 = segments[i+1];

		NxVec3 globalAnchor = NxVec3(0.075+i*0.15,boxStartHeight,0);
		NxVec3 globalAxis = NxVec3(0,1,0);
		NxRevoluteJoint* revJoint = CreateRevoluteJoint(seg1, seg2, globalAnchor, globalAxis);
		joints[i] = revJoint;


		float val = NxPi/180.0*70.0*cos(i*per_seg);

		targetAngle[i] = val;
		//targetAngle[i] = 0.0;
		//revJoint.setBreakable(False);
	}

	// //create actor with no shapes
	//NxShape* const *shape = pActor->getShapes();
	//NxBoxShape *boxShape = shape[0]->isBox();
	//assert(boxShape);
	//pActor->releaseShape(*boxShape);

	return pActor;
}


NxRevoluteJoint* CreateRevoluteJoint(NxActor* a0, NxActor* a1, NxVec3 globalAnchor, NxVec3 globalAxis)
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

NxActor* CreateSphere()
{
	// Set the sphere starting height to 3.5m so box starts off falling onto the ground
	NxReal sphereStartHeight = 3.5; 

	// Add a single-shape actor to the scene
	NxActorDesc actorDesc;
	NxBodyDesc bodyDesc;

	// The actor has one shape, a sphere, 1m on radius
	NxSphereShapeDesc sphereDesc;
	sphereDesc.radius		= 0.65f;
	sphereDesc.localPose.t	= NxVec3(0, 0, 0);

	actorDesc.shapes.pushBack(&sphereDesc);
	actorDesc.body			= &bodyDesc;
	actorDesc.density		= 10.0f;
	actorDesc.globalPose.t	= NxVec3(3.0f,sphereStartHeight,0);	
	return gScene->createActor(actorDesc);
}


NxActor* CreateCapsule()
{
	// Set the capsule starting height to 3.5m so box starts off falling onto the ground
	NxReal capsuleStartHeight = 3.5; 

	// Add a single-shape actor to the scene
	NxActorDesc actorDesc;
	NxBodyDesc bodyDesc;

	// The actor has one shape, a sphere, 1m on radius
	NxCapsuleShapeDesc capsuleDesc;
	capsuleDesc.radius		= 0.55f;
	capsuleDesc.height		= 0.75f;
	capsuleDesc.localPose.t = NxVec3(0, 0, 0);
	
	//Rotate capsule shape
	NxQuat quat(45, NxVec3(0, 0, 1));
	NxMat33 m33(quat);
	capsuleDesc.localPose.M = m33;

	actorDesc.shapes.pushBack(&capsuleDesc);
	actorDesc.body			= &bodyDesc;
	actorDesc.density		= 10.0f;
	actorDesc.globalPose.t	= NxVec3(6.0f,capsuleStartHeight,0);

	////Rotate actor
	//NxQuat quat1(45, NxVec3(1, 0, 0));
	//NxMat33 m331(quat1);
	//actorDesc.globalPose.M = m331;

	return gScene->createActor(actorDesc);
}

void InitializeHUD()
{
	bHardwareScene = (gScene->getSimType() == NX_SIMULATION_HW);

	hud.Clear();

	//// Add hardware/software to HUD
	//if (bHardwareScene)
	//    hud.AddDisplayString("Hardware Scene", 0.74f, 0.92f);
	//else
	//	hud.AddDisplayString("Software Scene", 0.74f, 0.92f);

	// Add pause to HUD
	if (bPause)  
		hud.AddDisplayString("Paused - Hit \"p\" to Unpause", 0.3f, 0.55f);
	else
		hud.AddDisplayString("", 0.0f, 0.0f);
}

void InitNx()
{
	// Initialize camera parameters
	gCameraAspectRatio	= 1.0f;
	gCameraPos			= NxVec3(3,6,0);
	gCameraForward		= NxVec3(0,-0.9,0.1);
	gCameraRight		= NxVec3(-1,0,0);
	
	// Create the physics SDK
    gPhysicsSDK = NxCreatePhysicsSDK(NX_PHYSICS_SDK_VERSION);
    if (!gPhysicsSDK)  return;

	// Set the physics parameters
	//gPhysicsSDK->setParameter(NX_SKIN_WIDTH, 0.0001);
	//gPhysicsSDK->setParameter(NX_CONTINUOUS_CD, 1);
	//gPhysicsSDK->setParameter(NX_CCD_EPSILON, 0.001f);
	gPhysicsSDK->setParameter(NX_SKIN_WIDTH, 0.0005f);

	// Set the debug visualization parameters
	gPhysicsSDK->setParameter(NX_VISUALIZATION_SCALE, 1);
	gPhysicsSDK->setParameter(NX_VISUALIZE_COLLISION_SHAPES, 1);
	gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES, 1);


    // Create the scene
    NxSceneDesc sceneDesc;
 	sceneDesc.simType				= NX_SIMULATION_HW;
    sceneDesc.gravity               = gDefaultGravity;
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
	defaultMaterial->setRestitution(0.5);
	defaultMaterial->setStaticFriction(0.1);
	defaultMaterial->setDynamicFriction(0.1);

	// Create the objects in the scene
	groundPlane		= CreateGroundPlane();

	//gSelectedActor = CreateSegment();
	//CreateSegment();
	//CreateSegment();
	//CreateSegment();
	//CreateSegment();
	CreateSnake();

	//CreateWall();
	
	 NxVec3 points[4] = {NxVec3(-10.0, 0.0, -0.2),
						NxVec3(10.0, 0.0, -0.2),
						NxVec3(10.0, 0.0, 0.2),
						NxVec3(-10.0, 0.0, 0.2)
						};

	CreateTriangleMeshWall(4, points);
	//CreateTriangleMeshWall(-0.3);
	//gSelectedActor	= CreateBox();
	//CreateSphere();
	//CreateCapsule();

	// Initialize HUD
	InitializeHUD();

	// Get the current time
	getElapsedTime();

	// Start the first frame of the simulation
	if (gScene)  StartPhysics();
}

void ReleaseNx()
{
    if (gScene)
	{
		GetPhysicsResults();  // Make sure to fetchResults() before shutting down
		gPhysicsSDK->releaseScene(*gScene);
	}
	if (gPhysicsSDK)  gPhysicsSDK->release();
}

void ResetNx()
{
	ReleaseNx();
	InitNx();
}

void StartPhysics()
{
	// Update the time step
	gDeltaTime = getElapsedTime();
	//printf("%f\n", gDeltaTime);

	// Fixed simulation step.... 
    //gScene->simulate(0.001);

	// Start collision and dynamics for delta time since the last frame
	gScene->simulate(0.01);
	//gScene->simulate(gDeltaTime);
	//printf("%f\n", gDeltaTime);
	gScene->flushStream();

	for ( int i = 0 ; i < NUMSEGS-1 ; i++ ) {
		printf("%f ", joints[i]->getAngle());
	}
	printf("\n");
}

void GetPhysicsResults()
{
	// Get results from gScene->simulate(gDeltaTime)
	while (!gScene->fetchResults(NX_RIGID_BODY_FINISHED, false));
}
/*
int main(int argc, char** argv)
{
	PrintControls();
    InitGlut(argc, argv);
    InitNx();
    glutMainLoop();
	ReleaseNx();
	return 0;
}
*/


int main(int argc, char** argv)
{

	NxQuat yRot(0.0, NxVec3(0.0,1.0,0.0));
	double quat_param[4] = {yRot.x, yRot.y, yRot.z, yRot.w};
	double pos[3] = {3.65, 0.04, 0.0};
	NxSnake *ns = new NxSnake(quat_param, pos, 40, 0.15, 0.2, 0.15);

	while(1) {
		ns->Step();

		//ns->setServo(20,NxPi * 90.0);
		//if (ns->getServo(20) == -10.0)
		//   ns->setServo(20,NxPi * 90.0);
		//else
		//   ns->setServo(20,-10.0);
		

		//for ( int i = 0 ; i < NUMSEGS-1 ; i++ ) {
		//	printf("%f ", ns->getServo(i));
		//}
		//printf("\n");

		// Homogeneous Transformation Matrix
		//NxF32 foo[16];
		//ns->getGlobalPose(9, foo);

		double x, y, z, w;

		ns->getGlobalPosition(9, &x, &y, &z);
		printf("%f %f %f\n", x, y, z);

		ns->getGlobalOrientationQuat(9, &x, &y, &z, &w);
		printf("%f %f %f %f\n", x, y, z, w);

		//return 1;
		//printf("%f %f %f %f\n", foo[0], foo[1], foo[2], foo[3]);
		//printf("%f %f %f %f\n", foo[4], foo[5], foo[6], foo[7]);
		//printf("%f %f %f %f\n", foo[8], foo[9], foo[10], foo[11]);
		//printf("%f %f %f %f\n\n", foo[12], foo[13], foo[14], foo[15]);




		gDeltaTime = getElapsedTime();
		//printf("%f\n", gDeltaTime);
	}

	delete ns;

	return 0;

}

