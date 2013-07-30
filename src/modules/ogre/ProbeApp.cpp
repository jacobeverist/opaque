/*
-----------------------------------------------------------------------------
Filename:    ProbeApp.cpp
-----------------------------------------------------------------------------

This source file is part of the
   ___                 __    __ _ _    _
  /___\__ _ _ __ ___  / / /\ \ (_) | _(_)
 //  // _` | '__/ _ \ \ \/  \/ / | |/ / |
/ \_// (_| | | |  __/  \  /\  /| |   <| |
\___/ \__, |_|  \___|   \/  \/ |_|_|\_\_|
      |___/
      Tutorial Framework
      http://www.ogre3d.org/tikiwiki/
-----------------------------------------------------------------------------
*/
#include "ProbeApp.h"
#include <string>
#include <iostream>
#include <sstream>

//-------------------------------------------------------------------------------------
ProbeApp::ProbeApp(int numSegs, double segLength, double segHeight, double segWidth) {

	this->numSegs = numSegs;
	this->segLength = segLength;
	this->segHeight = segHeight;
	this->segWidth = segWidth;
	
	this->numWalls = 0;

	// assume no more than 100 walls
	this->wallPoints = new double*[100];
	this->numWallPoints = new int[100];

	this->probeNodes = new Ogre::SceneNode*[numSegs];
	this->probeMaterials = new Ogre::Material*[numSegs];
	this->probeEntities = new Ogre::Entity*[numSegs];

	#ifdef _DEBUG
    mResourcesCfg = "resources_d.cfg";
    mPluginsCfg = "plugins_d.cfg";
	#else
    mResourcesCfg = "resources.cfg";
    mPluginsCfg = "plugins.cfg";
	#endif

    bool returnVal = setup();

}
//-------------------------------------------------------------------------------------
ProbeApp::~ProbeApp(void)
{
}

//-------------------------------------------------------------------------------------
void ProbeApp::createScene(void)
{
    // create your scene here :)

    // Create the SkyBox
    //mSceneMgr->setSkyBox(true, "Examples/CloudyNoonSkyBox");


    // Set the scene's ambient light
    mSceneMgr->setAmbientLight(Ogre::ColourValue(0.8f, 0.8f, 0.8f));

    // Create a Light and set its position
    //Ogre::Light* light = mSceneMgr->createLight("MainLight");
    //light->setPosition(0.0f, 0.0f, 1.0f);
    //light->setPosition(20.0f, 80.0f, 50.0f);
	//light->setCastShadows(true);
	
	// Give us some sky
	mSceneMgr->setSkyBox(true,"kk3d/DesertVII", 5000, true); // Examples/SpaceSkyBox",5000,True)

	// Position and orient the camera
	mCamera->setPosition(Ogre::Vector3( 2, 10, 0 ));
	//#self.camera.lookAt(-7,0.5,0)
	mCamera->setNearClipDistance(0.1);


	// Create a default plane to act as the ground
	Ogre::StaticGeometry* s = mSceneMgr->createStaticGeometry("StaticFloor");
	s->setRegionDimensions(Ogre::Vector3(160.0, 100.0, 160.0));

	// Set the region origin so the center is at 0 world
	s->setOrigin(Ogre::Vector3(0,0,0));


	// Use a load of meshes to represent the floor
	int i = 0;
	for (int z = -80 ; z < 80 ; z += 20) {
		for (int x = -80 ; x < 80 ; x += 20) {

			std::stringstream ss;
			ss << i;
			std::string name("Plane_" + ss.str());

			i += 1;

			Ogre::Entity* entity = mSceneMgr->createEntity(name, "plane.mesh");
			//entity->setQueryFlags(STATIC_GEOMETRY_QUERY_MASK);
			entity->setCastShadows(false);
			s->addEntity(entity, Ogre::Vector3(x,0,z));

		}
	}

	// add coordinate system markers 

	Ogre::Vector3 xPnt(5.0,0.0,0.0);
	Ogre::Vector3 zPnt(0.0,0.0,5.0);

	Ogre::Entity* xEntity = mSceneMgr->createEntity("X_ent0", "Cube.mesh");
	xEntity->setCastShadows(false);
	xEntity->setMaterialName("Red");

	Ogre::Entity* zEntity = mSceneMgr->createEntity("Z_ent0", "Cube.mesh");
	zEntity->setCastShadows(false);
	zEntity->setMaterialName("Green");

	Ogre::Vector3 size = Ogre::Vector3(0.1,0.1,0.1);
	s->addEntity(xEntity, xPnt, Ogre::Quaternion::IDENTITY, size);
	s->addEntity(zEntity, zPnt, Ogre::Quaternion::IDENTITY, size);

	s->build();



	Ogre::Quaternion yRot(Ogre::Degree(180.0), Ogre::Vector3::UNIT_Y);
	Ogre::Vector3 pos(1.65, 0.04, 0.0);

	//setupMyWorld(Ogre::Quaternion::IDENTITY, Ogre::Vector3::ZERO);
	//setupMyWorld(yRot, pos);
	buildProbe();


	double WLEN2 = 7.0;
	double wall1[6] = {-14.0, -0.2, -4.0, -0.2, -4.0 + WLEN2*cos(M_PI/3), -0.2 - WLEN2*sin(M_PI/3)};
	double wall2[6] = {-4.0 + WLEN2*cos(M_PI/3), 0.2 + WLEN2*sin(M_PI/3), -4.0, 0.2 ,-14.0, 0.2};
	double wall5[4] = {wall2[2*2], wall2[2*2+1], wall1[0], wall1[1]};
	//w1 = wall1[2]
	//w2 = wall2[0]

	double wall3[8] = {wall1[2*2] + 0.4*cos(M_PI/6), wall1[2*2+1] + 0.4*sin(M_PI/6), 0.4*cos(M_PI/6) - 4, 0.0, wall2[0] + 0.4*cos(M_PI/6), wall2[1] - 0.4*sin(M_PI/6), wall2[0], wall2[1]};
	double wall4[4] = {wall1[2*2], wall1[2*2+1], wall1[2*2] + 0.4*cos(M_PI/3-M_PI/2), wall1[2*2+1] - 0.4*sin(M_PI/3-M_PI/2) };

	//wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
	//wall4 = [w1, [w1[0] + 0.4*cos(pi/3-pi/2), w1[1] - 0.4*sin(pi/3-pi/2)]]


	//walls = [wall1, wall2, wall3, wall4, wall5]


	//double wall1[4] = {-14.0,-0.2,14.0,-0.2};
	//double wall2[4] = {14.0,0.2,-14.0,0.2};

	//addWall(3,wall1);
	//addWall(3,wall2);
	//addWall(2,wall5);
	//addWall(4,wall3);
	//addWall(2,wall4);

	//createWalls();

    // Create a SceneNode and attach the Entity to it
    //Ogre::SceneNode* headNode = mSceneMgr->getRootSceneNode()->createChildSceneNode("HeadNode");
    //headNode->attachObject(ogreHead);

	mCamera->setPosition(Ogre::Vector3(pos[0],12.0,pos[2]));

	Ogre::Quaternion oriQuat(Ogre::Radian(-M_PI/2.0), Ogre::Vector3::UNIT_X);
	mCamera->setOrientation(oriQuat);


	mCameraMan->setTopSpeed(8.0);

}

void ProbeApp::addWall(int numPoints, double *points) {

	this->wallPoints[this->numWalls] = new double[numPoints*2];
	for ( int i = 0 ; i < numPoints ; i++ ) {
		this->wallPoints[numWalls][2*i] = points[2*i];
		this->wallPoints[numWalls][2*i+1] = points[2*i+1];
	}

	this->numWallPoints[this->numWalls] = numPoints;

	this->numWalls++;


}

void ProbeApp::createWalls() {

	// create the static geometry
	Ogre::StaticGeometry* s = mSceneMgr->createStaticGeometry("corridors");

	// region dimensions (FIXME: what does this mean?)
	s->setRegionDimensions(Ogre::Vector3(160.0, 100.0, 160.0));

	// Set the region origin so the center is at 0 world
	s->setOrigin(Ogre::Vector3(0,0,0));

	for ( int wallNum = 0 ; wallNum < this->numWalls ; wallNum++) {




		// create the names we will need for each wall
		std::stringstream ss1;
		ss1 << wallNum;
		std::string objIntName("WallObject_Interior" + ss1.str());
		std::string objExtName("WallObject_Exterior" + ss1.str());
		std::string objTopName("WallObject_Top" + ss1.str());

		std::string entIntName("WallEntity_Interior" + ss1.str());
		std::string entExtName("WallEntity_Exterior" + ss1.str());
		std::string entTopName("WallEntity_Top" + ss1.str());

		std::string meshIntName("WallMesh_Interior" + ss1.str());
		std::string meshExtName("WallMesh_Exterior" + ss1.str());
		std::string meshTopName("WallMesh_Top" + ss1.str());
		
		// start with a ManualObject for the mesh
		Ogre::ManualObject* mo1_int = new Ogre::ManualObject(objIntName);
		Ogre::ManualObject* mo1_ext = new Ogre::ManualObject(objExtName);
		Ogre::ManualObject* mo1_top = new Ogre::ManualObject(objTopName);

		// create the wall1 mesh
		mo1_int->begin("WallMaterial_Int", Ogre::RenderOperation::OT_TRIANGLE_LIST);
		mo1_ext->begin("WallMaterial_Ext", Ogre::RenderOperation::OT_TRIANGLE_LIST);
		mo1_top->begin("WallMaterial_Top", Ogre::RenderOperation::OT_TRIANGLE_LIST);

		double *currPoints = this->wallPoints[wallNum];
		int numPoints = this->numWallPoints[wallNum];

		// create mesh for top view of wall
		for (int i = 0 ; i < numPoints-1 ; i++ ) {

			Ogre::Vector2 vec(currPoints[i*2]-currPoints[(i+1)*2], currPoints[i*2+1]-currPoints[(i+1)*2+1]);

			vec.normalise();

			Ogre::Vector2 orthoVec(vec.x*cos(M_PI/2) - vec.y*sin(M_PI/2), vec.x*sin(M_PI/2) + vec.y*cos(M_PI/2));

			mo1_top->position(Ogre::Vector3(currPoints[i*2], 0.2, currPoints[i*2+1]));
			mo1_top->position(Ogre::Vector3(currPoints[i*2] + 0.05*orthoVec[0], 0.2, currPoints[i*2+1] + 0.05*orthoVec[1]));
			mo1_top->position(Ogre::Vector3(currPoints[(i+1)*2], 0.2, currPoints[(i+1)*2+1]));
			mo1_top->position(Ogre::Vector3(currPoints[(i+1)*2] + 0.05*orthoVec[0], 0.2, currPoints[(i+1)*2+1] + 0.05*orthoVec[1]));

		}

		// create triangles out of vertices
		for ( int i = 0 ; i < numPoints-1 ; i++ ) {
			mo1_top->triangle(4*i+2, 4*i+1, 4*i);
			mo1_top->triangle(4*i+1, 4*i+2, 4*i+3);
			if (i < numPoints-2)
				mo1_top->triangle(4*i+2, 4*i+5, 4*i+3);
		}
	
		// end the number of vertices
		mo1_top->end();

		mo1_top->convertToMesh(meshTopName);

		Ogre::Entity* entity1_top = mSceneMgr->createEntity(entTopName, meshTopName);
		entity1_top->setCastShadows(false);

		s->addEntity(entity1_top, Ogre::Vector3(0.0,0.0,0.0));



		// build the triangle meshes of the inner and outer wall view
		for ( int i = 0 ; i < numPoints ; i++ ) {
			double pX = currPoints[i*2];
			double pZ = currPoints[i*2+1];

			mo1_int->position(pX,0.0,pZ);
			mo1_ext->position(pX,0.0,pZ);
			mo1_int->position(pX,0.2,pZ);
			mo1_ext->position(pX,0.2,pZ);
		}

		for ( int i = 0 ; i < numPoints-1 ; i++ ) {
			mo1_int->triangle(2*i+2, 2*i+1, 2*i);
			mo1_int->triangle(2*i+1, 2*i+2, 2*i+3);

			mo1_ext->triangle(2*i, 2*i+1, 2*i+2);
			mo1_ext->triangle(2*i+3, 2*i+2, 2*i+1);
		}

		// end the number of vertices
		mo1_int->end();
		mo1_ext->end();

		mo1_int->convertToMesh(meshIntName);
		mo1_ext->convertToMesh(meshExtName);

		Ogre::Entity* entity1_int = mSceneMgr->createEntity(entIntName, meshIntName);
		Ogre::Entity* entity1_ext = mSceneMgr->createEntity(entExtName, meshExtName);

		entity1_int->setCastShadows(false);
		entity1_ext->setCastShadows(false);

		s->addEntity(entity1_int, Ogre::Vector3(0.0,0.0,0.0));
		s->addEntity(entity1_ext, Ogre::Vector3(0.0,0.0,0.0));

	}

	// now build the StaticGeometry so it will be rendered
	s->build();

}

void ProbeApp::buildProbe() {

	//self.segMaterials = []

	for (int i = 0 ; i < this->numSegs ; i++ ) {

		Ogre::Vector3 worldPose(i*this->segLength, 0.3, 0.0);

		// Create the visual representation (the Ogre entity and scene node)
		std::stringstream ss;
		ss << i;
		std::string name("segment_" + ss.str());

		Ogre::Entity* entity = mSceneMgr->createEntity(name, "Cube.mesh");
		Ogre::SceneNode* node = mSceneMgr->getRootSceneNode()->createChildSceneNode(name);
		node->attachObject(entity);
		entity->setCastShadows(false);


		// create separate materials for each segment so we can change
		// them in real-time for visualization purposes 
		std::stringstream ss2;
		ss2 << i;
		std::string name2("custom_" + ss2.str());
		Ogre::MaterialPtr newMaterial = entity->getSubEntity(0)->getMaterial()->clone(name2);

		newMaterial->setAmbient(0.0,0.0,1.0);

		entity->getSubEntity(0)->setMaterialName(name2);

		// Pick a size
		Ogre::Vector3 size(this->segLength, this->segHeight, this->segWidth);
		node->setScale(size.x,size.y,size.z);
		
		//Ogre::Vector3* currPos = self->getWorldPose(i);
		//node->setPosition(currPos[0], currPos[1], currPos[2]);
		node->setPosition(worldPose);
		//node->setOrientation(self.getWorldOrientation(i));
		node->setOrientation(Ogre::Quaternion::IDENTITY);
		
		// TODO:  keep pointers to all of these objects 
		// entities, scene nodes, materials
		this->probeNodes[i] = node;

	}

}

void ProbeApp::updatePose(double *positions, double *quaternions) {

	for ( int i = 0 ; i < numSegs ; i++ ) {
		Ogre::Vector3 pos(positions[3*i], positions[3*i+1], positions[3*i+2]);
		Ogre::Quaternion quat(quaternions[4*i], quaternions[4*i+1], quaternions[4*i+2], quaternions[4*i+3]);

		probeNodes[i]->setPosition(pos);
		probeNodes[i]->setOrientation(quat);

	}

}

void ProbeApp::updateCamera(double *pos, double *quat) {

	double xAvg = pos[0];
	double yAvg = pos[1];
	
	Ogre::Vector3 prevPose = mCamera->getPosition();
	double xPrev = prevPose[0] + 1;
	double yPrev = prevPose[2];
	double zPrev = prevPose[1];

	Ogre::Vector3 newPose;
	
	newPose[0] = xPrev*0.99 + 0.01*xAvg;
	newPose[1] = yPrev*0.99 + 0.01*yAvg;
	newPose[2] = 12;


	mCamera->setPosition(Ogre::Vector3(newPose[0]-1,newPose[2],newPose[1]));

	//mCamera->setPosition(Ogre::Vector3(pos[0],pos[1],pos[2]));

	Ogre::Quaternion newQuat(Ogre::Radian(-M_PI/2.0), Ogre::Vector3::UNIT_X);
	mCamera->setOrientation(newQuat);

	//std::cout << "camera quaternion: " << newQuat.x << " " << newQuat.y << " " << newQuat.z << " " << newQuat.w << std::endl;

	//Ogre::Quaternion oriQuat(quat[0],quat[1],quat[2],quat[3]);
	//mCamera->setOrientation(oriQuat);
}



void ProbeApp::render() {
	mRoot->renderOneFrame();
}

void ProbeApp::shutdown() {
	destroyScene();
}

void ProbeApp::saveView(char *name) {
	std::string strName(name);
	mWindow->writeContentsToFile(strName);
}




#if OGRE_PLATFORM == OGRE_PLATFORM_WIN32
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if OGRE_PLATFORM == OGRE_PLATFORM_WIN32
    INT WINAPI WinMain( HINSTANCE hInst, HINSTANCE, LPSTR strCmdLine, INT )
#else
    int main(int argc, char *argv[])
#endif
    {
        // Create application object
        ProbeApp app(40, 0.15, 0.1, 0.15);

        try {
            //app.go();
			//while (1) {
				app.render();
			//}
        } catch( Ogre::Exception& e ) {
#if OGRE_PLATFORM == OGRE_PLATFORM_WIN32
            MessageBox( NULL, e.getFullDescription().c_str(), "An exception has occured!", MB_OK | MB_ICONERROR | MB_TASKMODAL);
#else
            std::cerr << "An exception has occured: " <<
                e.getFullDescription().c_str() << std::endl;
#endif
        }
	
		app.shutdown();

        return 0;
    }

#ifdef __cplusplus
}
#endif
