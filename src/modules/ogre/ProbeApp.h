/*
-----------------------------------------------------------------------------
Filename:    ProbeApp.h
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
#ifndef __ProbeApp_h_
#define __ProbeApp_h_

#include "BaseApplication.h"

class ProbeApp : public BaseApplication
{
public:
    ProbeApp(int numSegs, double segLength, double segHeight, double segWidth);
    virtual ~ProbeApp(void);

	void addWall(int numPoints, double *points);
	void createWalls();
	//void setupMyWorld(Ogre::Quaternion R, Ogre::Vector3 pos);


	void buildProbe();
	void updatePose(double *postions, double *quaternions);

	void updateCamera(double *pos, double *quat);


	void render();
	void shutdown();

protected:
    virtual void createScene(void);


private:
	int numWalls;
	double **wallPoints;
	int *numWallPoints;

	int numSegs;
	double segLength;
	double segHeight;
	double segWidth;

	Ogre::SceneNode **probeNodes;
	Ogre::Material **probeMaterials;
	Ogre::Entity **probeEntities;

};

#endif // #ifndef __ProbeApp_h_
