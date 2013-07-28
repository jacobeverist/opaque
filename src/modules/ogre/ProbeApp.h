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
    ProbeApp(void);
    virtual ~ProbeApp(void);

	void addWall(int numPoints, float *points);
	void createWalls();
	void setupMyWorld(Ogre::Quaternion R, Ogre::Vector3 pos);


protected:
    virtual void createScene(void);


private:
	int numWalls;
	float **wallPoints;
	int *numWallPoints;

};

#endif // #ifndef __ProbeApp_h_
