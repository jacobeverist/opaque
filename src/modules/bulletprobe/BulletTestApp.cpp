// BulletTestApp.cpp : Defines the entry point for the console application.
//

#ifdef _WIN32
#include "stdafx.h"
#endif

#include <iostream>
#include <btBulletDynamicsCommon.h>
#include "BulletSnake.h"


#ifdef _WIN32
int _tmain(int argc, _TCHAR* argv[])
#else
int main(int argc, char argv[])
#endif

{
	///std::cout << "Hello World!" << std::endl;

	/*
    btDefaultMotionState* fallMotionState =
            new btDefaultMotionState(btTransform(btQuaternion(0,0,0,1),btVector3(0,50,0)));
    btScalar mass = 1;
    btVector3 fallInertia(0,0,0);
    fallShape->calculateLocalInertia(mass,fallInertia);
    btRigidBody::btRigidBodyConstructionInfo fallRigidBodyCI(mass,fallMotionState,fallShape,fallInertia);
    btRigidBody* fallRigidBody = new btRigidBody(fallRigidBodyCI);
    dynamicsWorld->addRigidBody(fallRigidBody);
	*/

	double quat_param[4] = {0.0, 1.0, 0.0, 0.0};
	double pos[3] = {1.65, 0.04, 0.0};
	BulletSnake *ns = new BulletSnake(quat_param, pos, 40, 0.15, 0.2, 0.15, 0.7);

	// create the walls
	double wall1[] = {-8.0, -0.20000000000000001, 2.0, -0.20000000000000001, 3.5, -2.7980763912200928};
	double wall2[] = {4.5, 4.5301270484924316, 2.0, 0.20000000000000001, -8.0, 0.20000000000000001};
	double wall3[] = {3.8464102745056152, -2.598076343536377, 2.3464102745056152, 0.0, 4.846409797668457, 4.3301272392272949, 4.5, 4.5301270484924316};
	double wall4[] = {3.5, -2.7980763912200928, 5.1999998092651367, -5.7425627708435059, 5.546410083770752, -5.5425629615783691, 4.0464105606079102, -2.9444866180419922, 6.6444864273071289, -1.4444864988327026, 6.4444866180419922, -1.098076343536377, 3.8464102745056152, -2.598076343536377};

	ns->addWall(3, wall1);
	ns->addWall(3, wall2);
	ns->addWall(4, wall3);
	ns->addWall(7, wall4);
	ns->createWalls();




    for (int i=0 ; i<300 ; i++) {
		ns->Step();
    }


	delete ns;

	/*
    for (int i=0 ; i<300 ; i++) {
            dynamicsWorld->stepSimulation(1/60.f,10);

            btTransform trans;
            fallRigidBody->getMotionState()->getWorldTransform(trans);

            std::cout << "sphere height: " << trans.getOrigin().getY() << std::endl;
    }
	*/

	/*
    dynamicsWorld->removeRigidBody(fallRigidBody);
    delete fallRigidBody->getMotionState();
    delete fallRigidBody;

    dynamicsWorld->removeRigidBody(groundRigidBody);
    delete groundRigidBody->getMotionState();
    delete groundRigidBody;


    delete fallShape;

    delete groundShape;
	*/

	return 0;
}
