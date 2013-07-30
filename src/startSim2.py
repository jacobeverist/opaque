
import random
random.seed(0)

import os
import sys
from math import *
from copy import copy
sys.path.insert(1,"modules/")
sys.path.insert(1,"modules/nelmin")
sys.path.insert(1,"modules/bulletprobe")
sys.path.insert(1,"modules/medialaxis")
sys.path.insert(1,"modules/alphashape")
sys.path.insert(1,"modules/ogre")

#os.environ['PATH'] = 'c:\\opencv-2.4.6.0\\build\\x86\\vc9\\bin' + ';' + os.environ['PATH']
#os.environ['PATH'] = os.getcwd() + '\\binaries\\boost_1_53_0\\lib' + ';' + os.environ['PATH']

from bulletsnakeprobe import BulletSnakeProbe
from ogreprobe import ProbeApp


"""
numSegs = 40
segLength = 0.15
segHeight = 0.1
segWidth = 0.15
maxTorque = 1000.0
friction = 1.0
quat = [0.0, 1.0, 0.0, -1.62920684943e-07]
pos = [1.65, 0.04, 0.0]

bullet_snake = BulletSnakeProbe(quat, pos, numSegs, segLength, segHeight, segWidth, friction)


# Y - junction, test 1 				
WLEN2 = 7.0
wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN2*cos(pi/3), -0.2 - WLEN2*sin(pi/3)]]
wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
wall5 = [wall2[2],wall1[0]]
w1 = wall1[2]
w2 = wall2[0]

wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
wall4 = [w1, [w1[0] + 0.4*cos(pi/3-pi/2), w1[1] - 0.4*sin(pi/3-pi/2)]]

walls = [wall1, wall2, wall3, wall4, wall5]

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p
		
for wall in walls:
	bullet_snake.addWall(wall)

visProbe = ProbeApp(numSegs, segLength, segHeight, segWidth)

for wall in walls:
	visProbe.addWall(wall)

visProbe.createWalls()

while True:
	bullet_snake.frameStarted()

	positions = []
	quaternions = []
	for segI in range(numSegs):
		currPos = bullet_snake.getGlobalPosition(segI)
		quat_vars = bullet_snake.getGlobalOrientation(segI)

		positions.append(currPos)
		quaternions.append([quat_vars[3], quat_vars[0], quat_vars[1], quat_vars[2]])
	
	visProbe.updatePose(positions, quaternions)

	visProbe.render()

visProbe.shutdown()
"""

from opaque.robot.QuickProbe import QuickProbe
from opaque.robot.BulletProbe import BulletProbe

from opaque.TestNavigation import TestNavigation
from opaque.TestProcess import TestProcess
from opaque.TestUnit import TestUnit
from opaque.ControlError import *
from opaque.CmdDrawThings import CmdDrawThings
from opaque.DrawThings import DrawThings

import traceback 
import cProfile

sys.setrecursionlimit(10000)


" sub process cleanup code "
import opaque.maps.gen_icp as gen_icp
import atexit
def cleanup():
	if gen_icp.overlapPool != None:
		for p in gen_icp.overlapPool:
			p.terminate()
atexit.register(cleanup)


def createTest():
	#probe = QuickProbe(40,0.15,0.05,30.0*2.5)

	quat = [0.0, 1.0, 0.0, -1.62920684943e-07]
	pos = [1.65, 0.04, 0.0]
	probe = BulletProbe(quat,pos,40,0.15,0.1,0.15,1000.0,1.0)
	#probe = BulletProbe(quat,pos,40,0.15,0.1,0.15,1000.0,10.0)

	
	drawThings = DrawThings(probe)
	currControl = TestNavigation(probe, drawThings)
	
	probe.addControl(currControl)

	# Y - junction, test 1 				
	WLEN2 = 7.0
	wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN2*cos(pi/3), -0.2 - WLEN2*sin(pi/3)]]
	wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
	wall5 = [wall2[2],wall1[0]]
	w1 = wall1[2]
	w2 = wall2[0]

	wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
	wall4 = [w1, [w1[0] + 0.4*cos(pi/3-pi/2), w1[1] - 0.4*sin(pi/3-pi/2)]]

	walls = [wall1, wall2, wall3, wall4, wall5]

	for wall in walls:
		for i in range(len(wall)):
			p = copy(wall[i])
			p[0] += 6.0
			wall[i] = p
	
	probe.setWalls(walls)

	drawThings.setWalls(walls)

	return probe, drawThings

def runTest(probe, drawThings):

	while True:
		probe.frameStarted()
		
		adjustCamera( probe, drawThings)

		drawThings.updatePosture()
		drawThings.render()


def adjustCamera(probe, drawThings):

	pose = probe.getActualJointPose(9)
	#pose = self.probe.getActualJointPose(19)
	xAvg = pose[0]
	yAvg = pose[1]
	
	#prevPose = self.camera.getPosition()
	#xPrev = prevPose[0] + 1
	#yPrev = prevPose[2]
	#zPrev = prevPose[1]
	
	#newPose = [xPrev*0.99 + 0.01*xAvg, yPrev*0.99 + 0.01*yAvg, 12]
	newPose = [xAvg, yAvg, 12]


	#self.camera.setPosition(newPose[0]-1,newPose[2],newPose[1])

	#oriQuat = ogre.Quaternion(ogre.Radian(-pi/2.0), ogre.Vector3().UNIT_X)
	###self.camera.setOrientation(oriQuat)
	
	angQuat = [-0.707107, 0.0, 1.0, 0.707107]

	#drawThings.updateCamera([newPose[0]-1,newPose[2],newPose[1]], angQuat)
	drawThings.updateCamera(newPose, angQuat)
		
def startRun():
	try:
		" create the robot "
		" create the environment "
		" add control program "		
		probe, drawThings = createTest()
		
		runTest(probe, drawThings)
		
		" while loop "		
		#cProfile.run('runTest(probe)', 'test_prof')
	
	except ControlError as inst:
		print inst.value

	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]
	

if __name__ == '__main__':

	try:
		" create the robot "
		" create the environment "
		" add control program "		
		probe, drawThings = createTest()

		runTest(probe, drawThings)
		
		" while loop "		
		#cProfile.run('runTest(probe)', 'test_prof')
	
	except ControlError as inst:
		print inst.value

	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]

