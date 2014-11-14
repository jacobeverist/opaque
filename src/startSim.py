
from math import *
from copy import copy
import random
random.seed(0)

import os
import sys

relPath = os.path.dirname(os.path.realpath(__file__))
print relPath

sys.path.insert(1,relPath)
sys.path.insert(1,relPath + "/modules/")
sys.path.insert(1,relPath + "/modules/nelmin")
sys.path.insert(1,relPath + "/modules/bulletprobe")
sys.path.insert(1,relPath + "/modules/medialaxis")
sys.path.insert(1,relPath + "/modules/alphashape")
sys.path.insert(1,relPath + "/modules/ogre")

#os.environ['PATH'] = 'c:\\opencv-2.4.6.0\\build\\x86\\vc9\\bin' + ';' + os.environ['PATH']
#os.environ['PATH'] = os.getcwd() + '\\binaries\\boost_1_53_0\\lib' + ';' + os.environ['PATH']

import argparse


parser = argparse.ArgumentParser(description='DarkMapper Simulator')
parser.add_argument('--mapFile', type=str, help="text file containing environment map")
parser.add_argument('--restoreConfig', type=bool, default=False, help="restore ogre config from file instead of diplaying dialog")
parser.add_argument('--hideWindow', type=bool, default=False, help="hide Ogre window")
# maxNumPoses
parser.add_argument('--maxNumPoses', type=int, default=300, help="max number of poses")
# numPoseParticles
parser.add_argument('--numPoseParticles', type=int, default=40, help="number of pose particles")
# bloomFeature
parser.add_argument('--bloomFeature', type=bool, default=True, help="use bloom spatial features")
# bendFeature 
parser.add_argument('--bendFeature', type=bool, default=True, help="use bend spatial features")
# startOffset
parser.add_argument('--startOffset', type=float, default=0.0, help="probe start displacement offset")
args = parser.parse_args()



from bulletsnakeprobe import BulletSnakeProbe
from ogreprobe import ProbeApp

from opaque.robot.QuickProbe import QuickProbe
from opaque.robot.BulletProbe import BulletProbe

from opaque.TestNavigation import TestNavigation
from opaque.TestProcess import TestProcess
#from opaque.TestUnit import TestUnit
from opaque.ControlError import *
from opaque.CmdDrawThings import CmdDrawThings
from opaque.DrawThings import DrawThings

import traceback 
import cProfile

sys.setrecursionlimit(10000)

" sub process cleanup code "
import opaque.maps.gen_icp as gen_icp
import opaque.maps.MapProcess as MapProcess
import atexit
def cleanup():
	if gen_icp.overlapPool != None:
		for p in gen_icp.overlapPool:
			p.terminate()

	if len(MapProcess.pool_move) > 0:
		for p in MapProcess.pool_move:
			p.terminate()
	if len(MapProcess.pool_localize) > 0:
		for p in MapProcess.pool_localize:
			p.terminate()
	if len(MapProcess.pool_generate) > 0:
		for p in MapProcess.pool_generate:
			p.terminate()
	if len(MapProcess.pool_eval) > 0:
		for p in MapProcess.pool_eval:
			p.terminate()

			
atexit.register(cleanup)


def createTest():
	#probe = QuickProbe(40,0.15,0.05,30.0*2.5)


	numSegs = 40
	segLength = 0.15
	segHeight = 0.1
	segWidth = 0.15
	maxTorque = 1000.0
	friction = 1.0
	quat = [0.0, 1.0, 0.0, -1.62920684943e-07]
	pos = [1.65 + args.startOffset, 0.04, 0.0]
	#pos = [1.35, 0.04, 0.0]
	probe = BulletProbe(quat, pos, numSegs, segLength, segHeight, segWidth, maxTorque, friction)
	#probe = BulletProbe(quat,pos,40,0.15,0.1,0.15,1000.0,1.0)
	#probe = BulletProbe(quat,pos,40,0.15,0.1,0.15,1000.0,10.0)

	
	drawThings = DrawThings(probe, hideWindow = args.hideWindow, restoreConfig = args.restoreConfig)
	#drawThings = CmdDrawThings(probe.robotParam)

	# args.maxNumPoses
	# args.numPoseParticles
	# args.bloomFeature
	# args.bendFeature 
	currControl = TestNavigation(probe, drawThings, args)
	#maxNumPoses = args.maxNumPoses, numPoseParticles = args.numPoseParticles, bloomFeature = args.bloomFeature, bendFeature=args.bendFeature)
	
	probe.addControl(currControl)

	if args.mapFile != None:
		f = open(args.mapFile, 'r')
		saveStr = f.read()
		f.close()
		exec(saveStr)
	else:

		# junction test 
		targetMapFile = "testData/curveEnv/mapFile0004.txt"
		#targetMapFile = "testData/curveEnv/mapFile0010.txt"
		#targetMapFile = "testData/curveEnv/mapFile0002.txt"
		f = open(relPath + "/" + targetMapFile, 'r')
		str_f = f.read()
		walls = eval(str_f)
		wall5 = [walls[0][0], walls[0][-1]]
		walls.append(wall5)
	
	probe.setWalls(walls)

	drawThings.setWalls(walls)

	return probe, drawThings

def runTest(probe, drawThings):

	while True:
		probe.frameStarted()
		
		adjustCamera( probe, drawThings)

		#drawThings.render()


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

