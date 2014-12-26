
import random
random.seed(0)

import os
import sys

relPath = os.path.dirname(os.path.realpath(__file__))
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
# maxNumPoses
parser.add_argument('--maxNumPoses', type=int, default=300, help="max number of poses")
# numPoseParticles
parser.add_argument('--numPoseParticles', type=int, default=40, help="number of pose particles")

# bloomFeature
parser.add_argument('--bloomFeature', dest='bloomFeature', action='store_true')
parser.add_argument('--no-bloomFeature', dest='bloomFeature', action='store_false')
parser.set_defaults(bloomFeature=True)

# bendFeature
parser.add_argument('--bendFeature', dest='bendFeature', action='store_true')
parser.add_argument('--no-bendFeature', dest='bendFeature', action='store_false')
parser.set_defaults(bendFeature=False)


# startOffset
parser.add_argument('--startOffset', type=float, default=0.0, help="probe start displacement offset")
args = parser.parse_args()




from opaque.robot.QuickProbe import QuickProbe
from opaque.TestProcess import TestProcess
#from opaque.TestUnit import TestUnit
from opaque.ControlError import *
from opaque.CmdDrawThings import CmdDrawThings

import traceback 
import cProfile

sys.setrecursionlimit(10000)


" sub process cleanup code "
import opaque.maps.gen_icp as gen_icp
import opaque.maps.MapProcess as MapProcess
import opaque.maps.shoots as shoots
import opaque.maps.ParticleFilter as ParticleFilter
import atexit

def cleanup():

	if len(shoots.pool_branch) > 0:

		shoots.qin_branch.close()
		shoots.qin_branch.join_thread()

		shoots.qout_branch.close()
		qout_branch.cancel_join_thread()

		for p in shoots.pool_branch:
			p.terminate()

	if len(ParticleFilter.pool_posePart) > 0:

		ParticleFilter.qin_posePart.close()
		ParticleFilter.qin_posePart.join_thread()

		ParticleFilter.qout_posePart.close()
		ParticleFilter.qout_posePart.cancel_join_thread()

		for p in ParticleFilter.pool_posePart:
			p.terminate()

	if len(ParticleFilter.pool_dispPosePart2) > 0:

		ParticleFilter.qin_dispPosePart2.close()
		ParticleFilter.qin_dispPosePart2.join_thread()

		ParticleFilter.qout_dispPosePart2.close()
		ParticleFilter.qout_dispPosePart2.cancel_join_thread()

		for p in ParticleFilter.pool_dispPosePart2:
			p.terminate()

	if len(ParticleFilter.pool_localLandmark) > 0:

		ParticleFilter.qin_localLandmark.close()
		ParticleFilter.qin_localLandmark.join_thread()

		ParticleFilter.qout_localLandmark.close()
		ParticleFilter.qout_localLandmark.cancel_join_thread()

		for p in ParticleFilter.pool_localLandmark:
			p.terminate()

			
atexit.register(cleanup)


def createTest():
	probe = QuickProbe(40,0.15,0.05,30.0*2.5)
	
	drawThings = CmdDrawThings(probe.robotParam)
	currControl = TestProcess(probe, drawThings, args)
	#currControl = TestUnit(probe, drawThings)
	
	probe.addControl(currControl)

	return probe

def runTest(probe):

	while True:
		probe.frameStarted()


def startRun():
	try:
		" create the robot "
		" create the environment "
		" add control program "		
		probe = createTest()
		
		runTest(probe)
		
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
		probe = createTest()

		runTest(probe)
		
		" while loop "		
		#cProfile.run('runTest(probe)', 'test_prof')
	
	except ControlError as inst:
		print inst.value

	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]

