
import random
random.seed(0)

import os
import sys
sys.path.insert(1,"modules/")
sys.path.insert(1,"modules/nelmin")
sys.path.insert(1,"modules/bulletprobe")
sys.path.insert(1,"modules/medialaxis")
sys.path.insert(1,"modules/alphashape")
#sys.path.insert(1,"binaries/opencv")


os.environ['PATH'] = 'c:\\opencv-2.4.6.0\\build\\x86\\vc9\\bin' + ';' + os.environ['PATH']
os.environ['PATH'] = os.getcwd() + '\\binaries\\boost_1_53_0\\lib' + ';' + os.environ['PATH']


from opaque.robot.QuickProbe import QuickProbe
from opaque.TestProcess import TestProcess
from opaque.TestUnit import TestUnit
from opaque.ControlError import *
from opaque.CmdDrawThings import CmdDrawThings

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
	probe = QuickProbe(40,0.15,0.05,30.0*2.5)
	
	drawThings = CmdDrawThings(probe.robotParam)
	currControl = TestProcess(probe, drawThings)
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

