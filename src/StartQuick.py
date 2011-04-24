from opaque.robot.QuickProbe import QuickProbe
from opaque.TestProcess import TestProcess
from opaque.ControlError import *
from CmdDrawThings import CmdDrawThings

import traceback 
import cProfile
import sys

def createTest():
	probe = QuickProbe(40,0.15,0.05,30.0*2.5)
	
	drawThings = CmdDrawThings(probe.robotParam)
	currControl = TestProcess(probe, drawThings)
	
	probe.addControl(currControl)

	return probe

def runTest(probe):

	while True:
		probe.frameStarted()
	
if __name__ == '__main__':

	try:
		" create the robot "
		" create the environment "
		" add control program "		
		probe = createTest()

		
		#runTest(probe)
		
		
		" while loop "		
		cProfile.run('runTest(probe)', 'test_prof')
	
	except ControlError as inst:
		print inst.value

	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]

