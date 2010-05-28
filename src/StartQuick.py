from opaque.robot.QuickProbe import QuickProbe
#from opaque.TestModular import TestModular
from opaque.TestHoldPosition import TestHoldPosition
from opaque.TestHoldTransition import TestHoldTransition
from opaque.TestHoldTransition import TestHoldTransition
from opaque.TestHoldSlideTransition import TestHoldSlideTransition
from opaque.TestAnchorTransition import TestAnchorTransition
from opaque.TestFrontExtend import TestFrontExtend
from opaque.TestPokeWalls import TestPokeWalls
from opaque.TestAdaptiveStep import TestAdaptiveStep
from opaque.TestPathStep import TestPathStep
from opaque.TestMapGraph import TestMapGraph
from opaque.TestModular import TestModular
from opaque.TestTransform import TestTransform
from opaque.ControlError import *
from DrawThings import DrawThings

import traceback 
import cProfile
import sys

def createTest():
	probe = QuickProbe(40,0.15,0.05,30.0*2.5)
	
	drawThings = DrawThings(probe.robotParam)
	#currControl = TestHoldPosition(probe, drawThings)
	#currControl = TestHoldTransition(probe, drawThings)
	#currControl = TestHoldSlideTransition(probe, drawThings)
	#currControl = TestAnchorTransition(probe, drawThings)
	#currControl = TestFrontExtend(probe, drawThings)
	#currControl = TestPokeWalls(probe, drawThings)
	#currControl = TestAdaptiveStep(probe, drawThings)
	#currControl = TestPathStep(probe, drawThings)
	#currControl = TestMapGraph(probe, drawThings)
	#currControl = TestModular(probe, drawThings)
	currControl = TestTransform(probe, drawThings)
	
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
		
		
		" while loop "		
		cProfile.run('runTest(probe)', 'test_prof')
	
	except ControlError as inst:
		print inst.value

	except:
		traceback.print_exc()
		print "Exception:", sys.exc_info()[0]

