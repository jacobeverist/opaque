
#from opaque.common import *
#import opaque.environment as envs

from opaque.environment.MapGenerate import ForkEnv, PipeJunctions, Corridor, StraightJunction
from math import pi
import pylab

if __name__ == '__main__':
	
	'pipeLength', 'pipeSegMin', 'pipeSegMax', 'pipeAngMax', 'pipeWidthMin', 'pipeWidthMax', 'pipeWallSegLen', 'pipeEntranceLen', 'wallPerturbRange'
	
	maps = []
	for i in range(16):
		
		angle = i * 0.1 + 0.1

		#newmap = envs.Corridor(10.0, 0.1, 0.5, pi/8, 0.25, 0.7, 0.2, 5.0, 0.2)
		#newmap = ForkEnv(10.0, 0.1, 0.5, pi/16, 0.3, 0.4, 0.2, 5.0, 0.03)
		#newmap = ForkEnv(10.0, 0.1, 0.5, pi/512, 0.3, 0.4, 0.2, 5.0, 0.00)
		newmap = StraightJunction(3.0, 5.0, 0.4, angle)

		# (pipe_len, entr_len, pipe_width, turn_angle)
		
		newmap.writeToFile(".")
		maps.append(newmap)
		
		#pylab.clf()
		#newmap.drawCorridor()
		#newmap.draw()

		pylab.xlim(-10,10)
		pylab.ylim(-10,10)
		pylab.savefig("genMap_%04u.png" % i)
	
	#newmap.draw()
	
	