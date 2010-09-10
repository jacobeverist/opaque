
#from opaque.common import *
#import opaque.environment as envs

from opaque.environment.MapGenerate import ForkEnv, PipeJunctions, Corridor, StraightJunction, CrossJunction
from math import pi
import pylab

if __name__ == '__main__':
	
	'pipeLength', 'pipeSegMin', 'pipeSegMax', 'pipeAngMax', 'pipeWidthMin', 'pipeWidthMax', 'pipeWallSegLen', 'pipeEntranceLen', 'wallPerturbRange'
	
	maps = []
	for i in range(17):
	#for i in range(1):
		
		angle = i * pi/16 - pi/2
		#angle = i * 0.1 + 0.1

		#					pipe_len, seg_len_min, seg_len_max, seg_ang_max, pipe_min, pipe_max, feature_len, entr_len, perturb_range)
		#newmap = Corridor(10.0, 0.1, 0.5, pi/8, 0.4, 0.7, 0.2, 5.0, 0.2)
		#newmap = ForkEnv(10.0, 0.1, 0.5, pi/16, 0.3, 0.4, 0.2, 5.0, 0.03)
		#newmap = ForkEnv(10.0, 0.1, 0.5, pi/512, 0.3, 0.4, 0.2, 5.0, 0.00)

		
		#newmap = CrossJunction(3.0, 5.0, 0.4, pi/4)
		newmap = StraightJunction(3.0, 5.0, 0.4, angle)
		#newmap = PipeJunctions(5.0, 5.0, angle)
		
		# (pipe_len, entr_len, pipe_width, turn_angle)
		
		newmap.writeToFile(".")
		maps.append(newmap)
		
		pylab.clf()
		#newmap.drawCorridor()
		newmap.draw()

		pylab.xlim(-10,10)
		pylab.ylim(-10,10)
		pylab.savefig("genMap_%04u.png" % i)
	
	#newmap.draw()
	
	