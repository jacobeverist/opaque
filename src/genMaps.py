
#from opaque.common import *
#import opaque.environment as envs

from opaque.environment.MapGenerate import ForkEnv, PipeJunctions, Corridor, StraightJunction, CrossJunction
from math import pi
import pylab

if __name__ == '__main__':
	
	'pipeLength', 'pipeSegMin', 'pipeSegMax', 'pipeAngMax', 'pipeWidthMin', 'pipeWidthMax', 'pipeWallSegLen', 'pipeEntranceLen', 'wallPerturbRange'
	
	maps = []
	for i in range(17):
		
		
		# pipe_len, seg_len_min, seg_len_max, seg_ang_max, pipe_min, pipe_max, feature_len, entr_len, perturb_range)
		try:
			#newmap = Corridor(20.0, 0.1, 0.5, pi/4, 0.4, 0.7, 0.2, 10.0, 0.2)
			newmap = Corridor(20.0, 0.02, 0.2, pi/4, 0.4, 0.41, 0.2, 10.0, 0.1)
		except:
			pass
		else:
			#newmap = ForkEnv(10.0, 0.1, 0.5, pi/16, 0.3, 0.4, 0.2, 5.0, 0.03)
			#newmap = ForkEnv(10.0, 0.1, 0.5, pi/512, 0.3, 0.4, 0.2, 5.0, 0.00)
			#newmap = Corridor(10.0, 0.1, 0.5, pi/8, 0.4, 0.7, 0.2, 5.0, 0.2)
			#newmap = ForkEnv(10.0, 0.1, 0.5, pi/16, 0.3, 0.4, 0.2, 5.0, 0.03)
			#newmap = ForkEnv(10.0, 0.1, 0.5, pi/512, 0.3, 0.4, 0.2, 5.0, 0.00)
			#newmap = CrossJunction(3.0, 5.0, 0.4, pi/4)
			#newmap = StraightJunction(3.0, 5.0, 0.4, angle)
			#newmap = PipeJunctions(5.0, 5.0, angle)
	
			
			#angle = i * pi/16 - pi/2
			#newmap = StraightJunction(3.0, 5.0, 0.4, angle)
			#newmap = PipeJunctions(5.0, 5.0, angle)
			
			# (pipe_len, entr_len, pipe_width, turn_angle)
			
			newmap.writeToFile(".")
			maps.append(newmap)
			
			pylab.clf()
			newmap.drawCorridor()
			#newmap.draw()
	
			pylab.xlim(-15,15)
			pylab.ylim(-15,15)
			pylab.savefig("genMap_%04u.png" % (i+1))
	
	
