
from opaque.common import *
import opaque.environment as envs


if __name__ == '__main__':
	
	'pipeLength', 'pipeSegMin', 'pipeSegMax', 'pipeAngMax', 'pipeWidthMin', 'pipeWidthMax', 'pipeWallSegLen', 'pipeEntranceLen', 'wallPerturbRange'
	
	maps = []
	for i in range(20):
		newmap = envs.Corridor(10.0, 0.1, 0.5, pi/8, 0.25, 0.7, 0.2, 5.0, 0.2)
		newmap.writeToFile(".")
		#newmap = envs.ForkEnv(10.0, 0.1, 0.5, pi/16, 0.3, 0.4, 0.2, 5.0, 0.03)
		#newmap = envs.PipeJunctions(10.0, 5.0, pi/4)
		maps.append(newmap)
	
	#newmap.draw()
	
	