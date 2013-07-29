
from math import *
import ogreprobe

probe = ogreprobe.ProbeApp(40, 0.15, 0.1, 0.15)

WLEN2 = 7.0
wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN2*cos(pi/3), -0.2 - WLEN2*sin(pi/3)]]
wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
wall5 = [wall2[2],wall1[0]]
w1 = wall1[2]
w2 = wall2[0]
wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
wall4 = [w1, [w1[0] + 0.4*cos(pi/3-pi/2), w1[1] - 0.4*sin(pi/3-pi/2)]]

probe.addWall(wall1)
probe.addWall(wall2)
probe.addWall(wall5)
probe.addWall(wall3)
probe.addWall(wall4)

probe.createWalls()


while True:
	probe.render()


probe.shutdown()
