# 45 degree junction, left, test 6


PLEN1 = 10.0


mag1 = -0.2/-sin(pi/3+pi/2-pi/6)
x1 = mag1*cos(pi/3+pi/2-pi/6)

mag2 = 0.2/-sin(pi/3-pi/2-pi/6)
x2 = mag2*cos(pi/3-pi/2-pi/6)

wall1 = [[-14.0, -0.2], [-14.0+PLEN1 + x1, -0.2]]
wall2 = [[-14.0,  0.2], [-14.0+PLEN1 + x2,  0.2]]


WLEN0 = 6.0
#WLEN1 = 5.0 - 0.2*cos(pi/3)
#WLEN2 = 5.0 + 0.2*cos(pi/3)
WLEN1 = WLEN0 - 0.2*cos(pi/3)
WLEN2 = WLEN0 + 0.2*cos(pi/3)

#WLEN1 = WLEN2 = 5.0

wall3 = [wall1[-1], [wall1[-1][0] + WLEN0*cos(pi/3), wall1[-1][1] - WLEN0*sin(pi/3)]]
wall4 = [wall2[-1], [wall2[-1][0] + WLEN0*cos(pi/3), wall2[-1][1] - WLEN0*sin(pi/3)]]

wall5 = [wall3[-1], [wall3[-1][0] + WLEN2*cos(0), wall3[-1][1] - WLEN2*sin(0)]]
wall6 = [wall4[-1], [wall4[-1][0] + WLEN1*cos(0), wall4[-1][1] - WLEN1*sin(0)]]

wall2.reverse()
wall4.reverse()
wall6.reverse()

wall7 = [wall2[-1],wall1[0]]
wall8 = [wall5[-1],wall6[0]]

walls = [wall1, wall2,wall3, wall4, wall5, wall6, wall7, wall8]

"""
WLEN3 = 5.5

wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN2*cos(pi/3), -0.2 - WLEN2*sin(pi/3)]]
wall2 = [[-3.0, 0.2], [-14.0, 0.2]]

w1 = wall1[-1]
w2 = wall2[0]			

wall3 = [w1,
		[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)],
		[w1[0] + 0.4*cos(pi/6) - WLEN3*cos(pi/3), w1[1] + 0.4*sin(pi/6) + WLEN3*sin(pi/3)]]

wall5 = [wall2[-1],wall1[0]]
walls = [wall1, wall2, wall3, wall5]
"""

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

