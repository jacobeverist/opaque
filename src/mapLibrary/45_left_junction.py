# 45 degree junction, left, test 6

WLEN2 = 5.0
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

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

