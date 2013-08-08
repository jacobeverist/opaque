
# Y - junction, test 1 				
WLEN2 = 7.0
wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN2*cos(pi/3), -0.2 - WLEN2*sin(pi/3)]]
wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
wall5 = [wall2[2],wall1[0]]
w1 = wall1[-1]
w2 = wall2[0]

wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
wall4 = [w1, [w1[0] + 0.4*cos(pi/3-pi/2), w1[1] - 0.4*sin(pi/3-pi/2)]]

walls = [wall1, wall2, wall3, wall4, wall5]

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

