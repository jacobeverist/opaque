# L - junction, left, test 3

wall1 = [[-11.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [-3.6,0.2]]
wall2 = [[-11.0, 0.2], [-3.6,0.2]]
wall2.reverse()
wall5 = [wall2[-1],wall1[0]]
walls = [wall1, wall2, wall5]

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

