
# T - junction, from bottom, test 7
wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-5.2], [-3.6,-5.2]]
wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0, 5.2], [-3.6, 5.2], [-3.6, -5.2]]
wall2.reverse()
wall5 = [wall2[-1],wall1[0]]
walls = [wall1, wall2, wall5]

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

