

# T - junctions, from side, test 2
wall1 = [[-11.0, -0.2], [-4.0, -0.2], [-4,-5.2], [-3.6,-5.2], [-3.6,-0.2], [6.8,-0.2], [6.8,0.2]]
wall2 = [[-11.0, 0.2], [1.4, 0.2], [1.4,5.2], [1.8,5.2], [1.8,0.2], [6.8,0.2]]
wall2.reverse()

# close the entrance
wall5 = [wall2[-1],wall1[0]]
walls = [wall1, wall2, wall5]

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

