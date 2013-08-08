

# T - junctions, from side, test 2
wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4,-4.2], [-3.6,-4.2], [-3.6,-0.2], [2.0,-0.2], [2.0,0.2]]
wall2 = [[-14.0, 0.2], [-1.0, 0.2], [-1.0,3.6], [-0.6,3.6], [-0.6,0.2], [2.0,0.2]]
wall2.reverse()

# close the entrance
wall5 = [wall2[-1],wall1[0]]
walls = [wall1, wall2, wall5]

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

