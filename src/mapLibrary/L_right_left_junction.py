# L - junction, right-left test

wall1 = [[-14.0, -0.2], [-3.6,-0.2], [-3.6,0.2], [-3.6,5.8], [2.0,5.8]]
wall2 = [[-14.0, 0.2], [-4.0, 0.2], [-4.0,6.2], [2.0,6.2] ]
wall2.reverse()
wall3 = [wall2[0],wall1[-1]]
wall4 = [wall2[-1],wall1[0]]
walls = [wall1, wall2, wall3, wall4]



for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

