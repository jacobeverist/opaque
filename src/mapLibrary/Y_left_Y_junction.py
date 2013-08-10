

# triple-Y - multi-junction, test 1 				
WLEN = 5.0
WWID = 0.4
wall1 = [[-11.0, -WWID/2.0], [-4.0, -WWID/2.0], [-4.0 + WLEN*cos(-pi/3), -WWID/2.0 + WLEN*sin(-pi/3)]]
wall2 = [[-11.0, WWID/2.0], [-4.0, WWID/2.0], [-4.0 + WLEN*cos(pi/3), WWID/2.0 + WLEN*sin(pi/3)]]
#wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
w1 = wall1[2]
w2 = wall2[2]

wall3 = [[w1[0] + WWID*cos(pi/6), w1[1] + WWID*sin(pi/6)], [WWID*cos(pi/6) - 4, 0.0], [w2[0] + WWID*cos(pi/6), w2[1] - WWID*sin(pi/6)]]

wall4 = [w1, [w1[0] + WLEN*cos(-2*pi/3), w1[1] + WLEN*sin(-2*pi/3)]]			
wall5 = [w2, [w2[0] + WLEN*cos(2*pi/3), w2[1] + WLEN*sin(2*pi/3)]]

wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi/2), wall4[-1][1] + WWID*sin(-2*pi/3 + pi/2)])
wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi)])
wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3)])
wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi - pi/3 + pi/2.), wall4[-1][1] + WWID*sin(-2*pi/3 + pi - pi/3 + pi/2.)])
wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3 + pi)])

wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi/2), wall5[-1][1] + WWID*sin(2*pi/3 - pi/2)])
wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi)])
wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3)])
wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi + pi/3 - pi/2.), wall5[-1][1] + WWID*sin(2*pi/3 - pi + pi/3 - pi/2.)])
wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3 - pi)])

wall7 = [wall5[0], wall5[-1]]

#wall6 = [[-11.0, WWID/2.0],[-11.0, -WWID/2.0]]
wall6 = [wall2[0], wall1[0]]

wall2.reverse()
wall5.reverse()
wall7.reverse()

#walls = [wall1, wall2, wall3, wall4, wall5, wall6]
walls = [wall1, wall2, wall3, wall4, wall7, wall6]


for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] += 6.0
		wall[i] = p

