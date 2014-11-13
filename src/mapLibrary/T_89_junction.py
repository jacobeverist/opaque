from opaque.environment.MapGenerate import MapMaker


foo = MapMaker(8.0, 0.4)

segA, segB = foo.addSymJunction(0, 49*pi/100.0)

foo.addStraight(segA,6.0)
foo.addStraight(segB,6.0)

foo.addCap(segA)
foo.addCap(segB)

foo.done()

walls = foo.walls

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] -= 6.0
		wall[i] = p

