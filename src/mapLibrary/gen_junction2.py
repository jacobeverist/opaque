from opaque.environment.MapGenerate import MapMaker


foo = MapMaker(8.0, 0.4)

segA, segB = foo.addSymJunction(0, pi/3)

foo.addStraight(segA,3.0)
foo.addTurn(segA, pi/3)
foo.addStraight(segA,2.5)
segC, segD = foo.addSymJunction(segA, pi/3)


foo.addStraight(segC,4.0)
foo.addTurn(segC, -pi/3)
foo.addStraight(segC,2.5)
foo.addCap(segC)

foo.addStraight(segD,4.0)
foo.addTurn(segD, pi/3)
foo.addStraight(segD,2.5)
foo.addCap(segD)

foo.addStraight(segB,4.0)
foo.addTurn(segB, -pi/3)
foo.addStraight(segB,2.5)
foo.addCap(segB)

foo.done()

walls = foo.walls

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] -= 6.0
		wall[i] = p

