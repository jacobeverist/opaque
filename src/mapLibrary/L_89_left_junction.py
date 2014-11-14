from opaque.environment.MapGenerate import MapMaker


foo = MapMaker(8.0, 0.4)
foo.addTurn(0, -49.0*pi/100.0)
foo.addStraight(0,6.0)
foo.addCap(0)

foo.done()

walls = foo.walls

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] -= 6.0
		wall[i] = p

