from opaque.environment.MapGenerate import MapMaker


foo = MapMaker(8.0, 0.4)

#foo.addStraight(0,1.0)
segA, segB = foo.addSymJunction(0, pi/2)
#foo.addStraight(segA,5.0)
#foo.addStraight(segB,5.0)
#foo.addCap(segA)
#foo.addCap(segB)

#foo.done()


foo.addStraight(segA,3.0)
foo.addTurn(segA, -pi/3)
foo.addStraight(segA,2.5)
#foo.addTurn(segA, pi/3)
segC, segD = foo.addSymJunction(segA, pi/2)
foo.addStraight(segC,5.0)
foo.addStraight(segD,5.0)
foo.addCap(segC)
foo.addCap(segD)

#foo.addStraight(segA,5.0)
#foo.addCap(segA)

foo.addStraight(segB,2.5)
foo.addTurn(segB, -pi/3)
foo.addStraight(segB,2.5)
foo.addCap(segB)

"""
segC, segD = foo.addSymJunction(segB, pi/2)
foo.addStraight(segC,5.0)
foo.addStraight(segD,5.0)
foo.addCap(segC)
foo.addCap(segD)
"""

foo.done()

"""
foo.addTurn(0, pi/3)
foo.addStraight(0,1.0)
segA, segB = foo.addSymJunction(0, pi/3)

foo.addStraight(segA,1.0)
foo.addStraight(segB,1.0)

segC, segD = foo.addSymJunction(segA, pi/3)
foo.addStraight(segC,1.0)
foo.addStraight(segD,1.0)

foo.addCap(segB)
foo.addCap(segC)
foo.addCap(segD)
"""


walls = foo.walls

for wall in walls:
	for i in range(len(wall)):
		p = copy(wall[i])
		p[0] -= 6.0
		wall[i] = p

