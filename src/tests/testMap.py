from MapGenerate import *


#foo = MapMaker(4.0, 0.4)

foo = MapMaker(8.0, 0.4)

#foo.addStraight(0,1.0)
segA, segB = foo.addSymJunction(0, pi/2)
foo.addStraight(segA,0.1)
foo.addStraight(segB,0.1)
#foo.addCap(segA)
#foo.addCap(segB)

foo.done()

"""

#foo.addTurn(0, -pi/2)
foo.addStraight(0,1.0)
foo.addTurn(0, pi/3)
foo.addStraight(0,1.0)
#foo.addTurn(0, -pi/3)
#foo.addStraight(0,1.0)
#foo.addTurn(0, -pi/3)
#foo.addStraight(0,1.0)
#foo.addTurn(0, pi/3)
#foo.addStraight(0,1.0)

#segA, segB = foo.addSymJunction(0, pi/2-pi/16)
#segA, segB = foo.addSymJunction(0, pi/2)
segA, segB = foo.addSymJunction(0, pi/3)

foo.addStraight(segA,1.0)
foo.addStraight(segB,1.0)

segC, segD = foo.addSymJunction(segA, pi/3)
foo.addStraight(segC,1.0)
foo.addStraight(segD,1.0)

foo.addCap(segB)
foo.addCap(segC)
foo.addCap(segD)

foo.done()

"""

"""
foo.addTurn(0, pi/3)
foo.addStraight(0,1.0)
foo.addTurn(0, -pi/3)
foo.addStraight(0,1.0)
foo.addTurn(0, pi/4)
foo.addStraight(0,1.0)
foo.addTurn(0, -pi/2)
foo.addStraight(0,1.0)
foo.addTurn(0, pi/4)
foo.addStraight(0,1.0)
foo.addTurn(0, pi/4)
"""

foo.draw()
