#!/usr/bin/python
import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from copy import *
from math import *
from random import *
import graph
import pylab

import Image
from maps import *


if __name__ =="__main__":

	points = [[3.0146429223689015, -4.6360507807281364], [2.926984730994624, -4.2295197291351778], [3.037630958186071, -3.4562473086890697], [2.9074070191177088, -3.4197869808984453], [2.8333637708991524, -3.0208687236050422], [2.7794573796539588, -2.675103519280559], [2.7620356913597419, -2.3322100488592468], [2.6309896222468057, -2.3209970331147014], [2.5359292123289596, -2.0898562324227128], [2.5359292123289596, -2.0898562324227128]]
	
	points = [[-0.30290490388870239, 0.18356227874755859], [-0.26424437761306763, 0.38780254125595093], [-0.55323827266693115, 0.66320997476577759]]

	points2 = []
	for p in points:
		if p not in points2:
			points2.append(p)

	spline1 = SplineFit(points2, smooth = 0.1, kp = 2)
	
	xP = []
	yP = []
	for p in points2:
		print p
		xP.append(p[0])
		yP.append(p[1])
		
	pylab.scatter(xP,yP, color='0.0')
	pylab.plot(xP,yP, color='0.0')

	spline1.drawSpline()
				
	pylab.show()

	
	
	
	
	