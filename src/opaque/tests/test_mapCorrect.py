import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *

import Image
import ImageDraw
from copy import *
from math import *
from numpy import arange
from maps import *
#from AlphaMapCost import AlphaMapCost

	
if __name__ =="__main__":
	alphamap = AlphaMap()
	
	pass


