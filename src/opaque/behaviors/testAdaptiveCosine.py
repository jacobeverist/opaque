
from AdaptiveCosine import * 
from AdaptiveAnchorCurve import * 
import pylab
from math import *

"""
curve = AdaptiveCosine(4*pi, 1.0)

curve.setPeakAmp(0, 0.5)
curve.setPeakAmp(1, 0.5)
curve.setHeadLength(0.5)

# 3.65, -0.5

curve.draw()
pylab.show()
"""


curve = AdaptiveAnchorCurve(4*pi)

curve.setPeakAmp(0.5)
curve.draw()
pylab.show()
