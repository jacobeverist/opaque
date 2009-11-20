import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
    sys.path.append(dir)

from common import *
import robot

#snake = robot.SnakeProbe()