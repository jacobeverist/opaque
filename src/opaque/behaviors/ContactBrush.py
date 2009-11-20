import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *

class ContactBrush(Behavior):
	
	def __init__(self, probe, direction = False ):
		Behavior.__init__(self, probe)

		self.dir = direction

	def step(self):

		return False
	
	def getPoints():
		return []