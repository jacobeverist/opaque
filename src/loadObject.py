#!/usr/bin/python

from math import *
from copy import copy
import random
random.seed(0)

import os
import sys

relPath = os.path.dirname(os.path.realpath(__file__))
print relPath

sys.path.insert(1,relPath)
sys.path.insert(1,relPath + "/modules/")
sys.path.insert(1,relPath + "/modules/nelmin")
sys.path.insert(1,relPath + "/modules/bulletprobe")
sys.path.insert(1,relPath + "/modules/medialaxis")
sys.path.insert(1,relPath + "/modules/alphashape")
sys.path.insert(1,relPath + "/modules/ogre")

#os.environ['PATH'] = 'c:\\opencv-2.4.6.0\\build\\x86\\vc9\\bin' + ';' + os.environ['PATH']
#os.environ['PATH'] = os.getcwd() + '\\binaries\\boost_1_53_0\\lib' + ';' + os.environ['PATH']

import argparse

parser = argparse.ArgumentParser(description='DarkMapper Simulator')
parser.add_argument('--mapfile', type=str, help="text file containing environment map")
parser.add_argument('--restoreConfig', type=bool, default=False, help="restore ogre config from file instead of diplaying dialog")
parser.add_argument('--hideWindow', type=bool, default=False, help="hide Ogre window")
args = parser.parse_args()



from bulletsnakeprobe import BulletSnakeProbe
from ogreprobe import ProbeApp

from opaque.robot.QuickProbe import QuickProbe
from opaque.robot.BulletProbe import BulletProbe

from opaque.TestNavigation import TestNavigation
from opaque.TestProcess import TestProcess
#from opaque.TestUnit import TestUnit
from opaque.ControlError import *
from opaque.CmdDrawThings import CmdDrawThings
from opaque.DrawThings import DrawThings



import pickle

dirName = "testDir5"
pickleFile = 'map_0019.obj'


mapAlgorithm = None

with open(dirName + '/' + pickleFile, 'rb') as inputVal:
		
	mapAlgorithm = pickle.load(inputVal)


spatialFeatures = mapAlgorithm.poseData.spatialFeatures

for key, val in spatialFeatures.iteritems():

	resultObj = val[0]

	print "node", key, "has:"
	print "bloomPoint:", resultObj["bloomPoint"]
	print "archPoint:", resultObj["archPoint"]
	print "inflectionPoint:", resultObj["inflectionPoint"]

#for key, val in spatialFeatures.iteritems():

if False:

	key = 0
	val = spatialFeatures[key]

	print key
	#, val

	print len(val)

	print val[0].keys()

	#for item in val:
	#	print item


#print spatialFeatures[0].keys()


