#!/usr/bin/python

import sys, os, subprocess
import glob
import shlex
import time
import shutil
import argparse

parser = argparse.ArgumentParser(description='DarkMapper Batch')
parser.add_argument('--startOffset', type=float, default=0.0, help="probe start displacement offset")
parser.add_argument('--maxNumPoses', type=int, default=100, help="max number of poses")

# bendFeature
parser.add_argument('--bendFeature', dest='bendFeature', action='store_true')
parser.set_defaults(bendFeature=False)


args = parser.parse_args()


easyMaps = ['60_left_junction.py', '60_right_junction.py', 'L_89_left_junction.py', 'L_89_right_junction.py']

#hardMaps = ['cross_junction.py', 'T_side_junction.py', 'T_bottom_junction.py', 'Y_junction.py']
hardMaps = ['cross_junction.py', 'T_side_junction.py', 'T_89_junction.py', 'Y_junction.py']

comboMaps = ['Y_Y_junction.py', 'Y_right_Y_junction.py', 'Y_right_cross_junction.py', 'Y_left_Y_junction.py', 'T_side_T_junction.py', 'T_bottom_T_junction.py', 'L_right_left_junction.py', 'L_left_right_junction.py', '60_right_left_junction.py', '60_left_right_junction.py', 'cross_cross_junction.py', 'gen_function.py', 'gen_function2.py', 'gen_junction.py', 'gen_junction2.py']

#maps = [hardMaps[0],] + easyMaps + hardMaps[1:]
#maps = easyMaps
maps = hardMaps


"""
1. cycle through each test
2. run test with given map file
3. clean up and store test results
4. go to next test

"""

startOffset = args.startOffset
maxNumPoses = args.maxNumPoses

for mapFile in maps:

	" run the test "
	fileOut = open('out.txt', 'a')
	#command_line = "/usr/bin/python ../startSim.py --mapfile=../mapLibrary/%s --restoreConfig=True --hideWindow=True" % mapFile 
	if args.bendFeature:
		command_line = "/usr/bin/python ../startSim.py --hideWindow=True --restoreConfig=True --mapFile=../mapLibrary/%s --maxNumPoses=%d --startOffset=%1.2f --bendFeature" % (mapFile, maxNumPoses, startOffset)
	else:
		command_line = "/usr/bin/python ../startSim.py --hideWindow=True --restoreConfig=True --mapFile=../mapLibrary/%s --maxNumPoses=%d --startOffset=%f" % (mapFile, maxNumPoses, startOffset)

	
	parsed_args = shlex.split(command_line)
	print "args:", parsed_args
	#p = subprocess.Popen(args, stdout=fileOut, stderr=fileOut)

	subprocess.check_call(parsed_args, stdout=fileOut, stderr=fileOut)

	"paste results together"
	#os.system("python compileResults3.py")
		
	" mkdir "
	testDir = "results_%s" % mapFile
	os.mkdir(testDir)

	outFiles = glob.glob(os.path.join("*.out"))
	outFiles += glob.glob(os.path.join("*.err"))
	outFiles += glob.glob(os.path.join("*.txt"))
	outFiles += glob.glob(os.path.join("*.obj"))
	outFiles += glob.glob(os.path.join("*.png"))


	" move files "
	for outFile in outFiles:
		shutil.move(outFile, testDir)
		
	
