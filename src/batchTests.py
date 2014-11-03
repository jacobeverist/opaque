#!/usr/bin/python

import sys, os

maps = ['60_left_junction.py', '60_right_junction.py', 'cross_junction.py', 'L_left_junction.py', 'L_right_junction.py', 'T_side_junction.py', 'T_bottom_junction.py', 'Y_junction.py']

moreMaps = ['Y_Y_junction.py', 'Y_right_Y_junction.py', 'Y_right_cross_junction.py', 'Y_left_Y_junction.py', 'T_side_T_junction.py', 'T_bottom_T_junction.py', 'L_right_left_junction.py', 'L_left_right_junction.py', '60_right_left_junction.py', '60_left_right_junction.py', 'cross_cross_junction.py', 'gen_function.py', 'gen_function2.py']




# python startSim.py --mapfile=../mapLibrary/cross_junction.py --restoreConfig=True --hideWindow=True &> out.txt



"""
1. cycle through each test
2. run test with given map file
3. clean up and store test results
4. go to next test

"""

for mapFile in maps:

	" run the test "
	os.system("python ../startSim.py --mapfile=../mapLibrary/%s --restoreConfig=True --hideWindow=True &> out.txt" % mapFile)
		
	"paste results together"
	#os.system("python compileResults3.py")
		
	" mkdir "
	testDir = "results_%s" % mapFile
	os.system("mkdir " + testDir)
		
	" move files "
	os.system("mv *.png *.obj *.txt *.out *.err " + testDir )
	
	" clean up "
	pass

