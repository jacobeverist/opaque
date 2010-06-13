#!/usr/bin/python

import sys, os

"""
1. cycle through each test
2. run test with given map file
3. clean up and store test results
4. go to next test

"""

for i in range(0,17):
	
	for j in range(0,4):
		
		dist = 0.2*j + 0.5
		
		" run the test "
		os.system("python StartTest.py " + "completeFile%04u.txt " % i + str(dist))
		
		"paste results together"
		#os.system("python compileResults3.py")
		
		" mkdir "
		testDir = "results%02u_%f" % (i,dist)
		os.system("mkdir " + testDir)
		
		" move files "
		os.system("mv scene*.png " + testDir )
		os.system("mv poses*.txt " + testDir )
	
		" clean up "
		pass
	

