#!/usr/bin/python

import sys, os

"""
1. cycle through each test
2. run test with given map file
3. clean up and store test results
4. go to next test

"""

os.system("rm convexHull*.png")
os.system("rm gndLocal*.png")
os.system("rm localOcc*.png")
os.system("rm costFile*.txt")
os.system("rm estpose*.txt")
os.system("rm prof*.txt")
os.system("rm comparison*.png")

for i in range(1,8):

	" run the test "
	os.system("python StartSimulation.py " + "completeFile%04u.txt" % i)
	
	"paste results together"
	os.system("python compileResults3.py")
	
	" mkdir "
	testDir = "results%02u" % i
	os.system("mkdir " + testDir)
	
	" move files "
	os.system("mv convexHull*.png " + testDir )
	os.system("mv gndLocal*.png " + testDir)
	os.system("mv localOcc*.png " + testDir)
	os.system("mv costFile*.txt " + testDir)
	os.system("mv estpose*.txt " + testDir)
	os.system("mv prof*.txt " + testDir)
	os.system("mv comparison*.png " + testDir)
	
	" clean up "
	pass


