
# Single-Cluster Spectral Graph Partitioning for Robotics Applications
# Edwin Olson, Matthew Walter, Seth Teller, and John Leanord
# Proceedings of Robotics Science and Systems, 2005

# Written by Jacob Everist, 2011
# jacob.everist@gmail.com

import numpy
import random
from math import *
from numpy import *

def rand_vector(n):
	vec_n = numpy.matrix([[0.0] for i in range(n)],dtype=float)
	sum = 0.0

	for i in range(n):
		vec_n[i][0] = random.uniform(0,1)
		sum += vec_n[i][0]* vec_n[i][0]
	
	mag = sqrt(sum)
	
	for i in range(n):
		vec_n[i][0] /= mag
	
	return vec_n


" compute first two dominant eigenvectors from power method "
def dominantEigenvectors(A):
	N = A.shape[0]
	#print "N =", N
	K = 5
	#e = [numpy.matrix([],dtype=float) for i in range(N)]
	e = [numpy.matrix([],dtype=float) for i in range(3)]
	
	lmbda = [0.0 for i in range(N)]
	
	#for i in range(N):
	for i in range(2):
		e[i] = rand_vector(N)

		for iters in range(K):
			for j in range(i):
	
				alpha = numpy.dot(e[j].transpose(), e[i])
				alpha = alpha[0,0]
				e[i] = e[i] - alpha * e[j]	
		
			sum = 0.0
			for k in range(N):
				sum += e[i][k]* e[i][k]				
			mag = sqrt(sum)
		
			e[i] = e[i] / mag
			e[i] = A * e[i]
				
			sum = 0.0
			for k in range(N):
				sum += e[i][k]* e[i][k]				
			lmbda[i] = sqrt(sum)
		
		sum = 0.0
		for k in range(N):
			sum += e[i][k]* e[i][k]				
		mag = sqrt(sum)		
		
		e[i] = e[i] / mag
		
	return e, lmbda

def getIndicatorVector(e):

	" create tuples of eigenvector elements and associated hypotheses "
	tups = []
	domVector = e
	N = domVector.shape[0]
	#print "shape =", domVector.shape
	#print "N =", N
	for i in range(N):
		tups.append((domVector[i,0], i))
		
	tups.sort()
	tups.reverse()
	
	
	" maximize the dot product of v and w "
	v = numpy.matrix([[tups[i][0]] for i in range(len(tups))], dtype=float)
	w = numpy.matrix([[0.0] for i in range(len(tups))], dtype=float)
	wp = numpy.matrix([[0.0] for i in range(len(tups))], dtype=float)
	
	lastVal = 0.0
	for i in range(len(tups)):
		w[i,0] = 1
		val = numpy.dot(w.transpose(), v)
		#print "val:", val
		val = val[0,0] / sqrt(1+i)

		#print "val:", tups[i][0], "dot:", val
		
		if val > lastVal:
			lastVal = val
		else:
			w[i,0] = 0
			break

		
		#wp[tups[i][1],0] = 1	
		#val2 = (wp.transpose() * A * wp) / (numpy.dot(wp.transpose(),wp)[0,0])

	return w

if __name__ == '__main__':
		
	" disable the truncation "
	set_printoptions(threshold=nan)
	
	fp1 = open("../../testData/clustering/A_28.txt", 'r')
	A1 = eval(fp1.read())

	fp2 = open("../../testData/clustering/A_110.txt", 'r')
	A2 = eval(fp2.read())
	
	fp3 = open("../../testData/clustering/A_120.txt", 'r')
	A3 = eval(fp3.read())

	fp4 = open("../../testData/clustering/A_overlap.txt", 'r')
	A4 = eval(fp4.read())

	fp5 = open("../../testData/clustering/A_sensor.txt", 'r')
	A5 = eval(fp5.read())
	
	print A1.shape
	print A2.shape
	print A3.shape
	print A4.shape
	print A5.shape

		
	" do graph clustering of consistency matrix "
	w = []
	e = []
	for i in range(100):
		e, lmbda = dominantEigenvectors(A5)
		w = getIndicatorVector(e[1])
		if len(e) <= 1:
			break

		" threshold test l1 / l2"
		ratio = lmbda[0][0,0] / lmbda[1][0,0]
		print "ratio:", ratio
		if ratio >= 2.0:
			print e[0]
			print e[1]
			break
	
	print w
	"""
	e = dominantEigenvectors(A1)
	w = getIndicatorVector(e[0])
	print w

	e = dominantEigenvectors(A2)
	w = getIndicatorVector(e[0])
	print w
	

	e = dominantEigenvectors(A3)
	w = getIndicatorVector(e[0])
	print w

	e = dominantEigenvectors(A4)
	w = getIndicatorVector(e[0])
	print w
	
	e = dominantEigenvectors(A5)
	w = getIndicatorVector(e[0])
	print w
	"""

	#V = numpy.matrix([[0.0]*100]*1000)

	#print repr(A[0])
	#print V	
	
	#print "eigenvectors:", e[0]
	
	#e_w, e_v = numpy.linalg.eig(A)
	#print "eigenvectors2:", e_w
	
	
	#"FIXME: add in threshold test l1 / l2"
	
	#print "Power Method:"
	#print "lambda:", lmbda
	#for i in range(len(hypotheses)):
	#	print e[i]
