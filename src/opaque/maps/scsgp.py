
# Single-Cluster Spectral Graph Partitioning for Robotics Applications
# Edwin Olson, Matthew Walter, Seth Teller, and John Leanord
# Proceedings of Robotics Science and Systems, 2005

# Written by Jacob Everist, 2011
# jacob.everist@gmail.com

import numpy
import random
from math import *
from numpy import *
from scipy.sparse.linalg import eigsh, eigs



def poweig(A, x0, maxiter = 100, ztol= 1.0e-5, mode= 0, teststeps=1):
	"""
	Performs iterative power method for dominant eigenvalue.
	 A  - input matrix.
	 x0 - initial estimate vector.
	 maxiter - maximum iterations
	 ztol - zero comparison.
	 mode:
	   0 - divide by last nonzero element.
	   1 - unitize.
	Return value:
	 eigenvalue, eigenvector
	"""

	"""
	A = [[3,6],[1, 4]]
	x = [1,1]
	xlambda, x = poweig(A,x, mode =1, teststeps=1)
	print "poweig returns", xlambda, x
	"""
	
	m    = len(A)
	xi   = x0[:] 
 
	for n in range(maxiter):
		# matrix vector multiplication.
		xim1 = xi[:]
		for i in range(m):
			xi[i] = 0.0
			for j in range(m):
				xi[i] += A[i][j] * xim1[j]
		print n, xi
		if mode == 0:
			vlen = sqrt(sum([xi[k]**2 for k in range(m)]))
			xi = [xi[k] /vlen for k in range(m)]
		elif mode == 1:
			for k in range(m-1, -1, -1):
				c = abs(xi[k])
				if c > 1.0e-5:
					xi = [xi[k] /c for k in range(m)]
				break
		# early termination test.
		if n % teststeps == 0:
			S = sum([xi[k]-xim1[k] for k in range(m)])
			if abs(S) < ztol:
				break
		print n, xi
		
	# Compute Rayleigh quotient.
	numer = sum([xi[k] * xim1[k] for k in range(m)])
	denom = sum([xim1[k]**2 for k in range(m)])
	xlambda = numer/denom
	return xlambda, xi


def rand_vector(n):
	vec_n = numpy.matrix([[0.0] for i in range(n)],dtype=float64)
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
	e = [numpy.matrix([],dtype=float64) for i in range(2)]
	
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

def dominantEigenvectors2(A):
	N = A.shape[0]
	#print "N =", N
	K = 5
	#e = [numpy.matrix([],dtype=float) for i in range(N)]
	lmbda = [0.0 for i in range(N)]


	#lmbda2, e2 = eigs(A, k=2)
	#print lmbda2
	#print e2
	

	lmbda, e = eigsh(A, k=2)
	#print lmbda
	#print e

	e0 = []
	e1 = []
	l0 = 0.0
	l1 = 0.0
	
	if fabs(lmbda[0]) > fabs(lmbda[1]):
		l0 = lmbda[0]
		l1 = lmbda[1]
		for k in range(N):
			e0.append(e[k][0])
			e1.append(e[k][1])
	else:
		l0 = lmbda[1]
		l1 = lmbda[0]
		for k in range(N):
			e0.append(e[k][1])
			e1.append(e[k][0])

	lmbda = [0.0 for i in range(N)]
	lmbda[0] = numpy.matrix([l0])
	lmbda[1] = numpy.matrix([l1])
	
	#print lmbda2
	
	#print 
	#print e2
	
	e_mat = [numpy.matrix([0.0],dtype=float64) for i in range(2)]
	e0_mat = numpy.matrix([[e0[j]] for j in range(N)],dtype=float)
	e1_mat = numpy.matrix([[e1[j]] for j in range(N)],dtype=float)
	
	e_mat[0] = e0_mat
	e_mat[1] = e1_mat
	
	#print lmbda
	#print e_mat
	
	return e_mat, lmbda


def getIndicatorVector(e):

	" create tuples of eigenvector elements and associated hypotheses "
	tups = []
	domVector = e
	N = domVector.shape[0]
	#print "shape =", domVector.shape
	#print "N =", N
	for i in range(N):
		tups.append((fabs(domVector[i,0]), i))
		
	tups.sort()
	tups.reverse()
	
	
	" maximize the dot product of v and w "
	#v = numpy.matrix([[tups[i][0]] for i in range(len(tups))], dtype=float)
	v = numpy.matrix([[fabs(domVector[i,0])] for i in range(len(tups))], dtype=float)
	w = numpy.matrix([[0.0] for i in range(len(tups))], dtype=float)
	wp = numpy.matrix([[0.0] for i in range(len(tups))], dtype=float)
	
	indVec = [0 for i in range(len(tups))]
	
	#print "tups:", tups
	#print "vec:", domVector
	
	lastVal = 0.0
	for i in range(len(tups)):
		w[tups[i][1],0] = 1
		indVec[tups[i][1]] = 1
		#w[i,0] = 1
		
		for k in range(len(tups)):
			if indVec[k] == 1:
				wp[k,0] = 1.0 / sqrt(i+1)	
			else:
				wp[k,0] = 0.0	
		
		val = numpy.dot(wp.transpose(), v)
		#val = numpy.dot(w.transpose(), v)
		#print "val:", val
		#val = val[0,0] / (1+i)
		val = val[0,0]

		#print "w:", indVec
		#print "dot:", val
		
		if val > lastVal:
			lastVal = val
		else:
			w[tups[i][1],0] = 0
			indVec[tups[i][1]] = 0
			break

		
		#wp[tups[i][1],0] = 1	
		#val2 = (wp.transpose() * A * wp) / (numpy.dot(wp.transpose(),wp)[0,0])

	eigVec = [domVector[i,0] for i in range(len(tups))]
	uniVec = [wp[i,0] for i in range(len(tups))]

	#print
	#print "eigVec:", eigVec
	#print "w:", indVec
	#print "wp:", uniVec
	#print 
	
	return w

if __name__ == '__main__':
		
	" disable the truncation "
	set_printoptions(threshold=nan)
	
	N = 10
	
	" create the consistency matrix "
	A = matrix(N*N*[0.0], dtype=float)
	A.resize(N, N)
	
	" populate consistency matrix "			
	for i in range(N):
		for j in range(N):
			if [0,1,2].count(i) > 0 and [0,1,2].count(j) > 0:
				A[i,j] = 1.0
				A[j,i] = 1.0
			else:
				A[i,j] = 0.1
				A[j,i] = 0.1

	A[6,5] = 0.5
	A[5,6] = 0.5
	A[7,5] = 0.5
	A[5,7] = 0.5
	A[7,7] = 0.5
	A[6,7] = 0.5
	
	#print A
	print "A:"
	for i in range(N):
		printStr = ""
		for j in range(N):
			printStr += "%1.2f" % (A[i,j]) + " "
		
		print printStr

	" do graph clustering of consistency matrix "
	w = []
	e = []
	for i in range(100):

		e, lmbda = dominantEigenvectors2(A)

		print "lambda:", lmbda
		eigVec0 = [e[0][k,0] for k in range(N)]
		eigVec1 = [e[1][k,0] for k in range(N)]
		print "eigVec0:", eigVec0
		print "eigVec1:", eigVec1


		w0 = getIndicatorVector(e[0])
		w1 = getIndicatorVector(e[1])
		wVec0 = [w0[k,0] for k in range(N)]
		wVec1 = [w1[k,0] for k in range(N)]
		print "wVec0:", wVec0
		print "wVec1:", wVec1
		
		
		if len(e) <= 1:
			break

		" threshold test l1 / l2"
		ratio = fabs(lmbda[0][0,0]) / fabs(lmbda[1][0,0])
		print "ratio:", ratio
		if ratio >= 2.0:
			break
	
	
	exit()
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
