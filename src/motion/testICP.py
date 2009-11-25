#!/usr/bin/python

import unittest
import genICP
import math
from numpy import *

class TestSequenceFunctions(unittest.TestCase):

	def setUp(self):
		pass

	def testTangent(self):

		rotTanVec = genICP.findRotationTangent([0.3,1.5])
		self.assertEqual(rotTanVec,[-1.5,0.3])

	def testLocalNormal(self):


		a_data = []
		for i in range(20):
			a_data.append([1.0+i*0.1,1.0+i*0.1])

		b_data = []
		for i in range(20):
			b_data.append([0.3,1.0+i*0.1])

		norm1 = genICP.findLocalNormal([1.5,1.5],a_data)
		norm2 = genICP.findLocalNormal([0.3,1.5],b_data)
		
		self.assertAlmostEqual(norm1[0], -0.70710678118654757)
		self.assertAlmostEqual(norm1[1], -0.70710678118654757)
		self.assertAlmostEqual(norm2[0], 1.0)
		self.assertAlmostEqual(norm2[1], 0.0)

	def testVectorCovariance(self):
		
		x_var = 0.1
		y_var = 1.0

		cov1 = genICP.computeVectorCovariance([1.0,0.0],x_var,y_var)
		cov2 = genICP.computeVectorCovariance([0.0,1.0],x_var,y_var)
		cov3 = genICP.computeVectorCovariance([math.cos(math.pi/4),math.sin(math.pi/4)],x_var,y_var)
		cov4 = genICP.computeVectorCovariance([math.cos(math.pi/4),-math.sin(math.pi/4)],x_var,y_var)
		
		testCov1 = matrix([[ 0.1,  0. ], [ 0.,   1. ]])
		testCov2 = matrix([[ 1.,   0. ], [ 0.,   0.1]])
		testCov3 = matrix([[ 0.55,  -0.45], [ -0.45,  0.55]])
		testCov4 = matrix([[ 0.55, 0.45], [0.45,  0.55]])

		def compareMatrix(mat1,mat2):
			self.assertAlmostEqual(mat1[0,0],mat2[0,0])
			self.assertAlmostEqual(mat1[0,1],mat2[0,1])
			self.assertAlmostEqual(mat1[1,0],mat2[1,0])
			self.assertAlmostEqual(mat1[1,1],mat2[1,1])
			

		compareMatrix(cov1,testCov1)
		compareMatrix(cov2,testCov2)
		compareMatrix(cov3,testCov3)
		compareMatrix(cov4,testCov4)


	def testMatchError(self):

		Ca = matrix([	[0.1, 0.0],
				[0.0, 1.0]
				])

		Cb = matrix([	[0.2, 0.0],
				[0.0, 0.1]
				])
		 
		a = matrix([	[2.0],
				[0.5]
				])

		b = matrix([	[1.9],
				[0.4]
				])

		offset = [0.0, 0.0, 0.0]


		result = genICP.computeMatchError(offset, a, b, Ca, Cb)

		self.assertAlmostEqual(result,.0424242424242)

	def testFindClosestPoint(self):

		offset = [0.0,0.0,0.0]
		a_data = [[1.0,1.0],[1.1,1.1],[1.2,1.2],[1.3,1.31],[1.4,1.4],[1.51,1.5],[1.6,1.6]]
		b_data = [[0.3,1.0],[0.3,1.1],[0.3,1.2],[0.31,1.3],[0.3,1.4],[0.3,1.5],[0.3,1.6]]

		result, resultDist = genICP.findClosestPointInB(b_data, a_data[3], offset)

		self.assertAlmostEqual(result[0],0.31)
		self.assertAlmostEqual(result[1],1.3)
		self.assertAlmostEqual(resultDist, 0.990050503762)

		offset2 = [0.0,0.0,math.pi/4]

		result2, resultDist2 = genICP.findClosestPointInB(b_data, a_data[3], offset2)
		self.assertAlmostEqual(result2[0],0.3)
		self.assertAlmostEqual(result2[1],1.6)
		self.assertAlmostEqual(resultDist2, 0.393175284342)



if __name__ == '__main__':
    unittest.main()

