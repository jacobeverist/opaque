"""

Test using Expectation-Maximization (EM) to localize the robot and a single feature in a discrete 3-state world.

January 2014
Jacob Everist


Z = { 1 : 'detected feature', 0 : 'nothing' }
B = {'f' : 'forward', 'b' : 'backward'}

State space 
S = { 0 : 's0', 1 : 's1', 2 : 's2' }   # 3 possible locations for the robot 
F = { 0 : 'f0', 1 : 'f1', 2 : 'f2' }   # 3 possible locations for the feature, doesn't move

"""

import numpy as np
from copy import *


" Experience "
#Obs = (   0,   0,   1,   1,   0,   0,   1,   1,   0,   1)
#Acts = (     'f', 'f', 'f', 'b', 'b', 'f', 'f', 'b', 'f')

Obs = (   0,   0,   1,   1,   0,   1,   0,   1,   1,   1,   1,   1)
Acts = (     'f', 'f', 'f', 'b', 'f', 'b', 'f', 'f', 'f', 'f', 'f')


" current state probability distribution S x F " 
PI = []
pi_0 = np.zeros([3,3])
for i in range(3):
	for j in range(3):
		pi_0[i,j] = 1.0/9.0
PI.append( pi_0 )

" action model, state transition probabilities "
" row = s_t, column = s_t+1 "

" forward action "
PHI_f = np.matrix( [[0.2, 0.7, 0.1],
		[0.1, 0.2, 0.7],
		[0.01, 0.1, 0.89] ])

" backward action "
PHI_b = np.matrix( [[0.89, 0.1, 0.01],
		[0.7, 0.2, 0.1],
		[0.2, 0.7, 0.1] ])


" sensor model, rows = s0f0,s0f1,s0f2,s1f0... columns = 'nothing', 'detected feature' "
" P(z=1 | s1,f2 ) = THETA[1*3 + 2, 1] "

THETA_1 = np.matrix( [[0.9, 0.1, 0.1],
			[0.1, 0.9, 0.1],
			[0.1, 0.1, 0.9] ])

THETA_0 = np.matrix( [[0.1, 0.9, 0.9],
			[0.9, 0.1, 0.9],
			[0.9, 0.9, 0.1] ])


THETA = np.matrix( [[0.1, 0.9],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.1, 0.9],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.1, 0.9] ])

prevState = deepcopy(PI[0])
currZ = Obs[0]

for t in range(1,len(Obs)):

	newState = np.zeros([3,3])
	prevState = deepcopy(PI[t-1])

	currAct = Acts[t-1]
	prevZ = Obs[t-1]
	currZ = Obs[t]

	if currZ == 1:
		curr_THETA = THETA_1
	elif currZ == 0:
		curr_THETA = THETA_0

	if prevZ == 1:
		prev_THETA = THETA_1
	elif prevZ == 0:
		prev_THETA = THETA_0

	if currAct == 'f':
		denomSum = 0.0
		for s_i in range(3):  # s_i  (t+1)
			for f_j in range(3):  # f_j  (t+1)
				numerSum = 0.0
				for s_k in range(3):  # s_k  (t)
					for f_l in range(3):  # f_l  (t)

						# no chance that the action will change the location of the feature
						if f_j == f_l:
							numerSum += prevState[s_k,f_l] * prev_THETA[s_k,f_l] * PHI_f[s_k,s_i]
							denomSum += prevState[s_k,f_l] * prev_THETA[s_k,f_l] * PHI_f[s_k,s_i] * curr_THETA[s_i,f_j]

				newState[s_i,f_j] = curr_THETA[s_i,f_j] * numerSum

		for s_i in range(3):  # s_i  (t+1)
			for f_j in range(3):  # f_j  (t+1)
				newState[s_i,f_j] /= denomSum

		print "go forward"

	elif currAct == 'b':
		denomSum = 0.0
		for s_i in range(3):  # s_i  (t+1)
			for f_j in range(3):  # f_j  (t+1)
				numerSum = 0.0
				for s_k in range(3):  # s_k  (t)
					for f_l in range(3):  # f_l  (t)

						# no chance that the action will change the location of the feature
						if f_j == f_l:
							numerSum += prevState[s_k,f_l] * prev_THETA[s_k,f_l] * PHI_b[s_k,s_i]
							denomSum += prevState[s_k,f_l] * prev_THETA[s_k,f_l] * PHI_b[s_k,s_i] * curr_THETA[s_i,f_j]

				newState[s_i,f_j] = curr_THETA[s_i,f_j] * numerSum

		for s_i in range(3):  # s_i  (t+1)
			for f_j in range(3):  # f_j  (t+1)
				newState[s_i,f_j] /= denomSum

		print "go backward"
	else:
		raise


	probState = [0.0,0.0,0.0]
	probFeature = [0.0,0.0,0.0]

	for i in range(3):
		probState[i] = newState[i,0]+newState[i,1]+newState[i,2]
		probFeature[i] = newState[0,i]+newState[1,i]+newState[2,i]

	probSum = 0.0
	for i in range(3):  
		for j in range(3): 
			probSum += newState[i,j] 


	print currAct, currZ
	print newState
	#print probState
	#print probFeature
	PI.append(newState)

