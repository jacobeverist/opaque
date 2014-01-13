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


" current state probability distribution SxF " 
PI = []
#pi_0 =  np.matrix('0.11 0.11 0.11; 0.11 0.11 0.11 ; 0.11 0.11 0.11 ')
pi_0 = np.zeros([3,3])
for i in range(3):
	for j in range(3):
		pi_0[i,j] = 1.0/9.0
PI.append( pi_0 )

" action model, state transition probabilities "
" row = s_t, column = s_t+1 "

" forward action "
#PHI_f = np.matrix( [[0.2, 0.7, 0.1],
#					[0.1, 0.2, 0.7],
#					[0.1, 0.2, 0.7] ])

PHI_f = np.matrix( [[0.2, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1, 0.0, 0.0],
					[0.0, 0.2, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1, 0.0],
					[0.0, 0.0, 0.2, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1],

					[0.1, 0.0, 0.0, 0.2, 0.0, 0.0, 0.7, 0.0, 0.0],
					[0.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.0, 0.7, 0.0],
					[0.0, 0.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.0, 0.7],

					[0.01, 0.0, 0.0, 0.1, 0.0, 0.0, 0.89, 0.0, 0.0],
					[0.0, 0.01, 0.0, 0.0, 0.1, 0.0, 0.0, 0.89, 0.0],
					[0.0, 0.0, 0.01, 0.0, 0.0, 0.1, 0.0, 0.0, 0.89]])

#PHI_f = np.matrix( [[0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1],
#					[0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1],
#					[0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1],
#
#					[0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7],
#					[0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7],
#					[0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.7, 0.7, 0.7],
#
#					[0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.89, 0.89, 0.89],
#					[0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.89, 0.89, 0.89],
#					[0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.89, 0.89, 0.89]])
for i in range(9):
	newSum = 0.0
	for j in range(9):
		newSum += PHI_f[i,j]
	for j in range(9):
		PHI_f[i,j] = PHI_f[i,j] / newSum



" backward action "
#PHI_b = np.matrix( [[0.7, 0.2, 0.1],
#					[0.7, 0.2, 0.1],
#					[0.2, 0.7, 0.1] ])

PHI_b = np.matrix( [[0.89, 0.0, 0.0, 0.1, 0.0, 0.0, 0.01, 0.0, 0.0],
					[0.0, 0.89, 0.0, 0.0, 0.1, 0.0, 0.0, 0.01, 0.0],
					[0.0, 0.0, 0.89, 0.0, 0.0, 0.1, 0.0, 0.0, 0.01],

					[0.7, 0.0, 0.0, 0.2, 0.0, 0.0, 0.1, 0.0, 0.0],
					[0.0, 0.7, 0.0, 0.0, 0.2, 0.0, 0.0, 0.1, 0.0],
					[0.0, 0.0, 0.7, 0.0, 0.0, 0.2, 0.0, 0.0, 0.1],

					[0.2, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1, 0.0, 0.0],
					[0.0, 0.2, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1, 0.0],
					[0.0, 0.0, 0.2, 0.0, 0.0, 0.7, 0.0, 0.0, 0.1] ])

#PHI_b = np.matrix( [[0.89, 0.89, 0.89, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01],
#					[0.89, 0.89, 0.89, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01],
#					[0.89, 0.89, 0.89, 0.1, 0.1, 0.1, 0.01, 0.01, 0.01],
#
#					[0.7, 0.7, 0.7, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1],
#					[0.7, 0.7, 0.7, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1],
#					[0.7, 0.7, 0.7, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1],
#
#					[0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1],
#					[0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1],
#					[0.2, 0.2, 0.2, 0.7, 0.7, 0.7, 0.1, 0.1, 0.1] ])
for i in range(9):
	newSum = 0.0
	for j in range(9):
		newSum += PHI_b[i,j]
	for j in range(9):
		PHI_b[i,j] = PHI_b[i,j] / newSum


" sensor model, rows = s0f0,s0f1,s0f2,s1f0... columns = 'nothing', 'detected feature' "
" P(z=1 | s1,f2 ) = THETA[1*3 + 2, 1] "

#THETA = np.matrix( [[0.1, 0.9],
#					[1.0, 0.0],
#					[1.0, 0.0],
#					[1.0, 0.0],
#					[0.1, 0.9],
#					[1.0, 0.0],
#					[1.0, 0.0],
#					[1.0, 0.0],
#					[0.1, 0.9] ])

THETA = np.matrix( [[0.1, 0.9],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.1, 0.9],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.9, 0.1],
					[0.1, 0.9] ])

#pi_0 =  np.matrix('0.11 0.11 0.11; 0.11 0.11 0.11 ; 0.11 0.11 0.11 ')
#PI.append( pi_0 )

alphas = []
newAlpha = np.zeros([3,3])
prevState = deepcopy(PI[0])
currZ = Obs[0]
for i in range(3):
	for j in range(3):
		# S_i
		# F_j
		newAlpha[i,j] = prevState[i,j] * THETA[i*3+j,currZ]

alphas.append(newAlpha)




for t in range(1,len(Obs)):
	newAlpha = np.zeros([3,3])
	prevAlpha = deepcopy(alphas[t-1])

	newState = np.zeros([3,3])
	prevState = deepcopy(PI[t-1])

	currAct = Acts[t-1]
	prevZ = Obs[t-1]
	currZ = Obs[t]

	if currAct == 'f':
		denomSum = 0.0
		for i in range(3):  # s_i  (t+1)
			for j in range(3):  # f_j  (t+1)
				newAlphaSum = 0.0
				numerSum = 0.0
				for k in range(3):  # s_k  (t)
					for l in range(3):  # f_l  (t)
						newAlphaSum += prevAlpha[k,l] * PHI_f[k*3+l,i*3+j] * THETA[i*3+j, currZ]

						numerSum += prevState[k,l] * THETA[k*3+l, prevZ] * PHI_f[k*3+l,i*3+j]
						denomSum += prevState[k,l] * THETA[k*3+l, prevZ] * PHI_f[k*3+l,i*3+j] * THETA[i*3+j, currZ]

				newAlpha[i,j] = newAlphaSum 
				newState[i,j] = THETA[i*3+j, currZ] * numerSum

		for i in range(3):  # s_i  (t+1)
			for j in range(3):  # f_j  (t+1)
				newState[i,j] /= denomSum

		print "go forward"

	elif currAct == 'b':
		denomSum = 0.0
		for i in range(3):
			for j in range(3):
				newAlphaSum = 0.0
				numerSum = 0.0
				for k in range(3):
					for l in range(3):
						newAlphaSum += prevAlpha[k,l] * PHI_b[k*3+l,i*3+j] * THETA[i*3+j, currZ]

						numerSum += prevState[k,l] * THETA[k*3+l, prevZ] * PHI_b[k*3+l,i*3+j]
						denomSum += prevState[k,l] * THETA[k*3+l, prevZ] * PHI_b[k*3+l,i*3+j] * THETA[i*3+j, currZ]

				newAlpha[i,j] = newAlphaSum
				newState[i,j] = THETA[i*3+j, currZ] * numerSum

		for i in range(3):  # s_i  (t+1)
			for j in range(3):  # f_j  (t+1)
				newState[i,j] /= denomSum

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
	#print probState
	print probFeature
	#print probSum, newState
	#print newAlpha
	

	#print currAct, currZ, newAlpha, newAlpha[0,0]+newAlpha[0,1]+newAlpha[0,2], newAlpha[1,0]+newAlpha[1,1]+newAlpha[1,2], newAlpha[2,0]+newAlpha[2,1]+newAlpha[2,2] 
	alphas.append(newAlpha)
	PI.append(newState)



" beta initialized to 1.0 "
betas = []
newBeta = np.matrix('1.0 1.0 1.0; 1.0 1.0 1.0 ; 1.0 1.0 1.0 ')
betas.insert(0,newBeta)


times = range(1,len(Obs))
#times = range(len(Acts))
times.reverse()


for t in times:

	currAct = Acts[t-1]
	currZ = Obs[t]

	newBeta = np.zeros([3,3])
	prevBeta = betas[0]

	if currAct == 'f':
		for i in range(3):
			for j in range(3):
				newBetaSum = 0.0
				for k in range(3):
					for l in range(3):
						newBetaSum += PHI_f[i*3+j, k*3+l] * THETA[k*3+l, currZ] * prevBeta[k,l]

				newBeta[i,j] = newBetaSum

	elif currAct == 'b':
		for i in range(3):
			for j in range(3):
				newBetaSum = 0.0
				for k in range(3):
					for l in range(3):
						newBetaSum += PHI_b[i*3+j, k*3+l] * THETA[k*3+l, currZ] * prevBeta[k,l]

				newBeta[i,j] = newBetaSum

	else:
		raise


	betas.insert(0,newBeta)

" experience probability given model and action sequence "
P_O_AMC = 0.0
beta_0 = betas[0]
pi_0 = PI[0]
currZ = Obs[0]

for i in range(3):
	for j in range(3):
		P_O_AMC += pi_0[i,j] * THETA[i*3+j, currZ] * beta_0[i,j]

print len(alphas), len(betas), P_O_AMC


gammas = []
for t in range(len(alphas)):

	currAlpha = alphas[t]
	currBeta = betas[t]

	newGamma = np.zeros([3,3])
	for i in range(3):
		for j in range(3):
			newGamma[i,j] = currAlpha[i,j] * currBeta[i,j] / P_O_AMC

	#print newGamma

	gammas.append(newGamma)

epsilons = []
for t in range(1,len(alphas)):

	currAlpha = alphas[t-1]
	currBeta = betas[t]
	currAct = Acts[t-1]
	currZ = Obs[t]

	newEpsilon = np.zeros([9,9])

	if currAct == 'f':
		for i in range(3):
			for j in range(3):
				newEpsilonSum = 0.0
				for k in range(3):
					for l in range(3):
						newEpsilon[i*3+j, k*3+l] += currAlpha[i,j] * PHI_f[i*3+j, k*3+l] * THETA[k*3+l, currZ] * currBeta[k,l] / P_O_AMC

	elif currAct == 'b':
		for i in range(3):
			for j in range(3):
				for k in range(3):
					for l in range(3):
						newEpsilon[i*3+j, k*3+l] += currAlpha[i,j] * PHI_b[i*3+j, k*3+l] * THETA[k*3+l, currZ] * currBeta[k,l] / P_O_AMC

	else:
		raise

	epsilons.append(newEpsilon)

for t in range(0,len(gammas)):
	currGamma = gammas[t]

for t in range(1,len(alphas)):
	currEpsilon = epsilons[t-1]

new_PHI_f = np.zeros([9,9])
new_PHI_b = np.zeros([9,9])




