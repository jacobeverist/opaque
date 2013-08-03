
from alphamod import doAlpha

def testAlpha():

	" alpha shape circle radius "

	radius = 0.2
	perturbPoints = []

	#targetInputFile = "doAlphaInput_00000000.txt"
	targetInputFile = "doAlphaInput_00000252.txt"
	f = open(targetInputFile, 'r')
	str_f = f.read()
	exec(str_f)

	vertices = doAlpha(radius,perturbPoints)

	#print vertices

testAlpha()
#while True:
#	testAlpha()
