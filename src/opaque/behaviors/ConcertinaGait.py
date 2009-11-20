import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *


from SnakeControl import *
from copy import *
from math import *

import opaque.behaviors as behave
import opaque.pose as poses
import opaque.maps as maps

# super states
A = 0
B = 1
C = 2
D = 3
E = 4

class ConcertinaGait:

	# 3 Steps, go right path, go left path, go back to origin
	# compute path at each step

	def __init__(self, probe):
		self.probe = probe
		self.gaitTimer = 0
		self.childNodes = []
		self.refJoint = 19
		self.numSegs = self.probe.numSegs
		
		# frequency, number of cycles per unit snake length
		self.freq = 5.0

		# joint amplitude
		# working freq=5.0,joint_amp=70.0, c=7
		self.joint_amp = 70.0

		# angle increment per segment
		self.per_seg = 2*pi / (self.numSegs / self.freq)
		
		# distance between peaks along the snake from 0.0 to 1.0
		self.peak_dist = 0.90
		
		# amplification factor to make the gaussian taller and sharper
		self.amp_factor = 1.5
		
		# width of the Gaussian
		self.c_short = 4.0
		
		# period
		# distance between each gaussian peak determined by frequency
		self.period = self.numSegs * self.peak_dist

		
	def doBehavior(self):
		self.doCautiousGait(True, False)
		
	def doCautiousGait(self, is_moving, is_pause, dir = False):
		'''
		parameters for sinusoid
		1. freq
		2. joint_amp
		3. per_seg
		4. numSegs
		5. jointValue_i = joint_amp*cos(i*per_seg)

		rudimentary concertina locomotion
		c = gaussian width
		b = mean position of gaussian, along (0,(numSegs-1)*100)
		x = joint location
		'''

		self.gaitTimer += 1
		mask = [0.0 for i in range(self.numSegs-1)]

		isLoop = False

		if self.gaitTimer > 900:
			self.gaitTimer = 0
			isLoop = True

		# sliding position of the first gaussian
		if not dir:
			b = (self.gaitTimer)*4*2
		else:
			b = (900-self.gaitTimer)*4*2

		b = b / (100.0)

		# number and location of the gaussian peaks
		numPeaks = 2
		peaks = [0.0 for i in range(numPeaks)]
		peaksRef = [0.0 for i in range(numPeaks)]
		peakWidths = [5.0, 9.0]

		mid = numPeaks/2
		for i in range(numPeaks):
			peaks[i] = b + (i-mid)*period
			peaksRef[i] = b + (i-mid)*period

		for i in range(self.numSegs-1):

			# joint value
			val = self.joint_amp*cos(i*self.per_seg)

			if is_moving:
				x = i
				max_short = 0.0

				if i < peaksRef[0] or i > peaksRef[1] :
					max_short = 1.0

				else:
					# mask for a shorter widthed activity indicator
					amps_short = range(numPeaks)

					# active node gaussian peaks
					for j in range(numPeaks):
						amps_short[j] = exp(-(x-peaksRef[j])*(x-peaksRef[j])/(2*self.c_short*self.c_short))
						amps_short[j] *= self.amp_factor

					for j in range(numPeaks):
						if amps_short[j] > max_short:
							max_short = amps_short[j]

					if max_short > 1.0:
						max_short = 1.0

				max = 0.0

				if i < peaks[0] or i > peaks[1] :
					max = 1.0

				else:
					# gaussians for amplitude and activity indicator
					amps = range(numPeaks)

					# active node gaussian peaks
					for j in range(numPeaks):
						amps[j] = exp(-(x-peaks[j])*(x-peaks[j])/(2*peakWidths[j]*peakWidths[j]))
						amps[j] *= self.amp_factor

					for j in range(numPeaks):
						if amps[j] > max:
							max = amps[j]

					if max > 1.0:
						max = 1.0

				# mask for determining active reference points
				mask[i] = max_short

				# damp the signal according to sliding gaussian
				#if mask[i] >= 0.8:
				self.probe.setServo(i, val*max)
					

			else:
				self.probe.setServo(i, val)

		return isLoop, mask


