import os
import sys
dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if not dir in sys.path:
	sys.path.append(dir)

from common import *
from Behavior import *

class BlindConcertinaGait(Behavior):
	
	def __init__(self, probe, direction = False):
		Behavior.__init__(self, probe)
		
		self.numSegs = self.probe.numSegs

		self.setDirection(direction)

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

		# time divider
		self.setTimerAliasing(5)

		self.mask = [1.0 for i in range(self.numSegs-1)]

	def getMask(self):
		return self.mask
	
	def setDirection(self, direction):
		self.direction = direction

	def step(self):
		Behavior.step(self)
		'''
		parameters for sinusoid
		1. freq
		2. joint_amp
		3. per_seg
		4. numSegs
		5. jointValue_i = joint_amp*cos(i*per_seg)

		rudimentary concertina locomotion
		c = gaussian width
		b = mean position of gaussian, along (0,(numSegs-1)*1000)
		x = joint location
		'''
		
		if self.isStep():
		
			self.mask = [1.0 for i in range(self.numSegs-1)]
			self.resetJoints()
	
			isLoop = False
	
			if self.gaitTimer >= 900*self.timeInc:
				self.gaitTimer = 0
				isLoop = True
	
			# sliding position of the first gaussian
			if not self.direction:
				b = (self.gaitTimer)*4*2
			else:
				b = (900*self.timeInc-self.gaitTimer)*4*2
	
			b = b / (100.0*self.timeInc)
	
			# number and location of the gaussian peaks
			numPeaks = 2
			peaks = [0.0 for i in range(numPeaks)]
			peaksRef = [0.0 for i in range(numPeaks)]
			peakWidths = [5.0, 9.0]
	
			mid = numPeaks/2
			for i in range(numPeaks):
				peaks[i] = b + (i-mid)*self.period
				peaksRef[i] = b + (i-mid)*self.period
	
			for i in range(self.numSegs-1):
	
				# joint value
				val = self.joint_amp*cos(i*self.per_seg)
	
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
				self.mask[i] = max_short
	
				# damp the signal according to sliding gaussian
				#if mask[i] >= 0.8:
				#self.probe.setServo(i, val*max)
				self.setJoint(i, val*max)
				#self.joints[i] = val*max
	
			return isLoop
	
