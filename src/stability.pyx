cdef extern from "ValueStability.h":
	ctypedef struct c_ValueStability "ValueStability":

		int sampleCount, sampleIndex, sampleSize
		double sampleMean, sampleVar, varThresh
		int *samples

		void setThresh(double val)
		int isStable()
		double getMean()
		double getVar()
		int getSampleCount()
		void reset()
		int addData(double newValue)

	c_ValueStability *new_ValueStability "new ValueStability" (double thresh, int sample_size)
	void del_ValueStability "delete" (c_ValueStability *valStab)

cdef class ValueStability:
	cdef c_ValueStability *thisptr      # hold a C++ instance which we're wrapping
	def __cinit__(self, double thresh, int sample_size):
		self.thisptr = new_ValueStability(thresh, sample_size)
	def __dealloc__(self):
		del_ValueStability(self.thisptr)

	def setThresh(self, val):
		self.thisptr.setThresh(val)

	def isStable(self):
		return self.thisptr.isStable()

	def getMean(self):
		return self.thisptr.getMean()

	def getVar(self):
		return self.thisptr.getVar()

	def getSampleCount(self):
		return self.thisptr.getSampleCount()

	def reset(self):
		self.thisptr.reset()

	def addData(self, newValue):
		return self.thisptr.addData(newValue)
