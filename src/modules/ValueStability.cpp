#include "ValueStability.h"


ValueStability::ValueStability(double thresh, int sample_size) {

	sampleMean = 0.0;
	sampleVar = 0.0;
	sampleCount = 0;
	sampleIndex = 0;
	varThresh = thresh;
	//sampleSize = sample_size;
	sampleSize = SAMPLES;
	
	//samples = new int(sampleSize);

	for ( int i = 0 ; i < sampleSize ; i++ ) {
		samples[i] = 0.0;
	}
	
}

ValueStability::~ValueStability() {
	//delete samples;
}

void ValueStability::setThresh(double val) {
	varThresh = val;
}

int ValueStability::isStable() { 
	if (sampleCount < sampleSize) 
		return 0;

	if (sampleVar < varThresh)
		return 1;
	
	return 0;
}

double ValueStability::getMean() {
	return sampleMean;
}

double ValueStability::getVar() {
	return sampleVar;
}

int ValueStability::getSampleCount() {
	return sampleCount;
}

void ValueStability::reset() {
	sampleCount = 0;
	sampleIndex = 0;
	sampleMean = 0.0;
	sampleVar = 0.0;
}

int ValueStability::addData(double newValue) {

	// **** COMPUTE SAMPLE VARIANCE AND MEAN  ****
	//	
	// NOTE: Samples are taken at every time step and are used to determine
	// the stability of the probe position.  Low variance indicates that the
	// position is stable and it is safe to compute a feature point.
	//	

	// retrieve data sample
	samples[sampleIndex] = newValue;

	// increment the circular index
	sampleIndex += 1;
	sampleIndex = sampleIndex % sampleSize;

	if (sampleCount < sampleSize)
		sampleCount += 1;

	// compute the mean
	double sampleSum = 0.0;
	for (int i = 0 ; i < sampleCount ; i++ ) 
		sampleSum += samples[i];
	sampleMean = sampleSum / sampleCount;

	// compute the variance
	sampleSum = 0;
	for (int i = 0 ; i < sampleCount ; i++ ) 
		sampleSum += (samples[i] - sampleMean)*(samples[i] - sampleMean);
	sampleVar = sampleSum/sampleCount;

	//max = -1e200
	//min = 1e200
	//for i in range(0,sampleCount):
	//	if samples[i] > max:
	//		max = samples[i]
	//	if samples[i] < min:
	//		min = samples[i]

	// should have at least 10 values before we compute the variance
	if (sampleCount < sampleSize)
		// make extremely high for invalid variances
		sampleVar = 1e200;

	if (sampleVar < varThresh)
		return 1;

	return 0;
}

