
// number of servo position samples for stability determinatin
const int SAMPLES = 20;

class ValueStability {
public:
	int sampleCount, sampleIndex, sampleSize;
	double sampleMean, sampleVar, varThresh;
	int samples[SAMPLES];
    ValueStability(double thresh, int sample_size);
    ~ValueStability();

	void setThresh(double val);
	int isStable();
	double getMean();
	double getVar();
	int getSampleCount();
	void reset();
	int addData(double newValue);
};



