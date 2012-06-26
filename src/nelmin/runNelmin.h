#ifdef __cplusplus 
extern "C" {
void doTest(double *matchPairs, int numPairs, double *initGuess, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum); 

void doCost(double *matchPairs, int numPairs, double *initGuess, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset);
}

#else
void doTest(double *matchPairs, int numPairs, double *initGuess, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum); 
void doCost(double *matchPairs, int numPairs, double *initGuess, double *poses_1, double *poses_2, int numPoses, double *resultParam, double *resultSum, double *resultOffset);
#endif


