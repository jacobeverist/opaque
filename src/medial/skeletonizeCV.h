#include <stdio.h>
#include <math.h>
#include <OpenCV/cv.h>
#include <OpenCV/highgui.h>

int runAlg(int index, int sizeX, int sizeY, int numPoints, int thickSize, double *polyX, double *polyY, char *result);
//int runAlg(int index, int sizeX, int sizeY, int numPoints, double *polyX, double *polyY, char *result);
//int runAlg(int index, int sizeX, int sizeY, int numPoints, double *polyX, double *polyY, char *input, char *result);

int nays8(IplImage *im, int r, int c);
int connectivity(IplImage *im, int r, int c);
void deleteCB(IplImage *im, IplImage *tmp);
void stair(IplImage *im, IplImage *tmp, int dir);
void skeletonize(IplImage *src); 
