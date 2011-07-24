#include <stdio.h>
#include <OpenCV/cv.h>
#include <OpenCV/highgui.h>

int runAlg(int sizeX, int sizeY, char *input, char *result);
int nays8(IplImage *im, int r, int c);
int connectivity(IplImage *im, int r, int c);
void deleteCB(IplImage *im, IplImage *tmp);
void stair(IplImage *im, IplImage *tmp, int dir);
void skeletonize(IplImage *src); 
