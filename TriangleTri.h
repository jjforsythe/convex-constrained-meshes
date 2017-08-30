#ifndef TRIANGLETRI_H	
#define TRIANGLETRI_H	

#include <iostream>

using namespace std;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;

extern "C" {
#include "triangle.h"
}

class TriangleTri{
public:
	triangulateio otio;

	TriangleTri(double * outCloud, int * segCloud, int numPoints, int numSegments, char * ts);
	triangulateio* getOutput();
};

#endif