#include "TriangleTri.h"

//TriangleTri constructor
TriangleTri::TriangleTri(double * outCloud, int * segCloud, int numPoints, int numSegments, char * ts){
	triangulateio itio; //the input triangulateio struct

	//initialize the required variables in the input triangulateio struct
	itio.numberofpoints = numPoints;
	itio.numberofsegments = numSegments;
	itio.pointlist = (REAL *)malloc(itio.numberofpoints * 2 * sizeof(REAL));
	itio.pointmarkerlist = (int *)malloc(itio.numberofpoints * sizeof(int));
	itio.segmentlist = (int *)malloc(itio.numberofsegments * 2 * sizeof(int));

	//add the vertices and check for negative values

	itio.pointlist = outCloud;
	itio.segmentlist = segCloud;
	for (int i = 0; i < itio.numberofpoints; i++){
		itio.pointmarkerlist[i] = i + 2;
	}
	
	//other required initializations
	itio.numberofpointattributes = 0;
	itio.segmentmarkerlist = NULL;
	itio.numberofholes = 0;
	itio.numberofregions = 0;

	//initialize the required variables in the output triangulateio object
	otio.pointlist = NULL;
	otio.pointmarkerlist = NULL;
	otio.trianglelist = NULL;
	otio.segmentlist = NULL;
	otio.segmentmarkerlist = NULL;

	//call triangulate()
	triangulate(ts, &itio, &otio, NULL);
}

//returns the output triangulateio struct
triangulateio* TriangleTri::getOutput(){
	return &otio;
}