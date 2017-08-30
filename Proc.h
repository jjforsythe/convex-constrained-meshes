#ifndef PROC
#define PROC

#include <iostream>
using namespace std;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
using namespace CGAL;

#include <boost/optional/optional_io.hpp>

#include "Mesh.h"

typedef Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 P2;

class Proc{
public:
	Proc(double * outCloud, int n, Mat& image);

	void Line_Segment_Intersections();
	void Line_Line_Intersections();
	void Delete_Short_Segments();
	void Delete_Hanging_Vertices();

	void Delete_Isolated_Short_Segments();

	void Get_Disjoint_Segments();
	void Erase_Segment(int x);

	void Add_Vertex(double x, double y);

	double * PointsToDoubleArray();
	int * SegmentsToIntArray();

	bool OutsideBoundary(Segment_2<K>& seg);

	void shorten_segments_from_boundary();

	int getNumberOfPoints();

	int getNumberOfSegments();

	bool checkAngle(double angle);
	double getMinAngle(double angle);
	
	bool shorten_segment_ii(int a, int b, Point_2<K> pt, double angle, double reduction);
	bool shorten_segment_v1_ii(int a, int b, Point_2<K> pt, double angle, double reduction);

	void deal_with_floating_point_errors();
	Mat DrawSegments(int rows, int cols);
	Mat Proc::DrawSegmentsOnImage(Mat& imgX);

	Segment_2<K> getSegment(int i);
	Line_2<K> Proc::getLine(int i);

	Mat * img;
	double * points;
	vector<int> segments;
	vector<double> pts;
	int seg_number;
	vector<bool> fix_v;
	vector<double> extra_pts;

	int num_seg_pts;

	int getNumSegPts();
};

bool DecreasingInts(int const& p1, int const& p2);
bool IncreasingInts(int const& p1, int const& p2);

double minimum(double x, double y);
double modulusX(double x);


#endif