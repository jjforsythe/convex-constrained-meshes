#ifndef MESH_H	
#define MESH_H

//Define our parameters for epsilon and the minimum angle, and the scale factor for the output
//const double MIN_ANGLE = 0.1745329251; //10 degrees
//const double MIN_ANGLE = 0.2094395102; //12 degrees
//const double MIN_ANGLE = 0.2443460953; //14 degrees
//const double MIN_ANGLE = 0.2792526803; //16 degrees
//const double MIN_ANGLE = 0.3141592653; //18 degrees
const double MIN_ANGLE = 0.3490658504; //20 degrees
//extern double MIN_ANGLE;
//const double MIN_ANGLE = 1.047197551; //60 degrees
const double EPSILON = 3.0;
const double EPS = 0.0000000001;
const double DELTA = 0.01;
const int SCALE_FACTOR = 10;

#include<iostream>
using namespace std;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/intersections.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
using namespace CGAL;

typedef Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 P2;	

//typedef Exact_predicates_tag Itag;

extern "C" {
#include "triangle.h"
}

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
// ----------------------------------------------------------------------------

struct MyTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3f Normal;

	//couldn't get it to work with this double precision when calculating edge lengths
	//typedef OpenMesh::Vec3d Point;
	VertexTraits
	{
	public:
		//holds the number of constrained edges attached to a vertex
		int constedges = 0;
		bool fixed = false;
		Point2d point;
	};
	EdgeTraits
	{
	public:
		//identifies whether an edge is a subsegment of the input lines
		bool constrained = false;
		int segment_idx = -1;
		double length = -1;
		int index = -1;
	};
	FaceTraits
	{
		//index for the Polygon_Create
		int polyidx = -1;
	};
	VertexAttributes(OpenMesh::Attributes::Status);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
};

typedef OpenMesh::PolyMesh_ArrayKernelT<MyTraits>  MyMesh;
typedef MyMesh::Point MP;
typedef MyMesh::VertexHandle MVH;
typedef MyMesh::EdgeHandle MEH;
typedef MyMesh::HalfedgeHandle MHH;
typedef MyMesh::FaceHandle MFH;

const CvScalar Blue = CV_RGB(0, 0, 255);
const CvScalar Red = CV_RGB(255, 0, 0);
const CvScalar Green = CV_RGB(0, 255, 0);
const CvScalar Black = CV_RGB(0, 0, 0);

static int collapses = 0;

class Mesh {
public:


	static double f1s;
	static double f2s;

	MyMesh mesh;
	String filename;
	Mat iMesh;

	double faces_after_triangle, faces_after_merge, faces_after_modify;
	double f1, f2;
	void DrawMesh(String tag);
	void DrawMeshOnWhite(String tag);

	ofstream outfile;
	string outFolder;

	Mesh(triangulateio * tio, Mat image, String f, int n, string outputFolder, int pt_cut);

	string getOutfolder();
	
	void RefineMesh();

	void adjustTJuncs2();
	void resolveCloseEndpoints2();

	vector<MEH> Sort_Edges(vector<MEH>& edges);

	bool isCollapseOK(MHH hh);
	bool isCollapseToPointOK(MHH hh, OpenMesh::Vec3f);
	bool isVertexPositionValid(MVH vh, MHH hh);
	bool Mesh::isVertexPositionValid(MVH vh);	

	bool isFaceConvex(MFH mfh);
	
	bool checkForSegment(MVH vh, int seg_idx);

	bool isVertexFaceConvex(MVH mvh);
	
	int Mesh::Check_Every_Angle_Int();


	void printMeshData();
	Mat scaleImage(Mat I, int scale);

	bool Mesh::Check_Simple_Polygons();

};

//IndexValue class for merging faces
class IndexValue
{
public:
	int index;
	double value;
	IndexValue() {}
	IndexValue(int i, double v) { index = i; value = v; }
};

double Det(Point2d v1, Point2d v2);
void Merge_Faces(MyMesh& pmesh, double delta, double eps);

bool Find_Angle(MyMesh const& pmesh, MHH hh, double eps, double& angle);
bool Find_Angle(Point2d vec1, Point2d vec2, double eps, double& angle);
bool Find_Angle2(MyMesh const& pmesh, MHH hh, double eps, double& angle);

bool Increasing_Unconstrained(MyMesh& pmesh, double max_length,
	vector<MEH>& sorted_edges);

bool DecreasingValues(IndexValue const& p1, IndexValue const& p2);
bool IncreasingValues(IndexValue const& p1, IndexValue const& p2);

bool Merge_Faces_at_Edge(MyMesh& pmesh, MHH hh0, bool& convex0,
	bool& convex1, double delta, double eps);

void Modify_Faces(MyMesh& pmesh, double delta, double eps);

bool Triangular_Cut_Possible(MyMesh& pmesh, MHH hh, double delta, double eps);

bool Remove_Deg2Vertices(MyMesh& pmesh, double delta, double eps);

void addFaceEdges(MyMesh& mesh, MFH fh, vector<MEH>& sorted_edges);


#endif