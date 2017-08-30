#include <iostream>

#include <Windows.h>
using namespace std;

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace cv;

extern "C" {
#include "lsd.h"
}

#include "Proc.h"
#include "TriangleTri.h"
#include "Mesh.h"

//switch for running at the command line - or not
bool cl = true;
string filename;

//Q - no detail for triangle, V - some detailed info (add more Vs for more info)
char triSwitches[] = "pq20czQ";
string inputFolder = "Grayscale";
const int MAX_TESTS = 500;

double* ScanImageToDoubleArray(Mat& I);
vector<string> get_all_files_names_within_folder(string folder);

int main(int argc, char *argv[]){
	Mat image;	
	
	double times[6];

	string outputFolder = "Output/" + inputFolder + "/";

	//vector<string> v = get_all_files_names_within_folder("Input/" + inputFolder);
	
	double t = (double)getTickCount();

	//for (int i = 0; i < v.size() && i < MAX_TESTS; i++) {

		//filename = v[i];
		//cout << endl << "Image " << i + 1 << ": " << filename << endl;
		cout << endl << "Image: " << filename << endl;
		//Create opencv Mat to hold our image and read from file in grayscale format
		if (cl) {
			string f = argv[1];
			filename = "Input/" + f;
			image = imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
		}
		else {
			image = imread("Input/" + inputFolder + "/" + filename, CV_LOAD_IMAGE_GRAYSCALE);   // Read the file
		}

		if (!image.data)
		{
			cout << "No image data \n";
			std::system("pause");
			return -1;
		}

		//convert image from Mat type to double array - necessary input for LSD
		double * dImage = ScanImageToDoubleArray(image);

		//call lsd, dImage now holds our image information
		int n;
		double * out = lsd(&n, dImage, image.cols, image.rows);

		//t = ((double)getTickCount() - t) / getTickFrequency();
		//cout << "LSD time: " << t << endl;

		cout << "Number of lines detected by LSD: " << n << endl;

		//create a Proc instance - class to process the LSD output
		Proc p(out, n, image);

		double * pts = p.PointsToDoubleArray();
		int * segs = p.SegmentsToIntArray();

		int num_points = p.getNumberOfPoints();
		int num_segments = p.getNumberOfSegments();

		string s = filename;
		s.erase(s.find_last_of("."), string::npos);

		imwrite(outputFolder + s + "_LSD_refined" + "_" + to_string(p.seg_number) + ".png", p.DrawSegments(image.rows, image.cols));

		//cout << "No. of refined points: " << num_points << endl;
		//cout << "No. of refined segments: " << num_segments << endl;

		//Create a TriangleTri object which calls Triangle and stores the raw output
		TriangleTri tt(pts, segs, num_points, num_segments, triSwitches);

		//create a Mesh object from the raw output which will create, refine and draw/save the mesh
		Mesh m(tt.getOutput(), image, s, n, outputFolder, p.getNumSegPts());
	//}

	t = ((double)getTickCount() - t) / getTickFrequency();
	cout << "Total time: " << t << endl;

	//double rf1 = Mesh::f1s / v.size();
	//double rf2 = Mesh::f2s / v.size();

	//cout << "RF 1: " << rf1 << endl;
	//cout << "RF 2: " << rf2 << endl;

	imshow("Display window", image);
	waitKey(0);

	return 0;
}

vector<string> get_all_files_names_within_folder(string folder)
{
	vector<string> names;
	char search_path[200];
	sprintf(search_path, "%s/*.jpg", folder.c_str());
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(search_path, &fd);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			// read all (real) files in current folder
			// , delete '!' read other 2 default folder . and ..
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
				names.push_back(fd.cFileName);
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}

//function which scans over the image and ouputs the image as a double array
double* ScanImageToDoubleArray(Mat& I)
{
	// accept only char type matrices
	CV_Assert(I.depth() != sizeof(uchar));

	int channels = I.channels();

	int nRows = I.rows;
	int nCols = I.cols * channels;

	if (I.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;
	}

	double * dI = (double *)malloc(I.rows * I.cols * sizeof(double));

	int i, j;
	uchar* p;
	for (i = 0; i < nRows; ++i)
	{
		p = I.ptr<uchar>(i);
		for (j = 0; j < nCols; ++j)
		{
			dI[j + i*I.cols] = p[j];
		}
	}

	return dI;
}