#include "Proc.h"

/***CONSTRUCTOR***/

Proc::Proc(double * outCloud, int n, Mat& image){

	img = &image;

	//add the vertices into a double array where the segments are defined implicitly
	points = (double *)malloc((n + 2) * 4 * sizeof(double));

	ofstream outfile;
	outfile << fixed;
	outfile << setprecision(50);
	outfile.open("Output/lsd.txt");

	outfile << n << "\n";

	for (int i = 0; i < n; i++){
		points[4 * i] = outCloud[7 * i];
		points[4 * i + 1] = outCloud[7 * i + 1];
		points[4 * i + 2] = outCloud[7 * i + 2];
		points[4 * i + 3] = outCloud[7 * i + 3];
		//outfile << points[4 * i] << "\t" << points[4 * i + 1] << "\t" << points[4 * i + 2] << "\t" << points[4 * i + 3] << "\n";
	}	

	//add the vertices at the boundary
	points[4 * n] = 0;
	points[4 * n + 1] = 0;
	points[4 * n + 2] = img->cols;
	points[4 * n + 3] = 0;
	points[4 * n + 4] = img->cols;
	points[4 * n + 5] = img->rows;
	points[4 * n + 6] = 0;
	points[4 * n + 7] = img->rows;

	//points into the pts vector
	for (int i = 0; i < (n + 2) * 4 ; i+=2){
		Add_Vertex(points[i], points[i+1]);
	}

	//initialize segments - a vector to hold the vertex indices for the segments
	for (int i = 0; i < n; i++){
		segments.push_back(2 * i);
		segments.push_back(2 * i + 1);
	}
	//segments for boundary
	segments.push_back(2 * n);
	segments.push_back(2 * n + 1);
	segments.push_back(2 * n + 1);
	segments.push_back(2 * n + 2);
	segments.push_back(2 * n + 2);
	segments.push_back(2 * n + 3);
	segments.push_back(2 * n + 3);
	segments.push_back(2 * n);

	//holds total number of segments
	seg_number = segments.size() / 2;

	string of = "Output/";

	//output the LSD line segments
	//imwrite(of + "101070/101070" + "_LSD_output" + "_" + to_string(seg_number) + ".png", DrawSegments(image.rows, image.cols));
	//imwrite(of + "cameraman512/cameraman512" + "_LSD_outputX" + "_" + to_string(seg_number) + ".png", DrawSegmentsOnImage(image));

	shorten_segments_from_boundary();

	Line_Segment_Intersections();
	//imwrite(of + "24063/24063" + "_LSD_output1" + "_" + to_string(seg_number) + ".png", DrawSegments(image.rows, image.cols));
	Line_Line_Intersections();

	//cout << "#pts: " << pts.size() / 2 << endl;
	//imwrite(of + "24063/24063" + "_LSD_output2" + "_" + to_string(seg_number) + ".png", DrawSegments(image.rows, image.cols));
	Get_Disjoint_Segments();

	//cout << "#pts: " << pts.size() / 2 << endl;

	num_seg_pts = pts.size() / 2;
	for (int i = 0; i < 2 * num_seg_pts; i += 2) {
		outfile << pts[i] << "\t" << pts[i + 1] << "\n";
	}

	Delete_Isolated_Short_Segments();
	Delete_Short_Segments();
	Delete_Hanging_Vertices();
	deal_with_floating_point_errors();

	//cout << "Number of points before we get: " << pts.size() / 2 << endl;
	num_seg_pts = pts.size() / 2;

	for (int i = 0; i < 2*num_seg_pts; i+=2) {
		outfile << pts[i] << "\t" << pts[i + 1] << "\n";
 	}

	//cout << "#extra pts: " << extra_pts.size() << endl;

	for (int i = 0; i < extra_pts.size(); i++){
		pts.push_back(extra_pts[i]);
		outfile << extra_pts[i] << "\t";
	}

	outfile.close();

	seg_number = segments.size() / 2;
	cout << "Number of refined LSD lines: " << seg_number << endl;
}

/***PRIMARY FUNCTIONS***/

void Proc::shorten_segments_from_boundary() {

	double eps_cls = (EPSILON - 0.5) * (EPSILON - 0.5);
	bool deleteI;

	Segment_2<K> seg;
	Line_2<K> line;

	vector<Segment_2<K>> b_segs;
	vector<Line_2<K>> b_lines;
	double angle;
	Vector_2<K> v1, v2;

	b_segs.push_back(getSegment((segments.size() / 2) - 4));
	b_segs.push_back(getSegment((segments.size() / 2) - 3));
	b_segs.push_back(getSegment((segments.size() / 2) - 2));
	b_segs.push_back(getSegment((segments.size() / 2) - 1));

	b_lines.push_back(getLine((segments.size() / 2) - 4));
	b_lines.push_back(getLine((segments.size() / 2) - 3));
	b_lines.push_back(getLine((segments.size() / 2) - 2));
	b_lines.push_back(getLine((segments.size() / 2) - 1));
	
	for (int j = 0; j < 4; j++) {

		for (int i = 0; i < (segments.size() / 2) - 4; i++) {

			deleteI = false;

			seg = getSegment(i);
			line = getLine(i);

			//compare segment with each boundary segment



			if (squared_distance(b_segs[j], seg) <= eps_cls) {
				//if closer than eps - 0.5 we should find which vertex we shorten from

				auto pintersect = intersection(b_lines[j], line);

				if (pintersect) {
					if (const Line_2<K>* s = boost::get<Line_2<K>>(&*pintersect)) {
						// handle segment 
						// intersecting as a segment corresponds to being on the boundary, we do not consider input segments which could be colinear
					}
					else if (const P2 * ppp = boost::get<Point_2<K> >(&*pintersect)) {

						if (squared_distance(b_segs[j], seg.source()) < squared_distance(b_segs[j], seg.target())) {
							//seg's source is closer

							v1 = seg.source() - *ppp;
							v2 = b_segs[j].source() - *ppp;

							Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

							double s = sin(getMinAngle(angle));
							double reduction = (EPSILON - 0.5) / s;
							double reduction_sub = reduction - sqrt(v1.squared_length());

							//cout << "reduction: " << reduction << ", is finite: " << isfinite(reduction) << endl;
							//cout << "boundary seg: " << b_segs[j] << ", seg: " << seg << endl;

							if (isfinite(reduction)) {
								if (!shorten_segment_v1_ii(i, (segments.size() / 2) - (4 - j), *ppp, angle, reduction_sub)) {
									deleteI = true;
								}
							}
						}
						else {
							//seg's target is closer

							v1 = seg.target() - *ppp;
							v2 = b_segs[j].source() - *ppp;

							Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

							double s = sin(getMinAngle(angle));
							double reduction = (EPSILON - 0.5) / s;
							double reduction_sub = reduction - sqrt(v1.squared_length());

							//cout << "reduction: " << reduction << ", is finite: " << isfinite(reduction) << endl;
							//cout << "boundary seg: " << b_segs[j] << ", seg: " << seg << endl;

							if (isfinite(reduction)) {
								if (!shorten_segment_ii(i, (segments.size() / 2) - (4 - j), *ppp, angle, reduction_sub)) {
									deleteI = true;
								}
							}
						}

					}
				}
			}


			if (deleteI) {
				i--;
			}
		}

	}
}

//to find intersections on segments to be called before we make all segments disjoint
void Proc::Line_Segment_Intersections() {

	cout << "Performing line-segment intersections..." << endl;

	double eps_sq = EPSILON * EPSILON;
	vector<IndexValue> seg_listA, seg_listB;
	map<int, P2> vertices;
	vector<int> segs_to_add;
	Line_2<K> line;
	Segment_2<K> seg, seg1;
	Vector_2<K> v1, v2;
	double angle;
	bool Agood, Bgood;

	//iterate over all segments to compute intersections between its line and all other segments
	//we look to extend from both endpoints of every segment
	for (int i = 0; i < segments.size() / 2; i++) {

		//if 
		
		//set the current lsd segment
		line = getLine(i);
		seg1 = getSegment(i);

		Agood = true;
		Bgood = true;

		for (int j = 0; j < segments.size() / 2; j++) {
			if (j == i) continue;
			seg = getSegment(j);

			//if a vertex is already attached to 2 segments we will not extend from it
			
			//cout << seg1[0] << endl << seg[0] << endl << seg1[0] << endl << seg[1] << endl;

			//|| fix_v[segments[2 * i]]
			//|| fix_v[segments[2 * i + 1]]

			//workaround trick for computing whether points are fixed
			if ((seg1[0] == seg[0]) || (seg1[0] == seg[1])) {
				Agood = false;
			}

			if ((seg1[1] == seg[0]) || (seg1[1] == seg[1]) ) {
				Bgood = false;
			}

			// a work-around for the intersection function specifically in 23084.jpg at point 226.875, 225.113 where an incorrect intersection was found
			if (do_intersect(line, seg)) {

				auto pintersect = intersection(line, seg);

				if (pintersect) {
					if (const Segment_2<K>* s = boost::get<Segment_2<K>>(&*pintersect)) {
						// handle segment 
						// intersecting as a segment corresponds to being on the boundary, we do not consider input segments which could be colinear
					}
					else if (const P2 * ppp = boost::get<Point_2<K> >(&*pintersect)) {

						if ((seg[0] != seg1[0]) && (seg[0] != seg1[1]) && (seg[1] != seg1[0]) && (seg[1] != seg1[1])) {

							if ((!do_intersect(seg1, seg))) {

								Vector_2<K> vec = seg1.source() - *ppp;
								Vector_2<K> vecX = seg1.target() - *ppp;

								//if (vec.squared_length() < eps_sq && vec.squared_length() < vecX.squared_length()){
								if (vec.squared_length() < vecX.squared_length()) {
									//add to the lists
									seg_listA.push_back(IndexValue(j, vec.squared_length()));
									vertices.insert(make_pair(j, *ppp));
								}
								else {
									vec = seg1.target() - *ppp;

									//if (vec.squared_length() < eps_sq){
									seg_listB.push_back(IndexValue(j, vec.squared_length()));
									vertices.insert(make_pair(j, *ppp));
									//}
								}

							}
						}
						//shorten segments which are closer than min_distance here

						/*
						else if ((seg1[0] == seg[0]) || (seg1[0] == seg[1])){
						Agood = false;
						}
						else if ((seg1[1] == seg[0]) || (seg1[1] == seg[1])){
						Bgood = false;
						}
						*/
					}
				}
			}
		}

		//sort the lists here
		//A is for moving the first vertex, B is for moving the second
		sort(seg_listA.begin(), seg_listA.end(), IncreasingValues);
		sort(seg_listB.begin(), seg_listB.end(), IncreasingValues);

		//v1 is a vector along the line
		v1 = (seg1.source() - seg1.target());
		bool intersectedA = false;
		bool intersectedB = false;

		int indexA, indexB;

		bool deleteI = false;

		if (seg_listA.size() > 0 && Agood) {
			seg = getSegment(seg_listA[0].index);

			//code to check the angles
			v2 = (seg.point(0) - seg.point(1));

			Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

			indexA = seg_listA[0].index;
			P2 pt = vertices[indexA];

			//cout << vertices[indexA] << endl;

			Vector_2<K> vv1, vv2, vv3;

			vv1 = seg1.source() - pt;
			vv2 = seg.source() - pt;
			vv3 = seg.target() - pt;

			double vv2_angle, vv3_angle;
			Find_Angle(Point2d(vv1[0], vv1[1]), Point2d(vv2[0], vv2[1]), 0.0000000001, vv2_angle);
			Find_Angle(Point2d(vv1[0], vv1[1]), Point2d(vv3[0], vv3[1]), 0.0000000001, vv3_angle);

			if ((vv1.squared_length() < eps_sq) && checkAngle(angle)) {
				//cout << "A" << endl;
				pts[2 * segments[2 * i]] = pt.x();
				pts[2 * segments[2 * i] + 1] = pt.y();
				fix_v[segments[2 * i]] = true;

				segments.push_back(segments[2 * i]);
				segments.push_back(segments[2 * indexA + 1]);

				segments[2 * indexA + 1] = segments[2 * i];

				//cout << "seg i: " << getSegment(i) << ", segA: " << getSegment(indexA) << ", segC: " << getSegment((segments.size()/2)-1) << endl;

				intersectedA = true;
			}
			else if ((vv1.squared_length() < eps_sq) && (vv2.squared_length() < eps_sq) && ((minimum(vv2_angle, (2 * M_PI) - vv2_angle) < MIN_ANGLE))) {
				//cout << "B" << endl;
				pts[2 * segments[2 * i]] = pt.x();
				pts[2 * segments[2 * i] + 1] = pt.y();
				fix_v[segments[2 * i]] = true;

				segments[2 * indexA] = segments[2 * i];

				intersectedA = true;
			}
			else if ((vv1.squared_length() < eps_sq) && (vv3.squared_length() < eps_sq) && ((minimum(vv3_angle, (2 * M_PI) - vv3_angle) < MIN_ANGLE))) {
				//cout << "C" << endl;
				pts[2 * segments[2 * i]] = pt.x();
				pts[2 * segments[2 * i] + 1] = pt.y();
				fix_v[segments[2 * i]] = true;

				segments[2 * indexA + 1] = segments[2 * i];

				intersectedA = true;
			}
			else if (!checkAngle(angle) && (squared_distance(seg1.source(), seg) < eps_sq)) {
				//cout << "D" << endl;
				//if we cannot make an intersection we shorten the segment to epsilon away from the line

				//calculate the length to reduce segments by at the intersection

				//cout << "seg i: " << getSegment(i) << ", seg j: " << getSegment(seg_listA[0].index) << endl;
				double s = sin(getMinAngle(angle));
				double reduction = EPSILON / s;
				double reduction_sub = reduction - sqrt(vv1.squared_length());
				//note one of these subs creates a small unconstrained angle in 254033.jpg at 99.1, 200.8

				//cout << "SHORTENING SEGMENT FROM POINT: " << pt << "\t" << "REDUCTION: " << reduction_sub << endl;
				//cout << "SEGMENT: " << seg1.target() << "\t" << seg1.source() << endl;
				if (isfinite(reduction)) {
					if (!shorten_segment_v1_ii(i, seg_listA[0].index, pt, getMinAngle(angle), reduction_sub))
						deleteI = true;
				}
			}
		}
		if (seg_listB.size() > 0 && Bgood && !deleteI) {
			seg = getSegment(seg_listB[0].index);

			//cout << "#275B: " << getSegment(275) << endl;

			//code to check the angles
			v2 = (seg.point(0) - seg.point(1));

			Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

			indexB = seg_listB[0].index;
			P2 pt = vertices[indexB];

			//cout << vertices[indexB] << endl;

			Vector_2<K> vv1, vv2, vv3;

			vv1 = seg1.target() - pt;
			vv2 = seg.source() - pt;
			vv3 = seg.target() - pt;

			double vv2_angle, vv3_angle;
			Find_Angle(Point2d(vv1[0], vv1[1]), Point2d(vv2[0], vv2[1]), 0.0000000001, vv2_angle);
			Find_Angle(Point2d(vv1[0], vv1[1]), Point2d(vv3[0], vv3[1]), 0.0000000001, vv3_angle);

			if ((vv1.squared_length() < eps_sq) && checkAngle(angle)) {

				pts[2 * segments[2 * i + 1]] = pt.x();
				pts[2 * segments[2 * i + 1] + 1] = pt.y();
				fix_v[segments[2 * i + 1]] = true;

				segments.push_back(segments[2 * i + 1]);
				segments.push_back(segments[2 * indexB + 1]);

				segments[2 * indexB + 1] = segments[2 * i + 1];
				intersectedB = true;

			}
			else if ((vv1.squared_length() < eps_sq) && (vv2.squared_length() < eps_sq) && ((minimum(vv2_angle, (2 * M_PI) - vv2_angle) < MIN_ANGLE))) {
				pts[2 * segments[2 * i + 1]] = pt.x();
				pts[2 * segments[2 * i + 1] + 1] = pt.y();
				fix_v[segments[2 * i + 1]] = true;

				segments[2 * indexB] = segments[2 * i + 1];

				intersectedB = true;
			}
			else if ((vv1.squared_length() < eps_sq) && (vv3.squared_length() < eps_sq) && ((minimum(vv3_angle, (2 * M_PI) - vv3_angle) < MIN_ANGLE))) {
				pts[2 * segments[2 * i + 1]] = pt.x();
				pts[2 * segments[2 * i + 1] + 1] = pt.y();
				fix_v[segments[2 * i + 1]] = true;

				segments[2 * indexB + 1] = segments[2 * i + 1];

				intersectedB = true;
			}
			else if (!checkAngle(angle) && (squared_distance(seg1.target(), seg) < eps_sq)) {
				//if we cannot make an intersection we shorten the segment to epsilon away from the line

				double s = sin(getMinAngle(angle));
				double reduction = EPSILON / s;
				double reduction_sub = reduction - sqrt(vv1.squared_length());

				//cout << "SHORTENING SEGMENT FROM POINT: " << pt << "\t" << "REDUCTION: " << reduction_sub <<  endl;
				//cout << "SEGMENT: " << seg1.target() << "\t" << seg1.source() << endl;
				if (isfinite(reduction)) {
					if (!shorten_segment_ii(i, seg_listB[0].index, pt, getMinAngle(angle), reduction_sub))
						deleteI = true;
				}
			}

			//cout << "#275Be: " << getSegment(275) << endl;
		}

		if (deleteI) {
			i--;
		}

		seg_listA.clear();
		seg_listB.clear();
		vertices.clear();
	}

}

void Proc::Line_Line_Intersections() {

	cout << "Performing line-line intersections..." << endl;

	double eps_sq = EPSILON * EPSILON;

	vector<IndexValue> seg_listA, seg_listB;
	map<int, P2> vertices;

	Line_2<K> line1, line2;
	Segment_2<K> seg1, seg2;

	Vector_2<K> v1, v2;
	double angle;

	bool Agood, Bgood;

	for (int i = 0; i < (segments.size() / 2) - 1; i++) {
		//set the current lsd segment
		line1 = getLine(i);
		seg1 = getSegment(i);

		Agood = true;
		Bgood = true;

		for (int j = i + 1; j < segments.size() / 2; j++) {
			if (j == i) continue;
			line2 = getLine(j);
			seg2 = getSegment(j);

			//if a vertex is already attached to 2 segments we will not extend from it
			//cout << seg1[0] << endl << seg[0] << endl << seg1[0] << endl << seg[1] << endl;
			if ((seg1[0] == seg2[0]) || (seg1[0] == seg2[1])) {
				Agood = false;
			}

			if ((seg1[1] == seg2[0]) || (seg1[1] == seg2[1])) {
				Bgood = false;
			}

			if (fix_v[segments[2 * i]]) {
				Agood = false;
			}

			if (fix_v[segments[2 * i + 1]]) {
				Bgood = false;
			}

			auto pintersect = intersection(line1, line2);
			if (pintersect) {
				if (const Line_2<K>* s = boost::get<Line_2<K>>(&*pintersect)) {
					// handle segment 
					// intersecting as a line corresponds to being on the boundary, we do not consider input segments which could be colinear
				}
				else if (const P2 * ppp = boost::get<Point_2<K> >(&*pintersect)) {

					if ((seg1[0] != seg2[0]) && (seg1[0] != seg2[1]) && (seg1[1] != seg2[0]) && (seg1[1] != seg2[1])) {
						//impose the stronger condition that the line does not intersect with the second segment
						if ((!do_intersect(line1, seg2)) && !do_intersect(line2, seg1)) {

							Vector_2<K> vec = seg1.source() - *ppp;
							Vector_2<K> vec2 = seg2.source() - *ppp;
							Vector_2<K> vec3 = seg2.target() - *ppp;

							if (((vec2.squared_length() < eps_sq) && !fix_v[segments[2 * j]]) || ((vec3.squared_length() < eps_sq) && !fix_v[segments[2 * j + 1]])) {
								if (vec.squared_length() < eps_sq) {
									//add to the lists
									seg_listA.push_back(IndexValue(j, vec.squared_length()));
									vertices.insert(make_pair(j, *ppp));
								}
								else {
									vec = seg1.target() - *ppp;

									if (vec.squared_length() < eps_sq) {
										seg_listB.push_back(IndexValue(j, vec.squared_length()));
										vertices.insert(make_pair(j, *ppp));
									}
								}
							}
							//shorten segments which are closer than min_distance here...

						}
					}
					else if ((seg1[0] == seg2[0]) || (seg1[0] == seg2[1])) {
						Agood = false;
					}
					else if ((seg1[1] == seg2[0]) || (seg1[1] == seg2[1])) {
						Bgood = false;
					}
				}
			}
		}

		//sort the lists here
		//A is for moving the first vertex, B is for moving the second
		sort(seg_listA.begin(), seg_listA.end(), IncreasingValues);
		sort(seg_listB.begin(), seg_listB.end(), IncreasingValues);

		//this will still only look at the closest other endpoint but there may be the possiblility of others if the closest generates a small angle
		if (seg_listA.size() > 0 && Agood) {

			int index = seg_listA[0].index;
			P2 pt = vertices[index];

			seg2 = getSegment(seg_listA[0].index);

			//cout << "segments: " << seg1 << ", " << seg2 << endl;

			//this should generate the angle between the lines, and avoids missing big angle intersections
			v1 = (seg1.source() - pt);
			v2 = (seg2.source() - pt);

			Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

			Vector_2<K> vv1 = seg2.source() - pt;
			Vector_2<K> vv2 = seg2.target() - pt;

			Vector_2<K> vv3 = seg1.source() - pt;
			Vector_2<K> vv4 = seg1.target() - pt;

			//cout << pt << "\t" << angle << "\t" << getMinAngle(angle) << endl;

			if (minimum(angle, (2 * M_PI) - angle) > MIN_ANGLE) {
				//if (checkAngle(angle)) {
				pts[2 * segments[2 * i]] = pt.x();
				pts[2 * segments[2 * i] + 1] = pt.y();
				fix_v[segments[2 * i]] = true;

				//seg2's source is closer to the intersection point
				if (vv1.squared_length() <= vv2.squared_length()) {
					segments[2 * index] = segments[2 * i];
				}
				//else seg2's target is closer to the intersection point
				else {
					segments[2 * index + 1] = segments[2 * i];
				}

				//cout << "good intersection found at x: " << pt.x() << ", y: " << pt.y() << endl;
			}
			else {

				//using shorten_segment_v1_ii() here because shortening away from segment i's source

				//seg2's source is closer to the intersection point
				if (vv1.squared_length() <= vv2.squared_length()) {
					//shorten away from the source

					//cout << "I" << endl;
					double s = sin(getMinAngle(angle));
					double reduction = EPSILON / s;
					double reduction_sub = reduction - sqrt(vv3.squared_length());
					//double rss = reduction_sub * reduction_sub;
					shorten_segment_v1_ii(i, seg_listA[0].index, pt, getMinAngle(angle), reduction_sub);
				}
				//else seg2's target is closer to the intersection point
				else {
					//shorten away from the target

					//cout << "II" << endl;
					double s = sin(getMinAngle(angle));
					double reduction = EPSILON / s;
					double reduction_sub = reduction - sqrt(vv3.squared_length());
					//double rss = reduction_sub * reduction_sub;

					shorten_segment_v1_ii(i, seg_listA[0].index, pt, getMinAngle(angle), reduction_sub);
				}
			}
		}

		if (seg_listB.size() > 0 && Bgood) {

			int index = seg_listB[0].index;
			P2 pt = vertices[index];

			seg2 = getSegment(seg_listB[0].index);

			//cout << "segments: " << seg1 << ", " << seg2 << endl;

			//this should generate the angle between the lines, and avoids missing big angle intersections
			//should this be seg1.target()?
			v1 = (seg1.source() - pt);
			v2 = (seg2.source() - pt);

			Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

			Vector_2<K> vv1 = seg2.source() - pt;
			Vector_2<K> vv2 = seg2.target() - pt;

			Vector_2<K> vv3 = seg1.source() - pt;
			Vector_2<K> vv4 = seg1.target() - pt;
			//cout << pt << "\t" << angle << "\t" << getMinAngle(angle) << endl;

			if (minimum(angle, (2 * M_PI) - angle) > MIN_ANGLE) {
				//if (checkAngle(angle)) {
				pts[2 * segments[2 * i + 1]] = pt.x();
				pts[2 * segments[2 * i + 1] + 1] = pt.y();
				fix_v[segments[2 * i + 1]] = true;

				//seg2's source is closer to the intersection point
				if (vv1.squared_length() <= vv2.squared_length()) {
					segments[2 * index] = segments[2 * i + 1];
				}
				//else seg2's target is closer to the intersection point
				else {
					segments[2 * index + 1] = segments[2 * i + 1];
				}

				//cout << "good intersection found at x: " << pt.x() << ", y: " << pt.y() << endl;
			}
			else {

				//using shorten_segment_ii() here because shortening away from segment i's target

				//seg2's source is closer to the intersection point
				if (vv1.squared_length() <= vv2.squared_length()) {
					//shorten away from the source
					//cout << "IB" << endl;
					double s = sin(getMinAngle(angle));
					double reduction = EPSILON / s;
					double reduction_sub = reduction - sqrt(vv4.squared_length());
					//double rss = reduction_sub * reduction_sub;

					//cout << "seg1: " << seg1 << ", seg2: " << seg2 << endl;

					shorten_segment_ii(i, seg_listB[0].index, pt, getMinAngle(angle), reduction_sub);
				}
				//else seg2's target is closer to the intersection point
				else {
					//shorten away from the target
					//cout << "IIB" << endl;
					double s = sin(getMinAngle(angle));
					double reduction = EPSILON / s;
					double reduction_sub = reduction - sqrt(vv4.squared_length());
					//double rss = reduction_sub * reduction_sub;

					shorten_segment_ii(i, seg_listB[0].index, pt, getMinAngle(angle), reduction_sub);
				}
			}


		}

		seg_listA.clear();
		seg_listB.clear();
		vertices.clear();

	}
}

//Main function to create the PL complex (of vertices and edges)
void Proc::Get_Disjoint_Segments() {

	cout << "Getting disjoint segments..." << endl;

	int init_size = segments.size() / 2;
	int idx;

	Segment_2<K> seg_1, seg_2;

	//iterate over all segment intersections (takes N^2/2 time)
	for (int i = 0; i <(segments.size() / 2) - 1; i++) {
		//set the current lsd segment
		seg_1 = getSegment(i);

		for (int j = i + 1; j < (segments.size() / 2); j++) {
			seg_2 = getSegment(j);

			//get the intersection point of the two segments
			auto pintersect = intersection(seg_1, seg_2);

			//if the segments intersect (not at their endpoints) split up the intersection into segments
			if (pintersect) {
				if (const Segment_2<K>* s = boost::get<Segment_2<K>>(&*pintersect)) {
					// handle segment - we do not expect two colinear segments to intersect
					cout << "segment: " << seg_1.source() << "\t" << seg_1.target() << "\t" << seg_2.source() << "\t" << seg_2.target() << endl;
				}
				else if (const P2 * ppp = boost::get<Point_2<K> >(&*pintersect)) {

					//check that the intersection is not at the endpoints of either segment
					if ((seg_1[0] != seg_2[0]) && (seg_1[0] != seg_2[1]) && (seg_1[1] != seg_2[0]) && (seg_1[1] != seg_2[1])) {

						//cout << "segments: " << seg_1 << ", " << seg_2 << endl;

						//get vectors along the segments and compute the angle between these vectors
						double angle;
						Vector_2<K> v1 = (seg_1.point(0) - seg_1.point(1));
						Vector_2<K> v2 = (seg_2.point(0) - seg_2.point(1));

						Find_Angle(Point2d(v1[0], v1[1]), Point2d(v2[0], v2[1]), 0.0000000001, angle);

						//add the intersection point to our vector of vertices
						idx = pts.size() / 2;
						pts.push_back(ppp->x());
						pts.push_back(ppp->y());

						int temp_s = segments.size() / 2;

						//add to two new segments to our segment list and edit the two existing ones
						segments.push_back(segments[2 * i + 1]);
						segments.push_back(idx);

						segments.push_back(segments[2 * j + 1]);
						segments.push_back(idx);

						segments[2 * i + 1] = idx;

						segments[2 * j + 1] = idx;

						//deal with the situation we have a small angle at the intersection
						if (!checkAngle(angle)) {

							//get the lengths of the 4 segments at the intersection
							double i1 = getSegment(i).squared_length();
							double i2 = getSegment(temp_s).squared_length();
							double j1 = getSegment(j).squared_length();
							double j2 = getSegment(temp_s + 1).squared_length();

							//calculate the length to reduce segments by at the intersection
							double s = sin(getMinAngle(angle));
							double reduction = EPSILON / s;
							double reduction_sq = reduction * reduction;

							//vectors from the intersection point along each of the 4 segments
							Vector_2<K> vi1 = getSegment(i).source() - *ppp;
							Vector_2<K> vi2 = getSegment(temp_s).source() - *ppp;
							Vector_2<K> vj1 = getSegment(j).source() - *ppp;
							Vector_2<K> vj2 = getSegment(temp_s + 1).source() - *ppp;

							int im, jm;

							//find the longest subsegments of each original segment
							double ang;
							Vector_2<K> max_i, max_j;

							if (i1 > i2) {
								max_i = vi1;
								im = i;
							}

							else {
								max_i = vi2;
								im = temp_s;
							}

							if (j1 > j2) {
								max_j = vj1;
								jm = j;
							}
							else {
								max_j = vj2;
								jm = temp_s + 1;
							}

							//the case where neither segment lies outside the boundary -> neither segment is a boundary segment
							if (!OutsideBoundary(seg_1) && !OutsideBoundary(seg_2)) {

								Find_Angle(Point2d(max_i[0], max_i[1]), Point2d(max_j[0], max_j[1]), EPS, ang);

								//if the small angle is between the longer subsegments we shorten the shorten the shorter original segment otherwise we shorten the two shorter segments
								//some improvement possible if we do have short edges and so don't need to shorten the other by so much
								if (minimum(ang, (2 * M_PI) - ang) < MIN_ANGLE) {

									if ((i1 + i2) > (j1 + j2)) {

										//keep the 'i' segments
										segments[2 * i + 1] = segments[2 * temp_s + 1];

										//reduce the length of the second segment
										shorten_segment_ii(temp_s + 1, i, *ppp, minimum(ang, (2 * M_PI) - ang), reduction);

										segments.erase(segments.begin() + 2 * temp_s + 1);
										segments.erase(segments.begin() + 2 * temp_s);

										//reduce the length of j1
										shorten_segment_ii(j, i, *ppp, minimum(ang, (2 * M_PI) - ang), reduction);

									}

									else {

										//keep the 'j' segments
										segments[2 * j + 1] = segments[2 * (temp_s + 1) + 1];

										segments.erase(segments.begin() + 2 * (temp_s + 1) + 1);
										segments.erase(segments.begin() + 2 * (temp_s + 1));

										//reduce the length of the second segment
										shorten_segment_ii(temp_s, j, *ppp, minimum(ang, (2 * M_PI) - ang), reduction);

										//reduce the length of i1
										shorten_segment_ii(i, j, *ppp, minimum(ang, (2 * M_PI) - ang), reduction);

									}
								}
								else {
									if (j1 > j2) {
										//shorten_segment_ii(temp_s + 1, reduction);
										shorten_segment_ii(temp_s + 1, im, *ppp, getMinAngle(angle), reduction);
									}

									if (i1 > i2) {
										//shorten_segment_ii(temp_s, reduction);
										shorten_segment_ii(temp_s, jm, *ppp, getMinAngle(angle), reduction);
									}


									if (j1 < j2) {
										//shorten_segment_ii(j, reduction);
										shorten_segment_ii(j, im, *ppp, getMinAngle(angle), reduction);
									}

									if (i1 < i2) {
										//shorten_segment_ii(i, reduction);
										shorten_segment_ii(i, jm, *ppp, getMinAngle(angle), reduction);
									}
								}
							}
							else {
								//handle boundary terms
								//determine which of the original segments is intersecting the boundary and erase or shorten the subsegments
								//ERRORS HERE
								if (OutsideBoundary(seg_1)) {
									shorten_segment_ii(temp_s, jm, *ppp, getMinAngle(angle), reduction);

									shorten_segment_ii(i, jm, *ppp, getMinAngle(angle), reduction);

								}
								else if (OutsideBoundary(seg_2)) {

									shorten_segment_ii(temp_s + 1, im, *ppp, getMinAngle(angle), reduction);

									shorten_segment_ii(j, im, *ppp, getMinAngle(angle), reduction);

								}

							}
						}

						seg_1 = getSegment(i);

					}
				}
			}
		}
	}
	seg_number = segments.size() / 2;
}

//finds all segments less than epsilon and saves a vertex at the midpoint of each one
void Proc::Delete_Isolated_Short_Segments() {
	Segment_2<K> seg1, seg2;
	double eps_sq = EPSILON * EPSILON;

	for (int i = 0; i < segments.size() / 2; i++) {

		seg1 = getSegment(i);

		if ((seg1.squared_length() < eps_sq)) {
			bool isolated = true;

			for (int j = 0; j < segments.size() / 2; j++) {
				seg2 = getSegment(j);

				//cout << "segments: " << seg1 << ", " << seg2 << endl;
 
				if (seg1 != seg2) {
					if (seg1.source() == seg2.source() || seg1.source() == seg2.target() || seg1.target() == seg2.source() || seg1.target() == seg2.target()) {
						isolated = false;
					}
				}
			}

			if (isolated) {

				//still outputs in v3
				//cout << "SMALL ISOLATED SEGMENT\n";
				//cout << seg1.source() << "\t" << seg1.target() << endl;

				Vector_2<K> v = (Vector_2<K>(seg1.source()[0], seg1.source()[1]) + Vector_2<K>(seg1.target()[0], seg1.target()[1])) / 2;

				extra_pts.push_back(v[0]);
				extra_pts.push_back(v[1]);

			}

		}
	}
}

//deletes all segments shorter than EPSILON and any segments outside the boundary
void Proc::Delete_Short_Segments() {
	Segment_2<K> seg;
	double eps_sq = EPSILON * EPSILON;

	for (int i = 0; i < segments.size() / 2; i++) {

		seg = getSegment(i);

		if (OutsideBoundary(seg)) {
			segments.erase(segments.begin() + 2 * i + 1);
			segments.erase(segments.begin() + 2 * i);
			i--;
		}
		else if ((seg.squared_length() < eps_sq)) {
			Erase_Segment(i);
			i--;
		}
	}
}

//finds vertices which are not being used anymore and deletes them, no solution yet for deleting them as their position in the vector is important
void Proc::Delete_Hanging_Vertices() {
	vector<double> ptscopy = pts;

	for (int i = 0; i < segments.size() / 2; i++) {
		ptscopy[2 * segments[2 * i]] = -1;
		ptscopy[2 * segments[2 * i] + 1] = -1;
		ptscopy[2 * segments[2 * i + 1]] = -1;
		ptscopy[2 * segments[2 * i + 1] + 1] = -1;
	}

	vector<int> indices;
	for (int i = 0; i < ptscopy.size(); i++) {
		if ((ptscopy[i] != -1)) {

			if (i % 2 == 0)
				indices.push_back(i);
		}
	}

	sort(indices.begin(), indices.end(), DecreasingInts);

	for (int i = 0; i < indices.size(); i++) {
		pts.erase(pts.begin() + indices[i] + 1);
		pts.erase(pts.begin() + indices[i]);

		for (int j = 0; j < segments.size(); j++) {
			if (segments[j] >= (indices[i] / 2))
				segments[j]--;
		}
	}

}

/***SEGMENT SHORTENING***/

//takes a length of reduction and a point to reduce a segment away from, and the angle at the intersection
//reduces the segment away from the target()
bool Proc::shorten_segment_ii(int a, int b, Point_2<K> pt, double angle, double reduction){
	Point_2<K> p1, p2;
	int idx2;

	//cout << "reduction1: " << reduction << "\n";

	p1 = Point_2<K>(pts[2 * segments[2 * a]], pts[2 * segments[2 * a] + 1]);
	p2 = Point_2<K>(pts[2 * segments[2 * a + 1]], pts[2 * segments[2 * a + 1] + 1]);
	Vector_2<K> v = p1 - p2;
	//make v the length of the reduction
	v = (v * reduction) / sqrt(v.squared_length());
	//np = 'X'
	Point_2<K> np = p2 + v;
	//cout << "np: " << np << endl;

	//we must check our distance is from the segment not the infinite line
	Segment_2<K> seg2 = getSegment(b);

	double dist_x = sqrt(squared_distance(seg2, np));

	//cout << "dist_x: " << dist_x << endl;

	if (dist_x <= (EPSILON + EPS)){
		;//this is the right point
	}
	else{ 
		Point_2<K> cls_pt;
		if (sqrt(squared_distance(seg2.source(), np)) <= (dist_x + EPS)){
			cls_pt = seg2.source();
		}
		else{
			cls_pt = seg2.target();
		}

		double PY;
		Vector_2<K> AP = cls_pt - pt;
		double AP_length = sqrt(AP.squared_length());

		PY = AP_length*cos(angle) + sqrt((EPSILON * EPSILON) - (AP_length * AP_length)*(sin(angle)*sin(angle)));
		
		double red_tail = sqrt(squared_distance(pt, p2));
		reduction = PY - red_tail;
		/*
		cout << "PY: " << PY << endl;
		cout << "AP_length: " << AP_length << endl;
		cout << "angle: " << angle << endl;
		cout << "reduction: " << reduction << endl;
		cout << "SUM1: " << AP_length*cos(angle) << endl;
		cout << "SUM2: " << sqrt((EPSILON * EPSILON) - (AP_length * AP_length)*(sin(angle)*sin(angle))) << endl;
		cout << "SUM3: " << (EPSILON * EPSILON) - (AP_length * AP_length)*(sin(angle)*sin(angle)) << endl;
		cout << "distance: " << sqrt(squared_distance(getSegment(a), pt)) << endl;
		*/
		v = (v * reduction) / sqrt(v.squared_length());

		np = p2 + v;
	}

	if (getSegment(a).squared_length() <= (reduction * reduction) || getSegment(a).squared_length() <= EPSILON*EPSILON){
		//segments.erase(segments.begin() + 2 * a + 1);
		//segments.erase(segments.begin() + 2 * a);
		segments[2 * a] = segments[2 * a + 1];
 		return false;
	}

	
	//cout << getSegment(a).squared_length() << endl;

	idx2 = pts.size() / 2;
	Add_Vertex(np[0], np[1]);

	//pts.push_back(np[0]);
	//pts.push_back(np[1]);

	//cout << "New Point in 1: " << np << endl;
	//cout << "reduction: " << reduction << endl;


	segments[2 * (a) + 1] = idx2;

	return true;


}

//reduces the segment away from the source()
bool Proc::shorten_segment_v1_ii(int a, int b, Point_2<K> pt, double angle, double reduction){
	Point_2<K> p1, p2;
	int idx2;

	//cout << reduction << "\t";

	p1 = Point_2<K>(pts[2 * segments[2 * a]], pts[2 * segments[2 * a] + 1]);
	p2 = Point_2<K>(pts[2 * segments[2 * a + 1]], pts[2 * segments[2 * a + 1] + 1]);
	Vector_2<K> v = p2 - p1;

	//cout << "first v: " << v << endl;
	//make v the length of the reduction
	v = (v * reduction) / sqrt(v.squared_length());

	//cout << "second v: " << v << ", reduction" << reduction << endl;
	//np = 'X'
	Point_2<K> np = p1 + v;

	Segment_2<K> seg2 = getSegment(b);

	double dist_x = sqrt(squared_distance(seg2, np));

	if (dist_x <= (EPSILON + EPS)){
		;//this is the right point
	}
	else {
		Point_2<K> cls_pt;
		if (sqrt(squared_distance(seg2.source(), np)) <= (dist_x + EPS)) {
			cls_pt = seg2.source();
		}
		else {
			cls_pt = seg2.target();
		}

		double PY;
		Vector_2<K> AP = cls_pt - pt;
		double AP_length = sqrt(AP.squared_length());

		PY = AP_length*cos(angle) + sqrt((EPSILON * EPSILON) - (AP_length * AP_length)*(sin(angle)*sin(angle)));

		double red_tail = sqrt(squared_distance(pt, p1));
		reduction = PY - red_tail;

		v = (v * reduction) / sqrt(v.squared_length());

		np = p1 + v;

	}

	if (getSegment(a).squared_length() <= (reduction * reduction) || getSegment(a).squared_length() <= EPSILON*EPSILON){
		//segments.erase(segments.begin() + 2 * a + 1);
		//segments.erase(segments.begin() + 2 * a);
		segments[2 * a] = segments[2 * a + 1];
		return false;
	}

	idx2 = (pts.size() / 2);

	Add_Vertex(np[0], np[1]);
	//pts.push_back(np[0]);
	//pts.push_back(np[1]);

	//cout << idx2 << endl;
	
	segments[2 * (a)] = idx2;

	return true;
}

/***DATA STRUCTURE HELPER FUNCTIONS***/

//helper functions for CGAL types
Segment_2<K> Proc::getSegment(int i) {
	return Segment_2<K>(P2(pts[2 * segments[2 * i]], pts[2 * segments[2 * i] + 1]), P2(pts[2 * segments[2 * i + 1]], pts[2 * segments[2 * i + 1] + 1]));
}

Line_2<K> Proc::getLine(int i) {
	return Line_2<K>(P2(pts[2 * segments[2 * i]], pts[2 * segments[2 * i] + 1]), P2(pts[2 * segments[2 * i + 1]], pts[2 * segments[2 * i + 1] + 1]));
}

int Proc::getNumberOfPoints() {
	return pts.size() / 2;
}

int Proc::getNumberOfSegments() {
	return segments.size() / 2;
}

int Proc::getNumSegPts() {
	return num_seg_pts;
}

void Proc::Erase_Segment(int x) {
	Segment_2<K> seg = getSegment(x);

	//Vector_2<K> v = (Vector_2<K>(seg.source()[0], seg.source()[1]) + Vector_2<K>(seg.target()[0], seg.target()[1])) / 2;
	//extra_pts.push_back(v[0]);
	//extra_pts.push_back(v[1]);

	segments.erase(segments.begin() + 2 * x + 1);
	segments.erase(segments.begin() + 2 * x);
}

void Proc::Add_Vertex(double x, double y)
{
	pts.push_back(x);
	pts.push_back(y);
	fix_v.push_back(false);
}

double * Proc::PointsToDoubleArray() {
	double * pointlist = (double *)malloc(pts.size() * sizeof(double));

	for (int i = 0; i < pts.size(); i++) {
		pointlist[i] = pts[i];
	}

	return pointlist;
}

int * Proc::SegmentsToIntArray() {
	int * seglist = (int *)malloc(segments.size() * sizeof(int));

	for (int i = 0; i < segments.size(); i++) {
		seglist[i] = segments[i];
	}

	return seglist;
}

//deals with the floating point errors that arise mainly from CGAL::intersection()
void Proc::deal_with_floating_point_errors() {

	int max_x = img->cols;
	int max_y = img->rows;

	for (int i = 0; i < pts.size(); i += 2) {

		if ((modulusX(pts[i] - 0)) < EPS) {
			pts[i] = 0;
		}

		if ((modulusX(pts[i] - max_x)) < EPS) {

			pts[i] = max_x;
		}

		if ((modulusX(pts[i + 1] - 0)) < EPS) {
			pts[i + 1] = 0;
		}

		if ((modulusX(pts[i + 1] - max_y)) < EPS) {
			pts[i + 1] = max_y;
		}
	}
}

//returns true if a segment has an endpoint outside the boundary
bool Proc::OutsideBoundary(Segment_2<K>& seg) {

	if (seg.source()[0] < 0 || seg.source()[0] > img->cols ||
		seg.source()[1] < 0 || seg.source()[1] > img->rows ||
		seg.target()[0] < 0 || seg.target()[0] > img->cols ||
		seg.target()[1] < 0 || seg.target()[1] > img->rows) {
		return true;
	}
	else
		return false;
}

/***MATHEMATICAL HELPER FUNCTIONS***/

//returns the smallest positive value of 2 doubles
double minimum(double x, double y) {
	if (y < x && y >= 0) {
		return y;
	}
	else if (x >= 0)
		return x;
	else {
		cout << "Two negative angles computed in minimum()" << endl;
		return x;
	}
}

//returns |x|
double modulusX(double x) {
	if (x >= 0) {
		return x;
	}
	else return -x;
}

//checks whether the associated acute angle of an angle is greater than MIN_ANGLE
bool Proc::checkAngle(double angle) {

	if (getMinAngle(angle) > MIN_ANGLE)
		return true;
	else
		return false;

}

//takes an angle in (0, 2*M_PI) and returns an associated acute angle
double Proc::getMinAngle(double angle) {

	if (angle > M_PI) {
		angle = (2 * M_PI) - angle;
	}

	return minimum(angle, modulusX(M_PI - angle));

}

//helper functions to sort integers
bool DecreasingInts(int const& p1, int const& p2) {
	return p1 > p2;
}

bool IncreasingInts(int const& p1, int const& p2) {
	return p1 < p2;
}

/***DRAWING FUNCTIONS***/

Mat Proc::DrawSegments(int rows, int cols) {

	Mat img(SCALE_FACTOR * rows, SCALE_FACTOR * cols, CV_8UC3, Scalar::all(255));

	Segment_2<K> seg;
	int linewidth = 4;

	for (int i = 0; i < segments.size() / 2; i++) {
		seg = getSegment(i);

		line(img, Point2d(SCALE_FACTOR*seg.source()[0], SCALE_FACTOR*seg.source()[1]), Point2d(SCALE_FACTOR*seg.target()[0], SCALE_FACTOR*seg.target()[1]), Red, linewidth, CV_AA);
	}

	Mat res;
	resize(img, res, Size(0, 0), 0.2, 0.2, CV_INTER_AREA);

	return res;
}

Mat Proc::DrawSegmentsOnImage(Mat& imgX) {

	//Mat img(SCALE_FACTOR * img.rows, SCALE_FACTOR * img.cols, CV_8UC3, Scalar::all(255));

	Mat img, img2;
	resize(imgX, img2, Size(0, 0), 10, 10, CV_INTER_AREA);

	cvtColor(img2, img, COLOR_GRAY2RGB);

	Segment_2<K> seg;
	int linewidth = 4;

	for (int i = 0; i < segments.size() / 2; i++) {
		seg = getSegment(i);

		line(img, Point2d(SCALE_FACTOR*seg.source()[0], SCALE_FACTOR*seg.source()[1]), Point2d(SCALE_FACTOR*seg.target()[0], SCALE_FACTOR*seg.target()[1]), Red, linewidth, CV_AA);
	}

	Mat res;
	resize(img, res, Size(0, 0), 0.2, 0.2, CV_INTER_AREA);

	return res;
}
