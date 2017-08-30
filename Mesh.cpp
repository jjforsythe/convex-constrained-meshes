//MESH.CPP
//FOR EXTERNAL RELEASE

#include "Mesh.h"

//Mesh constructor - takes a triangulateio struct as input
Mesh::Mesh(triangulateio * tio, Mat image, String f, int n, string outputFolder, int pt_cut){
	filename = f;

	outFolder = outputFolder;

	/*
	outfile << fixed;
	outfile.open(outFolder + "faces_data_bsd.txt", ios::app);
	outfile << filename << "\t";

	outfile << n << "\t";
	*/

	//commented out for now
	cvtColor(image, iMesh, COLOR_GRAY2RGB);
	//iMesh = image;
	//add the vertices to the mesh
	map<int, MVH> vhandles;
	for (int i = 0; i < tio->numberofpoints; i++){
		auto vh = mesh.add_vertex(MP(tio->pointlist[2 * i], tio->pointlist[2 * i + 1], 0));

		//output the pointmarker trait
		if (tio->pointmarkerlist[i] >= pt_cut + 2) {
			//think this is showing isolated vertices, commented out for now (07/06/16)
			//cout << "Point marker: " << tio->pointmarkerlist[i] << endl;
			mesh.data(vh).fixed = true;
		}

		vhandles.insert(make_pair(i, vh));
	}

	//add the faces to the mesh
	vector<MVH> face_vhandles;
	for (int i = 0; i < tio->numberoftriangles; i++){
		face_vhandles.clear();
		for (int j = 0; j < 3; j++)
		{
			auto iv = vhandles.find(tio->trianglelist[3 * i + j]);
			face_vhandles.push_back(iv->second);
		}
		mesh.add_face(face_vhandles);
	}

	//cout << "Stats after Shewchuk triangulation:" << endl;
	//outfile << mesh.n_faces() << "\t";
	faces_after_triangle = mesh.n_faces();
	
	//iterate all vertices to set the '.point' trait
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
		mesh.data(*v_it).point = Point2d(mesh.point(*v_it)[0], mesh.point(*v_it)[1]);
	}
	
	//Set MyTraits: Triangle returns a list of all constrained segments in segmentlist[], each given by the two indices of the end vertices
	//we cycle over all outgoing halfedges attached to the first vertex of each constrained segment and check whether the half-edge points to the segment's second vertex
	//then we set the edge 'constrained' attribute to 'true' and increment 'constedges' for each vertex (constedges holds the number of constrained edges attached to a vertex)
	//I check against tester which should be incremented once for every constrained segment
	int tester = 0;
	for (int i = 0; i < tio->numberofsegments; i++){
		auto vh = mesh.vertex_handle(tio->segmentlist[2 * i]);
		//vertex is attached to constrained edge so we update constedges

		mesh.data(vh).constedges++;

		// circulate around the current vertex
		for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(vh); voh_it.is_valid(); ++voh_it)
		{
			auto vhhe = mesh.to_vertex_handle(*voh_it);
			if (vhhe.idx() == tio->segmentlist[2 * i + 1]){
				mesh.data(mesh.edge_handle(*voh_it)).constrained = true;
				mesh.data(mesh.edge_handle(*voh_it)).segment_idx = tio->segmentmarkerlist[i];

				//update constedges for the second vertex
				mesh.data(vhhe).constedges++;
				tester++;
			}
		}
	}

	int angle_check = Check_Every_Angle_Int();
	//outfile << angle_check << "\t";

	if (tester != tio->numberofsegments)
		cout << "Error in Mesh initialization" << endl;

	//send the number of triangles to cout
	cout << "Number of triangles: " << mesh.n_faces() << endl;

	//Save initial mesh
	//DrawMeshOnWhite("_triangulation");
	/*
	//iterate all vertices to set the '.point' trait
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
		mesh.data(*v_it).point = Point2d(mesh.point(*v_it)[0], mesh.point(*v_it)[1]);
	}
	

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++){
		mesh.data(*e_it).length = mesh.calc_edge_length(*e_it);
	}

	int ind = 1;
	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++){
		mesh.data(*e_it).index = ind;
		ind++;
	}
	*/

	double t = (double)getTickCount();
	//Run through the refinements to the mesh
	
	RefineMesh();
	//cout << "Stats after Refinement:" << endl;
	//outfile << mesh.n_faces() << "\t";
	
	angle_check = Check_Every_Angle_Int();
	//outfile << angle_check << "\t";
	
	//merge faces
	Merge_Faces(mesh, 0.01, 0.0000000001);
	//DrawMeshOnWhite("_simple_merge");
	//outfile << mesh.n_faces() << "\t";
	faces_after_merge = mesh.n_faces();

	angle_check = Check_Every_Angle_Int();
	//outfile << angle_check << "\t";

	/*
	ind = 1;
	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++){
		mesh.data(*e_it).index = ind;
		ind++;
	}
	*/
	//modify faces
	Modify_Faces(mesh, 0.01, 0.0000000001);
	
	//Remove_Deg2Vertices(mesh, 0.01, 0.0000000001);
	for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
		mesh.data(*v_it).point = Point2d(mesh.point(*v_it)[0], mesh.point(*v_it)[1]);
		mesh.data(*v_it).constedges = 0;
		for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(*v_it); voh_it.is_valid(); ++voh_it)
		{
			if (mesh.data(mesh.edge_handle(*voh_it)).constrained)
				mesh.data(*v_it).constedges++;
		}
	}

	for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++){
		mesh.data(*e_it).length = mesh.calc_edge_length(*e_it);
	}
	//DrawMeshOnWhite("_white_modify");
	//RefineMesh();
	

	//cout << "Stats after Modifying Faces:" << endl;
	//outfile << mesh.n_faces() << "\t";
	//iterate all vertices to set the '.point' trait
	
	Remove_Deg2Vertices(mesh, 0.01, 0.0000000001);
	faces_after_modify = mesh.n_faces();

	//DrawMesh("_ccm");
	//DrawMeshOnWhite("_polygonal_mesh");

	//DrawMeshOnWhite("_white_final");
	//DrawMesh("_final");

	cout << "Check_Simple_Polygons(): " << Check_Simple_Polygons() << endl;
	/*
	cout << "Writing to " << outputFolder + filename << ".off" << endl;
	if (!OpenMesh::IO::write_mesh(mesh, outputFolder  + filename + "_mesh" + ".off"))
	{
		std::cerr << "write error\n";
		exit(1);
	}
	*/
	
	//outfile.close();

	f1 = faces_after_triangle / faces_after_merge;
	f2 = faces_after_merge / faces_after_modify;

	f1s += f1;
	f2s += f2;

	cout << "CCM Faces: " << mesh.n_faces() << endl;

	angle_check = Check_Every_Angle_Int();

	//outfile << angle_check << "\n";

	if (angle_check == 0)
		cout << "All angles <= MIN_ANGLE" << endl;
	else
		cout << "Mesh contains angles > MIN_ANGLE" << endl;

	//cout << "Does Mesh satisfy MIN_ANGLE? " << angle_check << endl;
}

bool Mesh::Check_Simple_Polygons(){

	vector<vector<double>> vertices;

	bool all_good = true;

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++) {

		Polygon_2<K> p;
		
		for (auto fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); fv_it++) {
			
			OpenMesh::Vec3f v = mesh.point(*fv_it);
			//vector<double> temp_v;
			//temp_v.push_back(v[0]);
			//temp_v.push_back(v[1]);

			//vertices.push_back(temp_v);

			Point_2<K> point(v[0],v[1]);
			p.push_back(point);

		}
		if (!p.is_simple()) {
			cout << "Non-simple polygon error!!!" << endl;
			all_good = false;
		}
	}

	return all_good;
}


double Mesh::f1s = 0;
double Mesh::f2s = 0;
/*
Mesh::Mesh(){

}
*/
void Mesh::RefineMesh(){
	//normals are undefined until we call update_face_normals()
	mesh.update_face_normals();

	//adjustTJuncs2();

	resolveCloseEndpoints2();

}

void Mesh::adjustTJuncs2(){
	//iterate over all faces to find all faces with exactly 1 constrained edge where the opposite vertex is attached to exactly 1 constrained edge
	//then check for an intersection between the constrained segment on the face and the line from the constrained segment attached to the opposite vertex
	//if it exists we add a vertex at the intersection
	int numconstedges = 0;
	int vindexes[3];
	int eindexes[2];
	MEH meh1, meh2;
	MHH mhh, mhh2;
	MVH mvh;
	OpenMesh::Vec3f point1, point2, point3, point4;
	for (MyMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it){
		if (mesh.status(*f_it).deleted()) continue;
		//we only want to work with faces which are triangles here
		if (mesh.valence(*f_it) == 3){
			numconstedges = 0;
			for (MyMesh::FaceEdgeIter fe_it = mesh.fe_iter(*f_it); fe_it.is_valid(); ++fe_it){
				if (mesh.data(*fe_it).constrained){
					numconstedges++;
					meh1 = *fe_it;
				}
			}
			if (numconstedges == 1)
			{
				//get the indices of the edge's vertices
				mhh = mesh.halfedge_handle(meh1, 0);
				eindexes[0] = mesh.to_vertex_handle(mhh).idx();
				eindexes[1] = mesh.from_vertex_handle(mhh).idx();

				MyMesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it);
				for (int i = 0; i < 3; i++){
					if (fv_it.is_valid()){
						vindexes[i] = fv_it->idx();
						fv_it++;
					}
					else
						cout << "101 ERROR" << endl;
				}

				for (int j = 0; j < 3; j++){
					if ((vindexes[j] != eindexes[0]) && (vindexes[j] != eindexes[1])){
						//check there is only one constrained edge attached to the vertex
						mvh = mesh.vertex_handle(vindexes[j]);
						if (mesh.data((mvh)).constedges == 1){
							for (MyMesh::VertexEdgeIter ve_it = mesh.ve_iter(mvh); ve_it.is_valid(); ++ve_it){
								if (mesh.data(*ve_it).constrained){
									//we've found the constrained edge at the vertex
									meh2 = *ve_it;
								}
							}
							//see if line from the segment intersects with the segment on the triangle
							//now (hopefully) meh1 points to the segment on triangle, meh2 points to the segment attached to the opposite vertex
							point1 = mesh.point(mesh.vertex_handle(vindexes[(j + 1) % 3]));
							point2 = mesh.point(mesh.vertex_handle(vindexes[(j + 2) % 3]));

							mhh2 = mesh.halfedge_handle(meh2, 0);
							point3 = mesh.point(mesh.to_vertex_handle(mhh2));
							point4 = mesh.point(mesh.from_vertex_handle(mhh2));

							Segment_2<K> segment(P2(point1[0], point1[1]), P2(point2[0], point2[1]));
							Line_2<K> line(P2(point3[0], point3[1]), P2(point4[0], point4[1]));

							auto pintersect = intersection(line, segment);

							//P2 * ppp = boost::get<Point_2<K> >(&*pintersect);
							if (pintersect) {
								if (const Segment_2<K>* s = boost::get<Segment_2<K>>(&*pintersect)) {
									// handle segment
								}
								else {
									const P2 * ppp = boost::get<Point_2<K> >(&*pintersect);
									// handle point and test whether we can straighten out a line
									P2 p1(point1[0], point1[1]);
									P2 p2(point2[0], point2[1]);

									OpenMesh::Vec3f ppp1(ppp->x(), ppp->y(), 0);

									Comparison_result cr = compare_distance_to_point(*ppp, p1, p2);

									bool straightened = false;

									auto vhz = mesh.vertex_handle(vindexes[(j + 1) % 3]);
									auto vhz2 = mesh.vertex_handle(vindexes[(j + 2) % 3]);

									int edge_si = mesh.data(meh1).segment_idx;

									if ((cr == SMALLER) && (checkForSegment(vhz, edge_si))){
										//p1 closer than p2, so move point1 to ppp
										OpenMesh::Vec3f ppp1(ppp->x(), ppp->y(), 0);


										OpenMesh::Vec3f save_point = mesh.point(vhz);
										mesh.set_point(vhz, ppp1);
										mesh.data(vhz).point = Point2d(mesh.point(vhz)[0], mesh.point(vhz)[1]);

										if (isVertexPositionValid(vhz)){
											for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(vhz); voh_it.is_valid(); ++voh_it){
												if (mesh.to_vertex_handle(*voh_it).idx() == vindexes[j]){
													mesh.data(mesh.edge_handle(*voh_it)).constrained = true;
													mesh.data(mesh.to_vertex_handle(*voh_it)).constedges++;
													mesh.data(mesh.from_vertex_handle(*voh_it)).constedges++;
													mesh.data(vhz).fixed = true;
													//edge.length not set yet so do not need to update
												}
											}
											straightened = true;
										}
										else{
											mesh.set_point(vhz, save_point);
											mesh.data(vhz).point = Point2d(mesh.point(vhz)[0], mesh.point(vhz)[1]);
										}


									}
									else if ((cr == LARGER) && (checkForSegment(vhz2, edge_si))){
										//p2 closer than p1, so move point2 to ppp
										OpenMesh::Vec3f ppp2(ppp->x(), ppp->y(), 0);


										OpenMesh::Vec3f save_point = mesh.point(vhz2);
										mesh.set_point(vhz2, ppp2);
										mesh.data(vhz2).point = Point2d(mesh.point(vhz2)[0], mesh.point(vhz2)[1]);
										if (isVertexPositionValid(vhz2)){
											for (MyMesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(vhz2); voh_it.is_valid(); ++voh_it){
												if (mesh.to_vertex_handle(*voh_it).idx() == vindexes[j]){
													mesh.data(mesh.edge_handle(*voh_it)).constrained = true;
													mesh.data(mesh.to_vertex_handle(*voh_it)).constedges++;
													mesh.data(mesh.from_vertex_handle(*voh_it)).constedges++;
													mesh.data(vhz2).fixed = true;
													//edge.length not set yet so do not need to update
												}
											}
											straightened = true;
										}
										else{
											mesh.set_point(vhz2, save_point);
											mesh.data(vhz2).point = Point2d(mesh.point(vhz2)[0], mesh.point(vhz2)[1]);
										}
									}

									if (straightened) continue;

									//code for case we can't straighten out a line
									MVH vh = mesh.add_vertex(ppp1);

									//carry through data onto new edges
									int temp_seg_idx = mesh.data(meh1).segment_idx;

									mesh.split_edge(meh1, vh);
									mesh.data(vh).constedges = 3;
									mesh.data(vh).point = Point2d(mesh.point(vh)[0], mesh.point(vh)[1]);

									//halfedge handles for insert_edge
									MHH hh1, hh2;
									for (auto fh_it = mesh.fh_iter(*f_it); fh_it.is_valid(); fh_it++){
										if (mesh.to_vertex_handle(*fh_it) == vh){
											hh1 = *fh_it;
											hh2 = mesh.prev_halfedge_handle(*fh_it);

											//set data for the split edge
											MEH e1 = mesh.edge_handle(hh1);
											mesh.data(e1).constrained = true;
											mesh.data(e1).segment_idx = temp_seg_idx;

											MEH e2 = mesh.edge_handle(mesh.next_halfedge_handle(hh1));
											mesh.data(e2).constrained = true;
											mesh.data(e2).segment_idx = temp_seg_idx;
										}
									}
									MHH hh3 = mesh.insert_edge(hh1, hh2);
									MEH e3 = mesh.edge_handle(hh3);
									mesh.data(e3).constrained = true;
									mesh.data(e3).segment_idx = mesh.data(meh2).segment_idx;

									//now need to check angles, old face is *f_it new face can get from hh1
									MFH f1 = *f_it;
									MFH f2 = mesh.face_handle(hh1);

									bool small_angle = false;
									double angle;
									for (auto fh_it = mesh.fh_iter(f1); fh_it.is_valid(); fh_it++){
										Find_Angle(mesh, *fh_it, EPS, angle);
										if (angle < MIN_ANGLE){
											small_angle = true;
										}
									}
									for (auto fh_it = mesh.fh_iter(f2); fh_it.is_valid(); fh_it++){
										Find_Angle(mesh, *fh_it, EPS, angle);
										if (angle < MIN_ANGLE){

											small_angle = true;
										}
									}


									if (small_angle){

										MFH f3 = mesh.remove_edge(e3);

										for (auto fh_it = mesh.fh_iter(f3); fh_it.is_valid(); fh_it++){
											if (mesh.from_vertex_handle(*fh_it) == vh){
												if (!isCollapseOK(*fh_it)){
													cout << "Collapse failed!";
												}
											}
										}
										for (auto fh_it = mesh.fh_iter(f3); fh_it.is_valid(); fh_it++){
											double angle2;
											Find_Angle(mesh, *fh_it, EPS, angle2);
											if (angle2 < MIN_ANGLE)
												cout << "smaller angle detected" << endl;
										}

										mesh.garbage_collection(true, true, true);

									}

								}
							}
						}
					}
				}

			}
		}
	}


	cout << "adjustTJuncs completes!" << endl;
}

int Mesh::Check_Every_Angle_Int()
{

	double angle;
	int i;
	int small_angles = 0;

	for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); f_it++){
		i = 0;
		for (auto fh_it = mesh.fh_iter(*f_it); fh_it.is_valid(); fh_it++){
			Find_Angle(mesh, *fh_it, EPS, angle);

			if (angle < (MIN_ANGLE)){
				//uncomment for each small angle size and location
				//cout << "ANGLE TOO SMALL: " << angle << endl;
				//cout << "At point: " << mesh.data(mesh.to_vertex_handle(*fh_it)).point << endl;
				//cout << "face valence: " << mesh.valence(mesh.face_handle(*fh_it)) << endl;
				small_angles++;
			}
		}
	}

	return small_angles;
}

bool Mesh::isCollapseOK(MHH hh){

	if (!mesh.is_collapse_ok(hh)) return false;

	MVH v1 = mesh.to_vertex_handle(hh);
	MVH v2 = mesh.from_vertex_handle(hh);

	//save the original position of v2 then change its data
	OpenMesh::Vec3f p1 = mesh.point(v2);
	mesh.set_point(v2, mesh.point(v1));
	mesh.data(v2).point = Point2d(mesh.point(v1)[0], mesh.point(v1)[1]);
	bool failed = false;

	failed = !isVertexPositionValid(v2, hh);
	/*
	for (auto vf_it = mesh.vf_iter(v2); vf_it.is_valid(); vf_it++){
		mesh.update_normal(*vf_it);
		if (mesh.normal(*vf_it)[2] == -1){
			mesh.set_point(v2, p1);
			failed = true;
		}
	}
	*/
	//if we failed the test revert the position back to how it was...
	if (failed){
		mesh.set_point(v2, p1);
		mesh.data(v2).point = Point2d(p1[0], p1[1]);
		return false;
	}
	else{
		mesh.set_point(v2, p1);
		
		int c2 = mesh.data(v2).constedges;
		bool con = mesh.data(mesh.edge_handle(hh)).constrained;
		mesh.collapse(hh);
		if (con)
			mesh.data(v1).constedges += (c2 - 2);
		else
			mesh.data(v1).constedges += c2;
		//update the lengths for the moved edges which haven't been deleted
		for (auto ve_it = mesh.ve_iter(v1); ve_it.is_valid(); ve_it++){
			if (!mesh.status(*ve_it).deleted()){
				mesh.data(*ve_it).length = mesh.calc_edge_length(*ve_it);
			}
		}

		collapses++;
		return true;
	}
}

bool Mesh::isCollapseToPointOK(MHH hh, OpenMesh::Vec3f ptest){
	if (!mesh.is_collapse_ok(hh)) return false;

	MVH v1 = mesh.to_vertex_handle(hh);
	MVH v2 = mesh.from_vertex_handle(hh);

	MFH f1 = mesh.face_handle(hh);
	MFH f2 = mesh.face_handle(mesh.opposite_halfedge_handle(hh));

	OpenMesh::Vec3f p1 = mesh.point(v1);
	OpenMesh::Vec3f p2 = mesh.point(v2);
	OpenMesh::Vec3f vec1 = p1 - ptest;
	OpenMesh::Vec3f vec2 = p2 - ptest;
	
	bool failed = false;

	double eps_sq = EPSILON * EPSILON;
	Vector_2<K> vv1(vec1[0], vec1[1]);
	Vector_2<K> vv2(vec2[0], vec2[1]);

	if (vv1.squared_length() > eps_sq || vv2.squared_length() > eps_sq){
		failed = true;
	}

	//set points and point traits
	mesh.set_point(v1, ptest);
	mesh.set_point(v2, ptest);
	mesh.data(v1).point = Point2d(ptest[0], ptest[1]);
	mesh.data(v2).point = Point2d(ptest[0], ptest[1]);

	
	double angle;

	//if the faces around either vertex are flipped or have small angles we set failed to true
	if (!failed)
		failed = (!isVertexPositionValid(v1, hh) || !isVertexPositionValid(v2, hh));

	for (auto vf_it = mesh.vf_iter(v1); vf_it.is_valid(); vf_it++){
		if (((*vf_it != f1) && (*vf_it != f2)) || (*vf_it == f1 && mesh.valence(f1) != 3) || (*vf_it == f2 && mesh.valence(f2) != 3)){
			if (!isFaceConvex(*vf_it)){
				failed = true;
			}
		}
	}

	for (auto vf_it = mesh.vf_iter(v2); vf_it.is_valid(); vf_it++){
		if (((*vf_it != f1) && (*vf_it != f2)) || (*vf_it == f1 && mesh.valence(f1) != 3) || (*vf_it == f2 && mesh.valence(f2) != 3)){
			if (!isFaceConvex(*vf_it)){
				failed = true;
			}
		}
	}

	if (failed){

		mesh.set_point(v1, p1);
		mesh.set_point(v2, p2);
		mesh.data(v1).point = Point2d(p1[0], p1[1]);
		mesh.data(v2).point = Point2d(p2[0], p2[1]);
		return false;
	}

	else{
		//update constedges for the vertex being collapsed

		MVH vh = mesh.to_vertex_handle(hh);

		mesh.data(vh).constedges++;
		mesh.collapse(hh);
		//update the lengths for the moved edges which haven't been deleted
		for (auto ve_it = mesh.ve_iter(vh); ve_it.is_valid(); ve_it++){
			if (!mesh.status(*ve_it).deleted()){
				mesh.data(*ve_it).length = mesh.calc_edge_length(*ve_it);
			}
		}
		return true;
	}
	
}

bool Mesh::isFaceConvex(MFH mfh){
	double angle;
	
	for (auto fh_it = mesh.fh_iter(mfh); fh_it.is_valid(); fh_it++){
		Find_Angle(mesh, *fh_it, EPS, angle);
		if (angle > (M_PI + DELTA))
			return false;
	}
	return true;
}

bool Mesh::isVertexPositionValid(MVH vh, MHH hh){

	MFH f1 = mesh.face_handle(hh);
	MFH f2 = mesh.face_handle(mesh.opposite_halfedge_handle(hh));

	double angle;

	for (auto vf_it = mesh.vf_iter(vh); vf_it.is_valid(); vf_it++){
		if (mesh.status(*vf_it).deleted()) continue;
		mesh.update_normal(*vf_it);
		if (mesh.normal(*vf_it)[2] == -1){
			return false;
		}

		if ((*vf_it != f1) && (*vf_it != f2) || ((mesh.valence(f1) != 3) && (*vf_it == f1)) || ((mesh.valence(f2) != 3) && (*vf_it == f2))){
			for (auto fh_it = mesh.fh_iter(*vf_it); fh_it.is_valid(); fh_it++){
				if (mesh.valence(*vf_it) == 3){
					if (Find_Angle(mesh, *fh_it, EPS, angle)){
						//cout << mesh.data(mesh.to_vertex_handle(hh)).point << "::" << angle << endl;
						if (angle < MIN_ANGLE){
							return false;
						}
					}
				}
				else{
					if (Find_Angle2(mesh, *fh_it, EPS, angle)){
						//cout << mesh.data(mesh.to_vertex_handle(hh)).point << "::" << angle << endl;
						if (angle < MIN_ANGLE){
							return false;
						}
					}
				}
			}
		}
		//cout << "face done" << endl;
	}

	return true;
}

bool Mesh::isVertexPositionValid(MVH vh){

	double angle;

	for (auto vf_it = mesh.vf_iter(vh); vf_it.is_valid(); vf_it++){
		if (mesh.status(*vf_it).deleted()) continue;
		mesh.update_normal(*vf_it);
		if (mesh.normal(*vf_it)[2] == -1){
			return false;
		}

		
		for (auto fh_it = mesh.fh_iter(*vf_it); fh_it.is_valid(); fh_it++){
			if (Find_Angle(mesh, *fh_it, EPS, angle)){
				if (angle < MIN_ANGLE){
					return false;
				}
			}
		}
		
	}

	return true;
}

bool Mesh::checkForSegment(MVH vh, int seg_idx){
	int count = 0;
	for (auto ve_it = mesh.ve_iter(vh); ve_it.is_valid(); ve_it++){
		if (mesh.data(*ve_it).segment_idx == seg_idx)
			count++;
	}

	if (count == 2)
		return true;
	else
		return false;
}

void Mesh::resolveCloseEndpoints2(){
	int index = 0;
	int count = 0;
	//find all short (less than epsilon) unconstrained edges attached to exactly one constrained edge at each end
	//compute the intersection and calculate whether it is feasible to use this intersection
	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it){
		if (mesh.status(*e_it).deleted()) continue;

		//set the indices for the edges
		mesh.data(*e_it).index = index;
		index++;

		//set the lengths of the edges
		double l = mesh.calc_edge_length(*e_it);
		mesh.data(*e_it).length = l;

		MHH mhi = mesh.halfedge_handle(*e_it, 0);
		int d1, d2;
		MVH v1, v2;
		MEH e1, e2;

		OpenMesh::Vec3f point1, point2, point3, point4;
		l = 5 * EPSILON;
		if (l <= (2*EPSILON)){
		//if (l <= (EPSILON)){
			if (!mesh.data(*e_it).constrained){
				v1 = mesh.to_vertex_handle(mhi);
				v2 = mesh.from_vertex_handle(mhi);
				d1 = mesh.data(v1).constedges;
				d2 = mesh.data(v2).constedges;
				if (d1 == 1 && d2 == 1){

					//iterate around the first vertex to find the constrained edge
					for (MyMesh::VertexEdgeIter ve_it = mesh.ve_iter(v1); ve_it.is_valid(); ++ve_it){
						if (mesh.data(*ve_it).constrained){
							//we've found the constrained edge at the vertex
							e1 = *ve_it;
						}
					}

					//iterate around the second vertex to find the constrained edge
					for (MyMesh::VertexEdgeIter ve_it = mesh.ve_iter(v2); ve_it.is_valid(); ++ve_it){
						if (mesh.data(*ve_it).constrained){
							//we've found the constrained edge at the vertex
							e2 = *ve_it;
						}
					}

					MHH hh1 = mesh.halfedge_handle(e1, 0);
					MHH hh2 = mesh.halfedge_handle(e2, 0);

					point1 = mesh.point(mesh.to_vertex_handle(hh1));
					point2 = mesh.point(mesh.from_vertex_handle(hh1));

					point3 = mesh.point(mesh.to_vertex_handle(hh2));
					point4 = mesh.point(mesh.from_vertex_handle(hh2));

					Line_2<K> line1(P2(point1[0], point1[1]), P2(point2[0], point2[1]));
					Line_2<K> line2(P2(point3[0], point3[1]), P2(point4[0], point4[1]));

					auto pintersect = intersection(line1, line2);

					if (pintersect) {
						if (const Line_2<K>* s = boost::get<Line_2<K>>(&*pintersect)) {
							// handle line
						}
						else {
							//handle point
							const P2 * ppp = boost::get<Point_2<K> >(&*pintersect);

							if(isCollapseToPointOK(mhi, OpenMesh::Vec3f(ppp->x(), ppp->y(), 0)))
								count++;
						}
					}
				}
			}
		}
	}
	mesh.garbage_collection(true, true, true);
	//cout << "Number of resolved close endpoints: " << count << endl;
}

bool Mesh::isVertexFaceConvex(MVH mvh){
	MHH hh0, hh1;
	double angle0, angle1;
	double eps = 0.0000000001;
	for (auto voh_iter = mesh.voh_iter(mvh); voh_iter.is_valid(); voh_iter++){
		hh0 = *voh_iter;
		//hh0 = mesh.halfedge_handle(*ve_iter, 0);mesh.halfedge
		Find_Angle(mesh, hh0, eps, angle0);
		hh1 = mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(hh0));
		Find_Angle(mesh, hh1, eps, angle1);
		if ((angle0 + angle1) >= M_PI)
			return false;
	}
	return true;
}

void Mesh::DrawMesh(String tag){

	//Mat tempImg(SCALE_FACTOR * 512, SCALE_FACTOR * 512, CV_8UC3 );
	//iMesh.copyTo(tempImg);

	//scale up the image using function below by an integer value keeping exact pixels colours, ie. output will have SCALE_FACTOR^2 more pixels
	int scale = SCALE_FACTOR;
	//int scale = 1;
	//Mat bigImg = scaleImage(iMesh, scale);

	Mat cm = imread("Input/SPM Grayscale/cameraman512.jpg", CV_LOAD_IMAGE_COLOR);
	Mat bigImg = scaleImage(cm, scale);

	//Mat bigImg(SCALE_FACTOR * 512, SCALE_FACTOR * 512, CV_8UC3);

	int linewidth = 12;

	//draw all edges in red with constrained edges drawn with width 2, unconstrained edges with width 1
	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it){
		MyTraits::Point p1 = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0)));
		MyTraits::Point p2 = mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0)));

		if (mesh.data(*e_it).constrained){
			;//line(iMesh, Point2d(p1[0], p1[1]), Point2d(p2[0], p2[1]), Red, 1);
			//linewidth = 1;
		}
		else
			;//linewidth = 1;
		
		line(bigImg, Point2d(scale*p1[0], scale*p1[1]), Point2d(scale*p2[0], scale*p2[1]), Red, linewidth, CV_AA);
	}

	//draw all vertices in blue
	//int radius = 3;
	//for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
	//	circle(bigImg, Point2d(scale * mesh.point(*v_it)[0], scale * mesh.point(*v_it)[1]), radius, Blue, -1, 8);
	//}

	//save the output

	Mat outImg;
	resize(bigImg, outImg, Size(0, 0), 0.1, 0.1, CV_INTER_AREA);

	if (!imwrite(outFolder + filename + "_mesh_" + to_string(SCALE_FACTOR) + "x" + tag + ".png", outImg)){
		cout << "Error writing image" << endl;
	}

	//displays the original image with LSD lines drawn
	imshow("OpenMesh", iMesh);

	//tempImg.copyTo(iMesh);
}

string Mesh::getOutfolder(){
	return outFolder;
}

void Mesh::DrawMeshOnWhite(String tag){

	int radius = 1;
	int linewidth = 4;
	int ulinewidth = 4;
	if (tag == "_triangulation")
		ulinewidth = 2;
	int borderwidth = 15;

	Mat img(SCALE_FACTOR * iMesh.size().height, SCALE_FACTOR* iMesh.size().width, CV_8UC3, Scalar::all(255));
	//Mat img((SCALE_FACTOR * iMesh.size().height) + linewidth / 2, (SCALE_FACTOR* iMesh.size().width) + linewidth / 2, CV_8UC3, Scalar::all(255));

	//draw all edges constrained edges in red, unconstrained edges in blue
	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it){
		MyTraits::Point p1 = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0)));
		MyTraits::Point p2 = mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0)));

		if (!mesh.data(*e_it).constrained){
			line(img, Point2d(SCALE_FACTOR*p1[0], SCALE_FACTOR*p1[1]), Point2d(SCALE_FACTOR*p2[0], SCALE_FACTOR*p2[1]), Blue, ulinewidth, CV_AA);
		}
			
	}

	for (MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it){
		MyTraits::Point p1 = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*e_it, 0)));
		MyTraits::Point p2 = mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*e_it, 0)));

		if (mesh.data(*e_it).constrained){
			line(img, Point2d(SCALE_FACTOR*p1[0], SCALE_FACTOR*p1[1]), Point2d(SCALE_FACTOR*p2[0], SCALE_FACTOR*p2[1]), Red, linewidth, CV_AA);
		}
	}

	//draw all vertices in blue
	/*
	for (MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it){
		circle(img, Point2d(SCALE_FACTOR * mesh.point(*v_it)[0], SCALE_FACTOR * mesh.point(*v_it)[1]), radius, Black, -1);
	}
	*/

	//save the output
	//if (!imwrite("Output/" + filename + "_mesh_" + to_string(SCALE_FACTOR) + "x" + tag + ".png", img)){
	
	Mat res;
	resize(img, res, Size(0,0), 0.2, 0.2, CV_INTER_AREA);

	if (!imwrite(outFolder + filename + tag + "_" + to_string(mesh.n_faces()) + ".png", res)){
		cout << "Error writing image" << endl;
	}
}

void Mesh::printMeshData(){
	cout << "Number of OpenMesh vertices: " << mesh.n_vertices() << endl;
	cout << "Number of OpenMesh faces: " << mesh.n_faces() << endl;
	cout << "Number of OpenMesh edges: " << mesh.n_edges() << endl;
	cout << endl;

	outfile << mesh.n_faces() << "\t";
}

//function which scales up an image by an integer value
Mat Mesh::scaleImage(Mat I, int scale){
	Mat scaleImg = Mat(scale * I.size().width, scale * I.size().height, CV_8UC3);
	//tested continuous at scale = 10

	// accept only char type matrices
	CV_Assert(I.depth() == CV_8U);

	int channels = I.channels();

	int nRows = I.rows;
	int nCols = I.cols * channels;

	if (I.isContinuous())
	{
		nCols *= nRows;
		nRows = 1;
	}

	int i, j, k, l;
	int xx, yy, nx, ny, jj;
	uchar * p, *q;
	for (i = 0; i < nRows; ++i)
	{
		p = scaleImg.ptr<uchar>(i);
		q = I.ptr<uchar>(i);
		for (j = 0; j < nCols; j += 3)
		{
			xx = (j / 3) % I.cols;
			yy = ((j / 3) - xx) / I.cols;
			nx = scale * xx;
			ny = scale * yy;
			jj = ny* scaleImg.cols * 3 + nx * 3;
			for (k = 0; k < scale; k++){
				for (l = 0; l < scale; l++){
					p[jj + scaleImg.cols * 3 * k + 3 * l] = q[j];
					p[jj + scaleImg.cols * 3 * k + 3 * l + 1] = q[j + 1];
					p[jj + scaleImg.cols * 3 * k + 3 * l + 2] = q[j + 2];
				}
			}
		}
	}

	return scaleImg;
}

// Angle at the to_vertex of a half-edge hh in the face defined by hh
bool Find_Angle(MyMesh const& pmesh, MHH hh, double eps, double& angle)
{
	if (!pmesh.is_valid_handle(hh)) return false;
	MVH vh = pmesh.to_vertex_handle(hh);
	MVH vh_prev = pmesh.from_vertex_handle(hh);
	MHH hh_next = pmesh.next_halfedge_handle(hh);
	MVH vh_next = pmesh.to_vertex_handle(hh_next);
	Point2d vec1 = pmesh.data(vh_next).point - pmesh.data(vh).point;
	Point2d vec2 = pmesh.data(vh_prev).point - pmesh.data(vh).point;
	return Find_Angle(vec1, vec2, eps, angle);
}

bool Find_Angle2(MyMesh const& pmesh, MHH hh, double eps, double& angle)
{


	if (!pmesh.is_valid_handle(hh)) return false;
	MVH vh = pmesh.to_vertex_handle(hh);
	MVH vh_prev = pmesh.from_vertex_handle(hh);
	MHH hh_next = pmesh.next_halfedge_handle(hh);
	MVH vh_next = pmesh.to_vertex_handle(hh_next);
	Point2d vec1 = pmesh.data(vh_next).point - pmesh.data(vh).point;
	Point2d vec2 = pmesh.data(vh_prev).point - pmesh.data(vh).point;

	Point2d zero_vec(0, 0);
	if (vec2 == zero_vec){
		vh_prev = pmesh.from_vertex_handle(pmesh.prev_halfedge_handle(hh));
		vec2 = pmesh.data(vh_prev).point - pmesh.data(vh).point;
	}
	if (vec1 == zero_vec){
		hh_next = pmesh.next_halfedge_handle(hh_next);
		vec1 = pmesh.data(vh_next).point - pmesh.data(vh).point;
	}

	return Find_Angle(vec1, vec2, eps, angle);
}

//Determinant function
double Det(Point2d v1, Point2d v2){
	return v1.x * v2.y - v2.x * v1.y;
}

// Angle from vec1 to vec2 (from 0 to 2pi)
//anticlockwise angle from vec1 to vec2
bool Find_Angle(Point2d vec1, Point2d vec2, double eps, double& angle)
{
	Point2d zero_vec(0, 0);
	if (vec1 == zero_vec || vec2 == zero_vec)
		return false;
	double prod = norm(vec1) * norm(vec2);
	if (prod < eps * eps)
	{
		cout << "Error in Find_Angle: prod=" << prod << endl;\
		return false;
	}
	double det = Det(vec1, vec2);
	double dot = vec1.dot(vec2);
	if (abs(det) < eps)
	{
		if (dot < 0) angle = M_PI; else angle = 0;
		return true;
	}
	double cosine = dot / prod;
	if (cosine >= 1) { angle = 0; return true; }
	if (cosine <= -1) { angle = M_PI; return true; }
	if (det > 0) angle = acos(cosine);
	else angle = 2 * M_PI - acos(cosine);
	return true;
}


// Sort all unconstrained edge in the order of their increasing length
bool Increasing_Unconstrained(MyMesh& pmesh, double max_length,
	vector<MEH>& sorted_edges)
{
	int ind;
	double length;
	map<int, MEH> ehandles;
	vector<IndexValue> sorted_lengths;
	for (auto e_it = pmesh.edges_begin(); e_it != pmesh.edges_end(); ++e_it)
	{
		if (pmesh.data(*e_it).constrained) continue;
		length = pmesh.data(*e_it).length;
		//if (length > max_length) continue;
		ind = pmesh.data(*e_it).index;
		ehandles.insert(make_pair(ind, *e_it));
		sorted_lengths.push_back(IndexValue(ind, length));
	}
	sort(sorted_lengths.begin(), sorted_lengths.end(), IncreasingValues);
	for (size_t e = 0; e < sorted_lengths.size(); e++)
		sorted_edges.push_back(ehandles[sorted_lengths[e].index]);
	return true;
}

bool DecreasingValues(IndexValue const& p1, IndexValue const& p2) {
	return p1.value > p2.value;
}
bool IncreasingValues(IndexValue const& p1, IndexValue const& p2) {
	return p1.value < p2.value;
}

//***new functions***
// Edges (starting from longest) between faces are removed if a new face is convex
void Merge_Faces(MyMesh& pmesh, double delta, double eps)
{
	MEH eh;
	MVH vh0, vh1;
	MHH hh, hh0, hh1, hh0prev, hh1prev, hh0next, hh1next;
	Point2d p0, p1;
	bool convex0, convex1;
	vector<MEH> sorted_edges;
	Increasing_Unconstrained(pmesh, 1, sorted_edges);
	while (sorted_edges.size() > 0)
	{
		int e = (int)sorted_edges.size() - 1;
		eh = sorted_edges[e];
		sorted_edges.erase(sorted_edges.begin() + e);
		if (!pmesh.is_valid_handle(eh)) continue;
		if (!pmesh.is_simple_link(eh)) continue;
		hh0 = pmesh.halfedge_handle(eh, 0);
		Merge_Faces_at_Edge(pmesh, hh0, convex0, convex1, delta, eps
			); // merged by removing 1 edge
	}
	//added to make it all work
	pmesh.garbage_collection(true, true ,true);
}

bool Merge_Faces_at_Edge(MyMesh& pmesh, MHH hh0, bool& convex0,
	bool& convex1, double delta, double eps)
{
	double angle00, angle01, angle10, angle11;
	MHH hh0prev = pmesh.prev_halfedge_handle(hh0);
	MHH hh1 = pmesh.opposite_halfedge_handle(hh0);
	MHH hh1prev = pmesh.prev_halfedge_handle(hh1);
	// find angles of the edge with its 4 neighboring edges
	Find_Angle(pmesh, hh0, eps, angle01);
	Find_Angle(pmesh, hh0prev, eps, angle00);
	Find_Angle(pmesh, hh1, eps, angle11);
	Find_Angle(pmesh, hh1prev, eps, angle10);
	//cout << angle00 + angle11 << "  " << angle01 + angle10 << endl;
	convex0 = (angle00 + angle11 <= M_PI + delta); // convexity of hh0 tail
	convex1 = (angle01 + angle10 <= M_PI + delta); // convexity of hh0 head
	if (convex0 && convex1) // both new angles after merging are convex
	{
		pmesh.request_edge_status();
		pmesh.request_face_status();
		pmesh.remove_edge(pmesh.edge_handle(hh0));
		return true;
	}
	return false;
}

void Modify_Faces(MyMesh& pmesh, double delta, double eps)
{
	MEH eh;
	MVH vh0, vh1;
	MHH hh, hh0, hh1, hh0prev, hh1prev, hh0next, hh1next, hhn;
	bool convex0, convex1, cut0, cut1;
	vector<MEH> sorted_edges;
	Increasing_Unconstrained(pmesh, 1, sorted_edges);

	
	while (sorted_edges.size() > 0)
	{
		int e = (int)sorted_edges.size() - 1;
		eh = sorted_edges[e];
		sorted_edges.erase(sorted_edges.begin() + e);
		if (!pmesh.is_valid_handle(eh)) continue;
		if (pmesh.status(eh).deleted()) continue;
		if (!pmesh.is_simple_link(eh)) continue;
		hh0 = pmesh.halfedge_handle(eh, 0);
		hh1 = pmesh.halfedge_handle(eh, 1);
		
		hhn = pmesh.next_halfedge_handle(hh0);
		
		//try merging by removing one edge
		if (Merge_Faces_at_Edge(pmesh, hh0, convex0, convex1, delta, eps)) {
			addFaceEdges(pmesh, pmesh.face_handle(pmesh.prev_halfedge_handle(hhn)), sorted_edges);
			continue;
				}
		vh0 = pmesh.from_vertex_handle(hh0);
		vh1 = pmesh.to_vertex_handle(hh0);
		//errors here - have removed Constrained_Degree()
		cut0 = !convex0 && ((pmesh.valence(vh0) == 3) || (pmesh.valence(vh0) == 4)) && (pmesh.data(vh0).constedges == 0) && !(pmesh.data(vh0).fixed);
		if (cut0) cut0 = Triangular_Cut_Possible(pmesh, hh1, delta, eps);
		if (!convex0 && !cut0) continue; // impossible to straighten non - convex angle 0

		MHH hh1next = pmesh.next_halfedge_handle(hh1);
		MHH hh1test = pmesh.prev_halfedge_handle(hh1);
		MHH hh1_opp_prev = pmesh.prev_halfedge_handle(pmesh.opposite_halfedge_handle(hh1));
		MHH hh1_fin = pmesh.next_halfedge_handle(pmesh.opposite_halfedge_handle(hh1next));
		OpenMesh::Vec3f thepoint1;

		if (pmesh.valence(vh0) == 4){
			MVH v1 = pmesh.to_vertex_handle(hh1next);
			MVH v2 = pmesh.from_vertex_handle(hh1_opp_prev);

			auto p1 = pmesh.data(v1).point;
			auto p2 = pmesh.data(v2).point;

			Segment_2<K> line1(P2(p1.x, p1.y), P2(p2.x, p2.y));

			MVH v3 = vh0;
			MVH v4 = pmesh.to_vertex_handle(hh1_fin);

			auto p3 = pmesh.data(v3).point;
			auto p4 = pmesh.data(v4).point;

			Segment_2<K> line2(P2(p3.x, p3.y), P2(p4.x, p4.y));

			auto pintersect = intersection(line1, line2);

			if (pintersect){

				if (const Segment_2<K>* s = boost::get<Segment_2<K>>(&*pintersect)) {
					// handle line
				}
				else {
					const P2 * ppp = boost::get<Point_2<K> >(&*pintersect);
					thepoint1 = OpenMesh::Vec3f(ppp->x(), ppp->y(), 0);
				}

			}
			else
				cut0 = false;

		}

		if (!convex0 && !cut0) continue;

		//errors here - have removed Constrained_Degree()
		cut1 = !convex1 && ((pmesh.valence(vh1) == 3) || (pmesh.valence(vh1) == 4)) && (pmesh.data(vh1).constedges == 0) && !(pmesh.data(vh1).fixed);
		if (cut1) cut1 = Triangular_Cut_Possible(pmesh, hh0, delta, eps);
		if (!convex1 && !cut1) continue; // impossible to straighten non - convex angle 1

			
		// remove the edge
		MHH hh0next = pmesh.next_halfedge_handle(hh0);		
		MHH hh0test = pmesh.prev_halfedge_handle(hh0);		
		MHH hh0_opp_prev = pmesh.prev_halfedge_handle(pmesh.opposite_halfedge_handle(hh0));
		MHH hh0_fin = pmesh.next_halfedge_handle(pmesh.opposite_halfedge_handle(hh0next));
		OpenMesh::Vec3f thepoint2;

		if (pmesh.valence(vh1) == 4){
			MVH v1 = pmesh.to_vertex_handle(hh0next);
			MVH v2 = pmesh.from_vertex_handle(hh0_opp_prev);

			auto p1 = pmesh.data(v1).point;
			auto p2 = pmesh.data(v2).point;

			Segment_2<K> line1(P2(p1.x, p1.y), P2(p2.x, p2.y));

			MVH v3 = vh1;
			MVH v4 = pmesh.to_vertex_handle(hh0_fin);

			auto p3 = pmesh.data(v3).point;
			auto p4 = pmesh.data(v4).point;

			Segment_2<K> line2(P2(p3.x, p3.y), P2(p4.x, p4.y));

			auto pintersect = intersection(line1, line2);

			if (pintersect){

				if (const Segment_2<K>* s = boost::get<Segment_2<K>>(&*pintersect)) {
					// handle line
				}
				else {
					const P2 * ppp = boost::get<Point_2<K> >(&*pintersect);
					thepoint2 = OpenMesh::Vec3f(ppp->x(), ppp->y(), 0);
				}

			}
			else
				cut1 = false;
		}

		if (!convex1 && !cut1) continue;

		MFH f1, f2, f3;

		pmesh.request_edge_status();
		pmesh.request_face_status();
		pmesh.remove_edge(eh);
		// straighten angle 0 - at vh0
		if (!convex0)
		{
			pmesh.request_vertex_status();
			pmesh.request_edge_status();
			pmesh.request_face_status();

			MHH hh1check = pmesh.prev_halfedge_handle(hh1next);

			if (pmesh.valence(vh0) == 2){
				
				if (pmesh.is_collapse_ok(hh1next)){
					

					///sorted_edges.push_back(pmesh.edge_handle(hh1check));

					pmesh.collapse(hh1next);

					
				}
			}
			else if (pmesh.valence(vh0) == 3){
					pmesh.set_point(vh0, thepoint1);
					pmesh.data(vh0).point = Point2d(pmesh.point(vh0)[0], pmesh.point(vh0)[1]);			

					f3 = pmesh.face_handle(pmesh.opposite_halfedge_handle(hh1_opp_prev));
					addFaceEdges(pmesh, f3, sorted_edges);
			}

			f1 = pmesh.face_handle(pmesh.opposite_halfedge_handle(hh1check));
			addFaceEdges(pmesh, f1, sorted_edges);

			f2 = pmesh.face_handle(hh1check);
			addFaceEdges(pmesh, f2, sorted_edges);

		}
		// straighten angle 1 - at vh1
		if (!convex1)
		{
			pmesh.request_vertex_status();
			pmesh.request_edge_status();
			pmesh.request_face_status();

			MHH hh0check = pmesh.prev_halfedge_handle(hh0next);

			if (pmesh.valence(vh1) == 2){
				
				if (pmesh.is_collapse_ok(hh0next)){
					MHH hh0check = pmesh.prev_halfedge_handle(hh0next);

					//sorted_edges.push_back(pmesh.edge_handle(hh0check));

					pmesh.collapse(hh0next);

				}
			}
			else if (pmesh.valence(vh1) == 3){			

					pmesh.set_point(vh1, thepoint2);
					pmesh.data(vh1).point = Point2d(pmesh.point(vh1)[0], pmesh.point(vh1)[1]);
						
					f3 = pmesh.face_handle(pmesh.opposite_halfedge_handle(hh0_opp_prev));
					addFaceEdges(pmesh, f3, sorted_edges);
			}


			f1 = pmesh.face_handle(pmesh.opposite_halfedge_handle(hh0check));
			addFaceEdges(pmesh, f1, sorted_edges);

			f2 = pmesh.face_handle(hh0check);
			addFaceEdges(pmesh, f2, sorted_edges);

		}
	}
	pmesh.garbage_collection(true, true, true);
}

void addFaceEdges(MyMesh& mesh, MFH fh, vector<MEH>& sorted_edges){

	if (mesh.is_valid_handle(fh) && !mesh.status(fh).deleted()){
		for (auto fh_it = mesh.fh_iter(fh); fh_it.is_valid(); fh_it++){
			if (!mesh.data(mesh.edge_handle(*fh_it)).constrained)
				sorted_edges.push_back(mesh.edge_handle(*fh_it));
			
			if ((mesh.valence(mesh.to_vertex_handle(*fh_it)) == 3) || (mesh.valence(mesh.to_vertex_handle(*fh_it)) == 4)){
			//if ((mesh.valence(mesh.to_vertex_handle(*fh_it)) == 3) ){
				for (auto ve_it = mesh.ve_iter(mesh.to_vertex_handle(*fh_it)); ve_it.is_valid(); ve_it++){
					if (!mesh.data(*ve_it).constrained)
						sorted_edges.push_back(*ve_it);
				}
			}
		}
	}
	return;
}

bool Triangular_Cut_Possible(MyMesh& pmesh, MHH hh, double delta, double eps)
{
	double angle_prev, angle_next, ext_prev, ext_next, opp_prev, opp_next;
	MHH hh_next = pmesh.next_halfedge_handle(hh);
	MHH hh_cur = pmesh.prev_halfedge_handle(
		pmesh.opposite_halfedge_handle(hh));
	MHH hh_prev = pmesh.prev_halfedge_handle(hh_cur);
	MHH hh_prevoppnext = pmesh.prev_halfedge_handle(pmesh.opposite_halfedge_handle(hh_next));
	MHH hh_oppcur = pmesh.opposite_halfedge_handle(hh_cur);
	Find_Angle(pmesh, hh_next, eps, angle_next);
	if (angle_next > M_PI + delta) {
		cout << "\nError in Triangular_Cut_Possible: angle_next = " << angle_next; }
		MVH vh_next = pmesh.to_vertex_handle(hh_next);
		Point2d v_next = pmesh.data(vh_next).point;
		MVH vh = pmesh.to_vertex_handle(hh);
		Point2d v = pmesh.data(vh).point;
		MVH vh_prev = pmesh.to_vertex_handle(hh_prev);
		Point2d v_prev = pmesh.data(vh_prev).point;
		Find_Angle(v - v_next, v_prev - v_next, eps, ext_next);
		Find_Angle(pmesh, hh_prevoppnext, EPS, opp_next);
		if (angle_next + ext_next > M_PI + delta) return false; // the new angle would be non - convex
		if ((opp_next - ext_next) < MIN_ANGLE) return false; //a small angle would be created
		Find_Angle(pmesh, hh_prev, eps, angle_prev);
		if (angle_prev > M_PI + delta) cout << "\nError in Triangular_Cut_Possible: angle_prev = " << angle_prev;
		Find_Angle(v_next - v_prev, v - v_prev, eps, ext_prev);
		Find_Angle(pmesh, hh_oppcur, EPS, opp_prev);
		if ((opp_prev - ext_prev) < MIN_ANGLE) return false; //a small angle would be created
		if (angle_prev + ext_prev > M_PI + delta) return false; // the new angle would be non - convex
			return true;
	}

vector<MEH> Mesh::Sort_Edges(vector<MEH>& edges){
	double length;
	int ind;
	map<int, MEH> ehandles;
	vector<IndexValue> sorted_lengths;
	vector<MEH> sorted_edges;

	for (int i = 0; i < edges.size(); i++){
		length = mesh.data(edges[i]).length;
		ind = mesh.data(edges[i]).index;
		ehandles.insert(make_pair(ind, edges[i]));
		sorted_lengths.push_back(IndexValue(ind, length));
	}

	sort(sorted_lengths.begin(), sorted_lengths.end(), IncreasingValues);

	for (size_t e = 0; e < sorted_lengths.size(); e++)
		sorted_edges.push_back(ehandles[sorted_lengths[e].index]);

	return sorted_edges;
}

bool Remove_Deg2Vertices(MyMesh& pmesh, double delta, double eps)
{
	double angle;
	MHH hh_in, hh_out;
	for (auto v_it = pmesh.vertices_begin(); v_it !=
		pmesh.vertices_end(); ++v_it)
	{
		if (pmesh.valence(*v_it) != 2) continue; // only deg 2 vertices are removed
			auto hh_out_it = pmesh.voh_iter(*v_it);
		auto hh_in_it = hh_out_it; hh_in_it++;
		hh_in = pmesh.opposite_halfedge_handle(*hh_in_it);
		Find_Angle(pmesh, hh_in, eps, angle);
		if (angle < M_PI - delta || angle > M_PI + delta) continue;
		// only vertices on a straight line are removed
		bool constrained = pmesh.data(pmesh.edge_handle(*hh_out_it)
			).constrained;
		constrained = constrained || pmesh.data(pmesh.edge_handle(
			hh_in)).constrained;
		if (constrained) pmesh.data(pmesh.edge_handle(hh_in)
			).constrained = true;
		pmesh.data(pmesh.edge_handle(hh_in)).length += pmesh.data(
			pmesh.edge_handle(*hh_out_it)).length;
		pmesh.request_vertex_status();
		pmesh.request_edge_status();
		pmesh.request_face_status();
		pmesh.collapse(*hh_out_it);
	}
	pmesh.garbage_collection();
	return true;
}