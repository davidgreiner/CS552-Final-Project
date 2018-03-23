#pragma once
#include "icVector.h"
#include "ply.h"
#include "Triangle.h"
#include "Edge.h"
#include "Vertex.h"
#include "Corner.h"

class Surface
{
public:
	Surface();
	Surface(FILE *file);
	~Surface();

	void write_file(FILE *file);
	void initialize();
	void finalize();
	Triangle* find_common_edge(Triangle *f1, Vertex *v1, Vertex *v2);
	void create_edge(Vertex *v1, Vertex *v2);
	void create_edges();
	void vertex_to_tri_ptrs();
	Triangle* other_triangle(Edge *edge, Triangle *tri);
	void order_vertex_to_tri_ptrs(Vertex *v);
	int face_to_vertex_ref(Triangle *f, Vertex *v);
	void create_pointers();
	void calc_bounding_sphere();
	void calc_moments();
	void calc_edge_length();
	void calc_face_normals_and_area();
	void calc_corner_table();
	void test_corner_table();

	std::vector<Triangle*> triangles;
	std::vector<Edge*> edges;
	std::vector<Vertex*> vertices;
	std::vector<Corner*> corners;

private:
	int index;

	int max_tris;

	int orientation;
	int max_verts;

	int max_edges;

	icVector3 center;
	icVector3 gravity_center;
	std::vector<icVector3> moment_vertices;
	std::vector<icVector3> normal_vertices;
	double radius;
	double area;

	icVector3 min, max;

	int seed;

	PlyOtherProp *vert_other, *face_other;
};

