#pragma once

#include <vector>

class Edge
{
public:
	int index;
	Vertex *verts[2];
	int ntris;
	Triangle **tris;
	double length;

	Vertex *middle;
	double opposite_angles[2];
};