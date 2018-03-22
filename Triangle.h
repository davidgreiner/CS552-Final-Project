#pragma once

#include <vector>

#include "Vertex.h"
#include "Edge.h"

class Triangle
{
public:
	int index;
	int nverts;
	Vertex *verts[3];
	Edge *edges[3];

	double angle[3];
	float area;

	icVector3 normal;
	void *other_props;
	double angle[3];
	double area;
};