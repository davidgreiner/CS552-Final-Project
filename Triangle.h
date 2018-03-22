#pragma once

#include <vector>

#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "Vertex.h"
#include "Edge.h"

class Triangle
{
public:
	int index;
	int nverts;
	Vertex *verts[3];
	Edge *edges[3];

	icVector3 normal;
	void *other_props;
	double angle[3];
	double area;
};

#endif /* __TRIANGLE_H__ */