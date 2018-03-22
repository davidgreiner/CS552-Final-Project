#pragma once

#ifndef __EDGE_H__
#define __EDGE_H__

#include <vector>

class Triangle;
class Vertex;

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

#endif /* __EDGE_H__ */