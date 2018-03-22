#pragma once

#ifndef __VERTEX_H__
#define __VERTEX_H__

class Triangle;

class Vertex {
public:
	double x, y, z;
	int index;

	int ntris;
	Triangle **tris;
	int max_tris;

	icVector3 normal;
	void *other_props;

public:
	Vertex(double xx, double yy, double zz) { x = xx; y = yy; z = zz;}
};
#endif /* __VERTEX_H__ */