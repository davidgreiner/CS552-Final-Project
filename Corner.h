#pragma once

#ifndef __CORNER_H__
#define __CORNER_H__

class Vertex;
class Edge;
class Triangle;

class Corner {
public:
	int index;
	Vertex *vertex;
	Edge *edge;
	Triangle *triangle;
	Corner *n, *p, *o;

	double angle;
	bool equals(Corner *c) {
		return (c->index == this->index) && (c->vertex == this->vertex);
	}
};
#endif /* __CORNER_H__ */