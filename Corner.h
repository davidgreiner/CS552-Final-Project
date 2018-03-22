#pragma once

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