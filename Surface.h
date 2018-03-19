#pragma once
#include "Face.h"
#include "Edge.h"
#include "Vertex.h"

class Surface
{
public:
	Surface();
	~Surface();

	std::vector<Face*> faces;
	std::vector<Edge*> edges;
	std::vector<Vertex*> vertices;
};

