#pragma once

#include <vector>

#include "Vertex.h"

class Face
{
public:
	double area;
	std::vector<Vertex*> vertices;
};