#pragma once

#include <vector>

class Edge
{
public:
	std::vector<double> opposite_angles;
	std::vector<Vertex*> vertices;
};