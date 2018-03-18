#pragma once

#include <vector>
#include <complex>

class StripeGenerator
{
public:
	StripeGenerator();
	~StripeGenerator();

	void generatePattern(Surface surface, std::vector<float> edge_lengths, std::vector<std::complex<float>> unit_vectors, std::vector<float> target_line_frequency);

private:
	float vertexAngles(Surface surface, std::vector<float> edge_lengths);
	void edgeData(Surface surface, std::vector<float> edge_lengths, float theta, std::vector<std::complex<float>> unit_vectors, std::vector<float> target_line_frequency);
	void energyMatrix(Surface surface, std::vector<float> edge_lengths, std::vector<float> omega, std::vector<float> s);
	void massMatrix(Surface surface, std::vector<float> edge_lengths);
	void principleEigenvector(A, B);
	void textureCoordinates(Surface surface, Eigenvector, std::vector<float> omega, std::vector<float> s);
};

