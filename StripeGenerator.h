#pragma once

#include <vector>
#include <complex>
#include <Eigen/Dense>

#include "Surface.h"

using namespace Eigen;

class StripeGenerator
{
public:
	StripeGenerator();
	~StripeGenerator();

	void generatePattern(Surface* surface, VectorXd edge_lengths, VectorXcd unit_vectors, VectorXd target_line_frequency);

private:
	float vertexAngles(Surface* surface, VectorXd edge_lengths);
	void edgeData(Surface* surface, VectorXd edge_lengths, float theta, VectorXcd unit_vectors, VectorXd target_line_frequency);
	MatrixXd energyMatrix(Surface* surface, VectorXd edge_lengths, MatrixXd omega, MatrixXd s);
	MatrixXd massMatrix(Surface* surface, VectorXd edge_lengths);
	void principleEigenvector(MatrixXd energy_matrix, MatrixXd mass_matrix);
	void textureCoordinates(Surface* surface, Eigenvector, MatrixXd omega, MatrixXd s);
};

