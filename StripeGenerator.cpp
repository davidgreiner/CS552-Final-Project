#include "StripeGenerator.h"
#define _USE_MATH_DEFINES
#include <math.h>



StripeGenerator::StripeGenerator()
{
}


StripeGenerator::~StripeGenerator()
{
}

MatrixXd StripeGenerator::energyMatrix(Surface* surface, VectorXd edge_lengths, MatrixXd omega, MatrixXd s)
{
	int doubleVertices = 2 * surface->vertices.size();
	MatrixXd energyMatrix(doubleVertices, doubleVertices);
	energyMatrix = MatrixXd::Zero(doubleVertices, doubleVertices);

	for (auto edge : surface->edges) {
		std::vector<double> angles = edge->opposite_angles;
		double omega = 0.5 * (tan(M_PI_2 - angles[0]) + tan(M_PI_2 - angles[1]));
		for (auto vertex : edge->vertices) {
			energyMatrix(vertex->index, vertex->index) = energyMatrix(vertex->index, vertex->index) + omega;
		}
		if (s(edge->vertices[0]->index, edge->vertices[1]->index) >= 0) {
			energyMatrix(edge->vertices[0]->index, edge->vertices[1]->index) = -omega * exp(omega);
		}
		else {
			energyMatrix(edge->vertices[0]->index, edge->vertices[1]->index) = -omega * exp(omega);
		}
		energyMatrix(edge->vertices[1]->index, edge->vertices[0]->index) = energyMatrix(edge->vertices[0]->index, edge->vertices[1]->index);
	}

	return energyMatrix;
}

MatrixXd StripeGenerator::massMatrix(Surface* surface, VectorXd edge_lengths)
{
	int doubleVertices = 2 * surface->vertices.size();
	MatrixXd massMatrix(doubleVertices, doubleVertices);
	massMatrix = MatrixXd::Zero(doubleVertices, doubleVertices);

	for (auto face : surface->faces) {
		double area = face->area;
		for (auto vertex : face->vertices)
			massMatrix(vertex->index, vertex->index) = massMatrix(vertex->index, vertex->index) + area / 3.;
		}
	}

	return massMatrix;
}

VectorXd StripeGenerator::principleEigenvector(MatrixXd energy_matrix, MatrixXd mass_matrix)
{
	Eigen::MatrixXd L(energy_matrix.llt().matrixL());

	VectorXd x = VectorXd::Random(energy_matrix.size);

	for (int iter = 0; iter < 20; iter++)
	{
		x = mass_matrix * x;
		SimplicialLDLT<SparseMatrix<double>> cholesky;
		cholesky.compute(L);
		x = cholesky.solve(x);
	}

	return x;
}
