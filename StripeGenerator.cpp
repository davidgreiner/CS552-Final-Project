#include "StripeGenerator.h"



StripeGenerator::StripeGenerator()
{
}


StripeGenerator::~StripeGenerator()
{
}

MatrixXd StripeGenerator::energyMatrix(Surface* surface, VectorXd edge_lengths, MatrixXd omega, MatrixXd s)
{
	int doubleVertices = 2 * surface->numVertices();
	MatrixXd energyMatrix(doubleVertices, doubleVertices);
	energyMatrix = MatrixXd::Zero(doubleVertices, doubleVertices);

	return energyMatrix;
}

MatrixXd StripeGenerator::massMatrix(Surface* surface, VectorXd edge_lengths)
{
	int doubleVertices = 2 * surface->numVertices();
	MatrixXd massMatrix(doubleVertices, doubleVertices);
	massMatrix = MatrixXd::Zero(doubleVertices, doubleVertices);

	return massMatrix;
}

void StripeGenerator::principleEigenvector(MatrixXd energy_matrix, MatrixXd mass_matrix)
{
	Eigen::MatrixXd L(energy_matrix.llt().matrixL());

	VectorXd x = VectorXd::Random(energy_matrix.size);

}
