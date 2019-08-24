/*
 * MMBinaryReciprocityModel.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: duyvu
 */

#include "MMBinaryReciprocityModel.h"

using namespace std;

namespace mixrerg {

MMBinaryReciprocityModel::MMBinaryReciprocityModel() {
	// TODO Auto-generated constructor stub
	cout << "Creating MMBinaryReciprocityModel" << endl;
}

MMBinaryReciprocityModel::~MMBinaryReciprocityModel() {
	// TODO Auto-generated destructor stub
}

void MMBinaryReciprocityModel::runFixedPointEstimationEStep() {

	//cout << "MMBinaryReciprocityModel::runFixedPointEstimationEStep()" << endl;

	static DoubleMatrix A = DoubleMatrix(numOfVertices, numOfClasses);
	static DoubleMatrix s = DoubleMatrix(numOfVertices, numOfClasses);

	// Calculate the quadratic coefficients

	// Compute the norm term, i.e. \pi_kl^0
	DoubleMatrix logPi00 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi00, logPi00);
	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++) {
			A(i, k) = 0;
			for (int l = 0; l < numOfClasses; l++)
				A(i, k) += (logPi00(k, l) * (tauL(l) - tau(i, l)));
		}

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;

	// y_ij = 1, y_ji = 0 AND y_ij = 0, y_ji = 1
	DoubleMatrix logPi10 = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi10(k, l) = log(pi10(k, l) / pi00(k, l));
	for (i1_t i1 = stat10.begin1(); i1 != stat10.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++) {
					A(i, k) += tau(j, l) * logPi10(k, l);
					A(j, l) += tau(i, k) * logPi10(l, k);
				}
		}
	}

	// y_ij = 1, y_ji = 1
	DoubleMatrix logPi11 = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi11(k, l) = log(pi11(k, l) / pi00(k, l));
	for (i1_t i1 = stat11.begin1(); i1 != stat11.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi11(k, l);
		}
	}

	// Finalize by subtracting half of from 1 dividing tau_{ik}
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++) {
			// In theory, A(i, k) must be negative or 0.
			//if (A(i, k) > 0) // In reality, A(i, k) can be greater than 0 because of numerical precision.
			//	A(i, k) = 0; // Therefore, we cut it off to 0 in this case
			A(i, k) = 1 - A(i, k) / 2;
			A(i, k) /= tau(i, k);
		}

	// Calculate the linear coefficients
	DoubleVector logAlpha = DoubleVector(numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		logAlpha(k) = log(alpha(k));
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			s(i, k) = 1 + logAlpha(k) - log(tau(i, k));

	DoubleMatrix new_tau = DoubleMatrix(numOfVertices, numOfClasses);
	solveQP(A, s, new_tau, minTau);

	normalizeTau(new_tau, minTau);

	prevTau.swap(tau);
	tau.swap(new_tau);
}

bool MMBinaryReciprocityModel::isTauSignificantlyChanged(double _tauPrecision) {
	// since we only update tau one time in the variational E step
	return false;
}

}
