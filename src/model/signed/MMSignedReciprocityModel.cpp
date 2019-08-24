/*
 * MMSignedReciprocityModel.cpp
 *
 *  Created on: Jan 25, 2010
 *      Author: duyvu
 */

#include "MMSignedReciprocityModel.h"

using namespace std;

namespace mixrerg {

MMSignedReciprocityModel::MMSignedReciprocityModel() {
	// TODO Auto-generated constructor stub
	cout << "Creating MMSignedReciprocityModel" << endl;
}

MMSignedReciprocityModel::~MMSignedReciprocityModel() {
	// TODO Auto-generated destructor stub
}

void MMSignedReciprocityModel::runFixedPointEstimationEStep() {

	//cout << "MMSignedReciprocityModel::runFixedPointEstimationEStep()" << endl;

	static DoubleMatrix A = DoubleMatrix(numOfVertices, numOfClasses);
	static DoubleMatrix s = DoubleMatrix(numOfVertices, numOfClasses);

	// Calculate the quadratic coefficients

	// Compute the norm term, i.e. \pi_kl^0
	DoubleMatrix logPi5 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi5, logPi5);
	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++) {
			A(i, k) = 0;
			for (int l = 0; l < numOfClasses; l++)
				A(i, k) += (logPi5(k, l) * (tauL(l) - tau(i, l)));
		}

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;

	// y_ij = -1, y_ji = 0
	DoubleMatrix logPi1 = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi1(k, l) = log(pi1(k, l) / pi5(k, l));
	for (i1_t i1 = stat1a.begin1(); i1 != stat1a.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi1(k, l);
		}
	}

	// y_ij = 0, y_ji = -1
	for (i1_t i1 = stat1b.begin1(); i1 != stat1b.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi1(l, k);
		}
	}

	// y_ij = 1, y_ji = 0
	DoubleMatrix logPi2 = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi2(k, l) = log(pi2(k, l) / pi5(k, l));
	for (i1_t i1 = stat2a.begin1(); i1 != stat2a.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi2(k, l);
		}
	}

	// y_ij = 0, y_ji = 1
	for (i1_t i1 = stat2b.begin1(); i1 != stat2b.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi2(l, k);
		}
	}

	// y_ij = -1, y_ji = 1
	DoubleMatrix logPi3 = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi3(k, l) = log(pi3(k, l) / pi5(k, l));
	for (i1_t i1 = stat3a.begin1(); i1 != stat3a.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi3(k, l);
		}
	}

	// y_ij = 1, y_ji = -1
	for (i1_t i1 = stat3b.begin1(); i1 != stat3b.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi3(l, k);
		}
	}

	// y_ij = -1, y_ji = -1
	DoubleMatrix logPi4m = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi4m(k, l) = log(pi4m(k, l) / pi5(k, l));
	for (i1_t i1 = stat4m.begin1(); i1 != stat4m.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi4m(k, l);
		}
	}

	// y_ij = 1, y_ji = 1
	DoubleMatrix logPi4p = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			logPi4p(k, l) = log(pi4p(k, l) / pi5(k, l));
	for (i1_t i1 = stat4p.begin1(); i1 != stat4p.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					A(i, k) += tau(j, l) * logPi4p(k, l);
		}
	}

	// Finalize by subtracting half of from 1 dividing tau_{ik}
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++) {
			// In theory, A(i, k) must be negative or 0.
			if (A(i, k) > 0) { // In reality, A(i, k) can be greater than 0 because of numerical precision.
				A(i, k) = 0; // Therefore, we cut it off to 0 in this case
				cout << "A(" << i << ", " << k << ") is < 0." << endl;
			}
			A(i, k) = 1 - A(i, k) / 2;
			A(i, k) /= tau(i, k);
		}

	// calculate the linear coefficients
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

bool MMSignedReciprocityModel::isTauSignificantlyChanged(double _tauPrecision) {
	// since we only update tau one time in the variational E step
	return false;
}

}
