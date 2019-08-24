/*
 * BinaryReciprocityModel.cpp
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#include "BinaryReciprocityModel.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#define MATHLIB_STANDALONE
#include <Rmath.h>

using namespace std;
using namespace mixrerg;
using namespace boost::numeric::ublas;

namespace mixrerg {

BinaryReciprocityModel::BinaryReciprocityModel() {
	numOfPiParameterSets = 3;
	cout << "Creating BinaryReciprocityModel" << endl;
}

BinaryReciprocityModel::~BinaryReciprocityModel() {
}

void BinaryReciprocityModel::setModel(int _numOfClasses, double _minTau,
		double _minAlpha, double _minPi) {

	cout << "BinaryReciprocityModel::setModel" << endl;

	// Set up basic information
	cout << "Setting up basic information" << endl;
	numOfClasses = _numOfClasses;
	minTau = _minTau;
	minAlpha = _minAlpha;
	minPi = _minPi;

	// Set up temporary statistics
	cout << "Setting up statistics" << endl;
	// the generalized vector of coordinate vectors is faster for insert
	SparseGVOCoordinateV _stat10 = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (1, 0)
	SparseGVOCoordinateV _stat01 = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (0, 1)
	SparseGVOCoordinateV _stat11 = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (1, 1)
	SparseGVOCoordinateV _stat00 = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // != (0, 0)
	// also remove self edges
	for (int i = 0; i < numOfVertices; i++) {
		_stat00.insert_element(i, i, 1);
	}

	cout << "Calculating statistics " << endl;
	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = networkMatrix.begin1(); i1 != networkMatrix.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			int i = i2.index1();
			int j = i2.index2();
			short y_ij = *i2;
			short y_ji = networkMatrix(j, i);
			// (1,0): stat10 = stat01'
			if (y_ij == 1 && y_ji == 0) {
				_stat10.insert_element(i, j, 1);
				_stat00.insert_element(i, j, 1);
				_stat01.insert_element(j, i, 1);
				_stat00.insert_element(j, i, 1);
			}
			// (1,1): stat11 is symmetric
			if (i < j && y_ij == 1 && y_ji == 1) {
				_stat11.insert_element(i, j, 1);
				_stat00.insert_element(i, j, 1);
				_stat11.insert_element(j, i, 1);
				_stat00.insert_element(j, i, 1);
			}
		}
	}

	cout << "Converting to compressed matrices" << endl;
	// convert to the compressed matrix which is is faster for lookup and matrix operations
	stat10 = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat10.nnz()); // (1, 0)
	stat10.assign(_stat10);

	stat01 = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat01.nnz()); // (0, 1)
	stat01.assign(_stat01);

	stat11 = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat11.nnz()); // (1, 1)
	stat11.assign(_stat11);

	stat00 = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat00.nnz()); // != (0, 0)
	stat00.assign(_stat00);

	cout << "Creating parameter matrices" << endl;
	tau = DoubleMatrix(numOfVertices, numOfClasses);
	prevTau = DoubleMatrix(numOfVertices, numOfClasses);

	alpha = DoubleVector(numOfClasses);
	prevAlpha = DoubleVector(numOfClasses);

	pi10 = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi10 = DoubleMatrix(numOfClasses, numOfClasses);

	pi11 = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi11 = DoubleMatrix(numOfClasses, numOfClasses);

	pi00 = DoubleMatrix(numOfClasses, numOfClasses);
}

void BinaryReciprocityModel::resetParameters() {
	// set up parameters
	cout << "Resetting tau randomly" << endl;
	tau.clear();
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			tau(i, k) = runif(0, 1);
	normalizeTau(tau, minTau);

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi10.clear();
	pi11.clear();
	pi00.clear();
	runModelEstimationMStep();
}

void BinaryReciprocityModel::resetParameters(const UnsignedIntVector& Z) {
	// set up parameters
	cout << "Resetting tau using Z" << endl;
	tau.clear();
	for (int i = 0; i < numOfVertices; i++)
		tau(i, Z(i)) = 1;
	normalizeTau(tau, minTau);

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi10.clear();
	pi11.clear();
	pi00.clear();
	runModelEstimationMStep();
}

void BinaryReciprocityModel::resetParameters(const DoubleMatrix& initTau) {
	// set up parameters
	cout << "Resetting tau using initTau" << endl;
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			tau(i, k) = initTau(i, k);

	alpha.clear();
	pi10.clear();
	pi11.clear();
	pi00.clear();
	runModelEstimationMStep();
}

void BinaryReciprocityModel::updateTau(DoubleMatrix& new_tau,
		NetworkSparseMatrix& stat, const DoubleMatrix& tau,
		const DoubleMatrix& logPi, DoubleMatrix& temp) {

	temp.clear();

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = stat.begin1(); i1 != stat.end1(); ++i1)
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int l = 0; l < numOfClasses; l++)
				temp(i, l) += tau(j, l);
		}

	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			for (int l = 0; l < numOfClasses; l++)
				new_tau(i, k) += logPi(k, l) * temp(i, l);
}

void BinaryReciprocityModel::updateTauByNegativeReflection(
		DoubleMatrix& new_tau, NetworkSparseMatrix& stat,
		const DoubleMatrix& tau, const DoubleMatrix& logPi, DoubleMatrix& temp) {

	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	for (int i = 0; i < numOfVertices; i++)
		for (int l = 0; l < numOfClasses; l++)
			temp(i, l) = tauL(l);

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = stat.begin1(); i1 != stat.end1(); ++i1)
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int l = 0; l < numOfClasses; l++)
				temp(i, l) -= tau(j, l);
		}

	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			for (int l = 0; l < numOfClasses; l++)
				new_tau(i, k) += logPi(k, l) * temp(i, l);
}

void BinaryReciprocityModel::runFixedPointEstimationEStep() {

	DoubleMatrix new_tau = DoubleMatrix(numOfVertices, numOfClasses);
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			new_tau(i, k) = log(alpha(k));

	DoubleMatrix temp = DoubleMatrix(numOfVertices, numOfClasses);

	DoubleMatrix logPi10 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi10, logPi10);
	updateTau(new_tau, stat10, tau, logPi10, temp);

	DoubleMatrix logPi10_t = DoubleMatrix(numOfClasses, numOfClasses);
	logTransposedMatrix(pi10, logPi10_t);
	updateTau(new_tau, stat01, tau, logPi10_t, temp);

	DoubleMatrix logPi11 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi11, logPi11);
	updateTau(new_tau, stat11, tau, logPi11, temp);

	DoubleMatrix logPi00 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi00, logPi00);
	updateTauByNegativeReflection(new_tau, stat00, tau, logPi00, temp);

	normalizeLogTau2Tau(new_tau, minTau);

	prevTau.swap(tau);
	tau.swap(new_tau);
}

bool BinaryReciprocityModel::isTauSignificantlyChanged(double _tauPrecision) {
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			if (fabs(tau(i, k) - prevTau(i, k)) > _tauPrecision)
				return true;
	return false;
}

double BinaryReciprocityModel::getLargestTauChange() {
	double largestChange = 0;
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			if (fabs(tau(i, k) - prevTau(i, k)) > largestChange)
				largestChange = fabs(tau(i, k) - prevTau(i, k));
	return largestChange;
}

void BinaryReciprocityModel::updatePi(DoubleMatrix& pi,
		NetworkSparseMatrix& stat, const DoubleMatrix& tau,
		const DoubleMatrix& sumTaus) {
	// reset to 0
	pi.clear();

	// compute numerator sum
	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = stat.begin1(); i1 != stat.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int k = 0; k < numOfClasses; k++)
				for (int l = 0; l < numOfClasses; l++)
					pi(k, l) += tau(i, k) * tau(j, l);
		}
	}

	// divided by denominator sum
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			pi(k, l) /= sumTaus(k, l);
}

void BinaryReciprocityModel::runModelEstimationMStep() {

	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);

	prevAlpha.swap(alpha);
	alpha.assign(tauL);
	for (unsigned int k = 0; k < alpha.size(); k++)
		alpha(k) /= numOfVertices;
	normalizeVector(alpha, minAlpha);

	DoubleMatrix sumTaus = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			sumTaus(k, l) = 0;
			for (int i = 0; i < numOfVertices; i++)
				sumTaus(k, l) += tau(i, k) * (tauL(l) - tau(i, l));
		}

	prevPi10.swap(pi10);
	updatePi(pi10, stat10, tau, sumTaus);

	prevPi11.swap(pi11);
	updatePi(pi11, stat11, tau, sumTaus);

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l <= k; l++) {
			pi00(k, l) = 1 - pi10(k, l) - pi10(l, k) - pi11(k, l);
			pi00(l, k) = pi00(k, l);
		}

	// check min and normalize again
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			if (pi10(k, l) < minPi)
				pi10(k, l) = minPi;
			if (pi11(k, l) < minPi)
				pi11(k, l) = minPi;
			if (pi00(k, l) < minPi)
				pi00(k, l) = minPi;
		}
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l <= k; l++) {
			double norms = pi00(k, l) + pi10(k, l) + pi10(l, k) + pi11(k, l);
			pi10(k, l) = pi10(k, l) / norms;
			pi10(l, k) = pi10(l, k) / norms;
			pi11(k, l) = pi11(k, l) / norms;
			pi11(l, k) = pi11(k, l);
			pi00(k, l) = pi00(k, l) / norms;
			pi00(l, k) = pi00(k, l);
		}
}

bool BinaryReciprocityModel::areAlphaPiSignificantlyChanged(
		double _paramPrecision) {
	for (int k = 0; k < numOfClasses; k++)
		if (fabs(alpha(k) - prevAlpha(k)) > _paramPrecision)
			return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi10(k, l) - prevPi10(k, l)) > _paramPrecision)
				return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi11(k, l) - prevPi11(k, l)) > _paramPrecision)
				return true;

	return false;
}

double BinaryReciprocityModel::getCompleteLogLikelihood() {

	double CLL = 0;

	// tau * [ log(alpha) - log(tau) ]
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			CLL += tau(i, k) * log(alpha(k) / tau(i, k));
	//cout << "tau * [ log(alpha) - log(tau) ] " << CLL << endl;

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;

	// y_ij = 1 and y_ji = 0
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat10.begin1(); i1 != stat10.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += temp * log(pi10(k, l) / pi00(k, l));
		}
	//cout << "y_ij = 1 and y_ji = 0 " << CLL << endl;

	// y_ij = 1 and y_ji = 1
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat11.begin1(); i1 != stat11.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += .5 * temp * log(pi11(k, l) / pi00(k, l));
		}
	//cout << "y_ij = 1 and y_ji = 1 " << CLL << endl;

	// norm term
	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (int i = 0; i < numOfVertices; i++)
				temp += tau(i, k) * (tauL(l) - tau(i, l));
			CLL += .5 * temp * log(pi00(k, l));
		}
	//cout << "norm term " << CLL << endl;

	return CLL;
}

double BinaryReciprocityModel::getICL() {
	double CLL = getCompleteLogLikelihood();
	double ICL = CLL - 0.5 * (numOfClasses - 1.0) * log(numOfVertices) - 0.5
			* (numOfClasses / 2.0) * (1.0 + 3 * numOfClasses) * log(
			(numOfVertices * (numOfVertices - 1.0)) / 2.0);
	return ICL;
}

void BinaryReciprocityModel::copyParameters(DoubleMatrix& bestTau,
		DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs) {
	bestTau.assign(tau);
	bestAlpha.assign(alpha);
	bestPIs(0).assign(pi10);
	bestPIs(1).assign(pi11);
	bestPIs(2).assign(pi00);
}

void BinaryReciprocityModel::copyPreviousParameters(DoubleMatrix& _prevTau,
		DoubleVector& _prevAlpha, VectorOfDoubleMatrices& _prevPIs) {
	_prevTau.assign(prevTau);
	_prevAlpha.assign(prevAlpha);
	_prevPIs(0).assign(prevPi10);
	_prevPIs(1).assign(prevPi11);
}

void BinaryReciprocityModel::printPIs() {
	cout.precision(6);
	cout << "pi10: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi10(k, l) << "  ";
		cout << endl;
	}
	cout << "pi11: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi11(k, l) << "  ";
		cout << endl;
	}
}

void BinaryReciprocityModel::printTau() {
	cout << "Printing TAU: " << endl;
	cout.precision(6);
	for (int i = 0; i < numOfVertices; i++) {
		cout << i << ": ";
		for (int k = 0; k < numOfClasses; k++)
			cout << tau(i, k) << "  ";
		cout << endl;
	}
}

void BinaryReciprocityModel::printAlpha() {
	cout.precision(6);
	cout << "Printing ALPHA: " << endl;
	cout << alpha(0);
	for (int k = 1; k < numOfClasses; k++)
		cout << ",  " << alpha(k);
	cout << endl;
}

}
