/*
 * SignedReciprocityModel.cpp
 *
 *  Created on: Jan 17, 2010
 *      Author: duyvu
 */

#include "SignedReciprocityModel.h"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#define MATHLIB_STANDALONE
#include <Rmath.h>

using namespace std;
using namespace mixrerg;
using namespace boost::numeric::ublas;

namespace mixrerg {

SignedReciprocityModel::SignedReciprocityModel() {
	numOfPiParameterSets = 6;
	cout << "Creating SignedReciprocityModel" << endl;
}

SignedReciprocityModel::~SignedReciprocityModel() {
	// TODO Auto-generated destructor stub
}

void SignedReciprocityModel::setModel(int _numOfClasses, double _minTau,
		double _minAlpha, double _minPi) {

	cout << "SignedReciprocityModel::setModel" << endl;

	// set up basic information
	cout << "setting up basic information" << endl;
	numOfClasses = _numOfClasses;
	minTau = _minTau;
	minAlpha = _minAlpha;
	minPi = _minPi;

	// set up temp statistics
	cout << "setting up statistics" << endl;
	// the generalized vector of coordinate vectors is faster for insert
	SparseGVOCoordinateV _stat1a = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (-1, 0)
	SparseGVOCoordinateV _stat1b = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (0, -1)
	SparseGVOCoordinateV _stat2a = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (1, 0)
	SparseGVOCoordinateV _stat2b = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (0, 1)
	SparseGVOCoordinateV _stat3a = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (-1, 1)
	SparseGVOCoordinateV _stat3b = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (1, -1)
	SparseGVOCoordinateV _stat4m = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (-1, -1)
	SparseGVOCoordinateV _stat4p = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // (1, 1)
	SparseGVOCoordinateV _stat5 = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // != (0, 0)
	// also remove self-edge
	for (int i = 0; i < numOfVertices; i++) {
		_stat5.insert_element(i, i, 1);
	}
	SparseGVOCoordinateV _statMinus = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // y_ij = -1
	SparseGVOCoordinateV _statPlus = SparseGVOCoordinateV(numOfVertices,
			numOfVertices); // y_ij = +1

	cout
			<< "calculating statistics to generalized vectors of coordinate vectors"
			<< endl;
	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = networkMatrix.begin1(); i1 != networkMatrix.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			int i = i2.index1();
			int j = i2.index2();
			short y_ij = *i2;
			short y_ji = networkMatrix(j, i);

			// (-1,0): stat1b = stat1a'
			if (y_ij == -1 && y_ji == 0) {
				_stat1a.insert_element(i, j, 1);
				_stat1b.insert_element(j, i, 1);
				_stat5.insert_element(i, j, 1);
				_stat5.insert_element(j, i, 1);
			} else
			// (1,0): stat2b = stat2a'
			if (y_ij == 1 && y_ji == 0) {
				_stat2a.insert_element(i, j, 1);
				_stat2b.insert_element(j, i, 1);
				_stat5.insert_element(i, j, 1);
				_stat5.insert_element(j, i, 1);
			} else
			// (-1,1): stat3b = stat3a'
			if (y_ij == -1 && y_ji == 1) {
				_stat3a.insert_element(i, j, 1);
				_stat3b.insert_element(j, i, 1);
				_stat5.insert_element(i, j, 1);
				_stat5.insert_element(j, i, 1);
			}

			if (i < j) {
				// (-1,-1): stat4m
				if (y_ij == -1 && y_ji == -1) {
					_stat4m.insert_element(i, j, 1);
					_stat4m.insert_element(j, i, 1);
					_stat5.insert_element(i, j, 1);
					_stat5.insert_element(j, i, 1);
				} else
				// (1,1): stat4p
				if (y_ij == 1 && y_ji == 1) {
					_stat4p.insert_element(i, j, 1);
					_stat4p.insert_element(j, i, 1);
					_stat5.insert_element(i, j, 1);
					_stat5.insert_element(j, i, 1);
				}
			}

			if (y_ij == -1)
				_statMinus.insert_element(i, j, 1);
			else if (y_ij == 1)
				_statPlus.insert_element(i, j, 1);
		}
	}

	cout
			<< "copying statistics from generalized vectors of coordinate vectors to compressed matrices"
			<< endl;
	// the compressed matrix is faster for lookup and matrix operations
	stat1a = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat1a.nnz()); // (-1, 0)
	stat1a.assign(_stat1a);
	cout << "sum stat1a: " << sumAbsoluteNetworkSparseMatrix(stat1a) << endl;

	stat1b = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat1b.nnz()); // (0, -1)
	stat1b.assign(_stat1b);
	cout << "sum stat1b: " << sumAbsoluteNetworkSparseMatrix(stat1b) << endl;

	stat2a = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat2a.nnz()); // (1, 0)
	stat2a.assign(_stat2a);
	cout << "sum stat2a: " << sumAbsoluteNetworkSparseMatrix(stat2a) << endl;

	stat2b = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat2b.nnz()); // (0, 1)
	stat2b.assign(_stat2b);
	cout << "sum stat2b: " << sumAbsoluteNetworkSparseMatrix(stat2b) << endl;

	stat3a = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat3a.nnz()); // (-1, 1)
	stat3a.assign(_stat3a);
	cout << "sum stat3a: " << sumAbsoluteNetworkSparseMatrix(stat3a) << endl;

	stat3b = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat3b.nnz()); // (1, -1)
	stat3b.assign(_stat3b);
	cout << "sum stat3b: " << sumAbsoluteNetworkSparseMatrix(stat3b) << endl;

	stat4m = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat4m.nnz()); // (-1, -1)
	stat4m.assign(_stat4m);
	cout << "sum stat4m: " << sumAbsoluteNetworkSparseMatrix(stat4m) << endl;

	stat4p = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat4p.nnz()); // (1, 1)
	stat4p.assign(_stat4p);
	cout << "sum stat4p: " << sumAbsoluteNetworkSparseMatrix(stat4p) << endl;

	stat5 = NetworkSparseMatrix(numOfVertices, numOfVertices, _stat5.nnz()); // != (0, 0)
	stat5.assign(_stat5);
	cout << "sum stat5: " << sumAbsoluteNetworkSparseMatrix(stat5) << endl;

	statMinus = NetworkSparseMatrix(numOfVertices, numOfVertices,
			_statMinus.nnz()); // y_ij = -1
	statMinus.assign(_statMinus);
	cout << "sum statMinus: " << sumAbsoluteNetworkSparseMatrix(statMinus)
			<< endl;

	statPlus = NetworkSparseMatrix(numOfVertices, numOfVertices,
			_statPlus.nnz()); // y_ij = +1
	statPlus.assign(_statPlus);
	cout << "sum statPlus: " << sumAbsoluteNetworkSparseMatrix(statPlus)
			<< endl;

	cout << "Creating parameter matrices" << endl;
	tau = DoubleMatrix(numOfVertices, numOfClasses);
	prevTau = DoubleMatrix(numOfVertices, numOfClasses);

	alpha = DoubleVector(numOfClasses);
	prevAlpha = DoubleVector(numOfClasses);

	pi1 = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi1 = DoubleMatrix(numOfClasses, numOfClasses);

	pi2 = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi2 = DoubleMatrix(numOfClasses, numOfClasses);

	pi3 = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi3 = DoubleMatrix(numOfClasses, numOfClasses);

	pi4m = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi4m = DoubleMatrix(numOfClasses, numOfClasses);

	pi4p = DoubleMatrix(numOfClasses, numOfClasses);
	prevPi4p = DoubleMatrix(numOfClasses, numOfClasses);

	pi5 = DoubleMatrix(numOfClasses, numOfClasses);
}

void SignedReciprocityModel::resetParameters() {
	// set up parameters
	cout << "Resetting tau randomly" << endl;
	tau.clear();
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			tau(i, k) = runif(0, 1);
	normalizeTau(tau, minTau);

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi1.clear();
	pi2.clear();
	pi3.clear();
	pi4m.clear();
	pi4p.clear();
	pi5.clear();
	runModelEstimationMStep();
}

void SignedReciprocityModel::resetParameters(const UnsignedIntVector& Z) {
	// set up parameters
	cout << "Resetting tau using Z" << endl;
	tau.clear();
	for (int i = 0; i < numOfVertices; i++)
		tau(i, Z(i)) = 1;
	normalizeTau(tau, minTau);

	//printDoubleMatrix("tau after normalizing", tau);

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi1.clear();
	pi2.clear();
	pi3.clear();
	pi4m.clear();
	pi4p.clear();
	pi5.clear();
	runModelEstimationMStep();

	//printDoubleVector("alpha after estimation", alpha);
	//printDoubleMatrix("pi4p after estimation", pi4p);
}

void SignedReciprocityModel::resetParameters(const DoubleMatrix& initTau) {

	cout << "Resetting tau using initTau" << endl;
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			tau(i, k) = initTau(i, k);

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi1.clear();
	pi2.clear();
	pi3.clear();
	pi4m.clear();
	pi4p.clear();
	pi5.clear();
	runModelEstimationMStep();
}

void SignedReciprocityModel::updateTau(DoubleMatrix& new_tau,
		NetworkSparseMatrix& stat, const DoubleMatrix& tau,
		const DoubleMatrix& logPi, DoubleMatrix& temp) {

	temp.clear();

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = stat.begin1(); i1 != stat.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int l = 0; l < numOfClasses; l++) {
				temp(i, l) += tau(j, l);
			}
		}
	}

	for (int i = 0; i < numOfVertices; i++) {
		for (int k = 0; k < numOfClasses; k++) {
			for (int l = 0; l < numOfClasses; l++) {
				new_tau(i, k) += logPi(k, l) * temp(i, l);
			}
		}
	}

}

void SignedReciprocityModel::updateTauByNegativeReflection(
		DoubleMatrix& new_tau, NetworkSparseMatrix& stat,
		const DoubleMatrix& tau, const DoubleMatrix& logPi, DoubleMatrix& temp) {

	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	for (int i = 0; i < numOfVertices; i++)
		for (int l = 0; l < numOfClasses; l++)
			temp(i, l) = tauL(l);

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = stat.begin1(); i1 != stat.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
			unsigned int i = i2.index1();
			unsigned int j = i2.index2();
			for (int l = 0; l < numOfClasses; l++) {
				temp(i, l) -= tau(j, l);
			}
		}
	}

	for (int i = 0; i < numOfVertices; i++) {
		for (int k = 0; k < numOfClasses; k++) {
			for (int l = 0; l < numOfClasses; l++) {
				new_tau(i, k) += logPi(k, l) * temp(i, l);
			}
		}
	}

}

void SignedReciprocityModel::runFixedPointEstimationEStep() {

	//cout << "running FixedPointEstimationEStep..." << endl;

	DoubleMatrix new_tau = DoubleMatrix(numOfVertices, numOfClasses);
	//cout << "Initialize new tau" << endl;
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			new_tau(i, k) = log(alpha(k));

	DoubleMatrix temp = DoubleMatrix(numOfVertices, numOfClasses);

	//cout << "stat1a" << endl;
	DoubleMatrix logPi1 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi1, logPi1);
	updateTau(new_tau, stat1a, tau, logPi1, temp);

	//cout << "stat1b" << endl;
	DoubleMatrix logPi1_t = DoubleMatrix(numOfClasses, numOfClasses);
	logTransposedMatrix(pi1, logPi1_t);
	updateTau(new_tau, stat1b, tau, logPi1_t, temp);

	//cout << "stat2a" << endl;
	DoubleMatrix logPi2 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi2, logPi2);
	updateTau(new_tau, stat2a, tau, logPi2, temp);

	//cout << "stat2b" << endl;
	DoubleMatrix logPi2_t = DoubleMatrix(numOfClasses, numOfClasses);
	logTransposedMatrix(pi2, logPi2_t);
	updateTau(new_tau, stat2b, tau, logPi2_t, temp);

	//cout << "stat3a" << endl;
	DoubleMatrix logPi3 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi3, logPi3);
	updateTau(new_tau, stat3a, tau, logPi3, temp);

	//cout << "stat3b" << endl;
	DoubleMatrix logPi3_t = DoubleMatrix(numOfClasses, numOfClasses);
	logTransposedMatrix(pi3, logPi3_t);
	updateTau(new_tau, stat3b, tau, logPi3_t, temp);

	//cout << "stat4m" << endl;
	DoubleMatrix logPi4m = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi4m, logPi4m);
	updateTau(new_tau, stat4m, tau, logPi4m, temp);

	//cout << "stat4p" << endl;
	DoubleMatrix logPi4p = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi4p, logPi4p);
	updateTau(new_tau, stat4p, tau, logPi4p, temp);

	//cout << "stat5" << endl;
	DoubleMatrix logPi5 = DoubleMatrix(numOfClasses, numOfClasses);
	logMatrix(pi5, logPi5);
	updateTauByNegativeReflection(new_tau, stat5, tau, logPi5, temp);

	normalizeLogTau2Tau(new_tau, minTau);

	prevTau.swap(tau);
	tau.swap(new_tau);
}

bool SignedReciprocityModel::isTauSignificantlyChanged(double _tauPrecision) {

	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			if (fabs(tau(i, k) - prevTau(i, k)) > _tauPrecision)
				return true;

	return false;
}

double SignedReciprocityModel::getLargestTauChange() {
	double largestChange = 0;
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			if (fabs(tau(i, k) - prevTau(i, k)) > largestChange)
				largestChange = fabs(tau(i, k) - prevTau(i, k));
	return largestChange;
}

void SignedReciprocityModel::updatePi(DoubleMatrix& pi,
		NetworkSparseMatrix& stat, const DoubleMatrix& tau,
		const DoubleMatrix& sumTaus) {

	pi.clear();
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

	//	cout << "pi: " << endl;
	//	for (int k = 0; k < numOfClasses; k++) {
	//		for (int l = 0; l < numOfClasses; l++) {
	//			cout << pi(k, l) << "\t";
	//		}
	//		cout << endl;
	//	}
	//
	//	cout << "sumTaus: " << endl;
	//	for (int k = 0; k < numOfClasses; k++) {
	//		for (int l = 0; l < numOfClasses; l++) {
	//			cout << sumTaus(k, l) << "\t";
	//		}
	//		cout << endl;
	//	}

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			pi(k, l) /= sumTaus(k, l);

}

void SignedReciprocityModel::runModelEstimationMStep() {

	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	//printDoubleVector("tauL: ", tauL);

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
	//printDoubleMatrix("sumTaus: ", sumTaus);

	//cout << "Update pi1" << endl;
	prevPi1.swap(pi1);
	updatePi(pi1, stat1a, tau, sumTaus);

	//cout << "Update pi2" << endl;
	prevPi2.swap(pi2);
	updatePi(pi2, stat2a, tau, sumTaus);

	//cout << "Update pi3" << endl;
	prevPi3.swap(pi3);
	updatePi(pi3, stat3a, tau, sumTaus);

	//cout << "Update pi4m" << endl;
	prevPi4m.swap(pi4m);
	updatePi(pi4m, stat4m, tau, sumTaus);

	//cout << "Update pi4p" << endl;
	prevPi4p.swap(pi4p);
	updatePi(pi4p, stat4p, tau, sumTaus);

	//cout << "Update pi5" << endl;
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l <= k; l++) {
			pi5(k, l) = 1 - pi1(k, l) - pi1(l, k) - pi2(k, l) - pi2(l, k)
					- pi3(k, l) - pi3(l, k) - pi4m(k, l) - pi4p(k, l);
			pi5(l, k) = pi5(k, l);
		}
	// check min and normalize again
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			if (pi1(k, l) < minPi)
				pi1(k, l) = minPi;
			if (pi2(k, l) < minPi)
				pi2(k, l) = minPi;
			if (pi3(k, l) < minPi)
				pi3(k, l) = minPi;
			if (pi4m(k, l) < minPi)
				pi4m(k, l) = minPi;
			if (pi4p(k, l) < minPi)
				pi4p(k, l) = minPi;
			if (pi5(k, l) < minPi)
				pi5(k, l) = minPi;
		}
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l <= k; l++) {
			double norms = pi5(k, l) + pi1(k, l) + pi1(l, k) + pi2(k, l) + pi2(
					l, k) + pi3(k, l) + pi3(l, k) + pi4m(k, l) + pi4p(k, l);
			pi1(k, l) = pi1(k, l) / norms;
			pi1(l, k) = pi1(l, k) / norms;
			pi2(k, l) = pi2(k, l) / norms;
			pi2(l, k) = pi2(l, k) / norms;
			pi3(k, l) = pi3(k, l) / norms;
			pi3(l, k) = pi3(l, k) / norms;
			pi4m(k, l) = pi4m(k, l) / norms;
			pi4m(l, k) = pi4m(k, l);
			pi4p(k, l) = pi4p(k, l) / norms;
			pi4p(l, k) = pi4p(k, l);
			pi5(k, l) = pi5(k, l) / norms;
			pi5(l, k) = pi5(k, l);
		}

	//printDoubleMatrix("pi4m: ", pi4m);
	//printDoubleMatrix("pi4p: ", pi4p);
}

bool SignedReciprocityModel::areAlphaPiSignificantlyChanged(
		double _paramPrecision) {

	for (int k = 0; k < numOfClasses; k++)
		if (fabs(alpha(k) - prevAlpha(k)) > _paramPrecision)
			return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi1(k, l) - prevPi1(k, l)) > _paramPrecision)
				return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi2(k, l) - prevPi2(k, l)) > _paramPrecision)
				return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi3(k, l) - prevPi3(k, l)) > _paramPrecision)
				return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi4m(k, l) - prevPi4m(k, l)) > _paramPrecision)
				return true;

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++)
			if (fabs(pi4p(k, l) - prevPi4p(k, l)) > _paramPrecision)
				return true;

	return false;
}

double SignedReciprocityModel::getCompleteLogLikelihood() {

	double CLL = 0;

	// tau * [ log(alpha) - log(tau) ]
	for (int i = 0; i < numOfVertices; i++)
		for (int k = 0; k < numOfClasses; k++)
			CLL += tau(i, k) * log(alpha(k) / tau(i, k));
	//cout << "tau * [ log(alpha) - log(tau) ] " << CLL << endl;

	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;

	// y_ij = -1 and y_ji = 0
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat1a.begin1(); i1 != stat1a.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += temp * log(pi1(k, l) / pi5(k, l));
		}

	// y_ij = 1 and y_ji = 0
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat2a.begin1(); i1 != stat2a.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += temp * log(pi2(k, l) / pi5(k, l));
		}

	// y_ij = -1 and y_ji = 1
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat3a.begin1(); i1 != stat3a.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += temp * log(pi3(k, l) / pi5(k, l));
		}

	// y_ij = -1 and y_ji = -1
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat4m.begin1(); i1 != stat4m.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += .5 * temp * log(pi4m(k, l) / pi5(k, l));
		}

	// y_ij = 1 and y_ji = 1
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (i1_t i1 = stat4p.begin1(); i1 != stat4p.end1(); ++i1) {
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2) {
					unsigned int i = i2.index1();
					unsigned int j = i2.index2();
					temp += tau(i, k) * tau(j, l);
				}
			}
			CLL += .5 * temp * log(pi4p(k, l) / pi5(k, l));
		}

	// norm term
	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			double temp = 0;
			for (int i = 0; i < numOfVertices; i++)
				temp += tau(i, k) * (tauL(l) - tau(i, l));
			CLL += .5 * temp * log(pi5(k, l));
		}
	//cout << "norm term " << CLL << endl;

	return CLL;
}

double SignedReciprocityModel::getICL() {
	//	cout << "SignedReciprocityModel::getICL()" << endl;
	double CLL = getCompleteLogLikelihood();
	//	cout << CLL << endl;
	//	cout << log(numOfVertices) << endl;
	//	cout << numOfVertices << endl;
	//	cout << log((numOfVertices * (numOfVertices - 1.0)) / 2.0) << endl;
	double ICL = CLL - 0.5 * (numOfClasses - 1.0) * log(numOfVertices) - 0.5
			* numOfClasses * (1.0 + 4 * numOfClasses) * log(
			(numOfVertices * (numOfVertices - 1.0)) / 2.0);
	return ICL;
}

void SignedReciprocityModel::copyParameters(DoubleMatrix& bestTau,
		DoubleVector& bestAlpha, VectorOfDoubleMatrices& bestPIs) {
	bestTau.assign(tau);
	bestAlpha.assign(alpha);
	bestPIs(0).assign(pi1);
	bestPIs(1).assign(pi2);
	bestPIs(2).assign(pi3);
	bestPIs(3).assign(pi4m);
	bestPIs(4).assign(pi4p);
	bestPIs(5).assign(pi5);
}

void SignedReciprocityModel::copyPreviousParameters(DoubleMatrix& _prevTau,
		DoubleVector& _prevAlpha, VectorOfDoubleMatrices& _prevPIs) {
	_prevTau.assign(prevTau);
	_prevAlpha.assign(prevAlpha);
	_prevPIs(0).assign(prevPi1);
	_prevPIs(1).assign(prevPi2);
	_prevPIs(2).assign(prevPi3);
	_prevPIs(3).assign(prevPi4m);
	_prevPIs(4).assign(prevPi4p);
}

void SignedReciprocityModel::printPIs() {
	cout.precision(6);
	cout << "pi1: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi1(k, l) << "  ";
		cout << endl;
	}
	cout << "pi2: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi2(k, l) << "  ";
		cout << endl;
	}
	cout << "pi3: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi3(k, l) << "  ";
		cout << endl;
	}
	cout << "pi4m: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi4m(k, l) << "  ";
		cout << endl;
	}
	cout << "pi4p: " << endl;
	for (int k = 0; k < numOfClasses; k++) {
		for (int l = 0; l < numOfClasses; l++)
			cout << pi4p(k, l) << "  ";
		cout << endl;
	}
}

void SignedReciprocityModel::printTau() {
	cout << "Printing TAU: " << endl;
	cout.precision(6);
	for (int i = 0; i < numOfVertices; i++) {
		cout << i << ": ";
		for (int k = 0; k < numOfClasses; k++)
			cout << tau(i, k) << "  ";
		cout << endl;
	}
}

void SignedReciprocityModel::printAlpha() {
	cout.precision(6);
	cout << "Printing ALPHA: " << endl;
	cout << alpha(0);
	for (int k = 1; k < numOfClasses; k++)
		cout << ",  " << alpha(k);
	cout << endl;
}

}
