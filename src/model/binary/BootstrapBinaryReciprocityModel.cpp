/*
 * BootstrapBinaryReciprocityModel.cpp
 *
 *  Created on: Oct 6, 2011
 *      Author: duyvu
 */

#include "BootstrapBinaryReciprocityModel.h"

using namespace std;

namespace mixrerg {

BootstrapBinaryReciprocityModel::BootstrapBinaryReciprocityModel() {
	// TODO Auto-generated constructor stub
}

BootstrapBinaryReciprocityModel::~BootstrapBinaryReciprocityModel() {
	// TODO Auto-generated destructor stub
}

void BootstrapBinaryReciprocityModel::resetParameters(
		const UnsignedIntVector& Z) {

	cout << "Resetting tau using Z" << endl;
	tau.clear();
	for (int i = 0; i < numOfVertices; i++)
		tau(i, Z(i)) = 1;

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi10.clear();
	pi11.clear();
	pi00.clear();

	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);

	prevAlpha.swap(alpha);
	alpha.assign(tauL);
	for (unsigned int k = 0; k < alpha.size(); k++)
		alpha(k) /= numOfVertices;

	DoubleMatrix sumTaus = DoubleMatrix(numOfClasses, numOfClasses);
	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l < numOfClasses; l++) {
			sumTaus(k, l) = 0;
			for (int i = 0; i < numOfVertices; i++)
				sumTaus(k, l) += tau(i, k) * (tauL(l) - tau(i, l));
		}

	//cout << "Update pi10" << endl;
	prevPi10.swap(pi10);
	updatePi(pi10, stat10, tau, sumTaus);

	//cout << "Update pi11" << endl;
	prevPi11.swap(pi11);
	updatePi(pi11, stat11, tau, sumTaus);

	for (int k = 0; k < numOfClasses; k++)
		for (int l = 0; l <= k; l++) {
			pi00(k, l) = 1 - pi10(k, l) - pi10(l, k) - pi11(k, l);
			pi00(l, k) = pi00(k, l);
		}
}

bool BootstrapBinaryReciprocityModel::isTauSignificantlyChanged(
		double _tauPrecision) {
	return false;
}

void BootstrapBinaryReciprocityModel::runModelEstimationMStep() {
	return;
}

bool BootstrapBinaryReciprocityModel::areAlphaPiSignificantlyChanged(
		double _paramPrecision) {
	return false;
}

}
