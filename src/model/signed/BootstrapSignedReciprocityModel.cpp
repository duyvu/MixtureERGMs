/*
 * BootstrapSignedReciprocityModel.cpp
 *
 *  Created on: Oct 25, 2010
 *      Author: duyvu
 */

#include "BootstrapSignedReciprocityModel.h"

using namespace std;

namespace mixrerg {

BootstrapSignedReciprocityModel::BootstrapSignedReciprocityModel() {
	// TODO Auto-generated constructor stub

}

BootstrapSignedReciprocityModel::~BootstrapSignedReciprocityModel() {
	// TODO Auto-generated destructor stub
}

void BootstrapSignedReciprocityModel::resetParameters(
		const UnsignedIntVector& Z) {

	cout << "Resetting tau using Z" << endl;
	tau.clear();
	for (int i = 0; i < numOfVertices; i++)
		tau(i, Z(i)) = 1;

	cout << "Resetting model parameters" << endl;
	alpha.clear();
	pi1.clear();
	pi2.clear();
	pi3.clear();
	pi4m.clear();
	pi4p.clear();
	pi5.clear();

	DoubleVector tauL = DoubleVector(numOfClasses);
	sumDoubleMatrixByRow(tau, tauL);
	//printDoubleVector("tauL: ", tauL);

	prevAlpha.swap(alpha);
	alpha.assign(tauL);
	for (unsigned int k = 0; k < alpha.size(); k++)
		alpha(k) /= numOfVertices;
	//printDoubleVector("alpha: ", alpha);

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
}

bool BootstrapSignedReciprocityModel::isTauSignificantlyChanged(
		double _tauPrecision) {
	false;
}

void BootstrapSignedReciprocityModel::runModelEstimationMStep() {
	return;
}

bool BootstrapSignedReciprocityModel::areAlphaPiSignificantlyChanged(
		double _paramPrecision) {
	return false;
}

}
