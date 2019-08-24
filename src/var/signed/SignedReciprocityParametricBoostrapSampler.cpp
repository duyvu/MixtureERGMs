/*
 * SignedReciprocityParametricBoostrapSampler.cpp
 *
 *  Created on: Oct 4, 2010
 *      Author: duyvu
 */

#include "SignedReciprocityParametricBoostrapSampler.h"
#include "model/signed/SignedReciprocityModel.h"
#include "model/signed/MMSignedReciprocityModel.h"
#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>
#include <string>

#define MATHLIB_STANDALONE
#include <Rmath.h>

using namespace std;

namespace mixrerg {

SignedReciprocityParametricBoostrapSampler::SignedReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated constructor stub
}

SignedReciprocityParametricBoostrapSampler::~SignedReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated destructor stub
	for (unsigned int k = 0; k < dyadProbs.size1(); k++)
		for (unsigned int l = 0; l < dyadProbs.size2(); l++)
			if (dyadProbs(k, l) != NULL)
				delete dyadProbs(k, l);
}

unsigned int SignedReciprocityParametricBoostrapSampler::getNumOfThetas() {
	return 5;
}

void SignedReciprocityParametricBoostrapSampler::setParameters(
		unsigned int _numOfVertices, unsigned int _numOfClasses,
		const DoubleVector& _alpha, const VectorOfDoubleMatrices& _PIs) {

	numOfVertices = _numOfVertices;

	assert(_numOfClasses == _alpha.size());
	numOfClasses = _numOfClasses;
	alpha = DoubleVector(numOfClasses);
	alpha.assign(_alpha);

	//cout << "Done with alpha" << endl;

	assert(_PIs.size() == 6);
	numOfPiParameterSets = 6;
	PIs = VectorOfDoubleMatrices(6);
	for (unsigned int m = 0; m < 6; m++)
		PIs(m) = DoubleMatrix(numOfClasses, numOfClasses);
	PIs(0).assign(_PIs(0));
	PIs(1).assign(_PIs(1));
	PIs(2).assign(_PIs(2));
	PIs(3).assign(_PIs(3));
	PIs(4).assign(_PIs(4));
	PIs(5).assign(_PIs(5));

	//cout << "Done with PIs" << endl;

	dyadProbs = DoublePointerMatrix(numOfClasses, numOfClasses);
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++) {
			dyadProbs(k, l) = new double[9];
			// y_ij = [-1,0,1,0,-1,1,-1,1,0];
			// y_ji = [0,-1,0,1,1,-1,-1,1,0];
			// p0(k,l),p0(l,k),p1(k,l),p1(l,k),p2(k,l),p2(l,k),p3(k,l),p4(k,l),p5(k,l)
			dyadProbs(k, l)[0] = PIs(0)(k, l);
			dyadProbs(k, l)[1] = PIs(0)(l, k);
			dyadProbs(k, l)[2] = PIs(1)(k, l);
			dyadProbs(k, l)[3] = PIs(1)(l, k);
			dyadProbs(k, l)[4] = PIs(2)(k, l);
			dyadProbs(k, l)[5] = PIs(2)(l, k);
			dyadProbs(k, l)[6] = PIs(3)(k, l);
			dyadProbs(k, l)[7] = PIs(4)(k, l);
			dyadProbs(k, l)[8] = PIs(5)(k, l);
			// should validate sum of these probabilities = 1 here
		}
	normalizeDyadProbs();
	//cout << "Done with dyadProbs" << endl;
}

void SignedReciprocityParametricBoostrapSampler::generateSample(
		NetworkSparseMatrix& Y, UnsignedIntVector& Z) {

	// copy alpha to an pointer-based array
	double * mixingProbs = new double[numOfClasses];
	for (unsigned int k = 0; k < alpha.size(); k++)
		mixingProbs[k] = alpha(k);

	// generate the latent class membership
	Z.clear();
	int * z_i = new int[numOfClasses];
	for (unsigned int i = 0; i < numOfVertices; i++) {
		rmultinom(1, mixingProbs, numOfClasses, z_i);
		int selectedClass = 0;
		while (z_i[selectedClass] == 0)
			selectedClass++;
		Z(i) = selectedClass;
	}

	// generate the adjacency matrix
	SparseGVOCoordinateV tempY = SparseGVOCoordinateV(numOfVertices,
			numOfVertices);
	int * dyad_ij = new int[9];
	for (unsigned int i = 0; i < numOfVertices; i++)
		for (unsigned int j = 0; j < i; j++) {
			rmultinom(1, dyadProbs(Z(i), Z(j)), 9, dyad_ij);
			int selectedDyadValue = 0;
			while (dyad_ij[selectedDyadValue] == 0)
				selectedDyadValue++;
			// y_ij = [-1,0,1,0,-1,1,-1,1,0];
			// y_ji = [0,-1,0,1,1,-1,-1,1,0];
			switch (selectedDyadValue) {
			// do not set 0 value because the matrix is in the sparse format
			// if we set value 0, it causes troubles when we iterate through nonzero elements
			// even being setting a value 0, an element turns to be nonzero ;-(
			case 0:
				tempY(i, j) = -1;
				break;
			case 1:
				tempY(j, i) = -1;
				break;
			case 2:
				tempY(i, j) = 1;
				break;
			case 3:
				tempY(j, i) = 1;
				break;
			case 4:
				tempY(i, j) = -1;
				tempY(j, i) = 1;
				break;
			case 5:
				tempY(i, j) = 1;
				tempY(j, i) = -1;
				break;
			case 6:
				tempY(i, j) = -1;
				tempY(j, i) = -1;
				break;
			case 7:
				tempY(i, j) = 1;
				tempY(j, i) = 1;
				break;
			case 8:
				break;
			default:
				cout << "ERROR!!! Undefined dyad value!" << endl;
			}
		}
	Y.assign(tempY);

	// clean up some pointer-based vectors and matrices
	if (mixingProbs != NULL)
		delete mixingProbs;
	if (z_i != NULL)
		delete z_i;
	if (dyad_ij != NULL)
		delete dyad_ij;
}

void SignedReciprocityParametricBoostrapSampler::transformParameters(
		const VectorOfDoubleMatrices& PIs, VectorOfDoubleMatrices& Thetas) {
	Thetas.resize(5);
	Thetas(0) = DoubleMatrix(numOfClasses, numOfClasses);
	Thetas(1) = DoubleMatrix(numOfClasses, numOfClasses);
	Thetas(2) = DoubleMatrix(numOfClasses, numOfClasses);
	Thetas(3) = DoubleMatrix(numOfClasses, numOfClasses);
	Thetas(4) = DoubleMatrix(numOfClasses, numOfClasses);
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++) {
			Thetas(0)(k, l) = log(PIs(0)(k, l) / PIs(5)(k, l));
			Thetas(1)(k, l) = log(PIs(1)(k, l) / PIs(5)(k, l));
			Thetas(2)(k, l) = log((PIs(2)(k, l) * PIs(5)(k, l)) / (PIs(0)(k, l)
					* PIs(1)(l, k)));
			Thetas(3)(k, l) = log((PIs(3)(k, l) * PIs(5)(k, l)) / (PIs(0)(k, l)
					* PIs(0)(l, k)));
			Thetas(4)(k, l) = log((PIs(4)(k, l) * PIs(5)(k, l)) / (PIs(1)(k, l)
					* PIs(1)(l, k)));
		}
}

unsigned int SignedReciprocityParametricBoostrapSampler::getNumberOfClusteringCoefficients() {
	return 5;
}

void SignedReciprocityParametricBoostrapSampler::computeOverallClusteringCoefficients(
		DoubleVector& overallClusteringCoefficientSample,
		const DoubleVector& alpha, const VectorOfDoubleMatrices& PIs) {

	double C_fff = 0;
	double p_fff = 0;
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++)
			for (unsigned int m = 0; m < numOfClasses; m++) {
				C_fff = C_fff + alpha(k) * alpha(l) * alpha(m) * (PIs(1)(l, m)
						+ PIs(2)(m, l) + PIs(4)(l, m)) * (PIs(1)(k, l) + PIs(2)(l, k)
						+ PIs(4)(k, l)) * (PIs(1)(k, m) + PIs(2)(m, k) + PIs(4)(k, m));
				p_fff = p_fff + alpha(k) * alpha(l) * alpha(m) * (PIs(1)(l, m)
						+ PIs(2)(m, l) + PIs(4)(l, m)) * (PIs(1)(k, l) + PIs(2)(l, k)
						+ PIs(4)(k, l));
			}
	C_fff = C_fff / p_fff;
	overallClusteringCoefficientSample(0) = C_fff;

	double C_eef = 0;
	double p_eef = 0;
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++)
			for (unsigned int m = 0; m < numOfClasses; m++) {
				C_eef = C_eef + alpha(k) * alpha(l) * alpha(m) * (PIs(0)(l, m)
						+ PIs(2)(l, m) + PIs(3)(l, m)) * (PIs(0)(k, l) + PIs(2)(k, l)
						+ PIs(3)(k, l)) * (PIs(1)(k, m) + PIs(2)(m, k) + PIs(4)(k, m));
				p_eef = p_eef + alpha(k) * alpha(l) * alpha(m) * (PIs(0)(l, m)
						+ PIs(2)(l, m) + PIs(3)(l, m)) * (PIs(0)(k, l) + PIs(2)(k, l)
						+ PIs(3)(k, l));
			}
	C_eef = C_eef / p_eef;
	overallClusteringCoefficientSample(1) = C_eef;

	double C_fee = 0;
	double p_fee = 0;
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++)
			for (unsigned int m = 0; m < numOfClasses; m++) {
				C_fee = C_fee + alpha(k) * alpha(l) * alpha(m) * (PIs(1)(l, m)
						+ PIs(2)(m, l) + PIs(4)(l, m)) * (PIs(0)(k, l) + PIs(2)(k, l)
						+ PIs(3)(k, l)) * (PIs(0)(k, m) + PIs(2)(k, m) + PIs(3)(k, m));
				p_fee = p_fee + alpha(k) * alpha(l) * alpha(m) * (PIs(1)(l, m)
						+ PIs(2)(m, l) + PIs(4)(l, m)) * (PIs(0)(k, l) + PIs(2)(k, l)
						+ PIs(3)(k, l));
			}
	C_fee = C_fee / p_fee;
	overallClusteringCoefficientSample(2) = C_fee;

	double C_efe = 0;
	double p_efe = 0;
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++)
			for (unsigned int m = 0; m < numOfClasses; m++) {
				C_efe = C_efe + alpha(k) * alpha(l) * alpha(m) * (PIs(0)(l, m)
						+ PIs(2)(l, m) + PIs(3)(l, m)) * (PIs(1)(k, l) + PIs(2)(l, k)
						+ PIs(4)(k, l)) * (PIs(0)(k, m) + PIs(2)(k, m) + PIs(3)(k, m));
				p_efe = p_efe + alpha(k) * alpha(l) * alpha(m) * (PIs(0)(l, m)
						+ PIs(2)(l, m) + PIs(3)(l, m)) * (PIs(1)(k, l) + PIs(2)(l, k)
						+ PIs(4)(k, l));
			}
	C_efe = C_efe / p_efe;
	overallClusteringCoefficientSample(3) = C_efe;

	double C_aggreated = (C_fff * p_fff + C_eef * p_eef + C_fee * p_fee + C_efe
			* p_efe) / (p_fff + p_eef + p_fee + p_efe);
	overallClusteringCoefficientSample(4) = C_aggreated;
}

}
