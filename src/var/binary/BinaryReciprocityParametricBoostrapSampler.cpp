/*
 * BinaryReciprocityParametricBoostrapSampler.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: duyvu
 */

#include "BinaryReciprocityParametricBoostrapSampler.h"

#include "model/binary/BinaryReciprocityModel.h"
#include "model/binary/MMBinaryReciprocityModel.h"
#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>
#include <string>

#define MATHLIB_STANDALONE
#include <Rmath.h>

namespace mixrerg {

BinaryReciprocityParametricBoostrapSampler::BinaryReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated constructor stub

}

BinaryReciprocityParametricBoostrapSampler::~BinaryReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated destructor stub
	for (unsigned int k = 0; k < dyadProbs.size1(); k++)
		for (unsigned int l = 0; l < dyadProbs.size2(); l++)
			if (dyadProbs(k, l) != NULL)
				delete dyadProbs(k, l);
}

unsigned int BinaryReciprocityParametricBoostrapSampler::getNumOfThetas() {
	return 2; // 10 and 11
}

void BinaryReciprocityParametricBoostrapSampler::setParameters(
		unsigned int _numOfVertices, unsigned int _numOfClasses,
		const DoubleVector& _alpha, const VectorOfDoubleMatrices& _PIs) {

	numOfVertices = _numOfVertices;

	assert(_numOfClasses == _alpha.size());
	numOfClasses = _numOfClasses;
	alpha = DoubleVector(numOfClasses);
	alpha.assign(_alpha);

	//cout << "Done with alpha" << endl;
	assert(_PIs.size() == 3);
	numOfPiParameterSets = 3;
	PIs = VectorOfDoubleMatrices(3);
	for (unsigned int m = 0; m < 3; m++)
		PIs(m) = DoubleMatrix(numOfClasses, numOfClasses);
	PIs(0).assign(_PIs(0)); // 10
	PIs(1).assign(_PIs(1)); // 11
	PIs(2).assign(_PIs(2)); // 00

	//cout << "Done with PIs" << endl;

	dyadProbs = DoublePointerMatrix(numOfClasses, numOfClasses);
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++) {
			dyadProbs(k, l) = new double[4];
			// y_ij = [1,0,1,0];
			// y_ji = [0,1,1,0];
			// p0(k,l),p0(l,k),p1(k,l),p2(k,l)
			dyadProbs(k, l)[0] = PIs(0)(k, l);
			dyadProbs(k, l)[1] = PIs(0)(l, k);
			dyadProbs(k, l)[2] = PIs(1)(k, l);
			dyadProbs(k, l)[3] = PIs(2)(k, l);
			// should validate sum of these probabilities = 1 here
		}
	normalizeDyadProbs(); // only need for the efficient implementation

	//cout << "Done with dyadProbs" << endl;
}

void BinaryReciprocityParametricBoostrapSampler::generateSample(
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
	int * dyad_ij = new int[4];
	for (unsigned int i = 0; i < numOfVertices; i++)
		for (unsigned int j = 0; j < i; j++) {
			rmultinom(1, dyadProbs(Z(i), Z(j)), 4, dyad_ij);
			int selectedDyadValue = 0;
			while (dyad_ij[selectedDyadValue] == 0)
				selectedDyadValue++;
			// y_ij = [1,0,1,0];
			// y_ji = [0,1,1,0];
			switch (selectedDyadValue) {
			// do not set 0 value because the matrix is in the sparse format
			// if we set value 0, it causes troubles when we iterate through nonzero elements
			// even being setting a value 0, an element turns to be nonzero ;-(
			case 0:
				tempY(i, j) = 1;
				break;
			case 1:
				tempY(j, i) = 1;
				break;
			case 2:
				tempY(i, j) = 1;
				tempY(j, i) = 1;
				break;
			case 3:
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

void BinaryReciprocityParametricBoostrapSampler::transformParameters(
		const VectorOfDoubleMatrices& PIs, VectorOfDoubleMatrices& Thetas) {
	Thetas.resize(2);
	Thetas(0) = DoubleMatrix(numOfClasses, numOfClasses);
	Thetas(1) = DoubleMatrix(numOfClasses, numOfClasses);
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l <= k; l++) {
			Thetas(0)(k, l) = log(PIs(0)(k, l) / PIs(2)(k, l));
			Thetas(0)(l, k) = log(PIs(0)(l, k) / PIs(2)(k, l));
			Thetas(1)(k, l) = log(PIs(1)(k, l) / PIs(2)(k, l));
			Thetas(1)(l, k) = Thetas(1)(k, l);
		}
}

unsigned int BinaryReciprocityParametricBoostrapSampler::getNumberOfClusteringCoefficients() {
	return numOfClasses + 1;
}

void BinaryReciprocityParametricBoostrapSampler::computeOverallClusteringCoefficients(
		DoubleVector& overallClusteringCoefficientSample,
		const DoubleVector& alpha, const VectorOfDoubleMatrices& PIs) {

	for (unsigned int k = 0; k < numOfClasses; k++)
		overallClusteringCoefficientSample(k) = PIs(0)(k, k) + PIs(0)(k, k)
				+ PIs(1)(k, k);

	double triangles = 0;
	double twoPaths = 0;
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++)
			for (unsigned int m = 0; m < numOfClasses; m++) {
				triangles += alpha(k) * alpha(l) * alpha(m) * (PIs(0)(k, l)
						+ PIs(0)(l, k) + PIs(1)(k, l)) * (PIs(0)(l, m)
						+ PIs(0)(m, l) + PIs(1)(l, m)) * (PIs(0)(m, k)
						+ PIs(0)(k, m) + PIs(1)(m, k));
				twoPaths += alpha(k) * alpha(l) * alpha(m) * (PIs(0)(k, l)
						+ PIs(0)(l, k) + PIs(1)(k, l)) * (PIs(0)(l, m)
						+ PIs(0)(m, l) + PIs(1)(l, m));
			}
	overallClusteringCoefficientSample(numOfClasses) = triangles / twoPaths;
}

}
