/*
 * EfficientBinaryReciprocityParametricBoostrapSampler.cpp
 *
 *  Created on: Oct 18, 2011
 *      Author: duyvu
 */

#include "EfficientBinaryReciprocityParametricBoostrapSampler.h"

#define MATHLIB_STANDALONE
#include <Rmath.h>

namespace mixrerg {

EfficientBinaryReciprocityParametricBoostrapSampler::EfficientBinaryReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated constructor stub

}

EfficientBinaryReciprocityParametricBoostrapSampler::~EfficientBinaryReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated destructor stub
}

void EfficientBinaryReciprocityParametricBoostrapSampler::normalizeDyadProbs() {
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++) {
			double sum = 0;
			for (unsigned int t = 0; t < 3; t++)
				sum += dyadProbs(k, l)[t];
			for (unsigned int t = 0; t < 3; t++)
				dyadProbs(k, l)[t] /= sum;
		}
}

int EfficientBinaryReciprocityParametricBoostrapSampler::generateDyadicValue(
		unsigned int k, unsigned int l, int *& dyad_ij,
		SparseGVOCoordinateV & Y, unsigned int i, unsigned int j) {
	// draw a random dyad value
	rmultinom(1, dyadProbs(k, l), 3, dyad_ij);
	int selectedDyadValue = 0;
	while (dyad_ij[selectedDyadValue] == 0)
		selectedDyadValue++;

	// y_ij = [1,0,1,0];
	// y_ji = [0,1,1,0];
	switch (selectedDyadValue) {
	// do not set 0 value because the matrix is in the sparse format
	// if we set value 0, it causes troubles when we iterate through nonzero elements
	// even being setting a value 0, an element turns to be nonzero
	case 0:
		Y(i, j) = 1;
		break;
	case 1:
		Y(j, i) = 1;
		break;
	case 2:
		Y(i, j) = 1;
		Y(j, i) = 1;
		break;
	default:
		cout << "ERROR!!! Undefined dyad value!" << endl;
	}

	return selectedDyadValue;
}

void EfficientBinaryReciprocityParametricBoostrapSampler::generateSelfBlock(
		unsigned int k, SparseGVOCoordinateV& Y, unsigned int kStartIndex,
		unsigned int kN) {

	int * dyad_ij = new int[3];

	// draw the number of non-zero dyads
	long long size = ((long long) kN) * (kN - 1) / 2;
	unsigned int M = rbinom(size, 1 - dyadProbs(k, k)[3]);

	// randomly select each dyad then generate its value
	for (unsigned int edgeIndex = 0; edgeIndex < M; edgeIndex++) {
		unsigned int i;
		unsigned int j;
		// randomly select dyad ij
		do {
			i = kStartIndex + floor(runif(0, 1) * kN);
			j = kStartIndex + floor(runif(0, 1) * kN);
		} while (i == j || Y(i, j) != 0 || Y(j, i) != 0);
		// draw value for the selected dyad ij
		generateDyadicValue(k, k, dyad_ij, Y, i, j);
	}

	// clean up
	if (dyad_ij != NULL)
		delete dyad_ij;
}

void EfficientBinaryReciprocityParametricBoostrapSampler::generateCrossBlock(
		unsigned int k, unsigned int l, SparseGVOCoordinateV & Y,
		unsigned int kStartIndex, unsigned int kN, unsigned int lStartIndex,
		unsigned int lN) {

	int *dyad_ij = new int[3];

	// draw the number of non-zero dyads
	long long size = ((long long) kN) * lN;
	unsigned int M = rbinom(size, 1 - dyadProbs(k, l)[3]);

	// randomly select each dyad then generate its value
	for (unsigned int edgeIndex = 0; edgeIndex < M; edgeIndex++) {
		unsigned int i;
		unsigned int j;
		// randomly select dyad ij
		do {
			i = kStartIndex + floor(runif(0, 1) * kN);
			j = lStartIndex + floor(runif(0, 1) * lN);
		} while (Y(i, j) != 0 || Y(j, i) != 0);
		// draw value for the selected dyad ij
		generateDyadicValue(k, l, dyad_ij, Y, i, j);
	}

	// clean up
	if (dyad_ij != NULL)
		delete dyad_ij;
}

void EfficientBinaryReciprocityParametricBoostrapSampler::generateSample(
		NetworkSparseMatrix& Y, UnsignedIntVector& Z) {

	cout
			<< "EfficientBinaryReciprocityParametricBoostrapSampler::generateSample"
			<< endl;

	// copy alpha to a pointer-based array
	double * mixingProbs = new double[numOfClasses];
	for (unsigned int k = 0; k < alpha.size(); k++)
		mixingProbs[k] = alpha(k);

	// generate the latent class membership
	Z.clear();
	int * z_i = new int[numOfClasses];
	rmultinom(numOfVertices, mixingProbs, numOfClasses, z_i);
	UnsignedIntVector startIndices = UnsignedIntVector(numOfClasses);
	int startIndex = 0;
	for (unsigned int k = 0; k < numOfClasses; k++) {
		startIndices(k) = startIndex;
		unsigned int endIndex = startIndex + z_i[k];
		for (unsigned int i = startIndex; i < endIndex; i++)
			Z(i) = k;
		startIndex = endIndex;
	}

	// generate the adjacency matrix
	SparseGVOCoordinateV tempY = SparseGVOCoordinateV(numOfVertices,
			numOfVertices);
	for (unsigned int k = 0; k < numOfClasses; k++) {
		for (unsigned int l = 0; l < k; l++) {
			//cout << "generating cross block (" << k << ", " << l << ")" << endl;
			generateCrossBlock(k, l, tempY, startIndices(k), z_i[k],
					startIndices(l), z_i[l]);
		}
		//cout << "generating self block (" << k << ", " << k << ")" << endl;
		generateSelfBlock(k, tempY, startIndices(k), z_i[k]);
	}
	Y.assign(tempY);

	unsigned int edgeCount = 0;
	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = Y.begin1(); i1 != Y.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			edgeCount += *i2;
	}
	cout << "Number of nodes is " << numOfVertices << endl;
	cout << "Number of edges is " << edgeCount << endl;

	// clean up some pointer-based vectors and matrices
	if (mixingProbs != NULL)
		delete mixingProbs;
	if (z_i != NULL)
		delete z_i;
}

}
