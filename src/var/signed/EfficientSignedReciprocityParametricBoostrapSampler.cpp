/*
 * EfficientSignedReciprocityParametricBoostrapSampler.cpp
 *
 *  Created on: Oct 15, 2010
 *      Author: duyvu
 */

#include "EfficientSignedReciprocityParametricBoostrapSampler.h"

#define MATHLIB_STANDALONE
#include <Rmath.h>

namespace mixrerg {

EfficientSignedReciprocityParametricBoostrapSampler::EfficientSignedReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated constructor stub
}

EfficientSignedReciprocityParametricBoostrapSampler::~EfficientSignedReciprocityParametricBoostrapSampler() {
	// TODO Auto-generated destructor stub
}

void EfficientSignedReciprocityParametricBoostrapSampler::normalizeDyadProbs() {
	for (unsigned int k = 0; k < numOfClasses; k++)
		for (unsigned int l = 0; l < numOfClasses; l++) {
			double sum = 0;
			for (unsigned int t = 0; t < 8; t++)
				sum += dyadProbs(k, l)[t];
			for (unsigned int t = 0; t < 8; t++)
				dyadProbs(k, l)[t] /= sum;
		}
}

int EfficientSignedReciprocityParametricBoostrapSampler::generateDyadicValue(
		unsigned int k, unsigned int l, int *& dyad_ij,
		SparseGVOCoordinateV & Y, unsigned int i, unsigned int j) {

	//cout << "generateDyadicValue" << endl;

	rmultinom(1, dyadProbs(k, l), 8, dyad_ij);
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
		Y(i, j) = -1;
		break;
	case 1:
		Y(j, i) = -1;
		break;
	case 2:
		Y(i, j) = 1;
		break;
	case 3:
		Y(j, i) = 1;
		break;
	case 4:
		Y(i, j) = -1;
		Y(j, i) = 1;
		break;
	case 5:
		Y(i, j) = 1;
		Y(j, i) = -1;
		break;
	case 6:
		Y(i, j) = -1;
		Y(j, i) = -1;
		break;
	case 7:
		Y(i, j) = 1;
		Y(j, i) = 1;
		break;
	case 8:
		break;
	default:
		cout << "ERROR!!! Undefined dyad value!" << endl;
	}

	return selectedDyadValue;
}

void EfficientSignedReciprocityParametricBoostrapSampler::generateSelfBlock(
		unsigned int k, SparseGVOCoordinateV& Y, unsigned int kStartIndex,
		unsigned int kN) {

	//DoubleVector counts = DoubleVector(8);
	//counts.clear();

	//cout << "generateSelfBlock: " << k << endl;

	int * dyad_ij = new int[8];

	// draw the number of non-zero dyads
	//cout << "kN = " << kN << endl;
	long long size = ((long long) kN) * (kN - 1) / 2;
	//cout << "size = " << size << endl;
	unsigned int M = rbinom(size, 1 - dyadProbs(k, k)[8]);
	//cout << "M = " << M << endl;

	// randomly select each dyad then generate its value
	//cout << "kStartIndex = " << kStartIndex << endl;
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
		//int dyadicValue = generateDyadicValue(k, k, dyad_ij, Y, i, j);
		//counts(dyadicValue) += 1;
	}

	//	cout << "dyadic freq: " << endl;
	//	for (unsigned int k = 0; k < 8; k++)
	//		cout << counts(k) << "\t";
	//	cout << endl;

	if (dyad_ij != NULL)
		delete dyad_ij;
}

void EfficientSignedReciprocityParametricBoostrapSampler::generateCrossBlock(
		unsigned int k, unsigned int l, SparseGVOCoordinateV & Y,
		unsigned int kStartIndex, unsigned int kN, unsigned int lStartIndex,
		unsigned int lN) {

	//cout << "generateCrossBlock: " << k << "\t" << l << endl;

	int *dyad_ij = new int[8];

	// draw the number of non-zero dyads
	//cout << "kN = " << kN << endl;
	//cout << "lN = " << lN << endl;
	long long size = ((long long) kN) * lN;
	//cout << "size = " << size << endl;
	unsigned int M = rbinom(size, 1 - dyadProbs(k, l)[8]);
	//cout << "M = " << M << endl;

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

	if (dyad_ij != NULL)
		delete dyad_ij;
}

void EfficientSignedReciprocityParametricBoostrapSampler::generateSample(
		NetworkSparseMatrix& Y, UnsignedIntVector& Z) {

	cout
			<< "EfficientSignedReciprocityParametricBoostrapSampler::generateSample"
			<< endl;

	// copy alpha to an pointer-based array
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

	double edgeCount = 0;
	typedef NetworkSparseMatrix::iterator1 i1_t;
	typedef NetworkSparseMatrix::iterator2 i2_t;
	for (i1_t i1 = Y.begin1(); i1 != Y.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			edgeCount += abs(*i2);
	}
	cout << "Number of edges is " << edgeCount << endl;

	// clean up some pointer-based vectors and matrices
	if (mixingProbs != NULL)
		delete mixingProbs;
	if (z_i != NULL)
		delete z_i;
}

}
