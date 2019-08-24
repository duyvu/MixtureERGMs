/*
 * EfficientSignedReciprocityCoverageExperiment.h
 *
 *  Created on: Oct 27, 2010
 *      Author: duyvu
 */

#ifndef EFFICIENTSIGNEDRECIPROCITYCOVERAGEEXPERIMENT_H_
#define EFFICIENTSIGNEDRECIPROCITYCOVERAGEEXPERIMENT_H_

#include "SignedReciprocityCoverageExperiment.h"

namespace mixrerg {

class EfficientSignedReciprocityCoverageExperiment: public mixrerg::SignedReciprocityCoverageExperiment {
protected:
	int generateDyadicValue(unsigned int k, unsigned int l, int *& dyad_ij,
			SparseGVOCoordinateV & Y, unsigned int i, unsigned int j);
	void generateSelfBlock(unsigned int k, SparseGVOCoordinateV& Y,
			unsigned int kStartIndex, unsigned int kN);
	void generateCrossBlock(unsigned int k, unsigned int l,
			SparseGVOCoordinateV& Y, unsigned int kStartIndex, unsigned int kN,
			unsigned int lStartIndex, unsigned int lN);
public:
	EfficientSignedReciprocityCoverageExperiment();
	virtual ~EfficientSignedReciprocityCoverageExperiment();

	void normalizeDyadProbs();
	void generateSample(NetworkSparseMatrix& Y, UnsignedIntVector& Z);
};

}

#endif /* EFFICIENTSIGNEDRECIPROCITYCOVERAGEEXPERIMENT_H_ */
