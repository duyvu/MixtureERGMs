/*
 * EfficientSignedReciprocityParametricBoostrapSampler.h
 *
 *  Created on: Oct 15, 2010
 *      Author: duyvu
 */

#ifndef EFFICIENTSIGNEDRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_
#define EFFICIENTSIGNEDRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_

#include "SignedReciprocityParametricBoostrapSampler.h"

namespace mixrerg {

class EfficientSignedReciprocityParametricBoostrapSampler: public mixrerg::SignedReciprocityParametricBoostrapSampler {
protected:
	int generateDyadicValue(unsigned int k, unsigned int l, int *& dyad_ij,
			SparseGVOCoordinateV & Y, unsigned int i, unsigned int j);
	void generateSelfBlock(unsigned int k, SparseGVOCoordinateV& Y,
			unsigned int kStartIndex, unsigned int kN);
	void generateCrossBlock(unsigned int k, unsigned int l,
			SparseGVOCoordinateV& Y, unsigned int kStartIndex, unsigned int kN,
			unsigned int lStartIndex, unsigned int lN);
public:
	EfficientSignedReciprocityParametricBoostrapSampler();
	virtual ~EfficientSignedReciprocityParametricBoostrapSampler();

	void normalizeDyadProbs();
	void generateSample(NetworkSparseMatrix& Y, UnsignedIntVector& Z);
};

}

#endif /* EFFICIENTSIGNEDRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_ */
