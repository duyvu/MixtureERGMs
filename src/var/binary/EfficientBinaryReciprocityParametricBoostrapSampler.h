/*
 * EfficientBinaryReciprocityParametricBoostrapSampler.h
 *
 *  Created on: Oct 18, 2011
 *      Author: duyvu
 */

#ifndef EFFICIENTBINARYRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_
#define EFFICIENTBINARYRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_

#include "BinaryReciprocityParametricBoostrapSampler.h"

namespace mixrerg {

class EfficientBinaryReciprocityParametricBoostrapSampler: public mixrerg::BinaryReciprocityParametricBoostrapSampler {
protected:
	int generateDyadicValue(unsigned int k, unsigned int l, int *& dyad_ij,
			SparseGVOCoordinateV & Y, unsigned int i, unsigned int j);
	void generateSelfBlock(unsigned int k, SparseGVOCoordinateV& Y,
			unsigned int kStartIndex, unsigned int kN);
	void generateCrossBlock(unsigned int k, unsigned int l,
			SparseGVOCoordinateV& Y, unsigned int kStartIndex, unsigned int kN,
			unsigned int lStartIndex, unsigned int lN);
public:
	EfficientBinaryReciprocityParametricBoostrapSampler();
	virtual ~EfficientBinaryReciprocityParametricBoostrapSampler();

	void normalizeDyadProbs();
	void generateSample(NetworkSparseMatrix& Y, UnsignedIntVector& Z);
};

}

#endif /* EFFICIENTBINARYRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_ */
