/*
 * BinaryReciprocityParametricBoostrapSampler.h
 *
 *  Created on: Oct 18, 2011
 *      Author: duyvu
 */

#ifndef BINARYRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_
#define BINARYRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_

#include "var/ParametricBoostrapSampler.h"

namespace mixrerg {

class BinaryReciprocityParametricBoostrapSampler: public mixrerg::ParametricBoostrapSampler {
public:
	BinaryReciprocityParametricBoostrapSampler();
	virtual ~BinaryReciprocityParametricBoostrapSampler();

	unsigned int getNumOfThetas();

	void setParameters(unsigned int _numOfVertices, unsigned int _numOfClasses,
			const DoubleVector& _alpha, const VectorOfDoubleMatrices& _PIs);

	void generateSample(NetworkSparseMatrix& Y, UnsignedIntVector& Z);

	void transformParameters(const VectorOfDoubleMatrices& PIs,
			VectorOfDoubleMatrices& Thetas);

	unsigned int getNumberOfClusteringCoefficients();
	void computeOverallClusteringCoefficients(
			DoubleVector& overallClusteringCoefficientSample,
			const DoubleVector& alpha, const VectorOfDoubleMatrices& PIs);
};

}

#endif /* BINARYRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_ */
