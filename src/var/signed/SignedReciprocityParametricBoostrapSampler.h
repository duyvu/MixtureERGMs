/*
 * SignedReciprocityParametricBoostrapSampler.h
 *
 *  Created on: Oct 4, 2010
 *      Author: duyvu
 */

#ifndef SIGNEDRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_
#define SIGNEDRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_

#include "var/ParametricBoostrapSampler.h"

namespace mixrerg {

class SignedReciprocityParametricBoostrapSampler: public mixrerg::ParametricBoostrapSampler {
public:
	SignedReciprocityParametricBoostrapSampler();
	virtual ~SignedReciprocityParametricBoostrapSampler();

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

#endif /* SIGNEDRECIPROCITYPARAMETRICBOOSTRAPSAMPLER_H_ */
