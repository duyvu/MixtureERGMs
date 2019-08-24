/*
 * CoverageExperiment.h
 *
 *  Created on: Oct 13, 2010
 *      Author: duyvu
 */

#ifndef COVERAGEEXPERIMENT_H_
#define COVERAGEEXPERIMENT_H_

#include "em/DelayedDataEMEngine.h"
#include "ParametricBoostrapSampler.h"
#include "signed/SignedReciprocityParametricBoostrapSampler.h"

#include <string>

using namespace std;

namespace mixrerg {

class CoverageExperiment {
protected:
	unsigned int numOfVertices;
	unsigned int numOfClasses;
	DoubleVector alpha;
	unsigned int numOfPiParameterSets;
	VectorOfDoubleMatrices PIs;
	DoublePointerMatrix dyadProbs;

	DelayedDataEMEngine engine;
	ParametricBoostrapSampler* sampler;

	unsigned int coverageSampleSize;

	string emConfigurationFile;
	string bootstrapConfigurationFile;
	string inputDirectory;
	string outputDirectory;
	unsigned int outputPrecision;

	int seed1;
	int seed2;
public:
	CoverageExperiment();
	virtual ~CoverageExperiment();

	virtual void configure(const char *configurationFile, const char *delim)=0;
	virtual void normalizeDyadProbs()=0;
	virtual void
			generateSample(NetworkSparseMatrix & Y, UnsignedIntVector & Z)=0;

	void getCoverageSample(unsigned int coverageSampleIndex);
	void getCoverageSamples(unsigned int fromCoverageSampleIndex,
			unsigned int toCoverageSampleIndex);
	void getCoverageSamples();

	void getCoverageSampleUsingParameterFiles(unsigned int coverageSampleIndex);
	void getCoverageSamplesUsingParameterFiles(
			unsigned int fromCoverageSampleIndex,
			unsigned int toCoverageSampleIndex);
	void getCoverageSamplesUsingParameterFiles();

	int getSeed1() const {
		return seed1;
	}

	void setSeed1(int seed1) {
		this->seed1 = seed1;
	}

	int getSeed2() const {
		return seed2;
	}

	void setSeed2(int seed2) {
		this->seed2 = seed2;
	}

};

}

#endif /* COVERAGEEXPERIMENT_H_ */
