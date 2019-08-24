/*
 * SignedReciprocityCoverageExperiment.h
 *
 *  Created on: Oct 14, 2010
 *      Author: duyvu
 */

#ifndef SIGNEDRECIPROCITYCOVERAGEEXPERIMENT_H_
#define SIGNEDRECIPROCITYCOVERAGEEXPERIMENT_H_

#include "var/CoverageExperiment.h"

namespace mixrerg {

class SignedReciprocityCoverageExperiment: public mixrerg::CoverageExperiment {
public:
	SignedReciprocityCoverageExperiment();
	virtual ~SignedReciprocityCoverageExperiment();

	void configure(const char* configurationFile, const char* delim);
	void normalizeDyadProbs();
	void generateSample(NetworkSparseMatrix& Y, UnsignedIntVector& Z);
};

}

#endif /* SIGNEDRECIPROCITYCOVERAGEEXPERIMENT_H_ */
