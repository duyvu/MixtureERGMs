/*
 * BootstrapSignedReciprocityModel.h
 *
 *  Created on: Oct 25, 2010
 *      Author: duyvu
 */

#ifndef BOOTSTRAPSIGNEDRECIPROCITYMODEL_H_
#define BOOTSTRAPSIGNEDRECIPROCITYMODEL_H_

#include "SignedReciprocityModel.h"

namespace mixrerg {

class BootstrapSignedReciprocityModel: public mixrerg::SignedReciprocityModel {
public:
	BootstrapSignedReciprocityModel();
	virtual ~BootstrapSignedReciprocityModel();

	void resetParameters(const UnsignedIntVector& Z);
	bool isTauSignificantlyChanged(double _tauPrecision);
	void runModelEstimationMStep();
	bool areAlphaPiSignificantlyChanged(double _paramPrecision);
};

}

#endif /* BOOTSTRAPSIGNEDRECIPROCITYMODEL_H_ */
