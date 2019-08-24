/*
 * BootstrapBinaryReciprocityModel.h
 *
 *  Created on: Oct 6, 2011
 *      Author: duyvu
 */

#ifndef BOOTSTRAPBINARYRECIPROCITYMODEL_H_
#define BOOTSTRAPBINARYRECIPROCITYMODEL_H_

#include "BinaryReciprocityModel.h"

namespace mixrerg {

class BootstrapBinaryReciprocityModel: public mixrerg::BinaryReciprocityModel {
public:
	BootstrapBinaryReciprocityModel();
	virtual ~BootstrapBinaryReciprocityModel();

	void resetParameters(const UnsignedIntVector& Z);
	bool isTauSignificantlyChanged(double _tauPrecision);
	void runModelEstimationMStep();
	bool areAlphaPiSignificantlyChanged(double _paramPrecision);
};

}

#endif /* BOOTSTRAPBINARYRECIPROCITYMODEL_H_ */
