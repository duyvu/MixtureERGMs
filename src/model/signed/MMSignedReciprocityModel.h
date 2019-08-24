/*
 * MMSignedReciprocityModel.h
 *
 *  Created on: Jan 25, 2010
 *      Author: duyvu
 */

#ifndef MMSIGNEDRECIPROCITYMODEL_H_
#define MMSIGNEDRECIPROCITYMODEL_H_

#include "SignedReciprocityModel.h"

namespace mixrerg {

class MMSignedReciprocityModel: public mixrerg::SignedReciprocityModel {
public:
	MMSignedReciprocityModel();
	virtual ~MMSignedReciprocityModel();

	void runFixedPointEstimationEStep();
	bool isTauSignificantlyChanged(double _tauPrecision);
};

}

#endif /* MMSIGNEDRECIPROCITYMODEL_H_ */
