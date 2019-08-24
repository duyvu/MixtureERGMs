/*
 * MMBinaryReciprocityModel.h
 *
 *  Created on: Oct 6, 2011
 *      Author: duyvu
 */

#ifndef MMBINARYRECIPROCITYMODEL_H_
#define MMBINARYRECIPROCITYMODEL_H_

#include "BinaryReciprocityModel.h"

namespace mixrerg {

class MMBinaryReciprocityModel: public mixrerg::BinaryReciprocityModel {
public:
	MMBinaryReciprocityModel();
	virtual ~MMBinaryReciprocityModel();

	void runFixedPointEstimationEStep();
	bool isTauSignificantlyChanged(double _tauPrecision);
};

}
#endif /* MMBINARYRECIPROCITYMODEL_H_ */
