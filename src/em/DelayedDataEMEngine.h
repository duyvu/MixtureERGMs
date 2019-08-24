/*
 * DelayedDataEMEngine.h
 *
 *  Created on: Oct 14, 2010
 *      Author: duyvu
 */

#ifndef DELAYEDDATAEMENGINE_H_
#define DELAYEDDATAEMENGINE_H_

#include "EMEngine.h"

namespace mixrerg {

class DelayedDataEMEngine: public mixrerg::EMEngine {
public:
	DelayedDataEMEngine();
	virtual ~DelayedDataEMEngine();

	void configure(const char* configurationFile, const char* delim);
};

}

#endif /* DELAYEDDATAEMENGINE_H_ */
