/*
 * Configuration.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef BINARYCONFIGURATION_H_
#define BINARYCONFIGURATION_H_

#include "Configuration.h"

/**
 * @brief Prints a configuration in binary format.
 */

class BinaryConfiguration: public Configuration {
protected:
	virtual std::string _headers(llint step);
	virtual std::string _configuration(llint step);

public:
	BinaryConfiguration();
	virtual ~BinaryConfiguration();
};

#endif /* BINARYCONFIGURATION_H_ */
