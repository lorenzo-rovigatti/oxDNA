/*
 * MetalListFactory.h
 *
 *  Created for Metal backend
 */

#ifndef METALLISTFACTORY_H_
#define METALLISTFACTORY_H_

#include "MetalBaseList.h"

/**
 * @brief Static factory class for Metal lists
 */
class MetalListFactory {
public:
	MetalListFactory() = delete;
	virtual ~MetalListFactory() = delete;

	static MetalBaseList *make_list(input_file &inp);
};

#endif /* METALLISTFACTORY_H_ */
