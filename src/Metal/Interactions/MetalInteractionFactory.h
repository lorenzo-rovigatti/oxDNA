/*
 * MetalInteractionFactory.h
 *
 *  Created for Metal backend
 */

#ifndef METALINTERACTIONFACTORY_H_
#define METALINTERACTIONFACTORY_H_

#include "MetalBaseInteraction.h"

/**
 * @brief Static factory class for Metal interactions
 */
class MetalInteractionFactory {
public:
	MetalInteractionFactory() = delete;
	virtual ~MetalInteractionFactory() = delete;

	static MetalBaseInteraction *make_interaction(input_file &inp);
};

#endif /* METALINTERACTIONFACTORY_H_ */
