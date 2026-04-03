/*
 * ObservableFactory.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef OBSERVABLEFACTORY_H_
#define OBSERVABLEFACTORY_H_

#include "BaseObservable.h"

class ObservableOutput;

/**
 * @brief Static factory class. It produces observables.
 *
 * This class can be used to produce BaseObservable instances.
 */
class ObservableFactory {
public:
	ObservableFactory() = delete;
	virtual ~ObservableFactory() = delete;

	/**
	 * @brief Creates an observable given an input file containing its definition
	 *
	 * This method needs two input files as parameters: one is the observable input file. The other
	 * one is the simulation input file. The observable may need it.
	 * @param obs_inp observable input file (it is usually generated from a string)
	 * @param sim_inp simulation input file
	 * @return pointer to the newly created observable
	 */

	static ObservablePtr make_observable(input_file &obs_inp);

	static std::vector<std::shared_ptr<ObservableOutput>> make_observables(std::string prefix="");
};

#endif /* OBSERVABLEFACTORY_H_ */
