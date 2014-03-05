/*
 * ObservableFactory.h
 *
 *  Created on: Feb 12, 2013
 *      Author: rovigatti
 */

#ifndef OBSERVABLEFACTORY_H_
#define OBSERVABLEFACTORY_H_

#include "BaseObservable.h"
#include "../Utilities/parse_input/parse_input.h"

/**
 * @brief Static factory class. It produces observables.
 *
 * This class can be used to produce BaseObservable instances.
 */
class ObservableFactory {
private:
	/**
	 * @brief The constructor is private because our factories have just one single static method
	 */
	ObservableFactory();
public:
	virtual ~ObservableFactory();

	/**
	 * @brief Creates an observable given an input file containing its definition
	 *
	 * This method needs two input files as parameters: one is the observable input file. The other
	 * one is the simulation input file. The observable may need it.
	 * @param obs_inp observable input file (it is usually generated from a string)
	 * @param sim_inp simulation input file
	 * @return pointer to the newly created observable
	 */
	template<typename number>
	static BaseObservable<number> *make_observable(input_file &obs_inp, input_file &sim_inp);
};

#endif /* OBSERVABLEFACTORY_H_ */
