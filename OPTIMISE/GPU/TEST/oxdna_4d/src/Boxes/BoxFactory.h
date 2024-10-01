/*
 * BoxFactory.h
 *
 *  Created on: 17/mar/2013
 *      Author: lorenzo
 */

#ifndef BOXFACTORY_H_
#define BOXFACTORY_H_

#include "BaseBox.h"

/**
 * @brief Static factory class. Its only public method builds a {@link BaseBox}.
 *
 * @verbatim
 [list_type = cubic (Type of simulation box for CPU simulations.)]
 @endverbatim
 */
class BoxFactory {
public:
	BoxFactory() = delete;
	virtual ~BoxFactory() = delete;

	/**
	 * @brief Builds the box instance.
	 *
	 * @param inp
	 * @return a pointer to the newly built box
	 */

	static BoxPtr make_box(input_file &inp);
};

#endif /* BOXFACTORY_H_ */
