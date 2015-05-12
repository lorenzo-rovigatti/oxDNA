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
 * @brief Static factory class. Its only public method builds a {@link BaseBox list}.
 *
 * @verbatim
[list_type = cubic (Type of simulation box for CPU simulations.)]
@endverbatim
 */
class BoxFactory {
private:
	BoxFactory();
public:
	virtual ~BoxFactory();

	/**
	 * @brief Builds the box instance.
	 *
	 * @param inp
	 * @param N
	 * @param box
	 * @return a pointer to the newly built box
	 */
	template<typename number>
	static BaseBox<number> *make_box(input_file &inp);
};

#endif /* BOXFACTORY_H_ */
