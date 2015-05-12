/*
 * ListFactory.h
 *
 *  Created on: 05/nov/2013
 *      Author: lorenzo
 */

#ifndef LISTFACTORY_H_
#define LISTFACTORY_H_

#include "BaseList.h"

/**
 * @brief Static factory class. Its only public method builds a {@link BaseList list}.
 *
 * @verbatim
[list_type = verlet|cells|no (Type of neighbouring list to be used in CPU simulations. 'no' implies a O(N^2) computational complexity. Defaults to verlet.)]
@endverbatim
 */
class ListFactory {
private:
	ListFactory();
public:
	virtual ~ListFactory();

	/**
	 * @brief Builds the list instance.
	 *
	 * @param inp
	 * @param N
	 * @param box
	 * @return a pointer to the newly built list
	 */
	template<typename number>
	static BaseList<number> *make_list(input_file &inp, int &N, BaseBox<number> *box);
};

#endif /* LISTFACTORY_H_ */
