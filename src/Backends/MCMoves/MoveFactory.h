/*
 * MoveFactory.h
 *
 *  Created on: Jun 25, 2014
 *      Author: flavio
 */

#ifndef MOVEFACTORY_H_
#define MOVEACTORY_H_

#include "BaseMove.h"

/**
 * @brief Static factory class. Its only public method builds an {@link BaseMove move}.
 *
 * @verbatim
[type = rotation|traslation|possibly other as they get added (move to perform. No Defaults)]
@endverbatim
 */
class MoveFactory {
private:
	MoveFactory();
public:
	virtual ~MoveFactory();


	/**
	 * @brief Builds the interaction instance.
	 *
	 * @param inp input file
	 * @param sim_inp input file of the simulation
	 * @param Info pointer to a ConfingInfo object, which is the one that the move is going to use do do its magic
	 *
	 * @return a pointer to the newly built interaction
	 */
	template <typename number>
	static BaseMove<number> * make_move (input_file &inp, input_file &sim_inp, ConfigInfo<number> * Info);
};

#endif /* INTERACTIONFACTORY_H_ */
