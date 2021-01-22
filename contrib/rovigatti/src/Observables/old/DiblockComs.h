/*
 * Gyradius.h
 *
 *  Created on: Oct 30, 2013
 *      Author: Lorenzo
 */

#ifndef DIBLOCKCOMS_H_
#define DIBLOCKCOMS_H_

#include "Observables/BaseObservable.h"

/**
 * @brief Computes the A-A, A-B, B-B and intramolecular g(r) of a two-chain diblock copolymer system.
 *
 */

class DiblockComs: public BaseObservable {
protected:
	bool _only_intra;

public:
	DiblockComs();
	virtual ~DiblockComs();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
	virtual std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_DiblockComs() {
	return new DiblockComs();
}

#endif /* DIBLOCKCOMS_H_ */
