/*
 * PatchyBonds.h
 *
 *  Created on: 16/aug/2020
 *      Author: lorenzo
 */

#ifndef PATCHYBONDS_H_
#define PATCHYBONDS_H_

#include "Observables/Configurations/Configuration.h"

/**
 * @brief Output the bonding pattern of patchy systems.
 */
class PatchyBonds: public Configuration {
protected:
	bool _zero_based_indexes = false;

	virtual std::string _headers(llint step);

public:
	PatchyBonds();
	virtual ~PatchyBonds();

	void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string _configuration(llint step);
};

extern "C" BaseObservable *make_PatchyBonds() {
	return new PatchyBonds();
}

#endif /* PATCHYBONDS_H_ */
