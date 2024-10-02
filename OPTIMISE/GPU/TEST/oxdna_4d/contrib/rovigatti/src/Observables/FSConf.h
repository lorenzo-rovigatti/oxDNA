/*
 * FSConf.h
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#ifndef FSCONF_H_
#define FSCONF_H_

#include "Observables/Configurations/Configuration.h"

/**
 * @brief Prints a configuration in the FS format.
 *
 * @verbatim
 [in_box = <bool> (if true all the positions are brought back between -L/2 and L/2. Defaults to false)]
 @endverbatim
 */
class FSConf: public Configuration {
protected:
	int _N = -1;
	int _N_A = -1;
	int _N_B = -1;
	bool _in_box = false;
	bool _also_patch = false;
	bool _print_bonds = false;
	int _energy_term_id = -1;
	number _bond_threshold = -0.2;

	std::vector<std::map<int, int> > _bonds;

	virtual std::string _headers(llint step);
	virtual std::string _particle(BaseParticle *p);

public:
	FSConf();
	virtual ~FSConf();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();

	std::string _configuration(llint step);
};

extern "C" BaseObservable *make_FSConf() {
	return new FSConf();
}

#endif /* FSCONF_H_ */
