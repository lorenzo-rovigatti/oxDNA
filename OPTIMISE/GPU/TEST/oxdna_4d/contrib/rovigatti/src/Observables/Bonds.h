/*
 * Bonds.h
 *
 *  Created on: 19/jul/2021
 *      Author: lorenzo
 */

#ifndef BONDS_H_
#define BONDS_H_

#include "Observables/Configurations/Configuration.h"

/**
 * @brief Prints the list of bonds per particle.
 */
class Bonds: public BaseObservable {
protected:
	int _energy_term_id = -1;
	number _bond_threshold = -0.2;

	std::vector<std::map<int, int> > _bonds;

public:
	Bonds();
	virtual ~Bonds();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	void init();

	std::string get_output_string(llint step);
};

extern "C" BaseObservable *make() {
	return new Bonds();
}

#endif /* BONDS_H_ */
