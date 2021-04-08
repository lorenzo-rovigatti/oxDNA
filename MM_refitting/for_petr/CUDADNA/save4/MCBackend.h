/*
 * MCBackend.h
 *
 *  Created on: 25/nov/2010
 *      Author: lorenzo
 */

#ifndef MCBACKEND_H_
#define MCBACKEND_H_

#define MC_ENSEMBLE_NVT 0
#define MC_ENSEMBLE_NPT 1

#define MC_MOVES 5
#define MC_MOVE_TRANSLATION 0
#define MC_MOVE_ROTATION 1
#define MC_MOVE_VOLUME 2
#define MC_MOVE_CLUSTER_SIZE 3
#define MC_MOVE_NEW_BOND 4

#include "SimBackend.h"

template<typename number>
class MCBackend: public SimBackend<number> {
	friend class IOManager;
protected:
	bool _overlap;
	number _delta[MC_MOVES];
	long long int _tries[MC_MOVES];
	long long int _accepted[MC_MOVES];
	int _ensemble;
	int _check_energy_every;
	int _check_energy_counter;
	number _check_energy_threshold;

	void _get_number_settings(input_file &inp);
	virtual void _compute_energy() = 0;

public:
	MCBackend(IOManager *IO);
	virtual ~MCBackend();

	void get_settings(input_file &inp);
	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);

	void print_energy(llint curr_step);
	void print_conf(llint curr_step, bool reduced=false, bool only_last=false);
};

#endif /* MCBACKEND_H_ */
