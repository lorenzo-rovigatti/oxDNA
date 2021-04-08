/*
 * MDBackend.h
 *
 *  Created on: 23/nov/2010
 *      Author: lorenzo
 */

#ifndef MDBACKEND_H_
#define MDBACKEND_H_

#include "SimBackend.h"

using namespace std;

class IOManager;

template<typename number>
class MDBackend: public SimBackend<number> {
	friend class IOManager;
protected:
	number _dt;
	bool _refresh_velocities;

	// thermostat's parameters
	number _rescale_factor;
	int _newtonian_steps;
	int _thermostat;
	number _pt, _pr;
	number _diff_coeff;

	void _get_number_settings(input_file &inp);

	void _generate_vel();

public:
	enum {THERMOSTAT_NO = 0, THERMOSTAT_JOHN = 1, THERMOSTAT_REFRESH = 2};

	MDBackend(IOManager *IO);
	virtual ~MDBackend();

	void get_settings(input_file &inp);
	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);

	void print_energy(llint curr_step);
	void print_conf(llint curr_step, bool reduced=false, bool only_last=false);
};

#endif /* MDBACKEND_H_ */
