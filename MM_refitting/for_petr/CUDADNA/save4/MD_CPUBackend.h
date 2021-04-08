/*
 * MD_CPUBackend.h
 *
 *  Created on: 03/set/2010
 *      Author: lorenzo
 */

#ifndef MD_CPUBACKEND_H_
#define MD_CPUBACKEND_H_

#include "MDBackend.h"

template<typename number>
class MD_CPUBackend: public MDBackend<number> {
protected:
	LR_vector<number> _vcm;

	void _first_step(llint cur_step);
	number _excluded_volume(const LR_vector<number> &r, LR_vector<number> &force, number sigma, number rstar, number b, number rc);
	number _particle_particle_bonded_interaction(Particle<number> *p);
	number _particle_particle_interaction(Particle<number> *p, Particle<number> *q);
	void _compute_forces();
	void _second_step();
	void _activate_john_thermostat();
	void _activate_refresh_thermostat();

	void _rescale_velocities();

public:
	MD_CPUBackend(IOManager *IO);
	virtual ~MD_CPUBackend();

	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);

	void sim_step(llint cur_step);
	void activate_thermostat();
};

#endif /* MD_CPUBACKEND_H_ */
