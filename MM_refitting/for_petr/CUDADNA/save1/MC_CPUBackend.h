/*
 * MC_CPUBackend.h
 *
 *  Created on: 26/nov/2010
 *      Author: lorenzo
 */

#ifndef MC_CPUBACKEND_H_
#define MC_CPUBACKEND_H_

#include "MCBackend.h"

template<typename number>
class MC_CPUBackend: public MCBackend<number> {
protected:
	Particle<number> *_particles_old;

	inline number _excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc);
	//inline number _particle_particle_bonded_interaction(Particle<number> *p);
	//inline number _particle_particle_interaction(Particle<number> *p, Particle<number> *q);
	// inlining gives warning when inherited; compiler does not inline them anyways
	number _particle_particle_bonded_interaction(Particle<number> *p, bool update=false);
	number _particle_particle_interaction(Particle<number> *p, Particle<number> *q);
	void _compute_energy();
	inline void _translate_particle(Particle<number> *p);
	inline void _rotate_particle(Particle<number> *p);
	inline number _particle_energy(Particle<number> *p);

public:
	MC_CPUBackend(IOManager *IO);
	virtual ~MC_CPUBackend();

	//void init(ifstream &conf_input);
	void init(char conf_filename[256]);

	void sim_step(llint cur_step);
	void print_energy(llint curr_step);
};

#endif /* MC_CPUBACKEND_H_ */
