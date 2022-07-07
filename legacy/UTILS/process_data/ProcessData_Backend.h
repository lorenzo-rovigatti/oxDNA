/*
 * ProcessData_Backend.h
 *
 *  Created on: 03/set/2010
 *      Author:
 */

#ifndef PROC_CPUBACKEND_H_
#define PROC_CPUBACKEND_H_

#include "../../SimBackend.h"
using namespace std;

template<typename number>
class ProcessData_Backend: public SimBackend<number> {
protected:
	//LR_vector<number> _vcm;

	//void _first_step(llint cur_step);
	number _excluded_volume(const LR_vector<number> &r,
			LR_vector<number> &force, number sigma, number rstar, number b,
			number rc);
	number _particle_particle_bonded_interaction(Particle<number> *p);
	number _particle_particle_interaction(Particle<number> *p,
			Particle<number> *q);

	number mc_particle_particle_bonded_interaction(Particle<number> *p);
    number mc_particle_particle_interaction(Particle<number> *p,
				Particle<number> *q);
	void _compute_forces();
	//void _second_step();
	//void _activate_john_thermostat();
	//void _activate_refresh_thermostat();

	//void _rescale_velocities();

	 number _particle_particle_bonded_stacking(Particle<number> *p
			); //neighbortype is 3 or 5; returns 0 if neighbor is virtual
	 number _particle_particle_coaxial_stacking(Particle<number> *p,
			Particle<number> *q);
	 number _particle_particle_cross_stacking(Particle<number> *p,
			Particle<number> *q);
	 number _particle_particle_H_bond(Particle<number> *p, Particle<
			number> *q);
	 number _particle_particle_excluded_volume(Particle<number> *p,
			Particle<number> *q);
	 number _particle_particle_bonded_excluded_volume(
			Particle<number> *p);
	 number _excluded_volume(const LR_vector<number> &r, number sigma,
			number rstar, number b, number rc);
	 char   _particle_particle_H_bond_status(Particle<number> *p, Particle<number> *q);

	 number _particle_particle_bonded_FENE(Particle<number> *p);

public:
	ProcessData_Backend(IOManager *IO, int confs_to_skip = 0);
	virtual ~ProcessData_Backend();
	ostream& OutputBonds(
		 		ostream& out, int t = 0);

	//void init(ifstream &conf_input);
	void init (char conf_filename[256]);
	int init_first_state(ifstream& conf_input);
	int init_to_next_state(ifstream &conf_input); //this loads one configuration, returns the loaded iteration

	virtual void sim_step(llint) {}
	virtual void print_energy(llint) {}
	virtual void print_conf(llint, bool, bool) {}

	//void sim_step(llint cur_step);
	//void activate_thermostat();
};

#endif /* PROCESSDATA_BACKEND_H_ */
