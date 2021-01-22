/*
 * SRDThermostat.h
 *
 *  Created on: 05/dic/2013
 *      Author: lorenzo
 */

#ifndef SRDTHERMOSTAT_H_
#define SRDTHERMOSTAT_H_

#include "BaseThermostat.h"

/**
 * @brief Incapsulates "fluid particles" used by the {@link SRDThermostat SRD thermostat}.
 */

struct SRDParticle {
	LR_vector r;
	LR_vector v;
	LR_vector L;
	SRDParticle *next;
	int cell_index;
};

/**
 * @brief Incapsulates a single cell used by the {@link SRDThermostat SRD thermostat}.
 */

struct SRDCell {
	LR_vector centre;
	LR_vector P;
	LR_vector L;
	LR_vector L_spin;
	LR_vector PR;
	LR_vector LR;
	LR_vector LR_spin;
	LR_vector dLgn;
	number tot_mass;
	number tot_I;
	SRDParticle *head;
};

/**
 * @brief Incapsulates a stochastic rotation dynamics thermostat (see https://en.wikipedia.org/wiki/Multi-particle_collision_dynamics).
 */

class SRDThermostat: public BaseThermostat {
protected:
	int _N_per_cell;
	number _r_cell;
	int _N_cells;
	int _N_cells_side;
	SRDCell *_cells;
	int _N_particles;
	number _rescale_factor;
	/// The first _N_particles are solvent particles, the others are the regular particles
	SRDParticle *_srd_particles;
	/// Since this is a reference, srd *should* also work for variable boxes, should we ever implement an MD barostat
	BaseBox *_box;
	number _dt;
	int _apply_every;
	/// SRD particle mass
	number _m;

	bool _is_cuda;

	int _get_cell_index(LR_vector &r);

public:
	SRDThermostat(BaseBox * box);
	virtual ~SRDThermostat();

	void get_settings(input_file &inp);
	void init();
	void apply(std::vector<BaseParticle *> &particles, llint curr_step);
	void apply1(std::vector<BaseParticle *> &particles, llint curr_step);
	void apply2(std::vector<BaseParticle *> &particles, llint curr_step);
	void apply3(std::vector<BaseParticle *> &particles, llint curr_step);
};

#endif /* SRDTHERMOSTAT_H_ */
