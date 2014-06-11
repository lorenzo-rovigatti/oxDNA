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
template<typename number>
struct SRDParticle {
	LR_vector<number> r;
	LR_vector<number> v;
	LR_vector<number> L;
	SRDParticle<number> *next;
	int cell_index;
};

/**
 * @brief Incapsulates a single cell used by the {@link SRDThermostat SRD thermostat}.
 */
template<typename number>
struct SRDCell {
	LR_vector<number> centre;
	LR_vector<number> P;
	LR_vector<number> L;
	LR_vector<number> L_spin;
	LR_vector<number> PR;
	LR_vector<number> LR;
	LR_vector<number> LR_spin;
	LR_vector<number> dLgn;
	number tot_mass;
	number tot_I;
	SRDParticle<number> *head;
};

/**
 * @brief Incapsulates a stochastic rotation dynamics thermostat (see https://en.wikipedia.org/wiki/Multi-particle_collision_dynamics).
 */
template<typename number>
class SRDThermostat : public BaseThermostat<number> {
protected:
	int _N_per_cell;
	number _r_cell;
	int _N_cells;
	int _N_cells_side;
	SRDCell<number> *_cells;
	int _N_particles;
	number _rescale_factor;
	/// The first _N_particles are solvent particles, the otherl are the regular particles
	SRDParticle<number> *_srd_particles;
	/// Since this is a reference, srd *should* also work for variable boxes, should we ever implement an MD barostat
	number &_box_side;
	number _T;
	number _dt;
	int _apply_every;
	/// SRD particle mass
	number _m;

	bool _is_cuda;

	int _get_cell_index(LR_vector<number> &r);

public:
	SRDThermostat(number &box_side);
	virtual ~SRDThermostat();

	void get_settings(input_file &inp);
	void init(int N_part);
	void apply(BaseParticle<number> **particles, llint curr_step);
	void apply1(BaseParticle<number> **particles, llint curr_step);
	void apply2(BaseParticle<number> **particles, llint curr_step);
	void apply3(BaseParticle<number> **particles, llint curr_step);
};

#endif /* SRDTHERMOSTAT_H_ */
