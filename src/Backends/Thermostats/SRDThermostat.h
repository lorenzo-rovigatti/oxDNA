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
	SRDParticle<number> *next;
};

/**
 * @brief Incapsulates a single cell used by the {@link SRDThermostat SRD thermostat}.
 */
template<typename number>
struct SRDCell {
	LR_vector<number> v_avg;
	number tot_mass;
	SRDParticle<number> *head;
};

/**
 * @brief Incapsulates a stochastic rotation dynamics thermostat (see https://en.wikipedia.org/wiki/Multi-particle_collision_dynamics).
 */
template<typename number>
class SRDThermostat : public BaseThermostat<number> {
private:
	int _N_per_cell;
	number _r_cell;
	int _N_cells;
	int _N_cells_side;
	SRDCell<number> *_cells;
	int _N_particles;
	/// The first _N_particles are solvent particles, the other are the regular particles
	SRDParticle<number> *_srd_particles;
	/// Since this is a reference, srd *should* also work for variable boxes, should we ever implement an MD barostat
	number &_box_side;
	number _T;
	number _dt;
	int _apply_every;
	/// Collision angle
	number _alpha;
	/// SRD particle mass
	number _m;

	int _get_cell_index(LR_vector<number> &r);

public:
	SRDThermostat(number &box_side);
	virtual ~SRDThermostat();

	void get_settings(input_file &inp);
	void init(int N_part);
	void apply(BaseParticle<number> **particles, llint curr_step);
};

#endif /* SRDTHERMOSTAT_H_ */
