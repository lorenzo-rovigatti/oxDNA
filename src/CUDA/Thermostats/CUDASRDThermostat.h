/*
 * CUDASRDThermostat.h
 *
 *  Created on: Feb 15, 2013
 *      Author: rovigatti
 */

#ifndef CUDASRDTHERMOSTAT_H_
#define CUDASRDTHERMOSTAT_H_

#include "CUDABaseThermostat.h"
#include "../../Backends/Thermostats/SRDThermostat.h"

/**
 * @brief CUDA implementation of a {@link SRDThermostat SRD thermostat}.
 */
template<typename number, typename number4>
class CUDASRDThermostat: public CUDABaseThermostat<number, number4>, public SRDThermostat<number> {
protected:
	int *_d_cells;
	int *_d_counters_cells;
	bool *_d_cell_overflow;

	number4 *_d_poss;
	number4 *_d_vels;

	/// the fourth component of each element stores the mass of the particle
	number4 *_d_cells_dp;
	number4 *_d_reduced_cells_dp;
	int *_d_reduce_keys;
	int *_d_reduced_cells_keys;

	int _max_N_per_cell;
	int _N_tot;
	int _N_vec_size;

public:
	CUDASRDThermostat(number &box_side);
	virtual ~CUDASRDThermostat();

	virtual void get_settings(input_file &inp);
	virtual void init(int N);

	virtual void apply_cuda(number4 *d_poss, GPU_quat<number> *d_orientations, number4 *d_vels, number4 *d_Ls, llint curr_step);
	virtual bool would_activate(llint curr_step);
};

#endif /* CUDASRDTHERMOSTAT_H_ */
