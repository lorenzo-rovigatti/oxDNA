/*
 * CPAnalysis.h
 *
 *  Created on: Jul 19, 2016
 *      Author: rovigatti
 */

#ifndef CPANALYSIS_H_
#define CPANALYSIS_H_

#include "Observables/BaseObservable.h"
#include "Lists/Cells.h"

template<typename number>
class CPCells : public Cells<number> {
public:
	vector<number> cell_density;
	vector<number> coarse_grained_cell_density;

	CPCells(int &N, BaseBox<number> *box);
	virtual ~CPCells();

	virtual number assign_density_to_cells();
	virtual int get_N_particles_in_cell(int cell_idx);
};

/**
 * @brief Outputs the total kinetic energy of the system.
 */
template<typename number>
class CPAnalysis : public BaseObservable<number> {
protected:
	int _type;
	number _sigma;
	CPCells<number> *_cells;

public:
	CPAnalysis();
	virtual ~CPAnalysis();

	void get_settings (input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);

	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable<float> *make_float() { return new CPAnalysis<float>(); }
extern "C" BaseObservable<double> *make_double() { return new CPAnalysis<double>(); }

#endif /* CPANALYSIS_H_ */
