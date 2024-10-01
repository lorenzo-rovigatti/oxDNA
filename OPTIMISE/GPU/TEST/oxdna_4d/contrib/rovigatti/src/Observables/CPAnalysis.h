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

class CPCells: public Cells {
public:
	std::vector<number> cell_density;
	std::vector<number> coarse_grained_cell_density;

	CPCells(std::vector<BaseParticle *> &ps, BaseBox *box);
	CPCells() = delete;
	virtual ~CPCells();

	virtual number assign_density_to_cells();
	virtual int get_N_particles_in_cell(int cell_idx);
};

/**
 * @brief Outputs the total kinetic energy of the system.
 */

class CPAnalysis: public BaseObservable {
protected:
	int _type;
	number _sigma;
	CPCells *_cells;

public:
	CPAnalysis();
	virtual ~CPAnalysis();

	void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();

	std::string get_output_string(llint curr_step);
};

extern "C" BaseObservable *make_CPAnalysis() {
	return new CPAnalysis();
}

#endif /* CPANALYSIS_H_ */
