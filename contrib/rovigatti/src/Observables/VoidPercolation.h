/*
 * VoidPercolation.h
 *
 *  Created on: Sep 27, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef VOIDPERCOLATION_H_
#define VOIDPERCOLATION_H_

#include "Lists/Cells.h"

#include "Observables/BaseObservable.h"
#include "Interactions/DNAInteraction.h"
#include "Utilities/OrderParameters.h"

class VPCells: public Cells {
	vector<int> _clusters;
	vector<int> _sizes;
	map<int, LR_vector> _color_map;

	void _flip_neighs(BaseParticle *);
public:
	vector<int> csd;

	VPCells(std::vector<BaseParticle *> &ps, BaseBox *box);
	VPCells() = delete;
	virtual ~VPCells();

	void compute_csd();
	std::string get_coloured_mgl(number);
};

class VoidPercolation: public BaseObservable {
protected:
	number _probe_diameter;
	number _particle_diameter;
	number _rcut, _sqr_rcut;
	Cells *_cells;
	VPCells *_probe_cells;
	VPCells *_replica_probe_cells;
	std::vector<BaseParticle *> _probes;
	std::vector<BaseParticle *> _replica_probes;
	std::vector<BaseParticle *> _particles;
	BaseParticle *_probe;
	BoxPtr _replica_box;
	int _N, _N_replica;
	int _insertions;
	bool _print_mgl;

public:
	VoidPercolation();
	virtual ~VoidPercolation();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init();
};

extern "C" BaseObservable *make_VoidPercolation() {
	return new VoidPercolation();
}

#endif /* VOIDPERCOLATION_H_ */
