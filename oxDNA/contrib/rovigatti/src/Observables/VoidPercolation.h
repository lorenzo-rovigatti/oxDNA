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

template<typename number>
class VPCells : public Cells<number> {
	vector<int> _clusters;
	vector<int> _sizes;
	map<int, LR_vector<number> > _color_map;

	void _flip_neighs(BaseParticle<number> *);
public:
	vector<int> csd;

	VPCells(int &N, BaseBox<number> *box);
	virtual ~VPCells();

	void compute_csd();
	std::string get_coloured_mgl(number);
};

template<typename number>
class VoidPercolation : public BaseObservable<number>  {
protected:
	number _probe_diameter;
	number _particle_diameter;
	number _rcut, _sqr_rcut;
	Cells<number> *_cells;
	VPCells<number> *_probe_cells;
	VPCells<number> *_replica_probe_cells;
	BaseParticle<number> **_probes;
	BaseParticle<number> **_replica_probes;
	BaseParticle<number> **_particles;
	BaseParticle<number> *_probe;
	BaseBox<number> *_replica_box;
	int _N, _N_replica;
	int _insertions;
	bool _print_mgl;

public:
	VoidPercolation();
	virtual ~VoidPercolation();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new VoidPercolation<float>(); }
extern "C" BaseObservable<double> *make_double() { return new VoidPercolation<double>(); }

#endif /* VOIDPERCOLATION_H_ */
