/*
 * DeltaBackend.h
 *
 *  Created on: 20/set/2011
 *      Author: lorenzo
 */

#ifndef DELTABACKEND_H_
#define DELTABACKEND_H_

#include "../../MC_CPUBackend.h"
#include "Strand.h"

class DeltaBackend : public MC_CPUBackend<double> {
protected:
	LR_vector<double> _patches_on0[4];
	double _min_corona, _max_corona;
	Strand *_strands[2];
	double _initial_U;

	virtual void _compute_energy(bool all_interactions);

public:
	DeltaBackend(IOManager *IO);
	virtual ~DeltaBackend();

	virtual void set_confs_to_skip(int n) { _confs_to_skip = n; };
	//virtual double get_trial_volume() { return _box_side*_box_side*_box_side; }
	virtual double get_trial_volume();
	virtual void init(ifstream &conf_input);
	virtual double insertion_energy();
};

#endif /* DELTABACKEND_H_ */
