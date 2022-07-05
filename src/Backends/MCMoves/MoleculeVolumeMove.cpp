/**
 * @file    MoleculeMoleculeVolumeMove.cpp
 * @date    21/apr/2020
 * @author  lorenzo
 */

#include "MoleculeVolumeMove.h"

#include "../../Particles/Molecule.h"

/// traslation
MoleculeVolumeMove::MoleculeVolumeMove() {

}

MoleculeVolumeMove::~MoleculeVolumeMove() {

}

void MoleculeVolumeMove::init() {
	BaseMove::init();
	if(_restrict_to_type > 0) {
		OX_LOG(Logger::LOG_WARNING, "(MoleculeVolumeMove.cpp) Cant use MoleculeMoleculeVolumeMove with restrict_to_type. Ignoring");
	}
	OX_LOG(Logger::LOG_INFO, "(VolumeMove.cpp) MoleculeVolumeMove (isotropic = %d) initiated with T %g, delta %g, prob: %g", _isotropic, _T, _delta, prob);
}

void MoleculeVolumeMove::get_settings(input_file &inp, input_file &sim_inp) {
	BaseMove::get_settings(inp, sim_inp);

	getInputBool(&inp, "isotropic", &_isotropic, 0);
	getInputNumber(&inp, "delta", &_delta, 1);
	getInputNumber(&inp, "prob", &prob, 0);
	getInputNumber(&sim_inp, "P", &_P, 1);

	std::string tmps;
	if(getInputString(&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if(!tmps.compare("verlet")) {
			getInputNumber(&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}

void MoleculeVolumeMove::apply(llint curr_step) {
	// we increase the attempted count
	_attempted += 1;

	std::vector<BaseParticle*> &particles = _Info->particles();
	int N_molecules = _Info->molecules().size();

	LR_vector box_sides = _Info->box->box_sides();
	LR_vector old_box_sides = box_sides;

	number oldE;
	if(_compute_energy_before) {
		oldE = _Info->interaction->get_system_energy(_Info->particles(), _Info->lists);
	}
	else {
		oldE = (number) 0.f;
	}
	number oldV = _Info->box->V();

	if(_isotropic) {
		number dL = _delta * (drand48() - (number) 0.5);
		box_sides[0] += dL;
		box_sides[1] += dL;
		box_sides[2] += dL;
	}
	else {
		box_sides[0] += _delta * (drand48() - (number) 0.5);
		box_sides[1] += _delta * (drand48() - (number) 0.5);
		box_sides[2] += _delta * (drand48() - (number) 0.5);
	}

	_Info->box->init(box_sides[0], box_sides[1], box_sides[2]);
	number dExt = (number) 0.f;

	LR_vector shift_factor(
		box_sides[0] / old_box_sides[0] - 1.,
		box_sides[1] / old_box_sides[1] - 1.,
		box_sides[2] / old_box_sides[2] - 1.
	);

	static std::vector<LR_vector> shifts(N_molecules);
	// we account for possible changes to the number of molecules
	shifts.resize(N_molecules);
	int i = 0;
	for(auto mol : _Info->molecules()) {
		mol->update_com();
		shifts[i] = LR_vector(
			mol->com[0] * shift_factor[0],
			mol->com[1] * shift_factor[1],
			mol->com[2] * shift_factor[2]
		);
		i++;
	}

	for(auto p : particles) {
		dExt -= p->ext_potential;
	}

	for(auto p : particles) {
		p->pos += shifts[p->strand_id];
	}

	for(auto p : particles) {
		p->set_ext_potential(curr_step, _Info->box);
		dExt += p->ext_potential;
	}

	// this bit must be executed after updating the particles' positions
	if(!_Info->lists->is_updated()) {
		_Info->lists->global_update();
	}

	number newE = _Info->interaction->get_system_energy(_Info->particles(), _Info->lists);
	number dE = newE - oldE + dExt;
	number V = _Info->box->V();
	number dV = V - oldV;

	if(_Info->interaction->get_is_infinite() == false && exp(-(dE + _P * dV - (N_molecules + 1) * _T * log(V / oldV)) / _T) > drand48()) {
		_accepted++;
		if(curr_step < _equilibration_steps && _adjust_moves) {
			_delta *= _acc_fact;
		}
		CONFIG_INFO->notify(CONFIG_INFO->box->UPDATE_EVENT);
	}
	else {
		for(auto p : particles) {
			p->pos -= shifts[p->strand_id];
		}

		for(auto p : particles) {
			p->set_ext_potential(curr_step, _Info->box);
		}

		_Info->interaction->set_is_infinite(false);
		_Info->box->init(old_box_sides[0], old_box_sides[1], old_box_sides[2]);

		if(!_Info->lists->is_updated()) {
			_Info->lists->global_update();
		}

		if(curr_step < _equilibration_steps && _adjust_moves) {
			_delta /= _rej_fact;
		}
	}
	return;
}

void MoleculeVolumeMove::log_parameters() {
	BaseMove::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta %g", _delta);
}
