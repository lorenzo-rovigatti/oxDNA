/**
 * @file    VolumeMove.cpp
 * @date    30/may/2015
 * @author  flavio
 *
 *
 */

#include "VolumeMove.h"

/// traslation
VolumeMove::VolumeMove() {
	_verlet_skin = -1.f;
	_isotropic = true;
	_P = 0.;
	_delta = 0.;
}

VolumeMove::~VolumeMove() {

}

void VolumeMove::init() {
	BaseMove::init();
	_pos_old.resize(_Info->N());
	if(_restrict_to_type > 0) {
		OX_LOG(Logger::LOG_WARNING, "(VolumeMove.cpp) Cant use VolumeMove with restrict_to_type. Ignoring");
	}
	OX_LOG(Logger::LOG_INFO, "(VolumeMove.cpp) VolumeMove (isotropic = %d) initiated with T %g, delta %g, prob: %g", _isotropic, _T, _delta, prob);
}

void VolumeMove::get_settings(input_file &inp, input_file &sim_inp) {
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

void VolumeMove::apply(llint curr_step) {
	// we increase the attempted count
	_attempted += 1;

	std::vector<BaseParticle*> &particles = _Info->particles();
	int N = _Info->N();

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
		box_sides.x += dL;
		box_sides.y += dL;
		box_sides.z += dL;
	}
	else {
		box_sides.x += _delta * (drand48() - (number) 0.5);
		box_sides.y += _delta * (drand48() - (number) 0.5);
		box_sides.z += _delta * (drand48() - (number) 0.5);
	}

	_Info->box->init(box_sides[0], box_sides[1], box_sides[2]);
	number dExt = (number) 0.f;
	for(int k = 0; k < N; k++) {
		BaseParticle *p = particles[k];
		dExt -= p->ext_potential;
		_pos_old[k] = p->pos;
		p->pos.x *= box_sides[0] / old_box_sides[0];
		p->pos.y *= box_sides[1] / old_box_sides[1];
		p->pos.z *= box_sides[2] / old_box_sides[2];
		p->set_ext_potential(curr_step, _Info->box);
		dExt += p->ext_potential;
	}

	// this bit has to come after the update of particles' positions
	if(!_Info->lists->is_updated()) {
		_Info->lists->global_update();
	}

	number newE = _Info->interaction->get_system_energy(_Info->particles(), _Info->lists);
	number dE = newE - oldE + dExt;
	number V = _Info->box->V();
	number dV = V - oldV;

	if(_Info->interaction->get_is_infinite() == false && exp(-(dE + _P * dV - N * _T * log(V / oldV)) / _T) > drand48()) {
		_accepted++;
		if(curr_step < _equilibration_steps && _adjust_moves) {
			_delta *= _acc_fact;
		}
		CONFIG_INFO->notify(CONFIG_INFO->box->UPDATE_EVENT);
	}
	else {
		for(int k = 0; k < N; k++) {
			BaseParticle *p = particles[k];
			p->pos = _pos_old[k];
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

void VolumeMove::log_parameters() {
	BaseMove::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta %g", _delta);
}
