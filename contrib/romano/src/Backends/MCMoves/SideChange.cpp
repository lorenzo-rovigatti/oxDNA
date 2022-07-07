/**
 * @file    SideChange.cpp
 * @date    25/oct/2018
 * @author  flavio
 *
 *
 */

#include "SideChange.h"

SideChange::SideChange () {
	_delta = -1.f;
	_apply_x = false;
	_apply_y = false;
	_apply_z = false;
	_verlet_skin = -1.f;
	_P = 0.f;
}

SideChange::~SideChange () {

}

void SideChange::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "prob", &this->prob, 0);
	getInputNumber (&sim_inp, "P", &_P, 1);
	
	std::string tmps;
	getInputString (&inp, "axes", tmps, 1);
	
	int found = 0;
	if (!tmps.compare("x")) {
	   _apply_x = true;
	   found ++;
	}
	if (!tmps.compare("y")) {
		_apply_y = true;
		found ++;
	}	
	if (!tmps.compare("z")) {
		_apply_z = true;
		found ++;
	}
	if (!tmps.compare("xy")) {
		_apply_x = true;
		_apply_y = true;
		found ++;
	}
	if (!tmps.compare("xz")) {
		_apply_x = true;
		_apply_z = true;
		found ++;
	}
	if (!tmps.compare("yz")) {
		_apply_y = true;
		_apply_z = true;
		found ++;
	}
	if (found != 1) throw oxDNAException ("Error in SideChange move initialization: could not understand which axes to act upon.\n\tFound axes = %s, should be one of x, y, z, xy, xz, yz");

	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}

	// orthogonal box check?
}

void SideChange::init () {
	BaseMove::init();
	_pos_old.resize (_Info->N());
	if (this->_restrict_to_type > 0) OX_LOG(Logger::LOG_WARNING, "(SideChange.cpp) Cant use SideChange with restrict_to_type. Ignoring");
	OX_LOG(Logger::LOG_INFO, "(SideChange.cpp) SideChange initiated with T=%g, P=%g, delta=%g, prob=%g", this->_T, _P, _delta, this->prob);
	OX_LOG(Logger::LOG_INFO, "(SideChange.cpp)                      apply_x = %d, apply_y = %d, apply_z = %d (true=%d)", (int)_apply_x, (int)_apply_y, (int)_apply_z, (int)true);
}

void SideChange::apply (llint curr_step) {
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
	
	number mydelta = _delta * (drand48() - 0.5);
	if (_apply_x) box_sides.x += mydelta;
	if (_apply_y) box_sides.y += mydelta;
	if (_apply_z) box_sides.z += mydelta;

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

void SideChange::log_parameters() {
	BaseMove::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}

