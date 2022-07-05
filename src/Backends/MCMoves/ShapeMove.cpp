/**
 * @file    ShapeMove.cpp
 * @date    30/may/2015
 * @author  flavio
 *
 *
 */

#include "ShapeMove.h"

/// traslation

ShapeMove::ShapeMove () {
	_verlet_skin = -1.f;
}


ShapeMove::~ShapeMove () {

}


void ShapeMove::init () {
	BaseMove::init();
	_pos_old.resize (_Info->N());
	if (this->_restrict_to_type > 0) OX_LOG(Logger::LOG_WARNING, "(ShapeMove.cpp) Cant use ShapeMove with restrict_to_type. Ignoring");
	OX_LOG(Logger::LOG_INFO, "(ShapeMove.cpp) ShapeMove initiated with T %g, delta %g, prob: %g", this->_T, _delta, this->prob);
}


void ShapeMove::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "prob", &this->prob, 0);

	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}


void ShapeMove::apply (llint curr_step) {
	// we increase the attempted count
	this->_attempted += 1;

	std::vector<BaseParticle *> & particles = this->_Info->particles();
	int N = this->_Info->N();

	// select axis NOT to change
	int preserved_axis = ((int) lrand48()) % 3;
	int change_axis_1, change_axis_2;
	if (drand48() > 0.5) {
		change_axis_1 = (preserved_axis + 1) % 3;
		change_axis_2 = (preserved_axis + 2) % 3;
	}
	else {
		change_axis_1 = (preserved_axis + 2) % 3;
		change_axis_2 = (preserved_axis + 1) % 3;
	}

	//printf ("@#@@ preserving axis %d, changing %d and %d\n", preserved_axis, change_axis_1, change_axis_2);

	number dL = _delta * (drand48() - (number) 0.5);
	LR_vector box_sides = this->_Info->box->box_sides();
	LR_vector old_box_sides = this->_Info->box->box_sides();

	number oldE;
	if (this->_compute_energy_before) oldE = this->_Info->interaction->get_system_energy(this->_Info->particles(), this->_Info->lists);
	else oldE = (number) 0.f;
	if (this->_Info->interaction->get_is_infinite()) printf ("WHAAT\n");
	number oldV = (box_sides[0]) * (box_sides[1]) * (box_sides[2]);

	box_sides[change_axis_1] += dL;
	box_sides[change_axis_2] = oldV / (box_sides[preserved_axis] * box_sides[change_axis_1]);

	//printf ("@#@@ old sides: %g %g %g\n", old_box_sides[0], old_box_sides[1], old_box_sides[2]);
	//printf ("@#@@ old sides: %g %g %g\n", box_sides[0], box_sides[1], box_sides[2]);

	//this->_Info->interaction->box->init (box_sides[0], box_sides[1], box_sides[2]);
	this->_Info->box->init(box_sides.x, box_sides.y, box_sides.z);

	number dExt = (number) 0.f;
	for (int k = 0; k < N; k ++) {
		BaseParticle *p = particles[k];
		dExt -= -p->ext_potential;
		_pos_old[k] = p->pos;
		// p->pos *= (1. + dL / box_sides[0]);
		p->pos.x *= (box_sides.x / old_box_sides.x);
		p->pos.y *= (box_sides.y / old_box_sides.y);
		p->pos.z *= (box_sides.z / old_box_sides.z);

		p->set_ext_potential(curr_step, this->_Info->box);
		dExt += -p->ext_potential;
	}

	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	number newE = this->_Info->interaction->get_system_energy(this->_Info->particles(), this->_Info->lists);
	number dE = newE - oldE + dExt;
	number V = box_sides.x * box_sides.y * box_sides.z;
	number dV = V - oldV;

	if (fabs (dV) > 1.e-6) throw oxDNAException("Too high dV: %g\n", dV);

	if (this->_Info->interaction->get_is_infinite() == false && exp(- dE / this->_T) > drand48()) {
		this->_accepted ++;
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta *= this->_acc_fact;
		CONFIG_INFO->notify(CONFIG_INFO->box->UPDATE_EVENT);
	}
	else {
		//printf ("reject: dE = %g\n", dE);
		for (int k = 0; k < N; k ++) {
			BaseParticle *p = particles[k];
			p->pos = _pos_old[k];
			p->set_ext_potential(curr_step, this->_Info->box);
		}
		this->_Info->box->init(old_box_sides.x, old_box_sides.y, old_box_sides.z);
		this->_Info->interaction->set_is_infinite (false);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;

		//number xE  = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
		//if (fabs ((xE - oldE) / oldE) > 1.e-5) throw oxDNAException ("Look at this shit %g %g", xE, oldE);
	}

	return;
}


void ShapeMove::log_parameters() {
	BaseMove::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}
