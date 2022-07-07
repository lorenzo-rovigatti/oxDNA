/**
 * @file    CutVolume.cpp
 * @date    25/oct/2018
 * @author  flavio
 *
 *
 */

#include "CutVolume.h"

/// traslation
template<typename number>
CutVolume<number>::CutVolume () {
	_verlet_skin = -1.f;
	_P = 0.f;
}

template<typename number>
CutVolume<number>::~CutVolume () {

}

template<typename number>
void CutVolume<number>::init () {
	BaseMove<number>::init();
	_pos_old.resize (*this->_Info->N);
	if (this->_restrict_to_type > 0) OX_LOG(Logger::LOG_WARNING, "(CutVolume.cpp) Cant use CutVolume with restrict_to_type. Ignoring");
	OX_LOG(Logger::LOG_INFO, "(CutVolume.cpp) CutVolume initiated with T %g, P=%g, delta %g, prob: %g", this->_T, _P, _delta, this->prob);
}

template<typename number>
void CutVolume<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "prob", &this->prob, 0);
	getInputNumber (&sim_inp, "P", &_P, 1);

	std::string tmps;
	if (getInputString (&sim_inp, "list_type", tmps, 0) == KEY_FOUND) {
		if (!tmps.compare ("verlet")) {
			getInputNumber (&sim_inp, "verlet_skin", &_verlet_skin, 1);
		}
	}
}

template<typename number>
void CutVolume<number>::apply (llint curr_step) {
	// we increase the attempted count
	this->_attempted += 1;

	BaseParticle<number> ** particles = this->_Info->particles;
	int N = *(this->_Info->N);

	// select axis to cut/expand
	int affected_axis = ((int) lrand48()) % 3;

	//printf ("@#@@ preserving axis %d, changing %d and %d\n", preserved_axis, change_axis_1, change_axis_2);

	// chose position
	LR_vector<number> box_sides_old = this->_Info->box->box_sides();
	number cut_at = drand48() * box_sides_old[affected_axis];
	number oldV = (box_sides_old[0] * box_sides_old[1] * box_sides_old[2]);

	number oldE;
	if (this->_compute_energy_before) oldE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	else oldE = (number) 0.f;
	if (this->_Info->interaction->get_is_infinite()) printf ("WHAAT\n");

	bool expand = true;
	if (drand48() < 0.5) {
		expand = false;
	}

	LR_vector<number> box_sides = box_sides_old;
	if (expand) {
		box_sides[affected_axis] += _delta;
	}
	else {
		box_sides[affected_axis] -= _delta;
	}

	//this->_Info->interaction->box->init (box_sides[0], box_sides[1], box_sides[2]);
	this->_Info->box->init(box_sides.x, box_sides.y, box_sides.z);

	number dExt = (number) 0.f;
	bool reject = false;
	for (int k = 0; k < N; k ++) {
		BaseParticle<number> *p = particles[k];
		dExt -= -p->ext_potential;
		_pos_old[k] = p->pos;
		// p->pos *= (1. + dL / box_sides[0]);
		number aux_pos = p->pos[affected_axis];
		while (aux_pos > box_sides_old[affected_axis]) aux_pos -= box_sides_old[affected_axis];
		while (aux_pos < (number)0.f) aux_pos += box_sides_old[affected_axis];
		if (expand) {
			if (aux_pos > cut_at) {
				p->pos[affected_axis] = p->pos[affected_axis] + _delta;
			}
		}
		else {
			number hdelta = _delta / (number) 2.f;
			if (cut_at - hdelta < aux_pos && aux_pos < cut_at + hdelta) {
				reject = true;
			}
			if (aux_pos > cut_at) {
				p->pos[affected_axis] -= _delta;
			}
		}

		// if (reject) break; IMPORTANT: can't do this, since otherwise pos_old[k] won't all be set correctly

		p->set_ext_potential(curr_step, this->_Info->box);
		dExt += -p->ext_potential;
	}

	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	number newE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	number dE = newE - oldE + dExt;
	number V = box_sides.x * box_sides.y * box_sides.z;
	number dV = V - oldV;

	if (reject == false && this->_Info->interaction->get_is_infinite() == false && exp(-(dE + _P*dV - (N + 1)*this->_T*log(V/oldV)) / this->_T) > drand48()) {
		this->_accepted ++;
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta *= this->_acc_fact;
		//number xE  = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
		//if (fabs ((xE - newE) / newE) > 1.e-8) throw oxDNAException ("Look at this shit, JUST ACCEPTED %g %g", xE, oldE);
		//if (this->_Info->interaction->get_is_infinite()) throw oxDNAException ("No, no, no...!");
		// printf ("    accepted cut volume\n");
	}
	else {
		// printf ("    rejected cut volume: expand = %d, reject = %d, dE = %g, dV = %g\n", expand, reject, dE, _P*dV - (N + 1)*this->_T*log(V/oldV));
		//printf ("reject: dE = %g\n", dE);
		this->_Info->box->init(box_sides_old.x, box_sides_old.y, box_sides_old.z);
		for (int k = 0; k < N; k ++) {
			BaseParticle<number> *p = particles[k];
			p->pos = _pos_old[k];
			p->set_ext_potential(curr_step, this->_Info->box);
		}
		this->_Info->interaction->set_is_infinite (false);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}

	/*
	newE = this->_Info->interaction->get_system_energy(this->_Info->particles, *this->_Info->N, this->_Info->lists);
	if (this->_Info->interaction->get_is_infinite() == true) {
		printf ("here we go: just_accepted=%d, expand = %d, cut_at=%g, delta=%g\n", just_accepted, expand, cut_at, _delta);
		printf ("box_sides: %g %g %g\n", box_sides[0], box_sides[1], box_sides[2]);
		printf ("box_sides_old: %g %g %g\n", box_sides_old[0], box_sides_old[1], box_sides_old[2]);
		fflush (stdout);
		throw oxDNAException("NANA");
	}*/

	return;
}

template<typename number>
void CutVolume<number>::log_parameters() {
	BaseMove<number>::log_parameters();
	OX_LOG(Logger::LOG_INFO, "\tdelta = %g", _delta);
}

template class CutVolume<float>;
template class CutVolume<double>;
