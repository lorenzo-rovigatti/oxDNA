/*
 * StressAutocorrelation.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#include "StressAutocorrelation.h"

StressAutocorrelation::StressAutocorrelation() :
				BaseObservable() {

}

StressAutocorrelation::~StressAutocorrelation() {

}

void StressAutocorrelation::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	uint m = 2;
	uint p = 16;

	getInputUInt(&my_inp, "m", &m, 0);
	getInputUInt(&my_inp, "p", &p, 0);

	getInputDouble(&sim_inp, "dt", &_delta_t, 1);
	_delta_t *= _update_every;

	_sigma_xy = std::make_shared<Level>(m, p, 0);
	_sigma_yz = std::make_shared<Level>(m, p, 0);
	_sigma_zx = std::make_shared<Level>(m, p, 0);
	_N_xy = std::make_shared<Level>(m, p, 0);
	_N_yz = std::make_shared<Level>(m, p, 0);
	_N_xz = std::make_shared<Level>(m, p, 0);
}

void StressAutocorrelation::update_data(llint curr_step) {

	_config_info->interaction->begin_energy_computation();

	// pair_interaction will change these vectors, but we still need them in the next
	// first integration step. For this reason we copy and then restore their values
	// after the calculation
	_old_forces.resize(_config_info->N());
	_old_torques.resize(_config_info->N());
	for(int i = 0; i < _config_info->N(); i++) {
		_old_forces[i] = _config_info->particles()[i]->force;
		_old_torques[i] = _config_info->particles()[i]->torque;
	}

	StressTensor stress_tensor = { 0., 0., 0., 0., 0., 0. };
	double energy = 0.;

	for(auto p : _config_info->particles()) {
		std::vector<BaseParticle *> neighs = _config_info->lists->get_neigh_list(p);

		std::set<BaseParticle *> bonded_neighs;
		for(auto &pair : p->affected) {
			if(pair.first != p) {
				bonded_neighs.insert(pair.first);
			}
			if(pair.second != p) {
				bonded_neighs.insert(pair.second);
			}
		}

		neighs.insert(neighs.end(), bonded_neighs.begin(), bonded_neighs.end());

		for(auto q : neighs) {
			if(p->index > q->index) {
				LR_vector r = _config_info->box->min_image(p->pos, q->pos);
				_config_info->interaction->set_computed_r(r);

				p->force = LR_vector();
				energy += (double) _config_info->interaction->pair_interaction(p, q, false, true);

				stress_tensor[0] -= r.x * p->force.x;
				stress_tensor[1] -= r.y * p->force.y;
				stress_tensor[2] -= r.z * p->force.z;
				stress_tensor[3] -= r.x * p->force.y;
				stress_tensor[4] -= r.x * p->force.z;
				stress_tensor[5] -= r.y * p->force.z;
			}
		}
	}

	for(int i = 0; i < _config_info->N(); i++) {
		auto p = _config_info->particles()[i];

		p->force = _old_forces[i];
		p->torque = _old_torques[i];

		LR_vector vel = p->vel;

		stress_tensor[0] += SQR(vel.x);
		stress_tensor[1] += SQR(vel.y);
		stress_tensor[2] += SQR(vel.z);
		stress_tensor[3] += vel.x * vel.y;
		stress_tensor[4] += vel.x * vel.z;
		stress_tensor[5] += vel.y * vel.z;
	}

	double V = _config_info->box->V();
	for(auto &v : stress_tensor) {
		v /= V;
	}

	_sigma_xy->add_value(stress_tensor[3]);
	_sigma_yz->add_value(stress_tensor[5]);
	_sigma_zx->add_value(stress_tensor[4]);

	_N_xy->add_value(stress_tensor[0] - stress_tensor[1]);
	_N_xz->add_value(stress_tensor[0] - stress_tensor[2]);
	_N_yz->add_value(stress_tensor[1] - stress_tensor[2]);
}

std::string StressAutocorrelation::get_output_string(llint curr_step) {
	std::stringstream ss;

	std::vector<double> times;
	_sigma_xy->get_times(_delta_t, times);

	std::vector<double> acf_sigma_xy, acf_sigma_yz, acf_sigma_zx, acf_N_xy, acf_N_xz, acf_N_yz;
	_sigma_xy->get_acf(_delta_t, acf_sigma_xy);
	_sigma_yz->get_acf(_delta_t, acf_sigma_yz);
	_sigma_zx->get_acf(_delta_t, acf_sigma_zx);

	_N_xy->get_acf(_delta_t, acf_N_xy);
	_N_xz->get_acf(_delta_t, acf_N_xz);
	_N_yz->get_acf(_delta_t, acf_N_yz);

	double V = _config_info->box->V();
	double T = _config_info->temperature();
	for(uint i = 0; i < times.size(); i++) {
		double Gt = V / (5. * T) * (acf_sigma_xy[i] + acf_sigma_yz[i] + acf_sigma_zx[i]);
		Gt += V / (30. * T) * (acf_N_xy[i] + acf_N_xz[i] + acf_N_yz[i]);

		ss << times[i] << " " << Gt << std::endl;
	}

	return ss.str();
}
