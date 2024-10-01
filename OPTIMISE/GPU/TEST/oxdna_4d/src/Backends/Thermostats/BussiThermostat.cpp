/*
 * BussiThermostat.cpp
 *
 *  Created on: May 07, 2015
 *      Author: Flavio
 */

#include "BussiThermostat.h"
#include "../../Utilities/Utils.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Boxes/BaseBox.h"

BussiThermostat::BussiThermostat() :
				BaseThermostat() {
	_newtonian_steps = 0;
	_tau = -1;
	_exp_dt_tau = 1;
	_K_t = 0.;
	_K_r = 0.;
	_supports_shear = true;
}

BussiThermostat::~BussiThermostat() {

}

void BussiThermostat::get_settings(input_file &inp) {
	BaseThermostat::get_settings(inp);
	getInputInt(&inp, "newtonian_steps", &_newtonian_steps, 1);
	getInputInt(&inp, "bussi_tau", &_tau, 1);
	if(_newtonian_steps < 1) {
		throw oxDNAException("'newtonian_steps' must be > 0");
	}
}

void BussiThermostat::init() {
	BaseThermostat::init();

	_K_t = 0.5 * _current_translational_degrees_of_freedom() * _T;
	_K_r = 0.5 * _current_rotational_degrees_of_freedom() * _T;

	_exp_dt_tau = exp(-_newtonian_steps / (number) _tau);
}

void BussiThermostat::_update_K(number &K, int degrees_of_freedom) {
	// dynamics for the kinetic energy
	number K_target = 0.5 * degrees_of_freedom * _T;

	number rr = Utils::gaussian();
	number dK = (1.0 - _exp_dt_tau) * (K_target * (_sum_noises(degrees_of_freedom - 1) + rr * rr) / degrees_of_freedom - K) + 2.0 * rr * std::sqrt(K * K_target / degrees_of_freedom * (1.0 - _exp_dt_tau) * _exp_dt_tau);
	K += dK;
}

void BussiThermostat::apply(std::vector<BaseParticle *> &particles, llint curr_step) {
	if(curr_step % _newtonian_steps) {
		return;
	}

	LR_vector v_com;
	for(auto p: particles) {
		v_com += p->vel;
	}
	v_com /= CONFIG_INFO->N();

	// compute the total kinetic energy
	number K_now_t = 0.;
	number K_now_r = 0.;
	for(auto p: particles) {
		LR_vector v_rel = p->vel - v_com;
		if(_lees_edwards) {
			// we compute the instantaneous kinetic energy considering only two out of three dimensions, leaving out the flow direction x
			K_now_t += (SQR(v_rel.z) + SQR(v_rel.y)) * 3. / 4.;
		}
		else {
			K_now_t += 0.5 * (v_rel * v_rel);
		}
		K_now_r += 0.5 * (p->L * p->L);
	}

	_update_K(_K_t, _current_translational_degrees_of_freedom());
	_update_K(_K_r, _current_rotational_degrees_of_freedom());

	number rescale_factor_t = std::sqrt(_K_t / K_now_t);
	number rescale_factor_r = std::sqrt(_K_r / K_now_r);

	for(auto p: particles) {
		if(_lees_edwards) {
			number Ly = CONFIG_INFO->box->box_sides().y;
			number y_in_box = p->pos.y - floor(p->pos.y / Ly) * Ly - 0.5 * Ly;
			number flow_vx = y_in_box * _shear_rate;
			p->vel.x = (p->vel.x - flow_vx) * rescale_factor_t + flow_vx;
			p->vel.y = (p->vel.y - v_com.y) * rescale_factor_t + v_com.y;
			p->vel.z = (p->vel.z - v_com.z) * rescale_factor_t + v_com.z;
		}
		else {
			p->vel = (p->vel - v_com) * rescale_factor_t + v_com;
		}
		if(p->is_rigid_body()) {
			p->L *= rescale_factor_r;
		}
	}
}

// Bussi's methods
number BussiThermostat::_sum_noises(int nn) {
	number rr;
	if(nn == 0) {
		return 0.0;
	}
	else if(nn == 1) {
		rr = Utils::gaussian();
		return rr * rr;
	}
	else if(nn % 2 == 0) {
		return 2.0 * _gamdev(nn / 2);
	}
	else {
		rr = Utils::gaussian();
		return 2.0 * _gamdev((nn - 1) / 2) + SQR(rr);
	}
}

number BussiThermostat::_gamdev(int ia) {
	int j;
	number am, e, s, v1, v2, x, y;

	if(ia < 1) throw oxDNAException("BussiThermostat::_gamdev may not be called with an argument ia < 1");
	if(ia < 6) {
		x = 1.0;
		for(j = 1; j <= ia; j++)
			x *= drand48();
		x = -std::log(x);
	}
	else {
		do {
			do {
				do {
					v1 = drand48();
					v2 = 2.0 * drand48() - 1.0;
				} while(SQR(v1) + SQR(v2) > 1.0);
				y = v2 / v1;
				am = ia - 1;
				s = std::sqrt(2.0 * am + 1.0);
				x = s * y + am;
			} while(x <= 0.0);
			e = (1.0 + SQR(y)) * std::exp(am * std::log(x / am) - s * y);
		} while(drand48() > e);
	}
	return x;
}

int BussiThermostat::_current_translational_degrees_of_freedom() {
	return 3. * (CONFIG_INFO->N() - 1);
}

int BussiThermostat::_current_rotational_degrees_of_freedom() {
	return 3. * CONFIG_INFO->N();
}
