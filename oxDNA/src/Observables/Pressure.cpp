/*
 * Pressure.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: lorenzo
 */

#include "Pressure.h"

template<typename number>
Pressure<number>::Pressure() :
				BaseObservable<number>(),
				_with_stress_tensor(false),
				_spherical(false),
				_nonbonded_only(false),
				_PV_only(false),
				_P(0.),
				_P_norm(0.),
				_shear_rate(0.) {

}

template<typename number>
Pressure<number>::~Pressure() {

}

template<typename number>
void Pressure<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable<number>::get_settings(my_inp, sim_inp);

	char raw_T[256];
	getInputString(&sim_inp, "T", raw_T, 1);
	_T = Utils::get_temperature<number>(raw_T);

	getInputBool(&my_inp, "stress_tensor", &_with_stress_tensor, 0);
	getInputBool(&my_inp, "spherical_geometry", &_spherical, 0);
	getInputBool(&my_inp, "nonbonded_only", &_nonbonded_only, 0);
	getInputBool(&my_inp, "PV_only", &_PV_only, 0);

	bool lees_edwards = false;
	getInputBool(&sim_inp, "lees_edwards", &lees_edwards, 0);
	if(lees_edwards) getInputNumber(&sim_inp, "lees_edwards_shear_rate", &_shear_rate, 0);
}

template<typename number>
LR_vector<number> Pressure<number>::_get_com() {
	LR_vector<number> com;
	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		com += this->_config_info.box->get_abs_pos(p);
	}
	com /= *this->_config_info.N;
	return com;
}

template<typename number>
void Pressure<number>::update_pressure() {
	int N = *this->_config_info.N;
	vector<ParticlePair<number> > pairs = this->_config_info.lists->get_potential_interactions();

	LR_vector<number> com = (_spherical) ? _get_com() : LR_vector<number>();
	_P_norm = 0.;

	double virial = 0;
	_stress_tensor = LR_matrix<double>();
	double energy = 0.;
	// we loop over all the pairs in order to update the forces
	typename vector<ParticlePair<number> >::iterator it;
	for(it = pairs.begin(); it != pairs.end(); it++) {
		BaseParticle<number> *p = (*it).first;
		BaseParticle<number> *q = (*it).second;
		if(!(_nonbonded_only && p->is_bonded(q))) {
			LR_vector<number> r = this->_config_info.box->min_image(p->pos, q->pos);

			// pair_interaction will change these vectors, but we still need them in the next
			// first integration step. For this reason we copy and then restore their values
			// after the calculation
			LR_vector<number> old_p_force(p->force);
			LR_vector<number> old_q_force(q->force);
			LR_vector<number> old_p_torque(p->torque);
			LR_vector<number> old_q_torque(q->torque);

			p->force = q->force = p->torque = q->torque = LR_vector<number>();

			LR_vector<number> r_mutable(r);
			energy += (double) this->_config_info.interaction->pair_interaction(p, q, &r_mutable, true);

			_stress_tensor.v1.x -= r.x * p->force.x;
			_stress_tensor.v1.y -= r.x * p->force.y;
			_stress_tensor.v1.z -= r.x * p->force.z;
			_stress_tensor.v2.x -= r.y * p->force.x;
			_stress_tensor.v2.y -= r.y * p->force.y;
			_stress_tensor.v2.z -= r.y * p->force.z;
			_stress_tensor.v3.x -= r.z * p->force.x;
			_stress_tensor.v3.y -= r.z * p->force.y;
			_stress_tensor.v3.z -= r.z * p->force.z;

			virial -= (r * p->force);

			// see http://dx.doi.org/10.1080/102866202100002518a
			if(_spherical) {
				LR_vector<number> rp_rel = this->_config_info.box->min_image(p->pos, com);
				LR_vector<number> rq_rel = this->_config_info.box->min_image(q->pos, com);
				number l0 = sqrt(rp_rel.norm() - SQR(rp_rel*r)/r.norm());
				// check that the arguments of the square roots are non-negative (it does happen for finite-accuracy reasons)
				double delta_rp_sqr_l0_sqr = rp_rel.norm() - SQR(l0);
				double F0 = (delta_rp_sqr_l0_sqr < 0.) ? 0 : -l0*atan(sqrt(delta_rp_sqr_l0_sqr)/l0);
				double delta_rq_sqr_l0_sqr = rq_rel.norm() - SQR(l0);
				double F1 = (delta_rq_sqr_l0_sqr < 0.) ? 0 : -l0*atan(sqrt(delta_rq_sqr_l0_sqr)/l0);
				double Fpq_mod = -(p->force*r)/r.module();
				double rpq_mod = r.module();
				_P_norm += Fpq_mod*(rpq_mod + F1 - F0);
			}

			p->force = old_p_force;
			q->force = old_q_force;
			p->torque = old_p_torque;
			q->torque = old_q_torque;
		}
	}

	for(int i = 0; i < *this->_config_info.N; i++) {
		BaseParticle<number> *p = this->_config_info.particles[i];
		LR_vector<number> vel = p->vel;
		if(_shear_rate > 0.) {
			number Ly = CONFIG_INFO->box->box_sides().y;
			number y_in_box = p->pos.y - floor(p->pos.y/Ly)*Ly - 0.5*Ly;
			number flow_vx = y_in_box*_shear_rate;
			vel.x -= flow_vx;
		}
		_stress_tensor.v1.x += SQR(vel.x);
		_stress_tensor.v1.y += vel.x * vel.y;
		_stress_tensor.v1.z += vel.x * vel.z;
		_stress_tensor.v2.x += vel.y * vel.x;
		_stress_tensor.v2.y += SQR(vel.y);
		_stress_tensor.v2.z += vel.y * vel.z;
		_stress_tensor.v3.x += vel.z * vel.x;
		_stress_tensor.v3.y += vel.z * vel.y;
		_stress_tensor.v3.z += SQR(vel.z);

//		if(_spherical) {
//			LR_vector<number> rp_rel = this->_config_info.box->min_image(p->pos, com);
//			rp_rel.normalize();
//			_P_norm += 3*SQR(vel*rp_rel);
//			printf("%e %e\n", SQR(vel*rp_rel), vel.norm());
//		}
	}

	double V = (_PV_only) ? 1 : this->_config_info.box->V();
	_P = _T * (N / V) + virial / (3. * V);
	_stress_tensor /= 3. * V;
	_P_norm = _T * (N / V) + _P_norm / (3. * V);
}

template<typename number>
string Pressure<number>::get_output_string(llint curr_step) {
	update_pressure();

	string to_ret;

	if(_with_stress_tensor) {
		to_ret += Utils::sformat("% .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e % .8e", get_P(), _stress_tensor.v1.x, _stress_tensor.v1.y, _stress_tensor.v1.z, _stress_tensor.v2.x, _stress_tensor.v2.y, _stress_tensor.v2.z, _stress_tensor.v3.x, _stress_tensor.v3.y, _stress_tensor.v3.z);
	}
	else {
		to_ret += Utils::sformat("% .8e", get_P());
	}

	if(_spherical) {
		to_ret += Utils::sformat(" % .8e", _P_norm);
	}

	return to_ret;
}

template class Pressure<float> ;
template class Pressure<double> ;
