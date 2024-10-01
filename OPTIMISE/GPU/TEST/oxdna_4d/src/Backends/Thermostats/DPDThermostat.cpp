/*
 * DPDThermostat.cpp
 *
 *  Created on: Feb 16, 2013
 *      Author: Petr
 */

#include "DPDThermostat.h"

#include "../../Utilities/Utils.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Lists/BaseList.h"
#include "../../Interactions/BaseInteraction.h"

using std::vector;

DPDThermostat::DPDThermostat() :
				BaseThermostat() {
	_reduced_mass = 0.5;
	_exponent = 1;
	_supports_shear = true;
	_zeta = -1.;
	_rcut = -1.;
	_dt = -1.;
	_sqrt_dt = -1.;
}

DPDThermostat::~DPDThermostat() {

}

void DPDThermostat::get_settings(input_file &inp) {
	BaseThermostat::get_settings(inp);

	getInputNumber(&inp, "DPD_zeta", &_zeta, 1);
	getInputNumber(&inp, "DPD_rcut", &_rcut, 1);
	getInputNumber(&inp, "DPD_reduced_mass", &_reduced_mass, 0);
	getInputNumber(&inp, "DPD_exponent", &_exponent, 0);
	getInputNumber(&inp, "dt", &_dt, 1);
	_sqrt_dt = sqrt(_dt);

	OX_LOG(Logger::LOG_INFO, "DPD thermostat parameters: zeta = %g, rcut = %g, exponent =  %g", _zeta, _rcut, _exponent);
}

void DPDThermostat::init() {
	BaseThermostat::init();

	number inter_rcut = CONFIG_INFO->interaction->get_rcut();
	if(inter_rcut < _rcut) {
		throw oxDNAException("DPD thermostat: the DPD cut-off (%lf) may not be larger than the interaction's cutoff (%lf)", _rcut, inter_rcut);
	}
}

void DPDThermostat::apply(std::vector<BaseParticle *> &particles, llint curr_step) {
	vector<ParticlePair> pairs = CONFIG_INFO->lists->get_potential_interactions();
	BaseBox *box = CONFIG_INFO->box;

	for(typename vector<ParticlePair>::iterator it = pairs.begin(); it != pairs.end(); it++) {
		BaseParticle *p = it->first;
		BaseParticle *q = it->second;

		// here we use the method proposed by Peters, Europhys. Lett. 66 (3), pp. 311-317 (2004)
		LR_vector rpq = box->min_image(q, p);
		number rpq_mod = rpq.module();
		if(rpq_mod < _rcut) {
			rpq /= rpq_mod;
			number omega_d = pow(1. - rpq_mod / _rcut, _exponent);
			number a_pq = _zeta * omega_d;
			number b_pq = sqrt(2 * this->_T * a_pq * (1. - a_pq * _dt / (2. * _reduced_mass)));

			LR_vector vpq = p->vel - q->vel;
			LR_vector delta_v = rpq * (-a_pq * _dt * (vpq * rpq) + b_pq * Utils::gaussian() * _sqrt_dt);
//			LR_vector delta_v = rpq*(-a_pq*_dt*(vpq*rpq) + b_pq*sqrt(12.)*(drand48() - 0.5)*_sqrt_dt);

			p->vel += delta_v;
			q->vel -= delta_v;
		}
	}
}
