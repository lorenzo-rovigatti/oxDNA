/*
 * ConstantRateTorque.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "ConstantRateTorque.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"

template<typename number>
ConstantRateTorque<number>::ConstantRateTorque() : BaseForce<number> () {
	
}

template<typename number>
void ConstantRateTorque<number>::get_settings (input_file &inp) {
	getInputInt (&inp, "particle", &this->_particle, 1);
	
	getInputNumber(&inp, "stiff", &this->_stiff, 1);
	getInputNumber(&inp, "rate", &this->_rate, 1);
	getInputNumber(&inp, "base", &this->_F0, 1);

	std::string (strdir);
	double tmpf[3];
	int tmpi;
	getInputString (&inp, "axis", strdir, 1);
	tmpi = sscanf (strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf +2);
	if (tmpi != 3) throw oxDNAException ("Could not parse axis `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
	this->_axis = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	this->_axis.normalize();
	
	getInputString (&inp, "pos0", strdir, 1);
	tmpi = sscanf (strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf +2);
	if (tmpi != 3) throw oxDNAException ("Could not parse pos0 `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
	this->_pos0 = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
	
	getInputString (&inp, "center", strdir, 1);
	tmpi = sscanf (strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf +2);
	if (tmpi != 3) throw oxDNAException ("Could not parse center `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
	this->_center = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);

	if (getInputString (&inp, "mask", strdir, 0) == KEY_FOUND) {
		tmpi = sscanf (strdir.c_str(), "%lf, %lf, %lf", tmpf, tmpf + 1, tmpf + 2);
		if (tmpi != 3) throw oxDNAException("Could not parse mask `%s\' for ConstantRateTorque. Aborting", strdir.c_str());
		this->_mask = LR_vector<number> ((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]); 
	}
	else this->_mask = LR_vector<number> (0., 0., 0.);
}

template <typename number>
void ConstantRateTorque<number>::init (BaseParticle<number> ** particles, int N, number * box_side_ptr) {
	if (_particle >= N || N < -1) throw oxDNAException ("Trying to add a ConstantRateTorque on non-existent particle %d. Aborting", _particle);
	if (_particle != -1) {
		OX_LOG (Logger::LOG_INFO, "Adding ConstantRateTorque (F0==%g, rate=%g, pos0=%g,%g,%g, axis=%g,%g,%g, center=%g,%g,%g, mask=%g,%g,%g on particle %d", this->_F0, this->_pos0.x, this->_pos0.y, this->_pos0.z, this->_rate, this->_axis.x, this->_axis.y, this->_axis.z, this->_center.x, this->_center.y, this->_center.z, this->_mask.x, this->_mask.y, this->_mask.z, _particle);
		particles[_particle]->add_ext_force(this);
	}
	else { // force affects all particles
		OX_LOG (Logger::LOG_INFO, "Adding ConstantRateTorque (F0==%g, rate=%g, pos0=%g,%g,%g, axis=%g,%g,%g, center=%g,%g,%g, mask=%g,%g,%g on ALL particles", this->_F0, this->_pos0.x, this->_pos0.y, this->_pos0.z, this->_rate, this->_axis.x, this->_axis.y, this->_axis.z, this->_center.x, this->_center.y, this->_center.z, this->_mask.x, this->_mask.y, this->_mask.z);
		for (int i = 0; i < N; i ++) particles[i]->add_ext_force(this);
	}
}


template<typename number>
number ConstantRateTorque<number>::potential(llint step, LR_vector<number> &pos) {
	number t = (this->_F0 + this->_rate * (number) step);

	number sintheta = sin (t);
	number costheta = cos (t);
	number olcos = ((number) 1.) - costheta;

	number xyo = this->_axis.x * this->_axis.y * olcos;
	number xzo = this->_axis.x * this->_axis.z * olcos;
	number yzo = this->_axis.y * this->_axis.z * olcos;
	number xsin = this->_axis.x * sintheta;
	number ysin = this->_axis.y * sintheta;
	number zsin = this->_axis.z * sintheta;

	LR_matrix<number> R(this->_axis.x * this->_axis.x * olcos + costheta,
						xyo - zsin,
						xzo + ysin,

						xyo + zsin,
						this->_axis.y * this->_axis.y * olcos + costheta,
						yzo - xsin,

						xzo - ysin,
						yzo + xsin,
						this->_axis.z * this->_axis.z * olcos + costheta);

	LR_vector<number> postrap = R * (_pos0 - _center) + _center;

	// we "mask" the resulting vector;
	LR_vector<number> dr = pos - postrap;
	dr.x *= _mask.x;
	dr.y *= _mask.y;
	dr.z *= _mask.z;

	return (number) (0.5 * this->_stiff * (dr * dr));
}


template<typename number>
LR_vector<number> ConstantRateTorque<number>::value(llint step, LR_vector<number> &pos) {
	//
	// dobbiamo ruotare pos0 di (base + rate * t) radianti
	// intorno all'asse passante per il centro.
	// e quello ci da la pos della trappola;
	//
	number t = (this->_F0 + this->_rate * (number) step);

	number sintheta = sin (t);
	number costheta = cos (t);
	number olcos = ((number) 1.) - costheta;

	number xyo = this->_axis.x * this->_axis.y * olcos;
	number xzo = this->_axis.x * this->_axis.z * olcos;
	number yzo = this->_axis.y * this->_axis.z * olcos;
	number xsin = this->_axis.x * sintheta;
	number ysin = this->_axis.y * sintheta;
	number zsin = this->_axis.z * sintheta;

	LR_matrix<number> R(this->_axis.x * this->_axis.x * olcos + costheta,
						xyo - zsin,
						xzo + ysin,

						xyo + zsin,
						this->_axis.y * this->_axis.y * olcos + costheta,
						yzo - xsin,

						xzo - ysin,
						yzo + xsin,
						this->_axis.z * this->_axis.z * olcos + costheta);

	LR_vector<number> postrap = R * (_pos0 - _center) + _center;

	// we "mask" the resulting vector;
	number x = - this->_stiff * (pos.x - postrap.x) * _mask.x;
	number y = - this->_stiff * (pos.y - postrap.y) * _mask.y;
	number z = - this->_stiff * (pos.z - postrap.z) * _mask.z;

	return LR_vector<number>(x, y, z);
}

template class ConstantRateTorque<double>;
template class ConstantRateTorque<float>;
