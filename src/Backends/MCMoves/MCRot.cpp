/**
 * @file    MCRot.cpp
 * @date    30/apr/2014
 * @author  flavio
 *
 *
 */

#include "MCRot.h"

template<typename number>
MCRot<number>::MCRot (ConfigInfo<number> * Info) : BaseMove<number>(Info) {
	_orientation_old = LR_matrix<number> (0., 0., 0., 0., 0., 0., 0., 0., 0.); 
	_orientationT_old = LR_matrix<number> (0., 0., 0., 0., 0., 0., 0., 0., 0.); 
}

template<typename number>
MCRot<number>::~MCRot () {

}

template<typename number>
void MCRot<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta_rot", &_delta, 1);
	getInputNumber (&inp, "prob", &this->prob, 1);

	OX_LOG(Logger::LOG_INFO, "(MCRot.cpp) MCRot initiated with T %g, delta_tras %g, prob: %g", this->_T, _delta, this->prob);
}

template<typename number>
void MCRot<number>::apply (llint curr_step) {

	this->_attempted ++;

	int pi = (int) (drand48() * (*this->_Info->N));
	BaseParticle<number> *p = this->_Info->particles[pi];

	number delta_E = -this->particle_energy(p);
	p->set_ext_potential (curr_step, (*this->_Info->box_side));
	number delta_E_ext = -p->ext_potential;

	_orientation_old = p->orientation;
	_orientationT_old = p->orientationT;
				
	number t = (drand48() - (number)0.5f) * _delta;
	LR_vector<number> axis = Utils::get_random_vector<number>();
	
	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number)1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	p->orientation = p->orientation * R;
	p->orientationT = p->orientation.get_transpose();
	p->set_positions();
	
	this->_Info->lists->single_update(p);

	if(!this->_Info->lists->is_updated()) {
		this->_Info->lists->global_update();
	}

	delta_E += this->particle_energy(p);
	p->set_ext_potential(curr_step, *this->_Info->box_side);
	delta_E_ext += p->ext_potential;
	
	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48() )) {
		// move accepted
		// put here the adjustment of moves
		this->_accepted ++;
	}
	else {
		p->orientation = _orientation_old;
		p->orientationT = _orientationT_old;
		p->set_positions();
		this->_Info->lists->single_update(p);
		this->_Info->interaction->set_is_infinite(false);
		// adjust moves here as well...
	}
	
	return;
}

template class MCRot<float>;
template class MCRot<double>;

