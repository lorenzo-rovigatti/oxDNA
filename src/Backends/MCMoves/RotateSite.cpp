/**
 * @file    RotateSite.cpp
 * @date    23/mar/2016
 * @author  flavio
 *
 *
 */

#include "RotateSite.h"
#include "../../Particles/JordanParticle.h"

template<typename number>
RotateSite<number>::RotateSite (ConfigInfo<number> * Info) : BaseMove<number>(Info) {

}

template<typename number>
RotateSite<number>::~RotateSite () {

}

template<typename number>
void RotateSite<number>::get_settings (input_file &inp, input_file &sim_inp) {
	BaseMove<number>::get_settings (inp, sim_inp);

	getInputNumber (&inp, "delta", &_delta, 1);
	getInputNumber (&inp, "prob", &this->prob, 1);

	OX_LOG(Logger::LOG_INFO, "(RotateSite.cpp) RotateSite initiated with T %g, delta %g, prob: %g", this->_T, _delta, this->prob);
}

template<typename number>
void RotateSite<number>::apply (llint curr_step) {

	this->_attempted ++;

	int pi = (int) (drand48() * (*this->_Info->N));
	JordanParticle<number> *p = (JordanParticle<number> *)this->_Info->particles[pi];

	number delta_E = -this->particle_energy(p);
	delta_E -= p->int_potential();
	p->set_ext_potential (curr_step, (*this->_Info->box_side));
	number delta_E_ext = -p->ext_potential;

	// select site
	int i_patch = (int) (drand48() * (p->N_int_centers));
	LR_matrix<number> site_store = p->get_patch_rotation(i_patch);
	
	number t = drand48() * _delta;
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
	
	p->rotate_patch(i_patch, R);
	p->set_positions();
	
	if (p->is_rigid_body()) {
		this->_Info->lists->single_update(p);
		if(!this->_Info->lists->is_updated()) {
			this->_Info->lists->global_update();
		}
	}
	
	delta_E += this->particle_energy(p);
	delta_E += p->int_potential();
	p->set_ext_potential(curr_step, *this->_Info->box_side);
	delta_E_ext += p->ext_potential;
	
	// accept or reject?
	if (this->_Info->interaction->get_is_infinite() == false && ((delta_E + delta_E_ext) < 0 || exp(-(delta_E + delta_E_ext) / this->_T) > drand48() )) {
		// move accepted
		// put here the adjustment of moves
		this->_accepted ++;
		if (curr_step < this->_equilibration_steps && this->_adjust_moves) {
			_delta *= this->_acc_fact;
			if (_delta > M_PI) _delta = M_PI;
		}
	}
	else {
		p->set_patch_rotation(i_patch, site_store);
		p->set_positions();
	
		if (p->is_rigid_body()) {
			this->_Info->lists->single_update(p);
			if(!this->_Info->lists->is_updated()) {
				this->_Info->lists->global_update();
			}
		}
		this->_Info->interaction->set_is_infinite(false);

		if (curr_step < this->_equilibration_steps && this->_adjust_moves) _delta /= this->_rej_fact;
	}
	
	return;
}

template class RotateSite<float>;
template class RotateSite<double>;

