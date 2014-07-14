#include "DNAInteraction2.h"

#include "../Particles/DNANucleotide.h"


template<typename number>
DNA2Interaction<number>::DNA2Interaction() : DNAInteraction<number>() {
	this->_int_map[DEBYE_HUCKEL] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_debye_huckel;

	_debye_huckel_half_charged_ends = true;
	this->_grooving = true;
}

template<typename number>
number DNA2Interaction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

        if (r->norm() >= this->_sqr_rcut) return (number) 0.f;

	number energy = DNAInteraction<number>::pair_interaction_nonbonded(p, q, r, update_forces);
	energy += _debye_huckel(p, q, r, update_forces);

	return energy;
}

template<typename number>
void DNA2Interaction<number>::get_settings(input_file &inp) {
	DNAInteraction<number>::get_settings(inp);

	getInputNumber(&inp, "salt_concentration", &_salt_concentration, 1);
	OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel at salt concentration =  %g", this->_salt_concentration);

        getInputBool(&inp, "dh_half_charged_ends", &_debye_huckel_half_charged_ends, 0);
	OX_LOG(Logger::LOG_INFO,"dh_half_charged_ends = %s", _debye_huckel_half_charged_ends ? "true" : "false");

	// lambda-factor (the dh length at T = 300K, I = 1.0)
	float lambdafactor;
	if (getInputFloat(&inp, "dh_lambda", &lambdafactor, 0) == KEY_FOUND) {
		_debye_huckel_lambdafactor = (float) lambdafactor;
	} 
	else {
		_debye_huckel_lambdafactor = 0.3616455;
	}

	// the prefactor to the Debye-Huckel term
	float prefactor;
	if (getInputFloat(&inp, "dh_strength", &prefactor, 0) == KEY_FOUND) {
		_debye_huckel_prefactor = (float) prefactor;
	} 
	else {
		_debye_huckel_prefactor = 0.0543;
	}

	// notify the user that major-minor grooving is switched on
	// check whether it's set in the input file to avoid duplicate messages
	int tmp;
	if (this->_grooving && (getInputBoolAsInt(&inp, "major_minor_grooving", &tmp, 0) != KEY_FOUND)) OX_LOG(Logger::LOG_INFO, "Using different widths for major and minor grooves");
}

template<typename number>
void DNA2Interaction<number>::init() {
	DNAInteraction<number>::init();

	// We wish to normalise with respect to T=300K, I=1M. 300K=0.1 s.u. so divide this->_T by 0.1
	number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
	// RHIGH gives the distance at which the smoothing begins
	_debye_huckel_RHIGH = 3.0 * lambda;
	_minus_kappa = -1.0/lambda;
	
	// these are just for convenience for the smoothing parameter computation
	number x = _debye_huckel_RHIGH;
	number q = _debye_huckel_prefactor;
	number l = lambda;
	
	// compute the some smoothing parameters
	_debye_huckel_B = -(exp(-x/l) * q * q * (x + l)*(x+l) )/(-4.*x*x*x * l * l * q );
	_debye_huckel_RC = x*(q*x + 3. * q* l )/(q * (x+l));

	number debyecut;
	if (this->_grooving){
		debyecut = 2.0f * sqrt((POS_MM_BACK1)*(POS_MM_BACK1) + (POS_MM_BACK2)*(POS_MM_BACK2)) + _debye_huckel_RC;
	}
	else{
		debyecut =  2.0f * sqrt(SQR(POS_BACK)) + _debye_huckel_RC;
	}
	// the cutoff radius for the potential should be the larger of rcut and debyecut
	if (debyecut > this->_rcut){
		this->_rcut = debyecut;
		this->_sqr_rcut = debyecut*debyecut;
	}

	// NB lambda goes into the exponent for the D-H potential and is given by lambda = lambda_k * sqrt((T/300K)/(I/1M))
	OX_LOG(Logger::LOG_DEBUG,"Debye-Huckel parameters: Q=%f, lambda_0=%f, lambda=%f, r_high=%f, cutoff=%f", _debye_huckel_prefactor, _debye_huckel_lambdafactor, lambda, _debye_huckel_RHIGH, this->_rcut);
	OX_LOG(Logger::LOG_INFO,"The Debye length at this temperature and salt concentration is %f", lambda);

}

template<typename number>
number DNA2Interaction<number>::_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)
{
	number cut_factor = 1.0f;
	if(this->_are_bonded(p, q)) return (number) 0.f;
	// for each particle that is on a terminus, halve the charge
	if (this->_debye_huckel_half_charged_ends && (p->n3 == P_VIRTUAL || p->n5 == P_VIRTUAL)) cut_factor *= 0.5f;
	if (this->_debye_huckel_half_charged_ends && (q->n3 == P_VIRTUAL || q->n5 == P_VIRTUAL)) cut_factor *= 0.5f;
	
	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number energy = (number) 0.f;
	
	// Debye-Huckel energy
	if (rbackmod <_debye_huckel_RC) {
		if(rbackmod < _debye_huckel_RHIGH) {
			energy =   exp(rbackmod * _minus_kappa) * (  _debye_huckel_prefactor / rbackmod );
		}
		else {
			energy = _debye_huckel_B * SQR(rbackmod -  _debye_huckel_RC);
		}

		energy *= cut_factor; 

		if(update_forces && energy != 0.) {
			LR_vector<number> force(0.,0.,0.);
			LR_vector<number> torqueq(0.,0.,0.);
			LR_vector<number> torquep(0.,0.,0.);
			LR_vector<number> rbackdir = rback / rbackmod;
			
			if(rbackmod < _debye_huckel_RHIGH){
				force =  rbackdir *  ( -1.0f *  (_debye_huckel_prefactor * exp(_minus_kappa * rbackmod) ) *  (  _minus_kappa / rbackmod   - 1.0f /SQR(rbackmod) ) );
			}
			else {
				force = - rbackdir * (2.0f * _debye_huckel_B * (rbackmod -  _debye_huckel_RC) );
			}
			force *= cut_factor;
			
			torqueq = q->int_centers[DNANucleotide<number>::BACK].cross(force);
			torquep = -p->int_centers[DNANucleotide<number>::BACK].cross(force);
			
			p->force -= force;
			q->force += force;
			
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
			
		}
	}
	
	return energy;
}

template class DNA2Interaction<float>;
template class DNA2Interaction<double>;

