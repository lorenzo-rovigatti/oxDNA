/*
 * FakePressure.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: flavio
 */

#include "FakePressure.h"

template<typename number> FakePressure<number>::FakePressure() {
    _cells = NULL;
    _particles = NULL;
    _sigma = -1.f;
    _restrict_to_type = -1;
}

template<typename number>
FakePressure<number>::~FakePressure() {

}

template<typename number>
void FakePressure<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);

    OX_LOG (Logger::LOG_INFO, "(FakePressure.cpp) Using _sigma %g and restrict_to_type=%d", _sigma, _restrict_to_type);

}

template<typename number>
void FakePressure<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
    getInputNumber (&my_inp, "sigma", &_sigma, 1);
    getInputInt(&my_inp, "restrict_to_type", &_restrict_to_type, 1);
    //if (_restrict_to_type < 0 ) throw oxDNAException ("(FakePressure.cpp) Cannot run with restrict_to_type=%d. Needs to be positive")
	}

template<typename number>
std::string FakePressure<number>::get_output_string(llint curr_step) {

    std::string ret ("");

    int N_of_type = 0;
    // build particles and cells
    int N = *this->_config_info.N;
    _particles = new BaseParticle<number> * [N];
    for(int i = 0; i < N; i++) {
        _particles[i] = this->_config_info.particles[i];
        if (_particles[i]->type == _restrict_to_type) N_of_type ++;
    }

    number rcut = (number) 2.5f * _sigma;
    _cells = new Cells<number>(*this->_config_info.N, this->_config_info.box);
    _cells->init(_particles, rcut);

    std::vector<ParticlePair<number> > pairs = _cells->get_potential_interactions();
    typename std::vector<ParticlePair<number> >::iterator it;
    number pressure = (number) 0.f;
    for (it = pairs.begin(); it != pairs.end(); it ++) {
        BaseParticle<number> * p = (*it).first;
        BaseParticle<number> * q = (*it).second;

        if (p->type != _restrict_to_type || q->type != _restrict_to_type)
            continue;

        LR_vector<number> r = this->_config_info.box->min_image(p->pos, q->pos);
        number rnorm = r.norm();
        if (rnorm < (2.5*_sigma * 2.5 * _sigma)) {
            LR_vector<number> force_on_p = -((number)6.f * powf(_sigma, 6) / powf(rnorm, 4)) * r;
            //LR_vector<number> force_on_q = -force_on_p;
            pressure += -force_on_p * r; // should the minus be there?
        }
    }

    number volume = this->_config_info.box->V();

    //ret += Utils::sformat("%g ", pressure);

    // beta * P = rho + (beta / (3V)) * Sum(f_ji * r_ij)
    pressure = 1.f * (N_of_type / volume) + pressure / (3.f * volume);

    ret += Utils::sformat("%g", pressure);

    delete[] _particles;
    delete _cells;

    return ret;
}

template class FakePressure<float>;
template class FakePressure<double>;
