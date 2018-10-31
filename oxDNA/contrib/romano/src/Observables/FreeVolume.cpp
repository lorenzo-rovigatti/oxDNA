/*
 * FreeVolume.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: flavio
 */

#include "FreeVolume.h"

template<typename number> FreeVolume<number>::FreeVolume() {
    _restrict_to_type = -1;
    _sigma = 0.5f;
    _ntries = 0;
}

template<typename number>
FreeVolume<number>::~FreeVolume() {

}

template<typename number>
void FreeVolume<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);
    if (_ntries < 0) {
        _ntries = -_ntries * int(this->_config_info.box->V());
    }
    OX_LOG (Logger::LOG_INFO, "(FreeVolume.cpp) Using sigma=%g, ntries=%d (%g per unit volume) and restrict_to_type=%d", _sigma, _ntries, _ntries / this->_config_info.box->V(), _restrict_to_type);
}

template<typename number>
void FreeVolume<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
    getInputNumber (&my_inp, "sigma", &_sigma, 0);
    getInputInt(&my_inp, "restrict_to_type", &_restrict_to_type, 0);
    getInputInt (&my_inp, "ntries", &_ntries, 1);
}

template<typename number>
std::string FreeVolume<number>::get_output_string(llint curr_step) {

    std::string ret ("");

    LR_vector<number> box_sides = this->_config_info.box->box_sides();

    BaseParticle<number> * p;

    int N = *this->_config_info.N;
    int N_of_type = 0;
    for(int i = 0; i < N; i++) {
        p = this->_config_info.particles[i];
        if (_restrict_to_type < 0 || p->type == _restrict_to_type) {
            N_of_type ++;
        }
    }

    // build particles and cells
    BaseParticle<number> ** particles = new BaseParticle<number> * [N_of_type + 1]; // we ignore other types types; some of the memory won't be used
    int j = 0;
    for(int i = 0; i < N; i++) {
        p = this->_config_info.particles[i];
        if (_restrict_to_type <=0 || p->type == _restrict_to_type) {
            particles[j] = new BaseParticle<number> ();
            particles[j]->pos = p->pos;
            particles[j]->index = j;
            particles[j]->type = p->type;
            particles[j]->orientation = p->orientation;
            particles[j]->orientationT = p->orientationT;
            particles[j]->set_positions();
            j ++;
        }
    }

    // we build the ``test'' particle and add it to the array
    particles[N_of_type] = new BaseParticle<number> ();
    particles[N_of_type]->index = N_of_type;
    particles[N_of_type]->pos.x = drand48() * box_sides[0];
    particles[N_of_type]->pos.y = drand48() * box_sides[1];
    particles[N_of_type]->pos.z = drand48() * box_sides[2];
    if (_restrict_to_type >= 0) particles[N_of_type]->type = _restrict_to_type;
    else particles[N_of_type]->type = 1;

    p = particles[N_of_type];
    int real_N_of_type = N_of_type + 1;
    Cells<number> * cells = new Cells<number>(real_N_of_type, this->_config_info.box);
    cells->init(particles, 10 + 2.01 * _sigma);

    int ntries = 0;
    int nfull = 0;
    while (ntries < _ntries) {
        p->pos.x = box_sides[0] * drand48();
        p->pos.y = box_sides[1] * drand48();
        p->pos.z = box_sides[2] * drand48();
        cells->single_update(p);
        std::vector<BaseParticle<number> *> neighs = cells->get_complete_neigh_list(p);
        typename std::vector<BaseParticle<number> *>::iterator it;
        bool found = false;
        for (it = neighs.begin(); it != neighs.end() && found == false; it ++) {
            if (_restrict_to_type >= 0 && (*it)->type != _restrict_to_type) {
                //throw oxDNAException("fsf");
                continue;
            }
            //if (p->type != 1) throw oxDNAException("gg");
            //if ((*it)->type != 0) throw oxDNAException("hh %d", (*it)->index);
            //if ((*it)->index >= p->index) throw oxDNAException ("jj");
            this->_config_info.interaction->pair_interaction(*it, p);
            if (this->_config_info.interaction->get_is_infinite() == true) {
                found = true;
                this->_config_info.interaction->set_is_infinite(false);
            }
        }
        if (found) nfull ++;
        ntries ++;
    }

    int nempty = _ntries - nfull;
    number free_volume = (nempty / (double) _ntries) * this->_config_info.box->V();

    ret += Utils::sformat("%g %g", free_volume, (nempty / (double) _ntries));

    for (int i = 0; i <= N_of_type; i ++) delete particles[i];
    delete[] particles;
    delete cells;

    return ret;
}

template class FreeVolume<float>;
template class FreeVolume<double>;
