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
    OX_LOG (Logger::LOG_INFO, "(FreeVolume.cpp) Using sigma=%g, and restrict_to_type=%d", _sigma, _restrict_to_type);
}

template<typename number>
void FreeVolume<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
    getInputNumber (&my_inp, "sigma", &_sigma, 0);
    getInputInt(&my_inp, "restrict_to_type", &_restrict_to_type, 1);
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
        if (p->type == _restrict_to_type) {
            N_of_type ++;
        }
    }

    // build particles and cells
    BaseParticle<number> ** particles = new BaseParticle<number> * [N_of_type + 1]; // we ignore other types types; some of the memory won't be used
    int j = 0;
    for(int i = 0; i < N; i++) {
        p = this->_config_info.particles[i];
        if (p->type == _restrict_to_type) {
            particles[j] = new BaseParticle<number> ();
            particles[j]->pos = p->pos;
            particles[j]->index = j;
            particles[j]->type = _restrict_to_type;
            j ++;
        }
    }

    // we build the ``test'' particle and add it to the array
    particles[N_of_type] = new BaseParticle<number> ();
    particles[N_of_type]->index = N_of_type;
    particles[N_of_type]->pos.x = drand48() * box_sides[0];
    particles[N_of_type]->pos.y = drand48() * box_sides[1];
    particles[N_of_type]->pos.z = drand48() * box_sides[2];
    particles[N_of_type]->type = _restrict_to_type;

    p = particles[N_of_type];
    int real_N_of_type = N_of_type + 1;
    Cells<number> * cells = new Cells<number>(real_N_of_type, this->_config_info.box);
    cells->init(particles, 2.01 * _sigma);

    int ntries = 0;
    int nfull = 0;
    while (ntries < _ntries) {
        number min = 4.f * _sigma * _sigma + 1.;
        p->pos.x = box_sides[0] * drand48();
        p->pos.y = box_sides[1] * drand48();
        p->pos.z = box_sides[2] * drand48();
        cells->single_update(p);
        std::vector<BaseParticle<number> *> neighs = cells->get_complete_neigh_list(p);
        typename std::vector<BaseParticle<number> *>::iterator it;
        bool found = false;
        for (it = neighs.begin(); it != neighs.end() && found == false; it ++) {
            number mydist2 = this->_config_info.box->sqr_min_image_distance(p, *it);
            if ((*it)->type != _restrict_to_type) throw oxDNAException("fix things please 2");
            if (mydist2 < min) min = mydist2;
            if (min < 4.f * _sigma * _sigma) {
                nfull ++;
                found = true;
            }
        }
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
