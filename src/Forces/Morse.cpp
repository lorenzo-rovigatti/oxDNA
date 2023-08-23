/*
 * Morse.cpp
 *
 *  Created on: 18/oct/2011
 *      Author: Flavio 
 */

#include "Morse.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

Morse::Morse() :
				BaseForce() {
	_ref_id = -2;
	_particle = -2;
	_p_ptr = NULL;
	_r0 = -1.;
	PBC = false;
	_a = 0.f;
	_D = 0.f;
}

std::tuple<std::vector<int>, std::string> Morse::init(input_file &inp) {
	getInputInt(&inp, "particle", &_particle, 1);
	getInputInt(&inp, "ref_particle", &_ref_id, 1);
	getInputNumber(&inp, "r0", &_r0, 1);
	getInputNumber(&inp, "a", &_a, 1);
	getInputNumber(&inp, "D", &_D, 1);
	getInputBool(&inp, "PBC", &PBC, 0);

	int N = CONFIG_INFO->particles().size();
	if(_ref_id < 0 || _ref_id >= N) {
		throw oxDNAException("Invalid reference particle %d for Mutual Trap", _ref_id);
	}
	_p_ptr = CONFIG_INFO->particles()[_ref_id];

	if(_particle >= N || N < -1) {
		throw oxDNAException("Trying to add a Morse on non-existent particle %d. Aborting", _particle);
	}
	if(_particle == -1) {
		throw oxDNAException("Cannot apply Morse to all particles. Aborting");
	}

	std::string description = Utils::sformat("Morse (D=%g, a=%g, r0=%g, ref_particle=%d, PBC=%d)", _D, _a, _r0, _ref_id, PBC);

    // For Debugging
    // graph_data();

	return std::make_tuple(std::vector<int>{_particle}, description);
}

LR_vector Morse::_distance(LR_vector u, LR_vector v) {
    if(PBC) {
        return CONFIG_INFO->box->min_image(u, v);
    }
    else {
        return v - u;
    }
}

LR_vector Morse::value(llint step, LR_vector &pos) {   // Negative of Force
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
	return (dr / dr.module())*2*_a*_D*exp(-_a*dr.module())*(1- exp(-_a*dr.module()));
}

number Morse::potential(llint step, LR_vector &pos) {
	LR_vector dr = _distance(pos, CONFIG_INFO->box->get_abs_pos(_p_ptr));
	return _D*pow(1-exp(-_a*dr.module()), 2) - _D;
}

void Morse::graph_data(){
    double x = -0.5;
    double totaldist = 2.0;
    int npoints = 1000;
    double dt = totaldist/npoints;

    double xarray [npoints];
    double forcearray [npoints];
    double energyarray [npoints];

    for(int i = 0; i < npoints; i++){
        x += dt;
        xarray[i] = x;
        forcearray[i] = 2*_a*_D*exp(-_a*x)*(1- exp(-_a*x));
        energyarray[i] = _D*pow(1-exp(-_a*x), 2);
    }

    printf("Morse Potential X Values: [");
    for(int i = 0; i < npoints; i++){
        printf("%.4f,", xarray[i]);
    }
    printf("]\n");

    printf("Morse Force Values: [");
    for(int i = 0; i < npoints; i++){
        printf("%.4f,", forcearray[i]);
    }
    printf("]\n");

    printf("Morse Potential Values: [");
    for(int i = 0; i < npoints; i++){
        printf("%.4f,", energyarray[i]);
    }
    printf("]\n");
}
