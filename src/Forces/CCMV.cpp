/*
 * CCMV.cpp
 *
 *  Created on: 31/may/2024
 *      Author: Lorenzo
 */

#include "CCMV.h"
#include "../Utilities/oxDNAException.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"

CCMV::CCMV() :
				BaseForce() {

}

std::tuple<std::vector<int>, std::string> CCMV::init(input_file &inp) {
	BaseForce::init(inp);

	getInputNumber(&inp, "radius", &_radius, 1);

	getInputNumber(&inp, "A", &_V_amplt, 1);

	getInputNumber(&inp, "cutoff", &_cutoff, 1);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	std::string strdir;
	if(getInputString(&inp, "center", strdir, 0) == KEY_FOUND) {
		auto center_as_vector = Utils::split_to_numbers(strdir, ",");
		if(center_as_vector.size() == 3) {
			_center = LR_vector(center_as_vector[0], center_as_vector[1], center_as_vector[2]);
		}
		else {
			throw oxDNAException("Could not parse center %s in external forces file. Aborting", strdir.c_str());
		}
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "CCMV");
	std::string description = Utils::sformat("CCMV force (radius=%g, center=%g,%g,%g, A,B,C,D,E,F,G,H,I=%g,%g,%g,%g,%g,%g,%g,%g,%g, Ampl=%g, cutoff=%g)", _radius, _center.x, _center.y, _center.z, _A_Vittorio, _B_Vittorio, _C_Vittorio, _D_Vittorio, _E_Vittorio, _F_Vittorio, _G_Vittorio, _H_Vittorio, _I_Vittorio, _V_amplt, _cutoff);

	return std::make_tuple(particle_ids, description);
}

LR_vector CCMV::value(llint step, LR_vector &pos) {
	LR_vector dist = CONFIG_INFO->box->min_image(_center, pos);
	number dist_surface = _radius - dist.module();
	if(dist_surface < 0) {
		throw oxDNAException("Found a particle beyond the Yukawa sphere radius %lf. Check your initial configuration or your timestep.", _radius);
	}

	LR_vector force;
	if(dist_surface < _cutoff) {
		
		LR_vector direction = dist / dist.module(); 

		force = direction * _V_amplt * (exp(_A_Vittorio*dist_surface) * (
			_B_Vittorio +   
			dist_surface*(2*_C_Vittorio + _A_Vittorio*_B_Vittorio) +
			pow(dist_surface,2)*(3*_D_Vittorio + _A_Vittorio*_C_Vittorio) +
			pow(dist_surface,3)*(4*_E_Vittorio + _A_Vittorio*_D_Vittorio) +
			pow(dist_surface,4)*(5*_F_Vittorio + _A_Vittorio*_E_Vittorio) +
			pow(dist_surface,5)*(6*_G_Vittorio + _A_Vittorio*_F_Vittorio) +
			pow(dist_surface,6)*(7*_H_Vittorio + _A_Vittorio*_G_Vittorio) +
			pow(dist_surface,7)*(8*_I_Vittorio + _A_Vittorio*_H_Vittorio) +
			pow(dist_surface,8)*(9*_L_Vittorio + _A_Vittorio*_I_Vittorio) +
			pow(dist_surface,9)*(10*_M_Vittorio + _A_Vittorio*_L_Vittorio)+
			pow(dist_surface,10)*(_A_Vittorio*_M_Vittorio))
		);
	}

	return force;
}

number CCMV::potential(llint step, LR_vector &pos) {
	number sqr_r = CONFIG_INFO->box->sqr_min_image_distance(_center, pos);
	number dist_surface = _radius - sqrt(sqr_r);
	if(dist_surface < 0) {
		return 1e10;
	}

	number energy = 0.;
	if(dist_surface < _cutoff) {

		energy = _V_amplt * exp(_A_Vittorio * dist_surface) * (
			    _B_Vittorio * dist_surface +
			    _C_Vittorio * pow(dist_surface,2) +
			    _D_Vittorio * pow(dist_surface, 3) +
			    _E_Vittorio * pow(dist_surface, 4) +
				_F_Vittorio * pow(dist_surface, 5) +
				_G_Vittorio * pow(dist_surface, 6) +
				_H_Vittorio * pow(dist_surface, 7) +
				_I_Vittorio * pow(dist_surface, 8) +
				_L_Vittorio * pow(dist_surface, 9) +
				_M_Vittorio * pow(dist_surface, 10)
		);
	}
	return energy;
}
