/*
 * VectorAngle.cpp
 *
 *  Created on: Dec 9, 2014
 *      Author: Ferdinando Randisi
 */

#include "VectorAngle.h"
#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"

VectorAngle::VectorAngle() {

	_angle_index = 1;
	_first_particle_index = 0;
	_last_particle_index = -1;

	_print_local_details = true;
}

VectorAngle::~VectorAngle() {

}

void VectorAngle::init() {
	BaseObservable::init();

	std::vector<BaseParticle *> &p = _config_info->particles();
	const int N = _config_info->N();

	// check that _first_particle_index is in [0,N) and not the terminal particle
	if(_first_particle_index < 0 || _first_particle_index > N - 1) throw oxDNAException("VectorAngle: first_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.", N, _first_particle_index);
	if(p[_first_particle_index]->n5 == P_VIRTUAL) {
		throw oxDNAException("VectorAngle: first_particle_index must not be the index of the last particle of a strand, otherwise which particle should be last_particle_index referring to?");
	}

	if(_last_particle_index == -1) {
		_last_particle_index = _first_particle_index;
		do { // _last_particle_index will either be the index of the last of the chain or (if the chain is circular) _first_particle_index itself
			_last_particle_index = p[_last_particle_index]->n5->index;

		} while(p[_last_particle_index]->n5 != P_VIRTUAL && _last_particle_index != _first_particle_index);
	}
	// check that _last_particle_index is in [0,N)
	if(_last_particle_index < 0 || _last_particle_index > N - 1) throw oxDNAException("VectorAngle: last_particle_index must be greater or equal to 0 and less than the total number of particles (%d, in this simulation), but it is %d.", N, _last_particle_index);
	// check that _first_particle_index is different than _last_particle_index
	if(_first_particle_index == _last_particle_index) {
		throw oxDNAException("VectorAngle: first_particle_index can't be equal to last_particle_index. If you're talking about a circular molecule then just don't set last_particle_index.");
	}

	// if the molecule is not circular, the last bead is not considered since it just mirrors the behaviour of the first-but-last bead
	if(p[_last_particle_index]->n5 == P_VIRTUAL) {
		_last_particle_index = p[_last_particle_index]->n3->index;
	}
	// check that first and last particle are on the same strand.
	if(p[_first_particle_index]->strand_id != p[_last_particle_index]->strand_id) {
		throw oxDNAException("In observable vector_angle, the first particle (index %d) and the last particle (index %d) are not on the same strand. They're supposed to define a domain of topologically adjacent particles, but obviously they don't.");
	}
	// check that you can actually start from the first particle and get to the last one going forward along the chain.
	int test_index = p[_first_particle_index]->n5->index;
	while(p[test_index]->n5 != P_VIRTUAL && test_index != _first_particle_index && test_index != _last_particle_index) {
		test_index = p[test_index]->n5->index;
	}
	if(p[test_index]->n5 == P_VIRTUAL) {
		throw oxDNAException("In observable vector_angle, could not get from particle %d to particle %d by going forward.");
	}
	if(p[test_index]->n5->index == _first_particle_index) {
		throw oxDNAException("In observable vector_angle, could not get from particle %d to particle %d by going forward. This is very strange since they are both on the same strand as far as I know, so one of the developers (probably Ferdinando) messed something up. Please report the occurrence of this error to the developers.");
	}
}

void VectorAngle::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputInt(&my_inp, "first_particle_index", &_first_particle_index, 0);
	getInputInt(&my_inp, "last_particle_index", &_last_particle_index, 0);
	getInputBool(&my_inp, "print_local_details", &_print_local_details, 0);

	// get the angle to consider 
	if(getInputInt(&my_inp, "angle_index", &_angle_index, 0) == KEY_FOUND) {
		if(_angle_index < 0 || _angle_index > 4) throw oxDNAException("(VectorAngle.cpp) angle_index must be 0,1,2,3, or 4, but was set to %d.", _angle_index);
	}

}

std::string VectorAngle::get_output_string(llint curr_step) {
	int sign, number_of_values = 0;
	std::string result;
	std::vector<BaseParticle *> &p = _config_info->particles();
	double delta_omega = 0;
	number result_number = 0;
	LR_vector u, up, v, vp, f, fp;

	int i = _first_particle_index;
	do {
		if(_angle_index == 0) {
			// omega + gamma = twisting angle
			u = p[i]->orientationT.v1;
			up = p[i]->n5->orientationT.v1;

			f = p[i]->orientationT.v2;
			fp = p[i]->n5->orientationT.v2;

			v = p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;

			if(up * (f.cross(fp)) > 0) sign = 1;
			else sign = -1;
			delta_omega = LRACOS(( f*fp + v*vp )/(1 + u*up)) / (2 * PI);
			if(_print_local_details) {
				result += Utils::sformat("%14.14lf\t", delta_omega * sign) + ' ';
			}
			else {
				result_number += delta_omega * sign;
			}

		}
		else if(_angle_index == 1) {
			// u * u
			u = p[i]->orientationT.v1;
			up = p[i]->n5->orientationT.v1;

			if(_print_local_details) {
				result += Utils::sformat("%14.14lf\t", u * up) + ' ';
			}
			else {
				result_number += u * up;
			}
		}
		else if(_angle_index == 2) {
			// f * f
			f = p[i]->orientationT.v2;
			fp = p[i]->n5->orientationT.v2;

			if(_print_local_details) {
				result += Utils::sformat("%14.14lf\t", f * fp) + ' ';
			}
			else {
				result_number += f * fp;
			}
		}
		else if(_angle_index == 3) {
			// v * v
			v = p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;

			if(_print_local_details) {
				result += Utils::sformat("%14.14lf\t", v * vp) + ' ';
			}
			else {
				result_number += v * vp;
			}
		}
		else if(_angle_index == 4) {
			// Total Twist:
			// 1/(2*pi) * t_i . (v_i x v_i+1)
			v = p[i]->orientationT.v3;
			vp = p[i]->n5->orientationT.v3;
			u = p[i]->orientationT.v1;

			if(_print_local_details) {
				result += Utils::sformat("%14.14lf\t", v * vp) + ' ';
			}
			else {
				result_number += (1. / (2 * M_PI)) * (u * (v.cross(vp)));
			}

		}

		i = p[i]->n5->index;
		number_of_values++;
	} while(i != _last_particle_index);
	if(!_print_local_details) {
		if(_angle_index == 0 || _angle_index == 4) {
			result = Utils::sformat("%14.14lf", result_number);
		}
		else {
			result = Utils::sformat("%14.14lf", result_number / number_of_values);
		}
	}

	return result;
}
