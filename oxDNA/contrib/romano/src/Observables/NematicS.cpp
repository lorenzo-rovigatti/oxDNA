/*
 * NematicS.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: flavio
 */

#include "NematicS.h"

#include "../Utilities/Utils.h"
#include "../Interactions/TEPInteraction.h"

template<typename number> NematicS<number>::NematicS() {
 	_external_dir = false;
	_dir = LR_vector<number> (0.f, 0.f, 0.f);
	_axis_i = -1;
}

template<typename number>
NematicS<number>::~NematicS() {

}

template<typename number>
void NematicS<number>::init(ConfigInfo<number> &config_info) {
    BaseObservable<number>::init(config_info);
}

template<typename number>
void NematicS<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	std::string tmps;
	if (getInputString(&my_inp,"dir", tmps, 0) == KEY_FOUND) {
		_external_dir = true;
		std::vector<std::string> nos = Utils::split(tmps, ',');
		if (nos.size() != 3) throw oxDNAException("Could not parse dir in input_file for observable nematic_s");
		_dir = LR_vector<number> (atof(nos[0].c_str()), atof(nos[1].c_str()), atof(nos[2].c_str()));
		OX_LOG(Logger::LOG_INFO, "NematicS: Using external_dir = %g, %g, %g", _dir.x, _dir.y, _dir.z);
	}
	else {
		throw oxDNAException ("NematicS Not implemented yet without external_dir");
	}

	int tmpi;
	getInputInt(&my_inp, "axis", &tmpi, 1);

	if (tmpi < 0 || tmpi >=6) throw oxDNAException ("Wrond index in observable nematic_s");
	_axis_i = tmpi;
	OX_LOG (Logger::LOG_INFO, "Using axis %d", _axis_i);
}

template<typename number>
LR_vector<number> * NematicS<number>::_get_axis(BaseParticle<number> * p) {
	switch (_axis_i) {
		case NematicS::O_V1: return &(p->orientation.v1);
		case NematicS::O_V2: return &(p->orientation.v2);
		case NematicS::O_V3: return &(p->orientation.v3);
		case NematicS::OT_V1: return &(p->orientationT.v1);
		case NematicS::OT_V2: return &(p->orientationT.v2);
		case NematicS::OT_V3: return &(p->orientationT.v3);
		default: throw oxDNAException ("Should never get here....");
	}
	return (LR_vector<number> *) NULL;
}

template<typename number>
std::string NematicS<number>::get_output_string(llint curr_step) {
	BaseParticle<number> ** particles = this->_config_info.particles;
	int N = *(this->_config_info.N);

	number sum = 0.f;

	for (int i = 0; i < N; i ++) {
		BaseParticle<number> * p = particles[i];

		LR_vector<number> * v = _get_axis(p);

		number dot = (*v) * _dir;
		if (fabs(dot) >= 1.f) dot = copysign(1.f - FLT_EPSILON, dot);
		
		sum += 1.5f * (dot * dot) - 0.5f;
	}

	return Utils::sformat("% 10.6lf", (double) (sum / N));
}

template class NematicS<float>;
template class NematicS<double>;
