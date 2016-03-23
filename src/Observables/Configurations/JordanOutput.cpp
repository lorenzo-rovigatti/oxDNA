/*
 * JordanOutput.cpp
 *
 *  Created on: 23/mar/2016
 *      Author: Flavio 
 */

#include <sstream>
#include "JordanOutput.h"
#include "../../Particles/JordanParticle.h"

template<typename number>
JordanOutput<number>::JordanOutput() : Configuration<number>() {

}

template<typename number>
JordanOutput<number>::~JordanOutput() {

}

template<typename number>
void JordanOutput<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
}

template<typename number>
std::string JordanOutput<number>::_headers(llint step) {
	std::stringstream headers;
	
	return headers.str();
}

template<typename number>
std::string JordanOutput<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	
	JordanParticle<number> * me = reinterpret_cast<JordanParticle<number> *> (p);

	res << me->get_output_string();
	
	return res.str();
}


template<typename number>
std::string JordanOutput<number>::_configuration(llint step) {
	stringstream conf;
	for(set<int>::iterator it = this->_visible_particles.begin(); it != this->_visible_particles.end(); it ++) {
		if(it != this->_visible_particles.begin()) conf << endl;
		BaseParticle<number> *p = this->_config_info.particles[*it];
		conf << _particle(p);
	}
	return conf.str();
}

template class JordanOutput<float>;
template class JordanOutput<double>;

