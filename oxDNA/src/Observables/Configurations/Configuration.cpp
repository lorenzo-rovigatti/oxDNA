/*
 * Configuration.cpp
 *
 *  Created on: 08/ott/2013
 *      Author: lorenzo
 */

#include <sstream>
#include <typeinfo>
#include <map>

#include "Configuration.h"
#include "../../Particles/BaseParticle.h"

using namespace std;

template<typename number>
Configuration<number>::Configuration() {
	_back_in_box = false;
	_reduced = false;
}

template<typename number>
Configuration<number>::~Configuration() {

}

template<typename number>
void Configuration<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int bb;
	if(getInputBoolAsInt(&sim_inp, "back_in_box", &bb, 0) == KEY_FOUND) _back_in_box = (bool) bb;
	if(getInputBoolAsInt(&my_inp, "reduced", &bb, 0) == KEY_FOUND) _reduced = (bool) bb;

	string opt;
	bool show_on = false;
	if(getInputString(&my_inp, "show", opt, 0) == KEY_FOUND) {
		show_on = true;
		vector<string> spl = Utils::split(opt, ',');
		for(vector<string>::iterator it = spl.begin(); it != spl.end(); it++) _visible_particles.insert(atoi(it->c_str()));
	}

	if(getInputString(&my_inp, "hide", opt, 0) == KEY_FOUND) {
		if(show_on) OX_LOG(Logger::LOG_WARNING, "Observable '%s': Specifying both 'show' and 'hide' options does not make sense. I will not consider the 'hide' one.", typeid(*this).name());
		else {
			vector<string> spl = Utils::split(opt, ',');
			for(vector<string>::iterator it = spl.begin(); it != spl.end(); it++) _hidden_particles.insert(atoi(it->c_str()));
		}
	}
}

template<typename number>
void Configuration<number>::init(ConfigInfo<number> &config_info) {
   BaseObservable<number>::init(config_info);

   // if no 'show' options is found then the default behaviour is to show all the particles
   if(_visible_particles.size() == 0) for(int i = 0; i < *config_info.N; i++) _visible_particles.insert(i);
   // otherwise we have to check that all the particles the user wants to show exist
   else {
	   for(set<int>::iterator it = _visible_particles.begin(); it != _visible_particles.end(); it++) {
		   if(*it < 0 || *it >= *config_info.N) {
			   OX_LOG(Logger::LOG_WARNING, "Observable '%s': Particle '%d' cannot be shown because it does not exist", typeid(*this).name(), *it);
			   _visible_particles.erase(*it);
		   }
	   }
   }
   // if a 'hide' option is found then we remove these elements from the _visible_particles set
   if(_hidden_particles.size() > 0) {
	   for(set<int>::iterator it = _hidden_particles.begin(); it != _hidden_particles.end(); it++) {
		   if(_visible_particles.find(*it) == _visible_particles.end()) OX_LOG(Logger::LOG_WARNING, "Observable '%s': Particle '%d' cannot be hidden because it does not exist", typeid(*this).name(), *it);
		   else _visible_particles.erase(*it);
	   }
   }

   _tot_energy.init(config_info);
}

template<typename number>
string Configuration<number>::_headers(llint step) {
	stringstream headers;
	headers.precision(15);

	if(_reduced) {
		set<int> strands;
		for(set<int>::iterator it = _visible_particles.begin(); it != _visible_particles.end(); it++) {
			BaseParticle<number> *p = this->_config_info.particles[*it];
			strands.insert(p->strand_id);
		}

		headers << strands.size() << " " << this->_config_info.box->box_sides().x << " " << this->_config_info.box->box_sides().y << " " << this->_config_info.box->box_sides().z << endl;
	}
	else {
		number U = _tot_energy.get_U(step);
		number K = _tot_energy.get_K(step);
		
		headers << "t = " << step << endl;
		headers << "b = " << this->_config_info.box->box_sides().x << " " << this->_config_info.box->box_sides().y << " " << this->_config_info.box->box_sides().z << endl;
		headers << "E = " << U+K << " " << U << " " << K << endl;
	}

	return headers.str();
}

template<typename number>
string Configuration<number>::_particle(BaseParticle<number> *p) {
	stringstream conf;
	conf.precision(15);

	LR_vector<double> mypos;
	if (_back_in_box){
		mypos.x = p->pos.x - floor(p->pos.x / *this->_config_info.box_side)*(*this->_config_info.box_side);
		mypos.y = p->pos.y - floor(p->pos.y / *this->_config_info.box_side)*(*this->_config_info.box_side);
		mypos.z = p->pos.z - floor(p->pos.z / *this->_config_info.box_side)*(*this->_config_info.box_side);
	}
	else{
		mypos.x = p->get_abs_pos(this->_config_info.box->box_sides()).x;
		mypos.y = p->get_abs_pos(this->_config_info.box->box_sides()).y;
		mypos.z = p->get_abs_pos(this->_config_info.box->box_sides()).z;
	}
	LR_matrix<number> oT = p->orientation.get_transpose();
	conf << mypos.x << " " << mypos.y << " " << mypos.z << " ";
	conf << oT.v1.x << " " << oT.v1.y << " " << oT.v1.z << " ";
	conf << oT.v3.x << " " << oT.v3.y << " " << oT.v3.z << " ";
	conf << p->vel.x << " " << p->vel.y << " " << p->vel.z << " ";
	conf << p->L.x << " " << p->L.y << " " << p->L.z;

	return conf.str();
}

template<typename number>
string Configuration<number>::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	if(_reduced) {
		// prints only the centres of mass of the strands
		map<int, LR_vector<number> > coms;
		map<int, int> counters;

		for(set<int>::iterator it = _visible_particles.begin(); it != _visible_particles.end(); it++) {
			BaseParticle<number> *p = this->_config_info.particles[*it];

			coms[p->strand_id] += p->get_abs_pos(*this->_config_info.box_side);
			counters[p->strand_id]++;
		}

		for(typename map<int, LR_vector<number> >::iterator it = coms.begin(); it != coms.end(); it++) {
			LR_vector<number> com = it->second / (number) counters[it->first];

			if(it != coms.begin()) conf << endl;
			conf << com.x << " " << com.y << " " << com.z;
		}
	}
	else {
		// this is used to avoid printing empty lines
		bool empty = true;
		for(set<int>::iterator it = _visible_particles.begin(); it != _visible_particles.end(); it++) {
			if(it != _visible_particles.begin() && !empty) conf << endl;
			BaseParticle<number> *p = this->_config_info.particles[*it];
			string p_str = _particle(p);
			conf << p_str;
			empty = (p_str.size() == 0);
		}
	}

	return conf.str();
}

template<typename number>
string Configuration<number>::get_output_string(llint curr_step) {
	return _headers(curr_step) + _configuration(curr_step);
}

template<typename number>
void Configuration<number>::_fill_strands_cdm() {
	std::map<int, int> nin;
	_strands_cdm.clear();
	for(set<int>::iterator it = _visible_particles.begin(); it != _visible_particles.end(); it++) {
		BaseParticle<number> *p = this->_config_info.particles[*it];
		if (_strands_cdm.count(p->strand_id) == 0) {
			nin[p->strand_id] = 1;
			_strands_cdm[p->strand_id] = p->pos;
		}
		else {
			nin[p->strand_id] ++;
			_strands_cdm[p->strand_id] += p->pos;
		}
	}
	for (map<int, int>::iterator it = nin.begin(); it != nin.end(); it ++) {
		int strand_id = it->first;
		int nin_id = it->second;
		_strands_cdm[strand_id] /= (number) nin_id;
	}
}

template class Configuration<float>;
template class Configuration<double>;

