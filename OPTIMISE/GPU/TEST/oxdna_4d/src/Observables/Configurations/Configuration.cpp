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

Configuration::Configuration() {
	_back_in_box = false;
	_reduced = false;
	_only_type = -1;
	_print_momenta = true;
}

Configuration::~Configuration() {

}

void Configuration::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);

	getInputBool(&sim_inp, "back_in_box", &_back_in_box, 0);
	getInputBool(&my_inp, "back_in_box", &_back_in_box, 0);
	getInputBool(&my_inp, "reduced", &_reduced, 0);
	getInputInt(&my_inp, "only_type", &_only_type, 0);
	getInputBool(&my_inp, "print_momenta", &_print_momenta, 0);

	string opt;
	bool show_on = false;
	if(getInputString(&my_inp, "show", opt, 0) == KEY_FOUND) {
		show_on = true;
		vector<string> spl = Utils::split(opt, ',');
		for(vector<string>::iterator it = spl.begin(); it != spl.end(); it++)
			_visible_particles.insert(atoi(it->c_str()));
	}

	if(getInputString(&my_inp, "hide", opt, 0) == KEY_FOUND) {
		if(show_on) OX_LOG(Logger::LOG_WARNING, "Observable '%s': Specifying both 'show' and 'hide' options does not make sense. I will not consider the 'hide' one.", typeid(*this).name());
		else {
			vector<string> spl = Utils::split(opt, ',');
			for(vector<string>::iterator it = spl.begin(); it != spl.end(); it++) _hidden_particles.insert(atoi(it->c_str()));
		}
	}
}

void Configuration::init() {
	BaseObservable::init();

	// if no 'show' options is found then the default behaviour is to show all the particles
	if(_visible_particles.size() == 0) {
		for(int i = 0; i < _config_info->N(); i++) {
			_visible_particles.insert(i);
		}
	}
	// otherwise we have to check that all the particles the user wants to show exist
	else {
		for(set<int>::iterator it = _visible_particles.begin(); it != _visible_particles.end(); it++) {
			if(*it < 0 || *it >= _config_info->N()) {
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

	_tot_energy.init();
}

string Configuration::_headers(llint step) {
	stringstream headers;
	headers.precision(15);

	if(_reduced) {
		set<int> strands;
		for(auto p_idx : _visible_particles) {
			BaseParticle *p = _config_info->particles()[p_idx];
			strands.insert(p->strand_id);
		}

		headers << strands.size() << " " << _config_info->box->box_sides().x << " " << _config_info->box->box_sides().y << " " << _config_info->box->box_sides().z << endl;
	}
	else {
		number U = _tot_energy.get_U(step);
		number K = _tot_energy.get_K(step);

		headers << "t = " << step << endl;
		headers << "b = " << _config_info->box->box_sides().x << " " << _config_info->box->box_sides().y << " " << _config_info->box->box_sides().z << endl;
		headers << "E = " << U + K << " " << U << " " << K << endl;
	}

	return headers.str();
}

string Configuration::_particle(BaseParticle *p) {
	stringstream conf;
	conf.precision(15);

	LR_vector box_sides = _config_info->box->box_sides();

	LR_vector mypos;
	if(_back_in_box) {
		mypos.x = p->pos.x - floor(_strands_cdm[p->strand_id].x / box_sides.x) * box_sides.x;
		mypos.y = p->pos.y - floor(_strands_cdm[p->strand_id].y / box_sides.y) * box_sides.y;
		mypos.z = p->pos.z - floor(_strands_cdm[p->strand_id].z / box_sides.z) * box_sides.z;
	}
	else {
		LR_vector number_pos = _config_info->box->get_abs_pos(p);
		mypos.x = number_pos.x;
		mypos.y = number_pos.y;
		mypos.z = number_pos.z;
	}

	LR_matrix oT = p->orientation.get_transpose();
	conf << mypos.x << " " << mypos.y << " " << mypos.z << " ";
	conf << oT.v1.x << " " << oT.v1.y << " " << oT.v1.z << " ";
	conf << oT.v3.x << " " << oT.v3.y << " " << oT.v3.z;

	if(_print_momenta) {
		conf << " " << p->vel.x << " " << p->vel.y << " " << p->vel.z << " ";
		conf << p->L.x << " " << p->L.y << " " << p->L.z;
	}

	return conf.str();
}

string Configuration::_configuration(llint step) {
	stringstream conf;
	conf.precision(15);

	if(_reduced) {
		// prints only the centres of mass of the strands
		map<int, LR_vector> coms;
		map<int, int> counters;

		for(auto p_idx : _visible_particles) {
			BaseParticle *p = _config_info->particles()[p_idx];
			coms[p->strand_id] += _config_info->box->get_abs_pos(p);
			counters[p->strand_id]++;
		}

		for(typename map<int, LR_vector>::iterator it = coms.begin(); it != coms.end(); it++) {
			LR_vector com = it->second / (number) counters[it->first];

			if(it != coms.begin()) conf << endl;
			conf << com.x << " " << com.y << " " << com.z;
		}
	}
	else {
		if(_back_in_box) _fill_strands_cdm();
		// this is used to avoid printing empty lines
		bool empty = true;
		for(auto p_idx : _visible_particles) {
			BaseParticle *p = _config_info->particles()[p_idx];
			bool visible = (_only_type == -1 || p->type == _only_type);
			if(visible) {
				if(p_idx != *_visible_particles.begin() && !empty) conf << endl;
				string p_str = _particle(p);
				conf << p_str;
				empty = (p_str.size() == 0);
			}
		}
	}

	return conf.str();
}

string Configuration::get_output_string(llint curr_step) {
	return _headers(curr_step) + _configuration(curr_step);
}

void Configuration::_fill_strands_cdm() {
	std::map<int, int> nin;
	_strands_cdm.clear();
	for(auto p_idx : _visible_particles) {
		BaseParticle *p = _config_info->particles()[p_idx];
		if(_strands_cdm.count(p->strand_id) == 0) {
			nin[p->strand_id] = 1;
			_strands_cdm[p->strand_id] = p->pos;
		}
		else {
			nin[p->strand_id]++;
			_strands_cdm[p->strand_id] += p->pos;
		}
	}
	for(map<int, int>::iterator it = nin.begin(); it != nin.end(); it++) {
		int strand_id = it->first;
		int nin_id = it->second;
		_strands_cdm[strand_id] /= (number) nin_id;
	}
}

