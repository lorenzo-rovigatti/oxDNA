/*
 * TclOutput.cpp
 *
 *  Created on: 25/ott/2013
 *      Author: Flavio 
 */

#include <sstream>
#include "TclOutput.h"

template<typename number>
TclOutput<number>::TclOutput() : Configuration<number>() {
	_print_labels = false;
	_back_radius = 0.35;
	_base_radius = 0.25;
	_backbase_radius = 0.15;
	_backback_radius = 0.25;
	_resolution = 20;
	_ref_particle_id = -1;
	_ref_strand_id = -1;
	_back_in_box = true;
}

template<typename number>
TclOutput<number>::~TclOutput() {

}

template<typename number>
void TclOutput<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
	int tmp;
	if(getInputBoolAsInt(&my_inp, "back_in_box", &tmp, 0) == KEY_FOUND) _back_in_box = (bool)tmp;
	if(getInputBoolAsInt(&my_inp, "print_labels", &tmp, 0) == KEY_FOUND) _print_labels = (bool)tmp;
	if(getInputInt(&my_inp, "ref_particle", &tmp, 0) == KEY_FOUND) _ref_particle_id = tmp;
	if(getInputInt(&my_inp, "ref_strand", &tmp, 0) == KEY_FOUND) _ref_strand_id = tmp;
	if(getInputInt(&my_inp, "resolution", &tmp, 0) == KEY_FOUND) _resolution = tmp;
}

template<typename number>
std::string TclOutput<number>::_headers(llint step) {
	std::stringstream headers;
	
	headers << "color Display Background white" << endl;
	headers << "mol new" << endl;

	// we might want to be able to change these in the future...
	double _box_radius = 0.1;
	int _box_resolution = _resolution;
	
	number mybox = *this->_config_info.box_side;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" << -mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} {" << -mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} {" <<  mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} {" << -mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " << -mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " << -mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" <<  mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" <<  mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;
	headers << "graphics 0 cylinder {" << -mybox / 2. << " " <<  mybox / 2. << " " << -mybox / 2. << "} {" << -mybox / 2. << " " <<  mybox / 2. << " " <<  mybox / 2. << "} radius " << _box_radius << " resolution " << _box_resolution << " filled yes" << endl;

	return headers.str();
}

template<typename number>
std::string TclOutput<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	
	DNANucleotide<number> * me, * next;
	me = reinterpret_cast<DNANucleotide<number> *> (p);
	next = reinterpret_cast<DNANucleotide<number> *> (p->n5);
	
	LR_vector<number> zero (0., 0., 0.);
	if (_ref_particle_id >= 0 && this->_visible_particles.count(_ref_particle_id) == 1) zero = this->_config_info.particles[_ref_particle_id]->pos; 
	if (_ref_strand_id >= 0 && this->_strands_cdm.count(_ref_strand_id) == 1) zero = this->_strands_cdm[_ref_strand_id];
	
	// set the colour according to the strand id
	int colorid = p->strand_id;
	if (colorid >= 8) colorid ++;
	colorid = colorid % 33;
	
	res << "graphics 0 color " << colorid << endl;
	
	number mybox = *this->_config_info.box_side;
	LR_vector<number> my_strand_cdm = this->_strands_cdm[me->strand_id];
	LR_vector<number> origin (0., 0., 0.);
	origin = zero;
	LR_vector<number> back = (me->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + me->int_centers[0]; 
	LR_vector<number> base = (me->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + me->int_centers[2]; 
	
	if (_print_labels && p->n5 != P_VIRTUAL && p->n3 == P_VIRTUAL) res << "graphics 0 text {" << back.x << " " << back.y << " " << back.z << "} \"" << me->strand_id << "\" size 1.5 " << endl;
	
	// backbone
	res << "graphics 0 sphere {" << back.x << " " << back.y << " " << back.z << "} radius " << _back_radius << " resolution " << _resolution << endl;
	// base
	res << "graphics 0 sphere {" << base.x << " " << base.y << " " << base.z << "} radius " << _base_radius << " resolution " << _resolution << endl;
	// backbone-base cilinder
	res << "graphics 0 cylinder {" << base.x << " " << base.y << " " << base.z << "} {" << back.x << " " << back.y << " " << back.z << "} radius " << _backbase_radius << " resolution " << _resolution << " filled yes" << endl; 
	
	if (next != P_VIRTUAL) {
		LR_vector<number> nextback = (next->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + next->int_centers[0]; 
		res << "graphics 0 cylinder {" << nextback.x << " " << nextback.y << " " << nextback.z << "} {" << back.x << " " << back.y << " " << back.z << "} radius " << _backback_radius << " resolution " << _resolution << " filled yes" << endl; 
	}
	else {
		// we draw the cone
		if (me->n3 != P_VIRTUAL) {
			next = reinterpret_cast<DNANucleotide<number>  *> (me->n3);
			LR_vector<number> nextback = (next->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + next->int_centers[0]; 
			LR_vector<number> dir = back - nextback;
			LR_vector<number> start = back;
			LR_vector<number> end = back + dir * (2.75 / 3.);
			dir.normalize();
			res << "graphics 0 cone {" << start.x << " " << start.y << " " << start.z << "} {" << end.x << " " << end.y << " " << end.z << "} radius " << _back_radius - 0.01 << " resolution " << _resolution << endl; 
			//ret += "graphics 0 cone {%lf %lf %lf} {%lf %lf %lf} radius 0.35 resolution 20\n" % (start[0], start[1], start[2], end[0], end[1], end[2]);
		}
	}
	
	res << "graphics 0 color white";

	return res.str();
}


template<typename number>
std::string TclOutput<number>::_configuration(llint step) {
	stringstream conf;
	//conf.precision(15);
	if (_back_in_box) this->_fill_strands_cdm ();
	//return Configuration<number>::_configuration(step);
	for(set<int>::iterator it = this->_visible_particles.begin(); it != this->_visible_particles.end(); it++) {
		if(it != this->_visible_particles.begin()) conf << endl;
		BaseParticle<number> *p = this->_config_info.particles[*it];
		conf << _particle(p);
	}
	return conf.str();
}

template class TclOutput<float>;
template class TclOutput<double>;

