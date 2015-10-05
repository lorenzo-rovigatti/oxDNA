/*
 * TEPxyzOutput.h
 *
 *  Created on: 19/jun/2015
 *      Author: Ferdinando (after TcLOutput.cpp, written by Flavio, and TEPtclOutput.cpp, written by Ferdinando)
 */

#include <sstream>
#include "TEPxyzOutput.h"

template<typename number>
TEPxyzOutput<number>::TEPxyzOutput() : Configuration<number>() {
	_print_labels = false;
	
	_weak_bead_1_index = -1;
	_weak_bead_2_index = -1;
	
	_core_radius = 1;
	_side_radius = 0.7;
	_front_radius = 0.82;
	
	_side_shift=0.3;
	_front_shift=0.15;
	
	_resolution = 20;
	_ref_particle_id = -1;
	_ref_strand_id = -1;
	_back_in_box = true;

}

template<typename number>
TEPxyzOutput<number>::~TEPxyzOutput() {

}

template<typename number>
void TEPxyzOutput<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
	int tmp;
	if(getInputBoolAsInt(&my_inp, "back_in_box", &tmp, 0) == KEY_FOUND) _back_in_box = (bool)tmp;
	if(getInputBoolAsInt(&my_inp, "print_labels", &tmp, 0) == KEY_FOUND) _print_labels = (bool)tmp;
	if(getInputInt(&my_inp, "ref_particle", &tmp, 0) == KEY_FOUND) _ref_particle_id = tmp;
	if(getInputInt(&my_inp, "ref_strand", &tmp, 0) == KEY_FOUND) _ref_strand_id = tmp;
	if(getInputInt(&my_inp, "resolution", &tmp, 0) == KEY_FOUND) _resolution = tmp;
	
	getInputInt(&sim_inp,"TEP_weakened_bead_index",&_weak_bead_1_index,0);
	printf("%d weak\n",_weak_bead_1_index);
	getInputInt(&sim_inp,"TEP_weakened_bead_index2",&_weak_bead_2_index,0);
}

template<typename number>
std::string TEPxyzOutput<number>::_headers(llint step) {
	std::stringstream headers;
	
	const number mybox = *this->_config_info.box_side;
	const int N = *this->_config_info.N*3;

	headers << N << endl;
	headers << "#" << mybox <<" "<< mybox <<" "<< mybox << endl;
/* // Previous headers TODO remove them when the observable is finished
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

*/
	return headers.str();
}

template<typename number>
std::string TEPxyzOutput<number>::_particle(BaseParticle<number> *p) {
	std::stringstream res;
	
	TEPParticle<number> * me;
	me = reinterpret_cast<TEPParticle<number> *> (p);
	//next = reinterpret_cast<TEPParticle<number> *> (p->n5);
	
	LR_vector<number> zero (0., 0., 0.);
	if (_ref_particle_id >= 0 && this->_visible_particles.count(_ref_particle_id) == 1) zero = this->_config_info.particles[_ref_particle_id]->pos; 
	if (_ref_strand_id >= 0 && this->_strands_cdm.count(_ref_strand_id) == 1) zero = this->_strands_cdm[_ref_strand_id];
		
	// set the colour according to the strand id
	/* //this is used to set the color in tclOutput. I could do something similar with atom type should I ever need it.
	int colorid = p->strand_id;
	if (colorid >= 8) colorid ++;
	colorid = colorid % 33;
	res << "graphics 0 color " << colorid << endl;
	*/
	
	number mybox = *this->_config_info.box_side;
	LR_vector<number> my_strand_cdm = this->_strands_cdm[me->strand_id];
	LR_vector<number> origin (0., 0., 0.);
	origin = zero;

	LR_vector<number> core = (me->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox);
	LR_vector<number> side = (me->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + _side_shift*me->orientationT.v2; 
	LR_vector<number> front = (me->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + _front_shift*me->orientationT.v1; 

/*
//This was used to print labels in tcl. I don't think there's a way of doing it in xyz.	
	if (_print_labels && p->n5 != P_VIRTUAL && p->n3 == P_VIRTUAL) res << "graphics 0 text {" << core.x << " " << core.y << " " << core.z << "} \"" << me->strand_id << "\" size 1.5 " << endl;
*/	

	// core
	if (me->index == _weak_bead_1_index){
		res << "Po" <<" "<< core.x << " " << core.y << " " << core.z << endl;
	}
	else if (me->index == _weak_bead_2_index){
		res << "Fe" <<" "<< core.x << " " << core.y << " " << core.z << endl;
	}
	else{	
		res << "Cr" <<" "<< core.x << " " << core.y << " " << core.z << endl;
	}
	// base
	res << "H" <<" "<< side.x << " " << side.y << " " << side.z << endl;
	// front 
	res << "He" <<" "<< front.x << " " << front.y << " " << front.z ;
	
	
	//res << "graphics 0 color white";

	return res.str();
}


template<typename number>
std::string TEPxyzOutput<number>::_configuration(llint step) {
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

template class TEPxyzOutput<float>;
template class TEPxyzOutput<double>;

