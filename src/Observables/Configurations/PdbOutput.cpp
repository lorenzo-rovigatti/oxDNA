/*
 * PdbOutput.cpp
 *
 *  Created on: 18/jan/2014
 *      Author: C. Matek
 */

#include <sstream>
#include "PdbOutput.h"

const std::string strtypes[] = {"ALA","GLY","CYS","TYR","ARG","PHE","LYS","SER","PRO","VAL","ASN","ASP","CYX","HSP","HSD","MET","LEU"};

template<typename number>
PdbOutput<number>::PdbOutput() : Configuration<number>() {
	_backbase_radius = 0.15;
	_ref_particle_id = -1;
	_ref_strand_id = -1;
	_back_in_box = true;
	_model_nr=0;
}

template<typename number>
PdbOutput<number>::~PdbOutput() {

}

template<typename number>
void PdbOutput<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	Configuration<number>::get_settings(my_inp, sim_inp);
	int tmp;
	if(getInputBoolAsInt(&my_inp, "back_in_box", &tmp, 0) == KEY_FOUND) _back_in_box = (bool)tmp;
	if(getInputInt(&my_inp, "ref_particle", &tmp, 0) == KEY_FOUND) _ref_particle_id = tmp;
	if(getInputInt(&my_inp, "ref_strand", &tmp, 0) == KEY_FOUND) _ref_strand_id = tmp;
}

template<typename number>
std::string PdbOutput<number>::_headers(llint step) {
	std::stringstream headers;
	_model_nr+=1;
	char model_buf[264];
	sprintf(model_buf,"\nMODEL     %4i", (int) _model_nr);
	headers << "HEADER    frame t= "<< step<<model_buf<<" \nREMARK ## 0,0\n";

	return headers.str();
}

template<typename number>
std::string PdbOutput<number>::_particle(BaseParticle<number> *p,std::string strtype) {
	std::stringstream res;
	std::string _strtype = strtype;
	LR_matrix<number> I_b = LR_matrix<number>(LR_vector<number> (1.1,0,0),LR_vector<number> (0,1.5,0),LR_vector<number> (0,0,2.2));
	LR_matrix<number> I_l = p->orientation*I_b*p->orientationT;
	LR_matrix<number> anis = I_l*LR_matrix<number>(-1.,0.,0.,0.,-1.,0.,0.,0.,-1.);

	anis.v1.x = (I_l.v2.y+I_l.v3.z-I_l.v1.x)/2.0;
	anis.v2.y = (I_l.v1.x+I_l.v3.z-I_l.v2.y)/2.0;
	anis.v3.z = (I_l.v1.x+I_l.v2.y-I_l.v3.z)/2.0;

	number mybox = *this->_config_info.box_side;
	LR_vector<number> my_strand_cdm = this->_strands_cdm[p->strand_id];
	LR_vector<number> zero (0., 0., 0.);
	if(_ref_strand_id>=0 && this->_strands_cdm.count(_ref_strand_id) == 1) 
		zero = this->_config_info.particles[_ref_particle_id]->pos;

	if(_ref_particle_id>=0 && this->_visible_particles.count(_ref_particle_id) == 1) 
		zero = this->_config_info.particles[_ref_particle_id]->pos;
	else
		zero = this->_config_info.particles[0]->pos;

	LR_vector<number> origin(0., 0., 0.);
	origin = zero;
	LR_vector<number> back = (p->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + p->int_centers[0];
	LR_vector<number> base = (p->pos - my_strand_cdm) + my_strand_cdm.minimum_image (origin, mybox) + p->int_centers[2];
	base = base-(base-back)*_backbase_radius;
	
	char back_buf1[264];
	char back_buf2[264];
	char base_buf1[264];
	char base_buf2[264];
	char atom_id;
        switch(abs(p->btype%4)){
        case 0:
          atom_id='O';
	  break;
	case 1:
          atom_id='S';
          break;
        case 2:
          atom_id='K';
	  break;
	case 3:
          atom_id='P';
          break;
        default:
          //This should not really happen...                                                                                                          
          atom_id='H';
          break;
	}

	sprintf(back_buf1,"ATOM  %5d %4s %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",2*p->get_index()+1,"A",strtype.c_str(),'A',p->get_index()%10000,' ',back.x,back.y,back.z,1.0,7.895);
	sprintf(back_buf2,"ANISOU%5d %4s %3s %c%4d%c %7i%7i%7i%7i%7i%7i\n",2*p->get_index()+1,"A",strtype.c_str(),'A',p->get_index()%10000,' ',1000,1000,1000,0,0,0);

	sprintf(base_buf1,"ATOM  %5d %4c %3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",2*p->get_index()+2, atom_id, strtype.c_str(), 'C', p->get_index()%10000,' ',base.x,base.y,base.z,1.0,6.316);
	sprintf(base_buf2,"ANISOU%5d %4c %3s %c%4d%c %7i%7i%7i%7i%7i%7i",2*p->get_index()+2, atom_id, strtype.c_str(), 'C', p->get_index()%10000,' ', (int) (anis.v1.x*1000), (int) (anis.v2.y*1000), (int) (anis.v3.z*1000), (int) (anis.v1.y*1000), (int) (anis.v1.z*1000), (int) (anis.v2.z*1000));
	
	res<<back_buf1<<back_buf2<<base_buf1<<base_buf2;

	return res.str();
}


template<typename number>
std::string PdbOutput<number>::_configuration(llint step) {
	stringstream conf;
	//conf.precision(15);
	if (_back_in_box) this->_fill_strands_cdm ();

	for(set<int>::iterator it = this->_visible_particles.begin(); it != this->_visible_particles.end(); it++) {
		if(it != this->_visible_particles.begin()) conf << endl;
		BaseParticle<number> *p = this->_config_info.particles[*it];
		std::string cur_strtype = strtypes[(p->strand_id+1)%17];
		conf << _particle(p,cur_strtype);
	}
	conf << "\nREMARK  ######### \n\nTER \nENDMDL \n ";
	return conf.str();
}

template class PdbOutput<float>;
template class PdbOutput<double>;

