/*
 * CoaxVariables.cpp
 *
 *  Created on: Mar 11, 2014
 *      Author: Ben Snodin
 */

#include "CoaxVariables.h"
#include "../Particles/DNANucleotide.h"

#include <sstream>

template<typename number>
CoaxVariables<number>::CoaxVariables() {
	_dna_interaction = new DNAInteraction<number>();
}

template<typename number>
CoaxVariables<number>::~CoaxVariables() {
}

template<typename number>
void CoaxVariables<number>::get_settings(input_file &my_inp, input_file &sim_inp) {
	int tmp = 0;

	tmp = 0;
	getInputInt(&my_inp,"particle1_id", &tmp, 1);
	_particle1_id = tmp;

	tmp = 0;
	getInputInt(&my_inp,"particle2_id", &tmp, 1);
	_particle2_id = tmp;

	std::string inter_type;
	if (getInputString(&sim_inp, "interaction_type", inter_type, 0) == KEY_FOUND){
		if (inter_type.compare("DNA2") == 0) {
			_use_oxDNA2_coaxial_stacking = true;
		}
	}
	// initialise the DNAInteraction object
	_dna_interaction->get_settings(sim_inp);
	_dna_interaction->init();
}

template<typename number>
std::string CoaxVariables<number>::get_output_string(llint curr_step) {
	std::stringstream output_str;

	output_str << "#id1 id2 delta_r_coax within_cutoff theta1 theta4 theta5 theta6 cos(phi3) total_coax, t = " << curr_step << "\n";
	
	BaseParticle<number> *p;
	BaseParticle<number> *q;
	p = this->_config_info.particles[_particle1_id];
	q = this->_config_info.particles[_particle2_id];
	output_str << p->index << " " << q->index << " ";

	// Modified version of coaxial stacking calculation code
	LR_vector<number> r = q->pos.minimum_image(p->pos, *(this->_config_info.box_side));
	LR_vector<number> rstack = r + q->int_centers[DNANucleotide<number>::STACK] - p->int_centers[DNANucleotide<number>::STACK];
	number rstackmod = rstack.module();

	LR_vector<number> rstackdir = rstack / rstackmod;

	// particle axes according to Allen's paper
	LR_vector<number> &a1 = p->orientationT.v1;
	//LR_vector<number> &a2 = p->orientationT.v2;
	LR_vector<number> &a3 = p->orientationT.v3;
	LR_vector<number> &b1 = q->orientationT.v1;
	LR_vector<number> &b3 = q->orientationT.v3;

	// angles involved in the CXST interaction
	number cost1 = -a1 * b1;
	number cost4 =  a3 * b3;
	number cost5 =  a3 * rstackdir;
	number cost6 = -b3 * rstackdir;

	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the coaxial stacking interaction).
	LR_vector<number> rbackboneref = r + b1 * POS_BACK - a1 * POS_BACK;

	number rbackrefmod = rbackboneref.module();
	LR_vector<number> rbackbonerefdir = rbackboneref / rbackrefmod;
	number cosphi3 = rstackdir * (rbackbonerefdir.cross(a1));

	number energy = (number) 0.f;
	bool rstack_within_cutoff = false;
	if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		rstack_within_cutoff = true;
		// functions called at their relevant arguments
		number f2   = _dna_interaction->_f2(rstackmod, CXST_F2);
		number f4t1 = _dna_interaction->_custom_f4 (cost1, CXST_F4_THETA1);
		number f4t4 = _dna_interaction->_custom_f4 (cost4, CXST_F4_THETA4);
		number f4t5 = _dna_interaction->_custom_f4 (cost5, CXST_F4_THETA5) + _dna_interaction->_custom_f4 (-cost5, CXST_F4_THETA5);
		number f4t6 = _dna_interaction->_custom_f4 (cost6, CXST_F4_THETA6) + _dna_interaction->_custom_f4 (-cost6, CXST_F4_THETA6);

		if (_use_oxDNA2_coaxial_stacking){
			energy = f2 * f4t1 * f4t4 * f4t5 * f4t6;
		}
		else{
			number f5cosphi3 = _dna_interaction->_f5(cosphi3, CXST_F5_PHI3);
			energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3;
		}
	}
	
	number t1 = acos(cost1);
	number t4 = acos(cost4);
	number t5 = acos(cost5);
	number t6 = acos(cost6);

	// \delta r_{coax} is equal to dr_{stack}, which is called rstackmod here
	output_str << rstackmod << " ";
	// record whether we're within the cutoff radius for dr_{stack}
	if (rstack_within_cutoff) output_str << "True  ";
	else output_str << "False ";

	output_str << t1 << " " << t4 << " " << t5 << " " << t6 << " " << cosphi3 << " " << energy << "\n";

	return output_str.str();
}

template class CoaxVariables<float>;
template class CoaxVariables<double>;
