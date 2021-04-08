/*
 * MC_CPUBackend.cpp
 *
 *  Created on: 26/nov/2010
 *      Author: lorenzo
 */


#include "MC_CPUBackend.h"
#include "IOManager.h"

template<typename number>
MC_CPUBackend<number>::MC_CPUBackend(IOManager *IO) : MCBackend<number>(IO), _particles_old(NULL) {
	this->_is_CUDA_sim = false;
	// initialize the messages for the timings output
	this->_timer_msgs_number = 2;
	strncpy(this->_timer_msgs[0], "MC step", 256);
	strncpy(this->_timer_msgs[1], "Lists update", 256);
}

template<typename number>
MC_CPUBackend<number>::~MC_CPUBackend() {
	if (_particles_old != NULL) delete [] _particles_old;
	if(this->_N_updates > 0) divide_given_timing(&this->_timer, 1, this->_N / (double) this->_N_updates);
}

template<typename number>
//void MC_CPUBackend<number>::init(ifstream &conf_input) {
void MC_CPUBackend<number>::init(char conf_filename[256]) {
//	MCBackend<number>::init(conf_input);
	MCBackend<number>::init(conf_filename);

	_particles_old = new Particle<number>[this->_N];
	for(int i = 0; i < this->_N; i++) {
		_particles_old[i].index = i;
		_particles_old[i].type = this->_particles[i].type;
		_particles_old[i].init(this->_max_neigh);

		// this is needed for the first _compute_energy()
		this->_particles[i].set_positions();
		this->_particles[i].orientationT = this->_particles[i].orientation.get_transpose();
	}

	if(this->_delta[MC_MOVE_TRANSLATION] * sqrt(3) > this->_verlet_skin) this->_IO->die("verlet_skin must be > delta_translation times sqrt(3) (the maximum displacement)");

	this->_update_lists();
	_compute_energy();
	if(this->_overlap == true) this->_IO->die("There is an overlap in the initial configuration");
}

template<typename number>
inline number MC_CPUBackend<number>::_particle_particle_bonded_interaction(Particle<number> *p, bool update) {
	number energy = (number) 0;

	if(p->n3 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n3];
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = q->pos - p->pos;

		// FENE
		LR_vector<number> rback = r + q->pos_back - p->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;
		// we check wheter we ended up OUTSIDE of the FENE range
		if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
		    this->_overlap = true;
		    return (number) (1.e6 * this->_T);
		}
		energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

		// excluded volume

		// BASE-BASE
		LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
		energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);

		// P-BASE vs. Q-BACK
		rcenter = r + q->pos_back - p->pos_base;
		energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);

		// P-BACK vs. Q-BASE
		rcenter = r + q->pos_base - p->pos_back;
		energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);

		// STACKING
		LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 =-b3 * rstackdir;
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		// functions and their derivatives needed for energies and forces
		number f1     = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
		number f4t4   = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5   = this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6   = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		if (update) this->_U_stack += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	}

	if(p->n5 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n5];
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = p->pos - q->pos;

		// FENE
		LR_vector<number> rback = r + p->pos_back - q->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;
		// we check wheter we ended up OUTSIDE of the FENE range
		if (fabs(rbackr0) > FENE_DELTA - DBL_EPSILON) {
		    this->_overlap = true;
		    return (number) (1.e6 * this->_T);
		}
		energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

		// excluded volume

		// BASE-BASE
		LR_vector<number> rcenter = r + p->pos_base - q->pos_base;
		energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);

		// P-BASE vs. Q-BACK
		rcenter = r + p->pos_back - q->pos_base;
		energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);

		// P-BACK vs. Q-BASE
		rcenter = r + p->pos_base - q->pos_back;
		energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);

		// STACKING
		LR_vector<number> rstack = r + p->pos_stack - q->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 =-b3 * rstackdir;
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		number f1     = this->_interaction.f1(rstackmod, STCK_F1, p->type, q->type);
		number f4t4   = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5   = this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6   = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		if (update) this->_U_stack += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	}

	return energy;
}

template<typename number>
inline number MC_CPUBackend<number>::_excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc) {
	number rmod = r.module();
	number energy = 0;

	if(rmod < rc) {
		if(rmod > rstar) {
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
		}
		else {
			number lj_part = SQR(sigma / rmod) * SQR(sigma / rmod) * SQR(sigma / rmod);
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}

	return energy;
}

template<typename number>
inline number MC_CPUBackend<number>::_particle_particle_interaction(Particle<number> *p, Particle<number> *q) {
	// true if p and q are Watson-Crick pairs
	//bool is_pair = (q->type + p->type == 3);
	bool is_pair = (q->btype + p->btype == 3);

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	number energy = 0;

	// excluded volume

	// BASE-BASE
	LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
	energy += _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);

	// P-BASE vs. Q-BACK
	rcenter = r + q->pos_back - p->pos_base;
	energy += _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);

	// P-BACK vs. Q-BASE
	rcenter = r + q->pos_base - p->pos_back;
	energy += _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);

	// BACK-BACK
	rcenter = r + q->pos_back - p->pos_back;
	energy += _excluded_volume(rcenter, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);

	// HYDROGEN BONDING
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
	  	LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

	 	 // functions called at their relevant arguments
		number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
		number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[HYDR_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[HYDR_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[HYDR_F4_THETA3]);

		number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[HYDR_F4_THETA4]);
		number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[HYDR_F4_THETA7]);
		number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[HYDR_F4_THETA8]);

		number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += hb_energy;
		this->_U_hydr += hb_energy;
	}
	// END OF HYDROGEN BONDING

	// CROSS STACKING
	LR_vector<number> rcstack = rhydro;
	number rcstackmod = rhydromod;
	//number rcstackmod = rcstack.module();
	if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
	  	LR_vector<number> rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 =  a1 * rcstackdir;
		number cost4 =  a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 =  a3 * rcstackdir;

	 	 // functions called at their relevant arguments
		number f2   = this->_interaction.f2(rcstackmod, CRST_F2);
	  	number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[CRST_F4_THETA1]);
	  	number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[CRST_F4_THETA2]);
	  	number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[CRST_F4_THETA3]);
	  	number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[CRST_F4_THETA4]) + this->_interaction.query_mesh (-cost4, this->_interaction.mesh_f4[CRST_F4_THETA4]);
	  	number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[CRST_F4_THETA7]) + this->_interaction.query_mesh (-cost7, this->_interaction.mesh_f4[CRST_F4_THETA7]);
;
	  	number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[CRST_F4_THETA8]) + this->_interaction.query_mesh (-cost8, this->_interaction.mesh_f4[CRST_F4_THETA8]);
;

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += cstk_energy;
	}

	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
	number rstackmod = rstack.module();
	if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
	  	LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		number cost1 = -a1 * b1;
		number cost4 =  a3 * b3;
		number cost5 =  a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		LR_vector<number> rbackbone = r + q->pos_back - p->pos_back;
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));

	 	// functions called at their relevant arguments
		number f2   = this->_interaction.f2(rstackmod, CXST_F2);
	  	number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[CXST_F4_THETA1]);
	  	number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[CXST_F4_THETA4]);
	  	number f4t5 = this->_interaction.query_mesh (cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]) + this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]);
	  	number f4t6 = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]) + this->_interaction.query_mesh (-cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]);
		number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

		number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy += cxst_energy;
	}

	return energy;
}

template<typename number>
inline number MC_CPUBackend<number>::_particle_energy(Particle<number> *p) {
	this->_overlap = false;
	number res = _particle_particle_bonded_interaction(p);
	if(this->_overlap == true) return res;

	p->prepare_list();
	int neigh = p->next_neighbour();
	while(neigh != P_VIRTUAL) {
		res += _particle_particle_interaction(p, &this->_particles[neigh]);
		neigh = p->next_neighbour();

		if(this->_overlap == true) return res;
	}

	return res;
}

template<typename number>
void MC_CPUBackend<number>::_compute_energy() {
	int neigh;
	Particle<number> *p;

	this->_U = this->_U_hydr = this->_U_stack = (number) 0;
	for(int i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		this->_U += _particle_particle_bonded_interaction(p, true);

		p->prepare_list();
		neigh = p->next_neighbour();
		while(neigh != P_VIRTUAL) {
			this->_U += _particle_particle_interaction(p, &this->_particles[neigh]);
			neigh = p->next_neighbour();
		}
	}

	this->_U *= (number) 0.5;
	this->_U_hydr *= (number) 0.5;
	this->_U_stack *= (number) 0.5;
}

template<typename number>
inline void MC_CPUBackend<number>::_translate_particle(Particle<number> *p) {
	p->pos.x += (drand48() - (number)0.5)*this->_delta[MC_MOVE_TRANSLATION];
	p->pos.y += (drand48() - (number)0.5)*this->_delta[MC_MOVE_TRANSLATION];
	p->pos.z += (drand48() - (number)0.5)*this->_delta[MC_MOVE_TRANSLATION];

	if(p->pos_list.sqr_distance(p->pos) > this->_sqr_verlet_skin) this->_are_lists_old = true;
}

template<typename number>
inline void MC_CPUBackend<number>::_rotate_particle(Particle<number> *p) {
	number t = (drand48() - (number)0.5) * this->_delta[MC_MOVE_ROTATION];
	LR_vector<number> axis = Utils::get_random_vector<number>();

	number sintheta = sin(t);
	number costheta = cos(t);
	number olcos = ((number)1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix<number> R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
				xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
				xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	p->orientation = p->orientation * R;
}

template<typename number>
void MC_CPUBackend<number>::sim_step(llint curr_step) {
	LR_vector<number> tmp;

	get_time(&this->_timer, 0);

	for(int i = 0; i < this->_N; i++) {
		int pi = (int) (drand48() * this->_N);
		Particle<number> *p = &this->_particles[pi];

		int move = (drand48() < 0.5) ? MC_MOVE_TRANSLATION : MC_MOVE_ROTATION;

		this->_tries[move]++;
		number delta_E = -_particle_energy(p);
		p->set_ext_potential (curr_step);
		number delta_E_ext = -p->ext_potential;

		if(move == MC_MOVE_TRANSLATION) {
			_particles_old[pi].pos = p->pos;
			_translate_particle(p);
		}
		else {
			_particles_old[pi].orientation = p->orientation;
			_particles_old[pi].orientationT = p->orientationT;
			_rotate_particle(p);
			p->orientationT = p->orientation.get_transpose();
			p->set_positions();
		}

		get_time(&this->_timer, 2);
		if(this->_are_lists_old == true) {
			this->_update_lists();
			this->_N_updates++;
		}
		get_time(&this->_timer, 3);

		delta_E += _particle_energy(p);
		p->set_ext_potential(curr_step);
		delta_E_ext += p->ext_potential;

		// uncomment to check the energy at a given time step.
		// may be useful for debugging purposes
		//if (curr_step > 410000 && curr_step <= 420001)
		// printf("delta_E: %lf\n", (double)delta_E);

		if(!this->_overlap && ((delta_E + delta_E_ext) < 0 ||
				       exp(-(delta_E + delta_E_ext) / this->_T) > drand48())) {
			this->_accepted[move]++;
			this->_U += delta_E;
		}
		else {
			if(move == MC_MOVE_TRANSLATION) {
				p->pos = _particles_old[pi].pos;
			}
			else {
				p->orientation = _particles_old[pi].orientation;
				p->orientationT = _particles_old[pi].orientationT;
				p->set_positions();
			}
		}
		this->_overlap = false;
	}

	get_time(&this->_timer, 1);
	process_times(&this->_timer);
}

template<typename number>
void MC_CPUBackend<number>::print_energy(llint curr_step) {
	_compute_energy();
	MCBackend<number>::print_energy(curr_step);
}

template class MC_CPUBackend<float>;
template class MC_CPUBackend<double>;

