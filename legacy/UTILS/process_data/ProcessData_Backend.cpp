/*
 * MD_CPUBackend.cpp
 *
 *  Created on: 03/set/2010
 *      Author:
 */

#include "ProcessData_Backend.h"
#include "../../IOManager.h"

#define HBOND_CUTOFF (-0.5)

template<typename number>
ProcessData_Backend<number>::ProcessData_Backend(IOManager *IO,int confs_to_skip) :
	SimBackend<number> (IO) {
	this->_is_CUDA_sim = false;
	this->_confs_to_skip = confs_to_skip;
}

template<typename number>
ProcessData_Backend<number>::~ProcessData_Backend() {

}

template<typename number> ostream& ProcessData_Backend<number>::OutputBonds(
		ostream& out, int t) {

	number en_stack, en_coax, en_cstack, en_H_bond, en_exc, en_FENE, en_tot;
	char h_bond_status;
	number energy = 0; //total energy, for control purposes
	this->_update_lists(); //creates lists of interacting particles for each particle

	for (int i = 0; i < this->_N; i++) {
		this->_particles[i].set_positions();
		this->_particles[i].orientationT
				= this->_particles[i].orientation.get_transpose();
	}

	int neigh; //particle's interacting neighbors
	Particle<number> *p;

	out
			<< "# Particle_A  Particle_B Total_Interaction_E H_bond_Energy Stacking_Energy Cross_Stacking_Energy Coaxial_stacking_Energy Excluded_Vol_EN FENE_Spring H_bond_status "
			<< ", t=  " << t << endl;
	//cout << this->_N;
	number controlenergy = 0;
	number bondedenergy = 0;
	number intenergy = 0;
	//total energy in the system calculation
	for (int i = 0; i < this->_N; i++) {

		p = &this->_particles[i];
		en_tot = en_stack = en_coax = en_cstack = en_H_bond = en_exc = en_FENE
				= 0;
		h_bond_status = 'N';
        bondedenergy = 0;
		//first calculate all 3 and 5- neighbors on strands interactions
		int strand_neighbor = 0;

		if (p->n3 != P_VIRTUAL /* && i < p->n3 */) {

			strand_neighbor = 3;
			neigh = p->n3;
			en_stack = _particle_particle_bonded_stacking(p);

			en_exc = _particle_particle_bonded_excluded_volume(p);
			en_FENE = _particle_particle_bonded_FENE(p);
			en_tot = en_stack + en_exc + en_FENE;
			controlenergy += en_tot;
			bondedenergy += en_tot;

		}
		if (!(en_stack == 0 && en_coax == 0 && en_cstack == 0 && en_H_bond == 0
				&& en_exc == 0 && en_FENE == 0)) {
			if( i < neigh)
			{
				out << i << " " << neigh << " ";
			}
			else
			{
				out << neigh << " " << i << " ";
			}
			out << en_tot << " " << en_H_bond
								<< " " << en_stack << " " << en_cstack << " " << en_coax
								<< " " << en_exc << " " << en_FENE << " " << h_bond_status << '\n';
			//out << i << " " << neigh << " " << en_tot << " " << en_H_bond
		//			<< " " << en_stack << " " << en_cstack << " " << en_coax
			//		<< " " << en_exc << " " << en_FENE << '\n';
		}


		en_tot = en_stack = en_coax = en_cstack = en_H_bond = en_exc = en_FENE
				= 0;

	
#ifdef DBG
				if( bondedenergy != _particle_particle_bonded_interaction(p))
											cout << "Error in bonded particle " << i << " and " << neigh << ": " << bondedenergy - _particle_particle_bonded_interaction(p) << '\n';
#endif


		energy += mc_particle_particle_bonded_interaction(p); //this is for control


		p->prepare_list();
		neigh = p->next_neighbour();
		while (neigh != P_VIRTUAL) {
			Particle<number> *q = &this->_particles[neigh];
			energy += mc_particle_particle_interaction(p,
					&this->_particles[neigh]); //this is for control

			en_tot = en_stack = en_coax = en_cstack = en_H_bond = en_exc
					= en_FENE = 0;
			h_bond_status = 'N';
			if (i < neigh) //only in this case we will print out the interaction
			{
				en_coax = _particle_particle_coaxial_stacking(p, q);
				en_cstack = _particle_particle_cross_stacking(p, q);

				en_H_bond = _particle_particle_H_bond(p, q);
				if(en_H_bond < HBOND_CUTOFF)
					h_bond_status = 'H';
				else 
					h_bond_status = _particle_particle_H_bond_status(p,q);
				en_exc = _particle_particle_excluded_volume(p, q);
				en_tot = en_coax + en_cstack + en_H_bond + en_exc;
              controlenergy += en_tot;
				if (!(en_stack == 0 && en_coax == 0 && en_cstack == 0
						&& en_H_bond == 0 && en_exc == 0 && en_FENE == 0 && h_bond_status == 'N')) {
					out << i << " " << neigh << " " << en_tot << " "
							<< en_H_bond << " " << en_stack << " " << en_cstack
							<< " " << en_coax << " " << en_exc << " "
							<< en_FENE << ' ' << h_bond_status << '\n';
				}
#ifdef DBG
				if(en_tot != _particle_particle_interaction(p,
						&this->_particles[neigh]) )
											cout << "Error in particle " << i << " and " << neigh << ": " << en_tot - _particle_particle_interaction(p,
													&this->_particles[neigh]) << '\n';
#endif
			}
			neigh = p->next_neighbour();
		}

	}

	out << "Control energy: " << controlenergy / (this->_N) << '\n';
	this->_compute_forces();
	out << "Pot Energy should be " << this->_U / this->_N << '\n';
	return out;
}

//--------------------------------------------------------------------------------------------
template<typename number> char ProcessData_Backend<number>::_particle_particle_H_bond_status(
		Particle<number> *p, Particle<number> *q) { //return 1 if bonded, 0 if not bonded, according to old criteria (avg h-bond)

	number energy = 0;
	LR_vector<number> force , torquep, torqueq;
	char h_bond_status = 'N';
	bool is_pair = (q->type + p->type == 3);

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
		number rhydromod = rhydro.module();
		if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
			// vector, versor and magnitude of the base-base separation
		  	LR_vector<number> rhydrodir = rhydro / rhydromod;

			// particle axes according to Allen's paper
			LR_vector<number> a1 = p->orientationT.v1;
			LR_vector<number> a2 = p->orientationT.v2;
			LR_vector<number> a3 = p->orientationT.v3;
			LR_vector<number> b1 = q->orientationT.v1;
			LR_vector<number> b2 = q->orientationT.v2;
			LR_vector<number> b3 = q->orientationT.v3;

			// angles involved in the HB interaction
			  /*number t1 = LRACOS (-a1 * b1);
			number t2 = LRACOS (-b1 * rhydrodir);
			number t3 = LRACOS ( a1 * rhydrodir);
			number t4 = LRACOS ( a3 * b3);
			number t7 = LRACOS (-b3 * rhydrodir);
			number t8 = LRACOS ( a3 * rhydrodir);
			*/
			number cost1 = -a1 * b1;
			number cost2 = -b1 * rhydrodir;
			number cost3 =  a1 * rhydrodir;

			number cost4 = a3 * b3;
			number cost7 = -b3 * rhydrodir;
			number cost8 = a3 * rhydrodir;

		 	 // functions called at their relevant arguments
			//number f1   = this->_interaction.f1(rhydromod, HYDR_F1);
		  	/*
			number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
		  	number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
		  	number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
			*/
			number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
			number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[HYDR_F4_THETA1]);
			number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[HYDR_F4_THETA2]);
			number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[HYDR_F4_THETA3]);

			number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[HYDR_F4_THETA4]);
			number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[HYDR_F4_THETA7]);
			number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[HYDR_F4_THETA8]);
		  	/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
		  	number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
		  	number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

			number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
			
			if(1.077*hb_energy/f1 < -0.1)
			{
				h_bond_status = 'H';
			}

			energy += hb_energy;
			this->_U_hydr += hb_energy;

			
		}
		// END OF HYDROGEN BONDING
	return h_bond_status;
}


//-----------------------------------------------------------------------------------------------------
template<typename number> number ProcessData_Backend<number>::_particle_particle_H_bond(
		Particle<number> *p, Particle<number> *q) {

	number energy = 0;
	LR_vector<number> force , torquep, torqueq;

	bool is_pair = (q->btype + p->btype == 3);

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
		number rhydromod = rhydro.module();
		if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
			// vector, versor and magnitude of the base-base separation
		  	LR_vector<number> rhydrodir = rhydro / rhydromod;

			// particle axes according to Allen's paper
			LR_vector<number> a1 = p->orientationT.v1;
			LR_vector<number> a2 = p->orientationT.v2;
			LR_vector<number> a3 = p->orientationT.v3;
			LR_vector<number> b1 = q->orientationT.v1;
			LR_vector<number> b2 = q->orientationT.v2;
			LR_vector<number> b3 = q->orientationT.v3;

			// angles involved in the HB interaction
			  /*number t1 = LRACOS (-a1 * b1);
			number t2 = LRACOS (-b1 * rhydrodir);
			number t3 = LRACOS ( a1 * rhydrodir);
			number t4 = LRACOS ( a3 * b3);
			number t7 = LRACOS (-b3 * rhydrodir);
			number t8 = LRACOS ( a3 * rhydrodir);
			*/
			number cost1 = -a1 * b1;
			number cost2 = -b1 * rhydrodir;
			number cost3 =  a1 * rhydrodir;

			number cost4 = a3 * b3;
			number cost7 = -b3 * rhydrodir;
			number cost8 = a3 * rhydrodir;

		 	 // functions called at their relevant arguments
			//number f1   = this->_interaction.f1(rhydromod, HYDR_F1);
		  	/*
			number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
		  	number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
		  	number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
			*/
			number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
			number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[HYDR_F4_THETA1]);
			number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[HYDR_F4_THETA2]);
			number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[HYDR_F4_THETA3]);

			number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[HYDR_F4_THETA4]);
			number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[HYDR_F4_THETA7]);
			number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[HYDR_F4_THETA8]);
		  	/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
		  	number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
		  	number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

			number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
			energy += hb_energy;
			this->_U_hydr += hb_energy;

			// makes sense, since the above functions may return 0. exactly
			if (hb_energy != 0.) {
				// derivatives called at the relevant arguments
				number f1D      =  this->_interaction.f1D(rhydromod, HYDR_F1,q->type, p->type);
				/*number f4t1Dsin = -this->_interaction.f4Dsin(t1, HYDR_F4_THETA1);
				number f4t2Dsin = -this->_interaction.f4Dsin(t2, HYDR_F4_THETA2);
				number f4t3Dsin =  this->_interaction.f4Dsin(t3, HYDR_F4_THETA3);*/
				number f4t1Dsin = this->_interaction.query_meshD (cost1, this->_interaction.mesh_f4[HYDR_F4_THETA1]);
				number f4t2Dsin = this->_interaction.query_meshD (cost2, this->_interaction.mesh_f4[HYDR_F4_THETA2]);
				number f4t3Dsin = -this->_interaction.query_meshD (cost3, this->_interaction.mesh_f4[HYDR_F4_THETA3]);

				number f4t4Dsin = -this->_interaction.query_meshD (cost4, this->_interaction.mesh_f4[HYDR_F4_THETA4]);
				number f4t7Dsin = this->_interaction.query_meshD (cost7, this->_interaction.mesh_f4[HYDR_F4_THETA7]);
				number f4t8Dsin = -this->_interaction.query_meshD (cost8, this->_interaction.mesh_f4[HYDR_F4_THETA8]);
				/*number f4t4Dsin =  this->_interaction.f4Dsin(t4, HYDR_F4_THETA4);
				number f4t7Dsin = -this->_interaction.f4Dsin(t7, HYDR_F4_THETA7);
				number f4t8Dsin =  this->_interaction.f4Dsin(t8, HYDR_F4_THETA8);*/

				//f4t1 = 1.;
				//f4t2 = 1.;
				//f4t3 = 1.;
				//f4t4 = 1.;
				//f4t7 = 1.;
				//f4t8 = 1.;

				//f4t1Dsin = 0.;
				//f4t2Dsin = 0.;
				//f4t3Dsin = 0.;
				//f4t4Dsin = 0.;
				//f4t7Dsin = 0.;
				//f4t8Dsin = 0.;

				// RADIAL PART
				force = - rhydrodir * f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

				// TETA4; t4 = LRACOS (a3 * b3);
				LR_vector<number> dir = a3.cross(b3);
				number torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8 ;

				torquep -= dir * torquemod;
				torqueq += dir * torquemod;

				// TETA1; t1 = LRACOS (-a1 * b1);
				dir = a1.cross(b1);
				torquemod = - f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 ;

				torquep -= dir * torquemod;
				torqueq += dir * torquemod;

				// TETA2; t2 = LRACOS (-b1 * rhydrodir);
				number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
				//force += (b1 + rhydrodir * cos(t2)) / rhydromod * fact;
				force += (b1 + rhydrodir * cost2) / rhydromod * fact;
				dir = rhydrodir.cross(b1);
				//torquemod = - f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

				torqueq += -dir * fact;

				// TETA3; t3 = LRACOS (a1 * rhydrodir);
				//force += (a1 - rhydrodir * cos(t3)) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
				force += (a1 - rhydrodir * cost3) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

				LR_vector<number> t3dir = rhydrodir.cross(a1);
				torquemod = - f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

				torquep += t3dir * torquemod;

				// THETA7; t7 = LRACOS (-rhydrodir * b3);
				//force += (b3 + rhydrodir * cos(t7)) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;
				force += (b3 + rhydrodir * cost7) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

				LR_vector<number> t7dir = rhydrodir.cross(b3);
				torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

				torqueq += t7dir * torquemod;

				// THETA 8; t8 = LRACOS (rhydrodir * a3);
				//force +=  (a3 - rhydrodir * cos(t8)) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
				force +=  (a3 - rhydrodir * cost8) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

				LR_vector<number> t8dir = rhydrodir.cross(a3);
				torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

				torquep += t8dir * torquemod;

				// final update of forces and torques for HB
				p->force -= force;
				q->force += force;

				torquep -= p->pos_base.cross(force);
				torqueq += q->pos_base.cross(force);
			}
		}
		// END OF HYDROGEN BONDING
	return energy;
}

template<typename number> number ProcessData_Backend<number>::_particle_particle_excluded_volume(
		Particle<number> *p, Particle<number> *q) {
	bool is_pair = (q->type + p->type == 3);

		LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

		number energy = 0;

		// excluded volume

		// BASE-BASE
		LR_vector<number> force(0, 0, 0);
		LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
		energy += _excluded_volume(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
		LR_vector<number> torquep = -p->pos_base.cross(force);
		LR_vector<number> torqueq = q->pos_base.cross(force);

		p->force -= force;
		q->force += force;

		// P-BASE vs. Q-BACK
		rcenter = r + q->pos_back - p->pos_base;
		energy += _excluded_volume(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
		torquep += -p->pos_base.cross(force);
		torqueq += q->pos_back.cross(force);

		p->force -= force;
		q->force += force;

		// P-BACK vs. Q-BASE
		rcenter = r + q->pos_base - p->pos_back;
		energy += _excluded_volume(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
		torquep += -p->pos_back.cross(force);
		torqueq +=  q->pos_base.cross(force);

		p->force -= force;
		q->force += force;

		// BACK-BACK
		rcenter = r + q->pos_back - p->pos_back;
		energy += _excluded_volume(rcenter, force, EXCL_S1, EXCL_R1, EXCL_B1, EXCL_RC1);
		torquep += -p->pos_back.cross(force);
		torqueq +=  q->pos_back.cross(force);

		p->force -= force;
		q->force += force;
 return energy;
}

template<typename number> number ProcessData_Backend<number>::_particle_particle_bonded_stacking(
		Particle<number> *p) {

	number energy = (number) 0;
	LR_vector<number> force , torquep, torqueq;

		if(p->n3 == P_VIRTUAL) {
			return 0;
		}

		Particle<number> *q = &this->_particles[p->n3];
				LR_vector<number> a1 = p->orientationT.v1;
				LR_vector<number> a2 = p->orientationT.v2;
				LR_vector<number> a3 = p->orientationT.v3;
				LR_vector<number> b1 = q->orientationT.v1;
				LR_vector<number> b2 = q->orientationT.v2;
				LR_vector<number> b3 = q->orientationT.v3;
				LR_vector<number> r = q->pos - p->pos;

				LR_vector<number> rback = r + q->pos_back - p->pos_back;
						number rbackmod = rback.module();

	// STACKING
			LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
			number rstackmod = rstack.module();
			LR_vector<number> rstackdir = rstack / rstackmod;

			//number t4 = LRACOS( a3 * b3);
			//number t5 = LRACOS( a3 * rstackdir);
			//number t6 = LRACOS(-b3 * rstackdir);
			number cost4 = a3 * b3;
			number cost5 = a3 * rstackdir;
			number cost6 =-b3 * rstackdir;
			number cosphi1 = a2 * rback / rbackmod;
			number cosphi2 = b2 * rback / rbackmod;

			// functions and their derivatives needed for energies and forces
			//number f1     = this->_interaction.f1(rstackmod, STCK_F1);
			number f1     = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
			//number f4t4   = this->_interaction.f4(t4, STCK_F4_THETA4);
			//number f4t5   = this->_interaction.f4(PI - t5, STCK_F4_THETA5);
			//number f4t6   = this->_interaction.f4(t6, STCK_F4_THETA6);
			number f4t4   = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[STCK_F4_THETA4]);
			number f4t5   = this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[STCK_F4_THETA5]);
			number f4t6   = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[STCK_F4_THETA6]);
			number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
			number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

			// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
			// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
			//number f1D      = this->_interaction.f1D(rstackmod, STCK_F1);
			//number f4t4Dsin = this->_interaction.f4Dsin(t4, STCK_F4_THETA4);
			//number f4t5Dsin = this->_interaction.f4Dsin(PI - t5, STCK_F4_THETA5);
			//number f4t6Dsin = this->_interaction.f4Dsin(t6, STCK_F4_THETA6);
			//number f4t4Dsin = -this->_interaction.query_meshD (cost4, this->_interaction.mesh_f4[STCK_F4_THETA4]);
			//number f4t5Dsin = -this->_interaction.query_meshD (-cost5, this->_interaction.mesh_f4[STCK_F4_THETA5]);
			//number f4t6Dsin = -this->_interaction.query_meshD (cost6, this->_interaction.mesh_f4[STCK_F4_THETA6]);
			//number f5phi1D  = this->_interaction.f5D(cosphi1, STCK_F5_PHI1);
			//number f5phi2D  = this->_interaction.f5D(cosphi2, STCK_F5_PHI2);

			number stck_energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
			energy += stck_energy;

			// RADIAL
			//force = -rstackdir * f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

			// THETA 5
			//force += -(a3 - rstackdir * cos(t5)) / rstackmod * f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;
			//force += -(a3 - rstackdir * cost5) / rstackmod * f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;

			// THETA 6
			//force += -(b3 + rstackdir * cos(t6)) / rstackmod * f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;
			//force += -(b3 + rstackdir * cost6) / rstackmod * f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;

			// COS PHI 1
			// here particle p is referred to using the a while particle q is referred with the b
			//number gamma = POS_STACK - POS_BACK;
			//number rbackmodcub = rbackmod * rbackmod * rbackmod;

			//number ra2 = rstackdir*a2;
			//number ra1 = rstackdir*a1;
			//number rb1 = rstackdir*b1;
			/*number a2b1 = a2*b1;

			number parentesi = rstackmod*ra2 - a2b1*gamma;

			number dcosphi1dr    = (SQR(rstackmod)*ra2 - ra2*SQR(rbackmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*gamma + a2b1*(-ra1 + rb1)*SQR(gamma))/ rbackmodcub;
			number dcosphi1dra1  =  rstackmod * gamma * parentesi / rbackmodcub;
			number dcosphi1dra2  = -rstackmod / rbackmod;
			number dcosphi1drb1  = -rstackmod * gamma * parentesi / rbackmodcub;
			number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackmodcub;
			number dcosphi1da2b1 =  gamma / rbackmod;

			// this force part has a minus because of the definition of cos(phi1)
			// which is not consistent with the derivatives above (i.e. all the
			// derivatives should have a minus sign in front and the force below shouldn't)
			number force_part_phi1 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1D * f5phi2;

			force += -(rstackdir * dcosphi1dr +
					   ((a2 - rstackdir * ra2) * dcosphi1dra2 +
						(a1 - rstackdir * ra1) * dcosphi1dra1 +
						(b1 - rstackdir * rb1) * dcosphi1drb1) / rstackmod) * force_part_phi1;

			// COS PHI 2
			// here particle p -> b, particle q -> a
			ra2 = rstackdir*b2;
			ra1 = rstackdir*b1;
			rb1 = rstackdir*a1;
			a2b1 = b2*a1;

			parentesi = rstackmod*ra2 + a2b1*gamma;

			number dcosphi2dr    =  (parentesi * (rstackmod + (rb1 - ra1)*gamma) - ra2*SQR(rbackmod)) / rbackmodcub;
			number dcosphi2dra1  = -rstackmod*gamma*(rstackmod*ra2 + a2b1*gamma) / rbackmodcub;
			number dcosphi2dra2  = -rstackmod / rbackmod;
			number dcosphi2drb1  =  rstackmod*gamma* parentesi / rbackmodcub;
			number dcosphi2da1b1 = -SQR(gamma)* parentesi / rbackmodcub;
			number dcosphi2da2b1 = -gamma / rbackmod;

			// this force part has a minus because of the definition of cos(phi2) which is not consistent with the derivatives
			// above (i.e. all the derivatives should have a minus sign in front and the force below shouldn't)
			number force_part_phi2 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2D;

			force += - force_part_phi2 * (rstackdir * dcosphi2dr +
										  ((b2 - rstackdir * ra2) * dcosphi2dra2 +
										   (b1 - rstackdir * ra1) * dcosphi2dra1 +
										   (a1 - rstackdir * rb1) * dcosphi2drb1) / rstackmod);

			// Add the contribution from all the forces to the stored particles' forces
			p->force -= force;
			q->force += force;

			torquep -= p->pos_stack.cross(force);
			torqueq += q->pos_stack.cross(force);

			// handle the part on theta4
			LR_vector<number> t4dir = b3.cross(a3);
			number torquemod = f1 * f4t4Dsin * f4t5 * f4t6 * f5phi1 * f5phi2;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// handle the part on theta5
			LR_vector<number> t5dir = rstackdir.cross(a3);
			torquemod = -f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;

			torquep -= t5dir * torquemod;

			// handle the part on theta6
			LR_vector<number> t6dir = rstackdir.cross(b3);
			torquemod = f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;

			torqueq += t6dir * torquemod;

			// PHI 1
			torquep += rstackdir.cross(a2) * force_part_phi1 * dcosphi1dra2 +
					rstackdir.cross(a1) * force_part_phi1 * dcosphi1dra1;
			torqueq += rstackdir.cross(b1) * force_part_phi1 * dcosphi1drb1;

			LR_vector<number> puretorque = a2.cross(b1) * force_part_phi1 * dcosphi1da2b1 +
					a1.cross(b1) * force_part_phi1 * dcosphi1da1b1;

			torquep -= puretorque;
			torqueq += puretorque;

			// PHI 2
			torquep += rstackdir.cross(a1) * force_part_phi2 * dcosphi2drb1;
			torqueq += rstackdir.cross(b2) * force_part_phi2 * dcosphi2dra2 +
					rstackdir.cross(b1) * force_part_phi2 * dcosphi2dra1;

			puretorque = a1.cross(b2) * force_part_phi2 * dcosphi2da2b1 +
					a1.cross(b1) * force_part_phi2 * dcosphi2da1b1;

			torquep -= puretorque;
			torqueq += puretorque;

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
*/
	return energy;

}

template<typename number> number ProcessData_Backend<number>::_particle_particle_coaxial_stacking(
		Particle<number> *p, Particle<number> *q) {
	bool is_pair = (q->type + p->type == 3);
	LR_vector<number> force , torquep, torqueq;

			LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	number energy = 0;
	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
		number rstackmod = rstack.module();
		if(CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		  	LR_vector<number> rstackdir = rstack / rstackmod;

			// particle axes according to Allen's paper
			LR_vector<number> a1 = p->orientationT.v1;
			LR_vector<number> a2 = p->orientationT.v2;
			LR_vector<number> a3 = p->orientationT.v3;
			LR_vector<number> b1 = q->orientationT.v1;
			LR_vector<number> b2 = q->orientationT.v2;
			LR_vector<number> b3 = q->orientationT.v3;

			// angles involved in the CXST interaction
			/*
			number t1 = LRACOS (-a1 * b1);
			number t4 = LRACOS ( a3 * b3);
			number t5 = LRACOS ( a3 * rstackdir);
			number t6 = LRACOS (-b3 * rstackdir);
			*/
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
			/*
		  	number f4t1 = this->_interaction.f4(t1, CXST_F4_THETA1) + this->_interaction.f4(2 * PI - t1, CXST_F4_THETA1);
		  	number f4t4 = this->_interaction.f4(t4, CXST_F4_THETA4);
		  	number f4t5 = this->_interaction.f4(t5, CXST_F4_THETA5) + this->_interaction.f4(PI - t5, CXST_F4_THETA5);
		  	number f4t6 = this->_interaction.f4(t6, CXST_F4_THETA6) + this->_interaction.f4(PI - t6, CXST_F4_THETA6);
			*/
		  	number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[CXST_F4_THETA1]);
		  	number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[CXST_F4_THETA4]);
		  	number f4t5 = this->_interaction.query_mesh (cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]) + this->_interaction.query_mesh (-cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]);
		  	number f4t6 = this->_interaction.query_mesh (cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]) + this->_interaction.query_mesh (-cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]);
			number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

			number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
			energy += cxst_energy;

			/*
			// again, makes sense (see same line for HB)
			if(cxst_energy != 0.) {
				// derivatives called at the relevant arguments
				number f2D      =  this->_interaction.f2D(rstackmod, CXST_F2);
				/*
				number f4t1Dsin = -this->_interaction.f4Dsin(t1, CXST_F4_THETA1) + this->_interaction.f4Dsin(2 * PI - t1, CXST_F4_THETA1);
				number f4t4Dsin =  this->_interaction.f4Dsin(t4, CXST_F4_THETA4);
				number f4t5Dsin =  this->_interaction.f4Dsin(t5, CXST_F4_THETA5) - this->_interaction.f4Dsin(PI - t5, CXST_F4_THETA5);
				number f4t6Dsin = -this->_interaction.f4Dsin(t6, CXST_F4_THETA6) + this->_interaction.f4Dsin(PI - t6, CXST_F4_THETA6);
		  		number f4t1Dsin = this->_interaction.query_meshD (cost1, this->_interaction.mesh_f4[CXST_F4_THETA1]);
		  		number f4t4Dsin = -this->_interaction.query_meshD (cost4, this->_interaction.mesh_f4[CXST_F4_THETA4]);
		  		number f4t5Dsin = -this->_interaction.query_meshD (cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]) + this->_interaction.query_meshD (-cost5, this->_interaction.mesh_f4[CXST_F4_THETA5]);
		  		number f4t6Dsin =  this->_interaction.query_meshD (cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]) - this->_interaction.query_meshD (-cost6, this->_interaction.mesh_f4[CXST_F4_THETA6]);

				number f5Dcosphi3 = this->_interaction.f5D(cosphi3, CXST_F5_PHI3);

				//f4t1 = 1.;
				//f4t4 = 1.;
				//f4t5 = 1.;
				//f4t6 = 1.;
				//f5cosphi3 = 1.;

				//f4t1Dsin = 0.;
				//f4t4Dsin = 0.;
				//f4t5Dsin = 0.;
				//f4t6Dsin = 0.;
				//f5Dcosphi3 = 0.;
*/
				/*FILE * yo = fopen("yo.dat", "w");
				for (int kk = 0; kk<10000; kk++)
				  {
					fprintf(yo,"%lf %lf %lf\n", 0.0002 * (kk - 5000) * 2 * PI, this->_interaction.f5(0.0002 * (kk-5000) * 2 * PI, CXST_F5_PHI3), this->_interaction.f5D(0.0002 * (kk - 5000) * 2 * PI, CXST_F5_PHI3));
				  }
				fclose(yo);
				exit(-1);*/
/*
				// RADIAL PART
				force = - rstackdir * f2D * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

				// THETA1; t1 = LRACOS (-a1 * b1);
				LR_vector<number> dir = a1.cross(b1);
				number torquemod = - f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);

				torquep -= dir * torquemod;
				torqueq += dir * torquemod;

				// TETA4; t4 = LRACOS (a3 * b3);
				dir = a3.cross(b3);
				torquemod = - f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * SQR(f5cosphi3);

				torquep -= dir * torquemod;
				torqueq += dir * torquemod;

				// THETA5; t5 = LRACOS ( a3 * rstackdir);
				number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * SQR(f5cosphi3);
				//force += fact * (a3 - rstackdir * cos(t5)) / rstackmod;
				force += fact * (a3 - rstackdir * cost5) / rstackmod;
				dir = rstackdir.cross(a3);
				torquep -= dir * fact;

				// THETA6; t6 = LRACOS (-b3 * rstackdir);
				fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * SQR(f5cosphi3);
				//force += (b3 + rstackdir * cos(t6)) / rstackmod * fact;
				force += (b3 + rstackdir * cost6) / rstackmod * fact;
				dir = rstackdir.cross(b3);

				torqueq += - dir * fact;

				// Cosphi3 (qui son dolori...) (meno male che cosphi4 = cosphi3)
				// Definition used:
				// cosphi3 = gamma * (ra3 * a2b1 - ra2 * a3b1)/ rbackmod
				number gamma = POS_STACK - POS_BACK;
				number gammacub = gamma * gamma * gamma;
				number rbackmodcub = rbackmod * rbackmod * rbackmod;
				//number a1b1 = a1 * b1;
				number a2b1 = a2 * b1;
				number a3b1 = a3 * b1;
				number ra1 = rstackdir * a1;
				number ra2 = rstackdir * a2;
				number ra3 = rstackdir * a3;
				number rb1 = rstackdir * b1;

				number parentesi = (ra3 * a2b1 - ra2 * a3b1);

				number dcdr    = -gamma * parentesi * (gamma * (ra1 - rb1) + rstackmod) / rbackmodcub;
				number dcda1b1 =  gammacub * parentesi / rbackmodcub;
				number dcda2b1 =  gamma * ra3 / rbackmod;
				number dcda3b1 = -gamma * ra2 / rbackmod;
				number dcdra1  = -SQR(gamma) * parentesi * rstackmod / rbackmodcub;
				number dcdra2  = -gamma * a3b1 / rbackmod;
				number dcdra3  =  gamma * a2b1 / rbackmod;
				number dcdrb1  =  SQR(gamma) * parentesi * rstackmod / rbackmodcub;

				number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2 * f5cosphi3 * f5Dcosphi3;

				force += - force_c * (rstackdir * dcdr +
							 ((a1 - rstackdir * ra1) * dcdra1 +
							  (a2 - rstackdir * ra2) * dcdra2 +
							  (a3 - rstackdir * ra3) * dcdra3 +
							  (b1 - rstackdir * rb1) * dcdrb1) / rstackmod);

				torquep += force_c * (rstackdir.cross(a1) * dcdra1 +
							  rstackdir.cross(a2) * dcdra2 +
							  rstackdir.cross(a3) * dcdra3);
				torqueq += force_c * (rstackdir.cross(b1) * dcdrb1);

				LR_vector<number> puretorque = force_c * (a1.cross(b1) * dcda1b1 +
									  a2.cross(b1) * dcda2b1 +
									  a3.cross(b1) * dcda3b1);
				torquep -= puretorque;
				torqueq += puretorque;

				// final update of forces and torques for CXST
				p->force -= force;
				q->force += force;

				torquep -= p->pos_stack.cross(force);
				torqueq += q->pos_stack.cross(force);
			}
		}

		// Final update of the torques
		// total torques
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;

*/
		}
	return energy;
}

template<typename number> number ProcessData_Backend<number>::_particle_particle_cross_stacking(
		Particle<number> *p, Particle<number> *q) {
	LR_vector<number> force , torquep, torqueq;

	number energy = 0;
	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();
	// CROSS STACKING
	// CROSS STACKING
		LR_vector<number> rcstack = rhydro;
		//LR_vector<number> rcstack = r + q->pos_base - p->pos_base;
		number rcstackmod = rhydromod;
		//number rcstackmod = rcstack.module();
		if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
		  	LR_vector<number> rcstackdir = rcstack / rcstackmod;

			// particle axes according to Allen's paper
			LR_vector<number> a1 = p->orientationT.v1;
			LR_vector<number> a2 = p->orientationT.v2;
			LR_vector<number> a3 = p->orientationT.v3;
			LR_vector<number> b1 = q->orientationT.v1;
			LR_vector<number> b2 = q->orientationT.v2;
			LR_vector<number> b3 = q->orientationT.v3;

			// angles involved in the CRST interaction
			/*
			number t1 = LRACOS (-a1 * b1);
			number t2 = LRACOS (-b1 * rcstackdir);
			number t4 = LRACOS ( a3 * b3);
			number t3 = LRACOS ( a1 * rcstackdir);
			number t7 = LRACOS (-rcstackdir * b3);
			number t8 = LRACOS ( rcstackdir * a3);
			*/
			number cost1 = -a1 * b1;
			number cost2 = -b1 * rcstackdir;
			number cost3 =  a1 * rcstackdir;
			number cost4 =  a3 * b3;
			number cost7 = -b3 * rcstackdir;
			number cost8 =  a3 * rcstackdir;

		 	 // functions called at their relevant arguments
			number f2   = this->_interaction.f2(rcstackmod, CRST_F2);
			/*
		  	number f4t1 = this->_interaction.f4(t1, CRST_F4_THETA1);
		  	number f4t2 = this->_interaction.f4(t2, CRST_F4_THETA2);
		  	number f4t3 = this->_interaction.f4(t3, CRST_F4_THETA3);
		  	number f4t4 = this->_interaction.f4(t4, CRST_F4_THETA4) + this->_interaction.f4(PI - t4, CRST_F4_THETA4);
			*/
		  	number f4t1 = this->_interaction.query_mesh (cost1, this->_interaction.mesh_f4[CRST_F4_THETA1]);
		  	number f4t2 = this->_interaction.query_mesh (cost2, this->_interaction.mesh_f4[CRST_F4_THETA2]);
		  	number f4t3 = this->_interaction.query_mesh (cost3, this->_interaction.mesh_f4[CRST_F4_THETA3]);
		  	number f4t4 = this->_interaction.query_mesh (cost4, this->_interaction.mesh_f4[CRST_F4_THETA4]) + this->_interaction.query_mesh (-cost4, this->_interaction.mesh_f4[CRST_F4_THETA4]);
			/*
		  	number f4t7 = this->_interaction.f4(t7, CRST_F4_THETA7) + this->_interaction.f4(PI - t7, CRST_F4_THETA7);
		  	number f4t8 = this->_interaction.f4(t8, CRST_F4_THETA8) + this->_interaction.f4(PI - t8, CRST_F4_THETA8);
			*/
		  	number f4t7 = this->_interaction.query_mesh (cost7, this->_interaction.mesh_f4[CRST_F4_THETA7]) + this->_interaction.query_mesh (-cost7, this->_interaction.mesh_f4[CRST_F4_THETA7]);
	;
		  	number f4t8 = this->_interaction.query_mesh (cost8, this->_interaction.mesh_f4[CRST_F4_THETA8]) + this->_interaction.query_mesh (-cost8, this->_interaction.mesh_f4[CRST_F4_THETA8]);
	;

			number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
			energy += cstk_energy;

			// makes sense since the above f? can return exacly 0.
		}

	return energy;

}

template<typename number> number ProcessData_Backend<number>::_particle_particle_bonded_excluded_volume(
		Particle<number> *p) {

	  number energy =0;
	  LR_vector<number> force , torquep, torqueq;

	if(p->n3 == P_VIRTUAL) {
		return 0;
	}


			Particle<number> *q = &this->_particles[p->n3];
			LR_vector<number> a1 = p->orientationT.v1;
			LR_vector<number> a2 = p->orientationT.v2;
			LR_vector<number> a3 = p->orientationT.v3;
			LR_vector<number> b1 = q->orientationT.v1;
			LR_vector<number> b2 = q->orientationT.v2;
			LR_vector<number> b3 = q->orientationT.v3;
			LR_vector<number> r = q->pos - p->pos;


			// excluded volume

					// BASE-BASE
					force.x = force.y = force.z = (number) 0;
					LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
					energy += _excluded_volume(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2, EXCL_RC2);
					torquep -= p->pos_base.cross(force);
					torqueq += q->pos_base.cross(force);

					p->force -= force;
					q->force += force;

					// P-BASE vs. Q-BACK
					rcenter = r + q->pos_back - p->pos_base;
					energy += _excluded_volume(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3, EXCL_RC3);
					torquep -= p->pos_base.cross(force);
					torqueq += q->pos_back.cross(force);

					p->force -= force;
					q->force += force;

					// P-BACK vs. Q-BASE
					rcenter = r + q->pos_base - p->pos_back;
					energy += _excluded_volume(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4, EXCL_RC4);
					torquep -= p->pos_back.cross(force);
					torqueq += q->pos_base.cross(force);

					p->force -= force;
					q->force += force;



	return energy;
}

template<typename number> number ProcessData_Backend<number>::_particle_particle_bonded_FENE(
		Particle<number> *p) {
	number energy = (number) 0;

	if(p->n3 == P_VIRTUAL) {
		return 0;
	}


			Particle<number> *q = &this->_particles[p->n3];
			LR_vector<number> a1 = p->orientationT.v1;
			LR_vector<number> a2 = p->orientationT.v2;
			LR_vector<number> a3 = p->orientationT.v3;
			LR_vector<number> b1 = q->orientationT.v1;
			LR_vector<number> b2 = q->orientationT.v2;
			LR_vector<number> b3 = q->orientationT.v3;
			LR_vector<number> r = q->pos - p->pos;


			// FENE
			LR_vector<number> rback = r + q->pos_back - p->pos_back;
			number rbackmod = rback.module();
			number rbackr0 = rbackmod - FENE_R0;
			energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);
			LR_vector<number> force = rback * (-(FENE_EPS * rbackr0  / (FENE_DELTA2 - SQR(rbackr0))) / rbackmod);

			p->force -= force;
			q->force += force;

			LR_vector<number> torquep = -p->pos_back.cross(force);
			LR_vector<number> torqueq = q->pos_back.cross(force);

	return energy;
}

template<typename number>
inline number ProcessData_Backend<number>::_excluded_volume(const LR_vector<
		number> &r, number sigma, number rstar, number b, number rc) {
	number rmod = r.module();
	number energy = 0;

	if (rmod < rc) {
		if (rmod > rstar) {
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
		} else {
			number lj_part = SQR(sigma / rmod) * SQR(sigma / rmod)
					* SQR(sigma / rmod);
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}

	return energy;
}
/*
template<typename number>
void ProcessData_Backend<number>::_first_step(llint cur_step) {
	for (int i = 0; i < this->_N; i++) {
		this->_particles[i].vel += this->_particles[i].force * this->_dt
				* (number) 0.5;
		this->_particles[i].L += this->_particles[i].torque * this->_dt
				* (number) 0.5;

		// update of the orientation
		number norm = this->_particles[i].L.module();
		LR_vector<number> LVersor(this->_particles[i].L / norm);

		number sintheta = sin(this->_dt * norm);
		number costheta = cos(this->_dt * norm);
		number olcos = 1. - costheta;

		number xyo = LVersor[0] * LVersor[1] * olcos;
		number xzo = LVersor[0] * LVersor[2] * olcos;
		number yzo = LVersor[1] * LVersor[2] * olcos;
		number xsin = LVersor[0] * sintheta;
		number ysin = LVersor[1] * sintheta;
		number zsin = LVersor[2] * sintheta;

		LR_matrix<number> R(LVersor[0] * LVersor[0] * olcos + costheta, xyo
				- zsin, xzo + ysin, xyo + zsin, LVersor[1] * LVersor[1] * olcos
				+ costheta, yzo - xsin, xzo - ysin, yzo + xsin, LVersor[2]
				* LVersor[2] * olcos + costheta);

		this->_particles[i].orientation = this->_particles[i].orientation * R;
		this->_particles[i].pos += this->_particles[i].vel * this->_dt;
		// set back, base and stack positions
		this->_particles[i].set_positions();
		this->_particles[i].orientationT
				= this->_particles[i].orientation.get_transpose();

		this->_particles[i].set_initial_forces(cur_step);
		this->_particles[i].torque = LR_vector<number> ((number) 0, (number) 0,
				(number) 0);

		if (this->_particles[i].pos_list.sqr_distance(this->_particles[i].pos)
				> this->_sqr_verlet_skin)
			this->_are_lists_old = true;
	}
}
*/
template<typename number>
inline number ProcessData_Backend<number>::_particle_particle_bonded_interaction(
		Particle<number> *p) {
	number energy = (number) 0;

	if (p->n3 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n3];
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = q->pos - p->pos;

		// FENE
		LR_vector<number> rback = r + q->pos_back - p->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;
		energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);
		LR_vector<number> force = rback * (-(FENE_EPS * rbackr0 / (FENE_DELTA2
				- SQR(rbackr0))) / rbackmod);

		p->force -= force;
		q->force += force;

		LR_vector<number> torquep = -p->pos_back.cross(force);
		LR_vector<number> torqueq = q->pos_back.cross(force);

		// excluded volume

		// BASE-BASE
		force.x = force.y = force.z = (number) 0;
		LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
		energy += _excluded_volume(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2,
				EXCL_RC2);
		torquep -= p->pos_base.cross(force);
		torqueq += q->pos_base.cross(force);

		p->force -= force;
		q->force += force;

		// P-BASE vs. Q-BACK
		rcenter = r + q->pos_back - p->pos_base;
		energy += _excluded_volume(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3,
				EXCL_RC3);
		torquep -= p->pos_base.cross(force);
		torqueq += q->pos_back.cross(force);

		p->force -= force;
		q->force += force;

		// P-BACK vs. Q-BASE
		rcenter = r + q->pos_base - p->pos_back;
		energy += _excluded_volume(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4,
				EXCL_RC4);
		torquep -= p->pos_back.cross(force);
		torqueq += q->pos_base.cross(force);

		p->force -= force;
		q->force += force;

		// STACKING
		LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		//number t4 = LRACOS( a3 * b3);
		//number t5 = LRACOS( a3 * rstackdir);
		//number t6 = LRACOS(-b3 * rstackdir);
		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		// functions and their derivatives needed for energies and forces
		//number f1 = this->_interaction.f1(rstackmod, STCK_F1);
		number f1     = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
		//number f4t4   = this->_interaction.f4(t4, STCK_F4_THETA4);
		//number f4t5   = this->_interaction.f4(PI - t5, STCK_F4_THETA5);
		//number f4t6   = this->_interaction.f4(t6, STCK_F4_THETA6);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5 = this->_interaction.query_mesh(-cost5,
				this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6 = this->_interaction.query_mesh(cost6,
				this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D = this->_interaction.f1D(rstackmod, STCK_F1,q->type, p->type);
		//number f4t4Dsin = this->_interaction.f4Dsin(t4, STCK_F4_THETA4);
		//number f4t5Dsin = this->_interaction.f4Dsin(PI - t5, STCK_F4_THETA5);
		//number f4t6Dsin = this->_interaction.f4Dsin(t6, STCK_F4_THETA6);
		number f4t4Dsin = -this->_interaction.query_meshD(cost4,
				this->_interaction.mesh_f4[STCK_F4_THETA4]);
		number f4t5Dsin = -this->_interaction.query_meshD(-cost5,
				this->_interaction.mesh_f4[STCK_F4_THETA5]);
		number f4t6Dsin = -this->_interaction.query_meshD(cost6,
				this->_interaction.mesh_f4[STCK_F4_THETA6]);
		number f5phi1D = this->_interaction.f5D(cosphi1, STCK_F5_PHI1);
		number f5phi2D = this->_interaction.f5D(cosphi2, STCK_F5_PHI2);

		number stck_energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
		energy += stck_energy;

		// RADIAL
		force = -rstackdir * f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

		// THETA 5
		//force += -(a3 - rstackdir * cos(t5)) / rstackmod * f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;
		force += -(a3 - rstackdir * cost5) / rstackmod * f1 * f4t4 * f4t5Dsin
				* f4t6 * f5phi1 * f5phi2;

		// THETA 6
		//force += -(b3 + rstackdir * cos(t6)) / rstackmod * f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;
		force += -(b3 + rstackdir * cost6) / rstackmod * f1 * f4t4 * f4t5
				* f4t6Dsin * f5phi1 * f5phi2;

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		number gamma = POS_STACK - POS_BACK;
		number rbackmodcub = rbackmod * rbackmod * rbackmod;

		number ra2 = rstackdir * a2;
		number ra1 = rstackdir * a1;
		number rb1 = rstackdir * b1;
		number a2b1 = a2 * b1;

		number parentesi = rstackmod * ra2 - a2b1 * gamma;

		number dcosphi1dr = (SQR(rstackmod) * ra2 - ra2 * SQR(rbackmod)
				- rstackmod * (a2b1 + ra2 * (-ra1 + rb1)) * gamma + a2b1
				* (-ra1 + rb1) * SQR(gamma)) / rbackmodcub;
		number dcosphi1dra1 = rstackmod * gamma * parentesi / rbackmodcub;
		number dcosphi1dra2 = -rstackmod / rbackmod;
		number dcosphi1drb1 = -rstackmod * gamma * parentesi / rbackmodcub;
		number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackmodcub;
		number dcosphi1da2b1 = gamma / rbackmod;

		// this force part has a minus because of the definition of cos(phi1)
		// which is not consistent with the derivatives above (i.e. all the
		// derivatives should have a minus sign in front and the force below shouldn't)
		number force_part_phi1 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1D * f5phi2;

		force += -(rstackdir * dcosphi1dr + ((a2 - rstackdir * ra2)
				* dcosphi1dra2 + (a1 - rstackdir * ra1) * dcosphi1dra1 + (b1
				- rstackdir * rb1) * dcosphi1drb1) / rstackmod)
				* force_part_phi1;

		// COS PHI 2
		// here particle p -> b, particle q -> a
		ra2 = rstackdir * b2;
		ra1 = rstackdir * b1;
		rb1 = rstackdir * a1;
		a2b1 = b2 * a1;

		parentesi = rstackmod * ra2 + a2b1 * gamma;

		number dcosphi2dr = (parentesi * (rstackmod + (rb1 - ra1) * gamma)
				- ra2 * SQR(rbackmod)) / rbackmodcub;
		number dcosphi2dra1 = -rstackmod * gamma * (rstackmod * ra2 + a2b1
				* gamma) / rbackmodcub;
		number dcosphi2dra2 = -rstackmod / rbackmod;
		number dcosphi2drb1 = rstackmod * gamma * parentesi / rbackmodcub;
		number dcosphi2da1b1 = -SQR(gamma) * parentesi / rbackmodcub;
		number dcosphi2da2b1 = -gamma / rbackmod;

		// this force part has a minus because of the definition of cos(phi2) which is not consistent with the derivatives
		// above (i.e. all the derivatives should have a minus sign in front and the force below shouldn't)
		number force_part_phi2 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2D;

		force += -force_part_phi2 * (rstackdir * dcosphi2dr + ((b2 - rstackdir
				* ra2) * dcosphi2dra2 + (b1 - rstackdir * ra1) * dcosphi2dra1
				+ (a1 - rstackdir * rb1) * dcosphi2drb1) / rstackmod);

		// Add the contribution from all the forces to the stored particles' forces
		p->force -= force;
		q->force += force;

		torquep -= p->pos_stack.cross(force);
		torqueq += q->pos_stack.cross(force);

		// handle the part on theta4
		LR_vector<number> t4dir = b3.cross(a3);
		number torquemod = f1 * f4t4Dsin * f4t5 * f4t6 * f5phi1 * f5phi2;

		torquep -= t4dir * torquemod;
		torqueq += t4dir * torquemod;

		// handle the part on theta5
		LR_vector<number> t5dir = rstackdir.cross(a3);
		torquemod = -f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2;

		torquep -= t5dir * torquemod;

		// handle the part on theta6
		LR_vector<number> t6dir = rstackdir.cross(b3);
		torquemod = f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2;

		torqueq += t6dir * torquemod;

		// PHI 1
		torquep += rstackdir.cross(a2) * force_part_phi1 * dcosphi1dra2
				+ rstackdir.cross(a1) * force_part_phi1 * dcosphi1dra1;
		torqueq += rstackdir.cross(b1) * force_part_phi1 * dcosphi1drb1;

		LR_vector<number> puretorque = a2.cross(b1) * force_part_phi1
				* dcosphi1da2b1 + a1.cross(b1) * force_part_phi1
				* dcosphi1da1b1;

		torquep -= puretorque;
		torqueq += puretorque;

		// PHI 2
		torquep += rstackdir.cross(a1) * force_part_phi2 * dcosphi2drb1;
		torqueq += rstackdir.cross(b2) * force_part_phi2 * dcosphi2dra2
				+ rstackdir.cross(b1) * force_part_phi2 * dcosphi2dra1;

		puretorque = a1.cross(b2) * force_part_phi2 * dcosphi2da2b1 + a1.cross(
				b1) * force_part_phi2 * dcosphi2da1b1;

		torquep -= puretorque;
		torqueq += puretorque;

		// we need torques in the reference system of the particle
		p->torque += p->orientationT * torquep;
		q->torque += q->orientationT * torqueq;
	}

	return energy;
}

template<typename number>
inline number ProcessData_Backend<number>::_excluded_volume(const LR_vector<
		number> &r, LR_vector<number> &force, number sigma, number rstar,
		number b, number rc) {
	number rmod = r.module();
	number energy = 0;

	force.x = force.y = force.z = 0;
	if (rmod < rc) {
		if (rmod > rstar) {
			number rrc = rmod - rc;
			energy = EXCL_EPS * b * SQR(rrc);
			force = -r * 2 * EXCL_EPS * b * rrc / rmod;
		} else {
			number lj_part = SQR(sigma / rmod) * SQR(sigma / rmod)
					* SQR(sigma / rmod);
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			force = -r * 24 * EXCL_EPS * (lj_part - 2 * SQR(lj_part))
					/ SQR(rmod);
		}
	}

	return energy;
}

template<typename number>
inline number ProcessData_Backend<number>::_particle_particle_interaction(
		Particle<number> *p, Particle<number> *q) {
	// true if p and q are Watson-Crick pairs
	bool is_pair = (q->btype + p->btype == 3);

	LR_vector<number> r = q->pos.minimum_image(p->pos, this->_box_side);

	number energy = 0;

	// excluded volume

	// BASE-BASE
	LR_vector<number> force(0, 0, 0);
	LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
	energy += _excluded_volume(rcenter, force, EXCL_S2, EXCL_R2, EXCL_B2,
			EXCL_RC2);
	LR_vector<number> torquep = -p->pos_base.cross(force);
	LR_vector<number> torqueq = q->pos_base.cross(force);

	p->force -= force;
	q->force += force;

	// P-BASE vs. Q-BACK
	rcenter = r + q->pos_back - p->pos_base;
	energy += _excluded_volume(rcenter, force, EXCL_S3, EXCL_R3, EXCL_B3,
			EXCL_RC3);
	torquep += -p->pos_base.cross(force);
	torqueq += q->pos_back.cross(force);

	p->force -= force;
	q->force += force;

	// P-BACK vs. Q-BASE
	rcenter = r + q->pos_base - p->pos_back;
	energy += _excluded_volume(rcenter, force, EXCL_S4, EXCL_R4, EXCL_B4,
			EXCL_RC4);
	torquep += -p->pos_back.cross(force);
	torqueq += q->pos_base.cross(force);

	p->force -= force;
	q->force += force;

	// BACK-BACK
	rcenter = r + q->pos_back - p->pos_back;
	energy += _excluded_volume(rcenter, force, EXCL_S1, EXCL_R1, EXCL_B1,
			EXCL_RC1);
	torquep += -p->pos_back.cross(force);
	torqueq += q->pos_back.cross(force);

	p->force -= force;
	q->force += force;

	// HYDROGEN BONDING
	LR_vector<number> rhydro = r + q->pos_base - p->pos_base;
	number rhydromod = rhydro.module();
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
		LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		/*number t1 = LRACOS (-a1 * b1);
		 number t2 = LRACOS (-b1 * rhydrodir);
		 number t3 = LRACOS ( a1 * rhydrodir);
		 number t4 = LRACOS ( a3 * b3);
		 number t7 = LRACOS (-b3 * rhydrodir);
		 number t8 = LRACOS ( a3 * rhydrodir);
		 */
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 = a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

		// functions called at their relevant arguments
		//number f1 = this->_interaction.f1(rhydromod, HYDR_F1);
		number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
		/*
		 number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
		 number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
		 number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[HYDR_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh(cost2,
				this->_interaction.mesh_f4[HYDR_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh(cost3,
				this->_interaction.mesh_f4[HYDR_F4_THETA3]);

		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[HYDR_F4_THETA4]);
		number f4t7 = this->_interaction.query_mesh(cost7,
				this->_interaction.mesh_f4[HYDR_F4_THETA7]);
		number f4t8 = this->_interaction.query_mesh(cost8,
				this->_interaction.mesh_f4[HYDR_F4_THETA8]);
		/*number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
		 number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
		 number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);*/

		number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += hb_energy;
		this->_U_hydr += hb_energy;

		// makes sense, since the above functions may return 0. exactly
		if (hb_energy != 0.) {
			// derivatives called at the relevant arguments
			number f1D = this->_interaction.f1D(rhydromod, HYDR_F1,q->type, p->type);
			/*number f4t1Dsin = -this->_interaction.f4Dsin(t1, HYDR_F4_THETA1);
			 number f4t2Dsin = -this->_interaction.f4Dsin(t2, HYDR_F4_THETA2);
			 number f4t3Dsin =  this->_interaction.f4Dsin(t3, HYDR_F4_THETA3);*/
			number f4t1Dsin = this->_interaction.query_meshD(cost1,
					this->_interaction.mesh_f4[HYDR_F4_THETA1]);
			number f4t2Dsin = this->_interaction.query_meshD(cost2,
					this->_interaction.mesh_f4[HYDR_F4_THETA2]);
			number f4t3Dsin = -this->_interaction.query_meshD(cost3,
					this->_interaction.mesh_f4[HYDR_F4_THETA3]);

			number f4t4Dsin = -this->_interaction.query_meshD(cost4,
					this->_interaction.mesh_f4[HYDR_F4_THETA4]);
			number f4t7Dsin = this->_interaction.query_meshD(cost7,
					this->_interaction.mesh_f4[HYDR_F4_THETA7]);
			number f4t8Dsin = -this->_interaction.query_meshD(cost8,
					this->_interaction.mesh_f4[HYDR_F4_THETA8]);
			/*number f4t4Dsin =  this->_interaction.f4Dsin(t4, HYDR_F4_THETA4);
			 number f4t7Dsin = -this->_interaction.f4Dsin(t7, HYDR_F4_THETA7);
			 number f4t8Dsin =  this->_interaction.f4Dsin(t8, HYDR_F4_THETA8);*/

			//f4t1 = 1.;
			//f4t2 = 1.;
			//f4t3 = 1.;
			//f4t4 = 1.;
			//f4t7 = 1.;
			//f4t8 = 1.;

			//f4t1Dsin = 0.;
			//f4t2Dsin = 0.;
			//f4t3Dsin = 0.;
			//f4t4Dsin = 0.;
			//f4t7Dsin = 0.;
			//f4t8Dsin = 0.;

			// RADIAL PART
			force = -rhydrodir * f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> dir = a3.cross(b3);
			number torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7
					* f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = -f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			//force += (b1 + rhydrodir * cos(t2)) / rhydromod * fact;
			force += (b1 + rhydrodir * cost2) / rhydromod * fact;
			dir = rhydrodir.cross(b1);
			//torquemod = - f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += -dir * fact;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			//force += (a1 - rhydrodir * cos(t3)) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
			force += (a1 - rhydrodir * cost3) / rhydromod * f1 * f4t1 * f4t2
					* f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector<number> t3dir = rhydrodir.cross(a1);
			torquemod = -f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			//force += (b3 + rhydrodir * cos(t7)) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;
			force += (b3 + rhydrodir * cost7) / rhydromod * f1 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rhydrodir.cross(b3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			//force +=  (a3 - rhydrodir * cos(t8)) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
			force += (a3 - rhydrodir * cost8) / rhydromod * f1 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rhydrodir.cross(a3);
			torquemod = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->pos_base.cross(force);
			torqueq += q->pos_base.cross(force);
		}
	}
	// END OF HYDROGEN BONDING

	// CROSS STACKING
	LR_vector<number> rcstack = rhydro;
	//LR_vector<number> rcstack = r + q->pos_base - p->pos_base;
	number rcstackmod = rhydromod;
	//number rcstackmod = rcstack.module();
	if (CRST_RCLOW < rcstackmod && rcstackmod < CRST_RCHIGH) {
		LR_vector<number> rcstackdir = rcstack / rcstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		/*
		 number t1 = LRACOS (-a1 * b1);
		 number t2 = LRACOS (-b1 * rcstackdir);
		 number t4 = LRACOS ( a3 * b3);
		 number t3 = LRACOS ( a1 * rcstackdir);
		 number t7 = LRACOS (-rcstackdir * b3);
		 number t8 = LRACOS ( rcstackdir * a3);
		 */
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rcstackdir;
		number cost3 = a1 * rcstackdir;
		number cost4 = a3 * b3;
		number cost7 = -b3 * rcstackdir;
		number cost8 = a3 * rcstackdir;

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rcstackmod, CRST_F2);
		/*
		 number f4t1 = this->_interaction.f4(t1, CRST_F4_THETA1);
		 number f4t2 = this->_interaction.f4(t2, CRST_F4_THETA2);
		 number f4t3 = this->_interaction.f4(t3, CRST_F4_THETA3);
		 number f4t4 = this->_interaction.f4(t4, CRST_F4_THETA4) + this->_interaction.f4(PI - t4, CRST_F4_THETA4);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[CRST_F4_THETA1]);
		number f4t2 = this->_interaction.query_mesh(cost2,
				this->_interaction.mesh_f4[CRST_F4_THETA2]);
		number f4t3 = this->_interaction.query_mesh(cost3,
				this->_interaction.mesh_f4[CRST_F4_THETA3]);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[CRST_F4_THETA4])
				+ this->_interaction.query_mesh(-cost4,
						this->_interaction.mesh_f4[CRST_F4_THETA4]);
		/*
		 number f4t7 = this->_interaction.f4(t7, CRST_F4_THETA7) + this->_interaction.f4(PI - t7, CRST_F4_THETA7);
		 number f4t8 = this->_interaction.f4(t8, CRST_F4_THETA8) + this->_interaction.f4(PI - t8, CRST_F4_THETA8);
		 */
		number f4t7 = this->_interaction.query_mesh(cost7,
				this->_interaction.mesh_f4[CRST_F4_THETA7])
				+ this->_interaction.query_mesh(-cost7,
						this->_interaction.mesh_f4[CRST_F4_THETA7]);
		;
		number f4t8 = this->_interaction.query_mesh(cost8,
				this->_interaction.mesh_f4[CRST_F4_THETA8])
				+ this->_interaction.query_mesh(-cost8,
						this->_interaction.mesh_f4[CRST_F4_THETA8]);
		;

		number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
		energy += cstk_energy;

		// makes sense since the above f? can return exacly 0.
		if (cstk_energy != 0.) {
			// derivatives called at the relevant arguments
			number f2D = this->_interaction.f2D(rcstackmod, CRST_F2);
			/*number f4t1Dsin = -this->_interaction.f4Dsin(t1, CRST_F4_THETA1);
			 number f4t2Dsin = -this->_interaction.f4Dsin(t2, CRST_F4_THETA2);
			 number f4t3Dsin =  this->_interaction.f4Dsin(t3, CRST_F4_THETA3);
			 number f4t4Dsin =  this->_interaction.f4Dsin(t4, CRST_F4_THETA4) - this->_interaction.f4Dsin(PI - t4, CRST_F4_THETA4);*/
			number f4t1Dsin = this->_interaction.query_meshD(cost1,
					this->_interaction.mesh_f4[CRST_F4_THETA1]);
			number f4t2Dsin = this->_interaction.query_meshD(cost2,
					this->_interaction.mesh_f4[CRST_F4_THETA2]);
			number f4t3Dsin = -this->_interaction.query_meshD(cost3,
					this->_interaction.mesh_f4[CRST_F4_THETA3]);
			number f4t4Dsin = -this->_interaction.query_meshD(cost4,
					this->_interaction.mesh_f4[CRST_F4_THETA4])
					+ this->_interaction.query_meshD(-cost4,
							this->_interaction.mesh_f4[CRST_F4_THETA4]);
			/*
			 number f4t7Dsin = -this->_interaction.f4Dsin(t7, CRST_F4_THETA7) + this->_interaction.f4Dsin(PI - t7, CRST_F4_THETA7);
			 number f4t8Dsin =  this->_interaction.f4Dsin(t8, CRST_F4_THETA8) - this->_interaction.f4Dsin(PI - t8, CRST_F4_THETA8);
			 */
			number f4t7Dsin = this->_interaction.query_meshD(cost7,
					this->_interaction.mesh_f4[CRST_F4_THETA7])
					- this->_interaction.query_meshD(-cost7,
							this->_interaction.mesh_f4[CRST_F4_THETA7]);
			number f4t8Dsin = -this->_interaction.query_meshD(cost8,
					this->_interaction.mesh_f4[CRST_F4_THETA8])
					+ this->_interaction.query_meshD(-cost8,
							this->_interaction.mesh_f4[CRST_F4_THETA8]);

			//f4t1 = 1.;
			//f4t2 = 1.;
			//f4t3 = 1.;
			//f4t4 = 1.;
			//f4t7 = 1.;
			//f4t8 = 1.;

			//f4t1Dsin = 0.;
			//f4t2Dsin = 0.;
			//f4t3Dsin = 0.;
			//f4t4Dsin = 0.;
			//f4t7Dsin = 0.;
			//f4t8Dsin = 0.;

			/*
			 FILE * yo = fopen("yo.dat", "w");
			 for (int kk = 0; kk<1000.; kk++)
			 {
			 fprintf(yo,"%lf %lf %lf\n", 0.002 * (kk - 500) * 2 * PI, this->_interaction.f4(0.002 * (kk-500) * 2 * PI, CRST_F4_THETA3), this->_interaction.f4D(0.002 * (kk-500) * 2 * PI, CRST_F4_THETA3));
			 }
			 fclose(yo);
			 exit(-1);*/

			// RADIAL PART
			force = -rcstackdir * f2D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> t1dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7
					* f4t8;

			torquep -= t1dir * torquemod;
			torqueq += t1dir * torquemod;

			// TETA2; t2 = LRACOS (-b1 * rhydrodir);
			//force += (b1 + rcstackdir * cos(t2)) / rcstackmod * f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rcstackdir * cost2) / rcstackmod * f2 * f4t1
					* f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			LR_vector<number> t2dir = rcstackdir.cross(b1);
			torquemod = -f2 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;

			torqueq += t2dir * torquemod;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			//force += (a1 - rcstackdir * cos(t3)) / rcstackmod * f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;
			force += (a1 - rcstackdir * cost3) / rcstackmod * f2 * f4t1 * f4t2
					* f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector<number> t3dir = rcstackdir.cross(a1);
			torquemod = -f2 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> t4dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8;

			torquep -= t4dir * torquemod;
			torqueq += t4dir * torquemod;

			// THETA7; t7 = LRACOS (-rcsrackir * b3);
			//force += (b3 + rcstackdir * cos(t7)) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;
			force += (b3 + rcstackdir * cost7) / rcstackmod * f2 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rcstackdir.cross(b3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			// force +=  (a3 - rcstackdir * cos(t8)) / rcstackmod * f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;
			force += (a3 - rcstackdir * cost8) / rcstackmod * f2 * f4t1 * f4t2
					* f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rcstackdir.cross(a3);
			torquemod = -f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for CRST
			p->force -= force;
			q->force += force;

			torquep -= p->pos_base.cross(force);
			torqueq += q->pos_base.cross(force);
		}
	}

	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
	number rstackmod = rstack.module();
	if (CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		/*
		 number t1 = LRACOS (-a1 * b1);
		 number t4 = LRACOS ( a3 * b3);
		 number t5 = LRACOS ( a3 * rstackdir);
		 number t6 = LRACOS (-b3 * rstackdir);
		 */
		number cost1 = -a1 * b1;
		number cost4 = a3 * b3;
		number cost5 = a3 * rstackdir;
		number cost6 = -b3 * rstackdir;
		LR_vector<number> rbackbone = r + q->pos_back - p->pos_back;
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rstackmod, CXST_F2);
		/*
		 number f4t1 = this->_interaction.f4(t1, CXST_F4_THETA1) + this->_interaction.f4(2 * PI - t1, CXST_F4_THETA1);
		 number f4t4 = this->_interaction.f4(t4, CXST_F4_THETA4);
		 number f4t5 = this->_interaction.f4(t5, CXST_F4_THETA5) + this->_interaction.f4(PI - t5, CXST_F4_THETA5);
		 number f4t6 = this->_interaction.f4(t6, CXST_F4_THETA6) + this->_interaction.f4(PI - t6, CXST_F4_THETA6);
		 */
		number f4t1 = this->_interaction.query_mesh(cost1,
				this->_interaction.mesh_f4[CXST_F4_THETA1]);
		number f4t4 = this->_interaction.query_mesh(cost4,
				this->_interaction.mesh_f4[CXST_F4_THETA4]);
		number f4t5 = this->_interaction.query_mesh(cost5,
				this->_interaction.mesh_f4[CXST_F4_THETA5])
				+ this->_interaction.query_mesh(-cost5,
						this->_interaction.mesh_f4[CXST_F4_THETA5]);
		number f4t6 = this->_interaction.query_mesh(cost6,
				this->_interaction.mesh_f4[CXST_F4_THETA6])
				+ this->_interaction.query_mesh(-cost6,
						this->_interaction.mesh_f4[CXST_F4_THETA6]);
		number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

		number cxst_energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
		energy += cxst_energy;

		// again, makes sense (see same line for HB)
		if (cxst_energy != 0.) {
			// derivatives called at the relevant arguments
			number f2D = this->_interaction.f2D(rstackmod, CXST_F2);
			/*
			 number f4t1Dsin = -this->_interaction.f4Dsin(t1, CXST_F4_THETA1) + this->_interaction.f4Dsin(2 * PI - t1, CXST_F4_THETA1);
			 number f4t4Dsin =  this->_interaction.f4Dsin(t4, CXST_F4_THETA4);
			 number f4t5Dsin =  this->_interaction.f4Dsin(t5, CXST_F4_THETA5) - this->_interaction.f4Dsin(PI - t5, CXST_F4_THETA5);
			 number f4t6Dsin = -this->_interaction.f4Dsin(t6, CXST_F4_THETA6) + this->_interaction.f4Dsin(PI - t6, CXST_F4_THETA6);*/
			number f4t1Dsin = this->_interaction.query_meshD(cost1,
					this->_interaction.mesh_f4[CXST_F4_THETA1]);
			number f4t4Dsin = -this->_interaction.query_meshD(cost4,
					this->_interaction.mesh_f4[CXST_F4_THETA4]);
			number f4t5Dsin = -this->_interaction.query_meshD(cost5,
					this->_interaction.mesh_f4[CXST_F4_THETA5])
					+ this->_interaction.query_meshD(-cost5,
							this->_interaction.mesh_f4[CXST_F4_THETA5]);
			number f4t6Dsin = this->_interaction.query_meshD(cost6,
					this->_interaction.mesh_f4[CXST_F4_THETA6])
					- this->_interaction.query_meshD(-cost6,
							this->_interaction.mesh_f4[CXST_F4_THETA6]);

			number f5Dcosphi3 = this->_interaction.f5D(cosphi3, CXST_F5_PHI3);

			//f4t1 = 1.;
			//f4t4 = 1.;
			//f4t5 = 1.;
			//f4t6 = 1.;
			//f5cosphi3 = 1.;

			//f4t1Dsin = 0.;
			//f4t4Dsin = 0.;
			//f4t5Dsin = 0.;
			//f4t6Dsin = 0.;
			//f5Dcosphi3 = 0.;

			/*FILE * yo = fopen("yo.dat", "w");
			 for (int kk = 0; kk<10000; kk++)
			 {
			 fprintf(yo,"%lf %lf %lf\n", 0.0002 * (kk - 5000) * 2 * PI, this->_interaction.f5(0.0002 * (kk-5000) * 2 * PI, CXST_F5_PHI3), this->_interaction.f5D(0.0002 * (kk - 5000) * 2 * PI, CXST_F5_PHI3));
			 }
			 fclose(yo);
			 exit(-1);*/

			// RADIAL PART
			force = -rstackdir * f2D * f4t1 * f4t4 * f4t5 * f4t6
					* SQR(f5cosphi3);

			// THETA1; t1 = LRACOS (-a1 * b1);
			LR_vector<number> dir = a1.cross(b1);
			number torquemod = -f2 * f4t1Dsin * f4t4 * f4t5 * f4t6
					* SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// TETA4; t4 = LRACOS (a3 * b3);
			dir = a3.cross(b3);
			torquemod = -f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * SQR(f5cosphi3);

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA5; t5 = LRACOS ( a3 * rstackdir);
			number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * SQR(f5cosphi3);
			//force += fact * (a3 - rstackdir * cos(t5)) / rstackmod;
			force += fact * (a3 - rstackdir * cost5) / rstackmod;
			dir = rstackdir.cross(a3);
			torquep -= dir * fact;

			// THETA6; t6 = LRACOS (-b3 * rstackdir);
			fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * SQR(f5cosphi3);
			//force += (b3 + rstackdir * cos(t6)) / rstackmod * fact;
			force += (b3 + rstackdir * cost6) / rstackmod * fact;
			dir = rstackdir.cross(b3);

			torqueq += -dir * fact;

			// Cosphi3 (qui son dolori...) (meno male che cosphi4 = cosphi3)
			// Definition used:
			// cosphi3 = gamma * (ra3 * a2b1 - ra2 * a3b1)/ rbackmod
			number gamma = POS_STACK - POS_BACK;
			number gammacub = gamma * gamma * gamma;
			number rbackmodcub = rbackmod * rbackmod * rbackmod;
			//number a1b1 = a1 * b1;
			number a2b1 = a2 * b1;
			number a3b1 = a3 * b1;
			number ra1 = rstackdir * a1;
			number ra2 = rstackdir * a2;
			number ra3 = rstackdir * a3;
			number rb1 = rstackdir * b1;

			number parentesi = (ra3 * a2b1 - ra2 * a3b1);

			number dcdr = -gamma * parentesi
					* (gamma * (ra1 - rb1) + rstackmod) / rbackmodcub;
			number dcda1b1 = gammacub * parentesi / rbackmodcub;
			number dcda2b1 = gamma * ra3 / rbackmod;
			number dcda3b1 = -gamma * ra2 / rbackmod;
			number dcdra1 = -SQR(gamma) * parentesi * rstackmod / rbackmodcub;
			number dcdra2 = -gamma * a3b1 / rbackmod;
			number dcdra3 = gamma * a2b1 / rbackmod;
			number dcdrb1 = SQR(gamma) * parentesi * rstackmod / rbackmodcub;

			number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * 2 * f5cosphi3
					* f5Dcosphi3;

			force += -force_c * (rstackdir * dcdr + ((a1 - rstackdir * ra1)
					* dcdra1 + (a2 - rstackdir * ra2) * dcdra2 + (a3
					- rstackdir * ra3) * dcdra3 + (b1 - rstackdir * rb1)
					* dcdrb1) / rstackmod);

			torquep += force_c * (rstackdir.cross(a1) * dcdra1
					+ rstackdir.cross(a2) * dcdra2 + rstackdir.cross(a3)
					* dcdra3);
			torqueq += force_c * (rstackdir.cross(b1) * dcdrb1);

			LR_vector<number> puretorque = force_c * (a1.cross(b1) * dcda1b1
					+ a2.cross(b1) * dcda2b1 + a3.cross(b1) * dcda3b1);
			torquep -= puretorque;
			torqueq += puretorque;

			// final update of forces and torques for CXST
			p->force -= force;
			q->force += force;

			torquep -= p->pos_stack.cross(force);
			torqueq += q->pos_stack.cross(force);
		}
	}

	// Final update of the torques
	// total torques
	p->torque += p->orientationT * torquep;
	q->torque += q->orientationT * torqueq;

	return energy;
}

template<typename number>
void ProcessData_Backend<number>::_compute_forces() {
	int neigh;
	Particle<number> *p;

	this->_U = this->_U_hydr = (number) 0;
	for (int i = 0; i < this->_N; i++) {
		p = &this->_particles[i];
		this->_U += _particle_particle_bonded_interaction(p);

		p->prepare_list();
		neigh = p->next_neighbour();
		while (neigh != P_VIRTUAL) {
			this->_U += _particle_particle_interaction(p,
					&this->_particles[neigh]);
			neigh = p->next_neighbour();
		}
	}
}

/*
template<typename number>
void ProcessData_Backend<number>::_second_step() {
	this->_K = (number) 0;
	for (int i = 0; i < this->_N; i++) {
		this->_particles[i].vel += this->_particles[i].force * this->_dt
				* (number) 0.5;
		this->_particles[i].L += this->_particles[i].torque * this->_dt
				* (number) 0.5;

		this->_K += (this->_particles[i].vel.norm()
				+ this->_particles[i].L.norm()) * (number) 0.5;
	}
}

template<typename number>
void ProcessData_Backend<number>::_activate_john_thermostat() {
	for (int i = 0; i < this->_N; i++) {
		if (drand48() < this->_pt) {
			this->_particles[i].vel = LR_vector<number> (
					Utils::gaussian<number>(), Utils::gaussian<number>(),
					Utils::gaussian<number>()) * this->_rescale_factor;
		}
		if (drand48() < this->_pr) {
			this->_particles[i].L = LR_vector<number> (
					Utils::gaussian<number>(), Utils::gaussian<number>(),
					Utils::gaussian<number>()) * this->_rescale_factor;
		}
	}
}

template<typename number>
void ProcessData_Backend<number>::_activate_refresh_thermostat() {
	for (int i = 0; i < this->_N; i++) {
		this->_particles[i].vel = LR_vector<number> (Utils::gaussian<number>(),
				Utils::gaussian<number>(), Utils::gaussian<number>())
				* this->_rescale_factor;
		this->_particles[i].L = LR_vector<number> (Utils::gaussian<number>(),
				Utils::gaussian<number>(), Utils::gaussian<number>())
				* this->_rescale_factor;
	}
}

template<typename number>
void ProcessData_Backend<number>::_rescale_velocities() {
	_vcm = LR_vector<number> (0, 0, 0);
	for (int i = 0; i < this->_N; i++) {
		_vcm += this->_particles[i].vel;
	}
	_vcm /= this->_N;

	for (int i = 0; i < this->_N; i++) {
		this->_particles[i].vel -= _vcm;
	}
}

template<typename number>
void ProcessData_Backend<number>::sim_step(llint curr_step) {
	get_time(&this->_timer, 0);

	get_time(&this->_timer, 2);
	_first_step(curr_step);
	get_time(&this->_timer, 3);

	get_time(&this->_timer, 6);
	if (this->_are_lists_old == true) {
		this->_update_lists();
		this->_N_updates++;
	}
	get_time(&this->_timer, 7);

	get_time(&this->_timer, 8);
	_compute_forces();
	_second_step();
	get_time(&this->_timer, 9);

	get_time(&this->_timer, 10);
	if (this->_thermostat != this->THERMOSTAT_NO && (curr_step
			% this->_newtonian_steps == 0)) {
		if (this->_thermostat == this->THERMOSTAT_JOHN)
			_activate_john_thermostat();
		else
			_activate_refresh_thermostat();
	}
	get_time(&this->_timer, 11);

	get_time(&this->_timer, 1);

	process_times(&this->_timer);
}

*/
template<typename number>
inline number ProcessData_Backend<number>::mc_particle_particle_interaction(
		Particle<number> *p, Particle<number> *q) {
	// true if p and q are Watson-Crick pairs
	bool is_pair = (q->type + p->type == 3);

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
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		number t1 = LRACOS (-a1 * b1);
		number t2 = LRACOS (-b1 * rhydrodir);
		number t3 = LRACOS ( a1 * rhydrodir);
		number t4 = LRACOS ( a3 * b3);
		number t7 = LRACOS (-b3 * rhydrodir);
		number t8 = LRACOS ( a3 * rhydrodir);

		// functions called at their relevant arguments
		//number f1 = this->_interaction.f1(rhydromod, HYDR_F1);
		number f1   = this->_interaction.f1(rhydromod, HYDR_F1, q->type, p->type);
		number f4t1 = this->_interaction.f4(t1, HYDR_F4_THETA1);
		number f4t2 = this->_interaction.f4(t2, HYDR_F4_THETA2);
		number f4t3 = this->_interaction.f4(t3, HYDR_F4_THETA3);
		number f4t4 = this->_interaction.f4(t4, HYDR_F4_THETA4);
		number f4t7 = this->_interaction.f4(t7, HYDR_F4_THETA7);
		number f4t8 = this->_interaction.f4(t8, HYDR_F4_THETA8);

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
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CRST interaction
		number t1, t2, t3, t4, t7, t8;
		t1 = LRACOS (-a1 * b1);
		t2 = LRACOS (-b1 * rcstackdir);
		t3 = LRACOS ( a1 * rcstackdir);
		t4 = LRACOS ( a3 * b3);
		t7 = LRACOS (-rcstackdir * b3);
		t8 = LRACOS ( rcstackdir * a3);

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rcstackmod, CRST_F2);
		number f4t1 = this->_interaction.f4(t1, CRST_F4_THETA1);
		number f4t2 = this->_interaction.f4(t2, CRST_F4_THETA2);
		number f4t3 = this->_interaction.f4(t3, CRST_F4_THETA3);
		number f4t4 = this->_interaction.f4(t4, CRST_F4_THETA4)
				+ this->_interaction.f4(PI - t4, CRST_F4_THETA4);
		number f4t7 = this->_interaction.f4(t7, CRST_F4_THETA7)
				+ this->_interaction.f4(PI - t7, CRST_F4_THETA7);
		number f4t8 = this->_interaction.f4(t8, CRST_F4_THETA8)
				+ this->_interaction.f4(PI - t8, CRST_F4_THETA8);

		energy += f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	}

	// COAXIAL STACKING
	LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
	number rstackmod = rstack.module();
	if (CXST_RCLOW < rstackmod && rstackmod < CXST_RCHIGH) {
		LR_vector<number> rstackdir = rstack / rstackmod;

		// particle axes according to Allen's paper
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;

		// angles involved in the CXST interaction
		number t1 = LRACOS (-a1 * b1);
		number t4 = LRACOS ( a3 * b3);
		number t5 = LRACOS ( a3 * rstackdir);
		number t6 = LRACOS (-b3 * rstackdir);
		LR_vector<number> rbackbone = r + q->pos_back - p->pos_back;
		number rbackmod = rbackbone.module();
		LR_vector<number> rbackbonedir = rbackbone / rbackmod;
		number cosphi3 = rstackdir * (rbackbonedir.cross(a1));

		// functions called at their relevant arguments
		number f2 = this->_interaction.f2(rstackmod, CXST_F2);
		number f4t1 = this->_interaction.f4(t1, CXST_F4_THETA1)
				+ this->_interaction.f4(2 * PI - t1, CXST_F4_THETA1);
		number f4t4 = this->_interaction.f4(t4, CXST_F4_THETA4);
		number f4t5 = this->_interaction.f4(t5, CXST_F4_THETA5)
				+ this->_interaction.f4(PI - t5, CXST_F4_THETA5);
		number f4t6 = this->_interaction.f4(t6, CXST_F4_THETA6)
				+ this->_interaction.f4(PI - t6, CXST_F4_THETA6);
		number f5cosphi3 = this->_interaction.f5(cosphi3, CXST_F5_PHI3);

		energy += f2 * f4t1 * f4t4 * f4t5 * f4t6 * SQR(f5cosphi3);
	}

	return energy;
}





template<typename number>
inline number ProcessData_Backend<number>::mc_particle_particle_bonded_interaction(
		Particle<number> *p) {
	number energy = (number) 0;

	if (p->n3 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n3];
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = q->pos - p->pos;

		// FENE
		LR_vector<number> rback = r + q->pos_back - p->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;
		energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

		// excluded volume

		// BASE-BASE
		LR_vector<number> rcenter = r + q->pos_base - p->pos_base;
		energy
				+= _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2,
						EXCL_RC2);

		// P-BASE vs. Q-BACK
		rcenter = r + q->pos_back - p->pos_base;
		energy
				+= _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3,
						EXCL_RC3);

		// P-BACK vs. Q-BASE
		rcenter = r + q->pos_base - p->pos_back;
		energy
				+= _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4,
						EXCL_RC4);

		// STACKING
		LR_vector<number> rstack = r + q->pos_stack - p->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		number t4 = LRACOS( a3 * b3);
		number t5 = LRACOS( a3 * rstackdir);
		number t6 = LRACOS(-b3 * rstackdir);
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		// functions and their derivatives needed for energies and forces
		//number f1 = this->_interaction.f1(rstackmod, STCK_F1);
		number f1     = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
		number f4t4 = this->_interaction.f4(t4, STCK_F4_THETA4);
		number f4t5 = this->_interaction.f4(PI - t5, STCK_F4_THETA5);
		number f4t6 = this->_interaction.f4(t6, STCK_F4_THETA6);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	}

	if (p->n5 != P_VIRTUAL) {
		Particle<number> *q = &this->_particles[p->n5];
		LR_vector<number> a1 = p->orientationT.v1;
		LR_vector<number> a2 = p->orientationT.v2;
		LR_vector<number> a3 = p->orientationT.v3;
		LR_vector<number> b1 = q->orientationT.v1;
		LR_vector<number> b2 = q->orientationT.v2;
		LR_vector<number> b3 = q->orientationT.v3;
		LR_vector<number> r = p->pos - q->pos;

		// FENE
		LR_vector<number> rback = r + p->pos_back - q->pos_back;
		number rbackmod = rback.module();
		number rbackr0 = rbackmod - FENE_R0;
		energy += -FENE_EPS * 0.5 * log(1 - SQR(rbackr0) / FENE_DELTA2);

		// excluded volume

		// BASE-BASE
		LR_vector<number> rcenter = r + p->pos_base - q->pos_base;
		energy
				+= _excluded_volume(rcenter, EXCL_S2, EXCL_R2, EXCL_B2,
						EXCL_RC2);

		// P-BASE vs. Q-BACK
		rcenter = r + p->pos_back - q->pos_base;
		energy
				+= _excluded_volume(rcenter, EXCL_S3, EXCL_R3, EXCL_B3,
						EXCL_RC3);

		// P-BACK vs. Q-BASE
		rcenter = r + p->pos_base - q->pos_back;
		energy
				+= _excluded_volume(rcenter, EXCL_S4, EXCL_R4, EXCL_B4,
						EXCL_RC4);

		// STACKING
		LR_vector<number> rstack = r + p->pos_stack - q->pos_stack;
		number rstackmod = rstack.module();
		LR_vector<number> rstackdir = rstack / rstackmod;

		number t4 = LRACOS( a3 * b3);
		number t5 = LRACOS( a3 * rstackdir);
		number t6 = LRACOS(-b3 * rstackdir);
		number cosphi1 = a2 * rback / rbackmod;
		number cosphi2 = b2 * rback / rbackmod;

		//number f1 = this->_interaction.f1(rstackmod, STCK_F1);
		number f1     = this->_interaction.f1(rstackmod, STCK_F1, q->type, p->type);
		number f4t4 = this->_interaction.f4(t4, STCK_F4_THETA4);
		number f4t5 = this->_interaction.f4(PI - t5, STCK_F4_THETA5);
		number f4t6 = this->_interaction.f4(t6, STCK_F4_THETA6);
		number f5phi1 = this->_interaction.f5(cosphi1, STCK_F5_PHI1);
		number f5phi2 = this->_interaction.f5(cosphi2, STCK_F5_PHI2);

		energy += f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;
	}

	return energy;
}

template<typename number>
//void ProcessData_Backend<number>::init(ifstream &conf_input) {
void ProcessData_Backend<number>::init(char conf_filename[256]) {
	//SimBackend<number>::init(conf_input);
	SimBackend<number>::init(conf_filename);
}

template <typename number> int ProcessData_Backend<number>::init_first_state(ifstream& conf_input) {
	//init_timer(&(this->_timer),&(this->_timer_msgs_number),&(this->_timer_msgs));
	this->_interaction.init(this->_T);

	this->_read_topology();
	if (this->_external_forces) this->_read_external_forces ();

	int iteration = -1;
	conf_input.seekg(0);
	char line[512];
	conf_input.getline(line, 512);
	sscanf(line, "t = %d",&iteration);
	conf_input.getline(line, 512);
	double box_side;
	sscanf(line, "b = %lf %*f %*f", &box_side);
	this->_box_side = (number) box_side;
	conf_input.getline(line, 512);
	sscanf(line, "E = %*f %*f %*f");

	this->_sqr_verlet_skin = SQR(this->_verlet_skin);
	// we choose rcut as the max of the range interaction of excluded volume between backbones and
	// hydrogen bonding
	number rcutback = 2 * fabs(POS_BACK) + 2 * EXCL_RC1;
	number rcutbase = 2 * fabs(POS_BASE) + 2 * HYDR_RCHIGH;
	this->_rcut = fmax(rcutback, rcutbase);
	this->_sqr_rcut = SQR(this->_rcut);
	this->_sqr_rverlet = SQR(this->_rcut + this->_verlet_skin * 2);

	this->_N_cells_side = (int) floor(this->_box_side / sqrt(this->_sqr_rverlet));
	while(this->_N_cells_side > ceil(pow(2*this->_N, 1/3.)) && this->_N_cells_side > 3) {
		this->_N_cells_side--;
	}

	if(this->_N_cells_side < 3) this->_IO->die("N_cells_side (%d) must be > 2", this->_N_cells_side);

	this->_N_cells = this->_N_cells_side*this->_N_cells_side*this->_N_cells_side;
	this->_heads = new int[this->_N_cells];

	this->_IO->log(this->_IO->LOG_INFO, "N_cells: %d, N_cells_side: %d", this->_N_cells, this->_N_cells_side);

	this->_max_neigh = 0;

	int i = 0;
	while(!conf_input.eof() && i < this->_N) {
		conf_input >> this->_particles[i].pos.x >> this->_particles[i].pos.y >> this->_particles[i].pos.z;
		conf_input >> this->_particles[i].orientation.v1.x >> this->_particles[i].orientation.v1.y >> this->_particles[i].orientation.v1.z;
		conf_input >> this->_particles[i].orientation.v3.x >> this->_particles[i].orientation.v3.y >> this->_particles[i].orientation.v3.z;
		this->_particles[i].orientation.v1.normalize();
		this->_particles[i].orientation.v3.normalize();
		// orthonormalization
		this->_particles[i].orientation.v1 -= this->_particles[i].orientation.v3 * (this->_particles[i].orientation.v1*this->_particles[i].orientation.v3);
		this->_particles[i].orientation.v1.normalize();
		this->_particles[i].orientation.v2 = this->_particles[i].orientation.v3.cross(this->_particles[i].orientation.v1);
		this->_particles[i].orientation.transpone();

		conf_input >> this->_particles[i].vel.x >> this->_particles[i].vel.y >> this->_particles[i].vel.z;
		conf_input >> this->_particles[i].L.x >> this->_particles[i].L.y >> this->_particles[i].L.z;

		this->_particles[i].init(this->_max_neigh);

		i++;
	}

	if(i != this->_N) this->_IO->die("The number of lines found in configuration file (%d) doesn't match the parsed number of particles (%d)", i, this->_N);

	
	this->_IO->log(this->_IO->LOG_INFO, "N: %d",this->_N);

	this->_check_input_sanity();

	return iteration;
	

}

template <typename number> int ProcessData_Backend<number>::init_to_next_state(ifstream& conf_input)
		{

                int iteration = 0;
		char line[512];
		conf_input.getline(line, 512);
		conf_input.getline(line, 512);
		if(!conf_input.good()) 
		  return -1;
		sscanf(line, "t = %d",&iteration);
		conf_input.getline(line, 512);
		cout << line << endl;
		double box_side;
		sscanf(line, "b = %lf %*f %*f", &box_side);
		this->_box_side = (number) box_side;
		conf_input.getline(line, 512);
		sscanf(line, "E = %*f %*f %*f");
                if(!conf_input.good()) 
		  return -1;
               
		

		// if we are performing a CUDA simulation then we don't need CPU verlet lists

		int i = 0;
		while(!conf_input.eof() &&  i < this->_N) {
			conf_input >> this->_particles[i].pos.x >> this->_particles[i].pos.y >> this->_particles[i].pos.z;
			conf_input >> this->_particles[i].orientation.v1.x >> this->_particles[i].orientation.v1.y >> this->_particles[i].orientation.v1.z;
			conf_input >> this->_particles[i].orientation.v3.x >> this->_particles[i].orientation.v3.y >> this->_particles[i].orientation.v3.z;
			this->_particles[i].orientation.v1.normalize();
			this->_particles[i].orientation.v3.normalize();
			// orthonormalization
			this->_particles[i].orientation.v1 -= this->_particles[i].orientation.v3 * (this->_particles[i].orientation.v1*this->_particles[i].orientation.v3);
			this->_particles[i].orientation.v1.normalize();
			this->_particles[i].orientation.v2 = this->_particles[i].orientation.v3.cross(this->_particles[i].orientation.v1);
			this->_particles[i].orientation.transpone();

			conf_input >> this->_particles[i].vel.x >> this->_particles[i].vel.y >> this->_particles[i].vel.z;
			conf_input >> this->_particles[i].L.x >> this->_particles[i].L.y >> this->_particles[i].L.z;

			this->_particles[i].init(this->_max_neigh);

			i++;
		}
		if(i != this->_N)
		{
			return -1;
		}

		cout << "reading parts: " << i << endl;
		//this->_check_input_sanity();
                return iteration;

		}

template class ProcessData_Backend<float> ;
template class ProcessData_Backend<double> ;
