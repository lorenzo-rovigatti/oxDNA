#include "RNAInteraction2.h"

template<typename number>
RNA2Interaction<number>::RNA2Interaction() : RNAInteraction<number>() {

	this->_int_map[DEBYE_HUCKEL] = (number (RNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &RNA2Interaction<number>::_debye_huckel;
	_RNA_HYDR_MIS = 1;
    // log the interaction type
	OX_LOG(Logger::LOG_INFO,"Running modification of oxRNA with additional Debye-Huckel potential");
}



template<typename number>
number RNA2Interaction<number>::pair_interaction_nonbonded(
		BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r,
		bool update_forces) {

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		computed_r = this->_box->min_image(p->pos, q->pos);
		r = &computed_r;
	}

	// if the distance is beyond the cutoff, return 0
        if (r->norm() >= this->_sqr_rcut)
		return (number) 0.f;
	// compute the interaction energy as always ...
	number energy = RNAInteraction<number>::pair_interaction_nonbonded(p, q, r,
			update_forces);

	// ... and then add the debye_huckel energy
	energy += _debye_huckel(p, q, r, update_forces);

	return energy;
}

template<typename number>
void RNA2Interaction<number>::get_settings(input_file &inp) {
	RNAInteraction<number>::get_settings(inp);

	float salt;
	int mismatches;
	float prefactor; // this is the strength of the interaction
	float lambdafactor; //Lambda is _debye_huckel_LAMBDAFACTOR / salt^0.5
	float rh;
	int lambda_T_dependent = 0;
	int half_charged_ends = 0;
	//getInputString(&inp, "topology", this->_topology_filename, 1);

	// read salt concentration from file, or set it to the default value
	if (getInputFloat(&inp, "salt_concentration", &salt, 0) == KEY_FOUND) {
		this->_salt_concentration = (float) salt;
	} else {
		this->_salt_concentration = 1.0f;
	}

	// read lambda-factor (the dh length at I = 1.0), or set it to the default value
	if (getInputFloat(&inp, "dh_lambda", &lambdafactor, 0) == KEY_FOUND) {
			_debye_huckel_lambdafactor = (float) lambdafactor;
	} else {
			_debye_huckel_lambdafactor = 0.3667258;
	}


	if (getInputBoolAsInt(&inp, "lambda_T_dependent", &lambda_T_dependent, 0) == KEY_FOUND  && lambda_T_dependent) {
			OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel with temperature dependent lambda at T = %f",this->_T);
			_debye_huckel_lambdafactor *= sqrt(this->_T / 0.1f) ;
	} 

	// read the prefactor to the potential, or set it to the default value
	if (getInputFloat(&inp, "dh_strength", &prefactor, 0) == KEY_FOUND) {
				_debye_huckel_prefactor = (float) prefactor;
			} else {
				_debye_huckel_prefactor =  0.0858; // 0.510473f;  Q = 0.0858
	}
	// read the cutoff distance (todo - set the default value to something 
	// which makes sense, maybe that depends on lambdafactor and salt concentration
	number lambda = _debye_huckel_lambdafactor * sqrt(this->_T / 0.1f) / sqrt(_salt_concentration);
	if(getInputFloat(&inp, "debye_huckel_rhigh", &rh, 0) == KEY_FOUND) {
	                        _debye_huckel_RHIGH = (float) rh; 
	                } else {
	                        _debye_huckel_RHIGH = 3.0 * lambda;
	}
	// read the mismatch_repulsion flag (not implemented yet)
	if(getInputBoolAsInt(&inp, "mismatch_repulsion", &mismatches, 0) == KEY_FOUND) {
			this->_mismatch_repulsion = (bool) mismatches;
	}
	else {
			this->_mismatch_repulsion = false;
	}
        // whether to halve the charge on the terminus of the strand or not (default true)
        if(getInputBoolAsInt(&inp, "dh_half_charged_ends", &half_charged_ends, 0) == KEY_FOUND) {
        		this-> _debye_huckel_half_charged_ends = (bool) half_charged_ends;	
	}
	else {
			this-> _debye_huckel_half_charged_ends = true;
        }

	//log it 
	OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel at salt_concentration =  %g",this->_salt_concentration);



	if(this->_mismatch_repulsion)
	{
		float temp;
		if(getInputFloat(&inp, "mismatch_repulsion_strength", &temp, 0) == KEY_FOUND) {
					this->_RNA_HYDR_MIS = temp;
		}
		else
		{
			this->_RNA_HYDR_MIS = 1;
		}

	}



}

template<typename number> //initialise the interaction
void RNA2Interaction<number>::init() {

	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	//number rcutback = 2 * fabs(model->RNA_POS_BACK) + model->RNA_EXCL_RC1;
	//this->_T = T;

        //perform the initialisation as in RNAinteraction
	RNAInteraction<number>::init();

	//compute the DH length lambda
	number lambda = _debye_huckel_lambdafactor / sqrt(_salt_concentration);

	//_debye_huckel_RHIGH =  _debye_huckel_cutoff_factor * lambda - 0.1; //2.0 * lambda - 0.1;
	//_debye_huckel_Vrc  = (_debye_huckel_prefactor / (_debye_huckel_cutoff_factor * lambda))  * exp(- _debye_huckel_cutoff_factor ) ;
        _debye_huckel_Vrc = 0;
	_minus_kappa = -1.0/lambda;
	
	number x = _debye_huckel_RHIGH;
	number q = _debye_huckel_prefactor;
	number l = lambda;
	number V = 0.; //We don't want any offset for the potential
	
	// now the smoothening transition at 2 * lambda - 0.1
	_debye_huckel_B = -(exp(-x/l) * q * q * (x + l)*(x+l) )/(4.*x*x*x * l * l * (-q + exp(x/l) * V*x) );
	_debye_huckel_RC = x*(q*x + 3. * q* l - 2.0 * exp(x/l) * V * x * l )/(q * (x+l));


	number debyecut =  2. * sqrt(SQR(this->model->RNA_POS_BACK_a1)	+ SQR(this->model->RNA_POS_BACK_a2) + SQR(this->model->RNA_POS_BACK_a3)) + _debye_huckel_RC;
        //number debyecut = 2. * sqrt(SQR(POS_BACK) ) + _debye_huckel_RC;
	if (debyecut > this->_rcut)
	{
		this->_rcut = debyecut;
		this->_sqr_rcut = debyecut*debyecut;
	}
	// log the parameters of the Debye-Huckel

	OX_LOG(Logger::LOG_INFO,"DEBUGGING: rhigh is %g, Cutoff is %g, RC huckel is %g, B huckel is %g, V is %g, lambda is %g ",_debye_huckel_RHIGH,this->_rcut,_debye_huckel_RC, _debye_huckel_B,_debye_huckel_Vrc,lambda);
	OX_LOG(Logger::LOG_INFO,"DEBUGGING: dh_half_charged_ends = %s", _debye_huckel_half_charged_ends ? "true" : "false");




//
/*
	printf("# B = %g\n",_debye_huckel_B);
	printf("# Rc = %g\n",_debye_huckel_RC);
	printf("# Rh = %g\n",_debye_huckel_RHIGH);
	    for (double x = 0.01; x <= _debye_huckel_RC +0.02;  x += 0.02)
	    {
	    	//energy = _debye_huckel(this->partic)
	    	printf("%g %g \n",x,test_huckel(x));
	    }
	    exit(1);
*/

	if(_mismatch_repulsion)
	{
	 float temp = -1.0f  * _RNA_HYDR_MIS / this->model->RNA_HYDR_EPS;
	 this->F1_EPS[RNA_HYDR_F1][0][0] *= temp;
	 this->F1_SHIFT[RNA_HYDR_F1][0][0] *= temp;
	 OX_LOG(Logger::LOG_INFO,"Using mismatch repulsion with strength %f",_RNA_HYDR_MIS );

	}


}

//to be removed later, the following function is just for debugging
template<typename number>
number RNA2Interaction<number>::test_huckel(number rbackmod)
{
	number energy = 0;
	if (rbackmod <_debye_huckel_RC) {
			if(rbackmod < _debye_huckel_RHIGH)
			{
			  energy =   exp(rbackmod * _minus_kappa) * (  _debye_huckel_prefactor / rbackmod );// - _debye_huckel_Vrc;

			}
			//use the quadratic-smoothed potential at large distances
			else {
				energy = _debye_huckel_B * SQR(rbackmod -  _debye_huckel_RC);
			}
	    }
		return energy;
}



template<typename number>
number RNA2Interaction<number>::_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)
{
	number cut_factor = 1.0f;
	if(this->_are_bonded(p, q)) return (number) 0.f;
	// for each particle that is on a terminus, halve the charge
    if (  this->_debye_huckel_half_charged_ends && ( p->n3 == P_VIRTUAL || p->n5 == P_VIRTUAL)) cut_factor *= 0.5f;
    if (  this->_debye_huckel_half_charged_ends && ( q->n3 == P_VIRTUAL || q->n5 == P_VIRTUAL)) cut_factor *= 0.5f;
	

	LR_vector<number> rback = *r + q->int_centers[RNANucleotide<number>::BACK] - p->int_centers[RNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number energy = (number) 0.f;
	
	if (rbackmod <_debye_huckel_RC) {
	  if(rbackmod < _debye_huckel_RHIGH)
	    {
	      energy =   exp(rbackmod * _minus_kappa) * (  _debye_huckel_prefactor / rbackmod ) ;//- _debye_huckel_Vrc;
	    }
	  else {
	    energy = _debye_huckel_B * SQR(rbackmod -  _debye_huckel_RC);
	  }
	 energy *= cut_factor; 
	  //update forces are NOT TESTED !!! NEEDS DEBUGGING+careful tests!!!
	  if(update_forces && energy != 0.)
	    {
	      LR_vector<number> force(0.,0.,0.);
	      LR_vector<number> torqueq(0.,0.,0.);
	      LR_vector<number> torquep(0.,0.,0.);
	      LR_vector<number> rbackdir = rback / rbackmod;
	      // compute the new force
	      if(rbackmod < _debye_huckel_RHIGH)
		{
		  force =  rbackdir *  ( -1.0f *  (_debye_huckel_prefactor *
    	        			exp(_minus_kappa * rbackmod) ) *  (  _minus_kappa / rbackmod   - 1.0f /
			SQR(rbackmod) ) );
		}
	      else {
		force = - rbackdir * (2.0f * _debye_huckel_B * (rbackmod -  _debye_huckel_RC) );
	      }
		force *= cut_factor;
	      // computes the torque on q and p
	      torqueq = q->int_centers[RNANucleotide<number>::BACK].cross(force);
	      torquep = -p->int_centers[RNANucleotide<number>::BACK].cross(force);
	      // updates the forces and the torques.
	      p->force -= force;
	      q->force += force;
	      
	      p->torque += p->orientationT * torquep;
	      q->torque += q->orientationT * torqueq;
	      
	    }
	}
	
	//printf("This is huckel for %d %d and energy is %g, at distance %g \n",q->index,p->index,energy,rbackmod);
	return energy;
}



//repulsion potentials

template<typename number>
number RNA2Interaction<number>::_fX(number r, int type,int n3, int n5) {
	number val = (number) 0;
	if(r < this->F1_RCHIGH[type]) {
		if(r > this->F1_RHIGH[type]) {
			val = this->F1_EPS[type][n3][n5] * this->F1_BHIGH[type] * SQR(r - this->F1_RCHIGH[type]);
		}
		else if(r > this->F1_R0[type]) {
			number tmp = 1 - exp(-(r - this->F1_R0[type]) * this->F1_A[type]);
			val = this->F1_EPS[type][n3][n5] * SQR(tmp) - this->F1_SHIFT[type][n3][n5];
		}
                else  {
			val = -this->F1_SHIFT[type][n3][n5];
		}

	}

	return val;
}

template<typename number>
number RNA2Interaction<number>::_fXD(number r, int type,int n3, int n5) {
	number val = (number) 0;
	if(r < this->F1_RCHIGH[type]) {
		if(r > this->F1_RHIGH[type]) {
			val = 2 * this->F1_BHIGH[type] * (r - this->F1_RCHIGH[type]);
		}
		else if(r > this->F1_R0[type]) {
			number tmp = exp(-(r - this->F1_R0[type]) * this->F1_A[type]);
			val = 2 * (1 - tmp) * tmp * this->F1_A[type];
		}
		//else if(r > F1_RCLOW[type]) {
		//	val = 2 * F1_BLOW[type] * (r - F1_RCLOW[type]);
		//}
	}

	return this->F1_EPS[type][n3][n5] * val;
}





template<typename number>
number RNA2Interaction<number>::_hydrogen_bonding_repulsion(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(this->_are_bonded(p, q)) return (number) 0.f;

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	if(!this->_average) //allow for wobble bp
	{
		if( q->btype + p->btype == 4 && ( (q->type == N_T && p->type == N_G ) || (q->type == N_G && p->type == N_T)  ))
			is_pair = true;
	}
    if (is_pair)
    {
    	return 0;
    }

	LR_vector<number> rhydro = *r + q->int_centers[RNANucleotide<number>::BASE] - p->int_centers[RNANucleotide<number>::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;

	//for(double xx = 0.0; xx < model->RNA_HYDR_RCHIGH+0.02; xx += 0.01)
	//{
	//	printf("%f %f \n",xx,_fX(xx,RNA_HYDR_F1,0,0));
	//}
	//exit(1);

	if (rhydromod < this->model->RNA_HYDR_RCHIGH) {
		// vector, versor and magnitude of the base-base separation
		LR_vector<number> rhydrodir = rhydro / rhydromod;

		// particle axes according to Allen's paper
		LR_vector<number> &a1 = p->orientationT.v1;
		LR_vector<number> &a3 = p->orientationT.v3;
		LR_vector<number> &b1 = q->orientationT.v1;
		LR_vector<number> &b3 = q->orientationT.v3;

		// angles involved in the HB interaction
		number cost1 = -a1 * b1;
		number cost2 = -b1 * rhydrodir;
		number cost3 =  a1 * rhydrodir;

		number cost4 = a3 * b3;
		number cost7 = -b3 * rhydrodir;
		number cost8 = a3 * rhydrodir;

		 // functions called at their relevant arguments

		number f1 = this->_fX(rhydromod,RNA_HYDR_F1,0,0);

		number f4t1 = this->_query_mesh (cost1, this->_mesh_f4[RNA_HYDR_F4_THETA1]);
		number f4t2 = this->_query_mesh (cost2, this->_mesh_f4[RNA_HYDR_F4_THETA2]);
		number f4t3 = this->_query_mesh (cost3, this->_mesh_f4[RNA_HYDR_F4_THETA3]);

		number f4t4 = this->_query_mesh (cost4, this->_mesh_f4[RNA_HYDR_F4_THETA4]);
		number f4t7 = this->_query_mesh (cost7, this->_mesh_f4[RNA_HYDR_F4_THETA7]);
		number f4t8 = this->_query_mesh (cost8, this->_mesh_f4[RNA_HYDR_F4_THETA8]);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D;
            f1D  = this->_fXD(rhydromod,RNA_HYDR_F1,0,0);

			number f4t1Dsin = this->_query_meshD (cost1, this->_mesh_f4[RNA_HYDR_F4_THETA1]);
			number f4t2Dsin = this->_query_meshD (cost2, this->_mesh_f4[RNA_HYDR_F4_THETA2]);
			number f4t3Dsin = -this->_query_meshD (cost3, this->_mesh_f4[RNA_HYDR_F4_THETA3]);

			number f4t4Dsin = -this->_query_meshD (cost4, this->_mesh_f4[RNA_HYDR_F4_THETA4]);
			number f4t7Dsin = this->_query_meshD (cost7, this->_mesh_f4[RNA_HYDR_F4_THETA7]);
			number f4t8Dsin = -this->_query_meshD (cost8, this->_mesh_f4[RNA_HYDR_F4_THETA8]);

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
			force += (b1 + rhydrodir * cost2) / rhydromod * fact;
			dir = rhydrodir.cross(b1);

			torqueq += -dir * fact;

			// TETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) / rhydromod * f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			LR_vector<number> t3dir = rhydrodir.cross(a1);
			torquemod = - f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			LR_vector<number> t7dir = rhydrodir.cross(b3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force +=  (a3 - rhydrodir * cost8) / rhydromod * f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			LR_vector<number> t8dir = rhydrodir.cross(a3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[RNANucleotide<number>::BASE].cross(force);
			torqueq += q->int_centers[RNANucleotide<number>::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}


template<typename number>
number RNA2Interaction<number>::_hydrogen_bonding(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(this->_are_bonded(p, q)) return (number) 0.f;

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	if(!this->_average) //allow for wobble bp
	{
		if( q->btype + p->btype == 4 && ( (q->type == N_T && p->type == N_G ) || (q->type == N_G && p->type == N_T)  ))
			is_pair = true;
	}

	if(is_pair)
	{
		return this->RNAInteraction<number>::_hydrogen_bonding(p,q,r,update_forces);
	}
	else if (_mismatch_repulsion)
	{
		return _hydrogen_bonding_repulsion(p,q,r,update_forces);
	}
	else
	{
		return 0;
	}
}





template class RNA2Interaction<float>;
template class RNA2Interaction<double>;

