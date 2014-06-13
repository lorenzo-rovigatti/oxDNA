#include "DNAInteraction2.h"

#include "../Particles/DNANucleotide.h"


template<typename number> // empty constructor
DNA2Interaction<number>::DNA2Interaction() : DNAInteraction<number>() {
  // question - what is this?
	this->_int_map[DEBYE_HUCKEL] = (number (DNAInteraction<number>::*)(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)) &DNA2Interaction<number>::_debye_huckel;

	// log the interaction type
	OX_LOG(Logger::LOG_INFO,"Running modification of oxDNA with additional Debye-Huckel potential (second version)");
}



template<typename number>
number DNA2Interaction<number>::pair_interaction_nonbonded(
		BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r,
		bool update_forces) {

	LR_vector<number> computed_r(0, 0, 0);
	// if r == null, compute the distance taking into account periodic boundary conditions
	if(r == NULL) {
		computed_r = q->pos.minimum_image(p->pos, this->_box_side);
		r = &computed_r;
	}

	// if the distance is beyond the cutoff, return 0
        if (r->norm() >= this->_sqr_rcut)
		return (number) 0.f;
	// compute the interaction energy as always ...
	number energy = DNAInteraction<number>::pair_interaction_nonbonded(p, q, r,
			update_forces);

	// ... and then add the debye_huckel energy
	energy += _debye_huckel(p, q, r, update_forces);

	return energy;
}

// get the settings from the input file
template<typename number>
void DNA2Interaction<number>::get_settings(input_file &inp) {
	DNAInteraction<number>::get_settings(inp);

	float salt;
	int mismatches;
	float prefactor; // this is the strength of the interaction
	float lambdafactor; //Lambda is _debye_huckel_LAMBDAFACTOR / salt^0.5
	float rh;
	int lambda_T_dependent = 0;
        int charged_n3 = 0;
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
	} else {// change for RNA since the simulation units are different
			_debye_huckel_lambdafactor = 0.3616455;
	}


	if (getInputBoolAsInt(&inp, "lambda_T_dependent", &lambda_T_dependent, 0) == KEY_FOUND  && lambda_T_dependent) {
			OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel with temperature dependent lambda at T = %f",this->_T);
			_debye_huckel_lambdafactor *= sqrt(this->_T / 0.1f) ;
	} 

	// read the prefactor to the potential, or set it to the default value
	if (getInputFloat(&inp, "dh_strength", &prefactor, 0) == KEY_FOUND) {
				_debye_huckel_prefactor = (float) prefactor;
			} else {
				_debye_huckel_prefactor = 0.510473f;
	}
	// read the cutoff distance (todo - set the default value to something 
	// which makes sense, maybe that depends on lambdafactor and salt concentration
	if(getInputFloat(&inp, "debye_huckel_rhigh", &rh, 0) == KEY_FOUND) {
	                        _debye_huckel_RHIGH = (float) rh; 
	                } else {
	                        _debye_huckel_RHIGH = 3.0f; 
	}
	// read the mismatch_repulsion flag (not implemented yet)
	if(getInputBoolAsInt(&inp, "mismatch_repulsion", &mismatches, 0) == KEY_FOUND) {
			this->_mismatch_repulsion = (bool) mismatches;
	}
	else {
			this->_mismatch_repulsion = false;
	}
        // whether to put a charge on the 3 terminus or not
        if(getInputBoolAsInt(&inp, "dh_charged_n3", &charged_n3, 0) == KEY_FOUND) {
        		this-> _debye_huckel_charged_n3 = (bool) charged_n3;	
	}
	else {
			this-> _debye_huckel_charged_n3 = true;
        }

	//log it 
	OX_LOG(Logger::LOG_INFO,"Running Debye-Huckel at salt_concentration =  %g",this->_salt_concentration);


}

template<typename number> //initialise the interaction
void DNA2Interaction<number>::init() {
	// we choose rcut as the max of the range interaction of excluded
	// volume between backbones and hydrogen bonding
	//number rcutback = 2 * fabs(model->DNA_POS_BACK) + model->DNA_EXCL_RC1;
	//this->_T = T;

        //perform the initialisation as in DNAinteraction
	DNAInteraction<number>::init();

	//compute the DH length lambda
	number lambda = _debye_huckel_lambdafactor / sqrt(_salt_concentration);
	// upper limit of the unsmoothed potential
	//_debye_huckel_RHIGH = 2.0 * lambda - 0.1;this is set in the get_settings function
	// value of V_DH at r = 2 lambda
	_debye_huckel_Vrc  = 0; //|(_debye_huckel_prefactor / (2* lambda))  * exp(-2.0) ;
	_minus_kappa = -1.0/lambda;
	
	//some constants for quick referencing of the variables
	number x = _debye_huckel_RHIGH;
	number q = _debye_huckel_prefactor;
	number l = lambda;
	number V = 0.; //We don't want any offset for the potential
	
	// now the smoothening transition at 2 * lambda - 0.1
	// why is this so complicated? I got something much simpler
	// I would simply have
	// _debye_huckel_B = q*exp(-x/l)/(x*(x-_debye_huckel_RC)*(x-_debye_huckel_RC))
	// which to work properly should be moved after _debye_huckel_RC
	_debye_huckel_B = -(exp(-x/l) * q * q * (x + l)*(x+l) )/(4.*x*x*x * l * l * (-q + exp(x/l) * V*x) );
	// question: Why is this so complicated? I got something much simpler which only depends on r_h and only requires the continuity/differentiability of the potential
	// I would simply have
	// _debye_huckel_RC = x + 2*x/(1+x/l)
	_debye_huckel_RC = x*(q*x + 3. * q* l - 2.0 * exp(x/l) * V * x * l )/(q * (x+l));

	// define the cut-off of the debye-huckel interaction
	// question: what is POS_BACK?

        
	number debyecut =  2.0f * sqrt(SQR(POS_BACK)	) + _debye_huckel_RC;
	if (this->_grooving){
		debyecut = 2.0f * sqrt((POS_MM_BACK1)*(POS_MM_BACK1) + (POS_MM_BACK2)*(POS_MM_BACK2)) + _debye_huckel_RC;
	}
	// if debyecut > rcut, then debyecut is the new rcut
	if (debyecut > this->_rcut)
	{
		this->_rcut = debyecut;
		this->_sqr_rcut = debyecut*debyecut;
	}
	// log the parameters of the Debye-Huckel
	// question: why do we output Vrc?
	OX_LOG(Logger::LOG_INFO,"DEBUGGING: rhigh is %g, Cutoff is %g, RC huckel is %g, B huckel is %g V is %g ",_debye_huckel_RHIGH,this->_rcut,_debye_huckel_RC, _debye_huckel_B,_debye_huckel_Vrc);
	OX_LOG(Logger::LOG_INFO,"DEBUGGING: dh_charged_n3 = %s", _debye_huckel_charged_n3 ? "true" : "false");



/*
	    for (double x = 0.01; x <= _debye_huckel_RC +0.02;  x += 0.02)
	    {
	//    	//energy = _debye_huckel(this->partic)
	    	printf("%f %f \n",x,test_huckel(x));
	    }
	    exit(1);
*/
}

//to be removed later, the following function is just for debugging
template<typename number>
number DNA2Interaction<number>::test_huckel(number rbackmod)
{
	number energy = 0;
	//compute the interaction energy at a given distance
	if (rbackmod <_debye_huckel_RC) {
	  //use the actual potential at short distances
			if(rbackmod < _debye_huckel_RHIGH)
			{
			  energy =   exp(rbackmod * _minus_kappa) * (  _debye_huckel_prefactor / rbackmod ); //- _debye_huckel_Vrc;
			    //commented since we don't want to shift the potential down
			}
			//use the quadratic-smoothed potential at large distances
			else {
				energy = _debye_huckel_B * SQR(rbackmod -  _debye_huckel_RC);
			}
	    }
		return energy;
}



template<typename number>
number DNA2Interaction<number>::_debye_huckel(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces)
{
  // question: this checks if they are adjacent along the backbone? Why don't you want
  // DH in that case?
	if(this->_are_bonded(p, q)) return (number) 0.f;
        if ( ! this->_debye_huckel_charged_n3 && ( p->n3 == P_VIRTUAL || q->n3 == P_VIRTUAL)) return (number) 0.f;
	
	//compute the distance vector between the two particles
	//question: r is used to enforce periodic boundary conditions?
	LR_vector<number> rback = *r + q->int_centers[DNANucleotide<number>::BACK] - p->int_centers[DNANucleotide<number>::BACK];
	number rbackmod = rback.module();
	number energy = (number) 0.f;
	
	// calculate the DH energy (see above for details)
	if (rbackmod <_debye_huckel_RC) {
	  if(rbackmod < _debye_huckel_RHIGH)
	    {
	      energy =   exp(rbackmod * _minus_kappa) * (  _debye_huckel_prefactor / rbackmod ); //- _debye_huckel_Vrc;
	    }
	  else {
	    energy = _debye_huckel_B * SQR(rbackmod -  _debye_huckel_RC);
	  }
	  
	  //update forces are NOT TESTED !!! NEEDS DEBUGGING+careful tests!!!
	  if(update_forces && energy != 0.)
	    {
	      LR_vector<number> force(0.,0.,0.);
	      LR_vector<number> torqueq(0.,0.,0.);
	      LR_vector<number> torquep(0.,0.,0.);
	      LR_vector<number> rbackdir = rback / rbackmod;
	      // compute the new force
	      if(rbackmod < _debye_huckel_RHIGH)
		{// question: why is the force like this? I get a different result.
		  // I would set
		  //force = rbackdir*( 1/rbackmod - _minus_kappa)*_debye_huckel_prefactor*
		  // exp(rbackmod * _minus_kappa)/rbackmod
			force =  rbackdir *  ( -1.0f *  (_debye_huckel_prefactor *
			exp(_minus_kappa * rbackmod) ) *  (  _minus_kappa / rbackmod   - 1.0f /
			SQR(rbackmod) ) );
		}
	      else {
		force = - rbackdir * (2.0f * _debye_huckel_B * (rbackmod -  _debye_huckel_RC) );
	      }
	      // computes the torque on q and p
	      torqueq = q->int_centers[DNANucleotide<number>::BACK].cross(force);
	      torquep = -p->int_centers[DNANucleotide<number>::BACK].cross(force);
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

template class DNA2Interaction<float>;
template class DNA2Interaction<double>;

