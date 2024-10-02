#include "DNA2ModInteraction.h"

#include "../Particles/DNANucleotide.h"

template<typename number>
DNA2ModInteraction<number>::DNA2ModInteraction() : DNA2Interaction<number>() {
	_hb_multiplier_1 = 1.;
	_mod_stacking_roll_1 = 0.;
	_mod_stacking_r_roll_1 = 0.;
	_mod_stacking_tilt_1 = 0.;
	_mod_stacking_multiplier_1 = 1.;
	_hb_multiplier_2 = 1.;
	_mod_stacking_roll_2 = 0.;
	_mod_stacking_r_roll_2 = 0.;
	_mod_stacking_tilt_2 = 0.;
	_mod_stacking_multiplier_2 = 1.;
	_a_hb_multiplier = _a_stacking_roll = _a_stacking_r_roll = _a_stacking_tilt = _a_stacking_multiplier = NULL;
}
template<typename number>
DNA2ModInteraction<number>::~DNA2ModInteraction(){
	delete[] _a_hb_multiplier;
	delete[] _a_stacking_roll;
	delete[] _a_stacking_r_roll;
	delete[] _a_stacking_tilt;
	delete[] _a_stacking_multiplier;
}
template<typename number>
void DNA2Interaction<number>::init() {
	DNAInteraction<number>::init();
}

template<typename number>
std::vector<int> DNA2ModInteraction<number>::unsafeGetParticlesFromString(std::string particle_string, char const * identifier) {
	// The interaction doesn't have access to the particles array during initialisation phase. Therefore we can't use getParticlesFromString, and we have to use this version that doesn't check the input sanity. 
	
	// first remove all the spaces from the string, so that the parsing goes well.
	particle_string.erase(remove_if(particle_string.begin(), particle_string.end(), static_cast<int (*)(int)>( isspace )), particle_string.end());

	std::vector<std::string> temp = Utils::split (particle_string.c_str(), ',');
	std::vector<int> particles_index; //declare as empty

	for( std::vector<std::string>::size_type i = 0; i < temp.size(); i++){
		bool found_dash = temp[i].find('-') != std::string::npos;
		// if the string contains a dash, then it has to be interpreted as a list of particles
		// unless it's a negative number

		//if (found_dash && strcmp("-1",temp[i].c_str()) != 0 ){
		if (found_dash && '-'!= temp[i].c_str()[0] ){
			// get the two indices p0 and p1 and check they make sense
			std::vector<std::string> p0_p1_index = Utils::split(temp[i].c_str(),'-');

			int p[2]={0};
			// check whether the p0 and p1 keys can be understood, and set them
			for (int ii = 0; ii < 2; ii++){
				if ( Utils::is_integer(p0_p1_index[ii])){
					p[ii] = atoi(p0_p1_index[ii].c_str());
				}
				if ( ! Utils::is_integer(p0_p1_index[ii]))
					throw oxDNAException("In %s I couldn't interpret particle identifier \"%s\" used as a boundary particle.",identifier,p0_p1_index[ii].c_str());
			}
			// make sure that p[2] is greater than p[1]
			if (p[0] > p[1])
					throw oxDNAException("In %s particle identifier \"%s\" must be such that the first particle has a lower index than the second one when using function %s",identifier,temp[i].c_str(), __FUNCTION__);

			// add all the particles between p0 and p1 (extremes included)
			int j = p[0];
			bool found_p1 = false;
			do{
				particles_index.push_back(j);
				if (j == p[1]){
					found_p1 = true;
				}
				j++;

			} while( j != p[0] && !found_p1);
			// check that it hasn't got to either the end of the strand or back to p1
			if(!found_p1){
				throw oxDNAException("In %s I couldn't get from particle %d to particle %d.",identifier,p[0],p[1]);
			}

		}
		else if ( temp[i] == "all"){
			particles_index.push_back(-1);
		}
		else{//add it to the vector, and make sure that the identifier is not an unidentified string
			if ( temp[i] != "-1" && ! Utils::is_integer(temp[i])){
				throw oxDNAException("In %s I couldn't interpret particle identifier \"%s\".",identifier,temp[i].c_str());

			}
			int j = atoi(temp[i].c_str());

			particles_index.push_back(j);
		}

	}
	// check that if -1 is present then that's the only key - something must be wrong if you
	// specified -1 (all particles) and then some more particles.
	if (std::find(particles_index.begin(),particles_index.end(),-1) != particles_index.end() && particles_index.size()>1){
		throw oxDNAException("In %s there is more than one particle identifier, including -1 or \"all\". If either -1 or \"all\" are used as particle identifiers then they have to be the only one, as both translate to \"all the particles\". Dying badly.",identifier);
	}
	// check that no particle appears twice
	for( std::vector<int>::size_type i = 0; i < particles_index.size(); i++){
		for( std::vector<int>::size_type j = i+1; j < particles_index.size(); j++){
			if ( particles_index[i] == particles_index[j] ){
				throw oxDNAException("In %s particle index %d appears twice (both at position %d and at position %d), but each index can only appear once. Dying badly.",identifier,particles_index[i],i+1,j+1);
			}
		}
	}
	// finally return the vector.
	return particles_index;
}
template<typename number>
void DNA2ModInteraction<number>::get_settings(input_file &inp) {
	DNA2Interaction<number>::get_settings(inp);
 	//int index = 0;
	std::string temp;
	std::string backend="CPU";
	getInputString(&inp,"backend",backend,0);
	/*
	if (KEY_FOUND == getInputInt(&inp, "mod_nucleotide_group_1", &index, 0))
		_vec_group_1.push_back(index);
	*/
	if (KEY_FOUND == getInputNumber(&inp, "mod_stacking_angle_1", &_mod_stacking_roll_1,0))
		throw oxDNAException("Argoment mod_stacking_angle_1 now substituted by mod_stacking_angle_roll_1");
	if (KEY_FOUND == getInputNumber(&inp, "mod_stacking_angle_2", &_mod_stacking_roll_2,0))
		throw oxDNAException("Argoment mod_stacking_angle_2 now substituted by mod_stacking_angle_roll_2");
	if (KEY_FOUND == getInputString(&inp, "mod_nucleotide_group_1", temp, 0)){
		_vec_group_1 = unsafeGetParticlesFromString(temp,"DNA2_mod interaction");
	
		getInputNumber(&inp, "mod_hb_multiplier_1", &_hb_multiplier_1,0);
		getInputNumber(&inp, "mod_stacking_roll_1", &_mod_stacking_roll_1,0);
		if( getInputNumber(&inp, "mod_stacking_r_roll_1", &_mod_stacking_r_roll_1,0) == KEY_FOUND and backend == "CUDA")
			throw oxDNAException("mod_stacking_r_roll_1 has been disabled on CUDA to try and prevent Bus Error.");
		getInputNumber(&inp, "mod_stacking_tilt_1", &_mod_stacking_tilt_1,0);
		getInputNumber(&inp, "mod_stacking_multiplier_1", &_mod_stacking_multiplier_1,0);
		_mod_stacking_roll_1 *= M_PI / 180.;
		_mod_stacking_r_roll_1 *= M_PI / 180.;
		_mod_stacking_tilt_1 *= M_PI / 180.;
		printf(" mod stacking roll_1 %g\n",_mod_stacking_roll_1);
		printf(" mod stacking r_roll_1 %g\n",_mod_stacking_r_roll_1);
		printf(" mod stacking tilt_1 %g\n",_mod_stacking_tilt_1);
		printf(" mod stacking multiplier_1 %g\n",_mod_stacking_multiplier_1);
		printf(" mod hb multiplier 1 %g\n",_hb_multiplier_1);
	}
	if (KEY_FOUND == getInputString(&inp, "mod_nucleotide_group_2", temp, 0)){
		_vec_group_2 = unsafeGetParticlesFromString(temp,"DNA2_mod interaction");
	
		getInputNumber(&inp, "mod_hb_multiplier_2", &_hb_multiplier_2,0);
		getInputNumber(&inp, "mod_stacking_roll_2", &_mod_stacking_roll_2,0);
		if( getInputNumber(&inp, "mod_stacking_r_roll_2", &_mod_stacking_r_roll_2,0) == KEY_FOUND and backend == "CUDA")
			throw oxDNAException("mod_stacking_r_roll_2 has been disabled on CUDA to try and prevent Bus Error.");
		getInputNumber(&inp, "mod_stacking_tilt_2", &_mod_stacking_tilt_2,0);
		getInputNumber(&inp, "mod_stacking_multiplier_2", &_mod_stacking_multiplier_2,0);
		_mod_stacking_roll_2 *= M_PI / 180.;
		_mod_stacking_r_roll_2 *= M_PI / 180.;
		_mod_stacking_tilt_2 *= M_PI / 180.;
		printf(" mod stacking roll_2 %g\n",_mod_stacking_roll_2);
		printf(" mod stacking r_roll_2 %g\n",_mod_stacking_r_roll_2);
		printf(" mod stacking tilt_2 %g\n",_mod_stacking_tilt_2);
		printf(" mod stacking multiplier_2 %g\n",_mod_stacking_multiplier_2);
		printf(" mod hb multiplier 2 %g\n",_hb_multiplier_2);
	}
	// initialise the arrays
	int N = this->get_N_from_topology();
	_a_hb_multiplier = new number[N];
	_a_stacking_tilt = new number[N];
	_a_stacking_multiplier = new number[N];
	_a_stacking_roll = new number[N];
	_a_stacking_r_roll = new number[N];

	for (int i = 0; i < N; i++){
		if (_is_integer_in_group(i, _vec_group_1)){
			_a_hb_multiplier[i] = _hb_multiplier_1;
			_a_stacking_tilt[i] = _mod_stacking_tilt_1;
			_a_stacking_multiplier[i] = _mod_stacking_multiplier_1;
			_a_stacking_roll[i] = _mod_stacking_roll_1;
			_a_stacking_r_roll[i] = _mod_stacking_r_roll_1;
		}
		else if (_is_integer_in_group(i, _vec_group_2)){
			_a_hb_multiplier[i] = _hb_multiplier_2;
			_a_stacking_tilt[i] = _mod_stacking_tilt_2;
			_a_stacking_multiplier[i] = _mod_stacking_multiplier_2;
			_a_stacking_roll[i] = _mod_stacking_roll_2;
			_a_stacking_r_roll[i] = _mod_stacking_r_roll_2;
		}
		else{
			_a_hb_multiplier[i] = 1.;
			_a_stacking_tilt[i] = 0.;
			_a_stacking_multiplier[i] = 1.;
			_a_stacking_roll[i] = 0.;
			_a_stacking_r_roll[i] = 0.;
		}

	}

}

template<typename number>
number DNA2ModInteraction<number>::_stacking(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(!DNAInteraction<number>::_check_bonded_neighbour(&p, &q, r)) return (number) 0.f;

	LR_vector<number> &a1 = p->orientationT.v1;
	LR_vector<number> &a2 = p->orientationT.v2;
	LR_vector<number> &a3 = p->orientationT.v3;
	LR_vector<number> b1_1 = q->orientationT.v1;
	LR_vector<number> b2_1 = q->orientationT.v2;
	
	// The b3 vector is rolled and tilted as needed.
	LR_vector<number> b3_1 = get_rotated_T_v3(q);
	bool q_in_group_1 = _is_particle_in_group(q, _vec_group_1);
	bool q_in_group_2 = _is_particle_in_group(q, _vec_group_2);
	LR_vector<number> b3(0.,0.,0.);
	LR_vector<number> b2(0.,0.,0.);
	LR_vector<number> b1(0.,0.,0.);
	LR_matrix<number> R;
	number stack_roll=0, stack_r_roll=0, stack_tilt=0, hb_multi=1., stack_multi=1.;
	if (q_in_group_1){
		R = rotationMatrixAroundVersorByAngle(q->orientationT.v1, _mod_stacking_roll_1);
		R = rotationMatrixAroundVersorByAngle(q->orientationT.v2, _mod_stacking_tilt_1) * R;
		R = (R*q->orientation).get_transpose();
		b1 = R.v1;
		b2 = R.v2;
		b3 = R.v3;
		stack_roll=_mod_stacking_roll_1;
		stack_r_roll=_mod_stacking_r_roll_1;
		stack_tilt=_mod_stacking_tilt_1;
		hb_multi = _hb_multiplier_1;
		stack_multi = _mod_stacking_multiplier_1;
	}
	else if (q_in_group_2){
		R = rotationMatrixAroundVersorByAngle(q->orientationT.v1, _mod_stacking_roll_2);
		R = rotationMatrixAroundVersorByAngle(q->orientationT.v2, _mod_stacking_tilt_2) * R;
		R = (R*q->orientation).get_transpose();
		b1 = R.v1;
		b2 = R.v2;
		b3 = R.v3;
		stack_roll=_mod_stacking_roll_2;
		stack_r_roll=_mod_stacking_r_roll_2;
		stack_tilt=_mod_stacking_tilt_2;
		hb_multi = _hb_multiplier_2;
		stack_multi = _mod_stacking_multiplier_2;
		//b3 = R*q->orientationT.v3;
	}
	else b3 = q->orientationT.v3;
	// TODO: just use the array values for the angles every time this is necessary, and remove the choice structure above, as soon as it's clear that the following exception is never triggered. 
	// TODO: the same should be done in hydrogen_bonding
	if (stack_roll != _a_stacking_roll[q->index] or stack_r_roll != _a_stacking_r_roll[q->index] or stack_tilt != _a_stacking_tilt[q->index] or stack_multi != _a_stacking_multiplier[q->index] or hb_multi != _a_hb_multiplier[q->index]){
		throw oxDNAException("SBANGABANGA!!!!!\n");
	}

	number bcos = b3 * b3_1;
	number angle = LRACOS(bcos)*180./M_PI;
	number bcos_2 = b2 * b2_1;
	number angle_2 = LRACOS(bcos_2)*180./M_PI;
	number bcos_1 = b1 * b1_1;
	number angle_1 = LRACOS(bcos_1)*180./M_PI;
	if ((q_in_group_1 or q_in_group_2) and false){
		printf("p[%d].v1 = %g %g %g\nrotated.v1 = %g %g %g\nnew_rota.v1 = %g %g %g, angle %g, cos %g\n",q->index,q->orientationT.v1.x,q->orientationT.v1.y,q->orientationT.v1.z, b1_1.x, b1_1.y, b1_1.z, b1.x, b1.y, b1.z, angle_1, bcos_1);
		printf("p[%d].v2 = %g %g %g\nrotated.v2 = %g %g %g\nnew_rota.v2 = %g %g %g, angle %g, cos %g\n",q->index,q->orientationT.v2.x,q->orientationT.v2.y,q->orientationT.v2.z, b2_1.x, b2_1.y, b2_1.z, b2.x, b2.y, b2.z, angle_2, bcos_2);
		printf("p[%d].v3 = %g %g %g\nrotated.v3 = %g %g %g\nnew_rota.v3 = %g %g %g, angle %g, cos %g\n",q->index,q->orientationT.v3.x,q->orientationT.v3.y,q->orientationT.v3.z, b3_1.x, b3_1.y, b3_1.z, b3.x, b3.y, b3.z, angle, bcos);
		LR_matrix<number> m = q->orientationT;
		printf("orientationT:\n%g\t%g\t%g\n%g\t%g\t%g\n%g\t%g\t%g\n",m.v1.x,m.v1.y,m.v1.z, m.v2.x,m.v2.y,m.v2.z, m.v3.x,m.v3.y,m.v3.z);
		LR_matrix<number> rr = R*q->orientationT;
		printf("orientationT:\n%g\t%g\t%g\n%g\t%g\t%g\n%g\t%g\t%g\n",rr.v1.x,rr.v1.y,rr.v1.z, rr.v2.x,rr.v2.y,rr.v2.z, rr.v3.x,rr.v3.y,rr.v3.z);
		exit(1);
	}

	//LR_vector<number> &b3 = q->orientationT.v3;

	LR_vector<number> computed_r;
	if (r == NULL) {
		computed_r = q->pos - p->pos;
		r = &computed_r;
	}

	// This is the position the backbone would have with major-minor grooves the same width.
	// We need to do this to implement different major-minor groove widths because rback is
	// used as a reference point for things that have nothing to do with the actual backbone
	// position (in this case, the stacking interaction).
	LR_vector<number> rbackref = _computed_r + b1 * POS_BACK - a1 * POS_BACK;
	if (q_in_group_1) rbackref = rotateVectorAroundVersor(rbackref, q->orientationT.v1, _mod_stacking_r_roll_1);
	else if (q_in_group_2) rbackref = rotateVectorAroundVersor(rbackref, q->orientationT.v1, _mod_stacking_r_roll_2);
	number rbackrefmod = rbackref.module();

	LR_vector<number> rstack = _computed_r + q->int_centers[DNANucleotide<number>::STACK] - p->int_centers[DNANucleotide<number>::STACK];
	if (q_in_group_1) rstack = rotateVectorAroundVersor(rstack, q->orientationT.v1, _mod_stacking_r_roll_1);
	else if (q_in_group_2) rstack = rotateVectorAroundVersor(rstack, q->orientationT.v1, _mod_stacking_r_roll_2);
	number rstackmod = rstack.module();
	LR_vector<number> rstackdir = rstack / rstackmod;

	number cost4 =  a3 * b3;
	number cost5 =  a3 * rstackdir;
	number cost6 = -b3 * rstackdir;
	number cosphi1 = a2 * rbackref / rbackrefmod;
	number cosphi2 = b2 * rbackref / rbackrefmod;

	// functions and their derivatives needed for energies and forces
	number f1     = this->_f1(rstackmod, STCK_F1, q->type, p->type);
	number f4t4   = DNA2Interaction<number>::_custom_f4 (cost4, STCK_F4_THETA4);
	number f4t5   = DNA2Interaction<number>::_custom_f4 (-cost5, STCK_F4_THETA5);
	number f4t6   = DNA2Interaction<number>::_custom_f4 (cost6, STCK_F4_THETA6);
	number f5phi1 = this->_f5(cosphi1, STCK_F5_PHI1);
	number f5phi2 = this->_f5(cosphi2, STCK_F5_PHI2);

	number energy = f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2;

	if(update_forces && energy != (number) 0.f) {
		LR_vector<number> torquep(0, 0, 0);
		LR_vector<number> torqueq(0, 0, 0);

		// these are the derivatives of f4 with respect to t4, t5 and t6 over sin(t*). If t* is too small
		// we have numerical instabilities so we use the fact that sin(x) = x for x -> 0
		number f1D      = this->_f1D(rstackmod, STCK_F1, q->type, p->type);
		number f4t4Dsin = -DNA2Interaction<number>::_custom_f4D (cost4, STCK_F4_THETA4);
		number f4t5Dsin = -DNA2Interaction<number>::_custom_f4D (-cost5, STCK_F4_THETA5);
		number f4t6Dsin = -DNA2Interaction<number>::_custom_f4D (cost6, STCK_F4_THETA6);
		number f5phi1D  = this->_f5D(cosphi1, STCK_F5_PHI1);
		number f5phi2D  = this->_f5D(cosphi2, STCK_F5_PHI2);

		// RADIAL
		LR_vector<number> force = -rstackdir * (f1D * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2);

		// THETA 5
		force += -(a3 - rstackdir * cost5) * (f1 * f4t4 * f4t5Dsin * f4t6 * f5phi1 * f5phi2 / rstackmod);

		// THETA 6
		force += -(b3 + rstackdir * cost6) * (f1 * f4t4 * f4t5 * f4t6Dsin * f5phi1 * f5phi2 / rstackmod);

		// COS PHI 1
		// here particle p is referred to using the a while particle q is referred with the b
		number gamma = POS_STACK - POS_BACK;
		number rbackrefmodcub = rbackrefmod * rbackrefmod * rbackrefmod;

		number ra2 = rstackdir*a2;
		number ra1 = rstackdir*a1;
		number rb1 = rstackdir*b1;
		number a2b1 = a2*b1;

		number parentesi = rstackmod*ra2 - a2b1*gamma;

		number dcosphi1dr    = (SQR(rstackmod)*ra2 - ra2*SQR(rbackrefmod) - rstackmod*(a2b1 + ra2*(-ra1 + rb1))*gamma + a2b1*(-ra1 + rb1)*SQR(gamma))/ rbackrefmodcub;
		number dcosphi1dra1  =  rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi1dra2  = -rstackmod / rbackrefmod;
		number dcosphi1drb1  = -rstackmod * gamma * parentesi / rbackrefmodcub;
		number dcosphi1da1b1 = -SQR(gamma) * parentesi / rbackrefmodcub;
		number dcosphi1da2b1 =  gamma / rbackrefmod;

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

		number dcosphi2dr    =  (parentesi * (rstackmod + (rb1 - ra1)*gamma) - ra2*SQR(rbackrefmod)) / rbackrefmodcub;
		number dcosphi2dra1  = -rstackmod*gamma*(rstackmod*ra2 + a2b1*gamma) / rbackrefmodcub;
		number dcosphi2dra2  = -rstackmod / rbackrefmod;
		number dcosphi2drb1  =  rstackmod*gamma* parentesi / rbackrefmodcub;
		number dcosphi2da1b1 = -SQR(gamma)* parentesi / rbackrefmodcub;
		number dcosphi2da2b1 = -gamma / rbackrefmod;

		// this force part has a minus because of the definition of cos(phi2) which is not consistent with the derivatives
		// above (i.e. all the derivatives should have a minus sign in front and the force below shouldn't)
		number force_part_phi2 = -f1 * f4t4 * f4t5 * f4t6 * f5phi1 * f5phi2D;

		force += -force_part_phi2 * (rstackdir * dcosphi2dr +
									  ((b2 - rstackdir * ra2) * dcosphi2dra2 +
									   (b1 - rstackdir * ra1) * dcosphi2dra1 +
									   (a1 - rstackdir * rb1) * dcosphi2drb1) / rstackmod);

		// Add the contribution from all the forces to the stored particles' forces
		p->force -= force * stack_multi;
		q->force += force * stack_multi;

		torquep -= p->int_centers[DNANucleotide<number>::STACK].cross(force);
		torqueq += q->int_centers[DNANucleotide<number>::STACK].cross(force);

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
		p->torque += p->orientationT * torquep * stack_multi;
		q->torque += q->orientationT * torqueq * stack_multi;
	}

	return energy * stack_multi;
}

template<typename number>
number DNA2ModInteraction<number>::_hydrogen_bonding(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return (number) 0.f;

	// true if p and q are Watson-Crick-like pairs
	bool is_pair = (q->btype + p->btype == 3);
	number hb_multi = get_hb_multiplier(p,q);

	LR_vector<number> rhydro = _computed_r + q->int_centers[DNANucleotide<number>::BASE] - p->int_centers[DNANucleotide<number>::BASE];
	number rhydromod = rhydro.module();
	number energy = (number) 0.f;
	if (is_pair && HYDR_RCLOW < rhydromod && rhydromod < HYDR_RCHIGH) {
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
		number f1   = hb_multi * this->_f1(rhydromod, HYDR_F1, q->type, p->type);
		number f4t1 = DNA2Interaction<number>::_custom_f4 (cost1, HYDR_F4_THETA1);
		number f4t2 = DNA2Interaction<number>::_custom_f4 (cost2, HYDR_F4_THETA2);
		number f4t3 = DNA2Interaction<number>::_custom_f4 (cost3, HYDR_F4_THETA3);

		number f4t4 = DNA2Interaction<number>::_custom_f4 (cost4, HYDR_F4_THETA4);
		number f4t7 = DNA2Interaction<number>::_custom_f4 (cost7, HYDR_F4_THETA7);
		number f4t8 = DNA2Interaction<number>::_custom_f4 (cost8, HYDR_F4_THETA8);

		energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;

		// makes sense, since the above functions may return 0. exactly
		if(update_forces && energy != 0.) {
			LR_vector<number> force(0, 0, 0);
			LR_vector<number> torquep(0, 0, 0);
			LR_vector<number> torqueq(0, 0, 0);

			// derivatives called at the relevant arguments
			number f1D      =  hb_multi * this->_f1D(rhydromod, HYDR_F1, q->type, p->type);
			number f4t1Dsin = DNA2Interaction<number>::_custom_f4D (cost1, HYDR_F4_THETA1);
			number f4t2Dsin = DNA2Interaction<number>::_custom_f4D (cost2, HYDR_F4_THETA2);
			number f4t3Dsin = -DNA2Interaction<number>::_custom_f4D (cost3, HYDR_F4_THETA3);

			number f4t4Dsin = -DNA2Interaction<number>::_custom_f4D (cost4, HYDR_F4_THETA4);
			number f4t7Dsin = DNA2Interaction<number>::_custom_f4D (cost7, HYDR_F4_THETA7);
			number f4t8Dsin = -DNA2Interaction<number>::_custom_f4D (cost8, HYDR_F4_THETA8);

			// RADIAL PART
			force = - rhydrodir * (f1D * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8);

			// THETA4; t4 = LRACOS (a3 * b3);
			LR_vector<number> dir = a3.cross(b3);
			number torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4Dsin * f4t7 * f4t8 ;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA1; t1 = LRACOS (-a1 * b1);
			dir = a1.cross(b1);
			torquemod = - f1 * f4t1Dsin * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 ;

			torquep -= dir * torquemod;
			torqueq += dir * torquemod;

			// THETA2; t2 = LRACOS (-b1 * rhydrodir);
			number fact = f1 * f4t1 * f4t2Dsin * f4t3 * f4t4 * f4t7 * f4t8;
			force += (b1 + rhydrodir * cost2) * (fact / rhydromod);
			dir = rhydrodir.cross(b1);

			torqueq += -dir * fact;

			// THETA3; t3 = LRACOS (a1 * rhydrodir);
			force += (a1 - rhydrodir * cost3) * (f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8 / rhydromod);

			LR_vector<number> t3dir = rhydrodir.cross(a1);
			torquemod = - f1 * f4t1 * f4t2 * f4t3Dsin * f4t4 * f4t7 * f4t8;

			torquep += t3dir * torquemod;

			// THETA7; t7 = LRACOS (-rhydrodir * b3);
			force += (b3 + rhydrodir * cost7) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8 / rhydromod);

			LR_vector<number> t7dir = rhydrodir.cross(b3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7Dsin * f4t8;

			torqueq += t7dir * torquemod;

			// THETA 8; t8 = LRACOS (rhydrodir * a3);
			force +=  (a3 - rhydrodir * cost8) * (f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin / rhydromod);

			LR_vector<number> t8dir = rhydrodir.cross(a3);
			torquemod = - f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8Dsin;

			torquep += t8dir * torquemod;

			// final update of forces and torques for HB
			p->force -= force;
			q->force += force;

			torquep -= p->int_centers[DNANucleotide<number>::BASE].cross(force);
			torqueq += q->int_centers[DNANucleotide<number>::BASE].cross(force);

			// we need torques in the reference system of the particle
			p->torque += p->orientationT * torquep;
			q->torque += q->orientationT * torqueq;
		}
	}

	return energy;
}

extern "C" BaseInteraction<float> *make_DNA2ModInteraction_float() { return new DNA2ModInteraction<float>(); }
extern "C" BaseInteraction<double> *make_DNA2ModInteraction_double() { return new DNA2ModInteraction<double>(); }
template class DNA2ModInteraction<float>;
template class DNA2ModInteraction<double>;
