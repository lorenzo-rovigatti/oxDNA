#include "PSPInteraction.h"

PSPInteraction::PSPInteraction():BaseInteraction(){
    // Constructor
};

PSPInteraction::~PSPInteraction(){
    // Destructor
};

void PSPInteraction::get_settings(input_file &inp){
    BaseInteraction::get_settings(inp); // Get basic settings
	if(getInputString(&inp,"patchyRcut",temp,0)==KEY_FOUND){
		patchyRcut =stod(temp);
		patchyRcut2=SQR(patchyRcut);
		std::cout<<"New cutoff value for the patchy interaction = "<<patchyRcut<<std::endl;
	}
	if(getInputString(&inp,"rcut",temp,0)==KEY_FOUND){ //set the rcut from the input file
		_rcut=stod(temp);
        _sqr_rcut=sqrt(_rcut);
		std::cout<<"The new value of rcut = "<< _rcut<<std::endl;
	}else{
        _rcut=-1; // set the rcut after reading the topology. 2.5x largest radius
    }
}

void PSPInteraction::init(){
	number r8b10 = powf(patchyRcut, (number) 8.f) / patchyPowAlpha;
    patchyEcut = -1.001f * exp(-(number) 0.5f * r8b10 * patchyRcut2);

	//making sure the lock is 0
	#pragma omp parallel for collapse(2)
	for(i=0;i<PSPmaxParticles;i++){
		for(j=0;j<PSPmaxPatchOnParticle;j++){
			patchLock[i][j]=false;
		}
	}
}

void PSPInteraction::allocate_particles(std::vector<BaseParticle *> &particles){
    // Add particle to the system
    particles.resize(totPar);
	for(i=0;i<totPar;i++){
		particles[i] = new CCGParticle();
		particles[i]->index=i;
		particles[i]->type=0;
	}
}

void PSPInteraction::check_input_sanity(std::vector<BaseParticle *> &particles){
    // Check all the input file are correct.
}

// void PSPInteraction::coordinateConverter(){
// 	#pragma omp parallel for
// 	for(i=0;i<PSPmaxPatchColor;i++){
// 		if(patchType[i]==1){
// 			patches[i][2]*=
// 		}
// 	}
// }

void PSPInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles){
    // Read the top file
	std::ifstream topology(this->_topology_filename);
	if(!topology.good()) throw oxDNAException("Something wrong with the topology file: %s",this->_topology_filename);
	std::getline(topology,temp);
	std::stringstream head(temp);
	head >> totPar >> strads>> totParticleColor >> totPatchColor; //Read header here
	*N_strands = strads;

	allocate_particles(particles);
	i=0;j=0;
	int countP=0, countC=0;
	while(std::getline(topology,line)){
		if(line.empty()||line[0]=='#') continue; //skip empty line and comment
		if(line[0]=='i'){
			if(line[1]=='P'){ //these are Patchy informations
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					body>>temp;
					if(j==1) patchType[countP]=stoi(temp); //get the patch type.
					if(j>1 && j<6) patches[countP][j-2]=stof(temp);
					j++;
				}
				countP++;
				// continue;
			}
			if(line[1]=='C'){
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					body>>temp; //skip word "iC"
					// if(j==0) continue; //skip the first word
					if(j>0) particlePatches[countC][j]=stoi(temp);
					j++;
				}
				countC++;
				// continue;
			}
			continue;
		}
		auto *p = static_cast<CCGParticle *>(particles[i]);
		std::stringstream body(line);
		j=0;
		connections[i][0]=0; //setting the number of neighbours to 0 before the loop
		r0[i][0]=0; // all the eq radius are different
		k0[i][0]=0; // all the spring constant are different
		while(body.tellg()!=-1){
			body>>temp;
			if(j==0) p->type=stoi(temp);
			if(j==1){
				p->strand_id=stoi(temp); // not sure if needed, more like backward compatible
				particleTopology[i][0]=std::stoi(temp);
			};
			if(j==2){
				p->btype=stoi(temp); // particle color is also btype, this is redundant, to make it backward compatible
				particleTopology[i][1]=std::stoi(temp);
			};
			if(j==3){
				p->radius=(number)stof(temp);
				particleTopology[i][2]=std::stof(temp);
			}
			if(j>3){
				int connection = std::stoi(temp);
				if(body.tellg()==-1) throw oxDNAException("Missing color after connection");
				body>>temp;
				number bfact = std::stod(temp);
				if(body.tellg()==-1) throw oxDNAException("Missing r0 after color");
				body>>temp;
				number tempro = std::stod(temp);
				p->add_neighbour(particles[connection],bfact,tempro);
			}
			j++;
		}
		i++;
	};
	// coordinateConverter();

	//Setting up the matrix
	#pragma omp parallel for
	for(i=0;i<totPar;i++){
		auto *p = static_cast<CCGParticle *>(particles[i]);
		connections[i][0]=p->spring_neighbours.size();
		r0[i][0]=0; // all the eq radius are different
		k0[i][0]=0; // all the spring constant are different
		for(j=0;j<(int)p->spring_neighbours.size();j++){
			connections[i][j+1]=p->spring_neighbours[j];
			r0[i][j+1]=(float)p->ro[j];
			k0[i][j+1]=(float)p->Bfactor[j];
		};
	}
	setRcut(particles);
	std::cout<<"Finished reading topology file"<<std::endl;
	// Debugging
	// print2DArraytoFile("connections.ign",(int*)connections,totPar,PSPmaxNeighbour);
	// print2DArraytoFile("r0.ign",(float*)r0,totPar,PSPmaxNeighbour);
	// print2DArraytoFile("k0.ign",(float*)k0,totPar,PSPmaxNeighbour);
	// print2DArraytoFile("particleTopology.ign",(float*)particleTopology,totPar,3);
	// print2DArraytoFile("patches.ign",(float*)patches,PSPmaxPatchColor,5);
	// print2DArraytoFile("particlePatches.ign",(int*)particlePatches,PSPmaxParticleColor,PSPmaxPatchColor);
}

number PSPInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	if(pCG->has_bond(q)){
		return pair_interaction_bonded(p,q,compute_r,update_forces);
	}else{
		// this->strength=pCG->strength;
		return pair_interaction_nonbonded(p,q,compute_r,update_forces);
	}
}

number PSPInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	number energy=0;
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	auto *qCG = static_cast<CCGParticle*>(q);
	// std::cout<< "this is called "<< std::endl;
	// if(pCG->has_bond(q)){
		energy += (harmonics)?spring(pCG,qCG,compute_r,update_forces):fene(pCG,qCG,compute_r,update_forces);
		if(helixBubble){
			energy += bonded_twist(pCG, qCG, false, update_forces);
			energy += bonded_double_bending(pCG, qCG, false, update_forces);
			energy += bonded_alignment(pCG, qCG, false, update_forces);
		}
	// }
	return energy;
	return 0;
}

number PSPInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	// Non-bonded particle interaction
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	auto *qCG = static_cast<CCGParticle*>(q);
	number energy=0;
	// cout<< "this is called "<< std::endl;
	if(!pCG->has_bond(q)){
		// cout<<"Radius : " << pCG->radius << " " << qCG->radius << endl;
		energy+=exc_vol_nonbonded(pCG,qCG,compute_r,update_forces);
		energy+=patchyInteractionSimple(pCG,qCG,compute_r,update_forces);
	}
	return energy;
}

int PSPInteraction::bonded(int p, int q){
	for(i=1;i<=connections[p][0];i++){
		if(connections[p][i]==q) return i;
	}
	return -1;
}

number PSPInteraction::spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	// updating _computed_r is necessary.
	if(compute_r){ 
		_computed_r = _box->min_image(p->pos,q->pos);
		rmod=_computed_r.module();// calculate the minimum distance between p and q in preodic boundary
	}

	int index = bonded(p->index,q->index);
	if(index==-1) return 0; // not bonded

	number k = k0[p->index][index];
	number r = r0[p->index][index];
	number dist=rmod-r; // Distance between the particles - equilibrium distance
	// number energy = 0.25*k*SQR(dist);
	if (update_forces){
		LR_vector force = ((-k*dist))*(_computed_r/rmod);
		p->force-=force;
		q->force+=force;
	}
	return 0.25*k*SQR(dist);
}

number PSPInteraction::exc_vol_nonbonded(CCGParticle *p, CCGParticle *q, bool compute_r,bool update_forces){
	if(compute_r) _computed_r = _box->min_image(p->pos,q->pos);
	// if(p->index==0 && q->index==1) std::cout<< "Distance between particles = "<<_computed_r <<std::endl; ;
	number totalRadius = p->radius+q->radius;
	number sigma=totalRadius*patchySigma;
	number rstar=totalRadius*patchyRstar;
	number b = patchyB/SQR(totalRadius);
	number rc = patchyRc*totalRadius;
	number energy =0;
	LR_vector force;
	energy = linearRepulsion(patchyEpsilon,_computed_r,force,sigma,rstar,b,rc,update_forces);
	// energy = hardRepulsive(patchyEpsilon,_computed_r,force,sigma,rc,update_forces);
	if(update_forces){
		p->force-=force;
		q->force+=force;
	}
	// std::cout<<24.0*patchyEpsilon * (2*SQR(lj_part)-lj_part) / rnorm<<"\t"<<rnorm<<"\t"<<p->index<<"\t"<<q->index<<"\n";
	return energy;
}

number PSPInteraction::fene(CCGParticle *p, CCGParticle *q, bool compute_r,bool update_forces){ // For now this function is not being used 
	if(compute_r){ 
		_computed_r = _box->min_image(p->pos,q->pos);
		rmod=_computed_r.module();// calculate the minimum distance between p and q in preodic boundary
	}
	number k,r0;
	p->return_kro(q->index,&k,&r0);
	r0=rmod-r0; // Distance between the particles - equilibrium distance
	number energy = -0.5*k*log(1 - SQR(r0) / tepFeneDelta2);
	if(fabs(r0) > tepFeneDelta - DBL_EPSILON) {
		if(update_forces) {
			throw oxDNAException("(TEPInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(r0));
		}
		return (number) (1.e12); //Explode the simulation
	}
	number totalRadius = p->radius+q->radius;
	number sigma=totalRadius*patchySigma;
	number rstar=totalRadius*patchyRstar;
	number b = patchyB/SQR(totalRadius);
	number rc = patchyRc*totalRadius;
	LR_vector force;
	energy+=linearRepulsion(patchyEpsilon,_computed_r,force,sigma,rstar,b,rc,update_forces);
	if(update_forces) {
		force += _computed_r * (-(k * r0 / (tepFeneDelta2 - r0 * r0)) / rmod);
		p->force -= force;
		q->force += force;
	}
	return energy;
}

number PSPInteraction::bonded_double_bending(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces) {
	number energy = 0;
	LR_vector torque(0., 0., 0.);
	if(p->spring_neighbours.size()<2 || q->spring_neighbours.size()<2)return 0.; //return 0 for extreme top and bottom particles

	LR_vector &up = p->orientationT.v1;
	LR_vector &uq = q->orientationT.v1;
	LR_vector &vq = q->orientationT.v2;
	number &th_b_0 = p->th_b_0;
	number &beta_0 = p->beta_b_0;
	// the particle behind sees the particle ahead as if it were slightly tilted,
	// and tries to align with the tilted particle
	LR_vector uqq = rotateVectorAroundVersor(uq, vq, th_b_0);
	uqq = rotateVectorAroundVersor(uqq, uq, beta_0);
	number cosine = up * uqq;

	LR_vector sin_vector = up.cross(uqq);

	energy = (number) 0.f;

	const number tu = p->xu_bending;
	const number tk = p->xk_bending;
	const number kb1 = p->kb1;
	const number kb2 = p->kb2;

	number gu = cos(tu);
	number angle = acos(cosine);

	number A = (kb2 * sin(tk) - kb1 * sin(tu)) / (tk - tu);
	number B = kb1 * sin(tu);
	number C = kb1 * (1 - gu) - (A * SQR(tu) / 2. + tu * (B - tu * A));
	number D = cos(tk) + (A * SQR(tk) * 0.5 + tk * (B - tu * A) + C) / kb2;

	number g1 = (angle - tu) * A + B;
	number g_ = A * SQR(angle) / 2. + angle * (B - tu * A);
	number g = g_ + C;

// 	// unkinked bending regime
	if(acos(cosine) <tu) {
		if(update_forces) {
			torque = -_kb * kb1 * sin_vector;
			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
// 			// forces are not updated since this term of the potential only introduces a torque,
// 			// since it does not depend on r.
		}
		energy = _kb * (1 - cosine) * kb1;
	}		//intermediate regime
	else if(acos(cosine) < tk) {
		if(update_forces) {
			sin_vector.normalize();
			torque = -_kb * g1 * sin_vector;

			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
		}
		energy = _kb * g;
	}		// kinked bending regime - same as unkinked, but with an additive term to the energy and with a different bending.
	else {
		if(update_forces) {
			torque = -_kb * kb2 * sin_vector;
			p->torque -= p->orientationT * torque;
			q->torque += q->orientationT * torque;
// 			// forces are not updated since this term of the potential only introduces a torque,
// 			// since it does not depend on r.
		}
		// energy = _kb*( (1 - up*uq) * kb2 + _get_phi_bending(p->index) );
		energy = _kb * (D - cosine) * kb2;

	}
	return energy;

}

number PSPInteraction::bonded_twist(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces) {

	if(_kt == 0 || p->kt_pref == 0) return 0;// just return 0 if the k_t is 0.
//	return the twist part of the interaction energy.
	LR_vector &up = p->orientationT.v1;
	LR_vector &uq = q->orientationT.v1;
	LR_vector &fp = p->orientationT.v2;
	LR_vector &fq = q->orientationT.v2;
	LR_vector &vp = p->orientationT.v3;
	LR_vector &vq = q->orientationT.v3;

	LR_vector torque(0., 0., 0.);
	number energy = 0;
	number M = fp * fq + vp * vq;
	number L = 1 + up * uq;
	number cos_alpha_plus_gamma = M / L;

	// if it gets here, the twisting angle is too big
	if(-cos_alpha_plus_gamma >= _twist_b) {
		if(p->spring_neighbours.size()==2 || q->spring_neighbours.size()==2){
			// if it doesn't involve the terminal bead and we're using MD, exit with an error
			if(update_forces) {
				throw oxDNAException("(TEPInteraction.cpp) During the simulation, the twisting angle between bonded neighbors %d and %d exceeded acceptable values (cos_alpha_plus_gamma = %g)", p->index, q->index, cos_alpha_plus_gamma);
			}
			// if the forces are not needed, just log it.
			// OX_LOG(Logger::LOG_INFO,"the bond between the bead %d and the bead %d is too twisted! cos_alpha_plus_gamma = %g, average superhelical density = %g/N,_my_time1 = %lld,_my_time2 = %lld",p->get_index(),q->get_index(),cos_alpha_plus_gamma,_my_time1*(_o1_modulus/(2*PI)-_o2_modulus/(2*PI)),_my_time1, _my_time2);
		}
		energy = 1.e12;
	}
	else {
// 		// if it gets here, the angle is in the normal region
		if(-cos_alpha_plus_gamma <= _twist_a) {
			if(update_forces) {

				torque = -(_kt / L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq)) * p->kt_pref;

				p->torque -= p->orientationT * torque;
				q->torque += q->orientationT * torque;
				// forces are not updated since this term of the potential only introduces a torque,
				// as it only depends on r.
			}
			energy = _kt * (1 - cos_alpha_plus_gamma);
		}
		// if it gets here, the angle is in the intermediate region, where the behaviour is lorentzian
		else {
			// quadratic constants of the potential - here for sentimental reasons
			//number A = ( _twist_b - _twist_a)*( _twist_b - _twist_a)*( _twist_b - _twist_a )*0.5;	//A = (b - a)^3/2
			//number C = 0.5*(_twist_a - _twist_b) + _twist_a;// C = (a - b)/2 + a

			number A = (_twist_b - _twist_a) * (_twist_b - _twist_a) * (_twist_b - _twist_a) * (_twist_b - _twist_a) * (_twist_b - _twist_a) * 0.25;			//A = (b - a)^5/4
			number C = 0.25 * (_twist_a - _twist_b) + _twist_a;			// C = (a - b)/4 + a

			if(update_forces) {
				// the torque is made steeper by multiplication of the derivative of the potential. It should work.
				//torque.normalize();
				torque = -(_kt / L) * (fp.cross(fq) + vp.cross(vq) - cos_alpha_plus_gamma * up.cross(uq)) * p->kt_pref;				//this torque is the same as above
				// quadratic prefactor - here for sentimental reasons
				//torque *= 2*A/( (SQR(_twist_b+cos_alpha_plus_gamma))*(_twist_b+cos_alpha_plus_gamma) );
				torque *= 4 * A / ((SQR(_twist_b + cos_alpha_plus_gamma)) * (SQR(_twist_b + cos_alpha_plus_gamma)) * (_twist_b + cos_alpha_plus_gamma));

				p->torque -= p->orientationT * torque;
				q->torque += q->orientationT * torque;
			}
			energy = _kt * (1. + A / (SQR( -cos_alpha_plus_gamma - _twist_b) * SQR(-cos_alpha_plus_gamma - _twist_b)) + C);
		}
	}

	energy *= p->kt_pref;
	return energy;
}

number PSPInteraction::bonded_alignment(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces) {
//	return the alignment term of the interaction energy.

	LR_vector up, tp;
	// these particles are initialised to P_VIRTUAL to prevent gcc from complaining
	BaseParticle *backp = P_VIRTUAL, *frontp = P_VIRTUAL;
// //	make sure the particle q follows p and not the contrary

	if(q->index>p->index) {
		up = p->orientationT.v1;
		tp = q->pos - p->pos;
		backp = p;
		frontp = q;
	} //now this only happens for the last particle, so this means that the particle p has to be aligned with the vector between r_p - r_q
// 	else if(p == q->n5) { // subho my model already takes care of the problem
// 		//new bit
// 		up = p->orientationT.v1;
// 		tp = p->pos - q->pos;
// 		backp = p;
// 		frontp = q;

// 		//old bit - kept here for the records
// 		//up = q->orientationT.v1;
// 		//tp = p->pos - q->pos;
// 		//backp = q;
// 		//frontp = p;
// 		number tpm = tp.module();
// 		if(update_forces) {
// 			// prima che inizi e' 
// 			// LR_vector force = _ka*(up - tp*(up*tp)/SQR(tpm))/tpm;
// 			LR_vector force = -_ka * (up - tp * (up * tp) / SQR(tpm)) / tpm;
// 			backp->force -= force;
// 			frontp->force += force;
// 			//only the torque on p is updated, since this interaction term is basically a self-interaction
// 			//that keeps a particle's u vector aligned with  its tangent vector.
// 			// prima che inizi e' 
// 			//backp->torque -= backp->orientationT*((_ka*tp.cross(up))/tpm);
// 			backp->torque -= backp->orientationT * ((_ka * tp.cross(up)) / tpm);
// 			//LR_vector temp=((_ka*tp.cross(up))/tp.module());
// 			//printf("%lf %lf %lf %lf %lf %lf\n",temp.x,temp.y,temp.z, tp.cross(force).x, tp.cross(force).y,tp.cross(force).z);
// 		}
// 		return _ka * (1 - (up * tp) / tpm);
// 	}
// 	else {
// 		printf("CHe cazzo ci faccio qua io?\n");
// 		abort();
// 	}
	number tpm = tp.module();
	if(update_forces) {
		LR_vector force = _ka * (up - tp * (up * tp) / SQR(tpm)) / tpm;
		backp->force -= force;
		frontp->force += force;
		//only the torque on p is updated, since this interaction term is basically a self-interaction
		//that keeps a particle's u vector aligned with  its tangent vector.
		backp->torque -= backp->orientationT * ((_ka * tp.cross(up)) / tpm);
	}

	return _ka * (1 - (up * tp) / tpm);
}

// Rotate a vector around a versor (TODO: propose to add this to defs.h)
LR_vector PSPInteraction::rotateVectorAroundVersor(const LR_vector vector, const LR_vector versor, const number angle) {
	/* According to section 5.2 of this webpage http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ ,
	 // the image of rotating a vector (x,y,z) around a versor (u,v,w) is given by
	 //      / u(ux + vy + wz)(1 - cos th) + x cos th + (- wy + vz) sin th \
//      | v(ux + vy + wz)(1 - cos th) + y cos th + (+ wx - uz) sin th |
	 //      \ w(ux + vy + wz)(1 - cos th) + z cos th + (- vx + uy) sin th /
	 //      Note that the first parenthesis in every component contains the scalar product of the vectors,
	 //      and the last parenthesys in every component contains a component of the cross product versor ^ vector.
	 //      The derivation that they show on the webpage seems sound, and I've checked it with Mathematica in about//			 1 hour of ininterrupting swearing. */
	number costh = cos(angle);
	number sinth = sin(angle);
	number scalar = vector * versor;
	LR_vector cross = versor.cross(vector);
	return LR_vector(versor.x * scalar * (1. - costh) + vector.x * costh + cross.x * sinth, versor.y * scalar * (1. - costh) + vector.y * costh + cross.y * sinth, versor.z * scalar * (1. - costh) + vector.z * costh + cross.z * sinth);
}

number PSPInteraction::linearRepulsion(number patchyEpsilon, LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = patchyEpsilon * b * SQR(rrc);
			if(update_forces) force = -r * (2 * patchyEpsilon * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * patchyEpsilon * (lj_part - 2*SQR(lj_part)) /rnorm);
			// std::cout<<"Inner collision"<<std::endl;

		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}

number PSPInteraction::cubicRepulsion(number patchyEpsilon, LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces){
	rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = patchyEpsilon * b * SQR(SQR(rrc));
			if(update_forces) force = -r * (8 * patchyEpsilon * b * CUB(rrc) / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * EXCL_EPS * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * patchyEpsilon * (lj_part - 2*SQR(lj_part)) /rnorm);
			// std::cout<<"Inner collision"<<std::endl;

		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}


// Patchy Functions
number PSPInteraction::patchyInteractionSimple(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces){
	if(compute_r){
		_computed_r = _box->min_image(p->pos,q->pos);
		rnorm = _computed_r.norm();
	}
	if(p->btype==100||q->btype==100) return 0; // no color present ignore
	if(p->strand_id>=0 && p->strand_id==q->strand_id) return 0;
	// if(rnorm<patchyRcut) return 0;
	if(_computed_r.module()>patchyRcut+p->radius+q->radius) return 0;// if the particles are farther than the cutoff, ignore
	number energy=0;
	// int c = 0;
	LR_vector tmptorquep(0, 0, 0);
	LR_vector tmptorqueq(0, 0, 0);
	int pparticleColor = particleTopology[p->index][1];
	int qparticleColor = particleTopology[q->index][1];
	// LR_vector pforce(0, 0, 0),qforce(0, 0, 0),ptorque(0, 0, 0),qtorque(0, 0, 0);
	#pragma omp parallel for collapse(2) reduction(+:energy)
	for(int pi=0;pi<particlePatches[pparticleColor][0];pi++){
		for(int qi=0;qi<particlePatches[qparticleColor][0];qi++){
			int pPatch = particlePatches[pparticleColor][pi+1];
			int qPatch = particlePatches[qparticleColor][qi+1];

			if(patches[pPatch][0]+patches[qPatch][0]!=0) continue; // colors are not compilmentory
			if(patchLock[p->index][pi]|patchLock[q->index][qi]) continue; // patch is already locked
			
			LR_vector pPatchR = LR_vector(patches[pPatch][2],patches[pPatch][3],patches[pPatch][4]);
			LR_vector qPatchR = LR_vector(patches[qPatch][2],patches[qPatch][3],patches[qPatch][4]);

			LR_vector patchDist = _computed_r + qPatchR - pPatchR;
			number patchDist2 = patchDist.norm();
			if(patchDist2>patchyRcut2) continue; // patch is too far


			number r8b10 = patchDist2*patchDist2*patchDist2*patchDist2 / patchyPowAlpha;
			number exp_part = -1.001f * exp(-(number)0.5f * r8b10 * patchDist2);
			number energyIJ =  (patches[pPatch][1]+patches[qPatch][1]) * (exp_part - patchyEcut);

			energy += energyIJ;

			if(update_forces){
				// if(energyIJ<patchyLockCutOff){
				// 	patchLock[p->index][pi]=true;
				// 	patchLock[q->index][qi]=true;
				// }else{
				// 	patchLock[p->index][pi]=false;
				// 	patchLock[q->index][qi]=false;
				// }
				number forceMag = 5*exp_part*r8b10;
				LR_vector force = patchDist*(forceMag);
				LR_vector ptorque = p->orientationT*pPatchR.cross(force), qtorque = q->orientationT*qPatchR.cross(force);
				#pragma omp critical // Do this one by one, overwriting can be 
				{
				p->force-=force;
				q->force+=force;
				p->torque-=ptorque;
				q->torque+=qtorque;
				}
			}
		}
	}
	// for(int pi=0;pi)
	return energy;
}

number PSPInteraction::patchyInteraction2point(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces){
	number energy = (number) 0.f;
	number sqr_r = _computed_r.norm();
	if(sqr_r > _sqr_rcut) return (number) 0.f;
	for(uint pi = 0; pi < p->N_int_centers(); pi++) {
		LR_vector ppatch = p->int_centers[pi];
		for(uint pj = 0; pj < q->N_int_centers(); pj++) {
			LR_vector qpatch = q->int_centers[pj];

			LR_vector patch_dist = _computed_r + qpatch - ppatch;
			number dist = patch_dist.norm();
			if(dist < _sqr_patch_rcut) {
				number r_p = sqrt(dist);
				number exp_part = exp(_sigma_ss / (r_p - _rcut_ss));
				number tmp_energy = patchyEpsilon*_A_part * exp_part * (_B_part / SQR(dist) - 1.);

				energy += tmp_energy;

				number tb_energy = (r_p < _sigma_ss) ? 1 : -tmp_energy;

				PatchyBond p_bond(q, r_p, pi, pj, tb_energy);
				PatchyBond q_bond(p, r_p, pj, pi, tb_energy);

				if(update_forces) {
					number force_mod = patchyEpsilon * _A_part * exp_part * (4. * _B_part / (SQR(dist) * r_p)) + _sigma_ss * tmp_energy / SQR(r_p - _rcut_ss);
					LR_vector tmp_force = patch_dist * (force_mod / r_p);

					LR_vector p_torque = p->orientationT * ppatch.cross(tmp_force);
					LR_vector q_torque = q->orientationT * qpatch.cross(tmp_force);

					p->force -= tmp_force;
					q->force += tmp_force;

					p->torque -= p_torque;
					q->torque += q_torque;

					p_bond.force = tmp_force;
					p_bond.p_torque = p_torque;
					p_bond.q_torque = q_torque;

					q_bond.force = -tmp_force;
					q_bond.p_torque = -q_torque;
					q_bond.q_torque = -p_torque;
				}

				_particle_bonds(p).emplace_back(p_bond);
				_particle_bonds(q).emplace_back(q_bond);

				energy += patchyInteraction3point(p, p_bond, update_forces);
				energy += patchyInteraction3point(q, q_bond, update_forces);
			}
		}
	}
}

number PSPInteraction::patchyInteraction3point(CCGParticle *p, PatchyBond &new_bond, bool update_forces){
	number energy = 0.;

	number curr_energy = new_bond.energy;
	const auto &p_bonds = _particle_bonds(p);
	for(auto &other_bond : p_bonds) {
		if(other_bond.other != new_bond.other && other_bond.p_patch == new_bond.p_patch) {
			number other_energy = other_bond.energy;
			energy += _lambda * curr_energy * other_energy;

			if(update_forces) {
				if(new_bond.r_p > _sigma_ss) {
					BaseParticle *other = new_bond.other;
					number factor = -_lambda * other_energy;
					LR_vector tmp_force = factor * new_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * new_bond.p_torque;
					other->torque += factor * new_bond.q_torque;
				}
				if(other_bond.r_p > _sigma_ss) {
					BaseParticle *other = other_bond.other;

					number factor = -_lambda * curr_energy;
					LR_vector tmp_force = factor * other_bond.force;

					p->force -= tmp_force;
					other->force += tmp_force;

					p->torque -= factor * other_bond.p_torque;
					other->torque += factor * other_bond.q_torque;
				}
			}
		}
	}

	return energy;
}

number PSPInteraction::patchyInteractionColored(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces){
	return 0;
}

number PSPInteraction::patchyInteractionBubble(CCGParticle *p, CCGParticle *q, bool compute_r, bool update_forces){
	return 0;
}

void PSPInteraction::setRcut(std::vector<BaseParticle *> &particles){
	if(_rcut<0){
		// #pragma omp parallel for
		for(i=0;i<totPar;i++){
			auto *p = static_cast<CCGParticle*>(particles[i]);
			if(p->radius>_rcut) _rcut=p->radius;
		}
		_rcut+=patchyRcut;
		// _rcut+=patchyIntercept+patchySpacer;
		_rcut*=2.1f;
		// if(_rcut<patchyRcut) _rcut=patchyRcut; // This should not be the case
	}
	std::cout<<"Rcut is set to "<<_rcut<<std::endl;
	_sqr_rcut=SQR(_rcut);
}

//Debug function
void PSPInteraction::print2DArraytoFile(std::string filename, float* arr, int row, int col){
	std::ofstream file(filename);
	if(!file.good()) throw oxDNAException("Something wrong with the file: %s",filename);
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			file<<arr[i*col+j]<<" ";
		}
		file<<std::endl;
	}
	file.close();
}

void PSPInteraction::print2DArraytoFile(std::string filename, int* arr, int row, int col){
	std::ofstream file(filename);
	if(!file.good()) throw oxDNAException("Something wrong with the file: %s",filename);
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			file<<arr[i*col+j]<<" ";
		}
		file<<std::endl;
	}
	file.close();
}