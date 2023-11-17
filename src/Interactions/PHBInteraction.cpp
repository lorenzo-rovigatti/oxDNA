#include "PHBInteraction.h"
using namespace std;

PHBInteraction::PHBInteraction():BaseInteraction(){

};

PHBInteraction::~PHBInteraction(){

};
void PHBInteraction::get_settings(input_file &inp){
    BaseInteraction::get_settings(inp);
    if(getInputString(&inp,"rcut",temp,0)==KEY_FOUND){ //set the rcut from the input file
		_rcut=stod(temp);
        _sqr_rcut=sqrt(_rcut);
		std::cout<<"The new value of rcut = "<< _rcut<<std::endl;
	}else{
        _rcut=-1; // set the rcut after reading the topology. 2.5x largest radius
    }
};
void PHBInteraction::init(){};
void PHBInteraction::allocate_particles(std::vector<BaseParticle *> &particles){
	particles.resize(totPar);
	for(i=0;i<totPar;i++){
		particles[i] = new PHBParticle();
		particles[i]->index=i;
		particles[i]->type=0;
	}
}; 
void PHBInteraction::check_input_sanity(std::vector<BaseParticle *> &particles){};
number PHBInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    auto *pCG = dynamic_cast<PHBParticle*>(p);
	if(pCG->has_bond(q)){
		return pair_interaction_bonded(p,q,compute_r,update_forces);
	}else{
		// this->strength=pCG->strength;
		return pair_interaction_nonbonded(p,q,compute_r,update_forces);
	}
};
number PHBInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    number energy=0;
	auto *pCG = dynamic_cast<PHBParticle*>(p);
	auto *qCG = static_cast<PHBParticle*>(q);
	if(!pCG->has_bond(q)){
		energy += spring(pCG,qCG,compute_r,update_forces);
		energy += bonded_twist(pCG, qCG, false, update_forces);
		energy += bonded_double_bending(pCG, qCG, false, update_forces);
		energy += bonded_alignment(pCG, qCG, false, update_forces);
	}
	return energy;
};
number PHBInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    auto *pCG = dynamic_cast<PHBParticle*>(p);
	auto *qCG = static_cast<PHBParticle*>(q);
	number energy=0;
	if(pCG->has_bond(q)){
		energy+=exc_vol_nonbonded(p,q,compute_r,update_forces);
		energy+=patchy_interaction_notorsion(pCG,qCG,compute_r,update_forces);
	}
	return energy;
};

double PHBInteraction::exc_vol_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	if(compute_r) _computed_r = _box->min_image(p->pos,q->pos);
	// if(p->index==0 && q->index==1) std::cout<< "Distance between particles = "<<_computed_r <<std::endl; ;
	auto *pCCG = static_cast<PHBParticle*>(p);
	auto *qCCG = static_cast<PHBParticle*>(q);
	double totalRadius = pCCG->radius+qCCG->radius;
	double sigma=totalRadius*patchySigma;
	double rstar=totalRadius*patchyRstar;
	double b = patchyB/SQR(totalRadius);
	double rc = patchyRc*totalRadius;
	double energy =0;
	LR_vector force;
	energy = _repulsive_lj2(patchyEpsilon,_computed_r,force,sigma,rstar,b,rc,update_forces);
	if(update_forces){
		p->force-=force;
		q->force+=force;
	}
	// std::cout<<24.0*patchyEpsilon * (2*SQR(lj_part)-lj_part) / rnorm<<"\t"<<rnorm<<"\t"<<p->index<<"\t"<<q->index<<"\n";
	return energy;
}

number PHBInteraction::_repulsive_lj2(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			rmod = sqrt(rnorm);
			number rrc = rmod - rc;
			energy = prefactor * b * SQR(rrc);
			if(update_forces)
				force = -r * (2 * prefactor * b * rrc / rmod);
		}
		else {
			number tmp = SQR(sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			// the additive term was added by me to mimick Davide's implementation
			energy = 4 * prefactor * (SQR(lj_part) - lj_part) + prefactor;
			if(update_forces)
				force = -r * (24 * prefactor * (lj_part - 2 * SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0)
		force.x = force.y = force.z = (number) 0;

	return energy;
};

number PHBInteraction::fene(PHBParticle *p, PHBParticle *q, bool compute_r,bool update_forces){ // For now this function is not being used 
	if(compute_r){ 
		_computed_r = _box->min_image(p->pos,q->pos);
		rmod=_computed_r.module();// calculate the minimum distance between p and q in preodic boundary
	}
	double k,r0;
	p->return_kro(q->index,&k,&r0);
	r0=rmod-r0; // Distance between the particles - equilibrium distance
	double energy = -0.5*k*log(1 - SQR(r0) / tepFeneDelta2);
	if(fabs(r0) > tepFeneDelta - DBL_EPSILON) {
		if(update_forces && !_allow_broken_fene) {
			throw oxDNAException("(TEPInteraction.cpp) During the simulation, the distance between bonded neighbors %d and %d exceeded acceptable values (d = %lf)", p->index, q->index, fabs(r0));
		}
		return (number) (1.e12); //Explode the simulation
	}
	double totalRadius = p->radius+q->radius;
	double sigma=totalRadius*patchySigma;
	double rstar=totalRadius*patchyRstar;
	double b = patchyB/SQR(totalRadius);
	double rc = patchyRc*totalRadius;
	LR_vector force;
	energy+=_repulsive_lj2(patchyEpsilon,_computed_r,force,sigma,rstar,b,rc,update_forces);
	if(update_forces) {
		force += _computed_r * (-(k * r0 / (tepFeneDelta2 - r0 * r0)) / rmod);
		p->force -= force;
		q->force += force;
	}
	return energy;
}

number PHBInteraction::spring(PHBParticle *p, PHBParticle *q, bool compute_r,bool update_forces){
	// updating _computed_r is necessary.
	if(compute_r){ 
		_computed_r = _box->min_image(p->pos,q->pos);
		rmod=_computed_r.module();// calculate the minimum distance between p and q in preodic boundary
		// std::cout<<rmod<<"\n";
	}
	// _computed_r = _box->min_image(p->pos,q->pos);
	// rmod=_computed_r.module();// calculate the minimum distance between p and q in preodic boundary
	// auto *pCG = static_cast<CCGParticle*>(p);
	double k,r0;
	p->return_kro(q->index,&k,&r0);
	double dist=rmod-r0; // Distance between the particles - equilibrium distance
	double energy = 0.25*k*SQR(dist); // Energy = 1/2*k*x^2 but 1/2 for one particle and 1/2 for other
	// std::cout<<energy<<std::endl;
	if(update_forces){
		LR_vector force = ((-k*dist))*(_computed_r/rmod); //force = -k*(r_unit*dist) =-k(r-r0)
		p->force-= force+damp*p->vel;//substract force from p
		q->force+= force-damp*q->vel;//add force to q
	}
	return energy;
}

number PHBInteraction::bonded_double_bending(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces) {
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

number PHBInteraction::bonded_twist(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces) {

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

number PHBInteraction::bonded_alignment(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces) {
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


// number PHBInteraction::maxRadius(std::vector<PHBParticle*> &particles){
// 	number radius =0;
// 	for(i=0;i<totPar;i++){
// 		if(radius<particles[i]->radius) radius=particles[i]->radius;
// 	}
// 	return radius;
// };

void PHBInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles){
	//Variables that needs to be destroyed 
	std::vector<Patch> patches; //stores the initial patches
    Patch tempPatch; // temporary patch variable;
    std::vector<int> tempColors; 
	std::vector<std::vector<int>> particleColors;// store id corresponding to particle colors with patch info. 1 particleColor can have multiple patches
	
    // char line[2048];
	std::ifstream topology;
	topology.open(this->_topology_filename,std::ios::in);
	if(!topology.good()) throw oxDNAException("Problem with topology file: %s",this->_topology_filename);

	std::getline(topology,temp);
	std::cout<<temp<<std::endl;
	std::stringstream head(temp);
	head>>totPar>>strands; //saving header info
	*N_strands=strands; // This one is important don't forget

	allocate_particles(particles);
	i=0;
	int j=0;
	string line;
	while(std::getline(topology,line)){
		if(line.empty()||line[0]=='#') continue; //skip empty line and comment
		if(line[0]=='i'){
			if(line[1]=='P'){ //these are Patchy informations
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					body>>temp;
					// std::cout<< patches.size()<<std::endl;
					if(j==1 && stoi(temp)!= static_cast<int> (patches.size())) throw oxDNAException("Ids of patches starts with 0 and should be monotonically increasing with gap of 1. We all sometime face counting problems!!!!!!");
					if(j==2) tempPatch.color=stoi(temp);
					if(j==3) tempPatch.strength=stod(temp);
					if(j==4) tempPatch.position.x=stod(temp);
					if(j==5) tempPatch.position.y=stod(temp);
					if(j==6) tempPatch.position.z=stod(temp);
					if(j==7) tempPatch.a1static.x=stod(temp);
					if(j==8) tempPatch.a1static.y=stod(temp);
					if(j==9) tempPatch.a1static.z=stod(temp);
					if(j==10) tempPatch.a2static.x=stod(temp);
					if(j==11) tempPatch.a2static.y=stod(temp);
					if(j==12) tempPatch.a2static.z=stod(temp);
					j++;
				}
				patches.push_back(tempPatch);
				continue;
			}
			if(line[1]=='C'){//these are Color informations
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					body>>temp;
					if(j==1 && stoi(temp)!=static_cast<int> (particleColors.size())) throw oxDNAException("Ids of colors starts with 0 and should be monotonically increasing with gap of 1. We all sometime face counting problems!!!!!!");
					if(j>1){
						tempColors.push_back(stoi(temp));
					}
					j++;
				}
				particleColors.push_back(tempColors);
				continue;
			}
		}
		auto *q = static_cast< PHBParticle *>(particles[i]); //Start working with PHB particles
		std::stringstream body(line);
		j=0;
		while(body.tellg()!=-1){
			body>>temp;
			if(j==0) q->type=std::stoi(temp);
			if(j==1) q->strand_id=std::stoi(temp);
			if(j==2){
				q->btype=std::stoi(temp); // btype is color don't forget
				if(q->btype!=100) //if the color is not 100, it is non empty and do further operations
					for(uint id=0;id<particleColors[q->btype].size();id++)
						q->add_patch(patches[particleColors[q->btype][id]]);
			}
			if(j==3) q->radius=std::stod(temp);
			// if(q->type==-3){ /// These are helix particles
				// std::cout<<"Hexlix Particle"<<std::endl;
			if(j>3){
				int connection = std::stoi(temp);
				if(body.tellg()==-1) throw oxDNAException("Missing color after connection");
				body>>temp;
				double bfact = std::stod(temp);
				if(body.tellg()==-1) throw oxDNAException("Missing r0 after color");
				body>>temp;
				double tempro = std::stod(temp);
				q->add_neighbour(particles[connection],bfact,tempro);
			}
			// }
			j++;
		}
		i++;
	};
	// auto *q = static_cast< PHBParticle *>(particles[0]);
	// cout<<q->patches[0].a1static.x<<endl;
	std::cout<<"Successfully completed topology reading with total types of patches = "<<patches.size()<< " and types of colored particles = "<<particleColors.size()<<std::endl;

};


// Rotate a vector around a versor (TODO: propose to add this to defs.h)
LR_vector PHBInteraction::rotateVectorAroundVersor(const LR_vector vector, const LR_vector versor, const number angle) {
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


///////// Patchy scense //////////////

number PHBInteraction::patchy_interaction_notorsion(PHBParticle *p, PHBParticle *q, bool compute_r, bool update_forces){
	rnorm = _computed_r.norm();
	if(rnorm > this->patchyRcut2) return 0; // not within reach ignore
	if(p->btype==100||q->btype==100) return 0; // no color present ignore
	number energy=0;
	int c = 0;
	LR_vector tmptorquep(0, 0, 0);
	LR_vector tmptorqueq(0, 0, 0);
	for(uint pi=0;pi<p->patches.size();pi++){
		LR_vector ppatch = p->int_centers[pi];
		for(uint qi=0;qi<q->patches.size();qi++){
			// if(bondingAllowed(p->patches[pi],q->patches[qi]))
		}
	}

	return energy;
}

// number PatchyShapeInteraction::_patchy_interaction_notorsion(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
// 	LR_vector r=_computed_r;
// 	if(compute_r){ 
// 		_computed_r = _box->min_image(p->pos,q->pos);
// 		r=_computed_r;
// 	}

// 	number rnorm = r.norm();
// 	if(rnorm > this->_sqr_rcut) return (number) 0.f;

// 	number energy = (number) 0.f;



// 	PatchyShapeParticle *pp = static_cast<PatchyShapeParticle *>(p);
// 	PatchyShapeParticle *qq = static_cast<PatchyShapeParticle *>(q);

// 	int c = 0;
// 	LR_vector tmptorquep(0, 0, 0);
// 	LR_vector tmptorqueq(0, 0, 0);
// 	for(int pi = 0; pi < pp->N_patches; pi++) {
// 		LR_vector ppatch = p->int_centers[pi];

// 		for(int pj = 0; pj < qq->N_patches; pj++) {

// 			if(this->_bonding_allowed(pp,qq,pi,pj)  )
// 			{

// 				number K = pp->patches[pi].strength;
// 			    LR_vector qpatch = q->int_centers[pj];

// 			    LR_vector patch_dist = r + qpatch - ppatch;
// 			    number dist = patch_dist.norm();


// 			    if(dist < SQR(PATCHY_CUTOFF)) {

// 				    c++;
//                     number energy_ij = 0;

// 				    number r8b10 = dist*dist*dist*dist / _patch_pow_alpha;
// 				    number exp_part = -1.001f * exp(-(number)0.5f * r8b10 * dist);


// 				    number f1 =  K * (exp_part - _patch_E_cut);


// 				    energy_ij = f1;// * angular_part;
//                     energy += energy_ij;

//                     if(update_forces && this->_no_multipatch)
//                     {
//                      if (energy_ij < this->_lock_cutoff )
//                      {
//                     	qq->patches[pj].set_lock(p->index,pi,energy_ij);
//                         pp->patches[pi].set_lock(q->index,pj,energy_ij);
//                      }
//                      else
//                      {
//                     	qq->patches[pj].unlock();
//                     	pp->patches[pi].unlock();

//                      }

//                     }

// 				if(update_forces ) {
// 					number f1D =  (5 * exp_part * r8b10);
// 					LR_vector tmp_force = patch_dist * (f1D ); //patch_dist * (f1D * angular_part);
// 					LR_vector torqueq(0,0,0) ; //= dir;
//                     LR_vector torquep(0,0,0) ; //= dir;
// 					torquep += ppatch.cross(tmp_force);
// 					torqueq += qpatch.cross(tmp_force);


// 					p->torque -= p->orientationT * torquep;
// 					q->torque += q->orientationT * torqueq;

// 					p->force -= tmp_force;
// 					q->force += tmp_force;
// 				}
// 			   }
// 			}
// 		}
// 	}
// 	return energy;
// }