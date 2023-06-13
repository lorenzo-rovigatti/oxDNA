#include "CCGInteraction.h"
#include "../Particles/CCGParticle.h"

CCGInteraction::CCGInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(SPRING,spring);
	ADD_INTERACTION_TO_MAP(EXEVOL,exc_vol);
	ADD_INTERACTION_TO_MAP(PATCHY,patchy_interaction);
}

CCGInteraction::~CCGInteraction() {

}

// number CCGInteraction::ccg_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
// 	return 0.f;
// }

void CCGInteraction::get_settings(input_file &inp) {
	// OX_LOG(Logger::LOG_INFO,"Setting is being called");
	BaseInteraction::get_settings(inp);
	
	if(getInputString(&inp,"patchyAlpha",temp,0)==KEY_FOUND){
		patchyAlpha=stod(temp);
		OX_LOG(Logger::LOG_INFO,"New alpha value for patchy interaction = %d",patchyAlpha);
	}
	if(getInputString(&inp,"patchyCutoff",temp,0)==KEY_FOUND){
		patchyCutoff=stod(temp);
		OX_LOG(Logger::LOG_INFO,"New cutoff value for the patchy interaction = %d",patchyCutoff);
	}
	if(getInputString(&inp,"patchyRcut",temp,0)==KEY_FOUND){
		patchyRcut =stod(temp);
		OX_LOG(Logger::LOG_INFO,"New cutoff value for the patchy interaction = %d",patchyRcut);
	}
	if(getInputString(&inp,"rcut",temp,0)==KEY_FOUND){
		_rcut=stod(temp);
		_sqr_rcut = SQR(_rcut);
		OX_LOG(Logger::LOG_INFO,"New interaction radius cutoff = %d",_rcut);
	}else{
		_rcut= 1.2;
		_sqr_rcut=SQR(_rcut);
	}
}

void CCGInteraction::init() {
    OX_LOG(Logger::LOG_INFO,"CCG is initilized");
}

void CCGInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	OX_LOG(Logger::LOG_INFO,"Filling up the void with CCGParticles");
	particles.resize(totPar);//Confirming the particles vector has size to accomodate all the particles.
	for(i=0;i<totPar;i++){
		particles[i] = new CCGParticle();
		particles[i]->index=i;
		particles[i]->type=0;
	}
}

number CCGInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	if(pCG->has_bond(q)){
		return pair_interaction_bonded(p,q,compute_r,update_forces);
	}else{
		this->strength=pCG->strength;
		return pair_interaction_nonbonded(p,q,compute_r,update_forces);
	}
}

number CCGInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	// OX_DEBUG("Bonded interaction is called");
	if(update_forces){
		_computed_r = _box->min_image(p->pos,q->pos);
	}
	r=_computed_r; //remove this later
	this->rmod=_computed_r.module();
	number energy = spring(p,q,compute_r,update_forces);
	// number energy = debug(p,q,compute_r,update_forces);
	// energy+=exc_vol(p,q,compute_r,update_forces);
	return energy;
}

number CCGInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	// OX_DEBUG("Nonbonded interaction is called");
	// number energy = exc_vol(p,q,compute_r,update_forces);
	// energy+= patchy_interaction(p,q,compute_r,update_forces);
	return 0;
}

number CCGInteraction::spring(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	// r=_computed_r;
	// rmod=r.module();
	// std::cout <<p->index<<"\t"<<q->index<<"\t"<<p->pos<<"\t"<<q->pos<<"\t"<<_computed_r<<"\t"<<_box->min_image(p,q)<<std::endl;
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	double k,r0;
	pCG->return_kro(q->index,&k,&r0);
	// std::cout<< p->index <<"\t"<<q->index<<"\t"<<k<<"\t"<<r0<<"\n";
	double dist=rmod-r0; // Distance between the particles - equilibrium distance
	double energy = 0.5*k*SQR(dist); // Energy = 1/2*k*x^2
	// OX_DEBUG("For particle %d and %d Stiffness = %d and r0 = %d",p->index,q->index,k,r0);
	if(update_forces){
		LR_vector force = (-k*_computed_r*dist/rmod); //force = -k*(r_unit*dist) =-k(r-r0)
		// std::cout<< p->index << "\t"<<q->index<<"\t k="<<k<<"\t dist="<<dist<<"\t rmod = "<<rmod<<"\t force = "<<force.x<<","<<force.y<<","<<force.z<<"\n";
		p->force-= force;//substract force from p
		q->force+= force;//add force to q
	}
	return energy;
}

number CCGInteraction::debug(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	// pCG->return_bfactor(q->index);
	double k,r0;
	pCG->return_kro(q->index,&k,&r0);
	std::cout<< p->index <<"\t"<<q->index<<"\t"<<k<<"\t"<<r0<<"\n";
	return 0.f;
}

number CCGInteraction::patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	number energy =0.f;
	if(color_compatibility(p,q)){
		// double r2 = SQR(r.x)+SQR(r.y)+SQR(r.z); //r^2
		if(rmod<patchyCutoff){
			double r8b10 = pow(rmod,8)/pow(patchyAlpha,10);
			double expPart = -1.001*exp(-0.5*r8b10*rmod*rmod);
			energy=strength*(expPart-patchyCutoff); //patchy interaction potential

			if(update_forces){
				double f1D = 5 * expPart * r8b10;
				LR_vector force = r*f1D;
				p->force+=force;
				q->force-=force;
			}
		}
		//Locking algo; don't forgot to look into color compatibility and particle definition as well
		// auto *pCG = static_cast<CCGParticle*>(p);
		// auto *qCG = static_cast<CCGParticle*>(q);

		// if(update_forces && !pCG->multipatch && !qCG->multipatch){
		// 	if(energy<lockCutOff){
		// 		pCG->lockedTo = qCG->index;
		// 		qCG->lockedTo = pCG->index;
		// 	}
		// 	else{
		// 		pCG->lockedTo=0;
		// 		qCG->lockedTo=0;
		// 	}
		// }
	}
	return energy;
}

double CCGInteraction::exc_vol(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	auto *pCCG = static_cast<CCGParticle*>(p);
	auto *qCCG = static_cast<CCGParticle*>(q);
	double totalRadius = pCCG->radius+qCCG->radius;
	double sig= sigma*totalRadius;
	double rs = rstar*totalRadius;
	double bs = b/totalRadius;
	double rcs = rc*totalRadius;
	// double rmod = r.module();
	double energy = 0.f;
	if(rmod<rs){
		double factor = SQR(sig/rmod);
		double tmp = factor*factor*factor;
		energy=4*epsilon*(SQR(tmp)-tmp);
		if(update_forces){
			LR_vector force = -_computed_r*(24*epsilon*(tmp-2*SQR(tmp))/rmod);
			p->force-=force;
			q->force+=force;
		}
	}else if (rmod<rcs){
		double rrc = rmod-rcs;
		energy=epsilon*b*SQR(rrc);
		if(update_forces){
			LR_vector force = -_computed_r*(2*epsilon*bs*rrc/rmod);
			p->force-=force;
			q->force+=force;
		}
	}

	return energy;
}

bool CCGInteraction::color_compatibility(BaseParticle *p, BaseParticle *q){
	// auto *pCCG = static_cast<CCGParticle*>(p);
	if(this->strength==0) return false; //Prevent any further calculations if strength of p =0
	if(abs(p->btype)==100 || abs(q->btype)==100) return false; //100 and -100 are colorless particles
	// if(pCCG->lockedTo==q->btype) return true;
	// if(pCCG->lockedTo!=0) return false;
	if(abs(p->btype)<10 && abs(q->btype)<10){ //self-complementary domain 1 interact with 1 and so on till 9
		if(p->btype==q->btype){
			return true;
		}else{
			return false;
		}
	}else if(p->btype+q->btype ==0){ //only -5 interact with +5 
		return true;
	}else{
		return false;
	}
}

void CCGInteraction::read_topology(int *N_strands, std::vector<BaseParticle*> &particles) {
	OX_LOG(Logger::LOG_INFO,"Read Topology is being called.");
	char line[2048];
	std::ifstream topology;
	topology.open(this-> _topology_filename, std::ios::in);

	// There are other function to take care of this as well
	if(!topology.good())
		throw oxDNAException("Unable to open the topology file: %s. Can you please check the filename again ?",this->_topology_filename);
	
	topology.getline(line,2048); //Read the header from the topology file
	std::stringstream head(line);

	head>> totPar>>strands>>ccg>>ccg0>>noSpring>>noColor; //Saving header info
	if(ccg+ccg0 !=totPar) throw oxDNAException("Number of Colored + Color less particles is not equal to total Particles. Insanity caught in the header.");
	allocate_particles(particles);
	if(head.fail())
		throw oxDNAException("Please check the header, there seems to be something out of ordinary.");
	
	// if(version<currentVersion)
	// 	throw oxDNAException("Seems like the version of the file is an old one. Try looking into the documentation to upgrade it to newer one.");
	i=0,color=0;
	while(topology.getline(line,2048)){ //Read till the end of file.
		if (strlen(line)==0 || line[0]=='#') //Ignore empty line or line starting with # for comment
			continue;
		if(i>totPar) throw oxDNAException("Contains from molecules than Total Particles provided in the header. Mismatch of information between header and body.");
		particles[i]->index=i;
		particles[i]->strand_id=0;
		// particles[i]->type=20;

		// BaseParticle *p = particles[i];
		auto *q = dynamic_cast< CCGParticle *>(particles[i]);
		std::stringstream body(line);
		j=0;
		// connection=true;
		// bcall=true;
		while(body.tellg()!=-1){
			body>>temp;
			if(j==0){
				particleType=std::stoi(temp); // particle type = -2 for future references
			}else if(j==1){
				particles[i]->strand_id=std::stoi(temp); // id of the particle
			}else if(j==2){
				particles[i]->btype=std::stoi(temp); //color of the particle
			}else if(j==3){
				q->radius=std::stoi(temp); //radius of the particle
			}else{
				// if(connection){
				// 	q->add_neighbour(particles[std::stoi(temp)]); //all connected particle
				// 	connection=false;
				// }else if(bcall){
				// 	q->Bfactor.push_back(std::stod(temp)); //connected particle Bfactor
				// 	bcall=false;
				// }else{
				// 	q->ro.push_back(std::stod(temp)); //connected particle equilibrium spring distance r0
				// 	connection=true;
				// 	bcall=true;
				// };
				int connection = std::stoi(temp);
				body>>temp;
				double bfact = std::stod(temp);
				body>>temp;
				double tempro = std::stod(temp);
				q->add_neighbour(particles[connection],bfact,tempro);
			}
			j++;
		}
		// OX_DEBUG("Neighbours for %d is %d and %d",q->index,q->ro[0],q->ro[1]);
		// std::cout << "Spring Neighbours \n";
		// for(int p=0;p<q->spring_neighbours.size();p++){
		// 	std::cout <<q->index<<"\t"<<q->spring_neighbours[p]<<"\n";
		// }

		// std::cout<< "Stiffness \n";
		// for(int p=0;p<q->Bfactor.size();p++){
		// 	std::cout<< q->index<<"\t"<<q->Bfactor[p]<<"\n";
		// }

		// std::cout<<"Mean distance \n";
		// for(int p=0; p<q->ro.size();p++){
		// 	std::cout<<q->index<<"\t"<<q->ro[p]<<"\n";
		// }
		
		// std::cout<< " Color of the particles: " << q->btype<<std::endl;


		particles[i]->btype=color; //btype is used for coloring
		if(i==0)
			OX_LOG(Logger::LOG_INFO, "One loop successfully completed");
		i++;
	}
	*N_strands=strands;
};

void CCGInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}