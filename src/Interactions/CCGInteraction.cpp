#include "CCGInteraction.h"
#include "../Particles/CCGParticle.h"

CCGInteraction::CCGInteraction() :
				BaseInteraction() {
	// ADD_INTERACTION_TO_MAP(CCG, _ccg_interaction);
}

CCGInteraction::~CCGInteraction() {

}

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
}

void CCGInteraction::init() {
    OX_LOG(Logger::LOG_INFO,"CCG is initilized");
}

void CCGInteraction::allocate_particles(std::vector<BaseParticle*> &particles) {
	OX_LOG(Logger::LOG_INFO,"Filling up the void with CCGParticles");
	particles.resize(totPar);//Confirming the particles vector has size to accomodate all the particles.
	for(i=0;i<totPar;i++){
		particles[i] = new CCGParticle();
	}
}

number CCGInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	this->r=_computed_r;
	this->rmod=r.module();
	if(pCG->has_bond(q)){
		return pair_interaction_bonded(p,q,compute_r,update_forces);
	}else{
		this->strength=pCG->strength;
		return pair_interaction_nonbonded(p,q,compute_r,update_forces);
	}
}

number CCGInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = spring(p,q,compute_r,update_forces);
	energy+=exc_vol(p,q,compute_r,update_forces);
	return energy;
}

number CCGInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = exc_vol(p,q,compute_r,update_forces);
	return energy;
}

number CCGInteraction::spring(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	double k,r0;
	pCG->return_kro(q->index,&k,&r0);
	number dist=rmod-r0; // Distance between the particles - equilibrium distance
	number energy = 0.5*k*SQR(dist); // Energy = 1/2*k*x^2
	if(update_forces){
		LR_vector force = (-1.f*k*_computed_r*dist/rmod); //force = -k*(r_unit*dist) =-k(r-r0)
		p->force-= force;//substract force from p
		q->force+= force;//add force to q
	}
	return energy;
}

number CCGInteraction::patchy_interaction(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	number energy =0.f;
	if(color_compatibility(p,q)){
		double r2 = SQR(r.x)+SQR(r.y)+SQR(r.z);
		energy=strength*(-1.001*exp(-pow(r2,5)/pow(patchyAlpha,10))-patchyCutoff);
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
	if(this->strength==0) return false; //Prevent any further calculations if strength of p =0
	if(abs(p->btype)==100 || abs(q->btype)==100) return false; //100 and -100 are colorless particles 

	if(abs(p->btype)<10 && abs(q->type)<10){ //self-complementary domain 1 interact with 1 and so on till 9
		if(p->type==q->type){
			return true;
		}else{
			return false;
		}
	}else if(p->btype+q->type ==0){ //only -5 interact with +5 
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

	head >> version >> totPar>>strands>>ccg>>ccg0>>noSpring>>noColor; //Saving header info
	if(ccg+ccg0 !=totPar) throw oxDNAException("Number of Colored + Color less particles is not equal to total Particles. Insanity caught in the header.");
	allocate_particles(particles);
	if(head.fail())
		throw oxDNAException("Please check the header, there seems to be something out of ordinary.");
	
	if(version<currentVersion)
		throw oxDNAException("Seems like the version of the file is an old one. Try looking into the documentation to upgrade it to newer one.");
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
		connection=true;
		bcall=true;
		while(body.tellg()!=-1){
			body>>temp;
			if(j==0){
				particleType=std::stoi(temp); // particle type = -2 for future references
			}else if(j==1){
				particles[i]->type=std::stoi(temp); // id of the particle
			}else if(j==2){
				particles[i]->btype=std::stoi(temp); //color of the particle
			}else if(j==3){
				q->radius=std::stoi(temp); //radius of the particle
			}else{
				if(connection){
					q->spring_neighbours.push_back(std::stoi(temp)); //all connected particle
					connection=false;
				}else if(bcall){
					q->Bfactor.push_back(std::stod(temp)); //connected particle Bfactor
					bcall=false;
				}else{
					q->ro.push_back(std::stod(temp)); //connected particle equilibrium spring distance r0
					connection=true;
					bcall=true;
				};
			}
			j++;
		}
		particles[i]->btype=color; //btype is used for coloring
		if(i==0)
			OX_LOG(Logger::LOG_INFO, "One loop successfully completed");
		i++;
	}
	*N_strands=strands;
};

void CCGInteraction::check_input_sanity(std::vector<BaseParticle*> &particles) {

}