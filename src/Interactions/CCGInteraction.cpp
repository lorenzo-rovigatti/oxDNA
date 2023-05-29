#include "CCGInteraction.h"
#include "../Particles/CCGParticle.h"
#include "../Utilities/Utils.h"

CCGInteraction::CCGInteraction() :
				BaseInteraction() {
	// ADD_INTERACTION_TO_MAP(CCG, _ccg_interaction);
}

CCGInteraction::~CCGInteraction() {

}

void CCGInteraction::get_settings(input_file &inp) {
	OX_LOG(Logger::LOG_INFO,"Setting is being called");
	BaseInteraction::get_settings(inp);
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
	if(pCG->has_bond(q)){
		return pair_interaction_bonded(p,q,compute_r,update_forces);
	}else{
		return pair_interaction_nonbonded(p,q,compute_r,update_forces);
	}
}

number CCGInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = spring(p,q,compute_r,update_forces);
	return energy;
}

number CCGInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {

	return (number) 0.f;
}

number CCGInteraction::spring(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	auto *pCG = dynamic_cast<CCGParticle*>(p);
	double k,r0;
	pCG->return_kro(q->index,&k,&r0);
	number dist=_computed_r.module()-r0; // Distance between the particles - equilibrium distance
	number energy = 0.5*k*SQR(dist); // Energy = 1/2*k*x^2
	if(update_forces){
		LR_vector force = (-1.f*k*_computed_r*dist/_computed_r.module()); //force = -k*(r_unit*dist) =-k(r-r0)
		p->force-= force;//substract force from p
		q->force+= force;//add force to q
	}
	return energy;
}

number CCGInteraction::exc_vol(BaseParticle *p, BaseParticle *q, bool compute_r,bool update_forces){
	return 0.f;
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
	allocate_particles(particles);
	if(head.fail())
		throw oxDNAException("Please check the header, there seems to be something out of ordinary.");
	
	if(version<currentVersion)
		throw oxDNAException("Seems like the version of the file is an old one. Try looking into the documentation to upgrade it to newer one.");
	i=0,color=0;
	while(topology.getline(line,2048)){ //Read till the end of file.
		if (strlen(line)==0 || line[0]=='#') //Ignore empty line or line starting with # for comment
			continue;
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
				particleType=std::stoi(temp);
			}else if(j==1){
				particles[i]->type=std::stoi(temp);
			}else if(j==2){
				color=std::stoi(temp);
			}else{
				if(connection){
					q->spring_neighbours.push_back(std::stoi(temp));
					connection=false;
				}else if(bcall){
					q->Bfactor.push_back(std::stod(temp));
					bcall=false;
				}else{
					q->ro.push_back(std::stod(temp));
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