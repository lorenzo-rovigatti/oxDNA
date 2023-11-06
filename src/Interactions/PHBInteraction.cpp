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
    return 0;
};
number PHBInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    return 0;
};
number PHBInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
    return 0;
};

number PHBInteraction::_repulsive_lj2(number prefactor, const LR_vector &r, LR_vector &force, number sigma, number rstar, number b, number rc, bool update_forces) {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;
	if(rnorm < SQR(rc)) {
		if(rnorm > SQR(rstar)) {
			number rmod = sqrt(rnorm);
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

number PHBInteraction::spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	return 0;
}

number PHBInteraction::maxRadius(std::vector<PHBParticle*> &particles){
	number radius =0;
	for(i=0;i<totPar;i++){
		if(radius<particles[i]->radius) radius=particles[i]->radius;
	}
	return radius;
};

void PHBInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles){
    // char line[2048];
	std::ifstream topology;
	topology.open(this->_topology_filename,std::ios::in);
	if(!topology.good()) throw oxDNAException("Problem with topology file: %s",this->_topology_filename);

	std::getline(topology,temp);
	std::cout<<temp<<std::endl;
	std::stringstream head(temp);
	head>>totPar>>strands; //saving header info
	allocate_particles(particles);
	while(std::getline(topology,temp)){
		if(temp.empty()||temp[0]=='#') continue; //skip empty line and comment
		
	}

};