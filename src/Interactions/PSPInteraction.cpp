#include "PSPInteraction.h"

PSPInteraction::PSPInteraction():BaseInteraction(){
    // Constructor
};

PSPInteraction::~PSPInteraction(){
    // Destructor
};

void PSPInteraction::get_settings(input_file &inp){
    BaseInteraction::get_settings(inp); // Get basic settings
}

void PSPInteraction::init(){
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
	while(std::getline(topology,line)){
		if(line.empty()||line[0]=='#') continue; //skip empty line and comment
		if(line[0]=='i'){
			if(line[1]=='P'){ //these are Patchy informations
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					
				}
			}
		}
	}
}

number PSPInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	// Check bonded or non-bonded
	return 0;
}

number PSPInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	// Bonded particle interaction
	return 0;
}

number PSPInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	// Non-bonded particle interaction
	return 0;
}