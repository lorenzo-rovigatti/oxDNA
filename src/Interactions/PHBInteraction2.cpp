#include "PHBInteraction2.h"

PHBInteraction2::PHBInteraction2():BaseInteraction(){
    // Constructor
};

PHBInteraction2::~PHBInteraction2(){
    // Destructor
};

void PHBInteraction2::get_settings(input_file &inp){
    BaseInteraction::get_settings(inp); // Get basic settings
}

void PHBInteraction2::init(){
}

void PHBInteraction2::allocate_particles(std::vector<BaseParticle *> &particles){
    // Add particle to the system
    particles.resize(totPar);
	for(i=0;i<totPar;i++){
		particles[i] = new CCGParticle();
		particles[i]->index=i;
		particles[i]->type=0;
	}
}

void PHBInteraction2::check_input_sanity(std::vector<BaseParticle *> &particles){
    // Check all the input file are correct.
}

void PHBInteraction2::read_topology(int *N_strands, std::vector<BaseParticle *> &particles){
    // Read the top file
}