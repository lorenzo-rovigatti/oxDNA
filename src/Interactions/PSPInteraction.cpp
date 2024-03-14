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
		int countP=0;int countC=0;
		if(line[0]=='i'){
			if(line[1]=='P'){ //these are Patchy informations
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					body>>temp;
					if(j==0) continue; //skip the first word
					if(j==1 && stoi(temp) !=0) throw oxDNAException("Not supported yet");
					if(j>6) continue; //skip the rest
					patches[countP][j-2]=stoi(temp);
					j++;
				}
				countP++;
			}
			if(line[1]=='C'){
				std::stringstream body(line);
				j=0;
				while(body.tellg()!=-1){
					body>>temp; //skip word "iC"
					particlePatches[countC][j]=stoi(temp);
					j++;
				}
				countC++;
			}
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
				p->radius=stof(temp);
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

	//Setting up the matrix
	#pragma omp parallel for
	for(i=0;i<totPar;i++){
		auto *p = static_cast<CCGParticle *>(particles[i]);
		connections[i][0]=p->spring_neighbours.size();
		r0[i][0]=0; // all the eq radius are different
		k0[i][0]=0; // all the spring constant are different
		for(j=0;j<(int)p->spring_neighbours.size();j++){
			connections[i][j+1]=p->spring_neighbours[j];
			r0[i][j+1]=p->ro[j];
			k0[i][j+1]=p->Bfactor[j];
		};
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

number PSPInteraction::exeVolInt(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces){
	if(compute_r) _computed_r = _box->min_image(p->pos,q->pos);
	number totalRadius = particleTopology[p->index][2]+particleTopology[q->index][2];
	
	return 0;
}