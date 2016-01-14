/**
 * @file    ConstantRateForce.h
 * @date    18/oct/2011
 * @author  Flavio
 *          Modified by Ferdinando on 11/Dec/2015
 *
 */

#ifndef CONSTANTRATEFORCE_H_
#define CONSTANTRATEFORCE_H_

#include "BaseForce.h"
#include "../Utilities/Utils.h"

/**
 * @brief Force that increases with a constant (possibly 0) rate
 *
 * This is a simple force, where one or more particles are pulled with a force that 
 * is either constant or grows linearly with time in a fixed direction.
 * All units are oxDNA units (time, force = energy / distance, ...)
 * arguments:
@verbatim
particle = <int> (comma-separated list of indices of particles to apply the force to. -1 applies it to all particles. Entries separated by a dash "-" get expanded in a list of all the particles on a same strand comprised between the two indices. E.g., particle= 1,2,5-7 applies the force to 1,2,5,6,7 if 5 and 7 are on the same strand.)
F0 = <float> (Initial force.)
rate = <float> (growth rate of the force. It is [oxDNA energy units / (oxDNA distance units * (MD/MC) steps].)
@endverbatim
 */
template<typename number>
class ConstantRateForce : public BaseForce<number> {
private:
	std::string _particles_string;

public:
	ConstantRateForce();
	virtual ~ConstantRateForce();

	void get_settings (input_file &);
	void init (BaseParticle<number> **, int, number *);

	virtual LR_vector<number> value(llint step, LR_vector<number> &pos);
	virtual number potential (llint step, LR_vector<number> &pos);

	std::vector<int> getParticlesFromString(BaseParticle<number> **particles, int N, std::string particle_string);
	void assert_is_valid_particle(int n,int N);
};

template<typename number>
std::vector<int> ConstantRateForce<number>::getParticlesFromString(BaseParticle<number> **particles, int N,std::string particle_string){
	// first remove all the spaces from the string, so that the parsing goes well.
	particle_string.erase(remove_if(particle_string.begin(), particle_string.end(), static_cast<int(*)(int)>( isspace )), particle_string.end());

  std::vector<std::string> temp = Utils::split (particle_string.c_str(), ',');
	std::vector<int> particles_index;//declare as empty

	for( std::vector<std::string>::size_type i = 0; i < temp.size(); i++){
		bool found_dash = temp[i].find('-') != std::string::npos;
		// if the string contains a dash, then it has to be interpreted as a list of particles
		// unless it's a negative number
		//if (found_dash && strcmp("-1",temp[i].c_str()) != 0 ){
		if (found_dash && '-'!= temp[i].c_str()[0] ){
			// get the two indices p1 and p2 and check they make sense
			std::vector<std::string> p1_p2_index = Utils::split(temp[i].c_str(),'-');

			int p1 = atoi(p1_p2_index[0].c_str());
			int p2 = atoi(p1_p2_index[1].c_str());
			assert_is_valid_particle(p1,N);
			assert_is_valid_particle(p2,N);

			int j = p1;
			// add all the particles between p1 and p2 (extremes included)
			bool found_p2 = false;
			do{
				particles_index.push_back(j);
				if (j == p2){
					found_p2 = true;
				}
				if (particles[j]->n5 == P_VIRTUAL) break;
				j = particles[j]->n5->index;

			} while( j != p1 && !found_p2);
			// check that it hasn't got to either the end of the strand or back to p1
			if(!found_p2){
				throw oxDNAException("In force ConstantRateForce I couldn't get from particle %d to particle %d.",p1,p2);
			} 

		}	
		else{//just add it to the vector
			int j = atoi(temp[i].c_str());
			
			assert_is_valid_particle(j,N);
			particles_index.push_back(j);
		}
		
	}
	// check that if -1 is present then that's the only key - something must be wrong if you
  // specified -1 (all particles) and then some more particles.
	if (std::find(particles_index.begin(),particles_index.end(),-1) != particles_index.end() && particles_index.size()>1){
		throw oxDNAException("In force ConstantRateForce there are more than one particle index, including -1. If -1 is a particle index then it has to be the only one. Dying badly.");
	}
	// check that no particle appears twice
	for( std::vector<int>::size_type i = 0; i < particles_index.size(); i++){
		for( std::vector<int>::size_type j = i+1; j < particles_index.size(); j++){
			if ( particles_index[i] == particles_index[j] ){
				throw oxDNAException("In force ConstantRateForce particle index %d appears twice (both at position %d and at position %d), but each index can only appear once. Dying badly.",particles_index[i],i+1,j+1);
			}
		}
	}	
	// finally return the vector.
	return particles_index;
	

}

template<typename number>
void ConstantRateForce<number>::assert_is_valid_particle(int index, int N){
	if (index >= N || index < -1){
		throw oxDNAException ("Trying to add a ConstantRateForce on non-existent particle %d. Aborting", index);
	}
}

#endif /* CONSTANTRATEFORCE_H_ */
