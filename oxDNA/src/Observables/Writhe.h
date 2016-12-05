/*
 * Writhe.h
 *
 *  Created on: 16 Sep, 2015
 *      Author: Ferdinando Randisi
 */

#ifndef WRITHE_H_
#define WRITHE_H_

#include "BaseObservable.h"
/**
 * @brief Outputs the writhe of a strand, or the writhes of all the subdomains of a chain of the form [L,L+d], for all possible values of L and some d. It has been tested for the TEP model, but should in principle work
 for DNA as well.
 *
 * To use this observable, use type = writhe.
 *
 * This observable takes four optional arguments
 * @verbatim

first_particle_index = <int> (defaults to 0. index of the first particle on which to compute the angle with the next particle)

last_particle_index = <int> (defaults to the index of the first-but last bead in the same strand as the first particle. Therefore, if I have a strand of N beads, the last one will be the one with index N-2. This is because the last bead is atypical in the TEP model (e.g. it's aligned with the vector before it rather than the one in front of it.). index of the last particle of the chain on which to compute the angle. )

subdomain_size = <int> (defaults to the entire chain. If smaller, the writhe will be computed only on 
@endverbatim
 */

template<typename number>
class Writhe : public BaseObservable<number> {
private:
	// arguments
	int _first_particle_index;
	int _last_particle_index;
	int _subdomain_size;

	number ** _writhe_integrand_values;
	
	inline number _writheIntegrand(LR_vector<number> t,LR_vector<number> tp,LR_vector<number> r,LR_vector<number> rp);
public:
	Writhe();
	virtual ~Writhe();

	virtual void init(ConfigInfo<number> &config_info);
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	//virtual std::string get_output_stringOLD(llint curr_step);
	

	std::string get_output_string(llint curr_step);
};

template<typename number>
number Writhe<number>::_writheIntegrand(LR_vector<number> t,LR_vector<number> tp,LR_vector<number> r,LR_vector<number> rp){
	// (1/2pi) * t x t' . ( r - r')/ | r - r'|**3
	// The formula for writhe usually reports a prefactor of 1/4pi, but that's when both integration variables are left free. Since the double integral is symmetrical, we run the second integer with values greater than the first, hence the factor of 2 in the formula above.
	number rm = ( r - rp ).module();
	//return (1./2*M_PI)*((tp.cross(t))*( r - rp ))/(  rm*rm*rm );
	return ((t.cross(tp))*( r - rp ))/(2*M_PI*rm*rm*rm);
}

#endif /* WRITHE_H_ */
