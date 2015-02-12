/*
 * CoaxVariables.h
 *
 *  Created on: Mar 11, 2014
 *      Author: Ben Snodin
 */

#ifndef COAXVARIABLES_H_
#define COAXVARIABLES_H_

#include "BaseObservable.h"
// nasty hack to steal the protected _f2(), _f5() etc methods from the DNAInteraction class
#define protected public
#include "../Interactions/DNAInteraction.h"
#undef protected

/**
 * @brief Outputs the relevant angle and position variables for the coaxial stacking term between 2 specified particles in the oxDNA potential
 *
 * @verbatim
particle1_id = <int> (particle 1 id)
particle2_id = <int> (particle 2 id)
@endverbatim
 */
template<typename number>
class CoaxVariables: public BaseObservable<number> {
protected:
	int _particle1_id;
	int _particle2_id;
	DNAInteraction<number> *_dna_interaction;
	bool _use_oxDNA2_coaxial_stacking;
public:
	CoaxVariables();
	virtual ~CoaxVariables();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);
};

#endif /* COAXVARIABLES_H_ */
