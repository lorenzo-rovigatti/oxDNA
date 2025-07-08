/*
 * HBEnergy.h
 *
 *  Created on: Feb 14, 2013
 *      Author: petr
 */

#ifndef PLCLUSTER_H_
#define PLCLUSTER_H_

#define MCMOVE_CUSTOM
#include "../../../../src/Observables/BaseObservable.h"

/**
 * @brief Prints out all the interactions between all pairs of nucleotides with non-zero interaction energies (default) or between a pair of nucleotides (if specified)
 *
 * @verbatim
particle1_id = <int> (particle 1 id)
particle2_id = <int> (particle 2 id)
@endverbatim
 */
class PLCluster: public BaseObservable {
protected:

    bool _show_types ;
public:
	PLCluster();
	~PLCluster() override;

	std::string get_output_string(llint curr_step) override;

	void get_settings(input_file &my_inp, input_file &sim_inp) override;
};

#endif /* PAIRENERGY_H_ */
