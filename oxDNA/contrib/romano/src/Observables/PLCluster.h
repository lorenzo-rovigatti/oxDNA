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
template<typename number>
class PLCluster: public BaseObservable<number> {
protected:

    bool _show_types ;
public:
	PLCluster();
	virtual ~PLCluster();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
};

extern "C" BaseObservable<float> * make_PLCluster_float() { return new PLCluster<float>(); }
extern "C" BaseObservable<double> * make_PLCluster_double() { return new PLCluster<double>(); }


#endif /* PAIRENERGY_H_ */
