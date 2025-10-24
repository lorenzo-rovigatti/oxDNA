//
// Created by josh on 8/14/25.
//

#ifndef OXDNA_PLCLUSTERTOPOLOGY_H
#define OXDNA_PLCLUSTERTOPOLOGY_H

/*
 * HBEnergy.h
 *
 *  Created on: Feb 14, 2013
 *      Author: petr
 */

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
class PLClusterTopology: public BaseObservable {
protected:

    bool _show_types ;
public:
    PLClusterTopology();
    ~PLClusterTopology() override;

    std::string get_output_string(llint curr_step) override;

    void get_settings(input_file &my_inp, input_file &sim_inp) override;
};

extern "C" BaseObservable * make_PLClusterTopology() { return new PLClusterTopology(); }


#endif //OXDNA_PLCLUSTERTOPOLOGY_H
