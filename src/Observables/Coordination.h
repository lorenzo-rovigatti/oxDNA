/*
 * Coordination.h
 *
 * Created on: 10/17/2025
 *      Author: Lorenzo
*/

#ifndef COORDINATION_H
#define COORDINATION_H

#include "BaseObservable.h"

#include "../Forces/Metadynamics/meta_utils.h"

class Coordination: public BaseObservable {
public:
    Coordination();

    virtual ~Coordination();

    virtual void init();
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
    virtual std::string get_output_string(llint curr_step);

private:
    std::string _op_file;

    meta::CoordSettings _settings;
    std::vector<std::pair<BaseParticle *, BaseParticle *>> _all_pairs;
};

#endif /* COORDINATION_H */
