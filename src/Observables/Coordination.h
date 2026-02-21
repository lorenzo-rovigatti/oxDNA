/*
 * Coordination.h
 *
 * Created on: 10/17/2025
 *      Author: Lorenzo
*/

#ifndef COORDINATION_H
#define COORDINATION_H

#include "BaseObservable.h"

#include "../Utilities/OrderParameters.h"

class Coordination: public BaseObservable {
public:
    Coordination();

    virtual ~Coordination();

    virtual void init();
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
    virtual std::string get_output_string(llint curr_step);

private:
    std::string _op_file;
    number _d0 = 1.2; // optimal distance for coordination, defaults to 1.2 in internal units
    number _r0 = 0.5; // width of the switching function, defaults to 0.5 in internal units
    int _n = 6; // exponent of the switching function, defaults to 6
    std::vector<std::pair<BaseParticle *, BaseParticle *>> _all_pairs;
};

#endif /* COORDINATION_H */
