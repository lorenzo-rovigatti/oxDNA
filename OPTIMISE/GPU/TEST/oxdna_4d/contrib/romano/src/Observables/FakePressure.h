/*
 * FakePressure.h
 *
 *  Created on: Mar 10, 2016
 *      Author: flavio
 */

#ifndef FAKEPRESSURE_H_
#define FAKEPRESSURE_H_

#include "Observables/BaseObservable.h"
#include "Interactions/InteractionUtils.h"
#include "Lists/Cells.h"
#include <list>


/**
 * @brief Outputs the nematic order parameter and the three components of the
 * nematic director. It uses as a reference vector on each particle one of the
 * 6 orientation vectors, which can be specified in the input file. The
 * observable can also be used with a custom director: in this case, the
 * director won't be computed and the one provided will be used
 *
 * To use this observable, use: type = nematic_s
 *
 * This observable takes four optional arguments
 * @verbatim
ref_vec = <int> (between 0 and 5, corresponding to p->orientation->v1, v2, v3, p->orientationT.v1, v2, v3 in this order)
[custom_dir = <bool>] (default: False. To choose a pre-defined )
[dir = <float>, <float>, <float>] (mandatory if custom_dir is True. Set a pre-defined nematic director rather than computing one)
@endverbatim
 */

template<typename number>
class FakePressure : public BaseObservable<number> {
private:
	Cells<number> * _cells;
	BaseParticle<number> **_particles;
	number _sigma;
	int _restrict_to_type;

public:

	FakePressure();
	virtual ~FakePressure();

	virtual void init(ConfigInfo<number> &config_info);
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);

};

extern "C" BaseObservable<float> * make_FakePressure_float() { return new FakePressure<float>(); }
extern "C" BaseObservable<double> * make_FakePressure_double() { return new FakePressure<double>(); }

#endif /* FAKEPRESSURE_H_ */
