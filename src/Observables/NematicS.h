/*
 * NematicS.h
 *
 *  Created on: Mar 10, 2016
 *      Author: flavio
 */

#ifndef NEMATICS_H_
#define NEMATICS_H_

#include "BaseObservable.h"
/**
 * @brief Outputs the mematic order parameter. It uses as a reference vector one of the 6 orientation vectors, which can be specified in the input file. It can also be configured to use an pre-defined nematic axis rather than computing one.
 *
 * To use this observable, use: type = nematic_s
 *
 * This observable takes four optional arguments
 * @verbatim

ref_vec = <int> (between 0 and 5, corresponding to ->orientation->v1, v2, v3, p->orientationT.v1, v2, v3 in this order)
[dir = <float>, <float>, <float>] (set a pre-defined nemati axis rather than computing one)
@endverbatim
 */

template<typename number>
class NematicS : public BaseObservable<number> {
private:
	// arguments
	LR_vector<number> _dir;

	bool _external_dir;

	int _axis_i;

	LR_vector<number> * _get_axis(BaseParticle<number> *);

public:
	// we can use sites
	enum {
		O_V1 = 0,
		O_V2 = 1,
		O_V3 = 2,
		OT_V1 = 3,
		OT_V2 = 4,
		OT_V3 = 5
	};

	NematicS();
	virtual ~NematicS();

	virtual void init(ConfigInfo<number> &config_info);
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);

	void _compute_dir () { return; } // eventually, there will be the possibility to compute the directrion
};

#endif /* NEMATICS_H_ */
