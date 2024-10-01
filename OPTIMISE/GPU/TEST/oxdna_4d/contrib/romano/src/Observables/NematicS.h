/*
 * NematicS.h
 *
 *  Created on: Mar 10, 2016
 *      Author: flavio
 */

#ifndef NEMATICS_H_
#define NEMATICS_H_

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

class NematicS : public BaseObservable {
private:

	int _axis_i;

	Cells * _cells;
	BaseParticle **_particles;

	LR_vector * _get_axis(BaseParticle *);

	// cutoffs used for the clustering
	number _cluster_n_cutoff;
	number _cluster_angle_cutoff;
	number _cluster_distance_cutoff;
	number _rod_length;

public:
	// we choose which vector to use on the particle
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

	virtual void init();
	virtual void get_settings(input_file &my_inp, input_file &sim_inp);

	std::string get_output_string(llint curr_step);

	//std::vector<std::list<int> > build_clusters();
};

extern "C" BaseObservable * make_NematicS() { return new NematicS(); }

#endif /* NEMATICS_H_ */
