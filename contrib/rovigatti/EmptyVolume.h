/*
 * EmptyVolume.h
 *
 *  Created on: Jan 14, 2016
 *      Author: Lorenzo Rovigatti
 */

#ifndef EMPTYVOLUME_H_
#define EMPTYVOLUME_H_

#include "../Lists/Cells.h"

#include "BaseObservable.h"
#include "../Interactions/DNAInteraction.h"
#include "../Utilities/OrderParameters.h"

/**
 * @brief Prints four columns: construct_id_1, construct_id_2, base_pairs_between, minimum_image_distance
 */
template<typename number>
class EmptyVolume : public BaseObservable<number>  {
protected:
	number _probe_diameter;
	number _particle_diameter;
	number _rcut, _sqr_rcut;
	Cells<number> *_cells;
	BaseParticle<number> *_probe;
	BaseParticle<number> **_particles;
	int _N;
	int _tries;

public:
	EmptyVolume();
	virtual ~EmptyVolume();

	std::string get_output_string(llint curr_step);

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	virtual void init(ConfigInfo<number> &config_info);
};

extern "C" BaseObservable<float> *make_float() { return new EmptyVolume<float>(); }
extern "C" BaseObservable<double> *make_double() { return new EmptyVolume<double>(); }

#endif /* CONSTRUCTWISE_H_ */
