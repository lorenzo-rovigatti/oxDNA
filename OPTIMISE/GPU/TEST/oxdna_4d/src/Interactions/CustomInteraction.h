/*
 * CustomInteraction.h
 *
 *  Created on: 6 Mar 2014
 *      Author: lorenzo
 */

#ifndef CUSTOMINTERACTION_H_
#define CUSTOMINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Custom interaction that makes use of lookup tables to support both custom bonded and non-bonded interactions. Under active development.
 */

class CustomInteraction: public BaseInteraction {
protected:
	char _lt_filename[512];
	Mesh _non_bonded_mesh;
	Mesh _bonded_mesh;
	int _bonded_points;
	number _Ecut;

	struct base_function {
		int points;
		number *x, *fx, *dfx;
	};

	/**
	 * @brief Performs a linear interpolation on the x_data and fx_data array to compute f(x).
	 *
	 * @param x
	 * @param x_data
	 * @param fx_data
	 * @param points
	 * @return
	 */
	number _linear_interpolation(number x, number *x_data, number *fx_data, int points);

	/**
	 * @brief Performs a linear interpolation on the lookup tables provided by the last parameter to calculate the value of f(x).
	 *
	 * This method is required during the building of the actual lookup table. It expects the last parameter to be of type base_function.
	 *
	 * @param x
	 * @param par should be a (void *) base_function
	 * @return
	 */
	number _fx(number x, void *par);

	/**
	 * @brief Just like _fx(), but computes the derivative of f(x).
	 *
	 * See _fx().
	 *
	 * @param x
	 * @param par should be a (void *) base_function
	 * @return
	 */
	number _dfx(number x, void *par);
public:
	enum {
		BONDED = 0, NONBONDED = 1
	};

	CustomInteraction();
	virtual ~CustomInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

#endif /* CUSTOMINTERACTION_H_ */
