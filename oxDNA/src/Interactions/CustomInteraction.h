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
template <typename number>
class CustomInteraction: public BaseInteraction<number, CustomInteraction<number> > {
protected:
	char _lt_filename[512];
	Mesh<number> _bonded_mesh;
	int _bonded_points;

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
		BONDED = 0,
		NONBONDED = 1
	};

	CustomInteraction();
	virtual ~CustomInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(BaseParticle<number> **particles, int N);
	virtual void read_topology(int N, int *N_strands, BaseParticle<number> **particles);

	virtual number pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false);
	virtual number pair_interaction_term(int name, BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r=NULL, bool update_forces=false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, r, update_forces);
	}

	virtual void check_input_sanity(BaseParticle<number> **particles, int N);
};

#endif /* CUSTOMINTERACTION_H_ */
