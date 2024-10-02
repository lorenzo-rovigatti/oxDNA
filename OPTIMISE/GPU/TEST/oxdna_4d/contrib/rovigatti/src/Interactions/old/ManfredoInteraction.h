/*
 * ManfredoInteraction.h
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#ifndef MANFREDOINTERACTION_H_
#define MANFREDOINTERACTION_H_

#include "Interactions/BaseInteraction.h"
#include "Interactions/DNAInteraction.h"
#include "Particles/CustomParticle.h"
#include "Particles/DNANucleotide.h"

/**
 * @brief Handles the interaction between coarse-grained DNA tetramers.
 */
class ManfredoInteraction: public BaseInteraction<ManfredoInteraction> {
protected:
	enum {
		CENTRE, ARM, STICKY
	};

	// look up tables for intra-tetramer potentials
	char _intra_filename[3][512];
	Mesh _intra_mesh[3];
	int _intra_points[3];

	// look up tables for inter-tetramer potentials
	char _inter_filename[3][512];
	Mesh _inter_mesh[3];
	int _inter_points[3];

	number _T;

	// this is used only during the look-up table set up
	struct lt_data {
		int points;
		number *x, *fx, *dfx;
	};

	/**
	 * @brief Performs a linear interpolation on the x_data and fx_data array to compute f(x).
	 *
	 * @param x
	 * @param x_data
	 * @param fx_data
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

	void _build_lt(Mesh &mesh, int points, char *filename);

	int _N_tetramers;
	const int _N_per_tetramer;

	DNAInteraction _DNA_inter;

	// helper functions
	bool _any(int p_type, int q_type, int to_test) {
		return (p_type == to_test || q_type == to_test);
	}
	bool _both(int p_type, int q_type, int to_test) {
		return (p_type == to_test && q_type == to_test);
	}

	int _get_p_type(BaseParticle *p);
	int _get_inter_type(BaseParticle *p, BaseParticle *q);

	/**
	 * @brief Query the given mesh at the value x.
	 *
	 * This method is extended to support linear extrapolation for x-values which are out of range.
	 *
	 * @param x
	 * @param m
	 * @return
	 */
	number _query_mesh(number x, Mesh &m);
	/**
	 * @brief Query the derivative of the given mesh at the value x.
	 *
	 * This method is extended to support linear extrapolation for x-values which are out of range.
	 *
	 * @param x
	 * @param m
	 * @return
	 */
	number _query_meshD(number x, Mesh &m);
public:
	enum {
		CENTRE_ARM_INTRA, ARM_ARM_NEAR, ARM_ARM_FAR,
	};

	enum {
		CENTRE_ARM_INTER, ARM_ARM_INTER, CENTRE_CENTRE, CENTRE_STICKY, ARM_STICKY, STICKY_STICKY
	};

	ManfredoInteraction();
	virtual ~ManfredoInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) {
		return this->_pair_interaction_term_wrapper(this, name, p, q, compute_r, update_forces);
	}

	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles);
};

// Arms are composed by 7 particles. The last six, which constitute the sticky end, are regular
// DNA nucleotides. The first one, on the other hand, is a DNA nucleotide but it also interacts
// via user-defined potentials with tetramers' centres and other arms

class CustomArmParticle: public CustomParticle {
protected:
	LR_vector _principal_DNA_axis;
public:
	CustomArmParticle();
	virtual ~CustomArmParticle() {
	}
	;

	virtual bool is_rigid_body() {
		return true;
	}

	virtual bool is_bonded(BaseParticle *q) {
		if(q == this->n3 || q == this->n5) return true;
		return CustomParticle::is_bonded(q);
	}

	void set_positions();
};

extern "C" ManfredoInteraction *make_ManfredoInteractiont() {
	return new ManfredoInteraction();
}

#endif /* MANFREDOINTERACTION_H_ */
