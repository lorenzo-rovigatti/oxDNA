/*
 * JordanInteraction.h
 *
 *  Created on: 19/feb/2016
 *      Author: flavio
 */

#ifndef JORDANINTERACTION_H_
#define JORDANINTERACTION_H_

#include "BaseInteraction.h"

/**
 * @brief Manages the interaction between simple jordan particles (as described in http://jcp.aip.org/resource/1/jcpsa6/v131/i1/p014504_s1)
 *
 * This interaction is selected with
 * interaction_type = jordan
 *
 * @verbatim
 JORDAN_N_patches = <int> (number of patches)
 [JORDAN_s = <float> (sigma of the gaussian modulation, defaults to 0.3)]
 [JORDAN_m = <float> (exponent to the 2m-m lennard-jones part)]
 [JORDAN_phi = <float> (angle below the equator for the rest position of the patches, defaults to PI/6)]
 [JORDAN_int_k = <float> (stiffness of the internal spring, defaults to 0., i.e., free patches)]
 @endverbatim
 */

class JordanInteraction: public BaseInteraction {
protected:
	/// number of patches on each particle
	int _n_patches;

	/// sigma of the gaussian modulation
	number _s;

	/// sigma of simmetry breaking gaussian function
	number _SB_s;

	/// angle below the equator
	number _phi;

	/// "sensitivity" of the LJ contribution
	int _m;

	/// stiffness of the internal spring
	number _int_k;

	/// interaction
	number _jordan_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	/// number of particles with N3 patches and of particles in general
	int _my_N3, _my_N;

	std::string patches_file_in;
	std::string patches_file_out;

public:
	enum {
		JORDAN = 0
	};

	JordanInteraction();
	virtual ~JordanInteraction();

	virtual void get_settings(input_file &inp);
	virtual void init();

	virtual void allocate_particles(std::vector<BaseParticle *> &particles);

	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles);
};

#endif /* JORDANINTERACTION_H_ */
