/*
 * JordanParticle.h
 *
 *  Created on: 15/mar/2013
 *      Author: lorenzo
 */

#ifndef JORDANPARTICLE_H_
#define JORDANPARTICLE_H_

#include "BaseParticle.h"
#include "../Utilities/Utils.h"

/**
 * @brief Incapsulates a jordan particle with 2, 3, or 4 spherical patches. Used by JordanInteraction.
 */
template<typename number>
class JordanParticle : public BaseParticle<number> {
protected:
	LR_vector<number> * _base_patches; // patch equilibrium positions in the lab reference frame
	LR_matrix<number> * _patch_rotations; // matrixes for internal degrees of freedom of the patches

	void _set_base_patches(number phi);
	
	number _int_k; // stiffness of the internal spring

public:
	JordanParticle(int npatches, number phi, number int_k);
	virtual ~JordanParticle();

	void set_positions();
	
	LR_matrix<number> get_patch_rotation (int i) { return LR_matrix<number> (_patch_rotations[i].v1, _patch_rotations[i].v2, _patch_rotations[i].v3); }
	void set_patch_rotation (int i, LR_matrix<number> R) { _patch_rotations[i] = R; }
	void rotate_patch (int i, LR_matrix<number> R) { _patch_rotations[i] = R * _patch_rotations[i]; }
	number int_potential ();

	virtual bool is_rigid_body() {
		return true;
	}
	
	std::string get_output_string ();
};

#endif /* JORDANPARTICLE_H_ */
