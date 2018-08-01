/*
 * Icosahedron.h
 *
 *  Created on: 25/lug/2018
 *      Author: Flavio
 */

#ifndef ICOSAHEDRON_H_
#define ICOSAHEDRON_H_

#include "Particles/BaseParticle.h"

/**
 * @brief Incapsulates a patchy particle with 2, 3, or 4 spherical patches. Used by PatchyInteraction.
 */
template<typename number>
class Icosahedron : public BaseParticle<number> {
protected:
	/**
	 * number of patches
	 */
	int _N_patches;

	/**
	 * array containing the base vertexes
	 */
	LR_vector<number> * _vertexes; 

	/**
	 * array containing the indexes of the vertexes that have
	 * patches on them
	 * TODO: perhaps change to short int to save memory?
	 */
	int * _patch_indexes;

	/** 
	 * array storing the indexes of the vertexes that form each face
	 * TODO: perhaps change to short int to save memory?
	 */
	int * _faces;

	/**
	 * array storing to which faces each vertex belongs to
	 * TODO: perhaps change to short int to save memory?
	 */
	int * _close_faces; 

	void _set_vertexes();

public:
	Icosahedron(int N_patches);
	virtual ~Icosahedron();

	void set_positions();

	virtual bool is_rigid_body() {
		return true;
	}
};

#endif /* ICOSAHEDRON_H_ */
