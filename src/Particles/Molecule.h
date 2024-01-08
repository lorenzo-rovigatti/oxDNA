/*
 * Molecule.h
 *
 *  Created on: 21 apr 2020
 *      Author: lorenzo
 */

#ifndef SRC_PARTICLES_MOLECULE_H_
#define SRC_PARTICLES_MOLECULE_H_

#include "BaseParticle.h"

struct Molecule {
	Molecule();
	virtual ~Molecule();

	static void reset_id() {
		_current_id = 0;
	}

	std::vector<BaseParticle *> particles;
	LR_vector com;

	/**
	 * @brief Add a particle to the molecule
	 * @param p the particle to be added
	 */
	void add_particle(BaseParticle *p);

	/**
	 * @brief Reconstruct the coordinates of particles that belong to the molecule so that PBCs are not involved
	 *
	 * Particles that share molecular bonds should not be separated by a box size through PBC. This method attempts to
	 * reconstruct the positions of the particles so that they are all connected through short bonds.
	 */
	void normalise();

	/**
	 * @brief The number of particles that belong to the molecule.
	 */
	unsigned int N();

	/**
	 * @brief Compute the centre of mass of the molecule
	 */
	void update_com();

	int get_id() {
		return _id;
	}

	int get_topology_id() {
		return _id + 1;
	}

	bool shiftable();

private:
	static int _next_id() {
	   return _current_id++;
	}
	const int _id;
	static int _current_id;

	/// @brief true if the shiftable conditions should be re-evaluated
	bool _shiftable_dirty = false;
	bool _is_shiftable = true;
};

#endif /* SRC_PARTICLES_MOLECULE_H_ */
