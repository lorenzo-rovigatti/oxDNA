/**
 * @file    BaseParticle.h
 * @date    21/set/2210
 * @author   lorenzo
 *
 *
 */

#ifndef BASEPARTICLE_H_
#define BASEPARTICLE_H_

#include <cstring>
#include <cstdlib>
#include <cassert>

#include <stdio.h>

#include "../defs.h"
#include "../Boxes/BaseBox.h"
#include "../Forces/BaseForce.h"

class ParticlePair;

/**
 * @brief Base particle class. All particles must inherit from this class.
 */
class BaseParticle {
protected:
	void _check();

public:
	std::vector<BaseForce *> ext_forces;
	int next_particle;
	int strand_id;

	BaseParticle();
	virtual ~BaseParticle();

	std::vector<ParticlePair> affected;

	virtual void set_positions() {

	}

	/// number of boxes the particle has diffused in each direction
	int _pos_shift[3];

	virtual void copy_from(const BaseParticle &);
	inline void soft_copy_from(const BaseParticle * p) {
		pos = p->pos;
		orientation = p->orientation;
		orientationT = p->orientationT;
		en3 = p->en3;
		en5 = p->en5;
		esn3 = p->esn3;
		esn5 = p->esn5;
	}

	number en3, en5, esn3, esn5;
	bool inclust;

	void init();

	int get_index() const {
		return index;
	}

	/**
	 * @brief Add an external force.
	 *
	 * @param f
	 * @return true if the force was added, false otherwise
	 */
	bool add_ext_force(BaseForce *f);

	inline void set_initial_forces(llint step, BaseBox *box) {
		if(is_rigid_body()) {
			torque = LR_vector((number) 0.f, (number) 0.f, (number) 0.f);
		}
		force = LR_vector((number) 0.f, (number) 0.f, (number) 0.f);

		if(ext_forces.size() > 0) {
			LR_vector abs_pos = box->get_abs_pos(this);
			for(auto ext_force : ext_forces) {
				force += ext_force->value(step, abs_pos);
			}
		}
	}

	/**
	 * @brief Computes the interaction resulting from all the external forces acting on the particle. Stores the result in the ext_potential member.
	 *
	 * @param step current time step. Useful for forces that depend on time.
	 * @param box pointer to the box object
	 * @return true if the external potential was added, false otherwise
	 */
	inline void set_ext_potential(llint step, BaseBox *box) {
		if(ext_forces.size() > 0) {
			LR_vector abs_pos = box->get_abs_pos(this);
			ext_potential = (number) 0.;
			for(auto ext_force : ext_forces) {
				ext_potential += ext_force->potential(step, abs_pos);
			}
		}
	}

	/**
	 * @brief Checks whether q and the current particle are bonded neighbours (such as neighbouring particles on a DNA strand).
	 *
	 * @param q candidate bonded neighbour
	 * @return true if the current particle and q are bonded neighbours, false otherwise
	 */
	virtual bool is_bonded(BaseParticle *q) {
		return false;
	}

	/**
	 * @brief Defaults to false.
	 *
	 * @return true if the particle is a rigid body (i.e. orientational degrees of freedom are to be taken into account), false otherwise
	 */
	virtual bool is_rigid_body() {
		return false;
	}

	inline void set_pos_shift(int x, int y, int z) {
		_pos_shift[0] = x;
		_pos_shift[1] = y;
		_pos_shift[2] = z;
	}

	void get_pos_shift(int *arg) {
		arg[0] = _pos_shift[0];
		arg[1] = _pos_shift[1];
		arg[2] = _pos_shift[2];
	}

	virtual uint N_int_centers() {
		return int_centers.size();
	}

	/// Index of the particle. Usually it is a useful way of accessing arrays of particles
	int index;

	/// Particle type (the meaning of which depends on the chosen interaction type)
	int type;

	/// Needed for specific base pairing
	int btype;
	
	/// DNA bonded neighbours
	BaseParticle *n3, *n5;

	/// Torque exerted on the particle in its own reference system
	LR_vector torque;

	/// Force exerted on the particle
	LR_vector force;

	/// Total potential energy due to external forces
	number ext_potential;

	/// Particle position inside the box
	LR_vector pos;

	/// Positions of all interaction centers. This array must be initialized by child classes
	std::vector<LR_vector> int_centers;

	/// Angular momentum of the particle
	LR_vector L;

	/// Velocity of the particle
	LR_vector vel;

	/// BaseParticle orientational matrix
	LR_matrix orientation;

	/// transpose (= inverse) orientational matrix
	LR_matrix orientationT;
};

/*
 * helpers
 */

class ParticlePair {
public:
	BaseParticle *first;
	BaseParticle *second;

	ParticlePair(BaseParticle *p, BaseParticle *q) {
		if(p == q)
			throw oxDNAException("ParticlePair: p == q is not allowed");
		if(p->index < q->index) {
			first = p;
			second = q;
		}
		else {
			first = q;
			second = p;
		}
	}

	bool operator<(ParticlePair q) const {
		int p1 = first->index;
		int p2 = second->index;
		int q1 = q.first->index;
		int q2 = q.second->index;

		if(p1 == q1)
			return (p2 < q2);
		else
			return (p1 < q1);
	}
};

#endif /* BASEPARTICLE_H_ */
