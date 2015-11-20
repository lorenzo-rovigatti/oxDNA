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
#include "../Forces/BaseForce.h"

template <typename number> class ParticlePair;

/**
 * @brief Base particle class. All particles must inherit from this class.
 */
template<typename number>
class BaseParticle {
protected:
	int _max_neigh, _N_neigh;

	/// number of boxes the particle has diffused in each direction
	int _pos_shift[3];

	void _check();

public:
	int N_ext_forces;
	BaseForce<number> *ext_forces[MAX_EXT_FORCES];
	int next_particle;
	int strand_id;
	int N_int_centers;

	BaseParticle();
	virtual ~BaseParticle();

//	std::vector<std::pair<BaseParticle<number> *, BaseParticle<number> *> > affected;
	std::vector<ParticlePair<number> > affected;

	virtual void set_positions() { }

	virtual void copy_from(const BaseParticle<number> &);
	inline void soft_copy_from(const BaseParticle<number> * p) {
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
	void is_neighbour(const BaseParticle &p);

	int get_index() const { return index; }

	/**
	 * @brief Add an external force.
	 *
	 * @param f
	 * @return true if the force was added, false otherwise
	 */
	bool add_ext_force(BaseForce<number> *f);
	void set_initial_forces(llint step, number box) {
		LR_vector<number> abs_pos = get_abs_pos(box);
		this->force = LR_vector<number>(0, 0, 0);
		for(int i = 0; i < N_ext_forces; i++) {
			this->force += ext_forces[i]->value(step, abs_pos);
		}
	}

	/**
	 * @brief Add an external potential.
	 *
	 * @param f
	 * @return true if the external potential was added, false otherwise
	 */
	bool add_ext_potential (BaseForce<number> *f);

	/**
	 * @brief Add an external potential.
	 *
	 * @param step current time step. Useful for forces that depend on time.
	 * @param box the current box size. Useful for PBC-aware forces
	 * @return true if the external potential was added, false otherwise
	 */
	void set_ext_potential (llint step, number box) {
		LR_vector<number> abs_pos = get_abs_pos(box);
		this->ext_potential = (number) 0.;
		for(int i = 0; i < N_ext_forces; i++) {
			this->ext_potential += ext_forces[i]->potential(step, abs_pos);
		}
	}

	/**
	 * @brief Checks whether q and the current particle are bonded neighbours (such as neighbouring particles on a DNA strand).
	 *
	 * @param q candidate bonded neighbour
	 * @return true if the current particle and q are bonded neighbours, false otherwise
	 */
	virtual bool is_bonded(BaseParticle<number> *q) {
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

	/**
	 * @brief Shifts the particle's position and stores internally the shift. Used in the fix_diffusion procedure
	 *
	 * @param my_shift reference vector to put back in the box (i.e., strand c.o.m.)
	 * @param box
	 */
	inline void shift(LR_vector<number> &my_shift, LR_vector<number> &box_sides) {
		_pos_shift[0] += (int) floor(my_shift.x / box_sides.x);
		_pos_shift[1] += (int) floor(my_shift.y / box_sides.y);
		_pos_shift[2] += (int) floor(my_shift.z / box_sides.z);
		pos.x -= box_sides.x * floor(my_shift.x / box_sides.x);
		pos.y -= box_sides.y * floor(my_shift.y / box_sides.y);
		pos.z -= box_sides.z * floor(my_shift.z / box_sides.z);
	}
	
	inline void set_pos_shift (int x, int y, int z) {
		_pos_shift[0] = x;
		_pos_shift[1] = y;
		_pos_shift[2] = z;
	}
	
	void get_pos_shift (int * arg) {
		arg[0] = _pos_shift[0];
		arg[1] = _pos_shift[1];
		arg[2] = _pos_shift[2];
	}

	/// Returns the absolute position of the particle, useful for MSD and such
	LR_vector<number> get_abs_pos(number box) { return pos + box * LR_vector<number> ((number)_pos_shift[0], (number)_pos_shift[1], (number)_pos_shift[2]); }
	LR_vector<number> get_abs_pos(number box_x, number box_y, number box_z) { return pos + LR_vector<number> (box_x * (number)_pos_shift[0], box_y * (number)_pos_shift[1], box_z * (number)_pos_shift[2]); }
	LR_vector<number> get_abs_pos(LR_vector<number> box_sides) { return pos + LR_vector<number> (box_sides.x * (number)_pos_shift[0], box_sides.y * (number)_pos_shift[1], box_sides.z * (number)_pos_shift[2]); }

	/// Index of the particle. Usually it is a useful way of accessing arrays of particles
	int index;

	/// Particle type (the meaning of which depends on the chosen interaction type)
	int type;

	/// Needed for specific base pairing
	int btype;

	/// DNA bonded neighbours
	BaseParticle<number> *n3, *n5;

	LR_vector<number> pos_list;
	/// Torque exerted on the particle in its own reference system
	LR_vector<number> torque;

	/// Force exerted on the particle
	LR_vector<number> force;

	/// Total potential energy due to external forces
	number ext_potential;

	/// Particle position inside the box
	LR_vector<number> pos;

	/// Positions of all interaction centers. This array must be initialized by child classes
	LR_vector<number> *int_centers;

	/// Angular momentum of the particle
	LR_vector<number> L;

	/// Velocity of the particle
	LR_vector<number> vel;

	/// BaseParticle orientational matrix
	LR_matrix<number> orientation;

	/// transpose (= inverse) orientational matrix
	LR_matrix<number> orientationT;
};

/*
 * helpers
 */
template <typename number>
class ParticlePair {
public:
	BaseParticle<number> *first;
	BaseParticle<number> *second;

	ParticlePair (BaseParticle<number> *p, BaseParticle<number> *q) {
		if(p == q) throw oxDNAException("ParticlePair: p == q is not allowed");
		if (p->index < q->index) {
			first = p;
			second = q;
		}
		else {
			first = q;
			second = p;
		}
	}

	bool operator< (ParticlePair q) const {
		int p1 = first->index;
		int p2 = second->index;
		int q1 = q.first->index;
		int q2 = q.second->index;

		if(p1 == q1) return (p2 < q2);
		else return (p1 < q1);
	}
};

#endif /* BASEPARTICLE_H_ */
