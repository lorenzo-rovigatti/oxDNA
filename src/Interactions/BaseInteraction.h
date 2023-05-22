/*
 * BaseInteraction is a class that is meant "contain" all the interaction
 * types
 */

#ifndef BASE_INTERACTION_H
#define BASE_INTERACTION_H

#include "../defs.h"
#include "../Particles/BaseParticle.h"
#include "../Boxes/BaseBox.h"
#include "../Lists/BaseList.h"
#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"
#include "Mesh.h"

#include "../Lists/Cells.h"

#include <map>
#include <fstream>
#include <set>
#include <vector>
#include <functional>

#define ADD_INTERACTION_TO_MAP(index, member) {_interaction_map[index] = [this](BaseParticle *p, BaseParticle *q, bool compute_r, bool compute_forces) { return member(p, q, compute_r, compute_forces); };}

/**
 * @brief Base class for managing particle-particle interactions. It is an abstract class.
 */
class BaseInteraction {
private:

protected:
	BaseBox *_box;

	/// This is useful for "hard" potentials
	bool _is_infinite;

	/// This are needed to create initial configurations, and they are used by the generator functions
	number _energy_threshold, _temperature;
	/// If true, the generation of the initial configuration will try to take into account the fact that particles that are bonded neighbours should be close to each other
	bool _generate_consider_bonded_interactions;
	/// This controls the maximum at which bonded neighbours should be randomly placed to speed-up generation. Used by generator functions.
	number _generate_bonded_cutoff;

	number _rcut, _sqr_rcut;

	char _topology_filename[256];

	LR_vector _computed_r;

	StressTensor _stress_tensor;

	virtual void _update_stress_tensor(const LR_vector &r_p, const LR_vector &group_force);

	using energy_function = std::function<number(BaseParticle *, BaseParticle *, bool, bool)>;
	using interaction_map = std::map<int, energy_function>;
	interaction_map _interaction_map;

public:
	/**
	 * @brief Basic constructor. By default, it does not need anything.
	 */
	BaseInteraction();
	virtual ~BaseInteraction();

	/**
	 * @brief Computes all the contributions to the total potential energy and returns them in a map.
	 *
	 * @param particles
	 * @param N
	 * @return a map storing all the potential energy contributions.
	 */
	virtual std::map<int, number> get_system_energy_split(std::vector<BaseParticle *> &particles, BaseList *lists);

	virtual void set_box(BaseBox *box) {
		_box = box;
	}

	virtual void get_settings(input_file &inp);

	/**
	 * Initialization of class constants.
	 */
	virtual void init() = 0;

	/**
	 * @brief Handles particle allocation. Child classes must implement it.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void allocate_particles(std::vector<BaseParticle *> &particles) = 0;

	/**
	 * @brief Returns the maximum cut-off associated to the interaction
	 *
	 * @return Interaction cutoff
	 */
	virtual number get_rcut() const {
		return _rcut;
	}

	/**
	 * @brief Check whether the initial configuration makes sense.
	 *
	 * Since this operation is interaction-dependent, each interaction must provide its own implementation.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void check_input_sanity(std::vector<BaseParticle *> &particles) = 0;

	/**
	 * @brief Signals the interaction that an energy computation is about to begin.
	 *
	 * By default this method does nothing, but interactions inheriting from this interface may
	 * need to initialise or reset other data structures before computing the energy or the force acting on all particles.
	 *
	 */
	virtual void begin_energy_computation();

	/**
	 * @brief Signals the interaction that an energy and force computation is about to begin.
	 *
	 * By default this method resets the stress tensor data structures, calls begin_energy_computationa() and nothing else, but interactions inheriting from this interface may
	 * need to initialise or reset other data structures before computing the energy or the force acting on all particles.
	 *
	 */
	virtual void begin_energy_and_force_computation();

	/**
	 * @brief Returns true if the interaction computes the stress tensor internally, false otherwise.
	 *
	 * This method will return false if not overridden.
	 *
	 * @return
	 */
	virtual bool has_custom_stress_tensor() const {
		return false;
	}

	void reset_stress_tensor();

	void compute_standard_stress_tensor();

	StressTensor stress_tensor() const;

	void set_stress_tensor(StressTensor st) {
		_stress_tensor = st;
	}

	/**
	 * @brief Computes the total interaction between particles p and q.
	 *
	 * If r is not given or NULL, it is computed from scratch. It can optionally update forces and torques exerted on p and q.
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return pair-interaction energy
	 */
	virtual number pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) = 0;

	/**
	 * @brief Computes the bonded part interaction between particles p and q.
	 *
	 * See {\@link pair_interaction} for a description of the parameters.
	 */
	virtual number pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) = 0;

	/**
	 * @brief Computed the non bonded part of the interaction between particles p and q.
	 *
	 * See {\@link pair_interaction} for a description of the parameters.
	 */
	virtual number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) = 0;

	/**
	 * @brief Computes the requested term of the interaction energy between p and q.
	 *
	 * @param name identifier of the interaction method
	 * @param p
	 * @param q
	 * @param r
	 * @param update_forces
	 * @return
	 */
	virtual number pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false);

	/**
	 * @brief Returns the total potential energy of the system
	 *
	 * @param particles
	 * @param N
	 * @return
	 */
	virtual number get_system_energy(std::vector<BaseParticle *> &particles, BaseList *lists);

	/**
	 * @brief Returns the given energy term computed on all the system
	 *
	 * @param name
	 * @param particles
	 * @param N
	 * @return the required energy contribution
	 */
	virtual number get_system_energy_term(int name, std::vector<BaseParticle *> &particles, BaseList *lists);

	/**
	 * @brief Read the topology of the interactions. The defaut
	 * method is empty.
	 *
	 * This function is set up so that the topology of the interactions
	 * are read in the code. The different interactions might require
	 * different formats, so it is implemented here. The default
	 *
	 * @param N_strands Pointer to an int which will be filled with the
	 * number of "strands" (unbreakable clusters of connected particles)
	 * @param particles array of particles.
	 */
	virtual void read_topology(int *N_strands, std::vector<BaseParticle *> &particles);

	/**
	 * @brief Returns the number of particles, as written in the topology.
	 *
	 * @return number of particles
	 */
	virtual int get_N_from_topology();

	/**
	 * @brief Returns the state of the interaction
	 */
	bool get_is_infinite() {
		return _is_infinite;
	}

	/**
	 * @brief Sets _is_infinite
	 *
	 * @param arg bool
	 */
	void set_is_infinite(bool arg) {
		_is_infinite = arg;
	}

	/**
	 * @brief Sets the distance between particles to be used by the pair_interaction_* methods.
	 *
	 * @param r the new particle distance
	 */
	void set_computed_r(LR_vector &r) {
		_computed_r =r;
	}

	/**
	 * @brief Check whether the two particles overlaps. Used only in the generation of configurations. Can be overloaded.
	 *
	 * The default implementation does not allow two particles to have an
	 * interaction strength greater than 100. More subtle implementations
	 * may be needed in special cases.
	 *
	 * @param p
	 * @param q
	 * @return bool whether the two particles overlap
	 */
	virtual bool generate_random_configuration_overlap(BaseParticle *p, BaseParticle *q);

	/**
	 * @brief Generate an initial configuration. Can be overloaded.
	 *
	 * The default function creates a random configuration in the most simple way possible:
	 * puts particles in one at a time, and using {\@link generate_random_configuration_overlap}
	 * to test if two particles overlap.
	 *
	 * @param particles
	 * @param N
	 */
	virtual void generate_random_configuration(std::vector<BaseParticle *> &particles);
};

using InteractionPtr = std::shared_ptr<BaseInteraction>;

#endif /* BASE_INTERACTION_H */
