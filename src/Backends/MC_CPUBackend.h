/**
 * @file    MC_CPUBackend.h
 * @date    26/nov/2010
 * @author  lorenzo
 *
 *
 */

#ifndef MC_CPUBACKEND_H_
#define MC_CPUBACKEND_H_

#include "MCBackend.h"

 class ParticlePair;

/**
 * @brief Manages a MC simulation on CPU. It supports NVT and NPT simulations
 */

class MC_CPUBackend: public MCBackend {
protected:
	number _verlet_skin;
	std::vector<BaseParticle *> _particles_old;

	bool _enable_flip;

	std::shared_ptr<Timer> _timer_move;
	std::shared_ptr<Timer> _timer_box;
	std::shared_ptr<Timer> _timer_lists;

	number _target_box;
	number _box_tolerance;
	number _e_tolerance;

	inline number _excluded_volume(const LR_vector &r, number sigma, number rstar, number b, number rc);
	void _compute_energy();
	inline void _translate_particle(BaseParticle *p);
	inline void _rotate_particle(BaseParticle *p);
	inline number _particle_energy(BaseParticle *p, bool reuse=false);

	std::map<ParticlePair, number> _stored_bonded_interactions;
	std::map<ParticlePair, number> _stored_bonded_tmp;

public:
	MC_CPUBackend();
	virtual ~MC_CPUBackend();

	virtual void get_settings(input_file &inp);
	void init();

	void sim_step();
};

#endif /* MC_CPUBACKEND_H_ */
