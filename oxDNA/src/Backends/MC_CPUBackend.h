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

/**
 * @brief Manages a MC simulation on CPU. It supports NVT and NPT simulations
 */
template<typename number>
class MC_CPUBackend: public MCBackend<number> {
protected:
	number _verlet_skin;
	BaseParticle<number> **_particles_old;

	bool _enable_flip;

	Timer * _timer_move;
	Timer * _timer_box;
	Timer * _timer_lists;

	number _target_box;
	number _box_tolerance;
	number _e_tolerance;

	inline number _excluded_volume(const LR_vector<number> &r, number sigma, number rstar, number b, number rc);
	void _compute_energy();
	inline void _translate_particle(BaseParticle<number> *p);
	inline void _rotate_particle(BaseParticle<number> *p);
	inline number _particle_energy(BaseParticle<number> *p, bool reuse=false);

	std::map<ParticlePair<number>, number> _stored_bonded_interactions;
	std::map<ParticlePair<number>, number> _stored_bonded_tmp;
	//std::vector < std:pair< std::pair<BaseParticle<number> *, BaseParticle<number> *>, number > > _stored_bonded_interactions;

public:
	MC_CPUBackend();
	virtual ~MC_CPUBackend();

	virtual void get_settings(input_file &inp);
	void init();

	void sim_step(llint cur_step);
};

#endif /* MC_CPUBACKEND_H_ */
