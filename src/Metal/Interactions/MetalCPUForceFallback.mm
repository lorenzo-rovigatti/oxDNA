/**
 * @file    MetalCPUForceFallback.mm
 * @brief   CPU force fallback helper for Metal interactions
 */

#include "MetalCPUForceFallback.h"

#include "../../Interactions/BaseInteraction.h"
#include "../../Lists/BaseList.h"
#include "../../Particles/BaseParticle.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Utilities/oxDNAException.h"

#include <vector>

namespace {

inline LR_matrix _orientation_from_quaternion(const m_quat &q) {
	number sqx = (number) q.x * (number) q.x;
	number sqy = (number) q.y * (number) q.y;
	number sqz = (number) q.z * (number) q.z;
	number sqw = (number) q.w * (number) q.w;
	number xy = (number) q.x * (number) q.y;
	number xz = (number) q.x * (number) q.z;
	number xw = (number) q.x * (number) q.w;
	number yz = (number) q.y * (number) q.z;
	number yw = (number) q.y * (number) q.w;
	number zw = (number) q.z * (number) q.w;
	number norm = sqx + sqy + sqz + sqw;
	if(norm < 1e-12) {
		return LR_matrix((number) 1., (number) 0., (number) 0.,
		                 (number) 0., (number) 1., (number) 0.,
		                 (number) 0., (number) 0., (number) 1.);
	}
	number invs = 1.0 / norm;

	LR_matrix orientation;
	orientation.v1.x = (sqx - sqy - sqz + sqw) * invs;
	orientation.v1.y = 2 * (xy - zw) * invs;
	orientation.v1.z = 2 * (xz + yw) * invs;
	orientation.v2.x = 2 * (xy + zw) * invs;
	orientation.v2.y = (-sqx + sqy - sqz + sqw) * invs;
	orientation.v2.z = 2 * (yz - xw) * invs;
	orientation.v3.x = 2 * (xz - yw) * invs;
	orientation.v3.y = 2 * (yz + xw) * invs;
	orientation.v3.z = (-sqx - sqy + sqz + sqw) * invs;

	return orientation;
}

}

void MetalCPUForceFallback::compute(int N,
                                    id<MTLBuffer> d_poss,
                                    id<MTLBuffer> d_orientations,
                                    id<MTLBuffer> d_forces,
                                    id<MTLBuffer> d_torques,
                                    id<MTLBuffer> d_energies) {
	if(N <= 0) {
		return;
	}

	auto cfg = ConfigInfo::instance();
	BaseInteraction *interaction = cfg->interaction;
	BaseList *lists = cfg->lists;
	if(interaction == nullptr || lists == nullptr || cfg->box == nullptr) {
		throw oxDNAException("Metal CPU fallback requires ConfigInfo interaction, lists and box");
	}

	auto &particles = cfg->particles();
	if((int) particles.size() != N) {
		throw oxDNAException("Metal CPU fallback particle count mismatch (%d vs %zu)", N, particles.size());
	}

	std::vector<m_number4> h_poss(N);
	std::vector<m_quat> h_orientations(N);
	std::vector<m_number4> h_forces(N);
	std::vector<m_number4> h_torques(N);

	MetalUtils::copy_from_device<m_number4>(h_poss.data(), d_poss, N);
	MetalUtils::copy_from_device<m_quat>(h_orientations.data(), d_orientations, N);

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];

		p->pos.x = h_poss[i].x;
		p->pos.y = h_poss[i].y;
		p->pos.z = h_poss[i].z;

		p->orientation = _orientation_from_quaternion(h_orientations[i]);
		p->orientationT = p->orientation.get_transpose();
		p->set_positions();

		p->set_initial_forces(cfg->curr_step, cfg->box);
	}

	// Mirror the CPU MD backend: refresh the cell lists incrementally and only
	// rebuild the (expensive) Verlet neighbour lists when a particle has moved
	// farther than the Verlet skin. Previously this unconditionally forced a
	// full rebuild on every step, which made the fallback path slower than the
	// plain CPU backend.
	for(auto p : particles) {
		lists->single_update(p);
	}
	if(!lists->is_updated()) {
		lists->global_update();
	}
	interaction->begin_energy_and_force_computation();

	for(auto p : particles) {
		for(auto &pair : p->affected) {
			if(pair.first == p) {
				interaction->pair_interaction_bonded(pair.first, pair.second, true, true);
			}
		}

		for(auto q : lists->get_neigh_list(p)) {
			interaction->pair_interaction_nonbonded(p, q, true, true);
		}
	}

	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		h_forces[i].x = (m_number) p->force.x;
		h_forces[i].y = (m_number) p->force.y;
		h_forces[i].z = (m_number) p->force.z;
		h_forces[i].w = (m_number) 0.f;

		h_torques[i].x = (m_number) p->torque.x;
		h_torques[i].y = (m_number) p->torque.y;
		h_torques[i].z = (m_number) p->torque.z;
		h_torques[i].w = (m_number) 0.f;
	}

	MetalUtils::copy_to_device<m_number4>(d_forces, h_forces.data(), N);
	MetalUtils::copy_to_device<m_number4>(d_torques, h_torques.data(), N);

	if(d_energies != nil) {
		const size_t n_energies = [d_energies length] / sizeof(float);
		if(n_energies > 0) {
			std::vector<float> zeros(n_energies, 0.f);
			MetalUtils::copy_to_device<float>(d_energies, zeros.data(), n_energies);
		}
	}
}
