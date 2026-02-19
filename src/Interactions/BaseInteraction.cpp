/*
 * BaseInteraction.cpp
 *
 *  Created on: 11 ott 2019
 *      Author: lorenzo
 */

#include "BaseInteraction.h"

BaseInteraction::BaseInteraction() {
	_energy_threshold = (number) 100.f;
	_is_infinite = false;
	_box = NULL;
	_generate_consider_bonded_interactions = true;
	_generate_bonded_cutoff = 2.0;
	_temperature = 1.0;

	_rcut = 2.5;
	_sqr_rcut = SQR(_rcut);
}

BaseInteraction::~BaseInteraction() {

}

number BaseInteraction::pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		_computed_r = _box->min_image(p->pos, q->pos);
	}

	// _int_map.find(name) returns an iterator which points to the
	// requested function pointer. We dereference it, apply it to the
	// that, which is a pointer to an instance of a class which inherits
	// from BaseInteraction. We then pass in the required parameters and
	// return the interaction energy.
	typename interaction_map::iterator interaction = _interaction_map.find(name);
	if(interaction == _interaction_map.end()) {
		throw oxDNAException("%s, line %d: Interaction term '%d' not found", __FILE__, __LINE__, name);
	}

	return interaction->second(p, q, false, update_forces);
}

std::map<int, number> BaseInteraction::get_system_energy_split(std::vector<BaseParticle *> &particles, BaseList *lists) {
	begin_energy_computation();

	std::map<int, number> energy_map;

	for(auto it = _interaction_map.begin(); it != _interaction_map.end(); it++) {
		int name = it->first;
		energy_map[name] = (number) 0.f;
	}

	for(auto p: particles) {
		std::vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			if(p->index > q->index) {
				for(auto it = _interaction_map.begin(); it != _interaction_map.end(); it++) {
					int name = it->first;
					energy_map[name] += pair_interaction_term(name, p, q);
				}
			}
		}
	}

	return energy_map;
}

void BaseInteraction::get_settings(input_file &inp) {
	getInputString(&inp, "topology", _topology_filename, 1);
	getInputNumber(&inp, "energy_threshold", &_energy_threshold, 0);
	getInputNumber(&inp, "T", &_temperature, 1);
	getInputBool(&inp, "generate_consider_bonded_interactions", &_generate_consider_bonded_interactions, 0);
	if(_generate_consider_bonded_interactions) {
		getInputNumber(&inp, "generate_bonded_cutoff", &_generate_bonded_cutoff, 0);
		OX_LOG(Logger::LOG_INFO, "The generator will try to take into account bonded interactions by choosing distances between bonded neighbours no larger than %lf", _generate_bonded_cutoff);
	}
}

int BaseInteraction::get_N_from_topology() {
	std::ifstream topology;
	topology.open(_topology_filename, std::ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", _topology_filename);
	int ret;
	topology >> ret;
	topology.close();
	return ret;
}

void BaseInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	*N_strands = particles.size();
	allocate_particles(particles);
	int idx = 0;
	for(auto p : particles) {
		p->index = idx;
		p->type = 0;
		p->strand_id = idx;
		idx++;
	}
}

void BaseInteraction::begin_energy_computation() {

}

void BaseInteraction::begin_energy_and_force_computation() {
	if(has_custom_stress_tensor()) {
		reset_stress_tensor();
	}
	begin_energy_computation();
}

void BaseInteraction::_update_stress_tensor(const LR_vector &r_p, const LR_vector &group_force) {
	_stress_tensor[0] += r_p.x * group_force.x;
	_stress_tensor[1] += r_p.y * group_force.y;
	_stress_tensor[2] += r_p.z * group_force.z;
	_stress_tensor[3] += r_p.x * group_force.y;
	_stress_tensor[4] += r_p.x * group_force.z;
	_stress_tensor[5] += r_p.y * group_force.z;
}

void BaseInteraction::compute_standard_stress_tensor() {
	static std::vector<LR_vector> old_forces, old_torques;

	begin_energy_and_force_computation();

	// pair_interaction will change these vectors, but we still need them in the next
	// first integration step. For this reason we copy and then restore their values
	// after the calculation
	old_forces.resize(CONFIG_INFO->N());
	old_torques.resize(CONFIG_INFO->N());
	for(int i = 0; i < CONFIG_INFO->N(); i++) {
		old_forces[i] = CONFIG_INFO->particles()[i]->force;
		old_torques[i] = CONFIG_INFO->particles()[i]->torque;
	}

	StressTensor stress_tensor = { 0., 0., 0., 0., 0., 0. };
	double energy = 0.;

	for(auto p : CONFIG_INFO->particles()) {
		std::vector<BaseParticle *> neighs = CONFIG_INFO->lists->get_neigh_list(p);

		std::set<BaseParticle *> bonded_neighs;
		for(auto &pair : p->affected) {
			if(pair.first != p) {
				bonded_neighs.insert(pair.first);
			}
			if(pair.second != p) {
				bonded_neighs.insert(pair.second);
			}
		}

		neighs.insert(neighs.end(), bonded_neighs.begin(), bonded_neighs.end());

		for(auto q : neighs) {
			if(p->index > q->index) {
				LR_vector r = CONFIG_INFO->box->min_image(p->pos, q->pos);
				set_computed_r(r);

				p->force = LR_vector();
				energy += (double) pair_interaction(p, q, false, true);

				stress_tensor[0] -= r.x * p->force.x;
				stress_tensor[1] -= r.y * p->force.y;
				stress_tensor[2] -= r.z * p->force.z;
				stress_tensor[3] -= r.x * p->force.y;
				stress_tensor[4] -= r.x * p->force.z;
				stress_tensor[5] -= r.y * p->force.z;
			}
		}

		LR_vector &vel = p->vel;
		stress_tensor[0] += SQR(vel.x);
		stress_tensor[1] += SQR(vel.y);
		stress_tensor[2] += SQR(vel.z);
		stress_tensor[3] += vel.x * vel.y;
		stress_tensor[4] += vel.x * vel.z;
		stress_tensor[5] += vel.y * vel.z;
	}

	for(int i = 0; i < CONFIG_INFO->N(); i++) {
		auto p = CONFIG_INFO->particles()[i];

		p->force = old_forces[i];
		p->torque = old_torques[i];
	}

	_stress_tensor = stress_tensor;
}

void BaseInteraction::reset_stress_tensor() {
	std::fill(_stress_tensor.begin(), _stress_tensor.end(), 0.);
	for(auto p : CONFIG_INFO->particles()) {
		_stress_tensor[0] += SQR(p->vel.x);
		_stress_tensor[1] += SQR(p->vel.y);
		_stress_tensor[2] += SQR(p->vel.z);
		_stress_tensor[3] += p->vel.x * p->vel.y;
		_stress_tensor[4] += p->vel.x * p->vel.z;
		_stress_tensor[5] += p->vel.y * p->vel.z;
	}
}

StressTensor BaseInteraction::stress_tensor() const {
    StressTensor norm_st(_stress_tensor);
    for(auto &v : norm_st) {
        v /= CONFIG_INFO->box->V();
    }
    return norm_st;
}  // <-- THIS brace must be here, before the next function starts


StressTensor BaseInteraction::stress_tensor_subset(const std::vector<BaseParticle *>& subset,
                                                  const LR_vector& reference_point,
                                                  number volume,
                                                  bool use_min_image,
                                                  bool include_kinetic) const
	{
    StressTensor st = {0., 0., 0., 0., 0., 0.};

    if(volume <= (number)0.) return st;

    BaseBox *box = CONFIG_INFO->box;

    for(const auto *p : subset) {
        // Position relative to reference point
        LR_vector dr;
        if(use_min_image && box != NULL) {
            dr = box->min_image(p->pos, reference_point);
        } else {
            dr = p->pos - reference_point;
        }

        // Total force already computed on particle
        const LR_vector &F = p->force;

        // Virial part. Sign chosen to match the convention used in compute_standard_stress_tensor()
        // (which accumulates - r_ij \otimes F_ij).
        st[0] -= dr.x * F.x;
        st[1] -= dr.y * F.y;
        st[2] -= dr.z * F.z;
        st[3] -= dr.x * F.y;
        st[4] -= dr.x * F.z;
        st[5] -= dr.y * F.z;

        if(include_kinetic) {
            const LR_vector &v = p->vel;
            st[0] += SQR(v.x);
            st[1] += SQR(v.y);
            st[2] += SQR(v.z);
            st[3] += v.x * v.y;
            st[4] += v.x * v.z;
            st[5] += v.y * v.z;
        }
    }

    // Normalise by volume
    for(auto &x : st) x /= volume;
    return st;
}

StressTensor BaseInteraction::stress_tensor_subset_cylindrical(const std::vector<BaseParticle *>& subset,
                                                               const LR_vector& reference_point,
                                                               number volume,
                                                               bool use_min_image,
                                                               bool include_kinetic) const
{
    StressTensor st = {0., 0., 0., 0., 0., 0.}; // rr, tt, zz, rt, rz, tz
    if(volume <= (number)0.) return st;

    BaseBox *box = CONFIG_INFO->box;

    for(const auto *p : subset) {
        // Position relative to reference point (axis point)
        LR_vector dr;
        if(use_min_image && box != NULL) {
            dr = box->min_image(p->pos, reference_point);
        } else {
            dr = p->pos - reference_point;
        }

        const LR_vector &F = p->force;

        // Build local cylindrical basis around z-axis:
        // er = (dx,dy,0)/r; et = (-er_y, er_x, 0); ez = (0,0,1)
        const number rx = dr.x;
        const number ry = dr.y;
        const number r2 = rx*rx + ry*ry;

        // If a particle lies exactly on the axis, cylindrical basis is undefined.
        // Skip it (or you could fall back to Cartesian there).
        if(r2 <= (number)1e-20) continue;

        const number rinv = (number)1.0 / std::sqrt(r2);

        LR_vector er(rx * rinv, ry * rinv, (number)0.0);
        LR_vector et(-er.y, er.x, (number)0.0);
        LR_vector ez((number)0.0, (number)0.0, (number)1.0);

        // Projections of dr and F onto cylindrical basis
        const number r_r  = dr * er;
        const number r_t  = dr * et;
        const number r_z  = dr * ez;

        const number f_r  = F * er;
        const number f_t  = F * et;
        const number f_z  = F * ez;

        // Virial contribution: - (dr ⊗ F) expressed in cylindrical basis:
        // sigma_ab += - ( (dr·ea) * (F·eb) )
        st[0] -= r_r * f_r; // rr
        st[1] -= r_t * f_t; // tt
        st[2] -= r_z * f_z; // zz

        // Off-diagonals (keep same 6-entry convention)
        st[3] -= r_r * f_t; // rt
        st[4] -= r_r * f_z; // rz
        st[5] -= r_t * f_z; // tz

        if(include_kinetic) {
            const LR_vector &v = p->vel;
            const number v_r = v * er;
            const number v_t = v * et;
            const number v_z = v * ez;

            st[0] += v_r * v_r;
            st[1] += v_t * v_t;
            st[2] += v_z * v_z;
            st[3] += v_r * v_t;
            st[4] += v_r * v_z;
            st[5] += v_t * v_z;
        }
    }

    // Normalise by volume
    for(auto &x : st) x /= volume;

    return st;
}

number BaseInteraction::get_system_energy(std::vector<BaseParticle *> &particles, BaseList *lists) {
	begin_energy_computation();

	double energy = 0.;
	std::vector<ParticlePair> pairs = lists->get_potential_interactions();
	for(auto &pair : pairs) {
		BaseParticle *p = pair.first;
		BaseParticle *q = pair.second;
		energy += (double) pair_interaction(p, q);
		if(get_is_infinite()) {
			return energy;
		}
	}

	return (number) energy;
}

number BaseInteraction::get_system_energy_term(int name, std::vector<BaseParticle *> &particles, BaseList *lists) {
	begin_energy_computation();

	number energy = (number) 0.f;
	for(auto p : particles) {
		std::vector<BaseParticle *> neighs = lists->get_all_neighbours(p);

		for(unsigned int n = 0; n < neighs.size(); n++) {
			BaseParticle *q = neighs[n];
			if(p->index > q->index) energy += pair_interaction_term(name, p, q);
			if(get_is_infinite()) return energy;
		}
	}

	return energy;
}

bool BaseInteraction::generate_random_configuration_overlap(BaseParticle *p, BaseParticle *q) {
	_computed_r = _box->min_image(p, q);

	if(_computed_r.norm() >= _sqr_rcut) {
		return false;
	}

	number energy = pair_interaction_nonbonded(p, q, false, false);

	// in case we create an overlap, we reset the interaction state
	set_is_infinite(false);

	// if energy too large, reject
	return energy > _energy_threshold;
}

void BaseInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	Cells c(particles, _box);
	c.init(_rcut);

	for(auto p : particles) {
		p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y, drand48() * _box->box_sides().z);
	}

	c.global_update();

	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		bool same_strand = (_generate_consider_bonded_interactions && i > 0 && p->is_bonded(particles[i - 1]));

		bool inserted = false;
		do {
			if(same_strand) {
				p->pos = particles[i - 1]->pos + LR_vector((drand48() - 0.5), (drand48() - 0.5), (drand48() - 0.5)) * _generate_bonded_cutoff;
			}
			else {
				p->pos = LR_vector(drand48() * _box->box_sides().x, drand48() * _box->box_sides().y, drand48() * _box->box_sides().z);
			}
			// random orientation
			//p->orientation = Utils::get_random_rotation_matrix (2.*M_PI);
			p->orientation = Utils::get_random_rotation_matrix_from_angle(acos(2. * (drand48() - 0.5)));
			p->orientation.orthonormalize();
			p->orientationT = p->orientation.get_transpose();

			p->set_positions();
			p->set_ext_potential(0, _box);
			c.single_update(p);

			inserted = true;

			// we take into account the bonded neighbours
			for(unsigned int n = 0; n < p->affected.size(); n++) {
				BaseParticle *p1 = p->affected[n].first;
				BaseParticle *p2 = p->affected[n].second;
				if(p1->index <= p->index && p2->index <= p->index) {
					number e = pair_interaction_bonded(p1, p2);
					if(std::isnan(e) || e > _energy_threshold) inserted = false;
				}
			}

			// here we take into account the non-bonded interactions
			std::vector<BaseParticle *> neighs = c.get_complete_neigh_list(p);
			for(unsigned int n = 0; n < neighs.size(); n++) {
				BaseParticle *q = neighs[n];
				// particles with an index larger than p->index have not been inserted yet
				if(q->index < p->index && generate_random_configuration_overlap(p, q)) inserted = false;
			}

			// we take into account the external potential
			number boltzmann_factor = exp(-p->ext_potential / _temperature);
			if(std::isnan(p->ext_potential) || drand48() > boltzmann_factor) {
				inserted = false;
			}

		} while(!inserted);

		if(i > 0 && N > 10 && i % (N / 10) == 0) OX_LOG(Logger::LOG_INFO, "Inserted %d%% of the particles (%d/%d)", i*100/N, i, N);
	}
}