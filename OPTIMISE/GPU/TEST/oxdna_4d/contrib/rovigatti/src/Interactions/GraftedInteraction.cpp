/*
 * GraftedInteraction.cpp
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#include "GraftedInteraction.h"
#include <string>

using namespace std;

GraftedInteraction::GraftedInteraction() :
				BaseInteraction() {
	ADD_INTERACTION_TO_MAP(0, pair_interaction);

	_N_arms = _N_per_arm = -1;
	_colloid_n = -1;
	_walls = false;
}

GraftedInteraction::~GraftedInteraction() {

}

void GraftedInteraction::get_settings(input_file &inp) {
	BaseInteraction::get_settings(inp);
	_TSP_inter.get_settings(inp);

	getInputNumber(&inp, "GRF_alpha", &_alpha, 1);
	getInputNumber(&inp, "GRF_colloid_sigma", &_colloid_sigma, 1);
	getInputInt(&inp, "TSP_n", &_colloid_n, 1);
	getInputNumber(&inp, "TSP_rfene", &_colloid_rfene, 1);
	getInputString(&inp, "GRF_anchor_file", _anchor_file, 1);
	getInputBool(&inp, "GRF_walls", &_walls, 0);
	if(_walls) getInputNumber(&inp, "GRF_wall_distance", &_wall_distance, 1);

	bool fix_diffusion = true;
	getInputBool(&inp, "fix_diffusion", &fix_diffusion, 0);
	if(fix_diffusion) throw oxDNAException("GraftedInteraction is not compatible with fix_diffusion");
}

void GraftedInteraction::init() {
	_TSP_inter.init();

	_colloid_monomer_sqr_sigma = SQR(0.5 * (_colloid_sigma + 1.));
	_colloid_sqr_rfene = SQR(_colloid_rfene);

	number rep_rcut = pow(2., 1. / _colloid_n);
	_colloid_monomer_sqr_rep_rcut = _colloid_monomer_sqr_sigma * SQR(rep_rcut);
}

void GraftedInteraction::set_box(BaseBox *box) {
	BaseInteraction::set_box(box);
	_TSP_inter.set_box(box);

	if(_walls && box->box_sides().x < 2 * _wall_distance) throw oxDNAException("The box side (%lf) is too small to contain walls separated by %lf", box->box_sides().x, _wall_distance);
}

void GraftedInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
	// we need to parse the file storing the anchors' positions
	FILE *af = fopen(_anchor_file, "r");
	if(af == NULL) throw oxDNAException("The anchor file '%s' is not readable", _anchor_file);
	std::vector<LR_vector> anchors;
	double ax, ay, az;
	while(fscanf(af, "%lf %lf %lf\n", &ax, &ay, &az) == 3) {
		LR_vector apos(ax, ay, az);
		// we put the anchor point's position exactly on the colloid's surface
		apos.normalize();
		apos *= _colloid_sigma * 0.5;
		anchors.push_back(apos);
	}
	fclose(af);

	if((int) anchors.size() != _N_arms) throw oxDNAException("The anchor file '%s' contains only %d positions while, according to the topology, there should be %d", _anchor_file, anchors.size(), _N_arms);

	particles[0] = new Colloid(anchors);

	for(int i = 0; i < _N_arms * _N_per_arm; i++) {
		particles[i + 1] = new TSPParticle();
	}
}

void GraftedInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	int res = sscanf(line, "%*d %d %d\n", &_N_arms, &_N_per_arm);
	if(res != 2) throw oxDNAException("Incorrect topology");
	*N_strands = _N_arms + 1;
	if((_N_arms * _N_per_arm + 1) != N) {
		throw oxDNAException("Incoherent topology: the total number of particles should be equal to the number of arms (%d) times the number of particles in an arm (%d) plus one", _N_arms, _N_per_arm);
	}

	allocate_particles(particles);

	particles[0]->type = P_COLLOID;
	particles[0]->index = 0;
	particles[0]->strand_id = 0;

	int attractive_from = (int) round(_N_per_arm * (1. - _alpha));
	for(int i = 1; i < N; i++) {
		TSPParticle *p = (TSPParticle *) particles[i];
		p->index = i;
		p->strand_id = (i - 1) / _N_per_arm + 1;
		p->set_arm((i - 1) / _N_per_arm);

		int rel_index = (i - 1) % _N_per_arm;
		p->type = (rel_index >= attractive_from) ? P_B : P_A;

		p->add_bonded_neigh((TSPParticle *) particles[0]);
		if(rel_index != 0) {
			TSPParticle *p_n = (TSPParticle *) particles[i - 1];
			p->add_bonded_neigh(p_n);
		}
	}
}

number GraftedInteraction::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	else {
		return pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
}

number GraftedInteraction::_wall_interaction(BaseParticle *p, bool update_forces) {
	number energy = 0.;

	number z = this->_box->get_abs_pos(p).z;
	// we exert a constant force on monomers that are on the wrong side of the wall
	// this is useful to compress initial configurations
	if(z <= -(_wall_distance - 0.5) && update_forces) p->force.z += 10.;
	else if(z >= (_wall_distance - 0.5) && update_forces) p->force.z -= 10.;
	else {
		number dz = (z > 0.) ? _wall_distance - z : _wall_distance + z;
		energy += exp(-dz) / dz;
		if(update_forces) {
			number mod_force = energy * (1. + 1. / dz);
			p->force.z += (z > 0.) ? -mod_force : mod_force;
		}
	}

	return energy;
}

number GraftedInteraction::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	number energy = 0.f;
	if(q == P_VIRTUAL) {
		if(_walls) {
			energy += _wall_interaction(p, update_forces);
		}
		TSPParticle *TSPp = (TSPParticle *) p;
		for(typename set<TSPParticle *>::iterator it = TSPp->bonded_neighs.begin(); it != TSPp->bonded_neighs.end(); it++) {
			energy += pair_interaction_bonded(p, *it, true, update_forces);
		}
		return energy;
	}
	else if(p->type != P_COLLOID && q->type != P_COLLOID) return _TSP_inter.pair_interaction_bonded(p, q, compute_r, update_forces);

	// if we are here then we have to compute the interaction between a regular monomer and the colloid.
	// we have two possibilities: either the monomer is at the beginning of the chain or it's not.
	// if it's at the beginning of the chain then we also have to add the FENE
	BaseParticle *monomer = (p->type == P_COLLOID) ? q : p;
	BaseParticle *colloid = (p->type == P_COLLOID) ? p : q;
	if((monomer->index - 1) % _N_per_arm == 0) {
		int anchor_id = monomer->strand_id - 1;
		LR_vector apos_rel = colloid->int_centers[anchor_id];
		LR_vector apos = colloid->pos + apos_rel;
		LR_vector r_ma = this->_box->min_image(monomer->pos, apos);

		number sqr_r_ma = r_ma.norm();

		if(sqr_r_ma > _colloid_sqr_rfene) {
			if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r_ma), sqrt(_colloid_sqr_rfene));
			else return 1e10;
		}

		energy += -15. * _colloid_sqr_rfene * log(1. - sqr_r_ma / _colloid_sqr_rfene);

		if(update_forces) {
			LR_vector tmp_force = r_ma * (-30. * _colloid_sqr_rfene / (_colloid_sqr_rfene - sqr_r_ma));

			monomer->force -= tmp_force;
			colloid->force += tmp_force;

			colloid->torque += colloid->orientationT * apos_rel.cross(tmp_force);
		}
	}

	if(compute_r) {
		if(q != P_VIRTUAL && p != P_VIRTUAL) {
			_computed_r = q->pos - p->pos;
		}
	}

	number sqr_r = _computed_r.norm();
	// the repulsion is cut off at the minimum
	if(sqr_r < _colloid_monomer_sqr_rep_rcut) {
		number part = pow(_colloid_monomer_sqr_sigma / sqr_r, _colloid_n / 2.);
		energy += 4 * (part * (part - 1.)) + 1.;
		if(update_forces) {
			number force_mod = 4 * _colloid_n * part * (2 * part - 1) / sqr_r;
			p->force -= _computed_r * force_mod;
			q->force += _computed_r * force_mod;
		}
	}

	return energy;
}

number GraftedInteraction::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(p->type != P_COLLOID && q->type != P_COLLOID) {
		return _TSP_inter.pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}
	return 0.f;
}

void GraftedInteraction::check_input_sanity(std::vector<BaseParticle *> &particles) {

}

void GraftedInteraction::generate_random_configuration(std::vector<BaseParticle *> &particles) {
	BaseParticle *colloid = particles[0];
	colloid->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
	colloid->orientation.orthonormalize();
	colloid->orientationT = colloid->orientation.get_transpose();
	colloid->set_positions();

	colloid->pos = LR_vector(0, 0, 0);
	int p_index = 1;
	for(int na = 0; na < _N_arms; na++) {
		LR_vector apos_dist = colloid->int_centers[na];
		apos_dist.normalize();
		apos_dist *= 1.01;

		LR_vector current = colloid->pos + colloid->int_centers[na];
		int i = 0;
		do {
			BaseParticle *p = particles[p_index];
			current += apos_dist;
			p->pos = current;
			p->orientation = Utils::get_random_rotation_matrix_from_angle(M_PI);
			p->orientation.orthonormalize();

			i++;
			p_index++;
		} while(i < _N_per_arm);
	}
}

Colloid::Colloid(std::vector<LR_vector> &anchor_poss) :
				TSPParticle() {
	this->int_centers.resize(anchor_poss.size());

	_base_anchors = anchor_poss;
}

Colloid::~Colloid() {

}

void Colloid::set_positions() {
	for(uint i = 0; i < N_int_centers(); i++)
		this->int_centers[i] = this->orientation * _base_anchors[i];
}
