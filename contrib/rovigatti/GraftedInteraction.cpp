/*
 * GraftedInteraction.cpp
 *
 *  Created on: 10/feb/2013
 *      Author: lorenzo
 */

#include "GraftedInteraction.h"
#include <string>

using namespace std;

template<typename number>
GraftedInteraction<number>::GraftedInteraction() : BaseInteraction<number, GraftedInteraction<number> >() {
	this->_int_map[0] = &GraftedInteraction<number>::pair_interaction;
	_N_arms = _N_per_arm = -1;
	_colloid_n = -1;
	_walls = false;
}

template<typename number>
GraftedInteraction<number>::~GraftedInteraction() {

}

template<typename number>
void GraftedInteraction<number>::get_settings(input_file &inp) {
	IBaseInteraction<number>::get_settings(inp);
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


template<typename number>
void GraftedInteraction<number>::init() {
	_TSP_inter.init();

	_colloid_monomer_sqr_sigma = SQR(0.5*(_colloid_sigma + 1.));
	_colloid_sqr_rfene = SQR(_colloid_rfene);

	number rep_rcut = pow(2., 1./_colloid_n);
	_colloid_monomer_sqr_rep_rcut = _colloid_monomer_sqr_sigma * SQR(rep_rcut);
}

template<typename number>
void GraftedInteraction<number>::set_box_side(number box_side) {
	IBaseInteraction<number>::set_box_side(box_side);
	_TSP_inter.set_box_side(box_side);

	if(_walls && box_side < 2*_wall_distance) throw oxDNAException("The box side (%lf) is too small to contain walls separated by %lf", box_side, _wall_distance);
}

template<typename number>
void GraftedInteraction<number>::set_box(BaseBox<number> *box) {
	IBaseInteraction<number>::set_box(box);
	_TSP_inter.set_box(box);

	if(_walls && box->box_sides().x < 2*_wall_distance) throw oxDNAException("The box side (%lf) is too small to contain walls separated by %lf", box->box_sides().x, _wall_distance);
}

template<typename number>
void GraftedInteraction<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	// we need to parse the file storing the anchors' positions
	FILE *af = fopen(_anchor_file, "r");
	if(af == NULL) throw oxDNAException("The anchor file '%s' is not readable", _anchor_file);
	std::vector<LR_vector<number> > anchors;
	double ax, ay, az;
	while(fscanf(af, "%lf %lf %lf\n", &ax, &ay, &az) == 3) {
		LR_vector<number> apos(ax, ay, az);
		// we put the anchor point's position exactly on the colloid's surface
		apos.normalize();
		apos *= _colloid_sigma * 0.5;
		anchors.push_back(apos);
	}
	fclose(af);

	if((int)anchors.size() != _N_arms) throw oxDNAException("The anchor file '%s' contains only %d positions while, according to the topology, there should be %d", anchors.size(), _N_arms);

	particles[0] = new Colloid<number>(anchors);

	for(int i = 0; i < _N_arms*_N_per_arm; i++) {
		particles[i+1] = new TSPParticle<number>();
	}
}

template<typename number>
void GraftedInteraction<number>::read_topology(int N, int *N_strands, BaseParticle<number> **particles) {
	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[512];
	topology.getline(line, 512);
	topology.close();
	int res = sscanf(line, "%*d %d %d\n", &_N_arms, &_N_per_arm);
	if(res != 2) throw oxDNAException("Incorrect topology");
	*N_strands = _N_arms + 1;
	if((_N_arms*_N_per_arm+1) != N) throw oxDNAException("Incoherent topology: the total number of particles should be equal to the number of arms (%d) times the number of particles in an arm (%d) plus one", _N_arms, _N_per_arm);

	allocate_particles(particles, N);

	particles[0]->type = P_COLLOID;
	particles[0]->index = 0;
	particles[0]->strand_id = 0;

	int attractive_from = (int) round(_N_per_arm * (1. - _alpha));
	for(int i = 1; i < N; i ++) {
		TSPParticle<number> *p = (TSPParticle<number> *) particles[i];
		p->index = i;
		p->strand_id = (i-1) / _N_per_arm + 1;
		p->set_arm((i-1) / _N_per_arm);

		int rel_index = (i-1) % _N_per_arm;
		p->type = (rel_index >= attractive_from) ? P_B : P_A;

		p->add_bonded_neigh((TSPParticle<number> *) particles[0]);
		if(rel_index != 0) {
			TSPParticle<number> *p_n = (TSPParticle<number> *) particles[i-1];
			p->add_bonded_neigh(p_n);
		}
	}
}

template<typename number>
number GraftedInteraction<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->is_bonded(q)) return pair_interaction_bonded(p, q, r, update_forces);
	else return pair_interaction_nonbonded(p, q, r, update_forces);
}

template<typename number>
number GraftedInteraction<number>::_wall_interaction(BaseParticle<number> *p, bool update_forces) {
	number energy = 0.;

	number z = p->get_abs_pos(this->_box_side).z;
	// we exert a constant force on monomers that are on the wrong side of the wall
	// this is useful to compress initial configurations
	if(z <= -(_wall_distance-0.5) && update_forces) p->force.z += 10.;
	else if(z >= (_wall_distance-0.5) && update_forces) p->force.z -= 10.;
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

template<typename number>
number GraftedInteraction<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	number energy = 0.f;
	if(q == P_VIRTUAL) {
		if(_walls) energy += _wall_interaction(p, update_forces);
		TSPParticle<number> *TSPp = (TSPParticle<number> *) p;
		for(typename set<TSPParticle<number> *>::iterator it = TSPp->bonded_neighs.begin(); it != TSPp->bonded_neighs.end(); it++) {
			energy += pair_interaction_bonded(p, *it, r, update_forces);
		}
		return energy;
	}
	else if(p->type != P_COLLOID && q->type != P_COLLOID) return _TSP_inter.pair_interaction_bonded(p, q, r, update_forces);

	if(p->index < q->index && update_forces) return energy;

	// if we are here then we have to compute the interaction between a regular monomer and the colloid.
	// we have two possibilities: either the monomer is at the beginning of the chain or it's not.
	// if it's at the beginning of the chain then we also have to add the FENE
	BaseParticle<number> *monomer = (p->type == P_COLLOID) ? q : p;
	BaseParticle<number> *colloid = (p->type == P_COLLOID) ? p : q;
	if((monomer->index-1) % _N_per_arm == 0) {
		int anchor_id = monomer->strand_id - 1;
		LR_vector<number> apos_rel = colloid->int_centers[anchor_id];
		LR_vector<number> apos = colloid->pos + apos_rel;
		LR_vector<number> r_ma = apos.minimum_image(monomer->pos, this->_box_side);

		number sqr_r_ma = r_ma.norm();

		if(sqr_r_ma > _colloid_sqr_rfene) {
			if(update_forces) throw oxDNAException("The distance between particles %d and %d (%lf) exceeds the FENE distance (%lf)\n", p->index, q->index, sqrt(sqr_r_ma), sqrt(_colloid_sqr_rfene));
			else return 1e10;
		}

		energy += -15. * _colloid_sqr_rfene * log(1. - sqr_r_ma/_colloid_sqr_rfene);

		if(update_forces) {
			LR_vector<number> tmp_force = r_ma * (-30. * _colloid_sqr_rfene / (_colloid_sqr_rfene - sqr_r_ma));

			monomer->force -= tmp_force;
			colloid->force += tmp_force;

			colloid->torque += colloid->orientationT * apos_rel.cross(tmp_force);
		}
	}

	LR_vector<number> computed_r(0, 0, 0);
	if(r == NULL) {
		if (q != P_VIRTUAL && p != P_VIRTUAL) {
			computed_r = q->pos - p->pos;
			r = &computed_r;
		}
	}

	number sqr_r = r->norm();
	// the repulsion is cut off at the minimum
	if(sqr_r < _colloid_monomer_sqr_rep_rcut) {
		number part = pow(_colloid_monomer_sqr_sigma/sqr_r, _colloid_n/2.);
		energy += 4 * (part * (part - 1.)) + 1.;
		if(update_forces) {
			number force_mod = 4 * _colloid_n * part * (2*part - 1) / sqr_r;
			p->force -= *r * force_mod;
			q->force += *r * force_mod;
		}
	}

	return energy;
}

template<typename number>
number GraftedInteraction<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if(p->type != P_COLLOID && q->type != P_COLLOID) return _TSP_inter.pair_interaction_nonbonded(p, q, r, update_forces);
	return 0.f;
}

template<typename number>
void GraftedInteraction<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
void GraftedInteraction<number>::generate_random_configuration(BaseParticle<number> **particles, int N, number box_side) {
	BaseParticle<number> *colloid = particles[0];
	colloid->orientation = Utils::get_random_rotation_matrix_from_angle<number>(M_PI);
	colloid->orientation.orthonormalize();
	colloid->orientationT = colloid->orientation.get_transpose();
	colloid->set_positions();

	colloid->pos = LR_vector<number>(0, 0, 0);
	int p_index = 1;
	for(int na = 0; na < _N_arms; na++) {
		LR_vector<number> apos_dist = colloid->int_centers[na];
		apos_dist.normalize();
		apos_dist *= 1.01;

		LR_vector<number> current = colloid->pos + colloid->int_centers[na];
		int i = 0;
		do {
			BaseParticle<number> *p = particles[p_index];
			current += apos_dist;
			p->pos = current;
			p->orientation = Utils::get_random_rotation_matrix_from_angle<number>(M_PI);
			p->orientation.orthonormalize();

			i++;
			p_index++;
		} while(i < _N_per_arm);
	}
}

template<typename number>
Colloid<number>::Colloid(std::vector<LR_vector<number> > &anchor_poss) : TSPParticle<number>() {
	this->N_int_centers = anchor_poss.size();
	this->int_centers = new LR_vector<number>[this->N_int_centers];

	_base_anchors = anchor_poss;
}

template<typename number>
Colloid<number>::~Colloid() {

}

template<typename number>
void Colloid<number>::set_positions() {
	for(int i = 0; i < this->N_int_centers; i++) this->int_centers[i] = this->orientation * _base_anchors[i];
}

template class Colloid<float>;
template class Colloid<double>;

template class GraftedInteraction<float>;
template class GraftedInteraction<double>;
