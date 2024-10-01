/*
 * ChiralRodExplicit.cpp
 *
 *  Created on: 29/Oct/2013
 *      Author: Flavio
 */

#include "ChiralRodExplicit.h"

template<typename number>
ChiralRodExplicit<number>::ChiralRodExplicit() : BaseInteraction<number, ChiralRodExplicit<number> >() {
	this->_int_map[1] = &ChiralRodExplicit<number>::_chiral_pot;
	_length = (number) - 1.0f;
	_chiral_max = (number) 0.f;
	_chiral_min = (number) 0.f;
	_chiral_alpha = (number) 0.f;
	_chiral_delta = (number) 0.f;
	_chiral_interval = (number) 0.f;
	_sigma_solv = 0.2f;
	_N_solv = -1;
	_N_rods = -1;
	_togrow = false;
}

template<typename number>
ChiralRodExplicit<number>::~ChiralRodExplicit() {

}

template<typename number>
void ChiralRodExplicit<number>::get_settings(input_file &inp) {
	BaseInteraction<number>::get_settings(inp);

	std::string tmps;
	getInputString (&inp, "sim_type", tmps, 1);
	if (tmps.compare("MC") and tmps.compare("MC2")) throw oxDNAException ("Cannot run ChiralRodExplicit with MD");


	float tmpf;
	getInputFloat (&inp, "length", &tmpf, 1);
	_length = (number) tmpf;
	if (_length < 0.) throw oxDNAException ("Cannot run ChiralRod Interaction with negative length");

	getInputFloat (&inp, "chiral_delta", &tmpf, 1);
	_chiral_delta = (number) tmpf;

	getInputFloat (&inp, "chiral_alpha", &tmpf, 1);
	_chiral_alpha = (number) tmpf;

	getInputFloat (&inp, "chiral_interval", &tmpf, 1);
	_chiral_interval = (number) tmpf;

	getInputNumber (&inp, "sigma_solvent", &_sigma_solv, 1);

	getInputBool (&inp, "togrow", &_togrow, 0);
	// this assumes that the length of the cylinders is larger than
	// the delta by quite a lot...
	this->_rcut = 1.001 * (_length + (number)1.f);

}

template<typename number>
void ChiralRodExplicit<number>::init() {

	_chiral_min = _chiral_alpha - _chiral_interval;
	_chiral_max = _chiral_alpha + _chiral_interval;

	OX_LOG(Logger::LOG_INFO, "Initializing ChiralRodExplicit: length %g, chiral_delta = %g, chiral_alpha = %g, chiral_interval = %g (min:%g, max:%g)", _length, _chiral_delta, _chiral_alpha, _chiral_interval, _chiral_min, _chiral_max);

	this->_sqr_rcut = SQR(this->_rcut);
}

template<typename number>
void ChiralRodExplicit<number>::allocate_particles(BaseParticle<number> **particles, int N) {
	for(int i = 0; i < N; i++) particles[i] = new SpheroCylinder<number>(_length);
}

template<typename number>
void ChiralRodExplicit<number>::read_topology(int *N_strands, BaseParticle<number> **particles) {
	*N_strands = N;
	std::ifstream topology(this->_topology_filename, ios::in);
	if(!topology.good()) throw oxDNAException("Can't read topology file '%s'. Aborting", this->_topology_filename);
	char line[4096];
	topology.getline(line, 512);
	int my_N, my_N_rods, my_N_solv;
	sscanf(line, "%d %d %d\n", &my_N, &my_N_rods, &my_N_solv);

	_N_rods = my_N_rods;
	_N_solv = my_N_solv;

	if (my_N != my_N_rods + my_N_solv) throw oxDNAException ("Inconsistent topology file!");

	allocate_particles(particles);

	for (int i = 0; i < N; i ++) {
		SpheroCylinder<number> * p = static_cast<SpheroCylinder<number> *> (particles[i]);
		p->index = i;
		p->strand_id = i;
		if (i < _N_rods) p->type = 0;  // rods colloids
		else p->type = 1;              // solvent colloids

		if (p->type == 1) {
			if (_togrow) p->length = 1.e-6;
			else p->length = _sigma_solv;
		}
	}

	OX_LOG(Logger::LOG_INFO, "             ChiralRodExplicit: N_rods %d, N_solvent %d, sigma_solvent %g", _N_rods, _N_solv, _sigma_solv);

}

template<typename number>
number ChiralRodExplicit<number>::pair_interaction(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

template<typename number>
number ChiralRodExplicit<number>::pair_interaction_bonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	return (number) 0.f;
}

template<typename number>
number ChiralRodExplicit<number>::pair_interaction_nonbonded(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {

	LR_vector<number> my_r(0, 0, 0);
	SpheroCylinder<number> * pp = NULL, * qq = NULL;
	// the algorithm is not symmetric, so we have to compute things always in the same order
	if (p->index < q->index) {
		pp = static_cast< SpheroCylinder<number> *> (p);
		qq = static_cast< SpheroCylinder<number> *> (q);
		if (r == NULL) my_r = this->_box->min_image (pp, qq);
		else my_r = *r;
	}
	else {
		pp = static_cast< SpheroCylinder<number> *> (q);
		qq = static_cast< SpheroCylinder<number> *> (p);
		my_r = this->_box->min_image (pp, qq);
	}

	if (pp->index < _N_rods) {  // p is a rod
		if (qq->index < _N_rods) { // q is a rod
			return _chiral_pot (pp, qq, &my_r, update_forces);
		}
		else { // q is a solvent
			return _rod_solv (pp, qq, &my_r, update_forces);
		}
	}
	else { // p is a solvent
		if (qq->index < _N_rods) { // q is a rod
			return _rod_solv (pp, qq, &my_r, update_forces);
			throw oxDNAException ("should never happen %s %s", __FILE__, __LINE__);
		}
		else { // q is a solvent
			return _solv_solv (pp, qq, &my_r, update_forces);
		}
	}
}

template <typename number>
inline number ChiralRodExplicit<number>::_solv_solv(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	//if (r->norm() < ((number)4.f)*_sigma_solv*_sigma_solv) {
	number compare;
	if (_togrow == false) {
		compare = (number)2.f * _sigma_solv;
	}
	else {
		SpheroCylinder<number> * pp = static_cast<SpheroCylinder<number> *> (p);
		SpheroCylinder<number> * qq = static_cast<SpheroCylinder<number> *> (q);
		compare = pp->length + qq->length;
	}
	if (r->norm() < compare*compare) {
		//if (pp->length > 0.3) throw oxDNAException("oh yes");
		this->set_is_infinite(true);
		return (number) 1.e12;
	}
	else {
		return (number) 0.f;
	}
}

template <typename number>
inline number ChiralRodExplicit<number>::_rod_solv(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	// inside here, p is the cylinder and q the solvent, r goes from p to q;
	number rdotu = (_computed_r * p->orientation.v3);
	LR_vector<number> r_on_axis = rdotu * p->orientation.v3;
	LR_vector<number> r_off_axis = _computed_r - r_on_axis;

	if (fabs(rdotu) < (_length / (number)2.f + _sigma_solv)) {
		if (r_off_axis.norm() < (((number) 0.5f + _sigma_solv)*((number) 0.5f + _sigma_solv)) ) {
			this->set_is_infinite(true);
			return (number) 1.e12;
		}
	}

	return (number) 0.f;
}

template<typename number>
number ChiralRodExplicit<number>::_chiral_pot(BaseParticle<number> *p, BaseParticle<number> *q, LR_vector<number> *r, bool update_forces) {
	if (update_forces) throw oxDNAException ("No forces, figlio di ndrocchia");

	number rnorm = _computed_r.norm();
	if (rnorm > this->_sqr_rcut) return (number) 0.f;

	LR_vector<number> my_r = *r;

	// first we compute the distance between the two spherocylinders
	number hlength = _length / (number) 2.f;
	number drdotu1 = my_r * p->orientation.v3;
	number drdotu2 = my_r * q->orientation.v3;
	number u1dotu2 = p->orientation.v3 * q->orientation.v3;

	number mu, lambda;

	number cc = 1. - u1dotu2 * u1dotu2;
	if (cc < 1.e-12) {
		// parallel line segments
		lambda = drdotu1 / 2.;
		mu = -drdotu2 / 2.;
		if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
		if (fabs(mu) > hlength) mu = copysign(hlength, mu);
	}
	else {
		// line segments not parallel
		lambda = ( drdotu1 - u1dotu2 * drdotu2) / cc;
		mu     = (-drdotu2 + u1dotu2 * drdotu1) / cc;
		if (!(fabs (lambda) <= hlength && fabs (mu) <= hlength)) {
			number aux1 = fabs(lambda) - hlength;
			number aux2 = fabs(mu) - hlength;
			if (aux1 > aux2) {
				lambda = copysign (hlength, lambda);
				mu = lambda * u1dotu2 - drdotu2;
				if (fabs(mu) > hlength) mu = copysign (hlength, mu);
			}
			else {
				mu = copysign (hlength, mu);
				lambda = mu * u1dotu2 + drdotu1;
				if (fabs(lambda) > hlength) lambda = copysign(hlength, lambda);
			}
		}
	}
	LR_vector<number> dist = my_r - lambda * p->orientation.v3 + mu * q->orientation.v3;
	number dnorm = dist.norm();

	// early exit if spherocylinders too far away
	if (dnorm > ((number) 1.f + _chiral_delta)*((number) 1.f + _chiral_delta))
		return (number) 0.f;

	// we check if it goes through both rims
	//fabs(lambda + (dist * u1)) < hlength && fabs(mu + (-dist * u2)) < hlength
	bool through_rims = false;
	//if (fabs(lambda + (dist * u1)) < hlength && fabs(mu - (dist * u2)) < hlength) // the line below is equivalent
	if (fabs(mu * u1dotu2 + drdotu1) < hlength && fabs(lambda * u1dotu2 - drdotu2) < hlength)
		through_rims = true;

	if (through_rims) {
		// overlap between inner rims
		if (dnorm <= (number)1.f) {
			//printf ("overlap di tipo 1 %d %d\n", pp->index, qq->index);
			this->set_is_infinite(true);
			return (number) 1.e12;
		}

		// if we got here, it means that we have to handle the chirality.

		// first, we check that the projection of the v1 vector is positive for
		// both particles, for right handed helix
		LR_vector<number> pv3 = p->orientation.v3;
		LR_vector<number> qv3;
		if (u1dotu2 >= (number) 0.f)
			qv3 = q->orientation.v3;
		else
			qv3 = -q->orientation.v3;

		//handedness
		//if ((pp->orientation.v3.cross(qq->orientation.v3)) * dist < (number)0.f) {
		if (pv3.cross(qv3) * dist < 0.f) {
		//if ((pp->orientation.v1 * chiral_dir) > (number) 0.f && (qq->orientation.v1 * chiral_dir) > (number)0.f) {

			// the angle between pv3 and pq3 is what we are looking for
			number angle = LRACOS((pv3 * qv3));

			if (_chiral_min < angle && angle < _chiral_max) {
				return (number) 0.f;
			}
			else {
				// wrong angle
				// they are not happily chiral, so they are overlapping
				this->set_is_infinite(true);
				return (number) 1.e12;
			}
		}
		else { // wrong handedness
			// printf ("overlap di tipo 3 %d %d\n", pp->index, qq->index);
 			this->set_is_infinite(true);
			return (number) 1.e12;
		}
	}
	else {
		// actually compute the cylinder overlap
		number fact = 0.5 / (0.5 + _chiral_delta);
		bool cylinder_overlap = InteractionUtils::cylinder_overlap<number> (p, q, fact * my_r, fact * _length);
		if (cylinder_overlap == true) {
			// printf ("overlap di tipo 4 %d %d\n", pp->index, qq->index);
			this->set_is_infinite(true);
			return (number)1.e12;
		}
		else {
			return (number) 0.f;
		}
	}
}

template<typename number>
void ChiralRodExplicit<number>::check_input_sanity(BaseParticle<number> **particles, int N) {

}

template<typename number>
bool ChiralRodExplicit<number>::generate_random_configuration_overlap (BaseParticle<number> *p, BaseParticle<number> *q) {
	LR_vector<number> dr = this->_box->min_image (p->pos, q->pos);

	number energy = pair_interaction_nonbonded(p, q, &dr, false);

	/*
	if (this->_is_infinite == true) {
		this->set_is_infinite (false);
		dr = this->_box->min_image (q->pos, p->pos);
		_chiral_pot(q, p, &dr, false);
		if (this->_is_infinite == false) throw oxDNAException ("NOT SYMMETRIC");
	}
	*/

	this->set_is_infinite (false);
	return (energy > (number) this->_energy_threshold);
}

template class ChiralRodExplicit<float>;
template class ChiralRodExplicit<double>;
