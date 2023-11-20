/*
 * DRHInteraction_relax.cpp
 *
 *  Created on: May 25, 2023
 *      Author: Eryk R.
 * 
 * For relaxing initial configurations of systems containing DNA and RNA.
 */

#include <fstream>

#include "DRHInteraction_relax.h"
#include "../Particles/RNANucleotide.h"
#include "../Particles/DNANucleotide.h"

DRHInteraction_relax::DRHInteraction_relax() :
				DRHInteraction() {
	OX_LOG(Logger::LOG_INFO, "Using unphysical backbone (NA_relax interaction)");
	_constant_force = 0;
	_harmonic_force = 1;
}

DRHInteraction_relax::~DRHInteraction_relax() {
}

number DRHInteraction_relax::_backbone(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!DRHInteraction::DNA2Interaction::_check_bonded_neighbour(&p, &q, compute_r)) {
		return (number) 0.f;
	}

	//(excluding DNA-RNA bonding)
	if(_interaction_type(p,q) == 2) {
		return (number) 0.f;
	}

	if(compute_r) {
		_computed_r = q->pos - p->pos;
	}

	//different rbacks depending on DNA/RNA
	LR_vector rback;
	number rbackmod;
	number rbackr0;
	if(_interaction_type(p,q) == 0) {
		rback = _computed_r + q->int_centers[DNANucleotide::BACK] - p->int_centers[DNANucleotide::BACK];
		rbackmod = rback.module();
		rbackr0 = rbackmod - DNA2Interaction::_fene_r0;
	} 
	else if(_interaction_type(p,q) == 1){
		rback = _computed_r + q->int_centers[RNANucleotide::BACK] - p->int_centers[RNANucleotide::BACK];
		rbackmod = rback.module();
		rbackr0 = rbackmod - model->RNA_FENE_R0;
	}
	// since this is a bonded interaction, _interaction_type should return either 0 or 1, which means that
	// this else should never be reached. However, things may change (i.e. in the future this interaction)
	// could be extended to model chimeric strands, so we put this here to avoid issues
	else {
		throw oxDNAException("Bonded particles %d and %d are not of the same type: cannot relax the strand", p->index, q->index);
	}

	number energy = 0;
	if(_backbone_type == _constant_force) {
		energy = _backbone_k * fabs(rbackr0);
	}
	else if(_backbone_type == _harmonic_force) {
		energy = 0.5 * _backbone_k * rbackr0 * rbackr0;
	}

	if(update_forces) {
		LR_vector force;
		// this does not conserve energy very well -- but there is a discontinuity in the force, so...
		if(_backbone_type == _constant_force) {
			if(rbackr0 > 0) force = rback * (-_backbone_k / rbackmod);
			else force = rback * (_backbone_k / rbackmod);
		}
		// conserves energy about as well as the normal RNAInteraction
		else if(_backbone_type == _harmonic_force) force = rback * (-_backbone_k * rbackr0 / rbackmod);

		p->force -= force;
		q->force += force;

		// we need torques in the reference system of the particle
		if(_interaction_type(p,q) == 0) {
			p->torque -= p->orientationT * p->int_centers[DNANucleotide::BACK].cross(force);
			q->torque += q->orientationT * q->int_centers[DNANucleotide::BACK].cross(force);
		} else if(_interaction_type(p,q) == 1) {
			p->torque -= p->orientationT * p->int_centers[RNANucleotide::BACK].cross(force);
			q->torque += q->orientationT * q->int_centers[RNANucleotide::BACK].cross(force);
		}
	}

	return energy;
}

void DRHInteraction_relax::check_input_sanity(std::vector<BaseParticle *> &particles) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		BaseParticle *p = particles[i];
		if(p->n3 != P_VIRTUAL && p->n3->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n3 neighbor is %d, should be < N = %d)", i, p->n3->index, N);
		}
		if(p->n5 != P_VIRTUAL && p->n5->index >= N) {
			throw oxDNAException("Wrong topology for particle %d (n5 neighbor is %d, should be < N = %d)", i, p->n5->index, N);
		}

		//making sure there are no bonds between DNA and RNA particles
		if(p->n3 != P_VIRTUAL) {
			BaseParticle *q = p->n3;
			if(_interaction_type(p,q) == 2){
				throw oxDNAException("Bonds between DNA and RNA are not allowed (particles %d and %d)", i, p->n3->index);
			}
		}
		if(p->n5 != P_VIRTUAL) {
			BaseParticle *q = p->n5;
			if(_interaction_type(p,q) == 2){
				throw oxDNAException("Bonds between DNA and RNA are not allowed (particles %d and %d)", i, p->n5->index);
			}
		}
	}
}

void DRHInteraction_relax::get_settings(input_file &inp) {
	DRHInteraction::get_settings(inp);

	char tmps[256];
	getInputString(&inp, "relax_type", tmps, 1);
	if(strcmp(tmps, "constant_force") == 0) _backbone_type = _constant_force;
	else if(strcmp(tmps, "harmonic_force") == 0) _backbone_type = _harmonic_force;
	else throw oxDNAException("Error while parsing input file: relax_type '%s' not implemented; use constant_force or harmonic_force", tmps);

	float ftmp;
	if(getInputFloat(&inp, "relax_strength", &ftmp, 0) == KEY_FOUND) {
		_backbone_k = (number) ftmp;
		OX_LOG(Logger::LOG_INFO, "Using spring constant = %f for the RNA_relax interaction", _backbone_k);
	}
	else {
		if (_backbone_type == _harmonic_force) {
			_backbone_k = (number) 32.;
		}
		else {
			_backbone_k = (number) 1.;
		}
		OX_LOG(Logger::LOG_INFO, "Using default strength constant = %f for the RNA_relax interaction", _backbone_k);
	}
}
