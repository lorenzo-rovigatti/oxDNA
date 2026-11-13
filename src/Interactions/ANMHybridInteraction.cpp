/*
 * ANMHybridInteraction.cpp
 *
 *  Created Nov 2026
 *  Author: erik (from DNANMInteraction and RNANMInteraction by jonah and subhajit)
 */

#include "ANMHybridInteraction.h"

#include "DNA2Interaction.h"
#include "RNAInteraction2.h"
#include "DRHInteraction.h"
#include "../Particles/ANMParticle.h"

template<class NucleicInteraction>
ANMHybridInteraction<NucleicInteraction>::ANMHybridInteraction(bool torsion, const std::string &nucleic_type) :
				NucleicInteraction(),
				_anm(torsion),
				_nucleic_type(nucleic_type) {
	// the nucleic acid terms (0-7) have been registered by NucleicInteraction
	this->_interaction_map[SPRING] = [this](BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
		return this->_term_spring(p, q, compute_r, update_forces);
	};
	this->_interaction_map[PRO_EXC_VOL] = [this](BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
		return this->_term_protein_excluded_volume(p, q, compute_r, update_forces);
	};
	this->_interaction_map[PRO_NUC_EXC_VOL] = [this](BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
		return this->_term_protein_nucleotide_excluded_volume(p, q, compute_r, update_forces);
	};
	if(torsion) {
		this->_interaction_map[PRO_ANG_POT] = [this](BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
			return this->_term_angular(p, q, compute_r, update_forces);
		};
	}
}

template<class NucleicInteraction>
ANMHybridInteraction<NucleicInteraction>::~ANMHybridInteraction() {

}

template<class NucleicInteraction>
void ANMHybridInteraction<NucleicInteraction>::get_settings(input_file &inp) {
	NucleicInteraction::get_settings(inp);
	// the parameter file can be left out (or set to "none") for simulations without proteins, except with ANMT
	_anm.get_settings(inp, _anm.torsion());
}

template<class NucleicInteraction>
void ANMHybridInteraction<NucleicInteraction>::init() {
	NucleicInteraction::init();
}

template<class NucleicInteraction>
void ANMHybridInteraction<NucleicInteraction>::read_topology_counts() {
	_anm_topology = ANMTopology::read(this->_topology_filename, _nucleic_type);
}

template<class NucleicInteraction>
void ANMHybridInteraction<NucleicInteraction>::allocate_particles(std::vector<BaseParticle *> &particles) {
	if(_anm_topology.kinds.size() != particles.size()) {
		read_topology_counts();
		if(_anm_topology.kinds.size() != particles.size()) {
			throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");
		}
	}

	_is_amino_acid.resize(particles.size());
	for(uint i = 0; i < particles.size(); i++) {
		char kind = _anm_topology.kinds[i];
		_is_amino_acid[i] = (kind == 'P');
		particles[i] = (kind == 'P') ? new ANMParticle(_anm.torsion()) : _new_nucleotide(i, kind);
	}
}

template<class NucleicInteraction>
void ANMHybridInteraction<NucleicInteraction>::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
	read_topology_counts();
	if(_anm_topology.N != (int) particles.size()) {
		throw oxDNAException("Number of lines in the configuration file and number of particles in the topology files don't match. Aborting");
	}
	allocate_particles(particles);
	// with ANMT the amino acids are rigid bodies whose n3 and n5 are the previous and next amino acids
	ANMTopology::fill_particles(_anm_topology, particles, _anm.torsion());
	_anm.set_masses(particles);
	*N_strands = _anm_topology.N_strands;
	// protein bonds come from the parameter file only
	_anm.read_parameter_file(particles, _anm_topology.N_amino_acids);
	_anm.unlink_unbonded_neighbours(particles);
}

template<class NucleicInteraction>
void ANMHybridInteraction<NucleicInteraction>::check_input_sanity(std::vector<BaseParticle *> &particles) {
	// the checks of the nucleic acid interactions assume that all the particles are nucleotides
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::pair_interaction(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_is_protein(p) == _is_protein(q) && p->is_bonded(q)) {
		return pair_interaction_bonded(p, q, compute_r, update_forces);
	}
	return pair_interaction_nonbonded(p, q, compute_r, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::pair_interaction_bonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_is_protein(p)) {
		if(q != P_VIRTUAL && _is_protein(q)) {
			return (number) 0.f; // amino acids and nucleotides are never bonded
		}
		return NucleicInteraction::pair_interaction_bonded(p, q, compute_r, update_forces);
	}

	if(q == P_VIRTUAL || !_is_protein(q) || !p->is_bonded(q)) {
		return (number) 0.f;
	}

	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	number energy = _anm.spring(p, q, this->_computed_r, update_forces);
	energy += _anm.excluded_volume(p, q, this->_computed_r, update_forces);
	if(_anm.torsion() && q->index - p->index == 1) {
		energy += _anm.angular(p, q, this->_computed_r, update_forces);
	}
	return energy;
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	bool p_protein = _is_protein(p);
	bool q_protein = _is_protein(q);
	if(!p_protein && !q_protein) {
		return NucleicInteraction::pair_interaction_nonbonded(p, q, compute_r, update_forces);
	}

	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}

	if(p_protein && q_protein) {
		if(this->_computed_r.norm() >= _anm.protein_lj().sqr_rcut) {
			return (number) 0.f;
		}
		return _anm.excluded_volume(p, q, this->_computed_r, update_forces);
	}

	return _protein_nucleotide_excluded_volume(p, q, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::_protein_nucleotide_excluded_volume(BaseParticle *p, BaseParticle *q, bool update_forces) {
	if(this->_computed_r.norm() >= _anm.protein_nucleotide_sqr_rcut()) {
		return (number) 0.f;
	}
	// r is the position of the amino acid relative to the nucleotide
	if(_is_protein(p)) {
		return _anm.nucleotide_excluded_volume(p, q, -this->_computed_r, update_forces);
	}
	return _anm.nucleotide_excluded_volume(q, p, this->_computed_r, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::pair_interaction_term(int name, BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}
	// the nucleic acid terms are only defined between nucleotides
	if(name < SPRING && (_is_protein(p) || _is_protein(q))) {
		return (number) 0.f;
	}
	return BaseInteraction::pair_interaction_term(name, p, q, false, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::_term_spring(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_is_protein(p) || !_is_protein(q) || !p->is_bonded(q)) {
		return (number) 0.f;
	}
	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}
	return _anm.spring(p, q, this->_computed_r, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::_term_protein_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_is_protein(p) || !_is_protein(q)) {
		return (number) 0.f;
	}
	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}
	if(!p->is_bonded(q) && this->_computed_r.norm() >= _anm.protein_lj().sqr_rcut) {
		return (number) 0.f;
	}
	return _anm.excluded_volume(p, q, this->_computed_r, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::_term_protein_nucleotide_excluded_volume(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(_is_protein(p) == _is_protein(q)) {
		return (number) 0.f;
	}
	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}
	return _protein_nucleotide_excluded_volume(p, q, update_forces);
}

template<class NucleicInteraction>
number ANMHybridInteraction<NucleicInteraction>::_term_angular(BaseParticle *p, BaseParticle *q, bool compute_r, bool update_forces) {
	if(!_is_protein(p) || !_is_protein(q) || !p->is_bonded(q) || std::abs(q->index - p->index) != 1) {
		return (number) 0.f;
	}
	if(compute_r) {
		this->_computed_r = this->_box->min_image(p->pos, q->pos);
	}
	// the angular potential is defined for the pair (i, i + 1)
	if(p->index > q->index) {
		return _anm.angular(q, p, -this->_computed_r, update_forces);
	}
	return _anm.angular(p, q, this->_computed_r, update_forces);
}

template class ANMHybridInteraction<DNA2Interaction>;
template class ANMHybridInteraction<RNA2Interaction>;
template class ANMHybridInteraction<DRHInteraction>;
