/*
 * ANMModel.h
 *
 *  Created Nov 2026
 *  Author: erik (from the code of jonah and subhajit in DNANMInteraction, RNANMInteraction and ANMInteraction)
 *
 * The anisotropic network model (ANM) of proteins used by the AC/ACT and the hybrid DNANM/DNACT, RNANM/RNACT and
 * NANM/NACT interactions. It is not an interaction itself: the interactions own an ANMModel and call its potentials.
 *
 * Each amino acid is a single site. The potentials are:
 * - harmonic springs between the pairs listed in the parameter file (these pairs are the protein "bonds");
 * - a quartic-smoothed repulsive Lennard-Jones between all pairs of amino acids;
 * - (ANMT only, "torsion" models) bending and torsion terms between consecutive amino acids of the backbone, which make
 *   the amino acids rigid bodies with an orientation;
 * - a repulsive Lennard-Jones between amino acids and the backbone and base sites of nucleotides.
 */

#ifndef ANMMODEL_H_
#define ANMMODEL_H_

#include "../defs.h"
#include "../Utilities/parse_input/parse_input.h"

#include <array>
#include <functional>
#include <map>
#include <string>
#include <vector>

class BaseParticle;

class ANMModel {
public:
	struct Spring {
		double r0; // equilibrium distance
		double k; // stiffness
	};

	/// bending/torsion parameters of the backbone pair (i, i + 1), stored with the index of i
	struct Angular {
		std::array<double, 4> cosines; // equilibrium values of the four cosines (a0, b0, c0, d0)
		double k_bend, k_tor; // only used if the stiffnesses are read from the parameter file
	};

	/// parameters of a quartic-smoothed repulsive Lennard-Jones
	struct RepulsiveLJ {
		number sigma, rstar, b, rcut, sqr_rcut, stiffness;
	};

	explicit ANMModel(bool torsion);

	/**
	 * @brief Reads the protein-related options (parfile, massfile, bending_k, torsion_k).
	 *
	 * @param parfile_required throw if the parameter file is not set (otherwise "none" or a missing key mean "no proteins")
	 */
	void get_settings(input_file &inp, bool parfile_required);

	/**
	 * @brief Reads the parameter file and registers each spring as a bond between the two amino acids it connects.
	 *
	 * The particles must have been allocated and their indices set. Springs are stored with their smaller index first.
	 *
	 * @param N_amino_acids the number of amino acids in the topology, used for error messages only
	 */
	void read_parameter_file(std::vector<BaseParticle *> &particles, int N_amino_acids);

	/// parses the parameter file without touching the particles (used by the CUDA interactions)
	void parse_parameter_file(int N_particles);

	/// sets the mass of all the amino acids from the mass file (or the default masses)
	void set_masses(std::vector<BaseParticle *> &particles);

	bool torsion() const {
		return _torsion;
	}

	bool has_parameter_file() const {
		return !_parameter_file.empty();
	}

	const std::string &parameter_file() const {
		return _parameter_file;
	}

	/// true if the bending and torsion stiffnesses are read from the parameter file rather than set in the input
	bool per_pair_stiffness() const {
		return _per_pair_stiffness;
	}

	number global_k_bend() const {
		return _k_bend;
	}

	number global_k_tor() const {
		return _k_tor;
	}

	const std::map<std::pair<int, int>, Spring> &springs() const {
		return _springs;
	}

	const std::map<int, Angular> &angular_parameters() const {
		return _angular;
	}

	/// amino acid - amino acid excluded volume
	const RepulsiveLJ &protein_lj() const {
		return _protein_lj;
	}

	/// amino acid - nucleotide backbone excluded volume
	const RepulsiveLJ &backbone_lj() const {
		return _backbone_lj;
	}

	/// amino acid - nucleotide base excluded volume
	const RepulsiveLJ &base_lj() const {
		return _base_lj;
	}

	/// square of the distance between the centres of an amino acid and a nucleotide beyond which they do not interact
	number protein_nucleotide_sqr_rcut() const {
		return _protein_nucleotide_sqr_rcut;
	}

	/// harmonic spring between the bonded amino acids p and q (r = q - p)
	number spring(BaseParticle *p, BaseParticle *q, const LR_vector &r, bool update_forces);

	/// excluded volume between the amino acids p and q (r = q - p)
	number excluded_volume(BaseParticle *p, BaseParticle *q, const LR_vector &r, bool update_forces);

	/// bending and torsion between the consecutive amino acids p and q = p + 1 (r = q - p)
	number angular(BaseParticle *p, BaseParticle *q, const LR_vector &r, bool update_forces);

	/**
	 * @brief Excluded volume between an amino acid and the backbone and base sites of a nucleotide.
	 *
	 * @param r position of the amino acid relative to the centre of the nucleotide
	 */
	number nucleotide_excluded_volume(BaseParticle *protein, BaseParticle *nucleotide, const LR_vector &r, bool update_forces);

	/// the cosines in the parameter file must be in [-1, 1]
	static void check_cosines(const std::array<double, 4> &cosines);

private:
	bool _torsion;
	std::string _parameter_file;
	bool _per_pair_stiffness = false;
	number _k_bend = -1, _k_tor = -1;

	std::map<std::pair<int, int>, Spring> _springs;
	std::map<int, Angular> _angular;
	std::map<int, number> _masses;

	RepulsiveLJ _protein_lj, _backbone_lj, _base_lj;
	number _protein_nucleotide_sqr_rcut;

	void _load_massfile(const std::string &filename);
	void _parse_parameter_file(int N_particles, const std::function<void(int, int)> &on_spring);
	number _repulsive_lj(const RepulsiveLJ &lj, const LR_vector &r, LR_vector &force, bool update_forces) const;
};

#endif /* ANMMODEL_H_ */
