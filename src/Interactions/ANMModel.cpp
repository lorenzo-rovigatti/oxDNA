/*
 * ANMModel.cpp
 *
 *  Created Nov 2026
 *  Author: erik (from the code of jonah and subhajit in DNANMInteraction, RNANMInteraction and ANMInteraction)
 */

#include "ANMModel.h"

#include "../Particles/ANMParticle.h"
#include "../Utilities/oxDNAException.h"
#include "../Utilities/Logger.h"
#include "../model.h"

#include <algorithm>
#include <fstream>

ANMModel::ANMModel(bool torsion) :
				_torsion(torsion) {
	// amino acid - amino acid
	_protein_lj.sigma = 0.35f;
	_protein_lj.rstar = 0.349f;
	_protein_lj.b = 306484596.421f;
	_protein_lj.rcut = 0.352894;
	_protein_lj.sqr_rcut = 0.12454f;
	_protein_lj.stiffness = EXCL_EPS;

	// amino acid - nucleotide backbone
	_backbone_lj.sigma = 0.57f;
	_backbone_lj.rstar = 0.569f;
	_backbone_lj.b = 178699253.5f;
	_backbone_lj.rcut = 0.572934f;
	_backbone_lj.sqr_rcut = 0.3283f;
	_backbone_lj.stiffness = 1.0f;

	// amino acid - nucleotide base
	_base_lj.sigma = 0.36f;
	_base_lj.rstar = 0.359f;
	_base_lj.b = 296866090.f;
	_base_lj.rcut = 0.362897f;
	_base_lj.sqr_rcut = 0.1317f;
	_base_lj.stiffness = 1.0f;

	_protein_nucleotide_sqr_rcut = 3.0625f; // 1.75^2

	// default mass is 1 for everything
	for(int type = 0; type < 25; type++) {
		if(type != 4) {
			_masses[type] = 1.f;
		}
	}
}

void ANMModel::get_settings(input_file &inp, bool parfile_required) {
	std::string parfile;
	getInputString(&inp, "parfile", parfile, parfile_required ? 1 : 0);
	if(parfile != "none") {
		_parameter_file = parfile;
	}
	if(parfile_required && _parameter_file.empty()) {
		throw oxDNAException("This interaction requires a parameter file (parfile) that contains the protein springs");
	}

	if(_torsion) {
		float kb_tmp = -1.f;
		float kt_tmp = -1.f;
		getInputFloat(&inp, "bending_k", &kb_tmp, 0);
		getInputFloat(&inp, "torsion_k", &kt_tmp, 0);
		_k_bend = (number) kb_tmp;
		_k_tor = (number) kt_tmp;

		if((_k_bend < 0 && _k_tor >= 0) || (_k_tor < 0 && _k_bend >= 0)) {
			throw oxDNAException("Bending & Torsion Constants Must BOTH be set globally or in parameter file");
		}
		else if(_k_bend >= 0 && _k_tor >= 0) {
			OX_LOG(Logger::LOG_INFO, "Using Global kb/kt values set in Input File");
		}
		else {
			OX_LOG(Logger::LOG_INFO, "No Global kb/kt values set in Input File, Attempting to read from Parameter File");
			_per_pair_stiffness = true;
		}
	}

	std::string massfile;
	if(getInputString(&inp, "massfile", massfile, 0) == KEY_FOUND) {
		OX_LOG(Logger::LOG_INFO, "Using Provided Massfile");
		_load_massfile(massfile);
	}
	else {
		OX_LOG(Logger::LOG_INFO, "Using Default Masses");
	}
}

void ANMModel::_load_massfile(const std::string &filename) {
	std::ifstream mass_stream(filename, std::ios::in);
	if(!mass_stream.good()) {
		throw oxDNAException("Could Not Load Mass File, Aborting");
	}
	int masstypes;
	mass_stream >> masstypes;
	_masses.clear(); // remove the default masses
	int type;
	number mass;
	while(mass_stream >> type >> mass) {
		_masses[type] = mass;
	}
}

void ANMModel::set_masses(std::vector<BaseParticle *> &particles) {
	for(auto p : particles) {
		auto *aa = dynamic_cast<ANMParticle *>(p);
		if(aa != nullptr) {
			auto it = _masses.find(aa->type);
			if(it == _masses.end()) {
				throw oxDNAException("The mass file does not contain the mass of particle type %d", aa->type);
			}
			aa->mass = it->second;
			aa->invmass = 1.f / aa->mass;
		}
	}
}

void ANMModel::check_cosines(const std::array<double, 4> &cosines) {
	for(auto c : cosines) {
		if(c < -1.0 || c > 1.0) {
			throw oxDNAException("Cos of Angle in Parameter File not in Valid bounds");
		}
	}
}

void ANMModel::read_parameter_file(std::vector<BaseParticle *> &particles, int N_amino_acids) {
	if(_parameter_file.empty()) {
		if(N_amino_acids > 0) {
			throw oxDNAException("The topology contains %d amino acids, but no parameter file was given: set 'parfile' to the file that contains the protein springs", N_amino_acids);
		}
		OX_LOG(Logger::LOG_INFO, "Parfile: NONE, No protein parameters were filled");
		return;
	}

	_parse_parameter_file(particles.size(), [&particles](int key1, int key2) {
		auto *q = dynamic_cast<ANMParticle *>(particles[key1]);
		if(q == nullptr || dynamic_cast<ANMParticle *>(particles[key2]) == nullptr) {
			throw oxDNAException("The parameter file contains a spring between particles %d and %d, but only amino acids can be connected by springs", key1, key2);
		}
		if(key1 < key2) {
			q->add_bonded_neighbor(particles[key2]);
		}
	});
}

void ANMModel::unlink_unbonded_neighbours(std::vector<BaseParticle *> &particles) {
	for(auto p : particles) {
		auto *aa = dynamic_cast<ANMParticle *>(p);
		if(aa == nullptr) {
			continue;
		}
		if(aa->n3 != P_VIRTUAL && !aa->is_bonded(aa->n3)) {
			aa->n3 = P_VIRTUAL;
		}
		if(aa->n5 != P_VIRTUAL && !aa->is_bonded(aa->n5)) {
			aa->n5 = P_VIRTUAL;
		}
	}
}

void ANMModel::parse_parameter_file(int N_particles) {
	if(!_parameter_file.empty()) {
		_parse_parameter_file(N_particles, [](int, int) {});
	}
}

void ANMModel::_parse_parameter_file(int N_particles, const std::function<void(int, int)> &on_spring) {
	_springs.clear();
	_angular.clear();

	std::ifstream parameters(_parameter_file, std::ios::in);
	if(!parameters.good()) {
		throw oxDNAException("ParameterFile Could Not Be Opened");
	}

	// the first line is the number of amino acids. Spring indices are particle indices, so they are not bounded by it
	// when the nucleotides come first in the topology
	int N;
	parameters >> N;

	int key1, key2;
	char potswitch;
	double potential, dist;
	int spring_connection_num = 0;
	while(parameters >> key1 >> key2 >> dist >> potswitch >> potential) {
		if(key1 < 0 || key1 >= N_particles) throw oxDNAException("Invalid Particle ID %d in Parameter File (there are %d particles)", key1, N_particles);
		if(key2 < 0 || key2 >= N_particles) throw oxDNAException("Invalid Particle ID %d in Parameter File (there are %d particles)", key2, N_particles);
		if(dist < 0) throw oxDNAException("Invalid Eq Distance %d in Parameter File", dist);
		if(potswitch != 's') throw oxDNAException("Potential Type %c Not Supported", potswitch);
		if(potential < 0) throw oxDNAException("Spring Constant %f Not Supported", potential);

		on_spring(key1, key2);

		spring_connection_num += 1;

		if(_torsion && key2 - key1 == 1) {
			Angular ang;
			parameters >> ang.cosines[0] >> ang.cosines[1] >> ang.cosines[2] >> ang.cosines[3];
			check_cosines(ang.cosines);
			ang.k_bend = ang.k_tor = -1;
			if(_per_pair_stiffness) {
				number kb, kt;
				parameters >> kb >> kt;
				if(kb < 0 || kt < 0) {
					throw oxDNAException("Invalid pairwise kb/kt Value Declared in Parameter File");
				}
				ang.k_bend = kb;
				ang.k_tor = kt;
			}
			_angular[key1] = ang;
		}

		_springs[std::make_pair(key1, key2)] = Spring { dist, potential };
	}

	if(_torsion) {
		if(spring_connection_num == 1 && N > 2 && _per_pair_stiffness) {
			throw oxDNAException("Invalid Parameter File Format, pairwise kb and kt values incorrectly formatted or missing");
		}
		if(spring_connection_num == 1 && N > 2 && !_per_pair_stiffness) {
			throw oxDNAException("Invalid Parameter File Format, global kb and kt values have been set in Inputfile");
		}
	}
	else if(spring_connection_num == 1 && N > 2) {
		throw oxDNAException("Invalid Parameter File Format, cannot use a DNACT Parameter File");
	}
}

number ANMModel::_repulsive_lj(const RepulsiveLJ &lj, const LR_vector &r, LR_vector &force, bool update_forces) const {
	// this is a bit faster than calling r.norm()
	number rnorm = SQR(r.x) + SQR(r.y) + SQR(r.z);
	number energy = (number) 0;

	if(rnorm < SQR(lj.rcut)) {
		if(rnorm > SQR(lj.rstar)) {
			number rmod = sqrt(rnorm);
			number rrc = rmod - lj.rcut;
			energy = lj.stiffness * lj.b * SQR(SQR(rrc));
			if(update_forces) force = -r * (4 * lj.stiffness * lj.b * CUB(rrc) / rmod);
		}
		else {
			number tmp = SQR(lj.sigma) / rnorm;
			number lj_part = tmp * tmp * tmp;
			energy = 4 * lj.stiffness * (SQR(lj_part) - lj_part);
			if(update_forces) force = -r * (24 * lj.stiffness * (lj_part - 2 * SQR(lj_part)) / rnorm);
		}
	}

	if(update_forces && energy == (number) 0) force.x = force.y = force.z = (number) 0;

	return energy;
}

number ANMModel::spring(BaseParticle *p, BaseParticle *q, const LR_vector &r, bool update_forces) {
	auto it = _springs.find(std::make_pair(std::min(p->index, q->index), std::max(p->index, q->index)));
	if(it == _springs.end()) {
		return (number) 0.f;
	}
	const Spring &spring = it->second;

	number k = spring.k;
	number rinsta = r.module();
	number disp = rinsta - spring.r0;
	number energy = 0.5 * k * SQR(disp);

	if(update_forces) {
		LR_vector force(r);
		force *= (-1.0f * k) * disp / rinsta;

		p->force -= force;
		q->force += force;
	}

	return energy;
}

number ANMModel::excluded_volume(BaseParticle *p, BaseParticle *q, const LR_vector &r, bool update_forces) {
	LR_vector force(0, 0, 0);
	number energy = _repulsive_lj(_protein_lj, r, force, update_forces);

	if(update_forces) {
		p->force -= force;
		q->force += force;
		if(_torsion) {
			p->torque += p->orientationT * (r * 0.5).cross(-force);
			q->torque += q->orientationT * (-r * 0.5).cross(force);
		}
	}

	return energy;
}

number ANMModel::angular(BaseParticle *p, BaseParticle *q, const LR_vector &r, bool update_forces) {
	auto it = _angular.find(p->index);
	if(it == _angular.end()) {
		throw oxDNAException("The parameter file does not contain the bending and torsion parameters of the pair %d-%d", p->index, q->index);
	}
	const Angular &ang = it->second;
	const double &a0 = ang.cosines[0];
	const double &b0 = ang.cosines[1];
	const double &c0 = ang.cosines[2];
	const double &d0 = ang.cosines[3];

	number k_bend = _k_bend;
	number k_tor = _k_tor;
	if(_per_pair_stiffness) {
		k_bend = ang.k_bend;
		k_tor = ang.k_tor;
	}

	LR_vector rij_unit = r;
	rij_unit.normalize();
	LR_vector rji_unit = rij_unit * -1.f;

	LR_vector &a1 = p->orientationT.v1;
	LR_vector &b1 = q->orientationT.v1;
	LR_vector &a3 = p->orientationT.v3;
	LR_vector &b3 = q->orientationT.v3;

	double o1 = rij_unit * a1 - a0;
	double o2 = rji_unit * b1 - b0;
	double o3 = a1 * b1 - c0;
	double o4 = a3 * b3 - d0;

	// bending and torsion
	number energy = k_bend * 0.5f * (SQR(o1) + SQR(o2)) + k_tor * 0.5f * (SQR(o3) + SQR(o4));

	if(update_forces) {
		LR_vector force = -(rji_unit.cross(rij_unit.cross(a1)) * k_bend * o1 - rji_unit.cross(rij_unit.cross(b1)) * k_bend * o2) / r.module();

		p->force -= force;
		q->force += force;

		LR_vector ta = rij_unit.cross(a1) * o1 * k_bend;
		LR_vector tb = rji_unit.cross(b1) * o2 * k_bend;

		p->torque += p->orientationT * ta;
		q->torque += q->orientationT * tb;

		LR_vector torsion = a1.cross(b1) * o3 * k_tor + a3.cross(b3) * o4 * k_tor;

		p->torque += p->orientationT * -torsion;
		q->torque += q->orientationT * torsion;
	}

	return energy;
}

number ANMModel::nucleotide_excluded_volume(BaseParticle *protein, BaseParticle *nuc, const LR_vector &r, bool update_forces) {
	// the BACK and BASE sites have the same indices in DNANucleotide and RNANucleotide
	static const int BACK = 0;
	static const int BASE = 2;

	LR_vector force(0, 0, 0);
	LR_vector r_to_back = r - nuc->int_centers[BACK];
	LR_vector r_to_base = r - nuc->int_centers[BASE];

	LR_vector torquenuc(0, 0, 0);
	auto energy = (number) 0.f;

	if(r_to_back.norm() < _backbone_lj.sqr_rcut) {
		energy = _repulsive_lj(_backbone_lj, r_to_back, force, update_forces);
		if(update_forces) {
			torquenuc = nuc->int_centers[BACK].cross(-force);
			nuc->torque += nuc->orientationT * torquenuc;
			nuc->force -= force;
			protein->force += force;
			if(_torsion) {
				r_to_back.normalize();
				// point of contact on the amino acid relative to its centre
				protein->torque += protein->orientationT * (-r_to_back * _protein_lj.sigma * 0.5f).cross(force);
			}
		}
	}

	if(r_to_base.norm() < _base_lj.sqr_rcut) {
		energy += _repulsive_lj(_base_lj, r_to_base, force, update_forces);
		if(update_forces) {
			torquenuc = nuc->int_centers[BASE].cross(-force);
			nuc->torque += nuc->orientationT * torquenuc;
			nuc->force -= force;
			protein->force += force;
			if(_torsion) {
				r_to_base.normalize();
				protein->torque += protein->orientationT * (-r_to_base * _protein_lj.sigma * 0.5f).cross(force);
			}
		}
	}

	return energy;
}
