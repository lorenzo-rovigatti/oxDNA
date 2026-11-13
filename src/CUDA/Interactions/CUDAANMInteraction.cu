/*
 * CUDAANMInteraction.cu
 *
 *  Created Nov 2026
 *  Author: erik (replaces CUDADNANMInteraction by jonah)
 */

#include "CUDAANMInteraction.h"

#include "../../Particles/ANMParticle.h"
#include "../../Particles/RNANucleotide.h"
#include "../../Utilities/ConfigInfo.h"
#include "CUDA_DRH.cuh"
#include "CUDA_NA_constants.cuh"
#include "CUDA_ANM.cuh"

#include <algorithm>
#include <cstring>
#include <type_traits>

namespace {

// the parts of the nucleic acid model that a CPU interaction has. Only the overloads that take std::true_type access the
// parameters of the corresponding CPU classes, so they are not instantiated for interactions that do not have them
template<class I>
using has_dna = std::is_base_of<DNA2Interaction, I>;
template<class I>
using has_rna = std::is_base_of<RNA2Interaction, I>;
template<class I>
using has_hybrid = std::is_base_of<DRHInteraction, I>;

template<class I>
void init_dna(I &inter, int N, cuda_drh::HybridModel &hm, cuda_drh::KernelParameters &params, std::true_type) {
	NAConstants::dna(inter, N);
	NAConstants::dna_sites(hm);
	NAConstants::dna_kernel_parameters(inter, params);
}

template<class I>
void init_dna(I &, int, cuda_drh::HybridModel &, cuda_drh::KernelParameters &, std::false_type) {

}

template<class I>
void init_rna(I &inter, int N, cuda_drh::HybridModel &hm, cuda_drh::KernelParameters &params, std::true_type) {
	NAConstants::rna(inter, N);
	NAConstants::rna_sites(hm, inter.RNAInteraction::get_model());
	NAConstants::rna_kernel_parameters(inter, params);
}

template<class I>
void init_rna(I &, int, cuda_drh::HybridModel &, cuda_drh::KernelParameters &, std::false_type) {

}

template<class I>
void init_hybrid(I &inter, int N, const cuda_drh::HybridModel &hm, bool use_edge, int n_forces, std::true_type) {
	NAConstants::hybrid(inter, N, use_edge, n_forces);
}

template<class I>
void init_hybrid(I &, int N, const cuda_drh::HybridModel &hm, bool use_edge, int n_forces, std::false_type) {
	// only the positions of the interaction sites are used (by the amino acid-nucleotide excluded volume)
	NAConstants::copy_hybrid_model(hm, N, use_edge, n_forces);
}

cuda_anm::RepulsiveLJ to_device(const ANMModel::RepulsiveLJ &lj) {
	cuda_anm::RepulsiveLJ res;
	res.sigma = lj.sigma;
	res.rstar = lj.rstar;
	res.b = lj.b;
	res.rcut = lj.rcut;
	res.sqr_rcut = lj.sqr_rcut;
	res.stiffness = lj.stiffness;
	return res;
}

// the kernel parameters are filled by cuda_init and used by compute_forces. Each instantiation of the template has its own copy
template<class I>
cuda_drh::KernelParameters &kernel_parameters() {
	static cuda_drh::KernelParameters params;
	return params;
}

}

template<class CPUInteraction>
CUDAANMInteraction<CPUInteraction>::CUDAANMInteraction(bool torsion) :
				CUDABaseInteraction(),
				CPUInteraction(torsion) {
	_edge_compatible = true;
}

template<class CPUInteraction>
CUDAANMInteraction<CPUInteraction>::~CUDAANMInteraction() {
	if(_d_species != nullptr) CUDA_SAFE_CALL(cudaFree(_d_species));
	if(_d_spring_start != nullptr) CUDA_SAFE_CALL(cudaFree(_d_spring_start));
	if(_d_spring_neighbour != nullptr) CUDA_SAFE_CALL(cudaFree(_d_spring_neighbour));
	if(_d_spring_k != nullptr) CUDA_SAFE_CALL(cudaFree(_d_spring_k));
	if(_d_spring_r0 != nullptr) CUDA_SAFE_CALL(cudaFree(_d_spring_r0));
	if(_d_angular_cosines != nullptr) CUDA_SAFE_CALL(cudaFree(_d_angular_cosines));
	if(_d_angular_k != nullptr) CUDA_SAFE_CALL(cudaFree(_d_angular_k));
}

template<class CPUInteraction>
void CUDAANMInteraction<CPUInteraction>::get_settings(input_file &inp) {
	CPUInteraction::get_settings(inp);

	int sort_every = 0;
	getInputInt(&inp, "CUDA_sort_every", &sort_every, 0);
	if(sort_every > 0) {
		throw oxDNAException("The protein (ANM) interactions do not support CUDA_sort_every > 0 on CUDA");
	}
}

template<class CPUInteraction>
void CUDAANMInteraction<CPUInteraction>::_init_constants(int N) {
	cuda_drh::HybridModel hm;
	memset(&hm, 0, sizeof(hm));
	cuda_drh::KernelParameters &params = kernel_parameters<CPUInteraction>();
	memset(&params, 0, sizeof(params));
	params.average = true;
	params.hybrid_only = false;

	init_dna(*this, N, hm, params, has_dna<CPUInteraction>());
	init_rna(*this, N, hm, params, has_rna<CPUInteraction>());
	init_hybrid(*this, N, hm, _use_edge, _n_forces, has_hybrid<CPUInteraction>());

	const ANMModel &anm = this->_anm;
	cuda_anm::ANMParameters h_params;
	h_params.protein = to_device(anm.protein_lj());
	h_params.backbone = to_device(anm.backbone_lj());
	h_params.base = to_device(anm.base_lj());
	h_params.protein_nucleotide_sqr_rcut = anm.protein_nucleotide_sqr_rcut();
	h_params.torsion = anm.torsion();
	h_params.per_pair_stiffness = anm.per_pair_stiffness();
	h_params.k_bend = anm.global_k_bend();
	h_params.k_tor = anm.global_k_tor();
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_anm::anm_params, &h_params, sizeof(cuda_anm::ANMParameters)));

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_anm::MD_N, &N, sizeof(int)));
	int n_forces = _use_edge ? _n_forces : 1;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(cuda_anm::MD_n_forces, &n_forces, sizeof(int)));
}

template<class CPUInteraction>
void CUDAANMInteraction<CPUInteraction>::_init_arrays(int N) {
	// the topology (read by the CPU interaction) tells the species of each particle
	std::vector<BaseParticle *> &particles = CONFIG_INFO->particles();
	if((int) particles.size() != N) {
		throw oxDNAException("CUDAANMInteraction: the number of particles (%d) is different from N (%d)", (int) particles.size(), N);
	}
	std::vector<int> h_species(N);
	for(BaseParticle *p : particles) {
		if(dynamic_cast<ANMParticle *>(p) != nullptr) {
			h_species[p->index] = cuda_anm::SPECIES_PROTEIN;
		}
		else if(dynamic_cast<RNANucleotide *>(p) != nullptr) {
			h_species[p->index] = cuda_anm::SPECIES_RNA;
		}
		else {
			h_species[p->index] = cuda_anm::SPECIES_DNA;
		}
	}
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_species, sizeof(int) * N));
	CUDA_SAFE_CALL(cudaMemcpy(_d_species, h_species.data(), sizeof(int) * N, cudaMemcpyHostToDevice));

	// this object is not the one whose read_topology is called by the backend, so it parses the parameter file itself
	this->_anm.parse_parameter_file(N);
	const auto &springs = this->_anm.springs();

	// each spring is stored twice, once for each of its particles, sorted by the index of the other particle
	std::vector<std::vector<std::pair<int, const ANMModel::Spring *>>> per_particle(N);
	for(auto &spring : springs) {
		int i = spring.first.first;
		int j = spring.first.second;
		if(h_species[i] != cuda_anm::SPECIES_PROTEIN || h_species[j] != cuda_anm::SPECIES_PROTEIN) {
			throw oxDNAException("The parameter file contains a spring between particles %d and %d, but only amino acids can be connected by springs", i, j);
		}
		per_particle[i].emplace_back(j, &spring.second);
		per_particle[j].emplace_back(i, &spring.second);
	}

	int n_entries = 2 * springs.size();
	std::vector<int> h_start(N + 1);
	std::vector<int> h_neighbour(std::max(n_entries, 1));
	std::vector<c_number> h_k(std::max(n_entries, 1)), h_r0(std::max(n_entries, 1));
	int idx = 0;
	for(int i = 0; i < N; i++) {
		auto &neighs = per_particle[i];
		std::sort(neighs.begin(), neighs.end(), [](const std::pair<int, const ANMModel::Spring *> &a, const std::pair<int, const ANMModel::Spring *> &b) {
			return a.first < b.first;
		});
		h_start[i] = idx;
		for(auto &neigh : neighs) {
			h_neighbour[idx] = neigh.first;
			h_k[idx] = (c_number) neigh.second->k;
			h_r0[idx] = (c_number) neigh.second->r0;
			idx++;
		}
	}
	h_start[N] = idx;

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_spring_start, sizeof(int) * (N + 1)));
	CUDA_SAFE_CALL(cudaMemcpy(_d_spring_start, h_start.data(), sizeof(int) * (N + 1), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_spring_neighbour, sizeof(int) * h_neighbour.size()));
	CUDA_SAFE_CALL(cudaMemcpy(_d_spring_neighbour, h_neighbour.data(), sizeof(int) * h_neighbour.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&_d_spring_k, sizeof(c_number) * h_k.size()));
	CUDA_SAFE_CALL(cudaMemcpy(_d_spring_k, h_k.data(), sizeof(c_number) * h_k.size(), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&_d_spring_r0, sizeof(c_number) * h_r0.size()));
	CUDA_SAFE_CALL(cudaMemcpy(_d_spring_r0, h_r0.data(), sizeof(c_number) * h_r0.size(), cudaMemcpyHostToDevice));

	// bending and torsion parameters of the pair (i, i + 1), stored with the index of i
	std::vector<c_number> h_cosines(4 * N, (c_number) 0.f);
	std::vector<c_number> h_angular_k(2 * N, (c_number) 0.f);
	for(auto &ang : this->_anm.angular_parameters()) {
		int i = ang.first;
		const ANMModel::Angular &a = ang.second;
		for(int k = 0; k < 4; k++) {
			h_cosines[4 * i + k] = (c_number) a.cosines[k];
		}
		h_angular_k[2 * i] = a.k_bend;
		h_angular_k[2 * i + 1] = a.k_tor;
	}
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&_d_angular_cosines, sizeof(c_number) * 4 * N));
	CUDA_SAFE_CALL(cudaMemcpy(_d_angular_cosines, h_cosines.data(), sizeof(c_number) * 4 * N, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<c_number>(&_d_angular_k, sizeof(c_number) * 2 * N));
	CUDA_SAFE_CALL(cudaMemcpy(_d_angular_k, h_angular_k.data(), sizeof(c_number) * 2 * N, cudaMemcpyHostToDevice));

	OX_LOG(Logger::LOG_INFO, "ANM: %d springs on the GPU (%.2f MB)", (int) springs.size(),
			(sizeof(int) * (N + 1 + n_entries) + sizeof(c_number) * (2 * n_entries + 6 * N)) / (double) SQR(1024));
}

template<class CPUInteraction>
void CUDAANMInteraction<CPUInteraction>::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);

	if(_update_st) {
		throw oxDNAException("The protein (ANM) interactions do not support the calculation of the stress tensor on CUDA (CUDA_update_stress_tensor_every and CUDA_update_particle_stress_tensor)");
	}

	{
		ScopedLogMute mute("CUDA interaction init"); // the CPU interaction has already printed these messages
		CPUInteraction::init();
	}
	this->read_topology_counts();

	_init_constants(N);
	if(_d_species == nullptr) {
		_init_arrays(N);
	}
}

template<class CPUInteraction>
void CUDAANMInteraction<CPUInteraction>::_on_T_update() {
	cuda_init(_N);
}

template<class CPUInteraction>
void CUDAANMInteraction<CPUInteraction>::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	cuda_drh::KernelParameters &params = kernel_parameters<CPUInteraction>();

	if(!_use_edge) {
		cuda_anm::anm_forces
			<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_species, params,
				_d_spring_start, _d_spring_neighbour, _d_spring_k, _d_spring_r0, _d_angular_cosines, _d_angular_k, d_box);
		CUT_CHECK_ERROR("anm_forces error");

		_dump_forces_if_requested(d_poss, d_forces, d_torques);
		return;
	}

	if(_n_forces == 1) { // we can directly use d_forces and d_torques so that no sum is required
		cuda_anm::anm_forces_edge_nonbonded
			<<<(lists->N_edges - 1) / _launch_cfg.threads_per_block + 1, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, d_forces, d_torques, lists->d_edge_list, lists->N_edges, d_bonds, _d_species, params, d_box);
		CUT_CHECK_ERROR("anm_forces_edge_nonbonded error");
	}
	else { // sum required, somewhat slower
		cuda_anm::anm_forces_edge_nonbonded
			<<<(lists->N_edges - 1) / _launch_cfg.threads_per_block + 1, _launch_cfg.threads_per_block>>>
			(d_poss, d_orientations, _d_edge_forces, _d_edge_torques, lists->d_edge_list, lists->N_edges, d_bonds, _d_species, params, d_box);
		CUT_CHECK_ERROR("anm_forces_edge_nonbonded error");

		_sum_edge_forces_torques(d_forces, d_torques);
	}

	cuda_anm::anm_forces_edge_bonded
		<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>
		(d_poss, d_orientations, d_forces, d_torques, d_bonds, _d_species, params, _d_spring_start, _d_spring_neighbour, _d_spring_k,
			_d_spring_r0, _d_angular_cosines, _d_angular_k, d_box);
	CUT_CHECK_ERROR("anm_forces_edge_bonded error");

	_dump_forces_if_requested(d_poss, d_forces, d_torques);
}

template class CUDAANMInteraction<DNANMInteraction>;
template class CUDAANMInteraction<RNANMInteraction>;
template class CUDAANMInteraction<NANMInteraction>;
