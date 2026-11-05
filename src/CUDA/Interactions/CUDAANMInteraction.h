/*
 * CUDAANMInteraction.h
 *
 *  Created Nov 2026
 *  Author: erik (replaces CUDADNANMInteraction by jonah)
 */

#ifndef CUDAANMINTERACTION_H_
#define CUDAANMINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNANMInteraction.h"
#include "../../Interactions/RNANMInteraction.h"
#include "../../Interactions/NANMInteraction.h"

namespace {
struct NAConstants;
}

/**
 * @brief CUDA implementation of the protein (ANM/ANMT) + nucleic acid interactions: DNANM/DNACT, RNANM/RNACT and NANM/NACT.
 *
 * CPUInteraction is DNANMInteraction, RNANMInteraction or NANMInteraction. The nucleotides interact through the device code of the
 * NA interaction (CUDA_DRH.cuh), which contains oxDNA2, oxRNA2 and the DNA/RNA hybrid terms, and the amino acids through the code
 * of CUDA_ANM.cuh, which mirrors ANMModel.
 *
 * It requires use_edge = true and CUDA_sort_every = 0 (the default), and does not support the stress tensor.
 */
template<class CPUInteraction>
class CUDAANMInteraction: public CUDABaseInteraction, public CPUInteraction {
	friend struct NAConstants;

public:
	/// @param torsion use the ANMT protein model (DNACT, RNACT, NACT)
	explicit CUDAANMInteraction(bool torsion);
	virtual ~CUDAANMInteraction();

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() override {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) override;

protected:
	void _on_T_update() override;

	void _init_constants(int N);
	void _init_arrays(int N);

	/// species of each particle, by topology index: 0 DNA, 1 RNA, 2 amino acid
	int *_d_species = nullptr;

	// springs in compressed format, see anm_forces_edge_bonded()
	int *_d_spring_start = nullptr;
	int *_d_spring_neighbour = nullptr;
	c_number *_d_spring_k = nullptr;
	c_number *_d_spring_r0 = nullptr;

	// bending/torsion parameters (ANMT only)
	c_number *_d_angular_cosines = nullptr;
	c_number *_d_angular_k = nullptr;
};

using CUDADNANMInteraction = CUDAANMInteraction<DNANMInteraction>;
using CUDARNANMInteraction = CUDAANMInteraction<RNANMInteraction>;
using CUDANANMInteraction = CUDAANMInteraction<NANMInteraction>;

#endif /* CUDAANMINTERACTION_H_ */
