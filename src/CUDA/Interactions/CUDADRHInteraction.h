/*
 * CUDADRHInteraction.h
 */

#ifndef CUDADRHINTERACTION_H_
#define CUDADRHINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DRHInteraction.h"

/**
 * @brief CUDA implementation of the oxNA model (DNA2 + RNA2 + DNA/RNA hybrid interactions), as provided by DRHInteraction.
 *
 * It is selected with interaction_type = NA. Both DNA and RNA strands can be present in the same system.
 *
 * It is compatible with use_edge = true (forces are then computed by looping over the list of interacting pairs). The following options
 * are not supported and make the code exit with an error: the update of the stress tensor (CUDA_update_stress_tensor_every,
 * CUDA_update_particle_stress_tensor) and FFS.
 */
namespace {
struct NAConstants;
}

class CUDADRHInteraction: public CUDABaseInteraction, public DRHInteraction {
	friend struct NAConstants;

public:
	CUDADRHInteraction();
	virtual ~CUDADRHInteraction();

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() override {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) override;

protected:
	void _on_T_update() override;

	/// copies the parameters of the DNA part of the model (oxDNA2) to the GPU
	void _init_dna_constants(int N);
	/// copies the parameters of the RNA part of the model (oxRNA2) to the GPU
	void _init_rna_constants(int N);
	/// copies the parameters of the DNA/RNA hybrid terms to the GPU
	void _init_hybrid_constants(int N);

	/// device array whose i-th element is 1 if the i-th nucleotide of the topology is RNA and 0 if it is DNA
	int *_d_is_rna = nullptr;

	/// testing aid (CUDA_NA_dump_hybrid_only): only the DNA/RNA interactions are computed
	bool _dump_hybrid_only = false;
};

#endif /* CUDADRHINTERACTION_H_ */
