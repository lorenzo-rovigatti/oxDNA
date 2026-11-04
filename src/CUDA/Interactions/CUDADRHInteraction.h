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
 * The following options are not supported and make the code exit with an error: use_edge = true, the update of the stress tensor
 * (CUDA_update_stress_tensor_every, CUDA_update_particle_stress_tensor) and FFS.
 */
class CUDADRHInteraction: public CUDABaseInteraction, public DRHInteraction {
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

	/// writes forces, torques and energies to the file specified by CUDA_NA_dump_file (debugging and testing tool)
	void _dump(c_number4 *d_poss, c_number4 *d_forces, c_number4 *d_torques);

	/// device array whose i-th element is 1 if the i-th nucleotide of the topology is RNA and 0 if it is DNA
	int *_d_is_rna = nullptr;

	std::string _dump_file;
	int _dump_call = 0;
	bool _dump_hybrid_only = false;
	int _n_calls = 0;
};

#endif /* CUDADRHINTERACTION_H_ */
