/*
 * CUDADNA3Interaction.h
 *
 *  Created on: 13/may/25
 *      Author: lorenzo
 */

#ifndef CUDADNA3INTERACTION_H_
#define CUDADNA3INTERACTION_H_

#include "CUDABaseInteraction.h"
#include "../../Interactions/DNA3Interaction.h"

// We need PODs (plain old data) for constant memory
struct OxDNA3Params {
    static constexpr int D0 = TETRAMER_DIM_A;
    static constexpr int D1 = TETRAMER_DIM_B;
    static constexpr int D2 = TETRAMER_DIM_B;
    static constexpr int D3 = TETRAMER_DIM_A;
    static constexpr int total_size = D0 * D1 * D2 * D3;

    // Strides for row-major layout
    static constexpr int STR0 = D1 * D2 * D3;
    static constexpr int STR1 = D2 * D3;
    static constexpr int STR2 = D3;
    static constexpr int STR3 = 1;

    // Flat data storage
    c_number data[total_size];

    // Element access
    __host__ __device__
    c_number& operator()(int i, int j, int k, int l) {
        return data[i * STR0 + j * STR1 + k * STR2 + l * STR3];
    }

    __host__ __device__
    const c_number& operator()(int i, int j, int k, int l) const {
        return data[i * STR0 + j * STR1 + k * STR2 + l * STR3];
    }
};

/**
 * @brief CUDA implementation of the oxDNA3 model, as provided by DNA3Interaction.
 */
class CUDADNA3Interaction: public CUDABaseInteraction, public DNA3Interaction {
public:
	enum {
		DEBYE_HUCKEL = 7
	};
	CUDADNA3Interaction();
	virtual ~CUDADNA3Interaction();

    uint8_t *_d_particle_types = nullptr; // QUESTION: maybe it is accessed in a coalesced way... try with int
	int *_d_is_strand_end = nullptr;

	void get_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList*lists, c_number4 *d_poss, GPU_quat *d_qorientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox*d_box);

protected:
	std::vector<float> _convert_param_array(const MultiDimArray<TETRAMER_DIM_A, TETRAMER_DIM_B, TETRAMER_DIM_B, TETRAMER_DIM_A> *src, int N_arrays);
	void _on_T_update() override;
	void _init_strand_ends(LR_bonds *d_bonds);
};

#endif /* CUDADNA3INTERACTION_H_ */
