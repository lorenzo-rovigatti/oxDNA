/**
 * @file    MetalDNAInteraction.h
 * @date    Created for Metal backend
 * @author  antigravity
 */

#ifndef METALDNAINTERACTION_H_
#define METALDNAINTERACTION_H_

#include "MetalBaseInteraction.h"
#include "../../Interactions/DNAInteraction.h"

/**
 * @brief Metal implementation of the DNA interaction
 */
class MetalDNAInteraction: public MetalBaseInteraction, public DNAInteraction {
protected:
    id<MTLComputePipelineState> _dna_forces_pso;
    id<MTLComputePipelineState> _dna_forces_edge_nonbonded_pso;
    id<MTLComputePipelineState> _dna_forces_edge_bonded_pso;
    id<MTLComputePipelineState> _init_DNA_strand_ends_pso;

    // Debye-Huckel parameters
    bool _use_debye_huckel = false;
    m_number _salt_concentration = 0;
    bool _debye_huckel_half_charged_ends = false;
    m_number _debye_huckel_lambdafactor = 0;
    m_number _debye_huckel_prefactor = 0;
    m_number _debye_huckel_RHIGH = 0;
    m_number _debye_huckel_RC = 0;
    m_number _debye_huckel_B = 0;
    m_number _minus_kappa = 0;
    
    // Model flags
    bool _use_oxDNA2_coaxial_stacking = false;
    bool _use_oxDNA2_FENE = false;

    // Grooving parameters
    m_number _grooving_k_a = 0;
    m_number _grooving_k_rep = 0;
    m_number _grooving_r_a = 0;
    m_number _grooving_r_rep = 0;

    id<MTLBuffer> _d_is_strand_end;
    
    // Constant buffers for parameters
    // Constant buffers for parameters
    id<MTLBuffer> _d_dna_params;
    id<MTLBuffer> _d_init_args;
    id<MTLCommandQueue> _command_queue;
    
    void _init_strand_ends(id<MTLBuffer> d_bonds);
    void process_dna_force_kernel(id<MTLCommandBuffer> command_buffer, MetalBaseList *list, id<MTLBuffer> poss, id<MTLBuffer> orientations, id<MTLBuffer> forces, id<MTLBuffer> torques, id<MTLBuffer> bonds, id<MTLBuffer> metal_box, id<MTLBuffer> energies);
    
public:
	MetalDNAInteraction();
	virtual ~MetalDNAInteraction();

	void get_settings(input_file &inp) override;
    void check_input_sanity(std::vector<BaseParticle *> &particles) override;
    void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;
	void metal_init(int N, id<MTLDevice> device, id<MTLLibrary> library) override;
	m_number get_metal_rcut() override {
        return (m_number)this->_rcut;
    }

	void compute_forces(MetalBaseList *lists, id<MTLBuffer> d_poss, id<MTLBuffer> d_orientations,                                  id<MTLBuffer> forces,
                                  id<MTLBuffer> torques,
                                  id<MTLBuffer> bonds,
                                  id<MTLBuffer> metal_box,
                                  id<MTLBuffer> energies) override;

	void encode_forces(id<MTLCommandBuffer> command_buffer, MetalBaseList *lists, id<MTLBuffer> d_poss,
	                   id<MTLBuffer> d_orientations, id<MTLBuffer> forces, id<MTLBuffer> torques,
	                   id<MTLBuffer> bonds, id<MTLBuffer> metal_box, id<MTLBuffer> energies = nil) override;

protected:
    void _on_T_update() override;
};

#endif /* METALDNAINTERACTION_H_ */
