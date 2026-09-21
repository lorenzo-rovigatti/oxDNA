/**
 * @file    MD_MetalBackend.h
 * @brief   Molecular Dynamics backend for Metal GPU
 *
 * Metal equivalent of MD_CUDABackend
 */

#ifndef MD_METALBACKEND_H_
#define MD_METALBACKEND_H_

#include <Metal/Metal.h>
#include <Foundation/Foundation.h>

#include "MetalBaseBackend.h"
#include "../../Backends/MDBackend.h"
#include "../MetalUtils.h"
#include "../Lists/MetalBaseList.h"
#include "../Interactions/MetalBaseInteraction.h"
#include "../Thermostats/MetalBaseThermostat.h"

/**
 * @brief Numerical precision strategy for the Metal MD backend.
 *
 * Apple GPUs have no hardware double precision (`double` is a compile error in
 * Metal Shading Language), so the CUDA double/mixed kernels cannot be ported
 * directly. The tiers below are the achievable analogues:
 *
 *  - METAL_PREC_FLOAT     : everything in float32 on the GPU. Fastest.
 *  - METAL_PREC_MIXED     : GPU float force evaluation, but positions and the
 *                           velocity-Verlet integration use double-float (df64,
 *                           a pair of float32) arithmetic in the shader. Nearly
 *                           float speed, drift-resistant.
 *  - METAL_PREC_HARDMIXED : GPU float force evaluation, the velocity-Verlet
 *                           integration runs on the CPU in native double over
 *                           the shared (unified-memory) buffers, optionally
 *                           OpenMP-parallel. Most accurate; a per-step CPU pass.
 */
enum MetalPrecision {
    METAL_PREC_FLOAT = 0,
    METAL_PREC_MIXED,
    METAL_PREC_HARDMIXED
};

/**
 * @brief Manages a MD simulation on Apple GPU with Metal
 *
 * This class implements molecular dynamics simulation using the Metal API
 * for Apple Silicon (M-series) GPUs.
 */
class MD_MetalBackend : public MDBackend, public MetalBaseBackend {
protected:
    bool _use_edge;
    bool _any_rigid_body;
    bool _avoid_cpu_calculations;
    MetalPrecision _precision = METAL_PREC_FLOAT;

    /// Compute pipelines for different kernels
    id<MTLComputePipelineState> _first_step_pipeline;
    id<MTLComputePipelineState> _second_step_pipeline;
    id<MTLComputePipelineState> _forces_pipeline;
    id<MTLComputePipelineState> _zero_forces_pipeline;
    id<MTLComputePipelineState> _zero_torques_pipeline;
    id<MTLComputePipelineState> _update_angular_momenta_pipeline;
    id<MTLComputePipelineState> _update_orientations_pipeline;
    // df64 (double-float) integration kernels for backend_precision = mixed
    id<MTLComputePipelineState> _first_step_mixed_pipeline;
    id<MTLComputePipelineState> _second_step_mixed_pipeline;
    id<MTLComputePipelineState> _mixed_sync_vels_pipeline;

    /// Particle velocity and angular momentum buffers
    id<MTLBuffer> _d_vels;      // Linear velocities
    id<MTLBuffer> _d_Ls;        // Angular momenta
    id<MTLBuffer> _d_forces;    // Forces
    id<MTLBuffer> _d_torques;   // Torques
    // df64 (double-float) hi/lo mirrors of positions/velocities/momenta,
    // used by the `mixed` precision tier. _d_poss/_d_vels/_d_Ls remain the
    // float mirrors read by force kernels, lists and thermostats.
    id<MTLBuffer> _d_poss_hi;
    id<MTLBuffer> _d_poss_lo;
    id<MTLBuffer> _d_vels_hi;
    id<MTLBuffer> _d_vels_lo;
    id<MTLBuffer> _d_Ls_hi;
    id<MTLBuffer> _d_Ls_lo;

    m_number4 *_h_vels;
    m_number4 *_h_Ls;
    m_number4 *_h_forces;
    m_number4 *_h_torques;
    
    // Debug energies
    id<MTLBuffer> _d_energies;
    float *_h_energies;

    /// Molecular information for rigid bodies
    std::vector<int> _h_particles_to_mols;
    id<MTLBuffer> _d_particles_to_mols;
    id<MTLBuffer> _d_mol_sizes;
    id<MTLBuffer> _d_molecular_coms;

    /// Sorting buffers
    id<MTLBuffer> _d_buff_vels;
    id<MTLBuffer> _d_buff_Ls;
    id<MTLBuffer> _d_buff_particles_to_mols;

    /// Barostat statistics
    llint _barostat_attempts;
    llint _barostat_accepted;
    int _update_st_every;

    /// Energy output
    bool _print_energy;

    /// Error configuration output
    ObservableOutput *_obs_output_error_conf;
    std::string _error_conf_file;

    /// External forces
    id<MTLBuffer> _d_ext_forces;
    int _max_ext_forces;
    
    /// Metal components
    MetalBaseList *_metal_list;
    MetalBaseInteraction *_metal_interaction;
    MetalBaseThermostat *_metal_thermostat;

    /// Internal methods
    virtual void _gpu_to_host() override;
    virtual void _host_to_gpu() override;
    virtual void _apply_external_forces_changes();

    virtual void _sort_particles();
    virtual void _rescale_molecular_positions(m_number4 new_Ls, m_number4 old_Ls, bool is_reverse_move);
    virtual void _rescale_positions(m_number4 new_Ls, m_number4 old_Ls);

    virtual void _first_step();
    virtual void _encode_first_step(id<MTLCommandBuffer> command_buffer);
    virtual void _encode_first_step_mixed(id<MTLCommandBuffer> command_buffer);
    virtual void _encode_second_step_mixed(id<MTLCommandBuffer> command_buffer);
    virtual void _encode_mixed_sync_vels(id<MTLCommandBuffer> command_buffer);
    virtual void _split_mirrors_to_df();
    virtual void _apply_barostat();
    virtual void _forces_second_step();
    virtual void _set_external_forces();

    virtual void _thermalize();
    virtual void _update_stress_tensor();

    virtual void _init_metal_md_symbols();
    virtual void _create_compute_pipelines();
    virtual void _update_host_buffers_from_particles();
    virtual void _update_particles_from_host_buffers();
    virtual void _zero_force_and_torque_buffers();
    virtual void _encode_zero_force_and_torque(id<MTLCommandBuffer> command_buffer);
    virtual void _encode_second_step(id<MTLCommandBuffer> command_buffer);
    virtual void _sync_forces_torques_from_gpu();
    virtual void _sync_vels_Ls_from_gpu();

    /// true when any particle carries an external force (trap, wall, ...).
    bool _any_ext_forces = false;
    /// Native/mixed path: sync positions/orientations out of the GPU buffers,
    /// evaluate the CPU external forces and add them into _d_forces/_d_torques.
    virtual void _apply_cpu_external_forces_to_gpu();
    /// CPU-integration path: add the CPU external forces on top of whatever is
    /// already in p->force / p->torque (positions are already current on the CPU).
    virtual void _add_external_forces_to_particles();

public:
    MD_MetalBackend();
    virtual ~MD_MetalBackend();

    virtual void get_settings(input_file &inp) override;
    virtual void init() override;

    virtual void sim_step() override;

    virtual void apply_simulation_data_changes() override;
    virtual void apply_changes_to_simulation_data() override;
    
    // Debug method
    float *get_energies() { return _h_energies; }
};

#endif /* MD_METALBACKEND_H_ */
