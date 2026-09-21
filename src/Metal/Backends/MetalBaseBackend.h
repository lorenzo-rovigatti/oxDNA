/**
 * @file    MetalBaseBackend.h
 * @brief   Base backend for Metal GPU simulations
 *
 * Metal equivalent of CUDABaseBackend
 */

#ifndef METALBASEBACKEND_H_
#define METALBASEBACKEND_H_

#include <Metal/Metal.h>
#include <Foundation/Foundation.h>

#include "../metal_utils/MetalBox.h"
#include "../MetalUtils.h"
#include "../../Observables/BaseObservable.h"

/**
 * @brief Basic simulation backend on Metal. All Metal backends should inherit from this class
 *
 * This class provides basic Metal GPU facilities for molecular dynamics simulations.
 *
 * Configuration options:
 * - Metal_device: Device index to use (default: 0 for first Metal-capable device)
 * - Metal_sort_every: Sort particles every N steps for performance (default: 0, disabled)
 * - threads_per_threadgroup: Number of threads per threadgroup (default: 256)
 */
class MetalBaseBackend {
protected:
    /// Sort particles every _sort_every updates (0 = disabled)
    int _sort_every;

    /// Metal device and command infrastructure
    id<MTLDevice> _device;
    id<MTLCommandQueue> _command_queue;
    id<MTLLibrary> _library;

    /// Kernel configuration
    Metal_kernel_cfg _particles_kernel_cfg;

    /// Simulation box
    MetalBox _h_metal_box;
    id<MTLBuffer> _d_metal_box;

    /// Buffer sizes
    size_t _vec_size, _bonds_size, _orient_size;
    m_number _sqr_verlet_skin;

    /// Particle data buffers
    id<MTLBuffer> _d_poss;           // Positions (m_number4)
    id<MTLBuffer> _d_bonds;          // Bonds (MetalBonds)
    id<MTLBuffer> _d_orientations;   // Orientations (m_quat)
    id<MTLBuffer> _d_list_poss;      // Neighbor list positions

    /// Host-side arrays
    m_number4 *_h_poss;
    MetalBonds *_h_bonds;
    m_quat *_h_orientations;

    /// Sorting buffers
    id<MTLBuffer> _d_buff_poss;
    id<MTLBuffer> _d_buff_orientations;
    id<MTLBuffer> _d_buff_bonds;
    id<MTLBuffer> _d_hindex;
    id<MTLBuffer> _d_sorted_hindex;
    id<MTLBuffer> _d_inv_sorted_hindex;

    /// List update flag
    bool *_d_are_lists_old;

    /// Number of particles
    int _N;

    virtual void _host_to_gpu();
    virtual void _gpu_to_host();
    virtual void _sort_index();
    virtual void _init_metal_kernel_cfgs();

    /**
     * @brief Select Metal device (defaults to first available)
     */
    virtual void _choose_device();

public:
    MetalBaseBackend();
    virtual ~MetalBaseBackend();

    virtual void get_settings(input_file &inp);
    virtual void init_metal();

    /**
     * @brief Get the Metal device being used
     */
    id<MTLDevice> get_device() { return _device; }

    /**
     * @brief Get the command queue
     */
    id<MTLCommandQueue> get_command_queue() { return _command_queue; }
};

#endif /* METALBASEBACKEND_H_ */
