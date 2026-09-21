/**
 * @file    MetalBaseBackend.mm
 * @brief   Implementation of Metal base backend
 */

#include "MetalBaseBackend.h"
#include "../../Utilities/oxDNAException.h"

MetalBaseBackend::MetalBaseBackend() :
    _sort_every(0),
    _device(nil),
    _command_queue(nil),
    _library(nil),
    _d_metal_box(nil),
    _vec_size(0),
    _bonds_size(0),
    _orient_size(0),
    _sqr_verlet_skin(0),
    _d_poss(nil),
    _d_bonds(nil),
    _d_orientations(nil),
    _d_list_poss(nil),
    _h_poss(nullptr),
    _h_bonds(nullptr),
    _h_orientations(nullptr),
    _d_buff_poss(nil),
    _d_buff_orientations(nil),
    _d_buff_bonds(nil),
    _d_hindex(nil),
    _d_sorted_hindex(nil),
    _d_inv_sorted_hindex(nil),
    _d_are_lists_old(nullptr),
    _N(0) {
}

MetalBaseBackend::~MetalBaseBackend() {
    // Clean up Metal resources
    // Note: ARC (Automatic Reference Counting) handles most cleanup automatically
    // but we set to nil explicitly for clarity

    _d_poss = nil;
    _d_bonds = nil;
    _d_orientations = nil;
    _d_list_poss = nil;
    _d_metal_box = nil;
    _d_buff_poss = nil;
    _d_buff_orientations = nil;
    _d_buff_bonds = nil;
    _d_hindex = nil;
    _d_sorted_hindex = nil;
    _d_inv_sorted_hindex = nil;

    if(_h_poss) delete[] _h_poss;
    if(_h_bonds) delete[] _h_bonds;
    if(_h_orientations) delete[] _h_orientations;
    if(_d_are_lists_old) delete[] _d_are_lists_old;

    _command_queue = nil;
    _library = nil;
    _device = nil;
}

void MetalBaseBackend::get_settings(input_file &inp) {
    // Read Metal-specific settings
    getInputInt(&inp, "Metal_sort_every", &_sort_every, 0);

    int tmp;
    if(getInputInt(&inp, "threads_per_threadgroup", &tmp, 0) == KEY_FOUND) {
        if(tmp < 1) {
            throw oxDNAException("'threads_per_threadgroup' must be > 0");
        }
        _particles_kernel_cfg.threads_per_threadgroup = tmp;
    } else {
        // Default to 256 threads per threadgroup (good for Apple GPUs)
        _particles_kernel_cfg.threads_per_threadgroup = 256;
    }
}

void MetalBaseBackend::_choose_device() {
    @autoreleasepool {
        // Get default Metal device (typically the discrete GPU on systems with multiple GPUs)
        _device = MTLCreateSystemDefaultDevice();

        if (!_device) {
            throw oxDNAException("No Metal-capable device found. Metal requires Apple Silicon or compatible GPU.");
        }

        OX_LOG(Logger::LOG_INFO, "Using Metal device: %s", [[_device name] UTF8String]);
        OX_LOG(Logger::LOG_INFO, "Max threads per threadgroup: %lu", (unsigned long)[_device maxThreadsPerThreadgroup].width);
        OX_LOG(Logger::LOG_INFO, "Recommended max working set size: %.2f MB",
               [_device recommendedMaxWorkingSetSize] / 1048576.0);
    }
}

void MetalBaseBackend::init_metal() {
    @autoreleasepool {
        // Choose and initialize Metal device
        _choose_device();

        // Create command queue
        _command_queue = [_device newCommandQueue];
        METAL_CHECK_ERROR(_command_queue, "Failed to create Metal command queue");

        // Load default Metal library (compiled shaders)
        NSError *error = nil;
        _library = [_device newDefaultLibrary];
        if (!_library) {
            NSMutableArray<NSString *> *candidatePaths = [NSMutableArray array];

            NSArray<NSString *> *args = [[NSProcessInfo processInfo] arguments];
            if(args.count > 0) {
                NSString *execPath = args[0];
                if(execPath.length > 0) {
                    NSString *execDir = [execPath stringByDeletingLastPathComponent];
                    if(execDir.length > 0) {
                        [candidatePaths addObject:[execDir stringByAppendingPathComponent:@"shaders.metallib"]];
                    }
                }
            }

            NSString *cwdPath = [[[NSFileManager defaultManager] currentDirectoryPath] stringByAppendingPathComponent:@"shaders.metallib"];
            [candidatePaths addObject:cwdPath];

            NSFileManager *fm = [NSFileManager defaultManager];
            for(NSString *candidate in candidatePaths) {
                if([fm fileExistsAtPath:candidate]) {
                    _library = [_device newLibraryWithFile:candidate error:&error];
                    if(_library) {
                        OX_LOG(Logger::LOG_INFO, "Loaded Metal shader library from %s", [candidate UTF8String]);
                        break;
                    }
                }
            }
        }

        if(!_library) {
            const char *err = (error != nil) ? [[error localizedDescription] UTF8String] : "unknown error";
            throw oxDNAException("Failed to load Metal shader library (tried default/executable-dir/cwd): %s", err);
        }

        OX_LOG(Logger::LOG_INFO, "Metal backend initialized successfully");

        // Initialize kernel configurations
        _init_metal_kernel_cfgs();
    }
}

void MetalBaseBackend::_init_metal_kernel_cfgs() {
    // Calculate threadgroups based on particle count
    if(_N > 0) {
        uint32_t threads = _particles_kernel_cfg.threads_per_threadgroup;
        uint32_t threadgroups = (_N + threads - 1) / threads;

        _particles_kernel_cfg.threadgroups_per_grid = threadgroups;
        _particles_kernel_cfg.total_threads = threadgroups * threads;

        OX_LOG(Logger::LOG_INFO, "Kernel config: %u threadgroups × %u threads = %u total",
               threadgroups, threads, _particles_kernel_cfg.total_threads);
    }
}

void MetalBaseBackend::_host_to_gpu() {
    // Copy particle data from host to GPU
    if(_h_poss && _d_poss) {
        MetalUtils::copy_to_device<m_number4>(_d_poss, _h_poss, _N);
    }

    if(_h_bonds && _d_bonds) {
        MetalUtils::copy_to_device<MetalBonds>(_d_bonds, _h_bonds, _N);
    }

    if(_h_orientations && _d_orientations) {
        MetalUtils::copy_to_device<m_quat>(_d_orientations, _h_orientations, _N);
    }

    // Copy box information
    // Copy box information
    if(_d_metal_box) {
        MetalBox::BoxData data = _h_metal_box.get_box_data();
        MetalUtils::copy_to_device<MetalBox::BoxData>(_d_metal_box, &data, 1);
    }
}

void MetalBaseBackend::_gpu_to_host() {
    // Copy particle data from GPU to host
    if(_h_poss && _d_poss) {
        MetalUtils::copy_from_device<m_number4>(_h_poss, _d_poss, _N);
    }

    if(_h_bonds && _d_bonds) {
        MetalUtils::copy_from_device<MetalBonds>(_h_bonds, _d_bonds, _N);
    }

    if(_h_orientations && _d_orientations) {
        MetalUtils::copy_from_device<m_quat>(_h_orientations, _d_orientations, _N);
    }
}

void MetalBaseBackend::_sort_index() {
    // TODO: Implement particle sorting using Hilbert curve or similar
    // For now, this is a placeholder
    OX_LOG(Logger::LOG_DEBUG, "Particle sorting not yet implemented for Metal backend");
}
